/*
 *        File:        unitio_prmtop.c
 *
 ************************************************************************
 *                            LEAP                                      *
 *                                                                      *
 *                   Copyright (c) 1992, 1995                           *
 *           Regents of the University of California                    *
 *                     All Rights Reserved.                             *
 *                                                                      *
 *  This software provided pursuant to a license agreement containing   *
 *  restrictions on its disclosure, duplication, and use. This software *
 *  contains confidential and proprietary information, and may not be   *
 *  extracted or distributed, in whole or in part, for any purpose      *
 *  whatsoever, without the express written permission of the authors.  *
 *  This notice, and the associated author list, must be attached to    *
 *  all copies, or extracts, of this software. Any additional           *
 *  restrictions set forth in the license agreement also apply to this  *
 *  software.                                                           *
 ************************************************************************
 *
 *        Description:
 *                Input/output routines for PRMTOP and CRD files
 *                this has been separated out to make unitio.c smaller
 *
 */

#include        "basics.h"
#include        "classes.h"
#include        "avl.h"
#include        "defaults.h"
#include        "fortran.h"
#include        "mathop.h"
#include        "sort.h"
#include        "cmap.h"
#include        "unitio.h"
#include        "parser.h" // for PARMLIB GplAllParameters

#define INTFORMAT_C     "%8d"
#define INTFORMAT_F     "I8"
#define DBLFORMAT_C     "%16.8lE"
#define DBLFORMAT_F     "E16.8"
#define LBLFORMAT_C     "%-4s"
#define LBLFORMAT_F     "a4"

/* COMMENT_DIM: annotate a FLAG section with its dimension name and size.
   For 2D/strided arrays, add STRIDE=n manually after the macro.          */
#define COMMENT_DIM(name) \
    if (GDefaults.dPrmtopFormat>1.0) { \
    sprintf(sComment, "%%COMMENT DIMENSION:%s; SIZE:%d;", (#name), (name)); \
    FortranWriteString(sComment); }

/* COMMENT_SIZE: for arrays with no single dimension name (e.g. POINTERS) */
#define COMMENT_SIZE(size) \
    if (GDefaults.dPrmtopFormat>1.0) { \
    sprintf(sComment, "%%COMMENT SIZE:%d;", (size)); \
    FortranWriteString(sComment); }

#define COMMENT_DIM_UNITS(name,units) \
    if (GDefaults.dPrmtopFormat>1.0) { \
    sprintf(sComment, "%%COMMENT DIMENSION:%s; SIZE:%d; UNITS:%s;", (#name), (name), (units)); \
    FortranWriteString(sComment); }

#define COMMENT(str) \
    if (GDefaults.dPrmtopFormat>1.0) { \
    snprintf(sComment,sizeof(sComment), "%%COMMENT %s", (str)); \
    FortranWriteString(sComment); }

        /* RESTRAINTLOOP is used to loop over the RESTRAINTs */
        /* for adding constants to tables of constants */
#define        RESTRAINTLOOP( type, field, indexStart ) { \
int        ii, iiMax, jj = 0; \
    if ( (iiMax = iVarArrayElementCount( uUnit->vaRestraints )) ) { \
            srPRestraint = PVAI( uUnit->vaRestraints, SAVERESTRAINTt, 0 ); \
            for ( ii=0; ii<iiMax; ii++, srPRestraint++ ) { \
                if ( srPRestraint->iType == type ) { \
                            FortranWriteDouble( srPRestraint->field ); \
                            srPRestraint->iParmIndex = indexStart+jj; \
                            jj++; \
                } \
            } \
    } \
}
#define        bPERT_BOND(bp,a1,a2)        (bp && (bAtomFlagsSet(a1,ATOMPERTURB)\
                                || bAtomFlagsSet(a2,ATOMPERTURB)))
#define        bPERT_ANGLE(bp,a1,a2,a3) (bp && (bAtomFlagsSet(a1,ATOMPERTURB) \
                                || bAtomFlagsSet(a2,ATOMPERTURB)\
                                || bAtomFlagsSet(a3,ATOMPERTURB)))
#define        bPERT_TORSION(bp,a1,a2,a3,a4)        (bp && (bAtomFlagsSet(a1,ATOMPERTURB) \
                                || bAtomFlagsSet(a2,ATOMPERTURB)\
                                || bAtomFlagsSet(a3,ATOMPERTURB)\
                                || bAtomFlagsSet(a4,ATOMPERTURB)))


//---------------------------------------------------------------------------------------------
// BondAugmentationFound: test for the existence of bond augmentations.  Return 1 if they are
//                        found, which will ordeer them to be printed.
//
// Arguments:
//   uUnit:    the tleap Unit to save
//---------------------------------------------------------------------------------------------
int BondAugmentationFound(UNIT uUnit)
{
  int i, found;
  STRING sAtom1, sAtom2, sDesc;
  double dKb, dR0, dKpull, dRpull0, dKpress, dRpress0;

  found = 0;
  for (i = 0; i < iParmSetTotalBondParms(uUnit->psParameters); i++) {
    ParmSetBond(uUnit->psParameters, i, sAtom1, sAtom2, &dKb, &dR0, &dKpull, &dRpull0,
		&dKpress, &dRpress0, sDesc);
    if (fabs(dKpull) > 1.0e-8 || fabs(dKpress) > 1.0e-8) {
      found = 1;
      break;
    }
  }

  return found;
}

/*
 *        UnitIOFindAndCountMolecules
 *
 *        The caller must supply a VARARRAY whose elements are (int)s.
 *        uUnit->vaResidues
 *
 *        Loop through the residues list and count the number
 *        of molecules.
 *        For each molecule, count the number of ATOMs and place
 *        that count in a VARARRAY.  Return the VARARRAY, the number
 *        of molecules counted, and the index of the first solvent molecule.
 */
void
UnitIOFindAndCountMolecules(UNIT uUnit)
{
    SAVERESIDUEt *srPRes;
    int iResidues, iAtom, iCount;
    LOOP lSpanning;
    bool bSeenFirstSolvent = FALSE;
    ATOM aAtom;

    if (uUnit->vaAtomsPerMolecule) VarArrayDestroy(&(uUnit->vaAtomsPerMolecule));
    uUnit->vaAtomsPerMolecule = vaVarArrayCreate(sizeof(int));

    /* Clear the ATOMTOUCHED flag on all the ATOMs */

    ContainerResetAllAtomsFlags(CONTAINER_from(uUnit), ATOMTOUCHED);

    /* Get the first RESIDUE */

    srPRes = PVAI(uUnit->vaResidues, SAVERESIDUEt, 0);
    iResidues = iVarArrayElementCount(uUnit->vaResidues);

    /* Loop over all RESIDUES */
    int iMolecule = 0;
    for (int i=0; i < iResidues; i++, srPRes++) {

        /* Search for the next RESIDUE whose first ATOM has not */
        /* been touched */

        iAtom = srPRes->iAtomStartIndex - 1;
        if (iAtom < 0) {
            /* skip empty residue */
            continue;
        }
        aAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, iAtom)->aAtom;
        if (bAtomFlagsSet(aAtom, ATOMTOUCHED)) {
            continue;
        }

        /* Touch all of the ATOMs within the molecule that */
        /* contains the current RESIDUE */

        iMolecule++;
        iCount = 0;
        LOOPOVERALL(aAtom,SPANNINGTREE,aAtom,ATOM,lSpanning) {
            AtomSetFlags(aAtom, ATOMTOUCHED);
            iCount++;
            CONTAINER cParent = cContainerWithin(aAtom);
            if (iObjectType(cParent) == RESIDUEid)
                RESIDUE_from(cParent)->iMolecule = iMolecule;
        }

        if (!bSeenFirstSolvent) {
            if ( GDefaults.dPrmtopFormat > 1.0 ?
                     bResidueFlagsSet(srPRes->rResidue, RESIDUEBULKSOLVENT) :
                     cResidueType(srPRes->rResidue) == RESTYPESOLVENT ) {
                bSeenFirstSolvent = TRUE;
                uUnit->iFirstSolvent = iVarArrayElementCount(uUnit->vaAtomsPerMolecule);
            }
        }
        /* Add the molecule to the molecule ATOM count array */
        VarArrayAdd(uUnit->vaAtomsPerMolecule, (GENP) & iCount);
    }

    if (!bSeenFirstSolvent) {
        uUnit->iFirstSolvent = iVarArrayElementCount(uUnit->vaAtomsPerMolecule);
    }
}




/*
 *      zUnitIOBuildNonBondArrays
 *
 *        Author:        Christian Schafmeister (1991)
 *
 *      Build the two NON-BOND arrays.  One is NxN where N is the total
 *      number of types, and the other is Mx(M+1)/2 where M is the
 *      number of UNIQUE non-bond types.
 *      The NxN (vaNBIndexMatrix) array is square matrix containing
 *      integer indices+1 into the Mx(M+1)/2 (vaNBParameters) array
 *      which contains the NON-BOND A and C parameters.
 *      (vaPNBIndex) will return a pointer to a VarArray which contains
 *      indices into the (vaPNonBonds) array which contains the indices
 *      for the unique atom types (all those with unique parameters).
 */
void
zUnitIOBuildNonBondArrays(UNIT uUnit, VARARRAY * vaPNBIndexMatrix,
                          VARARRAY * vaPNBParameters,
                          VARARRAY * vaPNBIndex, VARARRAY * vaPNonBonds)
{
    VARARRAY vaNBIndex, vaNonBonds, vaPNBEdits;
    int i, j, iNonBonds, iNBIndices, iTemp, iI, iJ, iX, iY;
    int iIndex, iHBondIndex, iNBIndex, iElement, iHybridization;
    double dMass, dPolar, dE, dR, dE14, dR14, dA, dB, dEI, dRI, dEJ, dRJ;
    double dScreenF;
    STRING sType, sTypeI, sTypeJ, sDesc;
    NONBONDt *nbPFirst, *nbPCur, *nbPTemp, *nbPLast;

    /* Build the NON-BONDED arrays */
    /* First reduce the non-bond parameters to the absolute */
    /* minimum, by rolling up parameters with duplicate values */
    /* maintaining a many to one mapping into the */
    /* rolled up non-bond parameter array */


    vaNBIndex = vaVarArrayCreate(sizeof(int));
    vaNonBonds = vaVarArrayCreate(sizeof(NONBONDt));
    iNonBonds = iParmSetTotalAtomParms(uUnit->psParameters);
    iNBIndices = iNonBonds;
    VarArraySetSize(vaNonBonds, iNonBonds);
    VarArraySetSize(vaNBIndex, iNonBonds);
    vaPNBEdits = uUnit->psParameters->vaNBEdits;
    for (i = 0; i < iNonBonds; i++) {
        ParmSetAtom(uUnit->psParameters, i,
                    sType, &dMass, &dPolar, &dE, &dR, &dE14, &dR14,
                    &dScreenF, &iElement, &iHybridization, sDesc);
        PVAI(vaNonBonds, NONBONDt, i)->bCapableOfHBonding =
            bParmSetCapableOfHBonding(uUnit->psParameters, sType);
        PVAI(vaNonBonds, NONBONDt, i)->dE = dE;
        PVAI(vaNonBonds, NONBONDt, i)->dR = dR;
        PVAI(vaNonBonds, NONBONDt, i)->dE14 = dE14;
        PVAI(vaNonBonds, NONBONDt, i)->dR14 = dR14;
        strcpy(PVAI(vaNonBonds, NONBONDt, i)->sType, sType);
        *PVAI(vaNBIndex, int, i) = i;
    }

    /* Now roll up the equivalent NON-BOND parameters */

    if (!GDefaults.iCharmm) {
        nbPLast = PVAI(vaNonBonds, NONBONDt, iNonBonds - 1);
        for (iNBIndex = 0; iNBIndex < iVarArrayElementCount(vaNBIndex);
             iNBIndex++) {
            nbPFirst = PVAI(vaNonBonds, NONBONDt, 0);
            nbPCur =
                PVAI(vaNonBonds, NONBONDt, *PVAI(vaNBIndex, int, iNBIndex));
            if (!nbPCur->bCapableOfHBonding) {
                while (nbPFirst < nbPCur) {
                    if (nbPFirst->dE == nbPCur->dE &&
                            nbPFirst->dR == nbPCur->dR &&
                            nbPFirst->dE14 == nbPCur->dE14 &&
                            nbPFirst->dR14 == nbPCur->dR14 &&
			    CheckTypeNBEdit(nbPFirst->sType, vaPNBEdits) == 0 &&
			    CheckTypeNBEdit(nbPCur->sType, vaPNBEdits) == 0) {
                        for (nbPTemp = nbPCur; nbPTemp < nbPLast; nbPTemp++)
                            *nbPTemp = *(nbPTemp + 1);
                        iNonBonds--;
                        nbPLast--;

                        /* Update the indices into the NON-BOND array */
                        /* Making sure that indices into parameters that */
                        /* follow the one at j will be moved back on */

                        j = nbPFirst - PVAI(vaNonBonds, NONBONDt, 0);
                        *PVAI(vaNBIndex, int, iNBIndex) = j;
                        for (iTemp = iNBIndex + 1; iTemp < iNBIndices;
                             iTemp++) {
                            if (*PVAI(vaNBIndex, int, iTemp) > j)
                                 (*PVAI(vaNBIndex, int, iTemp))--;
                        }
                        break;
                    }
                    nbPFirst++;
                }
            }
        }
    }

    /* Now the first iNonBonds entries of the array vaNonBonds */
    /* contain unique NON-BOND parameters, and the array vaNBIndex */
    /* contains indices into the vaNonBonds array for every */
    /* NON-BOND type */
    /* Change the size of the vaNonBonds array to the number of */
    /* valid non-bond parameters */
    VarArraySetSize(vaNonBonds, iNonBonds);

    /* Build the NxN integer index array */

    *vaPNBIndexMatrix = vaVarArrayCreate(sizeof(int));
    VarArraySetSize((*vaPNBIndexMatrix), iNonBonds * iNonBonds);
    for (i = 0; i < iNBIndices; i++) {
        for (j = 0; j < iNBIndices; j++) {

            /* Calculate the position of the parameters for the */
            /* non-bond interaction i-j within the vaPNBParameters */
            /* array */

            iI = *PVAI(vaNBIndex, int, i);
            iJ = *PVAI(vaNBIndex, int, j);
            iX = MIN(iI, iJ);
            iY = MAX(iI, iJ);
            iIndex = iY * (iY + 1) / 2 + iX + 1;        /* +1 because they are FORTRAN */
            /* style arrays !!!!! */

            /* Check if there is an H-Bond parameter for this */
            /* interaction, if there is, make iIndex = -iHBondIndex */
            /***ve to signify that the index is into the HBOND tables */

            ParmSetAtom(uUnit->psParameters, i, sTypeI, &dMass, &dPolar,
                        &dE, &dR, &dE14, &dR14, &dScreenF, &iElement,
		       	&iHybridization, sDesc);
            ParmSetAtom(uUnit->psParameters, j, sTypeJ, &dMass, &dPolar,
                        &dE, &dR, &dE14, &dR14, &dScreenF, &iElement,
			&iHybridization, sDesc);
            iHBondIndex =
                iParmSetFindHBond(uUnit->psParameters, sTypeI, sTypeJ);
            if (iHBondIndex != PARM_NOT_FOUND) {
                VPWARN("Atom pair %s, %s is an H-bond pair\n",sTypeI,sTypeJ);
                iIndex = -(iHBondIndex + 1);
            }
            *PVAI(*vaPNBIndexMatrix, int, iI * iNonBonds + iJ) = iIndex;
        }
    }

    /* Now calculate the A,C parameters for all unique */
    /* NON-BOND interactions */

    *vaPNBParameters = vaVarArrayCreate(sizeof(NONBONDACt));
    VarArraySetSize((*vaPNBParameters), iNonBonds * (iNonBonds + 1) / 2);
    for (j = 0; j < iNonBonds; j++) {
        for (i = 0; i <= j; i++) {
            iX = i;
            iY = j;
            iIndex = iY * (iY + 1) / 2 + iX;
            dEI = PVAI(vaNonBonds, NONBONDt, i)->dE;
            dRI = PVAI(vaNonBonds, NONBONDt, i)->dR;
            dEJ = PVAI(vaNonBonds, NONBONDt, j)->dE;
            dRJ = PVAI(vaNonBonds, NONBONDt, j)->dR;
            MathOpConvertNonBondToAC(dEI, dRI, dEJ, dRJ, &dA, &dB);
	    CheckAgainstNBEdits(vaPNBEdits,
				PVAI(vaNonBonds, NONBONDt, i)->sType,
				PVAI(vaNonBonds, NONBONDt, j)->sType,
				&dA, &dB);
            PVAI(*vaPNBParameters, NONBONDACt, iIndex)->dA = dA;
            PVAI(*vaPNBParameters, NONBONDACt, iIndex)->dB = dB;
            if (GDefaults.iCharmm) {
                dEI = PVAI(vaNonBonds, NONBONDt, i)->dE14;
                dRI = PVAI(vaNonBonds, NONBONDt, i)->dR14;
                dEJ = PVAI(vaNonBonds, NONBONDt, j)->dE14;
                dRJ = PVAI(vaNonBonds, NONBONDt, j)->dR14;
                MathOpConvertNonBondToAC(dEI, dRI, dEJ, dRJ, &dA, &dB);
		CheckAgainstNBEdits(vaPNBEdits,
				    PVAI(vaNonBonds, NONBONDt, i)->sType,
				    PVAI(vaNonBonds, NONBONDt, j)->sType,
				    &dA, &dB);
                PVAI(*vaPNBParameters, NONBONDACt, iIndex)->dA14 = dA;
                PVAI(*vaPNBParameters, NONBONDACt, iIndex)->dB14 = dB;
            }
        }
    }

    /* Return the arrays that refer to the atoms */
    *vaPNonBonds = vaNonBonds;
    *vaPNBIndex = vaNBIndex;
    /* The other two arrays, vaPNBIndexMatrix and vaPNBParameters, */
    /* have already been assigned their return values. */

}



//---------------------------------------------------------------------------------------------
// zUnitIOSaveAmberParmFormat: save an Amber-format topology.
//
// Arguments:
//   uUnit:    the tleap Unit to save
//   fOut:     pointer to file that will be written
//   crdname:  name of the accompanying coordinates file
//   bPolar:   flag for the existence of polarization constants
//   bPert:    flag for the existence of (free energy) perturbation constants
//   bNetcdf:  flag to write a NetCDF coordinates file rather than the standard %12.7f format
//---------------------------------------------------------------------------------------------


void UnitIOSaveAmberParmFormat(UNIT uUnit, char *prmtopName,
                                bool bPolar, bool bPert, bool bNetcdf)
{
    int i, iMax, iIndex;
    LOOP lTemp, lSpan;
    ATOM aAtom, aA, aB, aC, aD;
    int iCount, iBondWith, iBondWithout, iAngleWith, iAngleWithout,
        iTorsionWith, iTorsionWithout, iNumExtra, nBondTypes;
    int iCountPerturbed, iCountBondPerturbed, iCountBondBoundary;
    int iCountAnglePerturbed, iCountAngleBoundary;
    int iCountTorsionPerturbed, iCountTorsionBoundary;
    VARARRAY vaExcludedAtoms, vaExcludedCount, vaNBIndexMatrix, vaNBParameters,
             vaNBIndex, vaNonBonds;
    SAVEBONDt *sbPBond;
    SAVEANGLEt *saPAngle;
    SAVEATOMt *saPAtom;
    SAVETORSIONt *stPTorsion;
    SAVERESTRAINTt *srPRestraint;
    double dMass, dPolar, dR, dKb, dR0, dKpull, dRpull0, dKpress, dRpress0,
           dKt, dT0, dTkub, dRkub, dKp, dP0, dB, dD, dTemp, dScEE, dScNB, dScreenF;
    atomNameStr sAtom1, sAtom2, sAtom3, sAtom4;
    typeStr sType1, sType2, sType;
    int iN, iAtoms, iMaxAtoms, iTemp, iAtom, iCalc14, iProper;
    int iElement, iHybridization, iStart;
    RESIDUE rRes;
    bool bFoundSome;
    char *cPTemp = NULL;
    double dX, dY, dZ, dEpsilon, dRStar, dEpsilon14, dRStar14;
    STRING sDesc;
    IX_REC eResEnt = {NULL};
    IX_DESC iResIx;
    char sVersionHeader[81], sComment[81];
    time_t tp;
    FILE *fOut = NULL;

    TurnOffDisplayerUpdates();

    if (bPolar && GDefaults.iIPOL <= 0) {
        VP0("  Conflict: polarizable prmtop can not have IPOL <= 0.\n");
        VP0("  Please change IPOL in frcmod/parmxx.dat or set default IPOL.\n");
        return;
    } else if (!bPolar && GDefaults.iIPOL > 0) {
        VP0("  Conflict: non-polarizable prmtop can not have IPOL > 0.\n");
        VP0("  Please change IPOL in frcmod/parmxx.dat or set default IPOL.\n");
        return;
    }

    /**** Build excluded atom list *********************************/
    MESSAGE("Building the excluded atom list\n");
    vaExcludedCount = vaVarArrayCreate(sizeof(int));
    vaExcludedAtoms = vaVarArrayCreate(sizeof(int));
    iCountPerturbed = 0;
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
        aAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom;
        if (bAtomFlagsSet(aAtom, ATOMPERTURB)) iCountPerturbed++;
        iCount = 0;
        bFoundSome = FALSE;
        iStart = iVarArrayElementCount(vaExcludedAtoms);
        LOOPOVERALL(aAtom,SPANNINGTREE,aA,ATOM,lSpan) {
            if (aA == aAtom) continue;
            // If the atom is more than three away from the first
            // atom then it is not in the excluded atom list
            if (iAtomBackCount(aA) >= 4) break;
            if (iContainerTempInt(aA) > iContainerTempInt(aAtom)) {
                VarArrayAdd(vaExcludedAtoms, (GENP) &iContainerTempInt(aA));
                bFoundSome = TRUE;
                iCount++;
            }
        }
        if (!bFoundSome) {
            iAtoms = 0;
            VarArrayAdd(vaExcludedAtoms, (GENP) &iAtoms);
            iCount++;
        } else {
            // Sort the part of the VARARRAY just added so that the
            // excluded ATOMs are in ascending order by index
            SortByInteger((GENP) PVAI(vaExcludedAtoms, int, iStart), iCount,
                          sizeof(int), (GENP) PVAI(vaExcludedAtoms, int, iStart), TRUE);
        }
        VarArrayAdd(vaExcludedCount, (GENP) &iCount);
    }

    /**** Mark atom tree classification ***************************/
    // Mark main chain atoms where possible, noting the 
    // number of atoms in the largest residue. keep
    // track of residues which can't be marked.
    VP0("Not Marking per-residue atom chain types.\n"); // FIXME: remove this
    create_index(&iResIx, IX_DUPKEYREC, IX_LEN_CSTRING);
    VP0("Marking per-residue atom chain types.\n");
    iMaxAtoms = 0;
    LOOPOVERALL(uUnit,RESIDUES,rRes,RESIDUE,lTemp) {
        iAtoms = iMarkMainChainAtoms(rRes, FALSE);
        if (iAtoms > 0) MarkSideChainAtoms(rRes);
        if (iAtoms < 0) {
            iAtoms = -iAtoms;
            // Couldn't mark main chains
            strcpy(eResEnt.key, rRes->cHeader.sName);
            if (add_key(&eResEnt, &iResIx) != IX_OK)
                DFATAL("add_key() residue chain\n");
        }
        if (iAtoms > iMaxAtoms) iMaxAtoms = iAtoms;
    }
    first_key(&iResIx);
    i = 1;
    while (next_key(&eResEnt, &iResIx) == IX_OK) {
        if (i) {
            VP0("  (Residues lacking connect0/connect1 - \n");
            VP0("   these don't have chain types marked:\n\n");
            VP0("\tres\ttotal affected\n\n");
            i = 0;
        }
        VP0("\t%s\t%d\n", eResEnt.key, eResEnt.count);
    }
    if (!i) VP0("  )\n");
    destroy_index(&iResIx);

    /******* Build the NON-BOND arrays ******/
    zUnitIOBuildNonBondArrays(uUnit, &vaNBIndexMatrix, &vaNBParameters,
                              &vaNBIndex, &vaNonBonds);

    /**** Fork: NetCDF or Fortran text ***************************/
    if (bNetcdf) {
        UnitIOSaveAmberParmNetcdf(prmtopName, uUnit, bPert, bPolar,
                                   vaExcludedAtoms, vaExcludedCount,
                                   vaNBIndexMatrix, vaNBParameters,
                                   vaNBIndex);
        VarArrayDestroy(&vaNBIndexMatrix);
        VarArrayDestroy(&vaNBParameters);
        VarArrayDestroy(&vaExcludedAtoms);
        VarArrayDestroy(&vaExcludedCount);
        VarArrayDestroy(&vaNBIndex);
        VarArrayDestroy(&vaNonBonds);
        return;
    }

  fOut = FOPENCOMPLAIN( prmtopName, "w" );
  if ( fOut == NULL ) {
    VPFATALEXIT("saveAmberParm: Could not open file: %s\n", prmtopName);
    return;
  }
  FortranFile(fOut);

#if 0
  // Turn on debugging of fortran format output file
  // by sticking comments into the file.
  FortranDebugOn();
#endif

  // coerce format to 1.0, 1.1 or 2.0
  if (GDefaults.dPrmtopFormat<1.0) GDefaults.dPrmtopFormat=1.0;
  else if (GDefaults.dPrmtopFormat>2.0) GDefaults.dPrmtopFormat=2.0;
  else if (GDefaults.dPrmtopFormat!=1.0 && GDefaults.dPrmtopFormat!=2.0) GDefaults.dPrmtopFormat=1.1;

  // -1- Save the title of the UNIT
  FortranDebug("-1-");
  MESSAGE("Saving the name of the UNIT\n");
  FortranFormat(1, "%-80s");
  time(&tp);
  strftime(sComment,sizeof(sComment), "DATE = %m/%d/%y  %H:%M:%S", localtime(&tp));
  snprintf(sVersionHeader,sizeof(sVersionHeader),"%%VERSION  VERSION_STAMP = V%08.3f  %.40s",
        GDefaults.dPrmtopFormat, sComment);
  FortranWriteString(sVersionHeader);
  FortranWriteString("%FLAG TITLE");
  FortranWriteString("%FORMAT(20" LBLFORMAT_F ")");
  FortranWriteString(sContainerName(uUnit));

/////////////////////////////////////////////////
  if (GDefaults.dPrmtopFormat>=2.0) {
      PARMSET psSet;
      int iLines=0;
      FortranFormat(1,"%-80s");
      FortranWriteString("%FLAG FORCE_FIELD_TYPE");
      ParmLibParmSetLoop(GplAllParameters);
      while ( bParmLibNextParmSet(GplAllParameters, &psSet) ) iLines += 2;
      COMMENT_SIZE(iLines)
      FortranWriteString("%FORMAT(A)");
      ParmLibParmSetLoop(GplAllParameters);
      while ( bParmLibNextParmSet(GplAllParameters, &psSet) ) {
          STRING sTemp;
          snprintf(sTemp,sizeof(sTemp),"Filename: %.1024s",psSet->sFname);
          FortranWriteString(sTemp);
          snprintf(sTemp,sizeof(sTemp),"Title: %.1024s",psSet->sTitle);
          FortranWriteString(sTemp);
      }
  }
/////////////////////////////////////////////////

  // NTOTAT->NATOM, NTOTRS->NRES, MUMBND->NUMBND, MUMANG->NUMANG NHB->NPHB
  int NATOM, NTYPES, NBONH, MBONA, NTHETH, MTHETA, NPHIH, MPHIA;
  int NHPARM, NPARM, NNB, NRES, NBONA, NTHETA, NPHIA;
  int NUMBND, NUMANG, NPTRA, NATYP, NPHB, IFPERT, NBPER, NGPER;
  int NDPER, MBPER, MGPER, MDPER, IFBOX, NMXRS, IFCAP, NUMEXTRA;

  NATOM=iVarArrayElementCount(uUnit->vaAtoms);

  NTYPES=iVarArrayElementCount(vaNonBonds);
  int NTTYP = NTYPES * (NTYPES + 1) / 2;  /* LJ parameter array size */

  // Count the number of bonds with hydrogens, and without
  iBondWith = 0;
  iBondWithout = 0;
  for (i = 0; i < iVarArrayElementCount(uUnit->vaBonds); i++) {
    sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, i);
    aA = PVAI(uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom1 - 1)->aAtom;
    aD = PVAI(uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom2 - 1)->aAtom;
    if (bPERT_BOND(bPert, aA, aD)) {
      continue;
    }
    if (iAtomElement(aA) == HYDROGEN || iAtomElement(aD) == HYDROGEN) {
      iBondWith++;
    }
    else {
      iBondWithout++;
    }
  }
  NBONH=iBondWith;
  MBONA=iBondWithout;

  // Count the number of angles with hydrogens, and without
  iAngleWith = 0;
  iAngleWithout = 0;
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAngles); i++) {
    saPAngle = PVAI(uUnit->vaAngles, SAVEANGLEt, i);
    aA = PVAI(uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom1 - 1)->aAtom;
    aB = PVAI(uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom2 - 1)->aAtom;
    aD = PVAI(uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom3 - 1)->aAtom;
    if (bPERT_ANGLE(bPert, aA, aB, aD)) {
      continue;
    }
    if (iAtomElement(aA) == HYDROGEN || iAtomElement(aB) == HYDROGEN ||
        iAtomElement(aD) == HYDROGEN) {
      iAngleWith++;
    }
    else {
      iAngleWithout++;
    }
  }
  NTHETH=iAngleWith;
  MTHETA=iAngleWithout;

  // Count the number of torsions with hydrogens, and without
  iTorsionWith = 0;
  iTorsionWithout = 0;
  for (i = 0; i < iVarArrayElementCount(uUnit->vaTorsions); i++) {
    stPTorsion = PVAI(uUnit->vaTorsions, SAVETORSIONt, i);
    aA = PVAI(uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom1 - 1)->aAtom;
    aB = PVAI(uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom2 - 1)->aAtom;
    aC = PVAI(uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom3 - 1)->aAtom;
    aD = PVAI(uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom4 - 1)->aAtom;
    if (bPERT_TORSION(bPert, aA, aB, aC, aD)) {
      continue;
    }
    if (iAtomElement(aA) == HYDROGEN || iAtomElement(aB) == HYDROGEN ||
        iAtomElement(aC) == HYDROGEN || iAtomElement(aD) == HYDROGEN) {
      iTorsionWith++;
    }
    else {
      iTorsionWithout++;
    }
  }
  NPHIH=iTorsionWith;
  MPHIA=iTorsionWithout;
  NHPARM=0;
  NPARM=0;

  // Write the number of excluded atoms
  NNB=iVarArrayElementCount(vaExcludedAtoms);
  NRES=iVarArrayElementCount(uUnit->vaResidues);

  // Write the number of bonds/angles/torsions without hydrogens
  // PLUS the number of RESTRAINT bonds/angles/torsions
  NBONA=iBondWithout + iUnitRestraintTypeCount(uUnit, RESTRAINTBOND);
  NTHETA=iAngleWithout + iUnitRestraintTypeCount(uUnit, RESTRAINTANGLE);
  NPHIA=iTorsionWithout + iUnitRestraintTypeCount(uUnit, RESTRAINTTORSION);

  // Write the number of unique bond types, angle types, and torsion types.  Add in the
  // number of RESTRAINT bonds/angles/torsion because they will have new parameters.

  NUMBND=iParmSetTotalBondParms(uUnit->psParameters) +
                  iUnitRestraintTypeCount(uUnit, RESTRAINTBOND);
  NUMANG=iParmSetTotalAngleParms(uUnit->psParameters) +
                  iUnitRestraintTypeCount(uUnit, RESTRAINTANGLE);
  NPTRA=iParmSetTotalTorsionParms(uUnit->psParameters) +
                  iParmSetTotalImproperParms(uUnit->psParameters) +
                  iUnitRestraintTypeCount(uUnit, RESTRAINTTORSION);

  // TODO - have different arrays for different restraint types
  if (iVarArrayElementCount(uUnit->vaRestraints)) {
    VP0(" Restraints:  Bond %d  Angle %d  Torsion %d\n",
        iUnitRestraintTypeCount(uUnit, RESTRAINTBOND),
        iUnitRestraintTypeCount(uUnit, RESTRAINTANGLE),
        iUnitRestraintTypeCount(uUnit, RESTRAINTTORSION));
  }
  else {
    VP0(" (no restraints)\n");
  }

  NATYP=iParmSetTotalAtomParms(uUnit->psParameters);
  NPHB=iParmSetTotalHBondParms(uUnit->psParameters);

  IFPERT = bPert ? 1 : 0;

  // Count the number of bonds to be perturbed, and those
  // across the perturbation/non-perturbed boundary
  iCountBondPerturbed = 0;
  iCountBondBoundary = 0;
  for (i = 0; i < iVarArrayElementCount(uUnit->vaBonds); i++) {
    sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, i);
    if ((sbPBond->fFlags & PERTURBED) != 0) {
      iCountBondPerturbed++;
      if ((sbPBond->fFlags & BOUNDARY) != 0) {
        MESSAGE("Boundary pert bond %d-%d\n", sbPBond->iAtom1, sbPBond->iAtom2);
        iCountBondBoundary++;
      }
    }
  }

  MESSAGE("Perturbed bonds: %d\n", iCountBondPerturbed);
  MESSAGE("Perturbed boundary bonds: %d\n", iCountBondBoundary);

  // Count the number of angles to be perturbed, and those on the boundary
  iCountAnglePerturbed = 0;
  iCountAngleBoundary = 0;
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAngles); i++) {
    saPAngle = PVAI(uUnit->vaAngles, SAVEANGLEt, i);
    if ((saPAngle->fFlags & PERTURBED) != 0) {
      iCountAnglePerturbed++;
    }
    if ((saPAngle->fFlags & BOUNDARY) != 0) {
      iCountAngleBoundary++;
    }
  }

  // Count the number of torsions and impropers to be perturbed and those on the boundary
  iCountTorsionPerturbed = 0;
  iCountTorsionBoundary = 0;
  for (i = 0; i < iVarArrayElementCount(uUnit->vaTorsions); i++) {
    stPTorsion = PVAI(uUnit->vaTorsions, SAVETORSIONt, i);
    if ((stPTorsion->fFlags & PERTURBED) != 0) {
      iCountTorsionPerturbed++;
    }
    if ((stPTorsion->fFlags & BOUNDARY) != 0) {
      iCountTorsionBoundary++;
    }
  }
  NBPER=iCountBondPerturbed;
  NGPER=iCountAnglePerturbed;
  NDPER=iCountTorsionPerturbed;
  MBPER=iCountBondPerturbed - iCountBondBoundary;
  MGPER=iCountAnglePerturbed - iCountAngleBoundary;
  MDPER=iCountTorsionPerturbed - iCountTorsionBoundary;

  // Save flag for periodic boundary conditions
  if (bUnitUseBox(uUnit)) IFBOX = bUnitBoxOct(uUnit) ? 2 : 1;
  else IFBOX = 0;

  // Save the number of atoms in the largest residue
  NMXRS=iMaxAtoms;

  // Save flag for cap information
  IFCAP = bUnitUseSolventCap(uUnit) ? 1 : 0;

  // NUMEXTRA
  iNumExtra = 0;
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
    cPTemp = sAtomType(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom);
    if (!strncmp(cPTemp, "EP", 2)) iNumExtra++;
  }
  NUMEXTRA=iNumExtra;

  // NCOPY was POINTER(32). C4Pairwise was POINTER(33)
  //NCOPY but not sure why it is not here. PMEMDCUDA reads it though.

  if (GDefaults.dPrmtopFormat < 2.0) {
#define POINTERS_SIZE 31
      // -2- Save control information
      FortranDebug("-2-");
      MESSAGE("Saving all the main control variables\n");
      FortranFormat(1, "%-80s");
      FortranWriteString("%FLAG POINTERS");
      COMMENT("MEMBERS:NATOM,NTYPES,NBONH,NBONA,NTHETH,MTHETA,NPHIH,MPHIA,NHPARM,")
      COMMENT("NPARM,NNB,NRES,MBONA,NTHETA,NPHIA,NUMBND,NUMANG,NPTRA,NATYP,NPHB,")
      COMMENT("IFPERT,NBPER,NGPER,NDPER,MBPER,MGPER,MDPER,IFBOX,NMXRS,IFCAP,NUMEXTRA;")
      COMMENT_SIZE(POINTERS_SIZE)
      FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
      FortranFormat(10, INTFORMAT_C);
      int pointers[POINTERS_SIZE]= {
              NATOM, NTYPES, NBONH, MBONA, NTHETH, MTHETA, NPHIH, MPHIA,
              NHPARM, NPARM, NNB, NRES, NBONA, NTHETA, NPHIA,
              NUMBND, NUMANG, NPTRA, NATYP, NPHB, IFPERT, NBPER, NGPER,
              NDPER, MBPER, MGPER, MDPER, IFBOX, NMXRS, IFCAP, NUMEXTRA };
      for (i=0; i<POINTERS_SIZE; i++) FortranWriteInt(pointers[i]);
      FortranEndLine();
  }

  // -3-  write out the names of the atoms
  FortranDebug("-3-");
  MESSAGE("Writing the names of the atoms\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG ATOM_NAME");
  COMMENT_DIM(NATOM)
  FortranWriteString("%FORMAT(20" LBLFORMAT_F ")");
  FortranFormat(20, LBLFORMAT_C);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
    cPTemp = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->sName;
    if (strlen(cPTemp) > 4) {
      cPTemp += (strlen(cPTemp) - 4);
    }
    FortranWriteString(cPTemp);
  }
  FortranEndLine();

  // -4- write out the atomic charges
  FortranDebug("-4-");
  MESSAGE("Writing the atomic charges\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG CHARGE");
  COMMENT_DIM_UNITS(NATOM, "sqrt(kcal*angstrom/mol)")
  COMMENT("NOTE:elementary charge multiplied by 18.2223")
  FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
  FortranFormat(5, DBLFORMAT_C);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
    FortranWriteDouble(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->dCharge * ELECTRONTOKCAL);
  }
  FortranEndLine();
  MESSAGE("Writing the atomic numbers\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG ATOMIC_NUMBER");
  COMMENT_DIM(NATOM)
  FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
  FortranFormat(10, INTFORMAT_C);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
    FortranWriteInt(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->iElement);
  }
  FortranEndLine();

  // -5- write out the atomic masses
  FortranDebug("-5-");
  MESSAGE("Writing the atomic masses\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG MASS");
  COMMENT_DIM_UNITS(NATOM, "g/mol")
  FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
  FortranFormat(5, DBLFORMAT_C);
  int iNumAtoms = iVarArrayElementCount(uUnit->vaAtoms);
  double *dGBRad, *dGBScreen;
  dGBRad = (double *)MALLOC(sizeof(double) * iNumAtoms);
  dGBScreen = (double *)MALLOC(sizeof(double) * iNumAtoms);
  for (i = 0; i < iNumAtoms; i++) {
    saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i);
    iIndex = iParmSetFindAtom(uUnit->psParameters, saPAtom->sType);
    ParmSetAtom(uUnit->psParameters, iIndex, sType, &dMass, &dPolar, &dEpsilon, &dRStar,
		&dEpsilon14, &dRStar14, &dScreenF, &iElement, &iHybridization, sDesc);
    FortranWriteDouble(dMass);
    dGBRad[i]=dGBRadiusForAtom(saPAtom, iElement, dMass, (i+1==iNumAtoms) );
    dGBScreen[i]=dGBScreenForElement(iElement);
  }
  FortranEndLine();

  // -6- write out the atomic types
  FortranDebug("-6-");
  MESSAGE("Writing the atomic types\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG ATOM_TYPE_INDEX");
  COMMENT_DIM(NATOM)
  FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
  FortranFormat(10, INTFORMAT_C);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
    iAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->iTypeIndex - 1;
    iTemp = *PVAI(vaNBIndex, int, iAtom);
    FortranWriteInt(iTemp + 1);
  }
  FortranEndLine();

  // -7- write out the starting index into the excluded atom list
  FortranDebug("-7-");
  MESSAGE("Writing the starting index into the excluded atom list\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG NUMBER_EXCLUDED_ATOMS");
  COMMENT_DIM(NATOM)
  FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
  FortranFormat(10, INTFORMAT_C);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
    FortranWriteInt(*PVAI(vaExcludedCount, int, i));
  }
  FortranEndLine();

  // -8- Write the index for the position of the non bond type of each type
  FortranDebug("-8-");
  MESSAGE("writing position of the non bond type of each type\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG NONBONDED_PARM_INDEX");
  COMMENT("DIMENSION:NTYPES,NTYPES;")
  COMMENT("DESC:Entry I,J = atom pair LJ parm index")
  COMMENT_SIZE(NTYPES*NTYPES)
  FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
  FortranFormat(10, INTFORMAT_C);
  for (i = 0; i < iVarArrayElementCount(vaNBIndexMatrix); i++) {
    FortranWriteInt(*PVAI(vaNBIndexMatrix, int, i));
  }
  FortranEndLine();

  // -9- Residue labels
  // Trim the string down to at most 3 characters by
  // taking the last three characters if it is too long
  FortranDebug("-9-");
  MESSAGE("Writing the residue labels\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG RESIDUE_LABEL");
  COMMENT_DIM(NRES)
  FortranWriteString("%FORMAT(20" LBLFORMAT_F ")");
  FortranFormat(20, LBLFORMAT_C);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaResidues); i++) {
    cPTemp = PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->sName;
    if (strlen(cPTemp) > 3) {
      cPTemp += (strlen(cPTemp) - 3);
    }
    FortranWriteString(cPTemp);
  }
  FortranEndLine();

  if (GDefaults.dPrmtopFormat>1.0) {
      MESSAGE("Writing the residue types\n");
      FortranFormat(1, "%-80s");
      FortranWriteString("%FLAG RESIDUE_TYPE");
      COMMENT_DIM(NRES)
      FortranWriteString("%FORMAT(80A1)");
      FortranFormat(80, "%1.1s");
      int iUnknownResTypes=0;
      for (i = 0; i < iVarArrayElementCount(uUnit->vaResidues); i++) {
        cPTemp = PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->sResidueType;
        if (*cPTemp == 0) { cPTemp[0]=RESTYPEUNDEFINED; cPTemp[1]=0; }
        if (*cPTemp == RESTYPEUNDEFINED) iUnknownResTypes++;
        FortranWriteString(cPTemp);
      }
      FortranEndLine();
      if (iUnknownResTypes>0) VPWARN("%d residues have undefined type\n",iUnknownResTypes);

      FortranFormat(1, "%-80s");
      FortranWriteString("%FLAG RESIDUE_MOLECULE");
      COMMENT_DIM(NRES)
      FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
      FortranFormat(10, INTFORMAT_C);
      for (i = 0; i < iVarArrayElementCount(uUnit->vaResidues); i++) {
        int j = PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->rResidue->iMolecule;
        FortranWriteInt(j);
      }
      FortranEndLine();
  }

  // -10- Pointer list for all the residues
  FortranDebug("-10-");


  // -10- Pointer list for all the residues
  FortranDebug("-10-");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG RESIDUE_POINTER");
  COMMENT_DIM(NRES)
  FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
  FortranFormat(10, INTFORMAT_C);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaResidues); i++) {
    FortranWriteInt(PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->iAtomStartIndex);
  }
  FortranEndLine();

  // -11- Force constants for bonds
  FortranDebug("-11-");
  MESSAGE("Writing bond force constants\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG BOND_FORCE_CONSTANT");
  COMMENT_DIM_UNITS(NUMBND, "kcal/mol/angstrom^2")
  FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
  FortranFormat(5, DBLFORMAT_C);
  for (i = 0; i < iParmSetTotalBondParms(uUnit->psParameters); i++) {
    ParmSetBond(uUnit->psParameters, i, sAtom1, sAtom2, &dKb, &dR0, &dKpull, &dRpull0,
                &dKpress, &dRpress0, sDesc);
    FortranWriteDouble(dKb);
  }

  // Write the RESTRAINT constants AND set the index
  // for where the interaction can find its constants
  RESTRAINTLOOP(RESTRAINTBOND, dKx, i + 1);
  FortranEndLine();

  // -12A- Equilibrium bond lengths
  FortranDebug("-12A-");
  MESSAGE("Writing equilibrium bond lengths\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG BOND_EQUIL_VALUE");
  COMMENT_DIM_UNITS(NUMBND, "angstrom")
  FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
  FortranFormat(5, DBLFORMAT_C);
  nBondTypes = iParmSetTotalBondParms(uUnit->psParameters);
  for (i = 0; i < nBondTypes; i++) {
    ParmSetBond(uUnit->psParameters, i, sAtom1, sAtom2, &dKb, &dR0, &dKpull, &dRpull0,
                &dKpress, &dRpress0, sDesc);
    FortranWriteDouble(dR0);
  }

  // Write the bond RESTRAINT constants AND set the index
  // for where the interaction can find its constants
  RESTRAINTLOOP(RESTRAINTBOND, dX0, i + 1);
  FortranEndLine();

  // If (and only if) bond augmentations exist, print them into the topology
  if (BondAugmentationFound(uUnit) == 1) {

    // -12B- Bond pulling adjustments: force constants
    FortranDebug("-12B-");
    MESSAGE("Writing bond pulling adjustments--force constants\n");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG BOND_STIFFNESS_PULL_ADJ");
    COMMENT_DIM_UNITS(NUMBND, "kcal/mol/angstrom^2")
    FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
    FortranFormat(5, DBLFORMAT_C);
    for (i = 0; i < nBondTypes; i++) {
      ParmSetBond(uUnit->psParameters, i, sAtom1, sAtom2, &dKb, &dR0, &dKpull, &dRpull0,
                  &dKpress, &dRpress0, sDesc);
      FortranWriteDouble(dKpull);
    }
    FortranEndLine();

    // -12C- Bond pulling adjustments: equilibrium lengths
    FortranDebug("-12C-");
    MESSAGE("Writing bond pulling adjustments--equilibria\n");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG BOND_EQUIL_PULL_ADJ");
    COMMENT_DIM_UNITS(NUMBND, "angstrom")
    FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
    FortranFormat(5, DBLFORMAT_C);
    for (i = 0; i < nBondTypes; i++) {
      ParmSetBond(uUnit->psParameters, i, sAtom1, sAtom2, &dKb, &dR0, &dKpull, &dRpull0,
                  &dKpress, &dRpress0, sDesc);
      FortranWriteDouble(dRpull0);
    }
    FortranEndLine();

    // -12D- Bond pulling adjustments: equilibrium lengths
    FortranDebug("-12D-");
    MESSAGE("Writing bond compression adjustments--force constants\n");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG BOND_STIFFNESS_PRESS_ADJ");
    COMMENT_DIM_UNITS(NUMBND, "kcal/mol/angstrom^2")
    FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
    FortranFormat(5, DBLFORMAT_C);
    for (i = 0; i < nBondTypes; i++) {
      ParmSetBond(uUnit->psParameters, i, sAtom1, sAtom2, &dKb, &dR0, &dKpull, &dRpull0,
                  &dKpress, &dRpress0, sDesc);
      FortranWriteDouble(dKpress);
    }
    FortranEndLine();

    // -12E- Bond pulling adjustments: equilibrium lengths
    FortranDebug("-12E-");
    MESSAGE("Writing bond compression adjustments--force constants\n");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG BOND_EQUIL_PRESS_ADJ");
    COMMENT_DIM_UNITS(NUMBND, "angstrom")
    FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
    FortranFormat(5, DBLFORMAT_C);
    for (i = 0; i < nBondTypes; i++) {
      ParmSetBond(uUnit->psParameters, i, sAtom1, sAtom2, &dKb, &dR0, &dKpull, &dRpull0,
                  &dKpress, &dRpress0, sDesc);
      FortranWriteDouble(dRpress0);
    }
    FortranEndLine();
  }

  // -13- Force constants for angles
  FortranDebug("-13-");
  MESSAGE("Writing angle force constants\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG ANGLE_FORCE_CONSTANT");
  COMMENT_DIM_UNITS(NUMANG, "kcal/mol/radian^2")
  FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
  FortranFormat(5, DBLFORMAT_C);
  for (i = 0; i < iParmSetTotalAngleParms(uUnit->psParameters); i++) {
    ParmSetAngle(uUnit->psParameters, i, sAtom1, sAtom2, sAtom3, &dKt, &dT0,
                 &dTkub, &dRkub, sDesc);
    FortranWriteDouble(dKt);
  }

  // Write the angle RESTRAINT constants AND set the index
  // for where the interaction can find its constants
  RESTRAINTLOOP(RESTRAINTANGLE, dKx, i + 1);
  FortranEndLine();

  // -14- Equilibrium angle values
  FortranDebug("-14-");

  MESSAGE("Writing equilibrium angle values\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG ANGLE_EQUIL_VALUE");
  COMMENT_DIM_UNITS(NUMANG, "radians")
  FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
  FortranFormat(5, DBLFORMAT_C);
  for (i = 0; i < iParmSetTotalAngleParms(uUnit->psParameters); i++) {
    ParmSetAngle(uUnit->psParameters, i, sAtom1, sAtom2, sAtom3, &dKt, &dT0,
                 &dTkub, &dRkub, sDesc);
    FortranWriteDouble(dT0);
  }

  // Write the angle RESTRAINT constants AND set the index
  // for where the interaction can find its constants
  RESTRAINTLOOP(RESTRAINTANGLE, dX0, i + 1);
  FortranEndLine();

  // -15- Force constants for torsions
  FortranDebug("-15-");
  MESSAGE("Writing torsional force constants\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG DIHEDRAL_FORCE_CONSTANT");
  COMMENT_DIM_UNITS(NPTRA, "kcal/mol")
  FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
  FortranFormat(5, DBLFORMAT_C);
  MESSAGE("There are %d torsions and %d impropers\n",
          iParmSetTotalTorsionParms(uUnit->psParameters),
          iParmSetTotalImproperParms(uUnit->psParameters));
  for (i = 0; i < iParmSetTotalTorsionParms(uUnit->psParameters); i++) {
    ParmSetTorsion(uUnit->psParameters, i, sAtom1, sAtom2, sAtom3, sAtom4,
                   &iN, &dKp, &dP0, &dScEE, &dScNB, sDesc);
    MESSAGE("Torsion %d  %s-%s-%s-%s %d %lf %lf\n", i, sAtom1, sAtom2,
            sAtom3, sAtom4, iN, dKp, dP0);
    FortranWriteDouble(dKp);
  }
  for (i = 0; i < iParmSetTotalImproperParms(uUnit->psParameters); i++) {
    ParmSetImproper(uUnit->psParameters, i, sAtom1, sAtom2, sAtom3, sAtom4,
                    &iN, &dKp, &dP0, sDesc);
    MESSAGE("Improper %d  %s-%s-%s-%s %d %lf %lf\n", i, sAtom1, sAtom2,
            sAtom3, sAtom4, iN, dKp, dP0);
    FortranWriteDouble(dKp);
  }

  // Write the torsion RESTRAINT constants AND set the index
  // for where the interaction can find its constants
  RESTRAINTLOOP(RESTRAINTTORSION, dKx, i + 1);
  FortranEndLine();

  // -16- Periodicity for the dihedral angles
  FortranDebug("-16-");
  MESSAGE("Writing periodicity of torsion interaction\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG DIHEDRAL_PERIODICITY");
  COMMENT_DIM(NPTRA)
  FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
  FortranFormat(5, DBLFORMAT_C);
  for (i = 0; i < iParmSetTotalTorsionParms(uUnit->psParameters); i++) {
    ParmSetTorsion(uUnit->psParameters, i, sAtom1, sAtom2, sAtom3, sAtom4,
                   &iN, &dKp, &dP0, &dScEE, &dScNB, sDesc);
    dTemp = iN;
    FortranWriteDouble(dTemp);
  }
  for (i = 0; i < iParmSetTotalImproperParms(uUnit->psParameters); i++) {
    ParmSetImproper(uUnit->psParameters, i, sAtom1, sAtom2, sAtom3, sAtom4,
                    &iN, &dKp, &dP0, sDesc);
    dTemp = iN;
    FortranWriteDouble(dTemp);
  }

  // Write the torsion RESTRAINT constants AND set the index
  // for where the interaction can find its constants
  RESTRAINTLOOP(RESTRAINTTORSION, dX0, i + 1);
  FortranEndLine();

  // -17- Phase for torsions
  FortranDebug("-17-");

  MESSAGE("Writing phase for torsion interactions\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG DIHEDRAL_PHASE");
  COMMENT_DIM_UNITS(NPTRA, "radians")
  FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
  FortranFormat(5, DBLFORMAT_C);
  for (i = 0; i < iParmSetTotalTorsionParms(uUnit->psParameters); i++) {
    ParmSetTorsion(uUnit->psParameters, i, sAtom1, sAtom2, sAtom3, sAtom4,
                   &iN, &dKp, &dP0, &dScEE, &dScNB, sDesc);
    FortranWriteDouble(dP0);
  }
  for (i = 0; i < iParmSetTotalImproperParms(uUnit->psParameters); i++) {
    ParmSetImproper(uUnit->psParameters, i, sAtom1, sAtom2, sAtom3, sAtom4,
                    &iN, &dKp, &dP0, sDesc);
    FortranWriteDouble(dP0);
  }

  // Write the torsion RESTRAINT constants AND set the index
  // for where the interaction can find its constants
  RESTRAINTLOOP(RESTRAINTTORSION, dN, i + 1);
  FortranEndLine();

  // -17B-
  FortranDebug("-17B-");
  MESSAGE("Writing SCEE_SCALE_FACTOR torsion\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG SCEE_SCALE_FACTOR");
  COMMENT_DIM(NPTRA)
  FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
  FortranFormat(5, DBLFORMAT_C);
  for (i = 0; i < iParmSetTotalTorsionParms(uUnit->psParameters); i++) {
    ParmSetTorsion(uUnit->psParameters, i, sAtom1, sAtom2, sAtom3, sAtom4,
                   &iN, &dKp, &dP0, &dScEE, &dScNB, sDesc);
    if (dScEE < 0.0) {
      dScEE = GDefaults.dSceeScaleFactor;
    }
    FortranWriteDouble(dScEE);
  }
  for (i = 0; i < iParmSetTotalImproperParms(uUnit->psParameters); i++) {
    FortranWriteDouble(0.0);
  }
  FortranEndLine();

  // -17C-
  FortranDebug("-17C-");
  MESSAGE("Writing SCNB_SCALE_FACTOR torsion\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG SCNB_SCALE_FACTOR");
  COMMENT_DIM(NPTRA)
  FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
  FortranFormat(5, DBLFORMAT_C);
  for (i = 0; i < iParmSetTotalTorsionParms(uUnit->psParameters); i++) {
    ParmSetTorsion(uUnit->psParameters, i, sAtom1, sAtom2, sAtom3, sAtom4,
                   &iN, &dKp, &dP0, &dScEE, &dScNB, sDesc);
    if (dScNB < 0.0) {
      dScNB = GDefaults.dScnbScaleFactor;
    }
    FortranWriteDouble(dScNB);
  }
  for (i = 0; i < iParmSetTotalImproperParms(uUnit->psParameters); i++) {
    FortranWriteDouble(0.0);
  }
  FortranEndLine();

  // -18- Not used, reserved for future use, uses NATYP
  // Corresponds to the AMBER SOLTY array
  FortranDebug("-18-");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG SOLTY");
  COMMENT_DIM(NATYP)
  FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
  FortranFormat(5, DBLFORMAT_C);
  for (i = 0; i < iParmSetTotalAtomParms(uUnit->psParameters); i++) {
    FortranWriteDouble(0.0);
  }
  FortranEndLine();

  // -19- Lennard jones r**12 term for all possible interactions
  // CN1 array
  FortranDebug("-19-");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG LENNARD_JONES_ACOEF");
  COMMENT("NOTE:NTTYP=NTYPES*(NTYPES+1)/2;")
  COMMENT_DIM_UNITS(NTTYP, "kcal/mol*angstrom^12")
  FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
  FortranFormat(5, DBLFORMAT_C);
  for (i = 0; i < iVarArrayElementCount(vaNBParameters); i++) {
    FortranWriteDouble(PVAI(vaNBParameters, NONBONDACt, i)->dA);
  }
  FortranEndLine();

  // -20- Lennard jones r**6 term for all possible interactions
  // CN2 array
  FortranDebug("-20-");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG LENNARD_JONES_BCOEF");
  COMMENT_DIM_UNITS(NTTYP, "kcal/mol*angstrom^6")
  FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
  FortranFormat(5, DBLFORMAT_C);
  for (i = 0; i < iVarArrayElementCount(vaNBParameters); i++) {
    FortranWriteDouble(PVAI(vaNBParameters, NONBONDACt, i)->dB);
  }
  FortranEndLine();

  // -21- Write the bond interactions that include hydrogen
  // Write the two indices into the atom table, then the index
  // into the interaction table
  FortranDebug("-21-");
  MESSAGE("Writing the bond interactions with hydrogens\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG BONDS_INC_HYDROGEN");
  COMMENT("STRIDE:3; MEMBERS:ATOM_I,ATOM_J,PARM_INDEX;")
  COMMENT_DIM(NBONH)
  FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
  FortranFormat(10, INTFORMAT_C);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaBonds); i++) {
    sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, i);
    aA = PVAI(uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom1 - 1)->aAtom;
    aD = PVAI(uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom2 - 1)->aAtom;
    if (bPERT_BOND(bPert, aA, aD)) {
      continue;
    }
    if (iAtomElement(aA) == HYDROGEN || iAtomElement(aD) == HYDROGEN) {
      FortranWriteInt(AMBERINDEX(sbPBond->iAtom1));
      FortranWriteInt(AMBERINDEX(sbPBond->iAtom2));
      FortranWriteInt(sbPBond->iParmIndex);
    }
  }
  FortranEndLine();

  // -22- Write the bond interactions that dont include hydrogen
  // Write the two indices into the atom table, then the index
  // into the interaction table
  FortranDebug("-22-");
  MESSAGE("Writing the bond interactions without hydrogens\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG BONDS_WITHOUT_HYDROGEN");
  COMMENT("STRIDE:3; MEMBERS:ATOM_I,ATOM_J,PARM_INDEX;")
  COMMENT_DIM(NBONA)
  FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
  FortranFormat(10, INTFORMAT_C);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaBonds); i++) {
    sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, i);
    aA = PVAI(uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom1 - 1)->aAtom;
    aD = PVAI(uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom2 - 1)->aAtom;
    if (bPERT_BOND(bPert, aA, aD)) {
      continue;
    }
    if (!(iAtomElement(aA) == HYDROGEN || iAtomElement(aD) == HYDROGEN)) {
      FortranWriteInt(AMBERINDEX(sbPBond->iAtom1));
      FortranWriteInt(AMBERINDEX(sbPBond->iAtom2));
      FortranWriteInt(sbPBond->iParmIndex);
    }
  }

  // Write out the (bond without H) RESTRAINT interactions
  // The iParmIndex field is set in RESTRAINTLOOP
  if ((iMax = iVarArrayElementCount(uUnit->vaRestraints))) {
    srPRestraint = PVAI(uUnit->vaRestraints, SAVERESTRAINTt, 0);
    for (i = 0; i < iMax; i++, srPRestraint++) {
      if (srPRestraint->iType == RESTRAINTBOND) {
        FortranWriteInt(AMBERINDEX(srPRestraint->iAtom1));
        FortranWriteInt(AMBERINDEX(srPRestraint->iAtom2));
        FortranWriteInt(srPRestraint->iParmIndex);
      }
    }
  }
  FortranEndLine();

  // -23- Write the angle interactions that include hydrogen
  // Write the three indices into the atom table, then the index
  // into the interaction table
  FortranDebug("-23-");

  MESSAGE("Writing the angle interactions with hydrogens\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG ANGLES_INC_HYDROGEN");
  COMMENT("STRIDE:4; MEMBERS:ATOM_I,ATOM_J,ATOM_K,PARM_INDEX;")
  COMMENT_DIM(NTHETH)
  FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
  FortranFormat(10, INTFORMAT_C);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAngles); i++) {
    saPAngle = PVAI(uUnit->vaAngles, SAVEANGLEt, i);
    aA = PVAI(uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom1 - 1)->aAtom;
    aB = PVAI(uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom2 - 1)->aAtom;
    aD = PVAI(uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom3 - 1)->aAtom;
    if (bPERT_ANGLE(bPert, aA, aB, aD)) {
      continue;
    }
    if (iAtomElement(aA) == HYDROGEN || iAtomElement(aB) == HYDROGEN ||
        iAtomElement(aD) == HYDROGEN) {
      FortranWriteInt(AMBERINDEX(saPAngle->iAtom1));
      FortranWriteInt(AMBERINDEX(saPAngle->iAtom2));
      FortranWriteInt(AMBERINDEX(saPAngle->iAtom3));
      FortranWriteInt(saPAngle->iParmIndex);
    }
  }
  FortranEndLine();

  // -24- Write the angle interactions that dont include hydrogen
  // Write the three indices into the atom table, then the index
  // into the interaction table
  FortranDebug("-24-");
  MESSAGE("Writing the angle interactions without hydrogens\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG ANGLES_WITHOUT_HYDROGEN");
  COMMENT("STRIDE:4; MEMBERS:ATOM_I,ATOM_J,ATOM_K,PARM_INDEX;")
  COMMENT_DIM(NTHETA)
  FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
  FortranFormat(10, INTFORMAT_C);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAngles); i++) {
    saPAngle = PVAI(uUnit->vaAngles, SAVEANGLEt, i);
    aA = PVAI(uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom1 - 1)->aAtom;
    aB = PVAI(uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom2 - 1)->aAtom;
    aD = PVAI(uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom3 - 1)->aAtom;
    if (bPERT_ANGLE(bPert, aA, aB, aD)) {
      continue;
    }
    if (!(iAtomElement(aA) == HYDROGEN || iAtomElement(aB) == HYDROGEN ||
          iAtomElement(aD) == HYDROGEN)) {
      FortranWriteInt(AMBERINDEX(saPAngle->iAtom1));
      FortranWriteInt(AMBERINDEX(saPAngle->iAtom2));
      FortranWriteInt(AMBERINDEX(saPAngle->iAtom3));
      FortranWriteInt(saPAngle->iParmIndex);
    }
  }

  // Write out the RESTRAINT interactions
  // The iParmIndex field is set in RESTRAINTLOOP
  if ((iMax = iVarArrayElementCount(uUnit->vaRestraints))) {
    srPRestraint = PVAI(uUnit->vaRestraints, SAVERESTRAINTt, 0);
    for (i = 0; i < iMax; i++, srPRestraint++) {
      if (srPRestraint->iType == RESTRAINTANGLE) {
        FortranWriteInt(AMBERINDEX(srPRestraint->iAtom1));
        FortranWriteInt(AMBERINDEX(srPRestraint->iAtom2));
        FortranWriteInt(AMBERINDEX(srPRestraint->iAtom3));
        FortranWriteInt(srPRestraint->iParmIndex);
      }
    }
  }
  FortranEndLine();

  // -25- Write the torsion interactions that include hydrogen
  // Write the three indices into the atom table, then the index
  // into the interaction table
  FortranDebug("-25-");
  MESSAGE("Writing the torsion interactions with hydrogens\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG DIHEDRALS_INC_HYDROGEN");
  COMMENT("STRIDE:5; MEMBERS:ATOM_I,ATOM_J,ATOM_K,ATOM_L,PARM_INDEX;")
  COMMENT_DIM(NPHIH)
  FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
  FortranFormat(10, INTFORMAT_C);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaTorsions); i++) {
    stPTorsion = PVAI(uUnit->vaTorsions, SAVETORSIONt, i);
    aA = PVAI(uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom1 - 1)->aAtom;
    aB = PVAI(uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom2 - 1)->aAtom;
    aC = PVAI(uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom3 - 1)->aAtom;
    aD = PVAI(uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom4 - 1)->aAtom;
    if (bPERT_TORSION(bPert, aA, aB, aC, aD)) {
      continue;
    }
    if (iAtomElement(aA) == HYDROGEN || iAtomElement(aB) == HYDROGEN ||
        iAtomElement(aC) == HYDROGEN || iAtomElement(aD) == HYDROGEN) {
      if ((AMBERINDEX(stPTorsion->iAtom3) == 0) || (AMBERINDEX(stPTorsion->iAtom4) == 0)) {
        MESSAGE("Had to turn torsion around to avoid K,L == 0\n");
        MESSAGE("Outer atoms: %s --- %s\n", sContainerName(aA), sContainerName(aD));
        MESSAGE("Old order %d %d %d %d\n", stPTorsion->iAtom1, stPTorsion->iAtom2,
                 stPTorsion->iAtom3, stPTorsion->iAtom4);
        SWAP(stPTorsion->iAtom1, stPTorsion->iAtom4, iTemp);
        SWAP(stPTorsion->iAtom2, stPTorsion->iAtom3, iTemp);
        MESSAGE("New order %d %d %d %d\n", stPTorsion->iAtom1, stPTorsion->iAtom2,
                 stPTorsion->iAtom3, stPTorsion->iAtom4);
      }
      if (stPTorsion->bProper) {
        iProper = 1;
      }
      else {
        iProper = -1;
      }
      if (stPTorsion->bCalc14) {
        iCalc14 = 1;
      }
      else {
        iCalc14 = -1;
      }
      if (GDefaults.iCharmm && iProper == -1) {
        FortranWriteInt(AMBERINDEX(stPTorsion->iAtom3));
        FortranWriteInt(AMBERINDEX(stPTorsion->iAtom2));
        FortranWriteInt(AMBERINDEX(stPTorsion->iAtom1) * iCalc14);
        FortranWriteInt(AMBERINDEX(stPTorsion->iAtom4) * iProper);
      }
      else {
        FortranWriteInt(AMBERINDEX(stPTorsion->iAtom1));
        FortranWriteInt(AMBERINDEX(stPTorsion->iAtom2));
        FortranWriteInt(AMBERINDEX(stPTorsion->iAtom3) * iCalc14);
        FortranWriteInt(AMBERINDEX(stPTorsion->iAtom4) * iProper);
      }
      FortranWriteInt(stPTorsion->iParmIndex);
    }
  }
  FortranEndLine();

  // -26- Write the torsion interactions that dont include hydrogen
  // Write the three indices into the atom table, then the index
  // into the interaction table
  FortranDebug("-26-");
  MESSAGE("Writing the torsion interactions without hydrogens\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG DIHEDRALS_WITHOUT_HYDROGEN");
  COMMENT("STRIDE:5; MEMBERS:ATOM_I,ATOM_J,ATOM_K,ATOM_L,PARM_INDEX;")
  COMMENT_DIM(NPHIA)
  FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
  FortranFormat(10, INTFORMAT_C);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaTorsions); i++) {
    stPTorsion = PVAI(uUnit->vaTorsions, SAVETORSIONt, i);
    aA = PVAI(uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom1 - 1)->aAtom;
    aB = PVAI(uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom2 - 1)->aAtom;
    aC = PVAI(uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom3 - 1)->aAtom;
    aD = PVAI(uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom4 - 1)->aAtom;
    if (bPERT_TORSION(bPert, aA, aB, aC, aD)) {
      continue;
    }
    if (!(iAtomElement(aA) == HYDROGEN || iAtomElement(aB) == HYDROGEN ||
          iAtomElement(aC) == HYDROGEN || iAtomElement(aD) == HYDROGEN)) {
      if ((AMBERINDEX(stPTorsion->iAtom3) == 0) || (AMBERINDEX(stPTorsion->iAtom4) == 0)) {
        MESSAGE("Had to turn torsion to avoid K,L == 0\n");
        MESSAGE("Outer atoms: %s --- %s\n", sContainerName(aA), sContainerName(aD));
        SWAP(stPTorsion->iAtom1, stPTorsion->iAtom4, iTemp);
        SWAP(stPTorsion->iAtom2, stPTorsion->iAtom3, iTemp);
      }
      if (stPTorsion->bCalc14) {
        iCalc14 = 1;
      }
      else {
        iCalc14 = -1;
      }
      if (stPTorsion->bProper) {
        iProper = 1;
      }
      else {
        iProper = -1;
      }
      if (GDefaults.iCharmm && iProper == -1) {
        FortranWriteInt(AMBERINDEX(stPTorsion->iAtom3));
        FortranWriteInt(AMBERINDEX(stPTorsion->iAtom2));
        FortranWriteInt(AMBERINDEX(stPTorsion->iAtom1) * iCalc14);
        FortranWriteInt(AMBERINDEX(stPTorsion->iAtom4) * iProper);
      }
      else {
        FortranWriteInt(AMBERINDEX(stPTorsion->iAtom1));
        FortranWriteInt(AMBERINDEX(stPTorsion->iAtom2));
        FortranWriteInt(AMBERINDEX(stPTorsion->iAtom3) * iCalc14);
        FortranWriteInt(AMBERINDEX(stPTorsion->iAtom4) * iProper);
      }
      FortranWriteInt(stPTorsion->iParmIndex);
    }
  }

  // Write out the RESTRAINT interactions
  // The iParmIndex field is set in RESTRAINTLOOP
  if ((iMax = iVarArrayElementCount(uUnit->vaRestraints))) {
    srPRestraint = PVAI(uUnit->vaRestraints, SAVERESTRAINTt, 0);
    for (i = 0; i < iMax; i++, srPRestraint++) {
      if (srPRestraint->iType == RESTRAINTTORSION) {
        if ((AMBERINDEX(srPRestraint->iAtom3) == 0) ||
            (AMBERINDEX(srPRestraint->iAtom4) == 0)) {
          MESSAGE("Had to turn RESTRAINT torsion around to avoid\n");
          MESSAGE("K,L == 0\n");
          SWAP(srPRestraint->iAtom1, srPRestraint->iAtom4, iTemp);
          SWAP(srPRestraint->iAtom2, srPRestraint->iAtom3, iTemp);
        }
        FortranWriteInt(AMBERINDEX(srPRestraint->iAtom1));
        FortranWriteInt(AMBERINDEX(srPRestraint->iAtom2));
        FortranWriteInt(AMBERINDEX(srPRestraint->iAtom3));
        FortranWriteInt(AMBERINDEX(srPRestraint->iAtom4));
        FortranWriteInt(srPRestraint->iParmIndex);
      }
    }
  }
  FortranEndLine();

  // -27- Write the excluded atom list
  FortranDebug("-27-");
  MESSAGE("Writing the excluded atom list\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG EXCLUDED_ATOMS_LIST");
  COMMENT_DIM(NNB)
  FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
  FortranFormat(10, INTFORMAT_C);
  for (i = 0; i < iVarArrayElementCount(vaExcludedAtoms); i++) {
    FortranWriteInt(*PVAI(vaExcludedAtoms, int, i));
  }
  FortranEndLine();

  // -28- Write the R^12 term for the Hydrogen bond equation
  FortranDebug("-28-");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG HBOND_ACOEF");
  COMMENT_DIM(NPHB)
  FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
  FortranFormat(5, DBLFORMAT_C);
  for (i = 0; i < iParmSetTotalHBondParms(uUnit->psParameters); i++) {
    ParmSetHBond(uUnit->psParameters, i, sType1, sType2, &dB, &dD, sDesc);
    FortranWriteDouble(dB);
  }
  FortranEndLine();

  // -29- Write the R^10 term for the Hydrogen bond equation
  FortranDebug("-29-");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG HBOND_BCOEF");
  COMMENT_DIM(NPHB)
  FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
  FortranFormat(5, DBLFORMAT_C);
  for (i = 0; i < iParmSetTotalHBondParms(uUnit->psParameters); i++) {
    ParmSetHBond(uUnit->psParameters, i, sType1, sType2, &dB, &dD, sDesc);
    FortranWriteDouble(dD);
  }
  FortranEndLine();

  // -30- No longer used, but stored
  FortranDebug("-30-");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG HBCUT");
  COMMENT_DIM(NPHB)
  FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
  FortranFormat(5, DBLFORMAT_C);
  for (i = 0; i < iParmSetTotalHBondParms(uUnit->psParameters); i++) {
    FortranWriteDouble(0.0);
  }
  FortranEndLine();

  // -31- List of atomic symbols
  FortranDebug("-31-");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG AMBER_ATOM_TYPE");
  COMMENT_DIM(NATOM)
  FortranWriteString("%FORMAT(20" LBLFORMAT_F ")");
  FortranFormat(20, LBLFORMAT_C);
  //FortranFormat(10, INTFORMAT_C);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
    FortranWriteString(sAtomType(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom));
    //FortranWriteInt(iAtomCoordination(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom)); //VeryNew
  }
  FortranEndLine();

  // -32- List of tree symbols
  FortranDebug("-32-");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG TREE_CHAIN_CLASSIFICATION");
  COMMENT_DIM(NATOM)
  COMMENT("DESC:mainchain tree class: M,E,S,B,3,4,5,6,X,BLA=unknown;")
  FortranWriteString("%FORMAT(20" LBLFORMAT_F ")");
  FortranFormat(20, LBLFORMAT_C);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
    aAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom;
    if (dAtomTemp(aAtom) == (double) 'M') {
      FortranWriteString("M  ");
    }
    else if (dAtomTemp(aAtom) == (double) 'E') {
      FortranWriteString("E  ");
    }
    else if (dAtomTemp(aAtom) == (double) 'S') {
      FortranWriteString("S  ");
    }
    else if (dAtomTemp(aAtom) == (double) 'B') {
      FortranWriteString("B  ");
    }
    else if (dAtomTemp(aAtom) == (double) '3') {
      FortranWriteString("3  ");
    }
    else if (dAtomTemp(aAtom) == (double) '4') {
      FortranWriteString("4  ");
    }
    else if (dAtomTemp(aAtom) == (double) '5') {
      FortranWriteString("5  ");
    }
    else if (dAtomTemp(aAtom) == (double) '6') {
     FortranWriteString("6  ");
    }
    else if (dAtomTemp(aAtom) == (double) 'X') {
     FortranWriteString("X  ");
    }
    else {
     FortranWriteString("BLA"); // nornally unknown type 'x'
    }
  }
  FortranEndLine();

  // -33- Tree Joining information !!!!!!! Add support for this !!!!!
  FortranDebug("-33-");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG JOIN_ARRAY");
  COMMENT_DIM(NATOM)
  FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
  FortranFormat(10, INTFORMAT_C);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
    FortranWriteInt(0);
  }
  FortranEndLine();

  // -34- Who knows, something to do with rotating atoms
  FortranDebug("-34-");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG IROTAT");
  COMMENT_DIM(NATOM)
  FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
  FortranFormat(10, INTFORMAT_C);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
    FortranWriteInt(0);
  }
  FortranEndLine();

  // -35A- The last residue before "solvent"
  // Number of molecules
  // Index of first molecule that is solvent
  if (bUnitUseBox(uUnit)) {
    FortranDebug("-35A-");

    // Find the index of the first solvent RESIDUE
    for (i = 0; i < iVarArrayElementCount(uUnit->vaResidues); i++) {
      if (PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->sResidueType[0] == RESTYPESOLVENT) {
        break;
      }
    }
    int IPTRES = i; // FIXME: first bulk solvent res?

    // Find the molecules and return the number of ATOMs in each
    // molecule, along with the index of the first solvent molecule
    //
    UnitIOFindAndCountMolecules( uUnit );
    int NSPM = iVarArrayElementCount(uUnit->vaAtomsPerMolecule);
    // +1 for FORTRAN indexing conversion
    int NSPSOL = uUnit->iFirstSolvent + 1;

    if (GDefaults.dPrmtopFormat<2.0) {
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG SOLVENT_POINTERS");
    COMMENT("MEMBERS:IPTRES,NSPM,NSPSOL;")
    COMMENT_SIZE(3)
    FortranWriteString("%FORMAT(3" INTFORMAT_F ")");
    FortranFormat(3, INTFORMAT_C);
    FortranWriteInt(IPTRES);
    FortranWriteInt(NSPM);
    FortranWriteInt(NSPSOL);
    FortranEndLine();
    } else {
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG IPTRES");
    FortranWriteString("%SIZE 1");
    FortranWriteString("%FORMAT(" INTFORMAT_F ")");
    FortranFormat(1, INTFORMAT_C);
    FortranWriteInt(IPTRES);

    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG NSPSOL");
    FortranWriteString("%SIZE 1");
    FortranWriteString("%FORMAT(" INTFORMAT_F ")");
    FortranFormat(1, INTFORMAT_C);
    FortranWriteInt(NSPSOL);
    FortranEndLine();
    }

    // -35B- The number of ATOMs in the Ith MOLECULE
    FortranDebug("-35B-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG ATOMS_PER_MOLECULE");
    COMMENT_DIM(NSPM)
    FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
    FortranFormat(10, INTFORMAT_C);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtomsPerMolecule); i++) {
      FortranWriteInt(*PVAI(uUnit->vaAtomsPerMolecule, int, i));
    }
    FortranEndLine();

    // -35C- BETA, (BOX(I), I=1,3 )
    FortranDebug("-35C-");
    FortranFormat(1, "%-80s");
    UnitGetBox(uUnit, &dX, &dY, &dZ);
    if (GDefaults.dPrmtopFormat >= 2.0 ||
         (dUnitBeta(uUnit)!=uUnit->dAlpha || dUnitBeta(uUnit)!=uUnit->dGamma) ) {
      FortranWriteString("%FLAG CELL_DIMENSIONS");
      COMMENT("SIZE: 6; UNITS:angstrom,degrees; MEMBERS:A,B,C,APHA,BETA,GAMMA;")
      FortranWriteString("%FORMAT(6E16.8)");
      FortranFormat(4, DBLFORMAT_C);
      FortranWriteDouble(dX);
      FortranWriteDouble(dY);
      FortranWriteDouble(dZ);
      FortranWriteDouble(uUnit->dAlpha / DEGTORAD);
      FortranWriteDouble(uUnit->dBeta  / DEGTORAD);
      FortranWriteDouble(uUnit->dGamma / DEGTORAD);
    } else {
      FortranWriteString("%FLAG BOX_DIMENSIONS");
      FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
      COMMENT("SIZE:4; UNITS:angstrom,degrees; MEMBERS:BETA,A,B,C;")
      FortranFormat(4, DBLFORMAT_C);
      FortranWriteDouble(dUnitBeta(uUnit) / DEGTORAD);
      FortranWriteDouble(dX);
      FortranWriteDouble(dY);
      FortranWriteDouble(dZ);
    }
    FortranEndLine();
  }

  // -35D- NATCAP
  if (bUnitUseSolventCap(uUnit)) {
    FortranDebug("-35D-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG CAP_INFO");
    COMMENT_SIZE(1)
    FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
    FortranFormat(1, INTFORMAT_C);
    FortranWriteInt(uUnit->iCapTempInt);
    FortranEndLine();

    // -35E- CUTCAP, XCAP, YCAP, ZCAP
    FortranDebug("-35E-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG CAP_INFO2");
    COMMENT("UNITS:angstrom; MEMBERS:CUTCAP,XCAP,YCAP,ZCAP;")
    COMMENT_SIZE(4)
    FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
    FortranFormat(4, DBLFORMAT_C);
    UnitGetSolventCap(uUnit, &dX, &dY, &dZ, &dR);
    FortranWriteDouble(dR);
    FortranWriteDouble(dX);
    FortranWriteDouble(dY);
    FortranWriteDouble(dZ);
    FortranEndLine();
  }

  // Write out the GB radii
  MESSAGE("Writing the GB radii\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG RADIUS_SET");
  COMMENT_SIZE(1)
  FortranWriteString("%FORMAT(1a80)");
  if (GDefaults.iGBparm>=0 && GDefaults.iGBparm<=8 && PBRadii_optionDesc[GDefaults.iGBparm])
      sprintf(sComment,"%s (%s)", PBRadii_optionDesc[GDefaults.iGBparm], PBRadii_options[GDefaults.iGBparm]);
  else
      sprintf(sComment,"Unknown radius set %d (programming error)",GDefaults.iGBparm );
  FortranWriteString(sComment);

  FortranWriteString("%FLAG RADII");
  COMMENT_DIM_UNITS(NATOM, "angstrom")
  COMMENT("DESC:GB radius;")
  FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
  FortranFormat(5, DBLFORMAT_C);
  saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, 0);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++,saPAtom++) FortranWriteDouble(dGBRad[i]);
  FortranEndLine();
  FREE(dGBRad);

  // Write out the GB screening parameters
  MESSAGE("Writing the GB screening parameters\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG SCREEN");
  COMMENT_DIM(NATOM)
  COMMENT("DESC:GB screening factor;")
  FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
  FortranFormat(5, DBLFORMAT_C);
  saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, 0);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++,saPAtom++) FortranWriteDouble(dGBScreen[i]);
  FortranEndLine();
  FREE(dGBScreen);

  // Write IPOL near the end of prmtop
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG IPOL");
  COMMENT_SIZE(1)
  COMMENT("DESC:Flag for polarizable force field;")
  FortranWriteString("%FORMAT(1I8)");
  FortranFormat(1, INTFORMAT_C);
  FortranWriteInt(GDefaults.iIPOL);
  FortranEndLine();

  // Write the perturbation information
  if (bPert) {

    // -36A- Bonds that are to be perturbed
    // Totally perturbed bonds first,
    // boundary second
    FortranDebug("-36A-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG PERT_BOND_ATOMS");
    COMMENT("STRIDE:2; MEMBERS:ATOM_I,ATOM_J;")
    COMMENT_DIM(NBPER)
    FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
    FortranFormat(10, INTFORMAT_C);
    if (iVarArrayElementCount(uUnit->vaBonds)) {
      sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, 0);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaBonds); i++, sbPBond++) {
      if (((sbPBond->fFlags & PERTURBED) != 0) && ((sbPBond->fFlags & BOUNDARY) == 0)) {
        FortranWriteInt(AMBERINDEX(sbPBond->iAtom1));
        FortranWriteInt(AMBERINDEX(sbPBond->iAtom2));
      }
    }
    }
    if (iVarArrayElementCount(uUnit->vaBonds)) {
      sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, 0);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaBonds); i++, sbPBond++) {
      if (((sbPBond->fFlags & PERTURBED) != 0) && ((sbPBond->fFlags & BOUNDARY) != 0)) {
        FortranWriteInt(AMBERINDEX(sbPBond->iAtom1));
        FortranWriteInt(AMBERINDEX(sbPBond->iAtom2));
      }
    }
    }
    FortranEndLine();

    // -36B- Index into bond interaction arrays
    FortranDebug("-36B-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG PERT_BOND_PARAMS");
    COMMENT_DIM(NBPER)
    FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
    FortranFormat(10, INTFORMAT_C);

    // First, LAMBDA = 0
    if (iVarArrayElementCount(uUnit->vaBonds)) {
      sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, 0);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaBonds); i++, sbPBond++) {
      if (((sbPBond->fFlags & PERTURBED) != 0) && ((sbPBond->fFlags & BOUNDARY) == 0)) {
        FortranWriteInt(sbPBond->iParmIndex);
      }
    }
    }
    if (iVarArrayElementCount(uUnit->vaBonds)) {
      sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, 0);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaBonds); i++, sbPBond++) {
      if (((sbPBond->fFlags & PERTURBED) != 0) && ((sbPBond->fFlags & BOUNDARY) != 0)) {
        FortranWriteInt(sbPBond->iParmIndex);
      }
    }
    }

    // Then LAMBDA = 1
    if (iVarArrayElementCount(uUnit->vaBonds)) {
      sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, 0);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaBonds); i++, sbPBond++) {
      if (((sbPBond->fFlags & PERTURBED) != 0) && ((sbPBond->fFlags & BOUNDARY) == 0)) {
        FortranWriteInt(sbPBond->iPertParmIndex);
      }
    }
    }
    if (iVarArrayElementCount(uUnit->vaBonds)) {
      sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, 0);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaBonds); i++, sbPBond++) {
      if (((sbPBond->fFlags & PERTURBED) != 0) && ((sbPBond->fFlags & BOUNDARY) != 0)) {
        FortranWriteInt(sbPBond->iPertParmIndex);
      }
    }
    }
    FortranEndLine();

    // -36C- Angles that are to be perturbed
    FortranDebug("-36C-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG PERT_ANGLE_ATOMS");
    COMMENT("STRIDE:3; MEMBERS:ATOM_I,ATOM_J,ATOM_K;")
    COMMENT_DIM(NGPER)
    FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
    FortranFormat(10, INTFORMAT_C);
    if (iVarArrayElementCount(uUnit->vaAngles)) {
      saPAngle = PVAI(uUnit->vaAngles, SAVEANGLEt, 0);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAngles); i++, saPAngle++) {
      if (((saPAngle->fFlags & PERTURBED) != 0) && ((saPAngle->fFlags & BOUNDARY) == 0)) {
        FortranWriteInt(AMBERINDEX(saPAngle->iAtom1));
        FortranWriteInt(AMBERINDEX(saPAngle->iAtom2));
        FortranWriteInt(AMBERINDEX(saPAngle->iAtom3));
      }
    }
    }
    if (iVarArrayElementCount(uUnit->vaAngles)) {
      saPAngle = PVAI(uUnit->vaAngles, SAVEANGLEt, 0);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAngles); i++, saPAngle++) {
      if (((saPAngle->fFlags & PERTURBED) != 0) && ((saPAngle->fFlags & BOUNDARY) != 0)) {
        FortranWriteInt(AMBERINDEX(saPAngle->iAtom1));
        FortranWriteInt(AMBERINDEX(saPAngle->iAtom2));
        FortranWriteInt(AMBERINDEX(saPAngle->iAtom3));
      }
    }
    }
    FortranEndLine();

    // -36D- Index into angle interaction arrays
    FortranDebug("-36D-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG PERT_ANGLE_PARAMS");
    COMMENT_DIM(NGPER)
    FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
    FortranFormat(10, INTFORMAT_C);

    // First LAMBDA = 0
    if (iVarArrayElementCount(uUnit->vaAngles)) {
      saPAngle = PVAI(uUnit->vaAngles, SAVEANGLEt, 0);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAngles); i++, saPAngle++) {
      if (((saPAngle->fFlags & PERTURBED) != 0) && ((saPAngle->fFlags & BOUNDARY) == 0)) {
        FortranWriteInt(saPAngle->iParmIndex);
      }
    }
    }
    if (iVarArrayElementCount(uUnit->vaAngles)) {
      saPAngle = PVAI(uUnit->vaAngles, SAVEANGLEt, 0);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAngles); i++, saPAngle++) {
      if (((saPAngle->fFlags & PERTURBED) != 0) && ((saPAngle->fFlags & BOUNDARY) != 0)) {
        FortranWriteInt(saPAngle->iParmIndex);
      }
    }
    }

    // Then LAMBDA = 1
    if (iVarArrayElementCount(uUnit->vaAngles)) {
      saPAngle = PVAI(uUnit->vaAngles, SAVEANGLEt, 0);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAngles); i++, saPAngle++) {
      if (((saPAngle->fFlags & PERTURBED) != 0) && ((saPAngle->fFlags & BOUNDARY) == 0)) {
        FortranWriteInt(saPAngle->iPertParmIndex);
      }
    }
    }
    if (iVarArrayElementCount(uUnit->vaAngles)) {
      saPAngle = PVAI(uUnit->vaAngles, SAVEANGLEt, 0);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAngles); i++, saPAngle++) {
      if (((saPAngle->fFlags & PERTURBED) != 0) && ((saPAngle->fFlags & BOUNDARY) != 0)) {
        FortranWriteInt(saPAngle->iPertParmIndex);
      }
    }
    }
    FortranEndLine();

    // -36E- Torsions that are to be perturbed
    FortranDebug("-36E-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG PERT_DIHEDRAL_ATOMS");
    COMMENT("STRIDE:4; MEMBERS:ATOM_I,ATOM_J,ATOM_K,ATOM_L;")
    COMMENT_DIM(NDPER)
    FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
    FortranFormat(10, INTFORMAT_C);
    if (iVarArrayElementCount(uUnit->vaTorsions)) {
      stPTorsion = PVAI(uUnit->vaTorsions, SAVETORSIONt, 0);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaTorsions); i++, stPTorsion++) {
      if (((stPTorsion->fFlags & PERTURBED) != 0) && ((stPTorsion->fFlags & BOUNDARY) == 0)) {
        if ((AMBERINDEX(stPTorsion->iAtom3) == 0) || (AMBERINDEX(stPTorsion->iAtom4) == 0)) {
          MESSAGE("Had to turn torsion around to avoid K,L == 0\n");
          MESSAGE("Outer atoms: %s --- %s\n", sContainerName(aA), sContainerName(aD));
          SWAP(stPTorsion->iAtom1, stPTorsion->iAtom4, iTemp);
          SWAP(stPTorsion->iAtom2, stPTorsion->iAtom3, iTemp);
        }
        if (stPTorsion->bProper) {
          iProper = 1;
        }
        else {
          iProper = -1;
        }
        if (stPTorsion->bCalc14) {
          iCalc14 = 1;
        }
        else {
          iCalc14 = -1;
	}
        FortranWriteInt(AMBERINDEX(stPTorsion->iAtom1));
        FortranWriteInt(AMBERINDEX(stPTorsion->iAtom2));
        FortranWriteInt(AMBERINDEX(stPTorsion->iAtom3) * iCalc14);
        FortranWriteInt(AMBERINDEX(stPTorsion->iAtom4) * iProper);
      }
    }
    }
    if (iVarArrayElementCount(uUnit->vaTorsions)) {
      stPTorsion = PVAI(uUnit->vaTorsions, SAVETORSIONt, 0);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaTorsions); i++, stPTorsion++) {
      if (((stPTorsion->fFlags & PERTURBED) != 0) && ((stPTorsion->fFlags & BOUNDARY) != 0)) {
        if ((AMBERINDEX(stPTorsion->iAtom3) == 0) || (AMBERINDEX(stPTorsion->iAtom4) == 0)) {
          MESSAGE("Had to turn torsion around to avoid K,L == 0\n");
          MESSAGE("Outer atoms: %s --- %s\n", sContainerName(aA), sContainerName(aD));
          SWAP(stPTorsion->iAtom1, stPTorsion->iAtom4, iTemp);
          SWAP(stPTorsion->iAtom2, stPTorsion->iAtom3, iTemp);
        }
        if (stPTorsion->bProper) {
          iProper = 1;
        }
        else {
          iProper = -1;
        }
        if (stPTorsion->bCalc14) {
          iCalc14 = 1;
        }
        else {
          iCalc14 = -1;
	}
        FortranWriteInt(AMBERINDEX(stPTorsion->iAtom1));
        FortranWriteInt(AMBERINDEX(stPTorsion->iAtom2));
        FortranWriteInt(AMBERINDEX(stPTorsion->iAtom3) * iCalc14);
        FortranWriteInt(AMBERINDEX(stPTorsion->iAtom4) * iProper);
      }
    }
    }
    FortranEndLine();

    // -36F- Index into torsion interaction arrays
    FortranDebug("-36F-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG PERT_DIHEDRAL_PARAMS");
    COMMENT_DIM(NDPER)
    FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
    FortranFormat(10, INTFORMAT_C);

    // First LAMBDA = 0
    if (iVarArrayElementCount(uUnit->vaTorsions)) {
      stPTorsion = PVAI(uUnit->vaTorsions, SAVETORSIONt, 0);
      for (i = 0; i < iVarArrayElementCount(uUnit->vaTorsions); i++, stPTorsion++) {
        if (((stPTorsion->fFlags & PERTURBED) != 0) && ((stPTorsion->fFlags & BOUNDARY) == 0)) {
          FortranWriteInt(stPTorsion->iParmIndex);
        }
      }
    }

    if (iVarArrayElementCount(uUnit->vaTorsions)) {
      stPTorsion = PVAI(uUnit->vaTorsions, SAVETORSIONt, 0);
      for (i = 0; i < iVarArrayElementCount(uUnit->vaTorsions); i++, stPTorsion++) {
        if (((stPTorsion->fFlags & PERTURBED) != 0) && ((stPTorsion->fFlags & BOUNDARY) != 0)) {
          FortranWriteInt(stPTorsion->iParmIndex);
        }
      }
    }

    // Then LAMBDA = 1
    if (iVarArrayElementCount(uUnit->vaTorsions)) {
      stPTorsion = PVAI(uUnit->vaTorsions, SAVETORSIONt, 0);
      for (i = 0; i < iVarArrayElementCount(uUnit->vaTorsions); i++, stPTorsion++) {
        if (((stPTorsion->fFlags & PERTURBED) != 0) && ((stPTorsion->fFlags & BOUNDARY) == 0)) {
          FortranWriteInt(stPTorsion->iPertParmIndex);
        }
      }
    }

    if (iVarArrayElementCount(uUnit->vaTorsions)) {
      stPTorsion = PVAI(uUnit->vaTorsions, SAVETORSIONt, 0);
      for (i = 0; i < iVarArrayElementCount(uUnit->vaTorsions); i++, stPTorsion++) {
        if (((stPTorsion->fFlags & PERTURBED) != 0) && ((stPTorsion->fFlags & BOUNDARY) != 0)) {
          FortranWriteInt(stPTorsion->iPertParmIndex);
        }
      }
    }
    FortranEndLine();

    // -36G- Residue labels at LAMBDA = 1
    // Just write the labels at LAMBDA = 0

    // Trim the string down to at most 3 characters by
    // taking the last three characters if it is too long
    FortranDebug("-36G-");
    MESSAGE("Writing the residue labels\n");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG PERT_RESIDUE_NAME");
    COMMENT_DIM(NRES)
    FortranWriteString("%FORMAT(20" LBLFORMAT_F ")");
    FortranFormat(20, LBLFORMAT_C);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaResidues); i++) {
      cPTemp = PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->sName;
      if (strlen(cPTemp) > 3) {
        cPTemp += (strlen(cPTemp) - 3);
      }
      FortranWriteString(cPTemp);
    }
    FortranEndLine();

    // -36H- Atom names at LAMBDA = 0
    FortranDebug("-36H-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG PERT_ATOM_NAME");
    COMMENT_DIM(NATOM)
    FortranWriteString("%FORMAT(20" LBLFORMAT_F ")");
    FortranFormat(20, LBLFORMAT_C);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
      cPTemp = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->sPertName;
      if (strlen(cPTemp) == 0) cPTemp = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->sName;
      if (strlen(cPTemp) > 4) cPTemp += (strlen(cPTemp) - 4);
      FortranWriteString(cPTemp);
    }
    FortranEndLine();

    // -36I- List of atomic symbols (atom types??????)
    FortranDebug("-36I-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG PERT_ATOM_SYMBOL");
    COMMENT_DIM(NATOM)
    FortranWriteString("%FORMAT(20" LBLFORMAT_F ")");
    FortranFormat(20, LBLFORMAT_C);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
      cPTemp = sAtomPertType(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom);
      if (strlen(cPTemp) == 0) cPTemp = sAtomType(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom);
      if (strlen(cPTemp) > 3) cPTemp += (strlen(cPTemp) - 3);
      FortranWriteString(cPTemp);
    }
    FortranEndLine();

    // -36J- Value of LAMBDA for each ATOM ?????????
    // TODO: Figure out what the hell this is
    FortranDebug("-36J-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG ALMPER");
    COMMENT_DIM(NATOM)
    FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
    FortranFormat(5, DBLFORMAT_C);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
      FortranWriteDouble(0.0);
    }
    FortranEndLine();

    // -36K- Flag to tell whether the atom is perturbed
    FortranDebug("-36K-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG IAPER");
    COMMENT_DIM(NATOM)
    FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
    FortranFormat(10, INTFORMAT_C);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
      FortranWriteInt(bAtomPerturbed(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom)?1:0);
    }
    FortranEndLine();

    // -36L- List of atom types - IACPER
    FortranDebug("-36L-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG PERT_ATOM_TYPE_INDEX");
    COMMENT_DIM(NATOM)
    FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
    FortranFormat(10, INTFORMAT_C);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
      iAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->iPertTypeIndex - 1;
      iTemp = *PVAI(vaNBIndex, int, iAtom);
      FortranWriteInt(iTemp + 1);
    }
    FortranEndLine();

    // -36M- Perturbed charges
    FortranDebug("-36M-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG PERT_CHARGE");
    COMMENT_DIM(NATOM)
    FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
    FortranFormat(5, DBLFORMAT_C);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
      aAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom;
      if (GDefaults.iGibbs) {
        if (bAtomPerturbed(aAtom)) {
          FortranWriteDouble(ELECTRONTOKCAL *
                             dAtomPertCharge(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom));
        }
        else {
          FortranWriteDouble(ELECTRONTOKCAL *
                             dAtomCharge(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom));
	}
      }
      else {
        FortranWriteDouble(ELECTRONTOKCAL *
                           (dAtomCharge(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom) +
                            dAtomPertCharge(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom)));
      }
    }
    FortranEndLine();
  }

  // Polarizabilities
  if (bPolar) {
    iCount = 0;
    iCountPerturbed = 0;
    iMax = iVarArrayElementCount(uUnit->vaAtoms);
    MESSAGE("Writing the atomic polarizabilities\n");

    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG POLARIZABILITY");
    COMMENT_DIM_UNITS(NATOM, "angstrom^3")
    FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
    FortranFormat(5, DBLFORMAT_C);
    saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, 0);
    for (i = 0; i < iMax; i++, saPAtom++) {
      iIndex = iParmSetFindAtom(uUnit->psParameters, saPAtom->sType);
      ParmSetAtom(uUnit->psParameters, iIndex, sType, &dMass, &dPolar, &dEpsilon, &dRStar,
                  &dEpsilon14, &dRStar14, &dScreenF, &iElement, &iHybridization, sDesc);
      if (dPolar == -1.0) {
        dPolar = 0.0;
        iCount++;
      }
      FortranWriteDouble(dPolar);
    }
    if (iCount > 0) {
      VP0("Total atoms with default polarization=0.0: %d of %d\n", iCount, iMax);
    }
    FortranEndLine();

    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG DIPOLE_DAMP_FACTOR");
    COMMENT_DIM(NATOM)
    FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
    FortranFormat(5, DBLFORMAT_C);
    saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, 0);
    for (i = 0; i < iMax; i++, saPAtom++) {
      iIndex = iParmSetFindAtom(uUnit->psParameters, saPAtom->sType);
      ParmSetAtom(uUnit->psParameters, iIndex, sType, &dMass, &dPolar, &dEpsilon, &dRStar,
                  &dEpsilon14, &dRStar14, &dScreenF, &iElement, &iHybridization, sDesc);
      if (dScreenF == 0.0) {
        dScreenF = GDefaults.dDipoleDampFactor;
      }
      FortranWriteDouble(dScreenF);
    }
    FortranEndLine();

    if (bPert) {
      int iPertTot = 0;
      FortranFormat(1, "%-80s");
      FortranWriteString("%FLAG PERT_POLARIZABILITY");
      COMMENT_DIM_UNITS(NATOM, "angstrom^3")
      FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
      FortranFormat(5, DBLFORMAT_C);
      saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, 0);
      for (i = 0; i < iMax; i++, saPAtom++) {
        bool bTmp = bAtomPerturbed(saPAtom->aAtom);
        if (bTmp) {
          iIndex = iParmSetFindAtom(uUnit->psParameters, saPAtom->sPertType);
          iPertTot++;
        }
        else {
          iIndex = iParmSetFindAtom(uUnit->psParameters, saPAtom->sType);
        }
        ParmSetAtom(uUnit->psParameters, iIndex, sType, &dMass, &dPolar, &dEpsilon, &dRStar,
                    &dEpsilon14, &dRStar14, &dScreenF, &iElement, &iHybridization, sDesc);
        if (dPolar == -1.0) {
          dPolar = 0.0;
          if (bTmp) {
            iCountPerturbed++;
	  }
        }
        FortranWriteDouble(dPolar);
      }
      FortranEndLine();
      if (iCountPerturbed > 0) {
        VP0("Total pert atoms with default polarization=0.0: %d of %d\n", iCountPerturbed,
             iPertTot);
      }
    }
  }

  //  Charmm-style parameters
  if (GDefaults.iCharmm) {

    // -19- Lennard jones r**12 term for all 14 interactions
    // CN114 array
    FortranDebug("-19-");
    FortranWriteString("%FLAG LENNARD_JONES_14_ACOEF");
    COMMENT_DIM_UNITS(NTTYP, "kcal/mol*angstrom^12")
    FortranFormat(5, DBLFORMAT_C);
    for (i = 0; i < iVarArrayElementCount(vaNBParameters); i++) {
      FortranWriteDouble(PVAI(vaNBParameters, NONBONDACt, i)->dA14);
    }
    FortranEndLine();

    // -20- Lennard jones r**6 term for all 14 interactions
    // CN214 array
    FortranDebug("-20-");
    FortranWriteString("%FLAG LENNARD_JONES_14_BCOEF");
    COMMENT_DIM_UNITS(NTTYP, "kcal/mol*angstrom^6")
    FortranFormat(5, DBLFORMAT_C);
    for (i = 0; i < iVarArrayElementCount(vaNBParameters); i++) {
      FortranWriteDouble(PVAI(vaNBParameters, NONBONDACt, i)->dB14);
    }
    FortranEndLine();

    // -13- Force constants for Urey-Bradley
    FortranDebug("-13-");
    FortranWriteString("%FLAG CHARMM_UREY_BRADLEY_FORCE_CONSTANT");
    COMMENT_DIM(NUMANG)
    FortranFormat(5, DBLFORMAT_C);
    for (i = 0; i < iParmSetTotalAngleParms(uUnit->psParameters); i++) {
      ParmSetAngle(uUnit->psParameters, i, sAtom1, sAtom2, sAtom3, &dKt,
                   &dT0, &dTkub, &dRkub, sDesc);
      FortranWriteDouble(dTkub);
    }
    FortranEndLine();

    // -14- Equilibrium distances for Urey-Bradley
    FortranDebug("-14-");
    FortranWriteString("%FLAG CHARMM_UREY_BRADLEY_EQUIL_VALUE");
    COMMENT_DIM(NUMANG)
    FortranFormat(5, DBLFORMAT_C);
    for (i = 0; i < iParmSetTotalAngleParms(uUnit->psParameters); i++) {
      ParmSetAngle(uUnit->psParameters, i, sAtom1, sAtom2, sAtom3, &dKt,
                   &dT0, &dTkub, &dRkub, sDesc);
      FortranWriteDouble(dRkub);
    }
    FortranEndLine();
  }

  if ( GDefaults.bPdbKeepChainId) {
    MESSAGE("Writing residue PDB ResId\n");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG RESIDUE_NUMBER");
    COMMENT_DIM(NRES)
    FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
    FortranFormat(10, INTFORMAT_C);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaResidues); i++) {
      FortranWriteInt( PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->iPdbResSeq );
    }
    FortranEndLine();

    MESSAGE("Writing residue PDB ChainId\n");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG RESIDUE_CHAINID");
    COMMENT_DIM(NRES)
    FortranWriteString("%FORMAT(20" LBLFORMAT_F ")");
    FortranFormat(20, LBLFORMAT_C);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaResidues); i++) {
      FortranWriteString(PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->sChainId);
    }
    FortranEndLine();
  }

  // CMAP parameters, Mengjuei Hsieh and Yong Duan
  SaveAmberParmCMAP(uUnit, fOut);

  fclose(fOut);

  // Clean up arrays and return
  VarArrayDestroy(&vaNBIndexMatrix);
  VarArrayDestroy(&vaNBParameters);
  VarArrayDestroy(&vaExcludedAtoms);
  VarArrayDestroy(&vaExcludedCount);
  VarArrayDestroy(&vaNBIndex);
  VarArrayDestroy(&vaNonBonds);

}

void UnitIOSaveAmberCoord(UNIT uUnit, char *crdName) {
  FILE *fCrd;
  VECTOR vPos;
  int i;

  // Open the coordinate file
  fCrd = FOPENCOMPLAIN(crdName, "w");
  if (fCrd == NULL) {
    VPFATAL("%s: Could not open file: %s\n", crdName);
  }

  FortranFile(fCrd);
  FortranFormat(1, "%s");
  FortranWriteString(sContainerName(uUnit));
  FortranEndLine();
  FortranFormat(1, "%6d");
  FortranWriteInt(iVarArrayElementCount(uUnit->vaAtoms));
  FortranEndLine();
  FortranFormat(6, "%12.7lf");
  if (bUnitUseBox(uUnit)) {
    double dX, dY, dZ;
    double dX2, dY2, dZ2;
    UnitGetBox(uUnit, &dX, &dY, &dZ);
    if (GDefaults.nocenter == 0) {
      dX2 = dX * 0.5;
      dY2 = dY * 0.5;
      dZ2 = dZ * 0.5;
    } else {
      dX2 = 0.0;
      dY2 = 0.0;
      dZ2 = 0.0;
    }

    // Shift box to Amber spot; later, add a cmd opt or environment
    // var to switch between 0,0,0 center (spasms) or corner
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
      vPos = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->vPos;
      FortranWriteDouble(dVX(&vPos) + dX2);
      FortranWriteDouble(dVY(&vPos) + dY2);
      FortranWriteDouble(dVZ(&vPos) + dZ2);
    }
    FortranEndLine();
    FortranWriteDouble(dX);
    FortranWriteDouble(dY);
    FortranWriteDouble(dZ);
    FortranWriteDouble(dUnitBeta(uUnit) / DEGTORAD);
    FortranWriteDouble(dUnitBeta(uUnit) / DEGTORAD);
    FortranWriteDouble(dUnitBeta(uUnit) / DEGTORAD);
    FortranEndLine();
  }
  else {
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
      vPos = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->vPos;
      FortranWriteDouble(dVX(&vPos));
      FortranWriteDouble(dVY(&vPos));
      FortranWriteDouble(dVZ(&vPos));
    }
    FortranEndLine();
  }
  fclose(fCrd);
}


/*
 *      UnitIOSaveAmberParmFormat_old
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Save the UNIT in the AMBER PARM file format.
 *      This requires that the UNIT tables be built and that
 *      the UNIT contain a parameter set.
 *      The iContainerTempInt(atom) should still return the index
 *      of the atom within the vaAtoms array.
 *      Atom coordinates are written to the file (fCrd).
 *
 *NOTE: This routine depends on the order of the RESIDUEs in
 *      vaResidues being such that solvent residues follow
 *      all other RESIDUEs.  I know that this is going
 *      against the philosophy that the data written to
 *      OFF files has NO implicit order, and outside of
 *      this program that is how they should be handled.
 *      But it was SO convenient to sort the RESIDUEs
 *      as they are put into the table that I could
 *      not resist.
 */
#define INTFORMAT       "%6d"
#define DBLFORMAT       "%16.8lE"
#define LBLFORMAT       "%-4s"
#undef ELECTRONTOKCAL
#define ELECTRONTOKCAL  18.2223

void
UnitIOSaveAmberParmFormat_old( UNIT uUnit, FILE *fOut, char *crdName, 
        bool bPolar, bool bPert)
{
int             i, iMax, iIndex;
LOOP            lTemp, lSpan;
ATOM            aAtom, aA, aB, aC, aD;
int             iCount, iBondWith, iBondWithout;
int             iAngleWith, iAngleWithout;
int             iTorsionWith, iTorsionWithout;
VARARRAY        vaExcludedAtoms, vaExcludedCount;
VARARRAY        vaNBIndexMatrix, vaNBParameters;
VARARRAY        vaNBIndex, vaNonBonds;
int             iCountPerturbed, iCountBondPerturbed, iCountBondBoundary;
int             iCountAnglePerturbed, iCountAngleBoundary;
int             iCountTorsionPerturbed, iCountTorsionBoundary;
int             iNumExtra;
SAVEBONDt       *sbPBond;
SAVEANGLEt      *saPAngle;
SAVEATOMt       *saPAtom;
SAVETORSIONt    *stPTorsion;
SAVERESTRAINTt  *srPRestraint;
double          dMass, dPolar, dR, dKb, dR0, dKt, dT0, dTkub, dRkub, dKp, dP0, dC, dD, dTemp;
 double		dScEE, dScNB, dScreenF, dKpull, dRpull0, dKpress, dRpress0;
STRING          sAtom1, sAtom2, sAtom3, sAtom4, sType1, sType2;
int             iN, iAtoms, iMaxAtoms, iTemp, iAtom, iCalc14, iProper;
int             iElement, iHybridization, iStart;
RESIDUE         rRes;
bool            bFoundSome;
VECTOR          vPos;
char            *cPTemp;
double          dX, dY, dZ, dEpsilon, dRStar, dEpsilon14, dRStar14;
STRING          sDesc, sType;
IX_REC          *ePResEnt;
IX_DESC         iResIx;

  // Open the coordinate file
  FILE *fCrd = FOPENCOMPLAIN( crdName, "w" );
  if ( fOut == NULL ) {
    VP0( "Could not open file: %s\n", crdName );
  }

                /* Build the excluded atom list */

    MESSAGE( "Building the excluded atom list\n" );
    vaExcludedCount = vaVarArrayCreate( sizeof(int) );
    vaExcludedAtoms = vaVarArrayCreate( sizeof(int) );

    iCountPerturbed = 0;
    for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
        aAtom = PVAI( uUnit->vaAtoms, SAVEATOMt, i )->aAtom;

        if ( bAtomFlagsSet( aAtom, ATOMPERTURB ) )
                iCountPerturbed++;

        lSpan = lLoop( (OBJEKT)aAtom, SPANNINGTREE );
        iCount = 0;
        bFoundSome = FALSE;
        iStart = iVarArrayElementCount( vaExcludedAtoms );
        while ( (aA = (ATOM)oNext(&lSpan)) ) {

            if ( aA == aAtom ) continue;
            
                /* If the atom is more than three away from the first atom */
                /* then it is not in the excluded atom list */

            if ( iAtomBackCount(aA) >= 4 ) break;

            if ( iContainerTempInt(aA) > iContainerTempInt(aAtom) ) { 
                VarArrayAdd( vaExcludedAtoms, (GENP)&iContainerTempInt(aA) );
                bFoundSome = TRUE;
                iCount++;
            }
        }
        if ( !bFoundSome ) {
            iAtoms = 0;
            VarArrayAdd( vaExcludedAtoms, (GENP)&iAtoms );
            iCount++;
        } else {

                /* Sort the part of the VARARRAY just added so that */
                /* the excluded ATOMs are in ascending order by index */

            SortByInteger( (GENP) PVAI( vaExcludedAtoms, int, iStart ),
                                iCount,
                                sizeof(int),
                               (GENP) PVAI( vaExcludedAtoms, int, iStart ),
                                TRUE );
        } 

        VarArrayAdd( vaExcludedCount, (GENP)&iCount );
    }

    /*
     *  mark main chain atoms where possible, noting the 
     *  number of atoms in the largest residue. keep
     *  track of residues which can't be marked.
     */
    VP0( "Not Marking per-residue atom chain types.\n" );
    iMaxAtoms = 0;
    create_index( &iResIx, 2, 0 );
    ePResEnt = (IX_REC *)MALLOC(sizeof(IX_REC)+8 );
    ePResEnt->recptr = NULL;    /* for Purify */

    VP0( "Marking per-residue atom chain types.\n" );

    iMaxAtoms = 0;
    lTemp = lLoop( (OBJEKT)uUnit, RESIDUES );
    while ( (rRes = (RESIDUE)oNext(&lTemp)) ) {
        iAtoms = iMarkMainChainAtoms( rRes, FALSE );
        if ( iAtoms > 0 )
                MarkSideChainAtoms( rRes );
        if ( iAtoms < 0 ) {
                iAtoms = -iAtoms;
                /*
                 *  couldn't mark main chains
                 */

                strcpy( ePResEnt->key, rRes->cHeader.sName );
                if ( add_key( ePResEnt, &iResIx ) != IX_OK )
                        DFATAL( "add_key() residue chain\n" );
        }
        if ( iAtoms > iMaxAtoms ) 
                iMaxAtoms = iAtoms;
    }
    /*
     *  print warnings
     */
    first_key( &iResIx );
    i = 1;
    while ( next_key( ePResEnt, &iResIx ) == IX_OK ) {
        if ( i ) {
                VP0( "  (Residues lacking connect0/connect1 - \n" );
                VP0( "   these don't have chain types marked:\n\n" );
                VP0( "\tres\ttotal affected\n\n" );
                i = 0;
        }
        VP0( "\t%s\t%d\n", ePResEnt->key, ePResEnt->count);
    }
    if (!i)
        VP0( "  )\n" );
    destroy_index( &iResIx );
    FREE( ePResEnt );

                /* Build the NON-BOND arrays that AMBER needs */

    zUnitIOBuildNonBondArrays( uUnit, &vaNBIndexMatrix, &vaNBParameters,
                                        &vaNBIndex, &vaNonBonds);
 

    FortranFile( fOut );

#if 0
        /*
         *---------------------------------------------------------
         *
         *      Turn on debugging of fortran format output file
         *      by sticking comments into the file.
         */

    FortranDebugOn();
#endif

    
        /* -1- Save the title of the UNIT */
    FortranDebug( "-1-" );
    MESSAGE( "Saving the name of the UNIT\n" );
    FortranFormat( 1, "%-80s" );
    FortranWriteString( sContainerName(uUnit) );

        /* -2- Save control information */
    FortranDebug( "-2-" );
    MESSAGE( "Saving all the main control variables\n" );
    FortranFormat( 12, INTFORMAT );

/*NTOTAT*/      
    FortranWriteInt( iVarArrayElementCount( uUnit->vaAtoms ) );
/*NTYPES*/      
    FortranWriteInt( iVarArrayElementCount( vaNonBonds ) );
    
        /* Count the number of bonds with hydrogens, and without */

    iBondWith = 0;
    iBondWithout = 0;
    for ( i=0; i<iVarArrayElementCount( uUnit->vaBonds ); i++ ) {
        sbPBond = PVAI( uUnit->vaBonds, SAVEBONDt, i );
        aA = PVAI( uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom1-1 )->aAtom;
        aD = PVAI( uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom2-1 )->aAtom;
        if ( bPERT_BOND(bPert,aA,aD) ) 
                continue;
        if ( iAtomElement(aA) == HYDROGEN ||
                iAtomElement(aD) == HYDROGEN ) iBondWith++;
        else iBondWithout++;
    }
/*NBONH*/       
    FortranWriteInt( iBondWith );
/*NBONA*/       
    FortranWriteInt( iBondWithout );

        /* Count the number of angles with hydrogens, and without */

    iAngleWith = 0;
    iAngleWithout = 0;
    for ( i=0; i<iVarArrayElementCount( uUnit->vaAngles ); i++ ) {
        saPAngle = PVAI( uUnit->vaAngles, SAVEANGLEt, i );
        aA = PVAI( uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom1-1 )->aAtom;
        aB = PVAI( uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom2-1 )->aAtom;
        aD = PVAI( uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom3-1 )->aAtom;
        if ( bPERT_ANGLE(bPert,aA,aB,aD) ) 
                continue;
        if ( iAtomElement(aA) == HYDROGEN 
                || iAtomElement(aB) == HYDROGEN 
                || iAtomElement(aD) == HYDROGEN )
            iAngleWith++;
        else
            iAngleWithout++;
    }
/*NTHETH*/      
    FortranWriteInt( iAngleWith );
/*NTHETA*/      
    FortranWriteInt( iAngleWithout );

        /* Count the number of torsions with hydrogens, and without */

    iTorsionWith = 0;
    iTorsionWithout = 0;
    for ( i=0; i<iVarArrayElementCount( uUnit->vaTorsions ); i++ ) {
        stPTorsion = PVAI( uUnit->vaTorsions, SAVETORSIONt, i );
        aA = PVAI( uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom1-1 )->aAtom;
        aB = PVAI( uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom2-1 )->aAtom;
        aC = PVAI( uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom3-1 )->aAtom;
        aD = PVAI( uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom4-1 )->aAtom;
        if ( bPERT_TORSION(bPert,aA,aB,aC,aD) ) 
            continue;
        if ( iAtomElement(aA) == HYDROGEN || 
                iAtomElement(aB) == HYDROGEN ||
                iAtomElement(aC) == HYDROGEN ||
                iAtomElement(aD) == HYDROGEN ) 
            iTorsionWith++;
        else 
            iTorsionWithout++;
    }
/*NPHIH*/       
    FortranWriteInt( iTorsionWith );
/*NPHIA*/       
    FortranWriteInt( iTorsionWithout );

/*JHPARM*/      
    FortranWriteInt( 0 );
/*JPARM*/       
    FortranWriteInt( 0 );

        /* Write the number of excluded atoms */

/*NEXT*/        
    FortranWriteInt( iVarArrayElementCount(vaExcludedAtoms) );
/*NTOTRS*/      
    FortranWriteInt( iVarArrayElementCount(uUnit->vaResidues) );

        /* Write the number of bonds/angles/torsions without hydrogens */
        /* PLUS the number of RESTRAINT bonds/angles/torsions */

/*MBONA*/       
    FortranWriteInt( iBondWithout+
                        iUnitRestraintTypeCount( uUnit, RESTRAINTBOND ) );
/*MTHETA*/      
    FortranWriteInt( iAngleWithout+
                        iUnitRestraintTypeCount( uUnit, RESTRAINTANGLE ) );
/*MPHIA*/       
    FortranWriteInt( iTorsionWithout+
                        iUnitRestraintTypeCount( uUnit, RESTRAINTTORSION ) );

        /* Write the number of unique bond types, angle types, torsion types */
        /* Add in the number of RESTRAINT bonds/angles/torsion because */
        /* they will have new parameters */

/*MUMBND*/      
    FortranWriteInt( iParmSetTotalBondParms(uUnit->psParameters)+
                        iUnitRestraintTypeCount( uUnit, RESTRAINTBOND ) );
/*MUMANG*/      
    FortranWriteInt( iParmSetTotalAngleParms(uUnit->psParameters)+
                        iUnitRestraintTypeCount( uUnit, RESTRAINTANGLE ) );
/*NPTRA*/       
    FortranWriteInt( iParmSetTotalTorsionParms(uUnit->psParameters) +
                     iParmSetTotalImproperParms(uUnit->psParameters) +
                     iUnitRestraintTypeCount( uUnit, RESTRAINTTORSION ) );

                /* TODO - have different arrays for different restraint types*/
    if ( iVarArrayElementCount( uUnit->vaRestraints ) )
        VP0( " Restraints:  Bond %d  Angle %d  Torsion %d\n",
                iUnitRestraintTypeCount( uUnit, RESTRAINTBOND ),
                iUnitRestraintTypeCount( uUnit, RESTRAINTANGLE ),
                iUnitRestraintTypeCount( uUnit, RESTRAINTTORSION ) );
    else
        VP0( " (no restraints)\n" );

        /* The next parameter corresponds to NATYP in AMBER */
        /* I don't know what it does, and Dave Spellmeyer says that */
        /* he only uses it to skip over the SOLTY array */

/*NATYP*/       
    FortranWriteInt( iParmSetTotalAtomParms(uUnit->psParameters) );
/*NHB*/         
    FortranWriteInt( iParmSetTotalHBondParms(uUnit->psParameters) );
/*IFPERT*/
    if ( bPert )
        FortranWriteInt( 1 );
    else
        FortranWriteInt( 0 );



        /* Count the number of bonds to be perturbed, and those across the */
        /* perturbation/non-perturbed boundary */

    iCountBondPerturbed = 0;
    iCountBondBoundary = 0;
    for ( i=0; i<iVarArrayElementCount( uUnit->vaBonds ); i++ ) {
        sbPBond = PVAI( uUnit->vaBonds, SAVEBONDt, i );
        if ( (sbPBond->fFlags&PERTURBED) != 0 ) {
            iCountBondPerturbed++;
            if ( (sbPBond->fFlags&BOUNDARY) != 0 ) {
                MESSAGE( "Boundary pert bond %d-%d\n", 
                                sbPBond->iAtom1, sbPBond->iAtom2 );
                iCountBondBoundary++;
            }
        }
    }

    MESSAGE( "Perturbed bonds: %d\n", iCountBondPerturbed );
    MESSAGE( "Perturbed boundary bonds: %d\n", iCountBondBoundary );

        /* Count the number of angles to be perturbed, and those on the */
        /* boundary */

    iCountAnglePerturbed = 0;
    iCountAngleBoundary = 0;
    for ( i=0; i<iVarArrayElementCount( uUnit->vaAngles ); i++ ) {
        saPAngle = PVAI( uUnit->vaAngles, SAVEANGLEt, i );
        if ( (saPAngle->fFlags&PERTURBED) != 0 ) iCountAnglePerturbed++;
        if ( (saPAngle->fFlags&BOUNDARY) != 0 ) iCountAngleBoundary++;
    }

        /* Count the number of torsions and impropers to be perturbed */
        /* and those on the boundary */

    iCountTorsionPerturbed = 0;
    iCountTorsionBoundary = 0;
    for ( i=0; i<iVarArrayElementCount( uUnit->vaTorsions ); i++ ) {
        stPTorsion = PVAI( uUnit->vaTorsions, SAVETORSIONt, i );
        if ( (stPTorsion->fFlags&PERTURBED) != 0 ) iCountTorsionPerturbed++;
        if ( (stPTorsion->fFlags&BOUNDARY) != 0 ) iCountTorsionBoundary++;
    }

    /*NBPER*/   
    FortranWriteInt( iCountBondPerturbed );
    /*NGPER*/   
    FortranWriteInt( iCountAnglePerturbed );
    /*NDPER*/   
    FortranWriteInt( iCountTorsionPerturbed );
    /*MBPER*/   
    FortranWriteInt( iCountBondPerturbed-iCountBondBoundary );
    /*MGPER*/   
    FortranWriteInt( iCountAnglePerturbed-iCountAngleBoundary );
    /*MDPER*/   
    FortranWriteInt( iCountTorsionPerturbed-iCountTorsionBoundary );

        /* Save flag for periodic boundary conditions */

    /*IFBOX*/
    if ( bUnitUseBox(uUnit) ) {
        if ( bUnitBoxOct(uUnit) )
                FortranWriteInt( 2 );
        else
                FortranWriteInt( 1 );
    } else
        FortranWriteInt( 0 );

        /* Save the number of atoms in the largest residue */

    /*NMXRS*/
    //printf("iMaxAtoms (2) %i\n",iMaxAtoms);
    FortranWriteInt( iMaxAtoms );

        /* Save flag for cap information */

    /*IFCAP*/
    if ( bUnitUseSolventCap(uUnit) )
        FortranWriteInt( 1 );
    else
        FortranWriteInt( 0 );

    /*NUMEXTRA*/
    iNumExtra = 0;
    for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
                cPTemp = sAtomType( PVAI(uUnit->vaAtoms, SAVEATOMt, i )->aAtom );
                if( !strncmp( cPTemp, "EP", 2 )) iNumExtra++;
    }
    FortranWriteInt( iNumExtra );

    FortranEndLine();


        /* -3-  write out the names of the atoms */
    FortranDebug( "-3-" );

    MESSAGE( "Writing the names of the atoms\n" );
    FortranFormat( 20, LBLFORMAT );
    for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
        cPTemp = PVAI( uUnit->vaAtoms, SAVEATOMt, i )->sName;
        if ( strlen( cPTemp ) > 4 ) cPTemp += ( strlen(cPTemp)-4 );
        FortranWriteString( cPTemp );
    }
    FortranEndLine();
    
        /* -4- write out the atomic charges */
    FortranDebug( "-4-" );
        
    MESSAGE( "Writing the atomic charges\n" );
    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
        FortranWriteDouble( PVAI( uUnit->vaAtoms, SAVEATOMt, i )->dCharge *
                                ELECTRONTOKCAL );
    }
    FortranEndLine();

        /* -5- write out the atomic masses */
    FortranDebug( "-5-" );

    MESSAGE( "Writing the atomic masses\n" );         
    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
        saPAtom = PVAI( uUnit->vaAtoms, SAVEATOMt, i );
        iIndex = iParmSetFindAtom( uUnit->psParameters, saPAtom->sType );
        ParmSetAtom( uUnit->psParameters, iIndex, sType,
                        &dMass, &dPolar, &dEpsilon, &dRStar, &dEpsilon14, &dRStar14, 
			&dScreenF,
            &iElement, &iHybridization, sDesc );
        FortranWriteDouble( dMass );
    }
    FortranEndLine();

        /* -6- write out the atomic types */
    FortranDebug( "-6-" );

    MESSAGE( "Writing the atomic types\n" );          
    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
        iAtom = PVAI( uUnit->vaAtoms, SAVEATOMt, i )->iTypeIndex-1;
        iTemp = *PVAI( vaNBIndex, int, iAtom );
        FortranWriteInt( iTemp+1 );
    }
    FortranEndLine();

        /* -7- write out the starting index into the excluded atom list */
    FortranDebug( "-7-" );

    MESSAGE( "Writing the starting index into the excluded atom list\n" );
    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
        FortranWriteInt( *PVAI( vaExcludedCount, int, i ) );
    }
    FortranEndLine();

        /* -8- Write the index for the position of the non bond type */
                /* of each type */
    FortranDebug( "-8-" );

    MESSAGE( "writing position of the non bond type of each type\n");
    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount(vaNBIndexMatrix); i++ ) {
        FortranWriteInt( *PVAI( vaNBIndexMatrix, int, i ) );
    }
    FortranEndLine();

        /* -9- Residue labels */
                /* Trim the string down to at most 3 characters by */
                /* taking the last three characters if it is too long */
    FortranDebug( "-9-" );

    MESSAGE( "Writing the residue labels\n" );
    FortranFormat( 20, LBLFORMAT );
    for ( i=0; i<iVarArrayElementCount(uUnit->vaResidues); i++ ) {
        cPTemp = PVAI( uUnit->vaResidues, SAVERESIDUEt, i )->sName;
        if ( strlen( cPTemp ) > 3 ) cPTemp += ( strlen(cPTemp)-3 );
        FortranWriteString( cPTemp );
    }
    FortranEndLine();

        /* -10- Pointer list for all the residues */
    FortranDebug( "-10-" );

    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount(uUnit->vaResidues); i++ ) {
        FortranWriteInt( PVAI( uUnit->vaResidues, 
                                SAVERESIDUEt, i )->iAtomStartIndex );
    }
    FortranEndLine();

        /* -11- Force constants for bonds */
    FortranDebug( "-11-" );

    MESSAGE( "Writing bond force constants\n" );
    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iParmSetTotalBondParms(uUnit->psParameters); i++ ) {
      ParmSetBond(uUnit->psParameters, i, sAtom1, sAtom2, &dKb, &dR0, &dKpull, &dRpull0,
		  &dKpress, &dRpress0, sDesc);
       FortranWriteDouble( dKb );
    }
                /* Write the RESTRAINT constants AND set the index */
                /* for where the interaction can find its constants */
    RESTRAINTLOOP( RESTRAINTBOND, dKx, i+1 );
    FortranEndLine();

        /* -12- Equilibrium bond lengths */
    FortranDebug( "-12-" );

    MESSAGE( "Writing equilibrium bond lengths\n" );
    FortranFormat( 5, DBLFORMAT );
    for (i = 0; i < iParmSetTotalBondParms(uUnit->psParameters); i++) {
      ParmSetBond(uUnit->psParameters, i, sAtom1, sAtom2, &dKb, &dR0, &dKpull, &dRpull0,
		  &dKpress, &dRpress0, sDesc);
      FortranWriteDouble( dR0 );
    }
                /* Write the bond RESTRAINT constants AND set the index */
                /* for where the interaction can find its constants */
    RESTRAINTLOOP( RESTRAINTBOND, dX0, i+1 );
    FortranEndLine();

        /* -13- Force constants for angles */
    FortranDebug( "-13-" );

    MESSAGE( "Writing angle force constants\n" );
    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iParmSetTotalAngleParms(uUnit->psParameters); i++ ) {
        ParmSetAngle( uUnit->psParameters, i, sAtom1, sAtom2, sAtom3,
                        &dKt, &dT0, &dTkub, &dRkub, sDesc );
        FortranWriteDouble( dKt );
    }
                /* Write the angle RESTRAINT constants AND set the index */
                /* for where the interaction can find its constants */
    RESTRAINTLOOP( RESTRAINTANGLE, dKx, i+1 );
    FortranEndLine();

        /* -14- Equilibrium angle values */
    FortranDebug( "-14-" );

    MESSAGE( "Writing equilibrium angle values\n" );
    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iParmSetTotalAngleParms(uUnit->psParameters); i++ ) {
        ParmSetAngle( uUnit->psParameters, i, sAtom1, sAtom2, sAtom3,
                        &dKt, &dT0, &dTkub, &dRkub, sDesc );
        FortranWriteDouble( dT0 );
    }
                /* Write the angle RESTRAINT constants AND set the index */
                /* for where the interaction can find its constants */
    RESTRAINTLOOP( RESTRAINTANGLE, dX0, i+1 );
    FortranEndLine();

        /* -15- Force constants for torsions */
    FortranDebug( "-15-" );

    MESSAGE( "Writing torsional force constants\n" );
    FortranFormat( 5, DBLFORMAT );
    MESSAGE( "There are %d torsions and %d impropers\n", 
                iParmSetTotalTorsionParms(uUnit->psParameters),
                iParmSetTotalImproperParms(uUnit->psParameters) );
    for ( i=0; i<iParmSetTotalTorsionParms(uUnit->psParameters); i++ ) {
        ParmSetTorsion( uUnit->psParameters, i, sAtom1, sAtom2, 
                        sAtom3, sAtom4,
                        &iN, &dKp, &dP0, &dScEE, &dScNB, sDesc );
        MESSAGE( "Torsion %d  %s-%s-%s-%s %d %lf %lf\n",
                        i, sAtom1, sAtom2, sAtom3, sAtom4,
                        iN, dKp, dP0 );
        FortranWriteDouble( dKp );
    }
    for ( i=0; i<iParmSetTotalImproperParms(uUnit->psParameters); i++ ) {
        ParmSetImproper( uUnit->psParameters, i, sAtom1, sAtom2, 
                        sAtom3, sAtom4,
                        &iN, &dKp, &dP0, sDesc );
        MESSAGE( "Improper %d  %s-%s-%s-%s %d %lf %lf\n",
                        i, sAtom1, sAtom2, sAtom3, sAtom4,
                        iN, dKp, dP0 );
        FortranWriteDouble( dKp );
    }
                /* Write the torsion RESTRAINT constants AND set the index */
                /* for where the interaction can find its constants */
    RESTRAINTLOOP( RESTRAINTTORSION, dKx, i+1 );
    FortranEndLine();

        /* -16- Division factor for the dihedral angles */
    FortranDebug( "-16-" );

    MESSAGE( "Writing multiplicity of torsion interaction\n" );
    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iParmSetTotalTorsionParms(uUnit->psParameters); i++ ) {
        ParmSetTorsion( uUnit->psParameters, i, sAtom1, sAtom2, 
                        sAtom3, sAtom4,
                        &iN, &dKp, &dP0, &dScEE, &dScNB, sDesc );
        dTemp = iN;
        FortranWriteDouble( dTemp );
    }
    for ( i=0; i<iParmSetTotalImproperParms(uUnit->psParameters); i++ ) {
        ParmSetImproper( uUnit->psParameters, i, sAtom1, sAtom2, 
                        sAtom3, sAtom4,
                        &iN, &dKp, &dP0, sDesc );
        dTemp = iN;
        FortranWriteDouble( dTemp );
    }
                /* Write the torsion RESTRAINT constants AND set the index */
                /* for where the interaction can find its constants */
    RESTRAINTLOOP( RESTRAINTTORSION, dX0, i+1 );
    FortranEndLine();

        /* -17- Phase for torsions */
    FortranDebug( "-17-" );

    MESSAGE( "Writing phase for torsion interactions\n" );
    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iParmSetTotalTorsionParms(uUnit->psParameters); i++ ) {
        ParmSetTorsion( uUnit->psParameters, i, sAtom1, sAtom2, 
                        sAtom3, sAtom4,
                        &iN, &dKp, &dP0, &dScEE, &dScNB, sDesc );
        FortranWriteDouble( dP0 );
    }
    for ( i=0; i<iParmSetTotalImproperParms(uUnit->psParameters); i++ ) {
        ParmSetImproper( uUnit->psParameters, i, sAtom1, sAtom2, 
                        sAtom3, sAtom4,
                        &iN, &dKp, &dP0, sDesc );
        FortranWriteDouble( dP0 );
    }
                /* Write the torsion RESTRAINT constants AND set the index */
                /* for where the interaction can find its constants */
    RESTRAINTLOOP( RESTRAINTTORSION, dN, i+1 );
    FortranEndLine();

        /* -18- Not used, reserved for future use, uses NATYP */
                /* Corresponds to the AMBER SOLTY array */
    FortranDebug( "-18-" );

    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iParmSetTotalAtomParms(uUnit->psParameters); i++ ) {
        FortranWriteDouble( 0.0 );
    }
    FortranEndLine();

        /* -19- Lennard jones r**12 term for all possible interactions */
                /* CN1 array */
    FortranDebug( "-19-" );

    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iVarArrayElementCount(vaNBParameters); i++ ) {
        FortranWriteDouble( PVAI( vaNBParameters, NONBONDACt, i )->dA );
    }
    FortranEndLine();
    
        /* -20- Lennard jones r**6 term for all possible interactions */
                /* CN2 array */
    FortranDebug( "-20-" );

    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iVarArrayElementCount(vaNBParameters); i++ ) {
        FortranWriteDouble( PVAI( vaNBParameters, NONBONDACt, i )->dB );
    }
    FortranEndLine();
 
        /* -21- Write the bond interactions that include hydrogen */
                /* Write the two indices into the atom table, then the index */
                /* into the interaction table */
    FortranDebug( "-21-" );

    MESSAGE( "Writing the bond interactions with hydrogens\n" );
    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount( uUnit->vaBonds ); i++ ) {
        sbPBond = PVAI( uUnit->vaBonds, SAVEBONDt, i );
        aA = PVAI( uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom1-1 )->aAtom;
        aD = PVAI( uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom2-1 )->aAtom;
        if ( bPERT_BOND(bPert,aA,aD) ) 
            continue;
        if ( iAtomElement(aA) == HYDROGEN ||
                iAtomElement(aD) == HYDROGEN ) {
            FortranWriteInt( AMBERINDEX(sbPBond->iAtom1) );
            FortranWriteInt( AMBERINDEX(sbPBond->iAtom2) );
            FortranWriteInt( sbPBond->iParmIndex );
        }
    }
    FortranEndLine();

        /* -22- Write the bond interactions that dont include hydrogen */
                /* Write the two indices into the atom table, then the index */
                /* into the interaction table */
    FortranDebug( "-22-" );

    MESSAGE( "Writing the bond interactions without hydrogens\n" );
    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount( uUnit->vaBonds ); i++ ) {
        sbPBond = PVAI( uUnit->vaBonds, SAVEBONDt, i );
        aA = PVAI( uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom1-1 )->aAtom;
        aD = PVAI( uUnit->vaAtoms, SAVEATOMt, sbPBond->iAtom2-1 )->aAtom;
        if ( bPERT_BOND(bPert,aA,aD) ) 
            continue;
        if ( !(iAtomElement(aA) == HYDROGEN ||
                iAtomElement(aD) == HYDROGEN) ) {
            FortranWriteInt( AMBERINDEX(sbPBond->iAtom1) );
            FortranWriteInt( AMBERINDEX(sbPBond->iAtom2) );
            FortranWriteInt( sbPBond->iParmIndex );
        }
    }
        /* Write out the (bond without H) RESTRAINT interactions */
        /* The iParmIndex field is set in RESTRAINTLOOP */
    if ( (iMax = iVarArrayElementCount( uUnit->vaRestraints )) ) {
        srPRestraint = PVAI( uUnit->vaRestraints, SAVERESTRAINTt, 0 );
        for ( i=0; i<iMax; i++, srPRestraint++ ) {
                if ( srPRestraint->iType == RESTRAINTBOND ) {
                        FortranWriteInt( AMBERINDEX(srPRestraint->iAtom1) );
                        FortranWriteInt( AMBERINDEX(srPRestraint->iAtom2) );
                        FortranWriteInt( srPRestraint->iParmIndex );
                }
        }
    }
    FortranEndLine();

        /* -23- Write the angle interactions that include hydrogen */
                /* Write the three indices into the atom table, then the index*/
                /* into the interaction table */
    FortranDebug( "-23-" );

    MESSAGE( "Writing the angle interactions with hydrogens\n" );
    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount( uUnit->vaAngles ); i++ ) {
        saPAngle = PVAI( uUnit->vaAngles, SAVEANGLEt, i );
        aA = PVAI( uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom1-1 )->aAtom;
        aB = PVAI( uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom2-1 )->aAtom;
        aD = PVAI( uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom3-1 )->aAtom;
        if ( bPERT_ANGLE(bPert,aA,aB,aD) ) 
            continue;
        if ( iAtomElement(aA) == HYDROGEN
                || iAtomElement(aB) == HYDROGEN
                || iAtomElement(aD) == HYDROGEN ) {
            FortranWriteInt( AMBERINDEX(saPAngle->iAtom1) );
            FortranWriteInt( AMBERINDEX(saPAngle->iAtom2) );
            FortranWriteInt( AMBERINDEX(saPAngle->iAtom3) );
            FortranWriteInt( saPAngle->iParmIndex );
        }
    }
    FortranEndLine();

        /* -24- Write the angle interactions that dont include hydrogen */
                /* Write the three indices into the atom table, then the index*/
                /* into the interaction table */
    FortranDebug( "-24-" );

    MESSAGE( "Writing the angle interactions without hydrogens\n" );
    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount( uUnit->vaAngles ); i++ ) {
        saPAngle = PVAI( uUnit->vaAngles, SAVEANGLEt, i );
        aA = PVAI( uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom1-1 )->aAtom;
        aB = PVAI( uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom2-1 )->aAtom;
        aD = PVAI( uUnit->vaAtoms, SAVEATOMt, saPAngle->iAtom3-1 )->aAtom;
        if ( bPERT_ANGLE(bPert,aA,aB,aD) ) 
            continue;
        if ( !(iAtomElement(aA) == HYDROGEN 
                || iAtomElement(aB) == HYDROGEN 
                || iAtomElement(aD) == HYDROGEN) ) {
            FortranWriteInt( AMBERINDEX(saPAngle->iAtom1) );
            FortranWriteInt( AMBERINDEX(saPAngle->iAtom2) );
            FortranWriteInt( AMBERINDEX(saPAngle->iAtom3) );
            FortranWriteInt( saPAngle->iParmIndex );
        }
    }
        /* Write out the RESTRAINT interactions */
        /* The iParmIndex field is set in RESTRAINTLOOP */
    if ( (iMax = iVarArrayElementCount( uUnit->vaRestraints )) ) {
        srPRestraint = PVAI( uUnit->vaRestraints, SAVERESTRAINTt, 0 );
        for ( i=0; i<iMax; i++, srPRestraint++ ) {
                if ( srPRestraint->iType == RESTRAINTANGLE ) {
                        FortranWriteInt( AMBERINDEX(srPRestraint->iAtom1) );
                        FortranWriteInt( AMBERINDEX(srPRestraint->iAtom2) );
                        FortranWriteInt( AMBERINDEX(srPRestraint->iAtom3) );
                        FortranWriteInt( srPRestraint->iParmIndex );
                }
        }
    }
    FortranEndLine();

        /* -25- Write the torsion interactions that include hydrogen */
               /* Write the three indices into the atom table, then the index*/
               /* into the interaction table */
    FortranDebug( "-25-" );

    MESSAGE( "Writing the torsion interactions with hydrogens\n" );
    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount( uUnit->vaTorsions ); i++ ) {
        stPTorsion = PVAI( uUnit->vaTorsions, SAVETORSIONt, i );
        aA = PVAI( uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom1-1 )->aAtom;
        aB = PVAI( uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom2-1 )->aAtom;
        aC = PVAI( uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom3-1 )->aAtom;
        aD = PVAI( uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom4-1 )->aAtom;
        if ( bPERT_TORSION(bPert,aA,aB,aC,aD) ) 
            continue;
        if ( iAtomElement(aA) == HYDROGEN 
                || iAtomElement(aB) == HYDROGEN 
                || iAtomElement(aC) == HYDROGEN
                || iAtomElement(aD) == HYDROGEN ) {
            if ( (AMBERINDEX(stPTorsion->iAtom3) == 0) ||
                 (AMBERINDEX(stPTorsion->iAtom4) == 0) ) {
                MESSAGE( "Had to turn torsion around to avoid K,L == 0\n" );
                MESSAGE( "Outer atoms: %s --- %s\n", 
                                sContainerName(aA), sContainerName(aD) );
                MESSAGE( "Old order %d %d %d %d\n",
                                stPTorsion->iAtom1,
                                stPTorsion->iAtom2,
                                stPTorsion->iAtom3,
                                stPTorsion->iAtom4 );
                SWAP( stPTorsion->iAtom1, stPTorsion->iAtom4, iTemp );
                SWAP( stPTorsion->iAtom2, stPTorsion->iAtom3, iTemp );
                MESSAGE( "New order %d %d %d %d\n",
                                stPTorsion->iAtom1,
                                stPTorsion->iAtom2,
                                stPTorsion->iAtom3,
                                stPTorsion->iAtom4 );
            }
            if ( stPTorsion->bProper )  iProper = 1;
            else                        iProper = -1;
            if ( stPTorsion->bCalc14 )  iCalc14 = 1;
            else                        iCalc14 = -1;
            if( GDefaults.iCharmm && iProper == -1 ){
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom3) );
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom2) );
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom1)*iCalc14 );
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom4)*iProper );
            } else {
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom1) );
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom2) );
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom3)*iCalc14 );
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom4)*iProper );
            }
            FortranWriteInt( stPTorsion->iParmIndex );
        }
    }
    FortranEndLine();

        /* -26- Write the torsion interactions that dont include hydrogen */
                /* Write the three indices into the atom table, then the index*/
                /* into the interaction table */
    FortranDebug( "-26-" );

    MESSAGE( "Writing the torsion interactions without hydrogens\n" );
    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount( uUnit->vaTorsions ); i++ ) {
        stPTorsion = PVAI( uUnit->vaTorsions, SAVETORSIONt, i );
        aA = PVAI( uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom1-1 )->aAtom;
        aB = PVAI( uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom2-1 )->aAtom;
        aC = PVAI( uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom3-1 )->aAtom;
        aD = PVAI( uUnit->vaAtoms, SAVEATOMt, stPTorsion->iAtom4-1 )->aAtom;
        if ( bPERT_TORSION(bPert,aA,aB,aC,aD) ) 
            continue;
        if ( !(iAtomElement(aA) == HYDROGEN 
                || iAtomElement(aB) == HYDROGEN 
                || iAtomElement(aC) == HYDROGEN
                || iAtomElement(aD) == HYDROGEN) ) {
            if ( (AMBERINDEX(stPTorsion->iAtom3) == 0) ||
                 (AMBERINDEX(stPTorsion->iAtom4) == 0) ) {
                MESSAGE( "Had to turn torsion to avoid K,L == 0\n" );
                MESSAGE( "Outer atoms: %s --- %s\n", 
                                sContainerName(aA), sContainerName(aD) );
                SWAP( stPTorsion->iAtom1, stPTorsion->iAtom4, iTemp );
                SWAP( stPTorsion->iAtom2, stPTorsion->iAtom3, iTemp );
            }
            if ( stPTorsion->bCalc14 )  iCalc14 = 1;
            else                        iCalc14 = -1;
            if ( stPTorsion->bProper )  iProper = 1;
            else                        iProper = -1;
            if( GDefaults.iCharmm && iProper == -1 ){
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom3) );
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom2) );
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom1)*iCalc14 );
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom4)*iProper );
            } else {
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom1) );
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom2) );
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom3)*iCalc14 );
              FortranWriteInt( AMBERINDEX(stPTorsion->iAtom4)*iProper );
            }
            FortranWriteInt( stPTorsion->iParmIndex );
        }
    }
        /* Write out the RESTRAINT interactions */
        /* The iParmIndex field is set in RESTRAINTLOOP */
    if ( (iMax = iVarArrayElementCount( uUnit->vaRestraints )) ) {
        srPRestraint = PVAI( uUnit->vaRestraints, SAVERESTRAINTt, 0 );
        for ( i=0; i<iMax; i++, srPRestraint++ ) {
            if ( srPRestraint->iType == RESTRAINTTORSION ) {
                if ( (AMBERINDEX(srPRestraint->iAtom3) == 0 ) ||
                     (AMBERINDEX(srPRestraint->iAtom4) == 0 ) ) {
                    MESSAGE( "Had to turn RESTRAINT torsion around to avoid\n" );
                    MESSAGE( "K,L == 0\n" );
                    SWAP( srPRestraint->iAtom1, srPRestraint->iAtom4, iTemp );
                    SWAP( srPRestraint->iAtom2, srPRestraint->iAtom3, iTemp );
                }
                FortranWriteInt( AMBERINDEX(srPRestraint->iAtom1) );
                FortranWriteInt( AMBERINDEX(srPRestraint->iAtom2) );
                FortranWriteInt( AMBERINDEX(srPRestraint->iAtom3) );
                FortranWriteInt( AMBERINDEX(srPRestraint->iAtom4) );
                FortranWriteInt( srPRestraint->iParmIndex );
            }
        }
    }
    FortranEndLine();

        /* -27- Write the excluded atom list */
    FortranDebug( "-27-" );

    MESSAGE( "Writing the excluded atom list\n" );
    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount( vaExcludedAtoms ); i++ ) {
        FortranWriteInt( *PVAI( vaExcludedAtoms, int, i ) );
    }
    FortranEndLine();

        /* -28- Write the R^12 term for the Hydrogen bond equation */
    FortranDebug( "-28-" );

    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iParmSetTotalHBondParms(uUnit->psParameters); i++ ) {
        ParmSetHBond( uUnit->psParameters, i, sType1, sType2, &dC, &dD, sDesc );
        FortranWriteDouble( dC );
    }
    FortranEndLine();

        /* -29- Write the R^10 term for the Hydrogen bond equation */
    FortranDebug( "-29-" );

    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iParmSetTotalHBondParms(uUnit->psParameters); i++ ) {
        ParmSetHBond( uUnit->psParameters, i, sType1, sType2, &dC, &dD, sDesc );
        FortranWriteDouble( dD );
    }
    FortranEndLine();

        /* -30- No longer used, but stored */
    FortranDebug( "-30-" );

    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iParmSetTotalHBondParms(uUnit->psParameters); i++ ) {
        FortranWriteDouble( 0.0 );
    }
    FortranEndLine();

        /* -31- List of atomic symbols */
    FortranDebug( "-31-" );

    FortranFormat( 20, LBLFORMAT );
    for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
        FortranWriteString( 
                sAtomType( PVAI(uUnit->vaAtoms, SAVEATOMt, i )->aAtom ) );
    }
    FortranEndLine();

        /* -32- List of tree symbols */
    FortranDebug( "-32-" );

    FortranFormat( 20, LBLFORMAT );
    for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
        aAtom = PVAI( uUnit->vaAtoms, SAVEATOMt, i )->aAtom;
        if ( dAtomTemp( aAtom ) == (double)'M' )
                FortranWriteString( "M  " );
        else if ( dAtomTemp( aAtom ) == (double)'E' )
                FortranWriteString( "E  " );
        else if ( dAtomTemp( aAtom ) == (double)'S' )
                FortranWriteString( "S  " );
        else if ( dAtomTemp( aAtom ) == (double)'B' )
                FortranWriteString( "B  " );
        else if ( dAtomTemp( aAtom ) == (double)'3' )
                FortranWriteString( "3  " );
        else if ( dAtomTemp( aAtom ) == (double)'4' )
                FortranWriteString( "4  " );
        else if ( dAtomTemp( aAtom ) == (double)'5' )
                FortranWriteString( "5  " );
        else if ( dAtomTemp( aAtom ) == (double)'6' )
                FortranWriteString( "6  " );
        else if ( dAtomTemp( aAtom ) == (double)'X' )
                FortranWriteString( "X  " );
        else
                FortranWriteString( "BLA" );
    }
    FortranEndLine();

        /* -33- Tree Joining information !!!!!!! Add support for this !!!!! */
    FortranDebug( "-33-" );

    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
        FortranWriteInt( 0 );
    }
    FortranEndLine();

        /* -34- Who knows, something to do with rotating atoms */
    FortranDebug( "-34-" );

    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
        FortranWriteInt( 0 );
    }
    FortranEndLine();

        /* -35A- The last residue before "solvent" */
                /* Number of molecules */
                /* Index of first molecule that is solvent */

    if ( bUnitUseBox(uUnit) ) {
        FortranDebug( "-35A-" );

                /* Find the index of the first solvent RESIDUE */

        for ( i=0; i<iVarArrayElementCount(uUnit->vaResidues); i++ ) {
            if ( PVAI(uUnit->vaResidues,SAVERESIDUEt,i)->sResidueType[0] ==
                RESTYPESOLVENT ) break;
        }
        int IPTRES = i;

        /* 
         *  Find the molecules and return the number of ATOMs in each 
         *  molecule, along with the index of the first solvent molecule
         */

        UnitIOFindAndCountMolecules( uUnit );
        int NSPM = iVarArrayElementCount(uUnit->vaAtomsPerMolecule);

        FortranFormat( 3, INTFORMAT );
        FortranWriteInt( IPTRES );
        FortranWriteInt( NSPM );
        FortranWriteInt( uUnit->iFirstSolvent+1 );     /* NSPSOL, FORTRAN index */

        FortranEndLine();

                /* -35B- The number of ATOMs in the Ith RESIDUE */

        FortranDebug( "-35B-" );
        FortranFormat( 12, INTFORMAT );
        for ( i=0; i<NSPM; i++ ) {
            FortranWriteInt( *PVAI(uUnit->vaAtomsPerMolecule,int,i) );
        }
        FortranEndLine();

                /* -35C- BETA, (BOX(I), I=1,3 ) */

        FortranDebug( "-35C-" );
        FortranFormat( 4, DBLFORMAT );
        FortranWriteDouble( dUnitBeta(uUnit)/DEGTORAD );
        UnitGetBox( uUnit, &dX, &dY, &dZ );
        FortranWriteDouble( dX );
        FortranWriteDouble( dY );
        FortranWriteDouble( dZ );
        FortranEndLine();
    }

        /* -35D- NATCAP */

    if ( bUnitUseSolventCap(uUnit) ) {
        FortranDebug( "-35D-" );
        FortranFormat( 1, INTFORMAT );
        FortranWriteInt( uUnit->iCapTempInt );
        FortranEndLine();
        
        /* -35E- CUTCAP, XCAP, YCAP, ZCAP */
        FortranDebug( "-35E-" );

        FortranFormat( 4, DBLFORMAT );
        UnitGetSolventCap( uUnit, &dX, &dY, &dZ, &dR );
        FortranWriteDouble( dR );
        FortranWriteDouble( dX );
        FortranWriteDouble( dY );
        FortranWriteDouble( dZ );
        FortranEndLine();
    }

                /* Write the perturbation information */

    if ( bPert ) {

                /* -36A- Bonds that are to be perturbed */
                        /* Totally perturbed bonds first, */
                        /* boundary second */
        FortranDebug( "-36A-" );

        FortranFormat( 12, INTFORMAT );
        if ( iVarArrayElementCount(uUnit->vaBonds) ) {
                sbPBond = PVAI( uUnit->vaBonds, SAVEBONDt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaBonds); i++, sbPBond++ ) {
            if ( ((sbPBond->fFlags&PERTURBED)!=0)
                 && ((sbPBond->fFlags&BOUNDARY)==0) ) {
                FortranWriteInt( AMBERINDEX(sbPBond->iAtom1) );
                FortranWriteInt( AMBERINDEX(sbPBond->iAtom2) );
            }
        }
                sbPBond = PVAI( uUnit->vaBonds, SAVEBONDt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaBonds); i++, sbPBond++ ) {
            if ( ((sbPBond->fFlags&PERTURBED)!=0)
                 && ((sbPBond->fFlags&BOUNDARY)!=0) ) {
                FortranWriteInt( AMBERINDEX(sbPBond->iAtom1) );
                FortranWriteInt( AMBERINDEX(sbPBond->iAtom2) );
            }
        }}
        FortranEndLine();

                /* -36B- Index into bond interaction arrays */
        FortranDebug( "-36B-" );

        FortranFormat( 12, INTFORMAT );


                        /* First LAMBDA = 0 */

        if ( iVarArrayElementCount(uUnit->vaBonds) ) {
                sbPBond = PVAI( uUnit->vaBonds, SAVEBONDt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaBonds); i++, sbPBond++ ) {
            if ( ((sbPBond->fFlags&PERTURBED)!=0)
                && ((sbPBond->fFlags&BOUNDARY)==0) ) {
                FortranWriteInt( sbPBond->iParmIndex );
            }
        }
                sbPBond = PVAI( uUnit->vaBonds, SAVEBONDt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaBonds); i++, sbPBond++ ) {
            if ( ((sbPBond->fFlags&PERTURBED)!=0)
                && ((sbPBond->fFlags&BOUNDARY)!=0) ) {
                FortranWriteInt( sbPBond->iParmIndex );
            }
        }

                        /* Then LAMBDA = 1 */

                sbPBond = PVAI( uUnit->vaBonds, SAVEBONDt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaBonds); i++, sbPBond++ ) {
            if ( ((sbPBond->fFlags&PERTURBED)!=0)
                && ((sbPBond->fFlags&BOUNDARY)==0) ) {
                FortranWriteInt( sbPBond->iPertParmIndex );
            }
        }
                sbPBond = PVAI( uUnit->vaBonds, SAVEBONDt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaBonds); i++, sbPBond++ ) {
            if ( ((sbPBond->fFlags&PERTURBED)!=0)
                && ((sbPBond->fFlags&BOUNDARY)!=0) ) {
                FortranWriteInt( sbPBond->iPertParmIndex );
            }
        }}
        FortranEndLine();
        

                /* -36C- Angles that are to be perturbed */
        FortranDebug( "-36C-" );

        FortranFormat( 12, INTFORMAT );
        if ( iVarArrayElementCount(uUnit->vaAngles) ) {
                saPAngle = PVAI( uUnit->vaAngles, SAVEANGLEt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAngles); i++, saPAngle++ ) {
            if ( ((saPAngle->fFlags&PERTURBED)!=0)
                && ((saPAngle->fFlags&BOUNDARY)==0) ) {
                FortranWriteInt( AMBERINDEX(saPAngle->iAtom1) );
                FortranWriteInt( AMBERINDEX(saPAngle->iAtom2) );
                FortranWriteInt( AMBERINDEX(saPAngle->iAtom3) );
            }
        }
                saPAngle = PVAI( uUnit->vaAngles, SAVEANGLEt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAngles); i++, saPAngle++ ) {
            if ( ((saPAngle->fFlags&PERTURBED)!=0)
                && ((saPAngle->fFlags&BOUNDARY)!=0) ) {
                FortranWriteInt( AMBERINDEX(saPAngle->iAtom1) );
                FortranWriteInt( AMBERINDEX(saPAngle->iAtom2) );
                FortranWriteInt( AMBERINDEX(saPAngle->iAtom3) );
            }
        }}
        FortranEndLine();

                /* -36D- Index into angle interaction arrays */
        FortranDebug( "-36D-" );

        FortranFormat( 12, INTFORMAT );


                        /* First LAMBDA = 0 */

        if ( iVarArrayElementCount(uUnit->vaAngles) ) {
                saPAngle = PVAI( uUnit->vaAngles, SAVEANGLEt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAngles); i++, saPAngle++ ) {
            if ( ((saPAngle->fFlags&PERTURBED)!=0)
                && ((saPAngle->fFlags&BOUNDARY)==0) ) {
                FortranWriteInt( saPAngle->iParmIndex );
            }
        }
                saPAngle = PVAI( uUnit->vaAngles, SAVEANGLEt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAngles); i++, saPAngle++ ) {
            if ( ((saPAngle->fFlags&PERTURBED)!=0)
                && ((saPAngle->fFlags&BOUNDARY)!=0) ) {
                FortranWriteInt( saPAngle->iParmIndex );
            }
        }

                        /* Then LAMBDA = 1 */

                saPAngle = PVAI( uUnit->vaAngles, SAVEANGLEt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAngles); i++, saPAngle++ ) {
            if ( ((saPAngle->fFlags&PERTURBED)!=0)
                && ((saPAngle->fFlags&BOUNDARY)==0) ) {
                FortranWriteInt( saPAngle->iPertParmIndex );
            }
        }
                saPAngle = PVAI( uUnit->vaAngles, SAVEANGLEt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAngles); i++, saPAngle++ ) {
            if ( ((saPAngle->fFlags&PERTURBED)!=0)
                && ((saPAngle->fFlags&BOUNDARY)!=0) ) {
                FortranWriteInt( saPAngle->iPertParmIndex );
            }
        }}
        FortranEndLine();

                /* -36E- Torsions that are to be perturbed */
        FortranDebug( "-36E-" );

        FortranFormat( 12, INTFORMAT );
        if ( iVarArrayElementCount(uUnit->vaTorsions) ) {
                stPTorsion = PVAI( uUnit->vaTorsions, SAVETORSIONt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaTorsions); 
                                                        i++, stPTorsion++ ) {
            if ( ((stPTorsion->fFlags&PERTURBED)!=0)
                && ((stPTorsion->fFlags&BOUNDARY)==0) ) {

                if ( (AMBERINDEX(stPTorsion->iAtom3) == 0) ||
                     (AMBERINDEX(stPTorsion->iAtom4) == 0) ) {
                    MESSAGE( "Had to turn torsion around to avoid K,L == 0\n" );
                    MESSAGE( "Outer atoms: %s --- %s\n", 
                                    sContainerName(aA), sContainerName(aD) );
                    SWAP( stPTorsion->iAtom1, stPTorsion->iAtom4, iTemp );
                    SWAP( stPTorsion->iAtom2, stPTorsion->iAtom3, iTemp );
                }
                if ( stPTorsion->bProper )  iProper = 1;
                else                        iProper = -1;
                if ( stPTorsion->bCalc14 )  iCalc14 = 1;
                else                        iCalc14 = -1;
                FortranWriteInt( AMBERINDEX(stPTorsion->iAtom1) );
                FortranWriteInt( AMBERINDEX(stPTorsion->iAtom2) );
                FortranWriteInt( AMBERINDEX(stPTorsion->iAtom3)*iCalc14 );
                FortranWriteInt( AMBERINDEX(stPTorsion->iAtom4)*iProper );
            }
        }
        stPTorsion = PVAI( uUnit->vaTorsions, SAVETORSIONt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaTorsions); 
                                                        i++, stPTorsion++ ) {
            if ( ((stPTorsion->fFlags&PERTURBED)!=0)
                && ((stPTorsion->fFlags&BOUNDARY)!=0) ) {

                if ( (AMBERINDEX(stPTorsion->iAtom3) == 0) ||
                     (AMBERINDEX(stPTorsion->iAtom4) == 0) ) {
                    MESSAGE( "Had to turn torsion around to avoid K,L == 0\n" );
                    MESSAGE( "Outer atoms: %s --- %s\n", 
                                    sContainerName(aA), sContainerName(aD) );
                    SWAP( stPTorsion->iAtom1, stPTorsion->iAtom4, iTemp );
                    SWAP( stPTorsion->iAtom2, stPTorsion->iAtom3, iTemp );
                }
                if ( stPTorsion->bProper )  iProper = 1;
                else                        iProper = -1;
                if ( stPTorsion->bCalc14 )  iCalc14 = 1;
                else                        iCalc14 = -1;
                FortranWriteInt( AMBERINDEX(stPTorsion->iAtom1) );
                FortranWriteInt( AMBERINDEX(stPTorsion->iAtom2) );
                FortranWriteInt( AMBERINDEX(stPTorsion->iAtom3)*iCalc14 );
                FortranWriteInt( AMBERINDEX(stPTorsion->iAtom4)*iProper );
            }
        }}
        FortranEndLine();

                /* -36F- Index into torsion interaction arrays */
        FortranDebug( "-36F-" );

        FortranFormat( 12, INTFORMAT );

                        /* First LAMBDA = 0 */

        if ( iVarArrayElementCount(uUnit->vaTorsions) ) {
                stPTorsion = PVAI( uUnit->vaTorsions, SAVETORSIONt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaTorsions);
                                                        i++, stPTorsion++ ) {
            if ( ((stPTorsion->fFlags&PERTURBED)!=0)
                && ((stPTorsion->fFlags&BOUNDARY)==0) ) {
                FortranWriteInt( stPTorsion->iParmIndex );
            }
        }
                stPTorsion = PVAI( uUnit->vaTorsions, SAVETORSIONt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaTorsions); 
                                                        i++, stPTorsion++ ) {
            if ( ((stPTorsion->fFlags&PERTURBED)!=0)
                && ((stPTorsion->fFlags&BOUNDARY)!=0) ) {
                FortranWriteInt( stPTorsion->iParmIndex );
            }
        }

                        /* Then LAMBDA = 1 */

                stPTorsion = PVAI( uUnit->vaTorsions, SAVETORSIONt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaTorsions); 
                                                        i++, stPTorsion++ ) {
            if ( ((stPTorsion->fFlags&PERTURBED)!=0)
                && ((stPTorsion->fFlags&BOUNDARY)==0) ) {
                FortranWriteInt( stPTorsion->iPertParmIndex );
            }
        }
                stPTorsion = PVAI( uUnit->vaTorsions, SAVETORSIONt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaTorsions);
                                                        i++, stPTorsion++ ) {
            if ( ((stPTorsion->fFlags&PERTURBED)!=0)
                && ((stPTorsion->fFlags&BOUNDARY)!=0) ) {
                FortranWriteInt( stPTorsion->iPertParmIndex );
            }
        }}
        FortranEndLine();

                /* -36G- Residue labels at LAMBDA = 1 */
                        /* Just write the labels at LAMBDA = 0 */

                /* Trim the string down to at most 3 characters by */
                /* taking the last three characters if it is too long */
        FortranDebug( "-36G-" );

        MESSAGE( "Writing the residue labels\n" );
        FortranFormat( 20, LBLFORMAT );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaResidues); i++ ) {
            cPTemp = PVAI( uUnit->vaResidues, SAVERESIDUEt, i )->sName;
            if ( strlen( cPTemp ) > 3 ) cPTemp += ( strlen(cPTemp)-3 );
            FortranWriteString( cPTemp );
        }
        FortranEndLine();

                /* -36H- Atom names at LAMBDA = 0 */
        FortranDebug( "-36H-" );

        FortranFormat( 20, LBLFORMAT );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
            cPTemp = PVAI( uUnit->vaAtoms, SAVEATOMt, i )->sPertName;
            if ( strlen( cPTemp ) == 0 )
                cPTemp = PVAI( uUnit->vaAtoms, SAVEATOMt, i )->sName;
            if ( strlen( cPTemp ) > 4 ) cPTemp += ( strlen(cPTemp)-4 );
            FortranWriteString( cPTemp );
        }
        FortranEndLine();

                /* -36I- List of atomic symbols (atom types??????) */
        FortranDebug( "-36I-" );

        FortranFormat( 20, LBLFORMAT );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
            cPTemp = sAtomPertType(PVAI(uUnit->vaAtoms,SAVEATOMt,i)->aAtom);
            if ( strlen(cPTemp) == 0 )
                cPTemp = sAtomType(PVAI(uUnit->vaAtoms,SAVEATOMt,i)->aAtom);
            if ( strlen( cPTemp ) > 3 ) cPTemp += ( strlen(cPTemp)-3 );
            FortranWriteString( cPTemp );
        }
        FortranEndLine();

                /* -36J- Value of LAMBDA for each ATOM ????????? */
                /* TODO: Figure out what the hell this is */
        FortranDebug( "-36J-" );

        FortranFormat( 5, DBLFORMAT );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
            FortranWriteDouble( 0.0 );
        }
        FortranEndLine();

                /* -36K- Flag to tell whether the atom is perturbed */
        FortranDebug( "-36K-" );

        FortranFormat( 12, INTFORMAT );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
            if ( bAtomPerturbed(PVAI(uUnit->vaAtoms,SAVEATOMt,i)->aAtom) ) 
                FortranWriteInt( 1 );
            else FortranWriteInt( 0 );
        }
        FortranEndLine();

                /* -36L- List of atom types - IACPER */
        FortranDebug( "-36L-" );

        FortranFormat( 12, INTFORMAT );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
            iAtom = PVAI( uUnit->vaAtoms, SAVEATOMt, i )->iPertTypeIndex-1;
            iTemp = *PVAI( vaNBIndex, int, iAtom );
            FortranWriteInt( iTemp+1 );
        }
        FortranEndLine();

                /* -36M- Perturbed charges */
        FortranDebug( "-36M-" );

        FortranFormat( 5, DBLFORMAT );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
            aAtom = PVAI(uUnit->vaAtoms,SAVEATOMt,i)->aAtom;
            if ( bAtomPerturbed( aAtom ) )
                FortranWriteDouble( ELECTRONTOKCAL *
                    dAtomPertCharge(PVAI(uUnit->vaAtoms,SAVEATOMt,i)->aAtom));
            else
                FortranWriteDouble( ELECTRONTOKCAL *
                    dAtomCharge(PVAI(uUnit->vaAtoms,SAVEATOMt,i)->aAtom));
        }
        FortranEndLine();

    }

    /*
     *  polarizabilities
     */
    if ( bPolar ) {
        iCount = 0;
        iCountPerturbed = 0;
        iMax = iVarArrayElementCount(uUnit->vaAtoms);
        MESSAGE( "Writing the atomic polarizabilities\n" );         
        FortranFormat( 5, DBLFORMAT );
        saPAtom = PVAI( uUnit->vaAtoms, SAVEATOMt, 0 );
        for ( i=0; i<iMax; i++, saPAtom++ ) {
            iIndex = iParmSetFindAtom( uUnit->psParameters, saPAtom->sType );
            ParmSetAtom( uUnit->psParameters, iIndex, sType,
                        &dMass, &dPolar, &dEpsilon, &dRStar, &dEpsilon14, &dRStar14,
			&dScreenF,
            &iElement, &iHybridization, sDesc );
            if ( dPolar == -1.0 ) {
                dPolar = 0.0;
                iCount++;
            }
            FortranWriteDouble( dPolar );
        }
        if ( iCount > 0 )
                VP0( "Total atoms with default polarization=0.0: %d of %d\n",
                                                        iCount, iMax);
                                                
	/*
        FortranEndLine();
	if (GDefaults.dDipoleDampFactor > 1.0) {
           FortranFormat(1, "%-80s");
           FortranWriteString("%FLAG DIPOLE_DAMP_FACTOR");
           FortranWriteString("%FORMAT(5E16.8)");
           FortranFormat(5, DBLFORMAT);
           for (i = 0; i < iMax; i++, saPAtom++) {
               iIndex = iParmSetFindAtom(uUnit->psParameters, saPAtom->sType);
               ParmSetAtom(uUnit->psParameters, iIndex, sType,
                           &dMass, &dPolar, &dEpsilon, &dRStar, &dEpsilon14,
                           &dRStar14, &dScreenF, &iElement, &iHybridization,
			   sDesc); 
	       if (dScreenF == 0.0) {
		  dScreenF = GDefaults.dDipoleDampFactor;
	       }
               FortranWriteDouble(dScreenF);
	   }                       
	}
	*/

        FortranEndLine();
        if ( bPert ) {
                int     iPertTot = 0;
                saPAtom = PVAI( uUnit->vaAtoms, SAVEATOMt, 0 );
                for ( i=0; i<iMax; i++, saPAtom++ ) {
                        bool    bTmp = bAtomPerturbed( saPAtom->aAtom );

                        if ( bTmp ) {
                                iIndex = iParmSetFindAtom( uUnit->psParameters, 
                                                        saPAtom->sPertType );
                                iPertTot++;
                        } else
                                iIndex = iParmSetFindAtom( uUnit->psParameters, 
                                                        saPAtom->sType );
                        ParmSetAtom( uUnit->psParameters, iIndex, sType,
                                        &dMass, &dPolar, &dEpsilon, &dRStar, &dEpsilon14,
                    &dRStar14, &dScreenF, &iElement, &iHybridization, sDesc );
                        if ( dPolar == -1.0 ) {
                                dPolar = 0.0;
                                if ( bTmp ) iCountPerturbed++;
                        }
                        FortranWriteDouble( dPolar );
                }
                if ( iCountPerturbed > 0 )
                    VP0( "Total pert atoms with default polarization=0.0: %d of %d\n",
                                        iCountPerturbed, iPertTot );
        }
    }

        /*  Charmm-style parameters  */

        if( GDefaults.iCharmm ){
        /* -19- Lennard jones r**12 term for all 14 interactions */
                /* CN114 array */
        FortranDebug( "-19-" );

        FortranFormat( 5, DBLFORMAT );
        for ( i=0; i<iVarArrayElementCount(vaNBParameters); i++ ) {
                FortranWriteDouble( PVAI( vaNBParameters, NONBONDACt, i )->dA14 );
        }
        FortranEndLine();
    
        /* -20- Lennard jones r**6 term for all 14 interactions */
                /* CN214 array */
        FortranDebug( "-20-" );

        FortranFormat( 5, DBLFORMAT );
        for ( i=0; i<iVarArrayElementCount(vaNBParameters); i++ ) {
                FortranWriteDouble( PVAI( vaNBParameters, NONBONDACt, i )->dB14 );
        }
        FortranEndLine();

        /* -13- Force constants for Urey-Bradley */
                FortranDebug( "-13-" );

                FortranFormat( 5, DBLFORMAT );
                for ( i=0; i<iParmSetTotalAngleParms(uUnit->psParameters); i++ ) {
                        ParmSetAngle( uUnit->psParameters, i, sAtom1, sAtom2, sAtom3,
                                &dKt, &dT0, &dTkub, &dRkub, sDesc );
                        FortranWriteDouble( dTkub );
                }
                FortranEndLine();

        /* -14- Equilibrium distances for Urey-Bradley*/
                FortranDebug( "-14-" );

                FortranFormat( 5, DBLFORMAT );
                for ( i=0; i<iParmSetTotalAngleParms(uUnit->psParameters); i++ ) {
                        ParmSetAngle( uUnit->psParameters, i, sAtom1, sAtom2, sAtom3,
                                                &dKt, &dT0, &dTkub, &dRkub, sDesc );
                        FortranWriteDouble( dRkub );
                }
                FortranEndLine();
 
        }


        /********************************************************/
        /* Write the coordinate file                            */
        /********************************************************/

    FortranFile( fCrd );

    FortranFormat( 1, "%s" );
    FortranWriteString( sContainerName( uUnit ) );
    FortranEndLine();

    FortranFormat( 1, "%5d" );
    FortranWriteInt( iVarArrayElementCount( uUnit->vaAtoms ) );
    FortranEndLine();

    FortranFormat( 6, "%12.7lf" );
    if ( bUnitUseBox(uUnit) ) {
        double  dX2, dY2, dZ2;

        UnitGetBox( uUnit, &dX, &dY, &dZ );
        dX2 = dX * 0.5;
        dY2 = dY * 0.5;
        dZ2 = dZ * 0.5;

        /*
         *  shift box to Amber spot; later, add a cmd opt or environment
         *      var to switch between 0,0,0 center (spasms) or corner
         */
        for ( i = 0; i<iVarArrayElementCount( uUnit->vaAtoms ); i++ ) {
            vPos = PVAI( uUnit->vaAtoms, SAVEATOMt, i )->vPos;
            FortranWriteDouble( dVX(&vPos) + dX2 );
            FortranWriteDouble( dVY(&vPos) + dY2 );
            FortranWriteDouble( dVZ(&vPos) + dZ2 );
        }
        FortranEndLine();
        FortranWriteDouble( dX );
        FortranWriteDouble( dY );
        FortranWriteDouble( dZ );
        FortranWriteDouble( dUnitBeta( uUnit ) / DEGTORAD );
        FortranWriteDouble( dUnitBeta( uUnit ) / DEGTORAD );
        FortranWriteDouble( dUnitBeta( uUnit ) / DEGTORAD );
        FortranEndLine();
    } else {
        for ( i = 0; i<iVarArrayElementCount( uUnit->vaAtoms ); i++ ) {
            vPos = PVAI( uUnit->vaAtoms, SAVEATOMt, i )->vPos;
            FortranWriteDouble( dVX(&vPos) );
            FortranWriteDouble( dVY(&vPos) );
            FortranWriteDouble( dVZ(&vPos) );
        }
        FortranEndLine();
    }

    VarArrayDestroy( &vaNBIndexMatrix );
    VarArrayDestroy( &vaNBParameters );
    VarArrayDestroy( &vaExcludedAtoms );
    VarArrayDestroy( &vaExcludedCount );
    VarArrayDestroy( &vaNBIndex );
    VarArrayDestroy( &vaNonBonds );
    fclose( fCrd );

}

