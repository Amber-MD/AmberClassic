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

#define AMBERINDEX(i)   3*(i-1)
#define INTFORMAT_C     "%8d"
#define INTFORMAT_F     "I8"
#define DBLFORMAT_C     "%16.8lE"
#define DBLFORMAT_F     "E16.8"
#define LBLFORMAT_C     "%-4s"
#define LBLFORMAT_F     "a4"

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

    ContainerResetAllAtomsFlags((CONTAINER) uUnit, ATOMTOUCHED);

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
        lSpanning = lLoop((OBJEKT) aAtom, SPANNINGTREE);
        FOREACH(aAtom, ATOM, lSpanning) {
            AtomSetFlags(aAtom, ATOMTOUCHED);
            iCount++;
            CONTAINER cParent = cContainerWithin(aAtom);
            if (iObjectType(cParent) == RESIDUEid)
                ((RESIDUE)cParent)->iMolecule = iMolecule;
        }

        if (!bSeenFirstSolvent) {
            if (cResidueType(srPRes->rResidue) == RESTYPESOLVENT) {
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
            if (iHBondIndex != PARM_NOT_FOUND)
                iIndex = -(iHBondIndex + 1);
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
    STRING sAtom1, sAtom2, sAtom3, sAtom4, sType1, sType2;
    int iN, iAtoms, iMaxAtoms, iTemp, iAtom, iCalc14, iProper;
    int iElement, iHybridization, iStart;
    RESIDUE rRes;
    bool bFoundSome;
    char *cPTemp = NULL;
    double dX, dY, dZ, dEpsilon, dRStar, dEpsilon14, dRStar14;
    STRING sDesc, sType;
    IX_REC eResEnt = {NULL};
    IX_DESC iResIx;
    char sVersionHeader[81];
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
        lSpan = lLoop((OBJEKT) aAtom, SPANNINGTREE);
        iCount = 0;
        bFoundSome = FALSE;
        iStart = iVarArrayElementCount(vaExcludedAtoms);
        while ((aA = (ATOM) oNext(&lSpan))) {
            if (aA == aAtom) continue;
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
            SortByInteger((GENP) PVAI(vaExcludedAtoms, int, iStart), iCount,
                          sizeof(int), (GENP) PVAI(vaExcludedAtoms, int, iStart), TRUE);
        }
        VarArrayAdd(vaExcludedCount, (GENP) &iCount);
    }

    /**** Mark atom tree classification ***************************/
    create_index(&iResIx, IX_DUPKEYREC, IX_LEN_CSTRING);
    VP0("Marking per-residue atom chain types.\n");
    iMaxAtoms = 0;
    lTemp = lLoop((OBJEKT) uUnit, RESIDUES);
    while ((rRes = (RESIDUE) oNext(&lTemp))) {
        int iAtoms = iMarkMainChainAtoms(rRes, 0);
        if (iAtoms > 0) MarkSideChainAtoms(rRes);
        if (iAtoms < 0) {
            iAtoms = -iAtoms;
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

  // -1- Save the title of the UNIT
  FortranDebug("-1-");
  MESSAGE("Saving the name of the UNIT\n");
  FortranFormat(1, "%-80s");
  time(&tp);
  strftime(sVersionHeader, 81,
           "%%VERSION  VERSION_STAMP = V0001.000  DATE = %m/%d/%y  %H:%M:%S",
           localtime(&tp));
  FortranWriteString(sVersionHeader);
  FortranWriteString("%FLAG TITLE");
  FortranWriteString("%FORMAT(20" LBLFORMAT_F ")");
  FortranWriteString(sContainerName(uUnit));

/////////////////////////////////////////////////
  FortranFormat(1,"%-80s");
  FortranWriteString("%FLAG FORCE_FIELD_TYPE");
  FortranWriteString("%FORMAT(I,A)");
  { 
      PARMSET psSet;
      int iLines=0;
      ParmLibParmSetLoop(GplAllParameters);
      while ( bParmLibNextParmSet(GplAllParameters, &psSet) ) iLines += 2;
      ParmLibParmSetLoop(GplAllParameters);
      while ( bParmLibNextParmSet(GplAllParameters, &psSet) ) {
          STRING sTemp;
          snprintf(sTemp,sizeof(sTemp),"%d Filename: %.1024s\n",iLines,psSet->sFname);
          FortranWriteString(sTemp);
          snprintf(sTemp,sizeof(sTemp),"%d Title: %.1024s\n",iLines,psSet->sTitle);
          FortranWriteString(sTemp);
      }
  }
/////////////////////////////////////////////////

  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG POINTERS");
  FortranWriteString("%FORMAT(10" INTFORMAT_F ")");

  // -2- Save control information
  FortranDebug("-2-");
  MESSAGE("Saving all the main control variables\n");
  FortranFormat(10, INTFORMAT_C);

  // NTOTAT
  FortranWriteInt(iVarArrayElementCount(uUnit->vaAtoms));

  // NTYPES
  FortranWriteInt(iVarArrayElementCount(vaNonBonds));

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

  // NBONH
  FortranWriteInt(iBondWith);

  // NBONA
  FortranWriteInt(iBondWithout);

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

  // NTHETH
  FortranWriteInt(iAngleWith);

  // NTHETA
  FortranWriteInt(iAngleWithout);

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

  // NPHIH
  FortranWriteInt(iTorsionWith);

  // NPHIA
  FortranWriteInt(iTorsionWithout);

  // JHPARM
  FortranWriteInt(0);

  // JPARM
  FortranWriteInt(0);

  // Write the number of excluded atoms

  // NEXT
  FortranWriteInt(iVarArrayElementCount(vaExcludedAtoms));

  // NTOTRS
  FortranWriteInt(iVarArrayElementCount(uUnit->vaResidues));

  // Write the number of bonds/angles/torsions without hydrogens
  // PLUS the number of RESTRAINT bonds/angles/torsions

  // MBONA
  FortranWriteInt(iBondWithout + iUnitRestraintTypeCount(uUnit, RESTRAINTBOND));

  // MTHETA
  FortranWriteInt(iAngleWithout + iUnitRestraintTypeCount(uUnit, RESTRAINTANGLE));
  // MPHIA
  FortranWriteInt(iTorsionWithout + iUnitRestraintTypeCount(uUnit, RESTRAINTTORSION));

  // Write the number of unique bond types, angle types, and torsion types.  Add in the
  // number of RESTRAINT bonds/angles/torsion because they will have new parameters.

  // MUMBND
  FortranWriteInt(iParmSetTotalBondParms(uUnit->psParameters) +
                  iUnitRestraintTypeCount(uUnit, RESTRAINTBOND));
  // MUMANG
  FortranWriteInt(iParmSetTotalAngleParms(uUnit->psParameters) +
                  iUnitRestraintTypeCount(uUnit, RESTRAINTANGLE));
  // NPTRA
  FortranWriteInt(iParmSetTotalTorsionParms(uUnit->psParameters) +
                  iParmSetTotalImproperParms(uUnit->psParameters) +
                  iUnitRestraintTypeCount(uUnit, RESTRAINTTORSION));

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

  // NATYP
  FortranWriteInt(iParmSetTotalAtomParms(uUnit->psParameters));

  // NHB
  FortranWriteInt(iParmSetTotalHBondParms(uUnit->psParameters));

  // IFPERT
  if (bPert) {
    FortranWriteInt(1);
  }
  else {
    FortranWriteInt(0);
  }

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

  // NBPER
  FortranWriteInt(iCountBondPerturbed);

  // NGPER
  FortranWriteInt(iCountAnglePerturbed);

  // NDPER
  FortranWriteInt(iCountTorsionPerturbed);

  // MBPER
  FortranWriteInt(iCountBondPerturbed - iCountBondBoundary);

  // MGPER
  FortranWriteInt(iCountAnglePerturbed - iCountAngleBoundary);

  // MDPER
  FortranWriteInt(iCountTorsionPerturbed - iCountTorsionBoundary);

  // Save flag for periodic boundary conditions
  // IFBOX
  if (bUnitUseBox(uUnit)) {
    if (bUnitBoxOct(uUnit)) {
      FortranWriteInt(2);
    }
    else {
      FortranWriteInt(1);
    }
  }
  else {
    FortranWriteInt(0);
  }

  // Save the number of atoms in the largest residue

  // NMXRS
  FortranWriteInt(iMaxAtoms);

  // Save flag for cap information

  // IFCAP
  if (bUnitUseSolventCap(uUnit)) {
    FortranWriteInt(1);
  }
  else {
    FortranWriteInt(0);
  }

  // NUMEXTRA
  iNumExtra = 0;
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
    cPTemp = sAtomType(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom);
    if (!strncmp(cPTemp, "EP", 2)) {
      iNumExtra++;
    }
  }
  FortranWriteInt(iNumExtra);

  //NCOPY but not sure why it is not here. PMEMDCUDA reads it though.
  //FortranWriteInt(0);
  FortranEndLine();

  // -3-  write out the names of the atoms
  FortranDebug("-3-");
  MESSAGE("Writing the names of the atoms\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG ATOM_NAME");
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
  FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
  FortranFormat(5, DBLFORMAT_C);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
    FortranWriteDouble(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->dCharge * ELECTRONTOKCAL);
  }
  FortranEndLine();
  MESSAGE("Writing the atomic numbers\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG ATOMIC_NUMBER");
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
  FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
  FortranFormat(5, DBLFORMAT_C);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
    saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i);
    iIndex = iParmSetFindAtom(uUnit->psParameters, saPAtom->sType);
    ParmSetAtom(uUnit->psParameters, iIndex, sType, &dMass, &dPolar, &dEpsilon, &dRStar,
		&dEpsilon14, &dRStar14, &dScreenF, &iElement, &iHybridization, sDesc);
    FortranWriteDouble(dMass);
  }
  FortranEndLine();

  // -6- write out the atomic types
  FortranDebug("-6-");
  MESSAGE("Writing the atomic types\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG ATOM_TYPE_INDEX");
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

  MESSAGE("Writing the residue types\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG RESIDUE_TYPE");
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


  // -10- Pointer list for all the residues
  FortranDebug("-10-");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG RESIDUE_POINTER");
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
    iTemp = i;

    // Find the molecules and return the number of ATOMs in each
    // molecule, along with the index of the first solvent molecule
    //
    UnitIOFindAndCountMolecules( uUnit );
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG SOLVENT_POINTERS");
    FortranWriteString("%FORMAT(3I8)");
    FortranFormat(3, INTFORMAT_C);
    FortranWriteInt(iTemp);
    FortranWriteInt(iVarArrayElementCount(uUnit->vaAtomsPerMolecule));

    // FORTRAN indexing conversion
    FortranWriteInt(uUnit->iFirstSolvent + 1);
    FortranEndLine();

    // -35B- The number of ATOMs in the Ith MOLECULE
    FortranDebug("-35B-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG ATOMS_PER_MOLECULE");
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
    if (dUnitBeta(uUnit)!=uUnit->dAlpha ||
      dUnitBeta(uUnit)!=uUnit->dGamma) {
      FortranWriteString("%FLAG CELL_DIMENSIONS");
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
    FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
    FortranFormat(1, INTFORMAT_C);
    FortranWriteInt(uUnit->iCapTempInt);
    FortranEndLine();

    // -35E- CUTCAP, XCAP, YCAP, ZCAP
    FortranDebug("-35E-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG CAP_INFO2");
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
  FortranWriteString("%FORMAT(1a80)");
  switch (GDefaults.iGBparm) {
    case 0:
      FortranWriteString("Bondi radii (bondi)");
      break;
    case 1:
      FortranWriteString("amber6 modified Bondi radii (amber6)");
      break;
    case 2:
      FortranWriteString("modified Bondi radii (mbondi)");
      break;
#if 0
    case 3:
      FortranWriteString("Huo and Kollman optimized radii (pbamber)");
      break;
#endif
    case 6:
      FortranWriteString("H(N)-modified Bondi radii (mbondi2)");
      break;
    case 7:
      FortranWriteString("Parse radii (parse)");
       break;
    case 8:
      FortranWriteString("ArgH and AspGluO modified Bondi2 radii (mbondi3)");
      break;
    default:
      FortranWriteString("Unknown radius set (leap needs to be modified!)");
      break;
  }
  FortranEndLine();
  FortranWriteString("%FLAG RADII");
  FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
  FortranFormat(5, DBLFORMAT_C);
  saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, 0);
  {
      int n = iVarArrayElementCount(uUnit->vaAtoms);
      for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
        iIndex = iParmSetFindAtom(uUnit->psParameters, saPAtom->sType);
        ParmSetAtom(uUnit->psParameters, iIndex, sType, &dMass, &dPolar,
                    &dEpsilon, &dRStar, &dEpsilon14, &dRStar14, &dScreenF,
                    &iElement, &iHybridization, sDesc);
        FortranWriteDouble(dGBRadiusForAtom(&saPAtom[i], iElement, dMass, (i+1==n) ));
      }
      FortranEndLine();
  }

  // Write out the GB screening parameters
  MESSAGE("Writing the GB screening parameters\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG SCREEN");
  FortranWriteString("%FORMAT(5" DBLFORMAT_F ")");
  FortranFormat(5, DBLFORMAT_C);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
    saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i);
    iIndex = iParmSetFindAtom(uUnit->psParameters, saPAtom->sType);
    ParmSetAtom(uUnit->psParameters, iIndex, sType, &dMass, &dPolar,
                &dEpsilon, &dRStar, &dEpsilon14, &dRStar14, &dScreenF,
                &iElement, &iHybridization, sDesc);
    FortranWriteDouble(dGBScreenForElement(iElement));
  }
  FortranEndLine();

  // Write IPOL near the end of prmtop
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG IPOL");
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
    FortranWriteString("%FORMAT(20" LBLFORMAT_F ")");
    FortranFormat(20, LBLFORMAT_C);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
      cPTemp = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->sPertName;
      if (strlen(cPTemp) == 0) {
        cPTemp = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->sName;
      }
      if (strlen(cPTemp) > 4) {
        cPTemp += (strlen(cPTemp) - 4);
      }
      FortranWriteString(cPTemp);
    }
    FortranEndLine();

    // -36I- List of atomic symbols (atom types??????)
    FortranDebug("-36I-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG PERT_ATOM_SYMBOL");
    FortranWriteString("%FORMAT(20" LBLFORMAT_F ")");
    FortranFormat(20, LBLFORMAT_C);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
      cPTemp = sAtomPertType(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom);
      if (strlen(cPTemp) == 0) {
        cPTemp = sAtomType(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom);
      }
      if (strlen(cPTemp) > 3) {
        cPTemp += (strlen(cPTemp) - 3);
      }
      FortranWriteString(cPTemp);
    }
    FortranEndLine();

    // -36J- Value of LAMBDA for each ATOM ?????????
    // TODO: Figure out what the hell this is
    FortranDebug("-36J-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG ALMPER");
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
    FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
    FortranFormat(10, INTFORMAT_C);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
      if (bAtomPerturbed(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom)) {
        FortranWriteInt(1);
      }
      else {
        FortranWriteInt(0);
      }
    }
    FortranEndLine();

    // -36L- List of atom types - IACPER
    FortranDebug("-36L-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG PERT_ATOM_TYPE_INDEX");
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
    FortranFormat(5, DBLFORMAT_C);
    for (i = 0; i < iVarArrayElementCount(vaNBParameters); i++) {
      FortranWriteDouble(PVAI(vaNBParameters, NONBONDACt, i)->dA14);
    }
    FortranEndLine();

    // -20- Lennard jones r**6 term for all 14 interactions
    // CN214 array
    FortranDebug("-20-");
    FortranFormat(5, DBLFORMAT_C);
    for (i = 0; i < iVarArrayElementCount(vaNBParameters); i++) {
      FortranWriteDouble(PVAI(vaNBParameters, NONBONDACt, i)->dB14);
    }
    FortranEndLine();

    // -13- Force constants for Urey-Bradley
    FortranDebug("-13-");
    FortranFormat(5, DBLFORMAT_C);
    for (i = 0; i < iParmSetTotalAngleParms(uUnit->psParameters); i++) {
      ParmSetAngle(uUnit->psParameters, i, sAtom1, sAtom2, sAtom3, &dKt,
                   &dT0, &dTkub, &dRkub, sDesc);
      FortranWriteDouble(dTkub);
    }
    FortranEndLine();

    // -14- Equilibrium distances for Urey-Bradley
    FortranDebug("-14-");
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
    FortranWriteString("%FORMAT(10" INTFORMAT_F ")");
    FortranFormat(10, INTFORMAT_C);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaResidues); i++) {
      FortranWriteInt( PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->iPdbResSeq );
    }
    FortranEndLine();

    MESSAGE("Writing residue PDB ChainId\n");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG RESIDUE_CHAINID");
    FortranWriteString("%FORMAT(20" LBLFORMAT_F ")");
    FortranFormat(20, LBLFORMAT_C);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaResidues); i++) {
      FortranWriteString(PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->sChainId);
    }
    FortranEndLine();
  }

  // CMAP parameters, Mengjuei Hsieh and Yong Duan
  SaveAmberParmCMAP(uUnit, fOut);

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

