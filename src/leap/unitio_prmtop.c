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
#ifdef BINTRAJ
#  include      "netcdf.h"
#endif

#define AMBERINDEX(i)   3*(i-1)
#define INTFORMAT       "%8d"
#define DBLFORMAT       "%16.8lE"
#define LBLFORMAT       "%-4s"
#define IDFORMAT       "%-8s"
#define ELECTRONTOKCAL  18.2223

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


typedef struct {
    SAVETORSIONt *tp;
} SAVETORSIONtp;

static void get4atoms(UNIT u, SAVETORSIONt * pt, SAVEATOMt * sa4[])
{
    sa4[0] = PVAI(u->vaAtoms, SAVEATOMt, pt->iAtom1 - 1);
    sa4[1] = PVAI(u->vaAtoms, SAVEATOMt, pt->iAtom2 - 1);
    sa4[2] = PVAI(u->vaAtoms, SAVEATOMt, pt->iAtom3 - 1);
    sa4[3] = PVAI(u->vaAtoms, SAVEATOMt, pt->iAtom4 - 1);
}

static int copyatoms(int atoms[], SAVETORSIONt * sa4, SAVETORSIONt * sb4,
                     int k1, int k2)
{
    int atma[4], atmb[4], i;

    if (k1 > 0) {
        atma[0] = AMBERINDEX(sa4->iAtom1);
        atma[1] = AMBERINDEX(sa4->iAtom2);
        atma[2] = AMBERINDEX(sa4->iAtom3);
        atma[3] = AMBERINDEX(sa4->iAtom4);
    } else {
        atma[3] = AMBERINDEX(sa4->iAtom1);
        atma[2] = AMBERINDEX(sa4->iAtom2);
        atma[1] = AMBERINDEX(sa4->iAtom3);
        atma[0] = AMBERINDEX(sa4->iAtom4);
    }

    if (k2 > 0) {
        atmb[0] = AMBERINDEX(sb4->iAtom1);
        atmb[1] = AMBERINDEX(sb4->iAtom2);
        atmb[2] = AMBERINDEX(sb4->iAtom3);
        atmb[3] = AMBERINDEX(sb4->iAtom4);
    } else {
        atmb[3] = AMBERINDEX(sb4->iAtom1);
        atmb[2] = AMBERINDEX(sb4->iAtom2);
        atmb[1] = AMBERINDEX(sb4->iAtom3);
        atmb[0] = AMBERINDEX(sb4->iAtom4);
    }
// check overlap
    for (i = 0; i < 3; i++)
        if (atmb[i] != atma[i + 1]) {
            VP0("Atom indices mismatch?? in copyatoms %i %i\n", atma[i + 1],
                 atmb[i]);
            return -1;
        } else {
            atoms[i] = atma[i];
            atoms[i + 2] = atmb[i + 1];
        }
    return 4;

}

static int cmpresname1(UNIT u, SAVEATOMt * sa4, WRD reslist[], int nres)
{
    int i;
    char *sname;

    sname = PVAI(u->vaResidues, SAVERESIDUEt, sa4->iResidueIndex - 1)->sName;
    for (i = 0; i < nres; i++) {
        if (strcmp(sname, reslist[i]) == 0)
            return (i + 1);
    }

    return -1;

}

static int cmpresname4(UNIT u, SAVEATOMt * sa4[], WRD reslist[], int nres)
{
    int i, j;

    for (j = 0; j < 4; j++) {
        if ((i = cmpresname1(u, sa4[j], reslist, nres)) > 0)
            return (i + 1);
    }

    return -1;

}

static int cmp4vs4(SAVEATOMt * sa4[], WRD atm4[])
{
    int i, l1, l2;
    l1 = 0;
    l2 = 0;
    for (i = 0; i < 4; i++)
        if (strcmp(sa4[i]->sName, atm4[i]) == 0)
            l1++;
    if (l1 == 4)
        return 4;
    for (i = 0; i < 4; i++)
        if (strcmp(sa4[3 - i]->sName, atm4[i]) == 0)
            l2++;
    if (l2 == 4)
        return -4;

    return 0;

}

static int cmp_residx(SAVEATOMt * sa4[], SAVEATOMt * sb4[], int *residx)
{
    int i, l1, l2;
    int idx0, idx1;

    idx0 = sa4[0]->iResidueIndex;
    idx1 = residx[0];
    l1 = 0;
    l2 = 0;
    for (i = 0; i < 4; i++) {
        if ((sa4[i]->iResidueIndex - idx0) == (residx[i] - idx1))
            l1++;
    }
    for (i = 0; i < 4; i++) {
        if ((sb4[i]->iResidueIndex - idx0) == (residx[i + 1] - idx1))
            l2++;
    }

    if (l1 == 4 && l2 == 4)
        return 1;

    return 0;

}


//---------------------------------------------------------------------------------------------
// BondAugmentationFound: test for the existence of bond augmentations.  Return 1 if they are
//                        found, which will ordeer them to be printed.
//
// Arguments:
//   uUnit:    the tleap Unit to save
//---------------------------------------------------------------------------------------------
static int BondAugmentationFound(UNIT uUnit)
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

//---------------------------------------------------------------------------------------------
//// C4PairwiseFound: test for the existence of adding C4 interaction.  Return 1 if they are
////                        found, which will ordeer them to be printed.
////
//// Arguments:
////   uUnit:    the tleap Unit to save
////---------------------------------------------------------------------------------------------
/*
static int C4PairwiseFound(UNIT uUnit)
{
  int i, found;

  found = 0;
  for (i = 0; i < iVarArrayElementCount(uUnit->vaC4Pairwise); i++) {
    if (PVAI(uUnit->vaC4Pairwise, SAVEC4Pairwiset, i)->daC4Pairwise > 0.1) {
      found = 1;
      break;
    }
  }

  return found;
}
*/

static void SaveAmberParmCMAP(UNIT uUnit, FILE * fOut)
{
    //
    // CMAP parameters, Mengjuei Hsieh and Yong Duan
    //
    int i, j, k, l;
    int mapid, maptypes;
//    int mapcount;
    int *mapflag, *mapidx;
    int iNumDIH;
    int nprospect, ires;
    SAVEATOMt *sa4[4], *sb4[4];
    SAVETORSIONt *stPTorsion2, *stPTorsion;
    SAVETORSIONtp *stdpt0;
    STRING sTmp;
    PHIPSI *phipsi;
    CMAPLST *cmaplstt;
    int k1, k2;
    CMAP *cmap;
    int maxmap;

    if (!GDefaults.iCMAP)
        return;
    if (mapnum <= 0)
        return;

    for (cmaplstt = cmaplst; cmaplstt->next != NULL; cmaplstt = cmaplstt->next) {
        cmap = cmaplstt->cmap;
    }
    iNumDIH = iVarArrayElementCount(uUnit->vaTorsions);
    if (iNumDIH > 0) {
        stPTorsion = PVAI(uUnit->vaTorsions, SAVETORSIONt, 0);
    } else {
        return;
    }

//    mapcount = 0;
    MALLOC(mapflag, int *, sizeof(int) * (mapnum + 1));
    MALLOC(mapidx, int *, sizeof(int) * (mapnum + 1));

    i = 0;
    for (i = 0; i <= mapnum; i++) {
        mapflag[i] = 0;
        mapidx[i] = 0;
    }

    maxmap = iNumDIH;
    MALLOC(phipsi, PHIPSI *, sizeof(PHIPSI) * (maxmap));
    MALLOC(stdpt0, SAVETORSIONtp *, sizeof(SAVETORSIONtp) * (iNumDIH));
    nprospect = 0;
//    mapcount = 0;
    // Loop over dihedral list
// pre-filter removes the irrelevant torsions first ...
    for (i = 0; i < iNumDIH; i++, stPTorsion++) {

//            if (stPTorsion->bCalc14 == 0) continue; //cycle

        get4atoms(uUnit, stPTorsion, sa4);

        for (cmaplstt = cmaplst; cmaplstt->next != NULL;
             cmaplstt = cmaplstt->next) {
            cmap = cmaplstt->cmap;
            if (cmpresname4(uUnit, sa4, cmap->reslist, cmap->nres) > 0) // potential match
            {
                if (abs(cmp4vs4(sa4, (WRD *) (&cmap->atmname[0]))) == 4 || abs(cmp4vs4(sa4, (WRD *) (&cmap->atmname[1]))) == 4) {       // keep in the list
                    stdpt0[nprospect].tp = stPTorsion;
                    nprospect++;
                    break;
                }
            } else if (cmap->termmap > 0) {
                if (cmpresname4(uUnit, sa4, cmap->creslist, cmap->nres) > 0) {
                    if (abs(cmp4vs4(sa4, (WRD *) (&cmap->catmname[0]))) == 4
                        || abs(cmp4vs4(sa4, (WRD *) (&cmap->catmname[1]))) ==
                        4) {
                        stdpt0[nprospect].tp = stPTorsion;
                        nprospect++;
                        break;
                    }
                } else if (cmpresname4(uUnit, sa4, cmap->nreslist, cmap->nres) >
                           0)
                    if (abs(cmp4vs4(sa4, (WRD *) (&cmap->natmname[0]))) == 4
                        || abs(cmp4vs4(sa4, (WRD *) (&cmap->natmname[1]))) ==
                        4) {
                        stdpt0[nprospect].tp = stPTorsion;
                        nprospect++;
                        break;
                    }
            }
        }
    }
//
    k = 0;
    int k1c, k2c, k1n, k2n;
    for (i = 0; i < nprospect; i++) {
        stPTorsion = stdpt0[i].tp;
        get4atoms(uUnit, stPTorsion, sa4);
        mapid = 0;
        for (cmaplstt = cmaplst; cmaplstt->next != NULL;
             cmaplstt = cmaplstt->next) {
            cmap = cmaplstt->cmap;
            mapid++;
            k1 = cmp4vs4(sa4, (WRD *) (&cmap->atmname[0]));     // check atmnames 0..3
            k1c = cmp4vs4(sa4, (WRD *) (&cmap->catmname[0]));   // check atmnames 0..3
            k1n = cmp4vs4(sa4, (WRD *) (&cmap->natmname[0]));   // check atmnames 0..3

            if (cmap->termmap == 0) {
                k1c = 0;
                k1n = 0;
            }

            if (abs(k1) == 4 || abs(k1c) == 4 || abs(k1n) == 4) {       // match the atom names
                int iresc=0, iresn=0;
                // "0" in residx[] marks "present residue"
                for (l = 0; l < 5 && cmap->residx[l] != 0; l++);
                if (l < 4) {
                    ires =
                        cmpresname1(uUnit, sa4[l], cmap->reslist, cmap->nres);
                    iresc =
                        cmpresname1(uUnit, sa4[l], cmap->creslist, cmap->nres);
                    iresn =
                        cmpresname1(uUnit, sa4[l], cmap->nreslist, cmap->nres);

                    if (abs(k1) != 4)
                        ires = 0;
                    if (abs(k1c) != 4)
                        iresc = 0;
                    if (abs(k1n) != 4)
                        iresn = 0;
                }

                if (ires > 0 || iresc > 0 || iresn > 0 ||       // found on the reslist of the cmap
                    l == 4) {   // or the "present residue" pointer is the last one.
                    for (j = 0; j < nprospect; j++) {
                        stPTorsion2 = stdpt0[j].tp;
                        get4atoms(uUnit, stPTorsion2, sb4);
                        if (l == 4) {
                            ires =
                                cmpresname1(uUnit, sb4[3], cmap->reslist,
                                            cmap->nres);
                            iresc =
                                cmpresname1(uUnit, sb4[3], cmap->creslist,
                                            cmap->nres);
                            iresn =
                                cmpresname1(uUnit, sb4[3], cmap->nreslist,
                                            cmap->nres);
                        }

                        if (cmap->termmap == 0) {
                            iresc = 0;
                            iresn = 0;
                        }

                        k2 = cmp4vs4(sb4, (WRD *) (&cmap->atmname[1])); // check atmnames 1..4
                        k2c = cmp4vs4(sb4, (WRD *) (&cmap->catmname[1]));       // check atmnames 1..4
                        k2n = cmp4vs4(sb4, (WRD *) (&cmap->natmname[1]));       // check atmnames 1..4

                        if (abs(k2) == 4 && abs(k1) == 4 && ires > 0)
                            // check residue indicies, 0: present, -1: residue before, +1: residue after
                            if (cmp_residx(sa4, sb4, cmap->residx))     // we found two torsions and the cmap
                            {
                                copyatoms(phipsi[k].atoms, stPTorsion,
                                          stPTorsion2, k1, k2);
                                phipsi[k].mapid = mapid;
                                // mark as used
                                mapflag[mapid - 1] = 1;
                                k++;
                            }

                        if (abs(k2c) == 4 && abs(k1c) == 4 && iresc > 0)
                            if (cmp_residx(sa4, sb4, cmap->cresidx)) {
                                copyatoms(phipsi[k].atoms, stPTorsion,
                                          stPTorsion2, k1c, k2c);
                                phipsi[k].mapid = mapid;
                                // mark as used
                                mapflag[mapid - 1] = 1;
                                k++;
                            }

                        if (abs(k2n) == 4 && abs(k1n) == 4 && iresn > 0)
                            if (cmp_residx(sa4, sb4, cmap->nresidx)) {
                                copyatoms(phipsi[k].atoms, stPTorsion,
                                          stPTorsion2, k1n, k2n);
                                phipsi[k].mapid = mapid;
                                // mark as used
                                mapflag[mapid - 1] = 1;
                                k++;
                            }
                    }
                }
            }
        }
    }

    int mk = 0;
    for (i = 0; i < k; i++) {
        int flag = 1;

// Remove the redundant torsions.

        int l;
        if (i > 0)
            for (l = 0; l < i; l++) {
                if (phipsi[i].atoms[0] == phipsi[l].atoms[0]
                    && phipsi[i].atoms[1] == phipsi[l].atoms[1]
                    && phipsi[i].atoms[2] == phipsi[l].atoms[2]
                    && phipsi[i].atoms[3] == phipsi[l].atoms[3]
                    && phipsi[i].atoms[4] == phipsi[l].atoms[4]
                    && phipsi[i].mapid == phipsi[l].mapid)
                    flag = -1;
            }
//
        for (j = 0; j < 5; j++)
            if (phipsi[i].atoms[j] < 0)
                flag = -1;
        if (phipsi[i].mapid <= 0)
            flag = -1;
        if (flag > 0)
            mk++;
    }

    //wmap
    maptypes = 0;

    for (i = 0; i < mapnum; i++) {
        if (mapflag[i] == 1) {
            mapidx[i] = maptypes;
            maptypes++;
        }
    }
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG CMAP_COUNT");
    //FortranWriteString("%COMMENT");
//        for (cmntt=cmnt0; cmntt->next != NULL; cmntt=cmntt->next){
//            sTmp[0]='\0';
//            strcat(sTmp,"%COMMENT");
//            strcat(sTmp,cmntt->record);
//            FortranWriteString(sTmp);
    //fprintf(fpout,"%%COMMENT%s",cmntt->record);
//        }
    FortranWriteString("%FORMAT(2I8)");
//        sprintf(sTmp,"%8d%8d",mapcount,mapnum);
//        sprintf(sTmp,"%8d%8d",k,maptypes);
    sprintf(sTmp, "%8d%8d", mk, maptypes);
    FortranWriteString(sTmp);
    //FortranEndLine();

    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG CMAP_RESOLUTION");
    FortranWriteString("%FORMAT(20I4)");
    FortranFormat(20, "%4d");
    mapid = 0;
    for (cmaplstt = cmaplst; cmaplstt->next != NULL; cmaplstt = cmaplstt->next) {
        cmap = cmaplstt->cmap;
        if (mapflag[mapid])
            FortranWriteInt(cmap->resolution);
        mapid++;
    }
    FortranEndLine();

    mapid = 0;
    for (cmaplstt = cmaplst; cmaplstt->next != NULL; cmaplstt = cmaplstt->next) {
        cmap = cmaplstt->cmap;
        int msize;
        if (mapflag[mapid] == 1) {
            sprintf(sTmp, "%%FLAG CMAP_PARAMETER_%02d", mapidx[mapid] + 1);
            FortranFormat(1, "%-80s");
            FortranWriteString(sTmp);
            sprintf(sTmp, "%%COMMENT  %s", cmap->title);
            FortranWriteString(sTmp);
            FortranWriteString("%FORMAT(8F9.5)");
            FortranFormat(8, "%9.5lf");
            msize = cmap->resolution * cmap->resolution;
            for (j = 0; j < msize; j += 8) {
                int l;
                for (l = j; l < msize && l < j + 8; l++) {
                    FortranWriteDouble(cmap->map[l]);
                }
            }
            FortranEndLine();
        }
        mapid++;
    }

    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG CMAP_INDEX");
    FortranWriteString("%FORMAT(6I8)");
    FortranFormat(6, "%8d");
    for (i = 0; i < k; i++) {
        int flag = 1;

// Remove the redundant torsions.

        int l;
        if (i > 0)
            for (l = 0; l < i; l++) {
                if (phipsi[i].atoms[0] == phipsi[l].atoms[0]
                    && phipsi[i].atoms[1] == phipsi[l].atoms[1]
                    && phipsi[i].atoms[2] == phipsi[l].atoms[2]
                    && phipsi[i].atoms[3] == phipsi[l].atoms[3]
                    && phipsi[i].atoms[4] == phipsi[l].atoms[4]
                    && phipsi[i].mapid == phipsi[l].mapid)
                    flag = -1;
            }
//
        for (j = 0; j < 5; j++)
            if (phipsi[i].atoms[j] < 0)
                flag = -1;
        if (phipsi[i].mapid <= 0)
            flag = -1;
        if (flag > 0) {
            for (j = 0; j < 5; j++)
                FortranWriteInt(phipsi[i].atoms[j] / 3 + 1);
            FortranWriteInt(mapidx[phipsi[i].mapid - 1] + 1);
        }
    }
    FortranEndLine();
    FREE(mapflag);
    FREE(mapidx);
    FREE(phipsi);
    FREE(stdpt0);
}



/*
 *        zUnitIOFindAndCountMolecules
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
static void
zUnitIOFindAndCountMolecules(UNIT uUnit)
{
    SAVERESIDUEt *srPRes;
    int i, iResidues, iAtom, iCount;
    LOOP lSpanning;
    BOOL bSeenFirstSolvent = FALSE;
    ATOM aAtom;

    if (uUnit->vaAtomsPerMolecule) VarArrayDestroy(&(uUnit->vaAtomsPerMolecule));
    uUnit->vaAtomsPerMolecule = vaVarArrayCreate(sizeof(int));

    /* Clear the ATOMTOUCHED flag on all the ATOMs */

    ContainerResetAllAtomsFlags((CONTAINER) uUnit, ATOMTOUCHED);

    /* Get the first RESIDUE */

    srPRes = PVAI(uUnit->vaResidues, SAVERESIDUEt, 0);
    iResidues = iVarArrayElementCount(uUnit->vaResidues);

    /* Loop over all RESIDUES */

    for (i = 0; i < iResidues; i++, srPRes++) {

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

        iCount = 0;
        lSpanning = lLoop((OBJEKT) aAtom, SPANNINGTREE);
        FOREACH(aAtom, ATOM, lSpanning) {
            AtomSetFlags(aAtom, ATOMTOUCHED);
            iCount++;
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
 *      CheckTypeNBEdit
 *
 *        Author:       David S. Cerutti (2013)
 *
 *      Check to see whether this atom type (which is imminently to be
 *      compacted in the nonbonded parameters array) is mentioned in any
 *      nonbonded pair potential adjustments.  If so, then it has to remain
 *      its own type.
 */
static int
CheckTypeNBEdit(typeStr sType, VARARRAY vaPNBEdits)
{
  int i;

  for (i = 0; i < vaPNBEdits->count; i++) {
    if (strcmp(sType, PVAI(vaPNBEdits, NBEDITt, i)->sType1) == 0 ||
	strcmp(sType, PVAI(vaPNBEdits, NBEDITt, i)->sType2) == 0) {
      return 1;
    }
  }

  return 0;
}

/*
 *      zCheckAgainstNBEdits
 *
 *        Author:       David S. Cerutti (2013)
 *
 *      Check to see whether an atom type, when paired with any other atom
 *      type, has any Lennard-Jones sigma and epsilon pairs that do not
 *      conform to the standard combining rules.
 */
static void
zCheckAgainstNBEdits(VARARRAY vaPNBEdits, typeStr tI, typeStr tJ,
		    double *dA, double *dC)
{
  int i;

  for (i = 0; i < vaPNBEdits->count; i++) {
    if ((strcmp(PVAI(vaPNBEdits, NBEDITt, i)->sType1, tI) == 0 &&
	 strcmp(PVAI(vaPNBEdits, NBEDITt, i)->sType2, tJ) == 0) ||
	(strcmp(PVAI(vaPNBEdits, NBEDITt, i)->sType1, tJ) == 0 &&
	 strcmp(PVAI(vaPNBEdits, NBEDITt, i)->sType2, tI) == 0)) {

      /* This is an edited pair interaction */
      MathOpConvertNonBondToAC(PVAI(vaPNBEdits, NBEDITt, i)->dEI,
			       PVAI(vaPNBEdits, NBEDITt, i)->dRI,
			       PVAI(vaPNBEdits, NBEDITt, i)->dEJ,
			       PVAI(vaPNBEdits, NBEDITt, i)->dRJ, dA, dC);
      break;
    }
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
static void
zUnitIOBuildNonBondArrays(UNIT uUnit, VARARRAY * vaPNBIndexMatrix,
                          VARARRAY * vaPNBParameters,
                          VARARRAY * vaPNBIndex, VARARRAY * vaPNonBonds,
                          char sA[8][16], char sB[8][16], double daC4Type[16], int iC4count ) //NewT
{
    VARARRAY vaNBIndex, vaNonBonds, vaPNBEdits;
    int i, j, iNonBonds, iNBIndices, iTemp, iI, iJ, iX, iY; //, iC4, jC4;
    int iIndex, iHBondIndex, iNBIndex, iElement, iHybridization;
    double dMass, dPolar, dE, dR, dE14, dR14, dA, dC, dEI, dRI, dEJ, dRJ;
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
                        for (nbPTemp = nbPCur; nbPTemp < nbPLast;
                             nbPTemp++) *nbPTemp = *(nbPTemp + 1);
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
            /* -ve to signify that the index is into the HBOND tables */

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
            MathOpConvertNonBondToAC(dEI, dRI, dEJ, dRJ, &dA, &dC);
	    zCheckAgainstNBEdits(vaPNBEdits,
				PVAI(vaNonBonds, NONBONDt, i)->sType,
				PVAI(vaNonBonds, NONBONDt, j)->sType,
				&dA, &dC);
            PVAI(*vaPNBParameters, NONBONDACt, iIndex)->dA = dA;
            PVAI(*vaPNBParameters, NONBONDACt, iIndex)->dC = dC;
            // Initialize all C4 to be zero NewT
            PVAI(*vaPNBParameters, NONBONDACt, iIndex)->d4 = 0.0;
            if (GDefaults.iCharmm) {
                dEI = PVAI(vaNonBonds, NONBONDt, i)->dE14;
                dRI = PVAI(vaNonBonds, NONBONDt, i)->dR14;
                dEJ = PVAI(vaNonBonds, NONBONDt, j)->dE14;
                dRJ = PVAI(vaNonBonds, NONBONDt, j)->dR14;
                MathOpConvertNonBondToAC(dEI, dRI, dEJ, dRJ, &dA, &dC);
		zCheckAgainstNBEdits(vaPNBEdits,
				    PVAI(vaNonBonds, NONBONDt, i)->sType,
				    PVAI(vaNonBonds, NONBONDt, j)->sType,
				    &dA, &dC);
                PVAI(*vaPNBParameters, NONBONDACt, iIndex)->dA14 = dA;
                PVAI(*vaPNBParameters, NONBONDACt, iIndex)->dC14 = dC;
            }
        }
    }
    // Only update C4 that matches both atom type names NewT
    for (int jj = 0; jj < iC4count; jj ++)
    {
	int rawA = iParmSetFindAtom(uUnit->psParameters, sA[jj]);
        int rawB = iParmSetFindAtom(uUnit->psParameters, sB[jj]);

        if (rawA == PARM_NOT_FOUND || rawB == PARM_NOT_FOUND) {
            VP0("C4Type: could not find type %s or %s in ParmSet\n", sA[jj], sB[jj]);
            continue;
        }

        // Map raw atom parm index -> rolled-up unique NB index (0-based)
        iI = *PVAI(vaNBIndex, int, rawA);
        iJ = *PVAI(vaNBIndex, int, rawB);

        iX = MIN(iI, iJ);
        iY = MAX(iI, iJ);

        // 0-based triangular index into vaPNBParameters
        int idx = iY * (iY + 1) / 2 + iX;

        PVAI(*vaPNBParameters, NONBONDACt, idx)->d4 = daC4Type[jj];
        //iC4 = -1;
        //jC4 = -1;
        //for (int ii = 0; ii < iVarArrayElementCount(uUnit->vaAtoms); ii++) {
        //    //VP0("Atom type is: %s\n", sAtomType(PVAI(uUnit->vaAtoms, SAVEATOMt, ii)->aAtom));
        //    //VP0("Target is: %s\n", sA[jj]);
        //    if (!strcmp(sAtomType(PVAI(uUnit->vaAtoms, SAVEATOMt, ii)->aAtom), sA[jj])) {
        //        iC4 = PVAI(uUnit->vaAtoms, SAVEATOMt, ii)->iTypeIndex;
        //        if (iC4 > iVarArrayElementCount(vaNonBonds)) {
        //                iC4 = iC4 - iVarArrayElementCount(vaNonBonds);
        //        }
	//	//VP0("iC4 matches!!! value is %i\n", iC4);
        //    }
        //    if (!strcmp(sAtomType(PVAI(uUnit->vaAtoms, SAVEATOMt, ii)->aAtom), sB[jj])) {
        //        jC4 = PVAI(uUnit->vaAtoms, SAVEATOMt, ii)->iTypeIndex;
        //        if (jC4 > iVarArrayElementCount(vaNonBonds)) {
        //                jC4 = jC4 - iVarArrayElementCount(vaNonBonds);
        //        }
	//	//VP0("jC4 matches!!! value is %i\n", jC4);
        //    }
        //}
        //if (iC4 != -1 && jC4 != -1) {
        //VP0("daC4Type value is %f\n", daC4Type[jj]);
        //PVAI(*vaPNBParameters, NONBONDACt, MAX(iC4,jC4)*(MAX(iC4,jC4)-1)/2+MIN(iC4,jC4)-1)->d4 = daC4Type[jj];
        //}
    }

    /* Return the arrays that refer to the atoms */
    *vaPNonBonds = vaNonBonds;
    *vaPNBIndex = vaNBIndex;
    /* The other two arrays, vaPNBIndexMatrix and vaPNBParameters, */
    /* have already been assigned their return values. */

}


/*
 *      zUnitIOSaveAmberNetcdf
 *
 *      Author: Robin Betz (2013)
 *
 *      Writes the coordinates in UNIT to a coordinate file in netCdf format.
 *      This is written based from NetcdfFile.cpp in cpptraj and AmberNetcdf.F90
 *      in pmemd/sander
 *
 *      Arguments:
 *          uUnit     - UNIT to save
 *          filename  - name of netcdf file to write
 */
static void
zUnitIOSaveAmberNetcdf(UNIT uUnit, char *filename)
{
#   ifndef BINTRAJ
    VPFATALEXIT("Built without NETCDF support. Rebuild with -DBINTRAJ\n");
#   else

    int ncid;                   // netcdf file handle
    int did_spatial, did_atom;       // dimension IDs
    int vid_spatial, vid_coord; // variable IDs
    int did_cell_spatial, did_cell_angular, did_label;
    int vid_cell_spatial, vid_cell_angular, vid_cell_length, vid_cell_angle,
        vid_time;
    int dimensionID[NC_MAX_VAR_DIMS];

    // Get number of atoms T_T
    int iAtomCount = iVarArrayElementCount(uUnit->vaAtoms);
    printf("There are %i atoms\n", iAtomCount);

    // Create the file
    if (nc_create(filename, NC_64BIT_OFFSET, &ncid) != NC_NOERR) {
        VPFATALEXIT("%s: Error creating file\n", filename);
    }
    // Spatial dimension and variable
    if (nc_def_dim(ncid, "spatial", 3, &did_spatial) != NC_NOERR) {
        VPFATALEXIT("%s: Error defining spatial dimension\n", filename);
    }
    dimensionID[0] = did_spatial;
    if (nc_def_var(ncid, "spatial", NC_CHAR, 1, dimensionID, &vid_spatial) !=
        NC_NOERR) {
        VPFATALEXIT("%s: Error defining spatial variable\n", filename);
    }
    // Atom dimension
    if (nc_def_dim(ncid, "atom", iAtomCount, &did_atom)
        != NC_NOERR) {
        VPFATALEXIT("%s: Error defining atom dimension\n", filename);
    }
    // Time dimension and variable
    if (nc_def_var(ncid, "time", NC_DOUBLE, 0, dimensionID, &vid_time) !=
        NC_NOERR) {
        VPFATALEXIT("%s: Error defining time variable\n", filename);
    }
    if (nc_put_att_text(ncid, vid_time, "units", 10, "picosecond") != NC_NOERR) {
        VPFATALEXIT("%s: Error setting time units to picosecond\n", filename);
    }
    dimensionID[0] = did_atom;
    dimensionID[1] = did_spatial;
    // Coord variable and attribute text
    if (nc_def_var(ncid, "coordinates", NC_DOUBLE, 2, dimensionID, &vid_coord)
        != NC_NOERR) {
        VPFATALEXIT("%s: Error defining coordinate variable\n", filename);
    }
    if (nc_put_att_text(ncid, vid_coord, "units", 8, "angstrom") != NC_NOERR) {
        VPFATALEXIT("%s: Error setting coordinate units to angstrom\n", filename);
    }
    // Define box if it exists
    if (bUnitUseBox(uUnit) == TRUE) {
        printf("Using the unit box\n");
        // Cell spatial
        if (nc_def_dim(ncid, "cell_spatial", 3, &did_cell_spatial) != NC_NOERR) {
            VPFATALEXIT("%s: Error defining cell spatial dimension\n", filename);
        }
        dimensionID[0] = did_cell_spatial;
        if (nc_def_var
            (ncid, "cell_spatial", NC_CHAR, 1, dimensionID, &vid_cell_spatial)
            != NC_NOERR) {
            VPFATALEXIT("%s: Error defining cell spatial variable\n", filename);
        }
        // Cell angular
        if (nc_def_dim(ncid, "label", 5, &did_label) != NC_NOERR) {
            VPFATALEXIT("%s: Error defining label dimension\n", filename);
        }
        if (nc_def_dim(ncid, "cell_angular", 3, &did_cell_angular) != NC_NOERR) {
            VPFATALEXIT("%s: Error defining cell angular dimension\n", filename);
        }
        dimensionID[0] = did_cell_angular;
        dimensionID[1] = did_label;
        if (nc_def_var
            (ncid, "cell_angular", NC_CHAR, 2, dimensionID, &vid_cell_angular)
            != NC_NOERR) {
            VPFATALEXIT("%s: Error defining cell angular variable\n", filename);
        }
        // Box dimensions
        dimensionID[0] = did_cell_spatial;
        if (nc_def_var
            (ncid, "cell_lengths", NC_DOUBLE, 1, dimensionID, &vid_cell_length)
            != NC_NOERR) {
            VPFATALEXIT("%s: Error defining cell lengths\n", filename);
        }
        if (nc_put_att_text(ncid, vid_cell_length, "units", 8, "angstrom") !=
            NC_NOERR) {
            VPFATALEXIT("%s: Error setting cell length units to angstrom\n",
                 filename);
        }
        dimensionID[0] = did_cell_angular;
        if (nc_def_var
            (ncid, "cell_angles", NC_DOUBLE, 1, dimensionID, &vid_cell_angle)
            != NC_NOERR) {
            VPFATALEXIT("%s: Error defining cell angles variable\n", filename);
        }
        if (nc_put_att_text(ncid, vid_cell_angle, "units", 6, "degree") !=
            NC_NOERR) {
            VPFATALEXIT("%s: Error setting cell angle units to degree\n", filename);
        }
    }
    // Conventions and file attributes
    if (nc_put_att_text
        (ncid, NC_GLOBAL, "title", strlen(sContainerName(uUnit)),
         sContainerName(uUnit)) != NC_NOERR) {
        VPFATALEXIT("%s: Error writing title\n", filename);
    }
    if (nc_put_att_text(ncid, NC_GLOBAL, "application", 5, "AMBER") != NC_NOERR) {
        VPFATALEXIT("%s: Error writing application string\n", filename);
    }
    if (nc_put_att_text(ncid, NC_GLOBAL, "program", 4, "leap") != NC_NOERR) {
        VPFATALEXIT("%s: Error writing program string\n", filename);
    }
    if (nc_put_att_text(ncid, NC_GLOBAL, "programVersion", 3, "1.0") !=
        NC_NOERR) {
        VPFATALEXIT("%s: Error writing program version string\n", filename);
    }
    if (nc_put_att_text(ncid, NC_GLOBAL, "Conventions", 12, "AMBERRESTART") !=
        NC_NOERR) {
        VPFATALEXIT("%s: Error writing conventions\n", filename);
    }
    if (nc_put_att_text(ncid, NC_GLOBAL, "ConventionVersion", 3, "1.0") !=
        NC_NOERR) {
        VPFATALEXIT("%s: Error writing conventions version\n", filename);
    }
    // Set fill mode and end definitions
    if (nc_set_fill(ncid, NC_NOFILL, dimensionID) != NC_NOERR) {
        VPFATALEXIT("%s: Error setting fill mode\n", filename);
    }
    if (nc_enddef(ncid) != NC_NOERR) {
        VPFATALEXIT("%s: NetCDF error on ending definitions\n", filename);
    }
    // Spatial dimension labels
    size_t start[2]={0}, count[2]={0};
    count[0] = 3; // 1 x abc
    const char *xyz = "xyz";
    if (nc_put_vara_text(ncid, vid_spatial, start, count, xyz) != NC_NOERR) {
        VPFATALEXIT("%s: Error writing spatial labels\n", filename);
    }
    if (bUnitUseBox(uUnit) == TRUE) {
        const char *abc = "abc";
        if (nc_put_vara_text(ncid, vid_cell_spatial, start, count, abc) !=
            NC_NOERR) {
            VPFATALEXIT("%s: Error writing cell spatial labels\n", filename);
        }
        count[1] = 5; // 3 x 5
        const char *abg = "alpha" "beta " "gamma";
        if (nc_put_vara_text(ncid, vid_cell_angular, start, count, abg) !=
            NC_NOERR) {
            VPFATALEXIT("%s: Error writing cell angular labels\n", filename);
        }
    }
// Create data array
    double *data;
    MALLOC(data, double *, 3 * iAtomCount * sizeof(double));
    int counter = 0;
    VECTOR vPos;

// Write time = 0
    double time = 0.0;
    if (nc_put_var_double(ncid, vid_time, &time) != NC_NOERR) {
        VPFATALEXIT("%s: Error writing start time\n", filename);
    }
// Calculate box shift if there's a box and we're centering
    double dX, dY, dZ;
    double dX2, dY2, dZ2;
    if (bUnitUseBox(uUnit) == TRUE && GDefaults.nocenter == 0) {
        UnitGetBox(uUnit, &dX, &dY, &dZ);
        dX2 = dX * 0.5;
        dY2 = dY * 0.5;
        dZ2 = dZ * 0.5;
    } else {
        dX2 = dY2 = dZ2 = 0.0;
    }

    // Fill the coordinate array with calculated shift (none if nobox)
    int i;
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
        vPos = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->vPos;
        data[counter] = dVX(&vPos) + dX2;
        ++counter;
        data[counter] = dVY(&vPos) + dY2;
        ++counter;
        data[counter] = dVZ(&vPos) + dZ2;
        ++counter;
    }

    // Write the coordinate array to the file
    start[0] = start[1] = 0;
    count[0] = iAtomCount;
    count[1] = 3;
    if (nc_put_vara_double(ncid, vid_coord, start, count, data) != NC_NOERR) {
        FREE(data);
        VPFATALEXIT("%s: Error writing coordinate data\n", filename);
    }
    // Write box lengths and angles if necessary
    if (bUnitUseBox(uUnit) == TRUE) {
        count[0] = 3;
        count[1] = 0;
        double lengths[3];
        lengths[0] = dX;
        lengths[1] = dY;
        lengths[2] = dZ;
        if (nc_put_vara_double(ncid, vid_cell_length, start, count, lengths) !=
            NC_NOERR) {
            FREE(data);
            VPFATALEXIT("%s: Error writing cell lengths\n", filename);
        }
        lengths[0] = lengths[1] = lengths[2] = dUnitBeta(uUnit) / DEGTORAD;
        if (nc_put_vara_double(ncid, vid_cell_angle, start, count, lengths) !=
            NC_NOERR) {
            FREE(data);
            VPFATALEXIT("%s: Error writing cell angles\n", filename);
        }
    }
    // Close the file
    if (nc_close(ncid) != NC_NOERR) {
        FREE(data);
        VPFATALEXIT("%s: Error closing file\n", filename);
    }
    FREE(data);
    printf("Successfully saved NetCDF inpcrd file \"%s\"\n", filename);
#   endif
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
void zUnitIOSaveAmberParmFormat(UNIT uUnit, FILE * fOut, char *crdName,
                                BOOL bPolar, BOOL bPert, BOOL bNetcdf,
                                char sA[8][16], char sB[8][16], double daC4Type[16], int iC4count ) //NewT
{
  int i, iMax, iIndex;
  LOOP lTemp, lSpan;
  ATOM aAtom, aAtomA, aA, aB, aC, aD;
  int iCount, iBondWith, iBondWithout, iAngleWith, iAngleWithout,
    iTorsionWith, iTorsionWithout, iNumExtra, iResidueIndex, nBondTypes;
  int iCountPerturbed, iCountBondPerturbed, iCountBondBoundary;
  int iCountAnglePerturbed, iCountAngleBoundary;
  int iCountTorsionPerturbed, iCountTorsionBoundary;
  VARARRAY vaExcludedAtoms, vaExcludedCount, vaNBIndexMatrix, vaNBParameters,
           vaNBIndex, vaNonBonds;
  SAVEBONDt *sbPBond;
  SAVEC4Pairwiset *scPC4Pairwise; //New
  SAVEANGLEt *saPAngle;
  SAVEATOMt *saPAtom;
  SAVETORSIONt *stPTorsion;
  SAVERESTRAINTt *srPRestraint;
  double dMass, dPolar, dR, dKb, dR0, dKpull, dRpull0, dKpress, dRpress0, dKt, dT0, dTkub,
         dRkub, dKp, dP0, dC, dD, dTemp, dScEE, dScNB, dScreenF;//, daC4Pairwise; //New dSceeScaleFactor not included
  STRING sAtom1, sAtom2, sAtom3, sAtom4, sType1, sType2;
  int iN, iAtoms, iMaxAtoms, iTemp, iAtom, iCalc14, iProper;
  int iElement, iHybridization, iStart;
  RESIDUE rRes;
  BOOL bFoundSome;
  VECTOR vPos;
  char *cPTemp=NULL;
  double dX, dY, dZ, dEpsilon, dRStar, dEpsilon14, dRStar14;
  STRING sDesc, sType;
  IX_REC eResEnt = {NULL};
  IX_DESC iResIx;
  char sVersionHeader[81];
  time_t tp;
  double dGBrad, dGBscreen;

  // Open the coordinate file
  FILE *fCrd = FOPENCOMPLAIN(crdName, "w");
  if (fCrd == NULL) {
    VP0("%s: Could not open file: %s\n", crdName);
  }

  // NOT allowed to save IPOL=0 for polarizable parm
  if (bPolar && GDefaults.iIPOL <= 0) {
    VP0("  Conflict: polarizable prmtop can not have IPOL <= 0.\n");
    VP0("  Please change IPOL in frcmod/parmxx.dat or set default IPOL.\n");
    return;
  }
  else if (!bPolar && GDefaults.iIPOL > 0) {
    VP0("  Conflict: non-polarizable prmtop can not have IPOL > 0.\n");
    VP0("  Please change IPOL in frcmod/parmxx.dat or set default IPOL.\n");
    return;
  }

  /*********** Build the excluded atom list *********/
  MESSAGE("Building the excluded atom list\n");
  vaExcludedCount = vaVarArrayCreate(sizeof(int));
  vaExcludedAtoms = vaVarArrayCreate(sizeof(int));

  iCountPerturbed = 0;
  // VP0("what iosaveparm returns %d\n", iVarArrayElementCount(uUnit->vaAtoms)); //NewTdebug
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
    aAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom;
    if (bAtomFlagsSet(aAtom, ATOMPERTURB)) {
      iCountPerturbed++;
    }
    lSpan = lLoop((OBJEKT) aAtom, SPANNINGTREE);
    iCount = 0;
    bFoundSome = FALSE;
    iStart = iVarArrayElementCount(vaExcludedAtoms);
    while ((aA = (ATOM) oNext(&lSpan))) {
      if (aA == aAtom) {
        continue;
      }

      // If the atom is more than three away from the first
      // atom then it is not in the excluded atom list

      if (iAtomBackCount(aA) >= 4) {
        break;
      }
      if (iContainerTempInt(aA) > iContainerTempInt(aAtom)) {
        VarArrayAdd(vaExcludedAtoms, (GENP) & iContainerTempInt(aA));
        bFoundSome = TRUE;
        iCount++;
      }
    }
    if (!bFoundSome) {
      iAtoms = 0;
      VarArrayAdd(vaExcludedAtoms, (GENP) & iAtoms);
      iCount++;
    }
    else {

      // Sort the part of the VARARRAY just added so that the
      // excluded ATOMs are in ascending order by index
      SortByInteger((GENP) PVAI(vaExcludedAtoms, int, iStart), iCount,
                    sizeof(int), (GENP) PVAI(vaExcludedAtoms, int, iStart), TRUE);
    }
    VarArrayAdd(vaExcludedCount, (GENP) & iCount);
  }

  /************* mark atom tree classification ********/
  // Mark main chain atoms where possible, noting the
  // number of atoms in the largest residue. keep
  // track of residues which can't be marked.
  // NOTE: Atom double dTemp is used to hold a character
  // tree flag, because
  //VP0("Not Marking per-residue atom chain types.\n");
  //iMaxAtoms = 0;

  create_index(&iResIx, IX_DUPKEYREC, IX_LEN_CSTRING);

  VP0("Marking per-residue atom chain types.\n");
  iMaxAtoms = 0;
  lTemp = lLoop((OBJEKT) uUnit, RESIDUES);
  while ((rRes = (RESIDUE) oNext(&lTemp))) {
    int iAtoms = iMarkMainChainAtoms(rRes, 0);
    if (iAtoms > 0) {
      MarkSideChains(rRes);
    }
    if (iAtoms < 0) {
      iAtoms = -iAtoms;

      // Couldn't mark main chains
      strcpy(eResEnt.key, rRes->cHeader.sName);
      if (add_key(&eResEnt, &iResIx) != IX_OK) {
        DFATAL("add_key() residue chain\n");
      }
    }
    if (iAtoms > iMaxAtoms) {
      iMaxAtoms = iAtoms;
    }
  }

  /***************** Print warnings *****************/
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
  if (!i) {
    VP0("  )\n");
  }
  destroy_index(&iResIx);

  /******* Build the NON-BOND arrays that AMBER needs ******/
  zUnitIOBuildNonBondArrays(uUnit, &vaNBIndexMatrix, &vaNBParameters,
                            &vaNBIndex, &vaNonBonds, sA, sB, daC4Type, iC4count ); //NewT
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
  FortranWriteString("%FORMAT(20a4)");
  FortranWriteString(sContainerName(uUnit));

  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG POINTERS");
  FortranWriteString("%FORMAT(10I8)");

  // -2- Save control information
  FortranDebug("-2-");
  MESSAGE("Saving all the main control variables\n");
  FortranFormat(10, INTFORMAT);

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

  //NCOPY but not sure why it is not here. PMEMDCUDA reads it though. Added for C4PairwiseCUDA
  //FortranWriteInt(0);

  // NUMC4Pairwise New2021
  if (iVarArrayElementCount(uUnit->vaC4Pairwise) > 0) {
    FortranWriteInt(0); // To make up the loss of NCOPY
    FortranWriteInt(iVarArrayElementCount(uUnit->vaC4Pairwise));
  }
  FortranEndLine();


  // -3-  write out the names of the atoms
  FortranDebug("-3-");
  MESSAGE("Writing the names of the atoms\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG ATOM_NAME");
  FortranWriteString("%FORMAT(20a4)");
  FortranFormat(20, LBLFORMAT);
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
  FortranWriteString("%FORMAT(5E16.8)");
  FortranFormat(5, DBLFORMAT);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
    FortranWriteDouble(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->dCharge * ELECTRONTOKCAL);
  }
  FortranEndLine();
  MESSAGE("Writing the atomic numbers\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG ATOMIC_NUMBER");
  FortranWriteString("%FORMAT(10I8)");
  FortranFormat(10, INTFORMAT);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
    FortranWriteInt(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->iElement);
  }
  FortranEndLine();

  // -5- write out the atomic masses
  FortranDebug("-5-");
  MESSAGE("Writing the atomic masses\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG MASS");
  FortranWriteString("%FORMAT(5E16.8)");
  FortranFormat(5, DBLFORMAT);
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
  FortranWriteString("%FORMAT(10I8)");
  FortranFormat(10, INTFORMAT);
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
  FortranWriteString("%FORMAT(10I8)");
  FortranFormat(10, INTFORMAT);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
    FortranWriteInt(*PVAI(vaExcludedCount, int, i));
  }
  FortranEndLine();

  // -8- Write the index for the position of the non bond type of each type
  FortranDebug("-8-");
  MESSAGE("writing position of the non bond type of each type\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG NONBONDED_PARM_INDEX");
  FortranWriteString("%FORMAT(10I8)");
  FortranFormat(10, INTFORMAT);
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
  FortranWriteString("%FORMAT(20a4)");
  FortranFormat(20, LBLFORMAT);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaResidues); i++) {
    cPTemp = PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->sName;
    if (strlen(cPTemp) > 3) {
      cPTemp += (strlen(cPTemp) - 3);
    }
    FortranWriteString(cPTemp);
  }
  FortranEndLine();

  // -10- Pointer list for all the residues
  FortranDebug("-10-");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG RESIDUE_POINTER");
  FortranWriteString("%FORMAT(10I8)");
  FortranFormat(10, INTFORMAT);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaResidues); i++) {
    FortranWriteInt(PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->iAtomStartIndex);
  }
  FortranEndLine();

  // -11- Force constants for bonds
  FortranDebug("-11-");
  MESSAGE("Writing bond force constants\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG BOND_FORCE_CONSTANT");
  FortranWriteString("%FORMAT(5E16.8)");
  FortranFormat(5, DBLFORMAT);
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
  FortranWriteString("%FORMAT(5E16.8)");
  FortranFormat(5, DBLFORMAT);
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
    FortranWriteString("%FORMAT(5E16.8)");
    FortranFormat(5, DBLFORMAT);
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
    FortranWriteString("%FORMAT(5E16.8)");
    FortranFormat(5, DBLFORMAT);
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
    FortranWriteString("%FORMAT(5E16.8)");
    FortranFormat(5, DBLFORMAT);
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
    FortranWriteString("%FORMAT(5E16.8)");
    FortranFormat(5, DBLFORMAT);
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
  FortranWriteString("%FORMAT(5E16.8)");
  FortranFormat(5, DBLFORMAT);
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
  FortranWriteString("%FORMAT(5E16.8)");
  FortranFormat(5, DBLFORMAT);
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
  FortranWriteString("%FORMAT(5E16.8)");
  FortranFormat(5, DBLFORMAT);
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
  FortranWriteString("%FORMAT(5E16.8)");
  FortranFormat(5, DBLFORMAT);
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
  FortranWriteString("%FORMAT(5E16.8)");
  FortranFormat(5, DBLFORMAT);
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
  FortranWriteString("%FORMAT(5E16.8)");
  FortranFormat(5, DBLFORMAT);
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
  FortranWriteString("%FORMAT(5E16.8)");
  FortranFormat(5, DBLFORMAT);
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
  FortranWriteString("%FORMAT(5E16.8)");
  FortranFormat(5, DBLFORMAT);
  for (i = 0; i < iParmSetTotalAtomParms(uUnit->psParameters); i++) {
    FortranWriteDouble(0.0);
  }
  FortranEndLine();

  // -19- Lennard jones r**12 term for all possible interactions
  // CN1 array
  FortranDebug("-19-");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG LENNARD_JONES_ACOEF");
  FortranWriteString("%FORMAT(5E16.8)");
  FortranFormat(5, DBLFORMAT);
  for (i = 0; i < iVarArrayElementCount(vaNBParameters); i++) {
    FortranWriteDouble(PVAI(vaNBParameters, NONBONDACt, i)->dA);
  }
  FortranEndLine();

  // -20- Lennard jones r**6 term for all possible interactions
  // CN2 array
  FortranDebug("-20-");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG LENNARD_JONES_BCOEF");
  FortranWriteString("%FORMAT(5E16.8)");
  FortranFormat(5, DBLFORMAT);
  for (i = 0; i < iVarArrayElementCount(vaNBParameters); i++) {
    FortranWriteDouble(PVAI(vaNBParameters, NONBONDACt, i)->dC);
  }
  FortranEndLine();

  // -21- Write the bond interactions that include hydrogen
  // Write the two indices into the atom table, then the index
  // into the interaction table
  FortranDebug("-21-");
  MESSAGE("Writing the bond interactions with hydrogens\n");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG BONDS_INC_HYDROGEN");
  FortranWriteString("%FORMAT(10I8)");
  FortranFormat(10, INTFORMAT);
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
  FortranWriteString("%FORMAT(10I8)");
  FortranFormat(10, INTFORMAT);
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
  FortranWriteString("%FORMAT(10I8)");
  FortranFormat(10, INTFORMAT);
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
  FortranWriteString("%FORMAT(10I8)");
  FortranFormat(10, INTFORMAT);
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
  FortranWriteString("%FORMAT(10I8)");
  FortranFormat(10, INTFORMAT);
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
  FortranWriteString("%FORMAT(10I8)");
  FortranFormat(10, INTFORMAT);
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
  FortranWriteString("%FORMAT(10I8)");
  FortranFormat(10, INTFORMAT);
  for (i = 0; i < iVarArrayElementCount(vaExcludedAtoms); i++) {
    FortranWriteInt(*PVAI(vaExcludedAtoms, int, i));
  }
  FortranEndLine();

  // -28- Write the R^12 term for the Hydrogen bond equation
  FortranDebug("-28-");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG HBOND_ACOEF");
  FortranWriteString("%FORMAT(5E16.8)");
  FortranFormat(5, DBLFORMAT);
  for (i = 0; i < iParmSetTotalHBondParms(uUnit->psParameters); i++) {
    ParmSetHBond(uUnit->psParameters, i, sType1, sType2, &dC, &dD, sDesc);
    FortranWriteDouble(dC);
  }
  FortranEndLine();

  // -29- Write the R^10 term for the Hydrogen bond equation
  FortranDebug("-29-");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG HBOND_BCOEF");
  FortranWriteString("%FORMAT(5E16.8)");
  FortranFormat(5, DBLFORMAT);
  for (i = 0; i < iParmSetTotalHBondParms(uUnit->psParameters); i++) {
    ParmSetHBond(uUnit->psParameters, i, sType1, sType2, &dC, &dD, sDesc);
    FortranWriteDouble(dD);
  }
  FortranEndLine();

  // -30- No longer used, but stored
  FortranDebug("-30-");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG HBCUT");
  FortranWriteString("%FORMAT(5E16.8)");
  FortranFormat(5, DBLFORMAT);
  for (i = 0; i < iParmSetTotalHBondParms(uUnit->psParameters); i++) {
    FortranWriteDouble(0.0);
  }
  FortranEndLine();

  // -31- List of atomic symbols
  FortranDebug("-31-");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG AMBER_ATOM_TYPE");
  FortranWriteString("%FORMAT(20a4)");
  FortranFormat(20, LBLFORMAT);
  //FortranFormat(10, INTFORMAT);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
    FortranWriteString(sAtomType(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom));
    //FortranWriteInt(iAtomCoordination(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom)); //VeryNew
  }
  FortranEndLine();

  // -32- List of tree symbols
  FortranDebug("-32-");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG TREE_CHAIN_CLASSIFICATION");
  FortranWriteString("%FORMAT(20a4)");
  FortranFormat(20, LBLFORMAT);
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
  FortranWriteString("%FORMAT(10I8)");
  FortranFormat(10, INTFORMAT);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
    FortranWriteInt(0);
  }
  FortranEndLine();

  // -34- Who knows, something to do with rotating atoms
  FortranDebug("-34-");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG IROTAT");
  FortranWriteString("%FORMAT(10I8)");
  FortranFormat(10, INTFORMAT);
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
    zUnitIOFindAndCountMolecules( uUnit );
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG SOLVENT_POINTERS");
    FortranWriteString("%FORMAT(3I8)");
    FortranFormat(3, INTFORMAT);
    FortranWriteInt(iTemp);
    FortranWriteInt(iVarArrayElementCount(uUnit->vaAtomsPerMolecule));

    // FORTRAN indexing conversion
    FortranWriteInt(uUnit->iFirstSolvent + 1);
    FortranEndLine();

    // -35B- The number of ATOMs in the Ith MOLECULE
    FortranDebug("-35B-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG ATOMS_PER_MOLECULE");
    FortranWriteString("%FORMAT(10I8)");
    FortranFormat(10, INTFORMAT);
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
      FortranFormat(4, DBLFORMAT);
      FortranWriteDouble(dX);
      FortranWriteDouble(dY);
      FortranWriteDouble(dZ);
      FortranWriteDouble(uUnit->dAlpha / DEGTORAD);
      FortranWriteDouble(uUnit->dBeta  / DEGTORAD);
      FortranWriteDouble(uUnit->dGamma / DEGTORAD);
    } else {
      FortranWriteString("%FLAG BOX_DIMENSIONS");
      FortranWriteString("%FORMAT(5E16.8)");
      FortranFormat(4, DBLFORMAT);
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
    FortranWriteString("%FORMAT(10I8)");
    FortranFormat(1, INTFORMAT);
    FortranWriteInt(uUnit->iCapTempInt);
    FortranEndLine();

    // -35E- CUTCAP, XCAP, YCAP, ZCAP
    FortranDebug("-35E-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG CAP_INFO2");
    FortranWriteString("%FORMAT(5E16.8)");
    FortranFormat(4, DBLFORMAT);
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
  FortranWriteString("%FORMAT(5E16.8)");
  FortranFormat(5, DBLFORMAT);
  iResidueIndex = -1;
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
    saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i);
    iIndex = iParmSetFindAtom(uUnit->psParameters, saPAtom->sType);
    ParmSetAtom(uUnit->psParameters, iIndex, sType, &dMass, &dPolar,
                &dEpsilon, &dRStar, &dEpsilon14, &dRStar14, &dScreenF,
                &iElement, &iHybridization, sDesc);
    if (GDefaults.iGBparm < 3 || GDefaults.iGBparm == 6 || GDefaults.iGBparm == 8) {

      // Bondi or modified Bondi radii
      switch (iElement) {
        case 1:
          dGBrad = 1.2;

          // Make the modifications that hydrogen radii
          // depend upon the atoms they are bonded to.
          // iGBparm=1 corresponds to Amber 6, JACS 122:2489 (2000);
          // iGBparm=2 adds the update of Biopolymers 56: 275 (2001)
          if (iAtomCoordination(saPAtom->aAtom) > 0) {

	    // For multiply bonded Hydrogen atoms use the first bond for
            // determining modified GB radii.  WAT contains multiply bonded
	    // Hydrogen atoms so do not emit a warning.
            aAtomA = aAtomBondedNeighbor(saPAtom->aAtom, 0);
            if (GDefaults.iGBparm == 1 || GDefaults.iGBparm == 2) {
              switch (aAtomA[0].iAtomicNumber) {
                case 6:
                  dGBrad = 1.3;
                  break;      // Carbon
                case 8:
                  dGBrad = 0.8;
                  break;      // Oxygen
                case 16:
                  dGBrad = 0.8;
                  break;      // Sulfur
                case 7:
                  if (GDefaults.iGBparm == 2) {
                    dGBrad = 1.3;
                  }
                  break;      // Nitrogen, mbondi
                case 1:        // Special case: water hydrogen
                  if ((sAtomType(aAtomA)[0] == 'H' || sAtomType(aAtomA)[0] == 'h') &&
                      (sAtomType(aAtomA)[1] == 'W' || sAtomType(aAtomA)[1] == 'w')) {
                                dGBrad = 0.8;
		  }
                  break;
              }
            }
	    else if (GDefaults.iGBparm == 6 || GDefaults.iGBparm == 8) {

              // Try Alexey's scheme
              if (aAtomA[0].iAtomicNumber == 7) {
                dGBrad = 1.3;
                if (GDefaults.iGBparm == 8) {

                  // Update residue as appropriate
                  if (saPAtom->iResidueIndex != iResidueIndex) {
                    iResidueIndex = saPAtom->iResidueIndex;
                    cPTemp = PVAI(uUnit->vaResidues, SAVERESIDUEt,
                                  iResidueIndex - 1)->sName;
                    if (strlen(cPTemp) > 3) {
                      cPTemp += (strlen(cPTemp) - 3);
		    }
                  }
		  // Adjust Arg HH and HE
                  if (!strcmp(cPTemp, "ARG") &&
		      !(strncmp(sAtomName(saPAtom->aAtom), "HH", 2) &&
                        strcmp(sAtomName(saPAtom->aAtom), "HE"))) {
                    dGBrad = 1.17;
                  }
                }
              }
            }
          }
          else {
            VPWARN("Unbonded Hydrogen atom %s in %s.\n"
                    " Cannot determine the requested GB radius for this atom.\n"
                    " Writing the unmodified Bondi GB radius.\n",
                    saPAtom->aAtom->cHeader.sName,
                    saPAtom->aAtom->cHeader.cContainedBy->sName);
          }
          break;
        case 6:

	  // Use the mass of the carbon atom. We are testing for
          // carbons here. C1 == CH, C2 == CH2, C3 == CH3. UA carbons
          // have a larger radius (2.2), so we want to make sure that
          // the C1, C2, and C3 atom types _really_ correspond to UA
          // UA carbons. C1 atoms should have a mass of 12.01 + 1.01,
          // C2 should be 12.01 + 2.02, and C3 should be 12.01 + 3.03.
          // This mneumonic will not work for 13C named "C1". This is
          // a (hopefully) temporary kludge.
	  if (strncmp(sType, "C1", 2) && strncmp(sType, "C2", 2) && strncmp(sType, "C3", 2)) {
            dGBrad = 1.7;
	  }
          else if (!strncmp(sType, "C1", 2) && dMass < 13.0) {
            dGBrad = 1.7;
	  }
	  else if (!strncmp(sType, "C2", 2) && dMass < 14.0) {
	    dGBrad = 1.7;
	  }
	  else if (!strncmp(sType, "C3", 2) && dMass < 15.0) {
	    dGBrad = 1.7;
	  }
          else {
            dGBrad = 2.2;
	  }
          break;
        case 7:
          dGBrad = 1.55;
          break;
        case 8:
          dGBrad = 1.5;
          if (GDefaults.iGBparm == 8) {

            // Update residue as appropriate
            if (saPAtom->iResidueIndex != iResidueIndex) {
              iResidueIndex = saPAtom->iResidueIndex;
              cPTemp = PVAI(uUnit->vaResidues, SAVERESIDUEt, iResidueIndex - 1)->sName;
              if (strlen(cPTemp) > 3) {
                cPTemp += (strlen(cPTemp) - 3);
	      }
            }
            // Adjust Asp OD and Glu OE, and terminal OXT
            if (!(strcmp(cPTemp, "ASP") || strncmp(sAtomName(saPAtom->aAtom), "OD", 2)) ||
                  !(strcmp(cPTemp, "AS4") || strncmp(sAtomName(saPAtom->aAtom), "OD", 2)) ||
                  !(strcmp(cPTemp, "GLU") || strncmp(sAtomName(saPAtom->aAtom), "OE", 2)) ||
		  !(strcmp(cPTemp, "GL4") || strncmp(sAtomName(saPAtom->aAtom), "OE", 2)) ||
                  (!strcmp(sAtomName(saPAtom->aAtom), "OXT") ||
                  (i + 1 < iVarArrayElementCount(uUnit->vaAtoms) &&
                   !strcmp(sAtomName(PVAI(uUnit->vaAtoms, SAVEATOMt, i + 1)->aAtom),
			   "OXT")))) {
	      dGBrad = 1.4;
	    }
          }
          break;
        case 9:
          dGBrad = 1.5;
          break;
        case 14:
          dGBrad = 2.1;
          break;
        case 15:
          dGBrad = 1.85;
          break;
        case 16:
          dGBrad = 1.8;
          break;
        case 17:
          dGBrad = 1.7;
          break;
        default:
          dGBrad = 1.5;
          break;
      }
    }
    else if (GDefaults.iGBparm == 3) {

      // Radii from Huo & Kollman
      switch (iElement) {
        case 1:
          dGBrad = 1.15;
          break;
        case 6:
          dGBrad = 1.85;
          break;
        case 7:
          dGBrad = 1.64;
          break;
        case 8:
          dGBrad = 1.53;
          break;
        case 9:
          dGBrad = 1.53;
          break;
        case 15:
          dGBrad = 2.02;
          break;
        case 16:
          dGBrad = 2.00;
          break;
        case 17:
          dGBrad = 1.97;
          break;
        case 35:
          dGBrad = 2.03;
          break;
        case 53:
          dGBrad = 2.10;
          break;
        default:
          dGBrad = 1.5;
          break;
      }
    }
    else if (GDefaults.iGBparm == 7) {

      // Parse radii
      switch (iElement) {
        case 1:
          dGBrad = 1.00;
          break;
        case 6:
          dGBrad = 1.70;
          break;
        case 7:
          dGBrad = 1.50;
          break;
        case 8:
          dGBrad = 1.40;
          break;
        case 16:
          dGBrad = 1.85;
          break;
        default:
          dGBrad = 1.50;
          break;
          // Radii from J. Phys. Chem. 1994, 98, 1978-1988
      }
    }
    FortranWriteDouble(dGBrad);
  }
  FortranEndLine();

  // Write out the GB screening parameters
  MESSAGE("Writing the GB screening parameters\n");
  dGBscreen = 0.0;
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG SCREEN");
  FortranWriteString("%FORMAT(5E16.8)");
  FortranFormat(5, DBLFORMAT);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
    saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i);
    iIndex = iParmSetFindAtom(uUnit->psParameters, saPAtom->sType);
    ParmSetAtom(uUnit->psParameters, iIndex, sType, &dMass, &dPolar,
                &dEpsilon, &dRStar, &dEpsilon14, &dRStar14, &dScreenF,
                &iElement, &iHybridization, sDesc);
    if (GDefaults.iGBparm < 4 || GDefaults.iGBparm == 6 || GDefaults.iGBparm == 8) {

      // For now, hardwire the Bondi radii
      switch (iElement) {
        case 1:
          dGBscreen = 0.85;
          break;
        case 6:
          dGBscreen = 0.72;
          break;
        case 7:
          dGBscreen = 0.79;
          break;
        case 8:
          dGBscreen = 0.85;
          break;
        case 9:
          dGBscreen = 0.88;
          break;
        case 15:
          dGBscreen = 0.86;
          break;
        case 16:
          dGBscreen = 0.96;
          break;
        default:
          dGBscreen = 0.8;
          break;          // or should fail??
      }
    }
    else if (GDefaults.iGBparm == 4) {    // param for Jayaram et al. 'GB'
      switch (iElement) {
        case 1:
          dGBscreen = 0.8461;
          break;
        case 6:
          dGBscreen = 0.9615;
          break;
        case 7:
          dGBscreen = 0.9343;
          break;
        case 8:
          dGBscreen = 1.0088;
          break;
        case 11:
          dGBscreen = 1.0000;
          break;
        case 12:
          dGBscreen = 1.0000;
          break;          // Set by HG
        case 15:
          dGBscreen = 1.0700;
          break;
        case 16:
          dGBscreen = 1.1733;
          break;
        default:
          dGBscreen = 0.8000;
          break;          // Set by HG
      }
    }
    else if (GDefaults.iGBparm == 5) {

      // Param for Jayaram et al. 'MGB'
      switch (iElement) {
        case 1:
          dGBscreen = 0.8846;
          break;
        case 6:
          dGBscreen = 0.9186;
          break;
        case 7:
          dGBscreen = 0.8733;
          break;
        case 8:
          dGBscreen = 0.8836;
          break;
        case 11:
          dGBscreen = 1.0000;
          break;
        case 12:
          dGBscreen = 1.0000;
          break;          // Set by HG
        case 15:
          dGBscreen = 0.9604;
          break;
        case 16:
          dGBscreen = 0.9323;
          break;
        default:
          dGBscreen = 0.8000;
          break;          // Set by HG
      }
    }
    FortranWriteDouble(dGBscreen);
  }
  FortranEndLine();

  // Write IPOL near the end of prmtop
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG IPOL");
  FortranWriteString("%FORMAT(1I8)");
  FortranFormat(1, INTFORMAT);
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
    FortranWriteString("%FORMAT(10I8)");
    FortranFormat(10, INTFORMAT);
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
    FortranWriteString("%FORMAT(10I8)");
    FortranFormat(10, INTFORMAT);

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
    FortranWriteString("%FORMAT(10I8)");
    FortranFormat(10, INTFORMAT);
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
    FortranWriteString("%FORMAT(10I8)");
    FortranFormat(10, INTFORMAT);

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
    FortranWriteString("%FORMAT(10I8)");
    FortranFormat(10, INTFORMAT);
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
    FortranWriteString("%FORMAT(10I8)");
    FortranFormat(10, INTFORMAT);

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
    FortranWriteString("%FORMAT(20a4)");
    FortranFormat(20, LBLFORMAT);
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
    FortranWriteString("%FORMAT(20a4)");
    FortranFormat(20, LBLFORMAT);
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
    FortranWriteString("%FORMAT(20a4)");
    FortranFormat(20, LBLFORMAT);
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
    FortranWriteString("%FORMAT(5E16.8)");
    FortranFormat(5, DBLFORMAT);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
      FortranWriteDouble(0.0);
    }
    FortranEndLine();

    // -36K- Flag to tell whether the atom is perturbed
    FortranDebug("-36K-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG IAPER");
    FortranWriteString("%FORMAT(10I8)");
    FortranFormat(10, INTFORMAT);
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
    FortranWriteString("%FORMAT(10I8)");
    FortranFormat(10, INTFORMAT);
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
    FortranWriteString("%FORMAT(5E16.8)");
    FortranFormat(5, DBLFORMAT);
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
    FortranWriteString("%FORMAT(5E16.8)");
    FortranFormat(5, DBLFORMAT);
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
    FortranWriteString("%FORMAT(5E16.8)");
    FortranFormat(5, DBLFORMAT);
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
      FortranWriteString("%FORMAT(5E16.8)");
      FortranFormat(5, DBLFORMAT);
      saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, 0);
      for (i = 0; i < iMax; i++, saPAtom++) {
        BOOL bTmp = bAtomPerturbed(saPAtom->aAtom);
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

  // NewT
  // -37- Lennard jones r**4 term for all possible interactions
  // CN1 array
  // Only print when tleap has a "addC4Type" command enabled
  if (iC4count) {
    FortranDebug("-37-");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG LENNARD_JONES_CCOEF");
    FortranWriteString("%FORMAT(5E16.8)");
    FortranFormat(5, DBLFORMAT);
    for (i = 0; i < iVarArrayElementCount(vaNBParameters); i++) {
      FortranWriteDouble(PVAI(vaNBParameters, NONBONDACt, i)->d4);
    }
    FortranEndLine();
  }

  // Add C4 output here.
  // -38- Write atom-specific C4 pairwise interaction
  // Write the two indices into the atom table, then the index
  // into the interaction table
  //VP0("vaC4Pairwise count is: %i\n", iVarArrayElementCount(uUnit->vaC4Pairwise));
  if (iVarArrayElementCount(uUnit->vaC4Pairwise) > 0) {
    FortranDebug("-38-");
    MESSAGE("Writing the atom-specific C4 pairwise indices\n");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG LENNARD_JONES_DCOEF");
    FortranWriteString("%FORMAT(3I8)");
    FortranFormat(10, INTFORMAT);
    //FortranWriteString("%FORMAT ATOM, ATOM, C4");
    //FortranFormat(5, DBLFORMAT);
    //VP0("vaC4Pairwise count is: %i\n", iVarArrayElementCount(uUnit->vaC4Pairwise));
    for (i = 0; i < iVarArrayElementCount(uUnit->vaC4Pairwise); i++) {
      scPC4Pairwise = PVAI(uUnit->vaC4Pairwise, SAVEC4Pairwiset, i);
      FortranFormat(10, INTFORMAT);
      FortranWriteInt(AMBERINDEX(scPC4Pairwise->iAtom1));
      FortranWriteInt(AMBERINDEX(scPC4Pairwise->iAtom2));
      FortranWriteInt(i+1);
      FortranEndLine();
    }

    FortranDebug("-39-");
    MESSAGE("Writing the atom-specific C4 pairwise values\n");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG LENNARD_JONES_DVALUE");
    FortranWriteString("%FORMAT(1E16.8)");
    FortranFormat(5, DBLFORMAT);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaC4Pairwise); i++) {
      scPC4Pairwise = PVAI(uUnit->vaC4Pairwise, SAVEC4Pairwiset, i);
      FortranWriteDouble(scPC4Pairwise->daC4Pairwise);
      FortranEndLine();
    }
  }
/*
  FortranDebug("-37-");
tring("%FORMAT(10I8)");
  FortranFormat(1, "%-80s");
  FortranWriteString("%FLAG LENNARD_JONES_DCOEF");
  FortranWriteString("%FORMAT(5E16.8)");
  //FortranFormat(20, LBLFORMAT);
  //FortranFormat(10, INTFORMAT);
  //FortranFormant(5, DBLFORMAT);
  for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
    //FortranWriteString(sAtomType(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom));
    if (iAtomC4Pairwise(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom) != 0) {
      saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i);
      //FortranWriteInt(iAtomC4Pairwise(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom)); //VeryNew
      FortranFormat(10, INTFORMAT);
      FortranWriteInt(saPAtom->iResidueIndex); //VeryNew
      FortranWriteInt(saPAtom->iSequence); //VeryNew
      FortranWriteInt(PVAI(uUnit->vaResidues, SAVERESIDUEt, saPAtom->iResidueIndex)->iAtomStartIndex + saPAtom->iSequence - 1);
      //FortranWriteInt(i);
      FortranFormat(20, LBLFORMAT);
      FortranWriteString("  ");
      FortranWriteString(sAtomName(saPAtom->aAtom));
      FortranEndLine();
      //TODO:Nested loop to print all related C4s if more than one,
      //then set self C4Coord to zero so it won't be re-visited by the other end.
      for (int j = 0; j < iAtomC4Pairwise(saPAtom->aAtom); j++) {
        FortranFormat(20, LBLFORMAT);
        FortranWriteString(sAtomType(saPAtom->aAtom->aaC4Pairwise[j]));
        FortranFormat(10, INTFORMAT);
        FortranWriteInt(saPAtom->aAtom->aaC4Pairwise[j]->iIndex);
        FortranFormat(5, DBLFORMAT);
        FortranWriteDouble(saPAtom->aAtom->daC4Pairwise[j]);
        FortranEndLine();
      }
      //FortranWriteInt(AMBERINDEX(uUnit->vaAtoms)); //VeryNew
    }
    //FortranWriteDouble(aAtomC4Pairwise(PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom, 0)); //New
  }
*/

  //  Charmm-style parameters
  if (GDefaults.iCharmm) {

    // -19- Lennard jones r**12 term for all 14 interactions
    // CN114 array
    FortranDebug("-19-");
    FortranFormat(5, DBLFORMAT);
    for (i = 0; i < iVarArrayElementCount(vaNBParameters); i++) {
      FortranWriteDouble(PVAI(vaNBParameters, NONBONDACt, i)->dA14);
    }
    FortranEndLine();

    // -20- Lennard jones r**6 term for all 14 interactions
    // CN214 array
    FortranDebug("-20-");
    FortranFormat(5, DBLFORMAT);
    for (i = 0; i < iVarArrayElementCount(vaNBParameters); i++) {
      FortranWriteDouble(PVAI(vaNBParameters, NONBONDACt, i)->dC14);
    }
    FortranEndLine();

    // -13- Force constants for Urey-Bradley
    FortranDebug("-13-");
    FortranFormat(5, DBLFORMAT);
    for (i = 0; i < iParmSetTotalAngleParms(uUnit->psParameters); i++) {
      ParmSetAngle(uUnit->psParameters, i, sAtom1, sAtom2, sAtom3, &dKt,
                   &dT0, &dTkub, &dRkub, sDesc);
      FortranWriteDouble(dTkub);
    }
    FortranEndLine();

    // -14- Equilibrium distances for Urey-Bradley
    FortranDebug("-14-");
    FortranFormat(5, DBLFORMAT);
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
    FortranWriteString("%FORMAT(10I8)");
    FortranFormat(10, "%8d");
    for (i = 0; i < iVarArrayElementCount(uUnit->vaResidues); i++) {
      FortranWriteInt( PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->iPdbResSeq );
    }
    FortranEndLine();

    MESSAGE("Writing residue PDB ChainId\n");
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG RESIDUE_CHAINID");
    FortranWriteString("%FORMAT(20a4)");
    FortranFormat(20, LBLFORMAT);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaResidues); i++) {
      FortranWriteString(PVAI(uUnit->vaResidues, SAVERESIDUEt, i)->sChainId);
    }
    FortranEndLine();
  }

  // CMAP parameters, Mengjuei Hsieh and Yong Duan
  SaveAmberParmCMAP(uUnit, fOut);

  // Write the coordinate file
  if (bNetcdf == TRUE) {
    zUnitIOSaveAmberNetcdf(uUnit, crdName);

    // Clean up arrays and return, netcdf file will be written later
    VarArrayDestroy(&vaNBIndexMatrix);
    VarArrayDestroy(&vaNBParameters);
    VarArrayDestroy(&vaExcludedAtoms);
    VarArrayDestroy(&vaExcludedCount);
    VarArrayDestroy(&vaNBIndex);
    VarArrayDestroy(&vaNonBonds);

    return;
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
  VarArrayDestroy(&vaNBIndexMatrix);
  VarArrayDestroy(&vaNBParameters);
  VarArrayDestroy(&vaExcludedAtoms);
  VarArrayDestroy(&vaExcludedCount);
  VarArrayDestroy(&vaNBIndex);
  VarArrayDestroy(&vaNonBonds);
  fclose(fCrd);
}

#if 0
/*
 *      zUnitIOSaveAmberParmFormat_old
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
 *
 *TODO: Add RESTRAINT code
 *TODO: Add CAP information
 */
#undef INTFORMAT
//#define AMBERINDEX(i)   3*(i-1)
#define INTFORMAT       "%6d"
//#define DBLFORMAT       "%16.8lE"
//#define LBLFORMAT       "%-4s"
//#define ELECTRONTOKCAL  18.2223

void
zUnitIOSaveAmberParmFormat_old( UNIT uUnit, FILE *fOut, char *crdName,
        BOOL bPolar, BOOL bPert, char sA[8][16], char sB[8][16], double daC4Type[16], int iC4count) //NewT
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
BOOL            bFoundSome;
VECTOR          vPos;
char            *cPTemp;
double          dX, dY, dZ, dEpsilon, dRStar, dEpsilon14, dRStar14;
STRING          sDesc, sType;
IX_REC          eResEnt = {NULL};
IX_DESC         iResIx;

    // Open the coordinate file
    FILE *fCrd = FOPENCOMPLAIN( crdName, "w" );
    if ( fOut == NULL ) {
        VP0("Could not open file: %s\n", crdName );
    }

                /* Build the excluded atom list */

    MESSAGE("Building the excluded atom list\n" );
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
    //VP0("Not Marking per-residue atom chain types.\n" );
    //iMaxAtoms = 0;

    create_index(&iResIx, IX_DUPKEYREC, IX_LEN_CSTRING);

    VP0("Marking per-residue atom chain types.\n" );

    iMaxAtoms = 0;
    lTemp = lLoop( (OBJEKT)uUnit, RESIDUES );
    while ( (rRes = (RESIDUE)oNext(&lTemp)) ) {
        int     iAtoms = iMarkMainChainAtoms( rRes, 0 );
        if ( iAtoms > 0 )
                zMarkSideChains( rRes );
        if ( iAtoms < 0 ) {
                iAtoms = -iAtoms;
                /*
                 *  couldn't mark main chains
                 */

                strcpy( eResEnt.key, rRes->cHeader.sName );
                if ( add_key( &eResEnt, &iResIx ) != IX_OK )
                        DFATAL("add_key() residue chain\n" );
        }
        if ( iAtoms > iMaxAtoms )
                iMaxAtoms = iAtoms;
    }
    /*
     *  print warnings
     */
    first_key( &iResIx );
    i = 1;
    while ( next_key( &eResEnt, &iResIx ) == IX_OK ) {
        if ( i ) {
                VP0("  (Residues lacking connect0/connect1 - \n" );
                VP0("   these don't have chain types marked:\n\n" );
                VP0("\tres\ttotal affected\n\n" );
                i = 0;
        }
        VP0("\t%s\t%d\n", eResEnt.key, eResEnt.count);
    }
    if (!i)
        VP0("  )\n" );
    destroy_index( &iResIx );

                /* Build the NON-BOND arrays that AMBER needs */

    zUnitIOBuildNonBondArrays( uUnit, &vaNBIndexMatrix, &vaNBParameters,
                                        &vaNBIndex, &vaNonBonds, sA, sB, daC4Type, iC4count ); //NewT


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
    MESSAGE("Saving the name of the UNIT\n" );
    FortranFormat( 1, "%-80s" );
    FortranWriteString( sContainerName(uUnit) );

        /* -2- Save control information */
    FortranDebug( "-2-" );
    MESSAGE("Saving all the main control variables\n" );
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
        VP0(" Restraints:  Bond %d  Angle %d  Torsion %d\n",
                iUnitRestraintTypeCount( uUnit, RESTRAINTBOND ),
                iUnitRestraintTypeCount( uUnit, RESTRAINTANGLE ),
                iUnitRestraintTypeCount( uUnit, RESTRAINTTORSION ) );
    else
        VP0(" (no restraints)\n" );

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
                MESSAGE("Boundary pert bond %d-%d\n",
                                sbPBond->iAtom1, sbPBond->iAtom2 );
                iCountBondBoundary++;
            }
        }
    }

    MESSAGE("Perturbed bonds: %d\n", iCountBondPerturbed );
    MESSAGE("Perturbed boundary bonds: %d\n", iCountBondBoundary );

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
    printf("iMaxAoms (2) %i\n",iMaxAtoms);
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

    MESSAGE("Writing the names of the atoms\n" );
    FortranFormat( 20, LBLFORMAT );
    for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
        cPTemp = PVAI( uUnit->vaAtoms, SAVEATOMt, i )->sName;
        if ( strlen( cPTemp ) > 4 ) cPTemp += ( strlen(cPTemp)-4 );
        FortranWriteString( cPTemp );
    }
    FortranEndLine();

        /* -4- write out the atomic charges */
    FortranDebug( "-4-" );

    MESSAGE("Writing the atomic charges\n" );
    FortranFormat( 5, DBLFORMAT );
    for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
        FortranWriteDouble( PVAI( uUnit->vaAtoms, SAVEATOMt, i )->dCharge *
                                ELECTRONTOKCAL );
    }
    FortranEndLine();

        /* -5- write out the atomic masses */
    FortranDebug( "-5-" );

    MESSAGE("Writing the atomic masses\n" );
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

    MESSAGE("Writing the atomic types\n" );
    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
        iAtom = PVAI( uUnit->vaAtoms, SAVEATOMt, i )->iTypeIndex-1;
        iTemp = *PVAI( vaNBIndex, int, iAtom );
        FortranWriteInt( iTemp+1 );
    }
    FortranEndLine();

        /* -7- write out the starting index into the excluded atom list */
    FortranDebug( "-7-" );

    MESSAGE("Writing the starting index into the excluded atom list\n" );
    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount(uUnit->vaAtoms); i++ ) {
        FortranWriteInt( *PVAI( vaExcludedCount, int, i ) );
    }
    FortranEndLine();

        /* -8- Write the index for the position of the non bond type */
                /* of each type */
    FortranDebug( "-8-" );

    MESSAGE("writing position of the non bond type of each type\n");
    FortranFormat( 12, INTFORMAT );
    for ( i=0; i<iVarArrayElementCount(vaNBIndexMatrix); i++ ) {
        FortranWriteInt( *PVAI( vaNBIndexMatrix, int, i ) );
    }
    FortranEndLine();

        /* -9- Residue labels */
                /* Trim the string down to at most 3 characters by */
                /* taking the last three characters if it is too long */
    FortranDebug( "-9-" );

    MESSAGE("Writing the residue labels\n" );
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

    MESSAGE("Writing bond force constants\n" );
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

    MESSAGE("Writing equilibrium bond lengths\n" );
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

    MESSAGE("Writing angle force constants\n" );
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

    MESSAGE("Writing equilibrium angle values\n" );
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

    MESSAGE("Writing torsional force constants\n" );
    FortranFormat( 5, DBLFORMAT );
    MESSAGE("There are %d torsions and %d impropers\n",
                iParmSetTotalTorsionParms(uUnit->psParameters),
                iParmSetTotalImproperParms(uUnit->psParameters) );
    for ( i=0; i<iParmSetTotalTorsionParms(uUnit->psParameters); i++ ) {
        ParmSetTorsion( uUnit->psParameters, i, sAtom1, sAtom2,
                        sAtom3, sAtom4,
                        &iN, &dKp, &dP0, &dScEE, &dScNB, sDesc );
        MESSAGE("Torsion %d  %s-%s-%s-%s %d %lf %lf\n",
                        i, sAtom1, sAtom2, sAtom3, sAtom4,
                        iN, dKp, dP0 );
        FortranWriteDouble( dKp );
    }
    for ( i=0; i<iParmSetTotalImproperParms(uUnit->psParameters); i++ ) {
        ParmSetImproper( uUnit->psParameters, i, sAtom1, sAtom2,
                        sAtom3, sAtom4,
                        &iN, &dKp, &dP0, sDesc );
        MESSAGE("Improper %d  %s-%s-%s-%s %d %lf %lf\n",
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

    MESSAGE("Writing multiplicity of torsion interaction\n" );
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

    MESSAGE("Writing phase for torsion interactions\n" );
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
        FortranWriteDouble( PVAI( vaNBParameters, NONBONDACt, i )->dC );
    }
    FortranEndLine();

        /* -21- Write the bond interactions that include hydrogen */
                /* Write the two indices into the atom table, then the index */
                /* into the interaction table */
    FortranDebug( "-21-" );

    MESSAGE("Writing the bond interactions with hydrogens\n" );
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

    MESSAGE("Writing the bond interactions without hydrogens\n" );
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

    MESSAGE("Writing the angle interactions with hydrogens\n" );
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

    MESSAGE("Writing the angle interactions without hydrogens\n" );
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

    MESSAGE("Writing the torsion interactions with hydrogens\n" );
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
                MESSAGE("Had to turn torsion around to avoid K,L == 0\n" );
                MESSAGE("Outer atoms: %s --- %s\n",
                                sContainerName(aA), sContainerName(aD) );
                MESSAGE("Old order %d %d %d %d\n",
                                stPTorsion->iAtom1,
                                stPTorsion->iAtom2,
                                stPTorsion->iAtom3,
                                stPTorsion->iAtom4 );
                SWAP( stPTorsion->iAtom1, stPTorsion->iAtom4, iTemp );
                SWAP( stPTorsion->iAtom2, stPTorsion->iAtom3, iTemp );
                MESSAGE("New order %d %d %d %d\n",
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

    MESSAGE("Writing the torsion interactions without hydrogens\n" );
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
                MESSAGE("Had to turn torsion to avoid K,L == 0\n" );
                MESSAGE("Outer atoms: %s --- %s\n",
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
                    MESSAGE("Had to turn RESTRAINT torsion around to avoid\n" );
                    MESSAGE("K,L == 0\n" );
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

    MESSAGE("Writing the excluded atom list\n" );
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
        iTemp = i;

        /*
         *  Find the molecules and return the number of ATOMs in each
         *  molecule, along with the index of the first solvent molecule
         */

        zUnitIOFindAndCountMolecules( uUnit );

        FortranFormat( 3, INTFORMAT );
        FortranWriteInt( iTemp );
        FortranWriteInt( iVarArrayElementCount(uUnit->vaAtomsPerMolecule) );
        FortranWriteInt( uUnit->iFirstSolvent+1 );     /* FORTRAN index */

        FortranEndLine();

                /* -35B- The number of ATOMs in the Ith RESIDUE */

        FortranDebug( "-35B-" );
        FortranFormat( 12, INTFORMAT );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAtomsPerMolecule); i++ ) {
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
        }
        if ( iVarArrayElementCount(uUnit->vaBonds) ) {
                sbPBond = PVAI( uUnit->vaBonds, SAVEBONDt, 0 );
            for ( i=0; i<iVarArrayElementCount(uUnit->vaBonds); i++, sbPBond++ ) {
                if ( ((sbPBond->fFlags&PERTURBED)!=0)
                     && ((sbPBond->fFlags&BOUNDARY)!=0) ) {
                    FortranWriteInt( AMBERINDEX(sbPBond->iAtom1) );
                    FortranWriteInt( AMBERINDEX(sbPBond->iAtom2) );
                }
            }
        }
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
        }
        if ( iVarArrayElementCount(uUnit->vaBonds) ) {
                sbPBond = PVAI( uUnit->vaBonds, SAVEBONDt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaBonds); i++, sbPBond++ ) {
            if ( ((sbPBond->fFlags&PERTURBED)!=0)
                && ((sbPBond->fFlags&BOUNDARY)!=0) ) {
                FortranWriteInt( sbPBond->iParmIndex );
            }
        }
        }

                        /* Then LAMBDA = 1 */

        if ( iVarArrayElementCount(uUnit->vaBonds) ) {
                sbPBond = PVAI( uUnit->vaBonds, SAVEBONDt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaBonds); i++, sbPBond++ ) {
            if ( ((sbPBond->fFlags&PERTURBED)!=0)
                && ((sbPBond->fFlags&BOUNDARY)==0) ) {
                FortranWriteInt( sbPBond->iPertParmIndex );
            }
        }
        }
        if ( iVarArrayElementCount(uUnit->vaBonds) ) {
                sbPBond = PVAI( uUnit->vaBonds, SAVEBONDt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaBonds); i++, sbPBond++ ) {
            if ( ((sbPBond->fFlags&PERTURBED)!=0)
                && ((sbPBond->fFlags&BOUNDARY)!=0) ) {
                FortranWriteInt( sbPBond->iPertParmIndex );
            }
        }
        }
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
        }
        if ( iVarArrayElementCount(uUnit->vaAngles) ) {
                saPAngle = PVAI( uUnit->vaAngles, SAVEANGLEt, 0 );
            for ( i=0; i<iVarArrayElementCount(uUnit->vaAngles); i++, saPAngle++ ) {
                if ( ((saPAngle->fFlags&PERTURBED)!=0)
                    && ((saPAngle->fFlags&BOUNDARY)!=0) ) {
                    FortranWriteInt( AMBERINDEX(saPAngle->iAtom1) );
                    FortranWriteInt( AMBERINDEX(saPAngle->iAtom2) );
                    FortranWriteInt( AMBERINDEX(saPAngle->iAtom3) );
                }
            }
        }
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
        }
        if ( iVarArrayElementCount(uUnit->vaAngles) ) {
                saPAngle = PVAI( uUnit->vaAngles, SAVEANGLEt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAngles); i++, saPAngle++ ) {
            if ( ((saPAngle->fFlags&PERTURBED)!=0)
                && ((saPAngle->fFlags&BOUNDARY)!=0) ) {
                FortranWriteInt( saPAngle->iParmIndex );
            }
        }
        }

                        /* Then LAMBDA = 1 */

        if ( iVarArrayElementCount(uUnit->vaAngles) ) {
                saPAngle = PVAI( uUnit->vaAngles, SAVEANGLEt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAngles); i++, saPAngle++ ) {
            if ( ((saPAngle->fFlags&PERTURBED)!=0)
                && ((saPAngle->fFlags&BOUNDARY)==0) ) {
                FortranWriteInt( saPAngle->iPertParmIndex );
            }
        }
        }
        if ( iVarArrayElementCount(uUnit->vaAngles) ) {
                saPAngle = PVAI( uUnit->vaAngles, SAVEANGLEt, 0 );
        for ( i=0; i<iVarArrayElementCount(uUnit->vaAngles); i++, saPAngle++ ) {
            if ( ((saPAngle->fFlags&PERTURBED)!=0)
                && ((saPAngle->fFlags&BOUNDARY)!=0) ) {
                FortranWriteInt( saPAngle->iPertParmIndex );
            }
        }
        }
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
                        MESSAGE("Had to turn torsion around to avoid K,L == 0\n" );
                        MESSAGE("Outer atoms: %s --- %s\n",
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
        }
        if ( iVarArrayElementCount(uUnit->vaTorsions) ) {
                stPTorsion = PVAI( uUnit->vaTorsions, SAVETORSIONt, 0 );
            for ( i=0; i<iVarArrayElementCount(uUnit->vaTorsions);
                                                        i++, stPTorsion++ ) {
                if ( ((stPTorsion->fFlags&PERTURBED)!=0)
                    && ((stPTorsion->fFlags&BOUNDARY)!=0) ) {

                    if ( (AMBERINDEX(stPTorsion->iAtom3) == 0) ||
                         (AMBERINDEX(stPTorsion->iAtom4) == 0) ) {
                        MESSAGE("Had to turn torsion around to avoid K,L == 0\n" );
                        MESSAGE("Outer atoms: %s --- %s\n",
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
        }
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
        }
        if ( iVarArrayElementCount(uUnit->vaTorsions) ) {
                stPTorsion = PVAI( uUnit->vaTorsions, SAVETORSIONt, 0 );
            for ( i=0; i<iVarArrayElementCount(uUnit->vaTorsions);
                                                            i++, stPTorsion++ ) {
                if ( ((stPTorsion->fFlags&PERTURBED)!=0)
                    && ((stPTorsion->fFlags&BOUNDARY)!=0) ) {
                    FortranWriteInt( stPTorsion->iParmIndex );
                }
            }
        }

                        /* Then LAMBDA = 1 */

        if ( iVarArrayElementCount(uUnit->vaTorsions) ) {
                stPTorsion = PVAI( uUnit->vaTorsions, SAVETORSIONt, 0 );
            for ( i=0; i<iVarArrayElementCount(uUnit->vaTorsions);
                                                            i++, stPTorsion++ ) {
                if ( ((stPTorsion->fFlags&PERTURBED)!=0)
                    && ((stPTorsion->fFlags&BOUNDARY)==0) ) {
                    FortranWriteInt( stPTorsion->iPertParmIndex );
                }
            }
        }
        if ( iVarArrayElementCount(uUnit->vaTorsions) ) {
                stPTorsion = PVAI( uUnit->vaTorsions, SAVETORSIONt, 0 );
            for ( i=0; i<iVarArrayElementCount(uUnit->vaTorsions);
                                                            i++, stPTorsion++ ) {
                if ( ((stPTorsion->fFlags&PERTURBED)!=0)
                    && ((stPTorsion->fFlags&BOUNDARY)!=0) ) {
                    FortranWriteInt( stPTorsion->iPertParmIndex );
                }
            }
        }
        FortranEndLine();

                /* -36G- Residue labels at LAMBDA = 1 */
                        /* Just write the labels at LAMBDA = 0 */

                /* Trim the string down to at most 3 characters by */
                /* taking the last three characters if it is too long */
        FortranDebug( "-36G-" );

        MESSAGE("Writing the residue labels\n" );
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
        MESSAGE("Writing the atomic polarizabilities\n" );
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
                VP0("Total atoms with default polarization=0.0: %d of %d\n",
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
                        BOOL    bTmp = bAtomPerturbed( saPAtom->aAtom );

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
                    VP0("Total pert atoms with default polarization=0.0: %d of %d\n",
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
                FortranWriteDouble( PVAI( vaNBParameters, NONBONDACt, i )->dC14 );
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
#endif
