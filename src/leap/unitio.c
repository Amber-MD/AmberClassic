/*
 *        File:        unitio.c
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
 *                                                                      *
 *     Designed by:    Christian Schafmeister                           *
 *     Author:         Christian Schafmeister                           *
 *                                                                      *
 *     VERSION: 1.0                                                     *
 *     Programmers:                                                     *
 *             Christian Schafmeister                                   *
 *             David Rivkin                                             *
 *                                                                      *
 *     Principal Investigator: Peter A. Kollman                         *
 *                                                                      *
 ************************************************************************
 *
 *        Description:
 *                Input/output routines for UNITs
 *                this has been separated out to make unit.c smaller
 *                and because there is ALOT of complex code here
 *                that doesn't have alot to do with the day to day operation
 *                of UNITs.
 *
 *
 *      Class:
 *              UNIT
 *      Superclass:
 *              CONTAINER, LOOP
 *
 *      Description:
 *              A UNIT is a subclass of CONTAINER.
 *              UNITS can contain molecules, residues, and atoms.
 *              UNITS are used to contain the entire molecular
 *              system plus other information about the system.
 *
 *        NOTE:        OFF files are used to write UNITs to files.
 *                In OFF files there is no implicit ordering
 *                of ANY DATA WHATSOEVER.  The PROGRAMMER is
 *                to assume that there is no order regardless
 *                of the ordering that the code in this file
 *                generates.
 *
 *        NOTE2:        When setting up tables for writing to OFF files,
 *                Indices are FORTRAN indices, where the first
 *                element has index=1.
 *
 */

/*      Modifications induced by the implementation of the savemol2 command
*       Christine Cezard (2007)
*       Universite de Picardie - Jules Verne
*       http://q4md-forcefieldtools.org
*       zbUnitIOIndexBondParameters and zUnitDoAtoms are now "extern functions"
*/

/*
 *      Modifications added for atom-spedific C4 interactions.
 *      Using similar algorithm as adding a bond
 *      Zhen Li (2020)
 *      Michigan State University
 */

#include        "basics.h"
#include        "classes.h"
#include        "restraint.h"
#include        "dictionary.h"
#include        "database.h"
#include        "parmLib.h"
#include        "avl.h"
#include        "defaults.h"
#include        "tools.h"
#include        "unitio.h"

int iFatal;


/*
static void debugtypes(UNIT uUnit, char *str)
{
    int i;
    ATOM aAtom;

    fprintf(stderr, "--------------------- %s\n", str);
    for (i = 0; i < iVarArrayElementCount(uUnit->vaAtoms); i++) {
        aAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i)->aAtom;
        fprintf(stderr, "\t%d\t%s\t%s\n",
                iAtomId(aAtom), sAtomName(aAtom), sAtomType(aAtom));
    }
    fprintf(stderr, "--------------------------------\n");
}
*/

static void
zSetTreeType(ATOM aAtom, int iCoordLeft)
{
    switch (iCoordLeft) {
    case 0:
        AtomSetTempDouble(aAtom, (double) 'E');
        break;
    case 1:
        AtomSetTempDouble(aAtom, (double) 'S');
        break;
    case 2:
        AtomSetTempDouble(aAtom, (double) 'B');
        break;
    case 3:
        AtomSetTempDouble(aAtom, (double) '3');
        break;
    case 4:
        VPWARN("  (%s: unexpected coordination 4 - %s)\n",
             sAtomName(aAtom), "using nonstandard tree type '4'");
        AtomSetTempDouble(aAtom, (double) '4');
        break;
    case 5:
        AtomSetTempDouble(aAtom, (double) '5');
        VPWARN("  (%s: unexpected coordination 5 - %s)\n",
             sAtomName(aAtom), "using nonstandard tree type '5'");
        break;
    case 6:
        AtomSetTempDouble(aAtom, (double) '6');
        VPWARN("  (%s: unexpected coordination 6 - %s)\n",
             sAtomName(aAtom), "using nonstandard tree type '6'");
        break;
    default:
        AtomSetTempDouble(aAtom, (double) 'X');
        VPWARN("  (%s: unexpected coordination %d - %s)\n",
             sAtomName(aAtom), iCoordLeft, "using tree type 'X'");
        break;
    }
}

/*
 *      zUnitLoadTables
 *
 *        Author:        Christian Schafmeister (1991)
 *
 *      Load the tables which can be used to construct the UNIT
 *      from a DATABASE.
 */
BOOL zbUnitIOLoadTables(UNIT uUnit, DATABASE db)
{
    int iSize, iCount, iAtomCount, iType;
    STRING sName;
    SAVEATOMt *saPAtom;
    SAVECONNECTIVITYt *scPConnectivity;
    SAVERESTRAINTt *srPRestraint;
    SAVERESIDUEt *srPResidue, *srPResTemp;
    SAVEMOLECULEt *smPMolecule;
    SAVEHIERARCHYt *shPHierarchy;
    SAVEGROUPSt *sgPGroupAtom;
    SAVEC4Pairwiset *scPC4Pairwise; // New2021
    int iBondCount, iRestraintCount, iSequence;
    BOOL bGotOne;
    SAVEBOXt sbBox;
    SAVECAPt scCap;
    int iTemp, iLen, i;


    DBPushPrefix(db, "unit.");

    if (!bDBGetValue(db, "name", &iLen, (GENP) sName, sizeof(sName))) {
        bGotOne = FALSE;
        goto DONE;
    }
    ContainerSetName(uUnit, sName); // FIXME check this is OK

    if (bDBGetValue(db, "description", &iLen, (GENP) sName, sizeof(sName))) {
        UnitSetDescription(uUnit, sName);
    }

    bDBGetValue(db, "childsequence", &iLen, (GENP) & iSequence, 0);
    ContainerSetNextChildsSequence(uUnit, iSequence);

    /* Construct an array of atoms with names and types */

    uUnit->vaAtoms = vaVarArrayCreate(sizeof(SAVEATOMt));
    bDBGetType(db, "atoms", &iType, &iAtomCount);
    VarArraySetSize((uUnit->vaAtoms), iAtomCount);
    saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, 0);
    iSize = sizeof(SAVEATOMt);
    bDBGetTable(db, "atoms", &iAtomCount,
                3, (char *) &(saPAtom->iTypeIndex), iSize,
                4, (char *) &(saPAtom->iResidueIndex), iSize,
                5, (char *) &(saPAtom->fFlags), iSize,
                6, (char *) &(saPAtom->iSequence), iSize,
                7, (char *) &(saPAtom->iElement), iSize,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                8, (char *) &(saPAtom->dCharge), iSize,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                1, (char *) saPAtom->sName, iSize,
                2, (char *) saPAtom->sType, iSize,
                0, NULL, 0, 0, NULL, 0, 0, NULL, 0);

    bDBGetTable(db, "atomspertinfo", &iAtomCount,
                3, (char *) &(saPAtom->iPertTypeIndex), iSize,
                4, (char *) &(saPAtom->iPertElement), iSize,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                5, (char *) &(saPAtom->dPertCharge), iSize,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                1, (char *) saPAtom->sPertName, iSize,
                2, (char *) saPAtom->sPertType, iSize,
                0, NULL, 0, 0, NULL, 0, 0, NULL, 0);

    /* Get the atom positions */

    bDBGetTable(db, "positions", &iAtomCount,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                1, (char *) &dVX(&(saPAtom->vPos)), sizeof(SAVEATOMt),
                2, (char *) &dVY(&(saPAtom->vPos)), sizeof(SAVEATOMt),
                3, (char *) &dVZ(&(saPAtom->vPos)), sizeof(SAVEATOMt),
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0, 0, NULL, 0, 0, NULL, 0, 0, NULL, 0);

    /* Get the atoms velocities */

    bDBGetTable(db, "velocities", &iAtomCount,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                1, (char *) &dVX(&(saPAtom->vVelocity)), sizeof(SAVEATOMt),
                2, (char *) &dVY(&(saPAtom->vVelocity)), sizeof(SAVEATOMt),
                3, (char *) &dVZ(&(saPAtom->vVelocity)), sizeof(SAVEATOMt),
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0, 0, NULL, 0, 0, NULL, 0, 0, NULL, 0);

    /* Load the bounding box information */

    bDBGetValue(db, "boundbox", &iLen, (GENP) & sbBox,
                sizeof(sbBox.dUseBox));
    UnitSetUseBox(uUnit, (sbBox.dUseBox > 0.0));
    UnitSetBeta(uUnit, sbBox.dBeta);
    UnitSetBox(uUnit, sbBox.dXWidth, sbBox.dYWidth, sbBox.dZWidth);
    ToolSanityCheckBox(uUnit);

    /* Load the cap information */

    bDBGetValue(db, "solventcap", &iLen, (GENP) & scCap,
                sizeof(scCap.dUseCap));
    UnitSetUseSolventCap(uUnit, (scCap.dUseCap > 0.0));
    UnitSetSolventCap(uUnit, scCap.dX, scCap.dY, scCap.dZ, scCap.dRadius);

    /* Load the connectivity information */

    uUnit->vaConnectivity = vaVarArrayCreate(sizeof(SAVECONNECTIVITYt));
    bDBGetType(db, "connectivity", &iType, &iBondCount);
    VarArraySetSize((uUnit->vaConnectivity), iBondCount);
    if (iBondCount) {
        scPConnectivity =
            PVAI(uUnit->vaConnectivity, SAVECONNECTIVITYt, 0);
        iSize = sizeof(SAVECONNECTIVITYt);
        bDBGetTable(db, "connectivity", &iBondCount,
                    1, (char *) &(scPConnectivity->iAtom1), iSize,
                    2, (char *) &(scPConnectivity->iAtom2), iSize,
                    3, (char *) &(scPConnectivity->fFlags), iSize,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0, 0, NULL, 0, 0, NULL, 0, 0, NULL, 0);

        /* Load the restraints information */
    }
    uUnit->vaRestraints = vaVarArrayCreate(sizeof(SAVERESTRAINTt));
    bDBGetType(db, "restraints", &iType, &iRestraintCount);
    VarArraySetSize((uUnit->vaRestraints), iRestraintCount);
    if (iRestraintCount) {
        srPRestraint = PVAI(uUnit->vaRestraints, SAVERESTRAINTt, 0);
        iSize = sizeof(SAVERESTRAINTt);
        bDBGetTable(db, "restraints", &iBondCount,
                    1, (char *) &(srPRestraint->iType), iSize,
                    2, (char *) &(srPRestraint->fFlags), iSize,
                    3, (char *) &(srPRestraint->iAtom1), iSize,
                    4, (char *) &(srPRestraint->iAtom2), iSize,
                    5, (char *) &(srPRestraint->iAtom3), iSize,
                    6, (char *) &(srPRestraint->iAtom4), iSize,
                    0, NULL, 0,
                    0, NULL, 0,
                    7, (char *) &(srPRestraint->dKx), iSize,
                    8, (char *) &(srPRestraint->dX0), iSize,
                    9, (char *) &(srPRestraint->dN), iSize,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0, 0, NULL, 0, 0, NULL, 0, 0, NULL, 0);
    }
    /* Load the UNIT connect atoms */

    uUnit->vaConnect = vaVarArrayCreate(sizeof(int));
    VarArraySetSize((uUnit->vaConnect), 2);
    bDBGetValue(db, "connect", &iLen,
                (GENP) PVAI(uUnit->vaConnect, int, 0), sizeof(int));

    /* Load an array of residues */

    uUnit->vaResidues = vaVarArrayCreate(sizeof(SAVERESIDUEt));
    bDBGetType(db, "residues", &iType, &iCount);
    VarArraySetSize((uUnit->vaResidues), iCount);
    srPResidue = PVAI(uUnit->vaResidues, SAVERESIDUEt, 0);
    iSize = sizeof(SAVERESIDUEt);
    if (MAXCONNECT != 6)
        DFATAL("MAXCONNECT has been changed, update UnitLoadTables");
    bDBGetTable(db, "residues", &iCount,
                2, (char *) &(srPResidue->iSequenceNumber), iSize,
                3, (char *) &(srPResidue->iNextChildSequence), iSize,
                4, (char *) &(srPResidue->iAtomStartIndex), iSize,
                6, (char *) &(srPResidue->iImagingAtomIndex), iSize,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                1, (char *) srPResidue->sName, iSize,
                5, (char *) srPResidue->sResidueType, iSize,
                0, NULL, 0, 0, NULL, 0, 0, NULL, 0);

    /* Load the RESIDUE PDB sequence numbers from an */
    /* extra table to remain compatible with old OFF files */

    if (bDBGetType(db, "residuesPdbSequenceNumber", &iType, &iTemp)) {
        bDBGetValue(db, "residuesPdbSequenceNumber", &iTemp,
                    (GENP) & (srPResidue->iPdbResSeq), iSize);
    } else {

        /* If no sequence numbers are defined then create some */
        srPResTemp = srPResidue;
        for (i = 0; i < iCount; i++) {
            //printf("make pdb seq %d\n",i+1);
            srPResTemp->iPdbResSeq = i + 1;
            srPResTemp++;
        }
    }

    /* Load the RESIDUE PDB chainID numbers from an */
    /* extra table to remain compatible with old OFF files */
    if (bDBGetType(db, "residuesPdbChainID", &iType, &iTemp)) {
        bDBGetValue(db, "residuesPdbChainID", &iTemp,
                    (GENP) & (srPResidue->sChainId), iSize);
    } else {
        /* ChainId not stored in OFF file, set to blank */
        for (i = 0,srPResTemp = srPResidue; i < iCount; i++,srPResTemp++)
            srPResTemp->sChainId[0]=0;
    }

    bDBGetTable(db, "residueconnect", &iCount,
                1, (char *) &(srPResidue->iaConnectIndex[0]), iSize,
                2, (char *) &(srPResidue->iaConnectIndex[1]), iSize,
                3, (char *) &(srPResidue->iaConnectIndex[2]), iSize,
                4, (char *) &(srPResidue->iaConnectIndex[3]), iSize,
                5, (char *) &(srPResidue->iaConnectIndex[4]), iSize,
                6, (char *) &(srPResidue->iaConnectIndex[5]), iSize,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0,
                0, NULL, 0, 0, NULL, 0, 0, NULL, 0, 0, NULL, 0);

    /* Construct an array of molecules */

    uUnit->vaMolecules = vaVarArrayCreate(sizeof(SAVEMOLECULEt));
    bDBGetType(db, "molecules", &iType, &iCount);
    VarArraySetSize((uUnit->vaMolecules), iCount);
    if (iCount) {
        smPMolecule = PVAI(uUnit->vaMolecules, SAVEMOLECULEt, 0);
        iSize = sizeof(SAVEMOLECULEt);
        bDBGetTable(db, "molecules", &iCount,
                    2, (char *) &(smPMolecule->iSequenceNumber), iSize,
                    3, (char *) &(smPMolecule->iNextChildSequence), iSize,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    1, (char *) smPMolecule->sName, iSize,
                    0, NULL, 0, 0, NULL, 0, 0, NULL, 0, 0, NULL, 0);
    }
    /* Construct an array which describes the hierarchy */

    uUnit->vaHierarchy = vaVarArrayCreate(sizeof(SAVEHIERARCHYt));
    bDBGetType(db, "hierarchy", &iType, &iCount);
    VarArraySetSize((uUnit->vaHierarchy), iCount);
    if (iCount) {
        shPHierarchy = PVAI(uUnit->vaHierarchy, SAVEHIERARCHYt, 0);
        iSize = sizeof(SAVEHIERARCHYt);
        bDBGetTable(db, "hierarchy", &iCount,
                    2, (char *) &(shPHierarchy->iAboveIndex), iSize,
                    4, (char *) &(shPHierarchy->iBelowIndex), iSize,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    0, NULL, 0,
                    1, (char *) shPHierarchy->sAboveType, iSize,
                    3, (char *) shPHierarchy->sBelowType, iSize,
                    0, NULL, 0, 0, NULL, 0, 0, NULL, 0);
    }

    /* New Load C4Pairwise interactions */
    uUnit->vaC4Pairwise = vaVarArrayCreate(sizeof(SAVEC4Pairwiset));
    bDBGetType(db, "C4Pairwise", &iType, &iCount);
    VarArraySetSize((uUnit->vaC4Pairwise), iCount);
    if (iCount) {
        scPC4Pairwise = PVAI(uUnit->vaC4Pairwise, SAVEC4Pairwiset, 0);
        iSize = sizeof(SAVEC4Pairwiset);
        bDBGetTable(db, "C4Pairwise", &iCount,
                   1, (char *) &(scPC4Pairwise->iAtom1), iSize,
                   2, (char *) &(scPC4Pairwise->iAtom2), iSize,
                   3, (char *) &(scPC4Pairwise->iParmIndex), iSize,
                   4, (char *) &(scPC4Pairwise->daC4Pairwise), iSize,
                   0, NULL, 0,
                   0, NULL, 0,
                   0, NULL, 0,
                   0, NULL, 0,
                   0, NULL, 0,
                   0, NULL, 0,
                   0, NULL, 0,
                   0, NULL, 0,
                   0, NULL, 0,
                   0, NULL, 0,
                   0, NULL, 0,
                   0, NULL, 0,
                   0, NULL, 0);
    }


    /* Load the atom groups */

    if (bDBGetType(db, "groupNames", &iType, &iCount)) {
        uUnit->vaGroupNames = vaVarArrayCreate(sizeof(STRING));
        uUnit->vaGroupAtoms = vaVarArrayCreate(sizeof(SAVEGROUPSt));
        VarArraySetSize((uUnit->vaGroupNames), iCount);
        if (iCount) {
            bDBGetValue(db, "groupNames", &iCount,
                        (GENP) PVAI(uUnit->vaGroupNames, STRING, 0),
                        sizeof(STRING));
            if (bDBGetType(db, "groupAtoms", &iType, &iCount)) {
                VarArraySetSize(uUnit->vaGroupAtoms, iCount);
                if (iCount) {
                    sgPGroupAtom = PVAI(uUnit->vaGroupAtoms, SAVEGROUPSt,
                                        0);
                    iSize = sizeof(SAVEGROUPSt);
                    bDBGetTable(db, "groupAtoms", &iCount,
                                1, (char *) &(sgPGroupAtom->iGroupIndex),
                                iSize, 2,
                                (char *) &(sgPGroupAtom->iIndexAtom),
                                iSize, 0, NULL, 0, 0, NULL, 0, 0, NULL, 0,
                                0, NULL, 0, 0, NULL, 0, 0, NULL, 0, 0,
                                NULL, 0, 0, NULL, 0, 0, NULL, 0, 0, NULL,
                                0, 0, NULL, 0, 0, NULL, 0, 0, NULL, 0, 0,
                                NULL, 0, 0, NULL, 0);
                }
            }
        }
    }

    /*
     *  DELETED: uUnit->psParameters = psParmSetLoad( db );
     *  since saved params confuse the issue of precedence.
     */

    bGotOne = TRUE;

  DONE:
    DBPopPrefix(db);
    return (bGotOne);
}

/*
 *      zUnitIOSaveTables
 *
 *        Author:        Christian Schafmeister (1991)
 *
 *      Save the tables which can be used to construct the UNIT
 *      from a DATABASE.
 */
void zUnitIOSaveTables(UNIT uUnit, DATABASE db)
{
    int iSize, iCount, iSequence;
    SAVEATOMt *saPAtom;
    SAVEBONDt *sbPBond;
    SAVEC4Pairwiset *scPC4Pairwise; //New
    SAVECONNECTIVITYt *scPConnectivity;
    SAVEANGLEt *saPAngle;
    SAVETORSIONt *stPTorsion;
    SAVERESTRAINTt *srPRestraint;
    SAVERESIDUEt *srPResidue;
    SAVEMOLECULEt *smPMolecule;
    SAVEHIERARCHYt *shPHierarchy;
    SAVEGROUPSt *sgPGroupAtom;
    SAVEBOXt sbBox;
    SAVECAPt scCap;

    if (!iVarArrayElementCount(uUnit->vaAtoms)) {
        VPFATALEXIT("\tUnit has no atoms!\n");
        return;
    }
    DBPushPrefix(db, "unit.");

    DBPutValue(db, "name", ENTRYSINGLE | ENTRYSTRING, 1,
               (GENP) sContainerName(uUnit), 0);

    if (sUnitDescription(uUnit)[0]!=0) {
        DBPutValue(db, "description", ENTRYSINGLE | ENTRYSTRING, 1,
               (GENP) sUnitDescription(uUnit), 0);
    }

    iSequence = iContainerNextChildsSequence(uUnit);
    DBPutValue(db, "childsequence", ENTRYSINGLE | ENTRYINTEGER, 1,
               (GENP) & iSequence, 0);

    /* Save the array of atoms with names and types */

    iCount = iVarArrayElementCount(uUnit->vaAtoms);
    saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, 0);
    iSize = sizeof(SAVEATOMt);
    DBPutTable(db, "atoms", iCount,
               3, "typex", (char *) &(saPAtom->iTypeIndex), iSize,
               4, "resx", (char *) &(saPAtom->iResidueIndex), iSize,
               5, "flags", (char *) &(saPAtom->fFlags), iSize,
               6, "seq", (char *) &(saPAtom->iSequence), iSize,
               7, "elmnt", (char *) &(saPAtom->iElement), iSize,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               8, "chg", (char *) &(saPAtom->dCharge), iSize,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               1, "name", (char *) saPAtom->sName, iSize,
               2, "type", (char *) saPAtom->sType, iSize,
               0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0);

    DBPutTable(db, "atomspertinfo", iCount,
               3, "ptypex", (char *) &(saPAtom->iPertTypeIndex), iSize,
               4, "pelmnt", (char *) &(saPAtom->iPertElement), iSize,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               5, "pchg", (char *) &(saPAtom->dPertCharge), iSize,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               1, "pname", (char *) saPAtom->sPertName, iSize,
               2, "ptype", (char *) saPAtom->sPertType, iSize,
               0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0);

    /* Get the atom positions */

    DBPutTable(db, "positions", iCount,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               1, "x", (char *) &dVX(&(saPAtom->vPos)), sizeof(SAVEATOMt),
               2, "y", (char *) &dVY(&(saPAtom->vPos)), sizeof(SAVEATOMt),
               3, "z", (char *) &dVZ(&(saPAtom->vPos)), sizeof(SAVEATOMt),
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0);

    /* Get the atom velocities */

    DBPutTable(db, "velocities", iCount,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               0, NULL, NULL, 0,
               1, "x", (char *) &dVX(&(saPAtom->vVelocity)),
               sizeof(SAVEATOMt), 2, "y",
               (char *) &dVY(&(saPAtom->vVelocity)), sizeof(SAVEATOMt), 3,
               "z", (char *) &dVZ(&(saPAtom->vVelocity)),
               sizeof(SAVEATOMt), 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0,
               NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
               NULL, 0);

    /* Save the bounding box information */

    ToolSanityCheckBox(uUnit);
    sbBox.dBeta = dUnitBeta(uUnit);
    UnitGetBox(uUnit, &(sbBox.dXWidth),
               &(sbBox.dYWidth), &(sbBox.dZWidth));
    if (bUnitUseBox(uUnit))
        sbBox.dUseBox = 1.0;
    else
        sbBox.dUseBox = -1.0;
    DBPutValue(db, "boundbox", ENTRYDOUBLE | ENTRYARRAY,
               sizeof(sbBox) / sizeof(sbBox.dUseBox),        /* # of doubles */
               (GENP) & sbBox, sizeof(sbBox.dUseBox));

    /* Save the cap information */

    if (bUnitUseSolventCap(uUnit))
        scCap.dUseCap = 1.0;
    else
        scCap.dUseCap = -1.0;
    UnitGetSolventCap(uUnit, &(scCap.dX), &(scCap.dY), &(scCap.dZ),
                      &(scCap.dRadius));
    DBPutValue(db, "solventcap", ENTRYDOUBLE | ENTRYARRAY,
               sizeof(scCap) / sizeof(scCap.dUseCap),
               (GENP) & scCap, sizeof(scCap.dUseCap));

    /* Save the connectivity information */

    if ((iCount = iVarArrayElementCount(uUnit->vaConnectivity))) {
        scPConnectivity =
            PVAI(uUnit->vaConnectivity, SAVECONNECTIVITYt, 0);
        iSize = sizeof(SAVECONNECTIVITYt);
        DBPutTable(db, "connectivity", iCount,
                   1, "atom1x", (char *) &(scPConnectivity->iAtom1), iSize,
                   2, "atom2x", (char *) &(scPConnectivity->iAtom2), iSize,
                   3, "flags", (char *) &(scPConnectivity->fFlags), iSize,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0);
    }
    /* Save the restraints information */

    if ((iCount = iVarArrayElementCount(uUnit->vaRestraints))) {
        srPRestraint = PVAI(uUnit->vaRestraints, SAVERESTRAINTt, 0);
        iSize = sizeof(SAVERESTRAINTt);
        DBPutTable(db, "restraints", iCount,
                   1, "type", (char *) &(srPRestraint->iType), iSize,
                   2, "flags", (char *) &(srPRestraint->fFlags), iSize,
                   3, "atom1x", (char *) &(srPRestraint->iAtom1), iSize,
                   4, "atom2x", (char *) &(srPRestraint->iAtom2), iSize,
                   5, "atom3x", (char *) &(srPRestraint->iAtom3), iSize,
                   6, "atom4x", (char *) &(srPRestraint->iAtom4), iSize,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0,
                   7, "kx", (char *) &(srPRestraint->dKx), iSize,
                   8, "x0", (char *) &(srPRestraint->dX0), iSize,
                   9, "n", (char *) &(srPRestraint->dN), iSize,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0,
                   0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0);
    }
    /* Save the UNIT connect atoms */

    DBPutValue(db, "connect", ENTRYARRAY | ENTRYINTEGER, 2,
               (GENP) PVAI(uUnit->vaConnect, int, 0), sizeof(int));

    /* Save the array of residues */

    if ((iCount = iVarArrayElementCount(uUnit->vaResidues))) {
        srPResidue = PVAI(uUnit->vaResidues, SAVERESIDUEt, 0);
        iSize = sizeof(SAVERESIDUEt);
        if (MAXCONNECT != 6)
            DFATAL("MAXCONNECT has been changed, update UnitSaveTables");

        DBPutTable(db, "residues", iCount,
                   2, "seq", (char *) &(srPResidue->iSequenceNumber),
                   iSize, 3, "childseq",
                   (char *) &(srPResidue->iNextChildSequence), iSize, 4,
                   "startatomx", (char *) &(srPResidue->iAtomStartIndex),
                   iSize, 6, "imagingx",
                   (char *) &(srPResidue->iImagingAtomIndex), iSize, 0,
                   NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0,
                   NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0,
                   NULL, NULL, 0, 0, NULL, NULL, 0, 1, "name",
                   (char *) srPResidue->sName, iSize, 5, "restype",
                   (char *) srPResidue->sResidueType, iSize, 0, NULL, NULL,
                   0, 0, NULL, NULL, 0, 0, NULL, NULL, 0);

        if (iVarArrayElementCount(uUnit->vaConnect)) {
            DBPutValue(db, "connect", ENTRYARRAY | ENTRYINTEGER, 2,
                       (GENP) PVAI(uUnit->vaConnect, int, 0), sizeof(int));
        }
        DBPutValue(db, "residuesPdbSequenceNumber",
                   ENTRYARRAY | ENTRYINTEGER, iCount,
                   (GENP) & srPResidue->iPdbResSeq, iSize);

        DBPutTable(db, "residueconnect", iCount,
                   1, "c1x", (char *) &(srPResidue->iaConnectIndex[0]),
                   iSize, 2, "c2x",
                   (char *) &(srPResidue->iaConnectIndex[1]), iSize, 3,
                   "c3x", (char *) &(srPResidue->iaConnectIndex[2]), iSize,
                   4, "c4x", (char *) &(srPResidue->iaConnectIndex[3]),
                   iSize, 5, "c5x",
                   (char *) &(srPResidue->iaConnectIndex[4]), iSize, 6,
                   "c6x", (char *) &(srPResidue->iaConnectIndex[5]), iSize,
                   0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0,
                   NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0,
                   NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0,
                   NULL, NULL, 0, 0, NULL, NULL, 0);
    }
    /* Save an array of molecules */

    if ((iCount = iVarArrayElementCount(uUnit->vaMolecules))) {
        smPMolecule = PVAI(uUnit->vaMolecules, SAVEMOLECULEt, 0);
        iSize = sizeof(SAVEMOLECULEt);
        DBPutTable(db, "molecules", iCount,
                   2, "seqnum", (char *) &(smPMolecule->iSequenceNumber),
                   iSize, 3, "childseq",
                   (char *) &(smPMolecule->iNextChildSequence), iSize, 0,
                   NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0,
                   NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0,
                   NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0,
                   NULL, NULL, 0, 1, "name", (char *) smPMolecule->sName,
                   iSize, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0);
    }
    /* Save an array which describes the hierarchy */

    if ((iCount = iVarArrayElementCount(uUnit->vaHierarchy))) {
        shPHierarchy = PVAI(uUnit->vaHierarchy, SAVEHIERARCHYt, 0);
        iSize = sizeof(SAVEHIERARCHYt);
        DBPutTable(db, "hierarchy", iCount,
                   2, "abovex", (char *) &(shPHierarchy->iAboveIndex),
                   iSize, 4, "belowx",
                   (char *) &(shPHierarchy->iBelowIndex), iSize, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 1, "abovetype",
                   (char *) shPHierarchy->sAboveType, iSize, 3,
                   "belowtype", (char *) shPHierarchy->sBelowType, iSize,
                   0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0);
    }
    /* Save the group information */

    if (uUnit->vaGroupNames) {
        DBPutValue(db, "groupNames", ENTRYSTRING | ENTRYARRAY,
                   iVarArrayElementCount(uUnit->vaGroupNames),
                   (GENP) PVAI(uUnit->vaGroupNames, STRING, 0),
                   iVarArrayElementSize(uUnit->vaGroupNames));
        sgPGroupAtom = PVAI(uUnit->vaGroupAtoms, SAVEGROUPSt, 0);
        iSize = iVarArrayElementSize(uUnit->vaGroupAtoms);
        DBPutTable(db, "groupAtoms",
                   iVarArrayElementCount(uUnit->vaGroupAtoms),
                   1, "groupIndex", (char *) &(sgPGroupAtom->iGroupIndex),
                   iSize, 2, "atomIndex",
                   (char *) &(sgPGroupAtom->iIndexAtom), iSize, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0);
    }

    /*
     *  DELETED: if ( uUnit->psParameters != NULL )
     *               ParmSetSave( uUnit->psParameters, db );
     *  since saving these confuses the precedence of params
     */

    /* Save the bond information */

    if ((iCount = iVarArrayElementCount(uUnit->vaBonds))) {
        sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, 0);
        iSize = sizeof(SAVEBONDt);
        DBPutTable(db, "bonds", iCount,
                   1, "atom1x", (char *) &(sbPBond->iAtom1), iSize,
                   2, "atom2x", (char *) &(sbPBond->iAtom2), iSize,
                   3, "parmx", (char *) &(sbPBond->iParmIndex), iSize,
                   4, "pertparmx", (char *) &(sbPBond->iPertParmIndex),
                   iSize, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0);
    }

    // New save the C4 information

    if ((iCount = iVarArrayElementCount(uUnit->vaC4Pairwise))) {
        scPC4Pairwise = PVAI(uUnit->vaC4Pairwise, SAVEC4Pairwiset, 0);
        iSize = sizeof(SAVEC4Pairwiset);
        DBPutTable(db, "C4Pairwise", iCount,
                   1, "atom1x", (char *) &(scPC4Pairwise->iAtom1), iSize,
                   2, "atom2x", (char *) &(scPC4Pairwise->iAtom2), iSize,
                   3, "parmx", (char *) &(scPC4Pairwise->iParmIndex), iSize,
                   4, "daC4Pairwise", (char *) &(scPC4Pairwise->daC4Pairwise), iSize,
                   0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0);
    }


    /* Save the angles information */

    if ((iCount = iVarArrayElementCount(uUnit->vaAngles))) {

        saPAngle = PVAI(uUnit->vaAngles, SAVEANGLEt, 0);
        iSize = sizeof(SAVEANGLEt);
        DBPutTable(db, "angles", iCount,
                   1, "atom1x", (char *) &(saPAngle->iAtom1), iSize,
                   2, "atom2x", (char *) &(saPAngle->iAtom2), iSize,
                   3, "atom3x", (char *) &(saPAngle->iAtom3), iSize,
                   4, "parmx", (char *) &(saPAngle->iParmIndex), iSize,
                   5, "pertparmx", (char *) &(saPAngle->iPertParmIndex),
                   iSize, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0);
    }


    /* Save the torsions information */

    if ((iCount = iVarArrayElementCount(uUnit->vaTorsions))) {

        stPTorsion = PVAI(uUnit->vaTorsions, SAVETORSIONt, 0);
        iSize = sizeof(SAVETORSIONt);
        DBPutTable(db, "torsions", iCount,
                   1, "atom1x", (char *) &(stPTorsion->iAtom1), iSize,
                   2, "atom2x", (char *) &(stPTorsion->iAtom2), iSize,
                   3, "atom3x", (char *) &(stPTorsion->iAtom3), iSize,
                   4, "atom4x", (char *) &(stPTorsion->iAtom4), iSize,
                   5, "parmx", (char *) &(stPTorsion->iParmIndex), iSize,
                   6, "pertparmx", (char *) &(stPTorsion->iPertParmIndex),
                   iSize, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0, 0, NULL,
                   NULL, 0, 0, NULL, NULL, 0, 0, NULL, NULL, 0);
    }

    DBPopPrefix(db);
}





/*
 *      zUnitIOTableAddAtom
 *
 *        Author:        Christian Schafmeister (1991)
 *
 *      Add an atom to the UNITs tables.
 *        Return FALSE in (*bPFailed) if the atom could not be added,
 *        this can happen when the type is unknown.
 */
static void
zUnitIOTableAddAtom(UNIT uUnit, ATOM aAtom, int i, PARMLIB plParameters,
                    BOOL * bPFailed, BOOL bPert)
{
    SAVEATOMt *saPAtom;
    int iIndex, iElement, iHybridization, iTemp;
    double dMass, dPolar, dDepth, dRStar, dDepth14, dRStar14, dScreenF;
    STRING sType, sDesc, sTemp;
    PARMSET psTemp;

    *bPFailed = FALSE;

    /* Define the ATOM index number in the SAVEATOMt array */

    ContainerSetTempInt(aAtom, i + 1);
    saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, i);


    saPAtom->aAtom = aAtom;
    REF(aAtom);

    strcpy(saPAtom->sName, sContainerName(aAtom));
    if (strlen(sAtomPertName(aAtom)) != 0) {
        strcpy(saPAtom->sPertName, sAtomPertName(aAtom));
    } else {
        MESSAGE(" no pert atom name set for %s - using nonpert\n",
                 sContainerFullDescriptor((CONTAINER) aAtom, sTemp));
        strcpy(saPAtom->sPertName, sContainerName(aAtom));
    }
    strcpy(saPAtom->sType, sAtomType(aAtom));
    if (strlen(sAtomPertType(aAtom)) != 0) {
        strcpy(saPAtom->sPertType, sAtomPertType(aAtom));
    } else {
        MESSAGE(" no pert atom type set for %s - using nonpert\n",
                 sContainerFullDescriptor((CONTAINER) aAtom, sTemp));
        strcpy(saPAtom->sPertType, sAtomType(aAtom));
    }
    saPAtom->dCharge = dAtomCharge(aAtom);
    saPAtom->dPertCharge = dAtomPertCharge(aAtom);
    saPAtom->iSequence = iContainerSequence(aAtom);
    saPAtom->iElement = iAtomElement(aAtom);
    saPAtom->iPertElement = iAtomPertElement(aAtom);

    /* If the atom is contained by a residue then set */
    /* which residue it is, otherwise set it to zero */

    if (iObjectType(cContainerWithin(aAtom)) == RESIDUEid)
        saPAtom->iResidueIndex =
            iContainerTempInt(cContainerWithin(aAtom));
    else
        saPAtom->iResidueIndex = 0;
    saPAtom->vPos = vAtomPosition(aAtom);
    saPAtom->vVelocity = vAtomVelocity(aAtom);
    saPAtom->fFlags = fAtomFlags(aAtom);
    saPAtom->iTypeIndex = 0;
    saPAtom->iPertTypeIndex = 0;

    /* If parameters should be generated then lookup the atom */
    /* in the PARMLIB, and add it to the UNITs PARMSET */
    /* then point the atoms TypeIndex and PertTypeIndex's to */
    /* the entry in the PARMSET */

    if (strlen(sAtomType(aAtom)) == 0) {
        /*
         *  may not matter, since UnitCheck catches
         *      this for saving parm, & not a problem
         *      for saveoff
         */
        /* try to mess things up.. bPFailed seems to be ignored */
        iIndex = -9999;
    } else if (plParameters != NULL) {
        iIndex = iParmSetFindAtom(uUnit->psParameters, sAtomType(aAtom));
        if (iIndex == PARM_NOT_FOUND) {

            PARMLIB_LOOP(plParameters, psTemp,
                         (iTemp = iParmSetFindAtom(psTemp,
                                                   sAtomType(aAtom))));
            if (iTemp != PARM_NOT_FOUND) {
                ParmSetAtom(psTemp, iTemp, sType,
                            &dMass, &dPolar, &dDepth, &dRStar,
			    &dDepth14, &dRStar14, &dScreenF, &iElement,
			    &iHybridization, sDesc);
                iIndex =
                    iParmSetAddAtom(uUnit->psParameters, sAtomType(aAtom),
                                    dMass, dPolar, dDepth, dRStar,
                                    dDepth14, dRStar14, dScreenF, iElement,
                                    iHybridization, sDesc);
            } else {
                iIndex = 0;
                VPERROR("For atom (%s) could not find vdW (or other) "
                        "parameters for type (%s)\n",
                     sContainerFullDescriptor((CONTAINER) aAtom, sTemp),
                     sAtomType(aAtom));
                *bPFailed = TRUE;
            }
        } else {
            ParmSetAtom(uUnit->psParameters, iIndex, sType,
                        &dMass, &dPolar, &dDepth, &dRStar, &dDepth14,
                        &dRStar14, &dScreenF, &iElement, &iHybridization, sDesc);
        }
        saPAtom->iTypeIndex = iIndex + 1;
        if(iElement != NOELEMENT)
            saPAtom->iElement = iElement;

        if (bPert && bAtomFlagsSet(aAtom, ATOMPERTURB)) {
            iIndex = iParmSetFindAtom(uUnit->psParameters,
                                      saPAtom->sPertType);
            if (iIndex == PARM_NOT_FOUND) {
                iTemp = PARM_NOT_FOUND;
                PARMLIB_LOOP(plParameters, psTemp,
                             (iTemp = iParmSetFindAtom(psTemp,
                                                       saPAtom->
                                                       sPertType)));
                if (iTemp != PARM_NOT_FOUND) {
                    ParmSetAtom(psTemp, iTemp, sType,
                                &dMass, &dPolar, &dDepth, &dRStar,
                                &dDepth14, &dRStar14, &dScreenF, &iElement,
                                &iHybridization, sDesc);
                    iIndex =
                        iParmSetAddAtom(uUnit->psParameters,
                                        saPAtom->sPertType, dMass, dPolar,
                                        dDepth, dRStar, dDepth14, dRStar14,
					dScreenF,
                                        iElement, iHybridization, sDesc);
                } else {
                    iIndex = 0;
                    VPERROR("For atom (%s) of type (%s) could not find "
                        "perturbed type (%s)\n",
                        sContainerFullDescriptor((CONTAINER) aAtom, sTemp),
                        sType, saPAtom->sPertType );
                    *bPFailed = TRUE;
                }
            } else {
                ParmSetAtom(uUnit->psParameters, iIndex, sType,
                            &dMass, &dPolar, &dDepth, &dRStar, &dDepth14,
                            &dRStar14, &dScreenF, &iElement, &iHybridization,
			    sDesc);

            }
        }
        if (bPert)
            saPAtom->iPertTypeIndex = iIndex + 1;
    }

}


/*
 *  avl tree related declarations - the 1-4 ones are here
 *        for sharing between routines, the impropers to keep
 *        the 1-4's company.
 */
static IX_REC e14, *ePImp = NULL;
static IX_DESC scr14_index, improper_index;
static int *Pint1, *Pint2;


/*
 *      zUnitIOSetCalc14Flags
 *
 *        Author:        Christian Schafmeister (1991)
 *
 *      This routine sets the flags that determine whether
 *      or not 1-4 interactions are calculated for the two
 *      atoms across the torsion.  It does this
 *      by checking the bond tables, and angle tables to see
 *      if the two outer atoms in the torsion are represented
 *      as having a bond or bond angle between them, if this is
 *      TRUE then there is no 1-4 interaction.  Also check
 *      to make sure that no 1-4 interaction has already been
 *      set to be calculated.
 *
 *        The *bPCalc14, *bPCalcPert14 flags are used to determine
 *        if the Calc14 flag has already been set for this torsion.
 *        For each torsion, the first time this routine is called,
 *        both must be set to TRUE.  The first term that has an
 *        interaction will have it's Calc14 flag set.
 *
 */
static void
zUnitIOSetCalc14Flags(SAVETORSIONt * stPTorsion, BOOL * bPCalc14,
                      BOOL * bPCalcPert14)
{
    BOOL bCheck;

    /*
     *  order the 1st, last pointers by atom #
     */
    if (stPTorsion->iAtom1 < stPTorsion->iAtom4) {
        *Pint1 = stPTorsion->iAtom1;
        *Pint2 = stPTorsion->iAtom4;
    } else {
        *Pint1 = stPTorsion->iAtom4;
        *Pint2 = stPTorsion->iAtom1;
    }

    stPTorsion->bCalc14 = FALSE;
    stPTorsion->bPertCalc14 = FALSE;

    /* Check if we have to check torsions to see if the 1-4 has been set */

    bCheck = FALSE;
    if (stPTorsion->iParmIndex != 0 && *bPCalc14)
        bCheck = TRUE;
    if (stPTorsion->iPertParmIndex != 0 && *bPCalcPert14)
        bCheck = TRUE;

    if (!bCheck)
        return;

    switch (find_key(&e14, &scr14_index,1)) {
    case IX_OK:
        if (e14.recptr != NULL) {
            /*
             *  bond or angle
             */
            return;
        }
        /*
         *  duplicates a previous torsion's 1-4
         */
        if (stPTorsion->iParmIndex != 0)
            *bPCalc14 = FALSE;

        if (stPTorsion->iPertParmIndex != 0)
            *bPCalcPert14 = FALSE;
        break;
    case IX_FAIL:

        /*
         *  a new torsion - add to index
         */
        e14.recptr = (void *) NULL;
        if (!add_key(&e14, &scr14_index))
            DFATAL("1-4 torsions: cannot add key %d %d\n",
                    *Pint1, *Pint2);
        break;
    default:
        DFATAL("unexpected index error\n");
    }
    if (stPTorsion->iParmIndex != 0) {
        stPTorsion->bCalc14 = *bPCalc14;
        *bPCalc14 = FALSE;
    }
    if (stPTorsion->iPertParmIndex != 0) {
        stPTorsion->bPertCalc14 = *bPCalcPert14;
        *bPCalcPert14 = FALSE;
    }

}




/*
 *        zUnitIndexBondParameters
 *
 *        Author:        Christian Schafmeister (1991)
 *
 *        For all of the bonds, search for their parameters within the
 *        UNITs PARMSET.  If they are found then set the index to the entry,
 *        otherwise add the bond parameter to the PARMSET and set the index.
 *        If either of the ATOMs in the bond are to be perturbed then
 *        do the same with the perturbation parameter.
 *
 *        Return TRUE if there was a problem generating parameters.
 */
BOOL
zbUnitIOIndexBondParameters(PARMLIB plLib, UNIT uUnit, BOOL bPert)
{
    int iCount, iIndex, iTemp;
    LOOP lTemp;
    SAVEBONDt *sbPBond;
    ATOM aAtom1, aAtom2;
    BOOL bFailedGeneratingParameters;
    double dKb, dR0, dKpull, dRpull0, dKpress, dRpress0;
    STRING sAtom1, sAtom2, sDesc;
    PARMSET psTemp;
#ifdef  DEBUG
    STRING sTemp1, sTemp2;
#endif


    VPTRACEENTER("zbUnitIOIndexBondParameters" );
    bFailedGeneratingParameters = FALSE;

    if (uUnit->vaBonds != NULL) {
        VP0("Rebuilding bond parameters.\n");
        VarArrayDestroy(&(uUnit->vaBonds));
    } else
        VP0("Building bond parameters.\n");

    uUnit->vaBonds = vaVarArrayCreate(sizeof(SAVEBONDt));
    iCount = 0;
    lTemp = lLoop((OBJEKT) uUnit, BONDS);
    while (oNext(&lTemp) != NULL)
        iCount++;
    VarArraySetSize((uUnit->vaBonds), iCount);
    if (iCount) {
        lTemp = lLoop((OBJEKT) uUnit, BONDS);
        sbPBond = PVAI(uUnit->vaBonds, SAVEBONDt, 0);
        if (sbPBond == NULL)
            DFATAL(" ?? null\n");
        for (; oNext(&lTemp) != NULL; sbPBond++) {
            LoopGetBond(&lTemp, &aAtom1, &aAtom2);
            sbPBond->iAtom1 = iContainerTempInt(aAtom1);
            sbPBond->iAtom2 = iContainerTempInt(aAtom2);
            sbPBond->iParmIndex = 0;
            sbPBond->iPertParmIndex = 0;
            sbPBond->fFlags = 0;
            strcpy(sAtom1, sAtomType(aAtom1));
            strcpy(sAtom2, sAtomType(aAtom2));
            VPTRACE("Searching for bond parameter for atom type %s with name %s,\n"
                    "       atomic number %i, Id %i, and Index %i \n"
                    "       at position %8.3f,%8.3f,%8.3f \n",
                    sAtom1, sContainerName( aAtom1), iAtomElement(aAtom1),
                    iAtomId(aAtom1), iAtomIndex(aAtom1), vAtomPosition(aAtom1).dX,
                    vAtomPosition(aAtom1).dY, vAtomPosition(aAtom1).dZ );
            iIndex = iParmSetFindBond(uUnit->psParameters, sAtom1, sAtom2);
            if (iIndex == PARM_NOT_FOUND) {
                iTemp = PARM_NOT_FOUND;
                PARMLIB_LOOP(plLib, psTemp,
                        (iTemp = iParmSetFindBond(psTemp, sAtom1, sAtom2)));
                if (iTemp != PARM_NOT_FOUND) {
                    ParmSetBond(psTemp, iTemp, sAtom1, sAtom2, &dKb, &dR0,
                            &dKpull, &dRpull0, &dKpress, &dRpress0, sDesc);
                    iIndex = iParmSetAddBond(uUnit->psParameters, sAtom1, sAtom2,
                            dKb, dR0, dKpull, dRpull0, dKpress, dRpress0, sDesc);
                } else {
                    bFailedGeneratingParameters = TRUE;
                    iIndex = 0;
                    VECTOR vPos1 = vAtomPosition(aAtom1);
                    VECTOR vPos2 = vAtomPosition(aAtom2);
                    RESIDUE rRes1 = (RESIDUE)cContainerWithin(aAtom1);
                    RESIDUE rRes2 = (RESIDUE)cContainerWithin(aAtom2);
                    VPERROR("Could not find bond parameter for atom types: %s - %s\n"
                            "        for atom %s.%s:%d.%s at position %8.3f,%8.3f,%8.3f \n"
                            "        and atom %s.%s:%d.%s at position %8.3f,%8.3f,%8.3f.\n",
                            sAtom1, sAtom2,
                            sContainerName(rRes1),sResidueChainId(rRes1),iResiduePdbSequence(rRes1),sContainerName(aAtom1),
                            vPos1.dX, vPos1.dY, vPos1.dZ,
                            sContainerName(rRes2),sResidueChainId(rRes2),iResiduePdbSequence(rRes2),sContainerName(aAtom2),
                            vPos2.dX, vPos2.dY, vPos2.dZ );
                }
            }
            sbPBond->iParmIndex = iIndex + 1;

            if (bPert &&
                (bAtomFlagsSet(aAtom1, ATOMPERTURB) ||
                 bAtomFlagsSet(aAtom2, ATOMPERTURB))) {

                /* Note that the bond is perturbed and whether or */
                /* not it is on the boundary between perturbed and */
                /* non-perturbed */
                MESSAGE("Pert interaction between: %s-%s\n",
                         sContainerFullDescriptor((CONTAINER) aAtom1,
                                                  sTemp1),
                         sContainerFullDescriptor((CONTAINER) aAtom2,
                                                  sTemp2));
                MESSAGE("Indexes = %d-%d\n", sbPBond->iAtom1,
                         sbPBond->iAtom2);
                sbPBond->fFlags |= PERTURBED;
                if (!(bAtomFlagsSet(aAtom1, ATOMPERTURB) &&
                      bAtomFlagsSet(aAtom2, ATOMPERTURB))) {
                    sbPBond->fFlags |= BOUNDARY;
                    MESSAGE("-Boundary\n");
                }

                if (bAtomFlagsSet(aAtom1, ATOMPERTURB)) {
                    if (strlen(sAtomPertType(aAtom1)) != 0)
                        strcpy(sAtom1, sAtomPertType(aAtom1));
                    else
                        strcpy(sAtom1, sAtomType(aAtom1));
                } else
                    strcpy(sAtom1, sAtomType(aAtom1));
                if (bAtomFlagsSet(aAtom2, ATOMPERTURB)) {
                    if (strlen(sAtomPertType(aAtom2)) != 0)
                        strcpy(sAtom2, sAtomPertType(aAtom2));
                    else
                        strcpy(sAtom2, sAtomType(aAtom2));
                } else
                    strcpy(sAtom2, sAtomType(aAtom2));
                iIndex = iParmSetFindBond(uUnit->psParameters,
                                          sAtom1, sAtom2);
                if (iIndex == PARM_NOT_FOUND) {
                    iTemp = PARM_NOT_FOUND;
                    PARMLIB_LOOP(plLib, psTemp,
                                 (iTemp = iParmSetFindBond(psTemp,
                                                           sAtom1,
                                                           sAtom2)));
                    if (iTemp != PARM_NOT_FOUND) {
                        ParmSetBond(psTemp, iTemp, sAtom1, sAtom2,
                                    &dKb, &dR0, &dKpull, &dRpull0, &dKpress, &dRpress0, sDesc);
                        iIndex = iParmSetAddBond(uUnit->psParameters, sAtom1, sAtom2, dKb, dR0,
                                                 dKpull, dRpull0, dKpress, dRpress0, sDesc);
                    } else {
                        bFailedGeneratingParameters = TRUE;
                        iIndex = 0;
                        VPERROR("No bond parameter for: %s - %s\n",
                                sAtom1, sAtom2);
                    }
                }
                sbPBond->iPertParmIndex = iIndex + 1;
            }
        }
    }

    VPTRACEEXIT("zbUnitIOIndexBondParameters" );
    return (bFailedGeneratingParameters);
}

/*
 *        zUnitIndexC4Pairwise
 *
 *        New feature (2021)
 *
 *        For all of the C4 interactions, search for their parameters
 *        within the UNITs PARMSET.  If they are
 *        found then set the index to the entry, otherwise
 *        add the C4 parameter to the PARMSET and set the index.
 *
 *        Return TRUE if there was a problem generating parameters.
 */

BOOL
zbUnitIOIndexC4Pairwise(UNIT uUnit)
{
    int iCount, iIndex;
    LOOP lTemp;
    SAVEC4Pairwiset *scPC4Pairwise;
    ATOM aAtom1, aAtom2;
    double daC4Pairwise; // New
    BOOL bFailedGeneratingParameters;
    STRING sAtom1, sAtom2, sDesc;

    daC4Pairwise = 0.0; // New
    bFailedGeneratingParameters = FALSE;

    if (uUnit->vaC4Pairwise != NULL) {
        // VP0("Rebuilding C4 parameters.\n"); // C4PairwiseDebug
        VarArrayDestroy(&(uUnit->vaC4Pairwise));
    }
    // else
        // VP0("Building C4 parameters.\n"); // C4PairwiseDebug

    uUnit->vaC4Pairwise = vaVarArrayCreate(sizeof(SAVEC4Pairwiset));
    iCount = 0;
    lTemp = lLoop((OBJEKT) uUnit, C4Pairwise);
    while (oNext(&lTemp) != NULL)
        iCount++;
    //VP0("iCount is: %i\n", iCount);
    VarArraySetSize((uUnit->vaC4Pairwise), iCount);
    if (iCount) {
        lTemp = lLoop((OBJEKT) uUnit, C4Pairwise);
        scPC4Pairwise = PVAI(uUnit->vaC4Pairwise, SAVEC4Pairwiset, 0);
        if (scPC4Pairwise == NULL)
            DFATAL(" ?? null\n");
        iIndex = 0; //FIXME why is scPC4Pairwise->iParmIndex set twice below? -- JMK
        for (; oNext(&lTemp) != NULL; scPC4Pairwise++) {
            LoopGetC4Pairwise(&lTemp, &aAtom1, &aAtom2, &daC4Pairwise);
            scPC4Pairwise->iAtom1 = iContainerTempInt(aAtom1);
            scPC4Pairwise->iAtom2 = iContainerTempInt(aAtom2);
            scPC4Pairwise->iParmIndex = 0;
            scPC4Pairwise->daC4Pairwise = daC4Pairwise;
            strcpy(sAtom1, sAtomType(aAtom1));
            strcpy(sAtom2, sAtomType(aAtom2));
            iIndex = iParmSetAddC4Pairwise(uUnit->psParameters, sAtom1, sAtom2, daC4Pairwise, sDesc);
        }
            scPC4Pairwise->iParmIndex = iIndex + 1;
    }
    return (bFailedGeneratingParameters);
}

/*
 *        zUnitIndexAngleParameters
 *
 *        Author:        Christian Schafmeister (1991)
 *
 *        For all of the angles, search for their parameters
 *        within the UNITs PARMSET.  If they are
 *        found then set the index to the entry, otherwise
 *        add the angle parameter to the PARMSET and set the index.
 *        If any of the ATOMs in the angle  are to be perturbed then
 *        do the same with the perturbation parameter.
 *
 *        Return TRUE if there was a problem generating parameters.
 */
static BOOL
zbUnitIOIndexAngleParameters(PARMLIB plLib, UNIT uUnit, BOOL bPert)
{
    LOOP lTemp;
    SAVEANGLEt saAngle;
    ATOM aAtom1, aAtom2, aAtom3;
    int iIndex, iTemp;
    BOOL bFailedGeneratingParameters;
    STRING sAtom1, sAtom2, sAtom3, sDesc;
    PARMSET psTemp;
    double dKt, dT0, dTkub, dRkub;

    bFailedGeneratingParameters = FALSE;

    /* Now generate the ANGLE table */

    if (uUnit->vaAngles != NULL) {
        VP0("Rebuilding angle parameters.\n");
        VarArrayDestroy(&(uUnit->vaAngles));
    } else
        VP0("Building angle parameters.\n");

    uUnit->vaAngles = vaVarArrayCreate(sizeof(SAVEANGLEt));

    lTemp = lLoop((OBJEKT) uUnit, ANGLES);
    while (oNext(&lTemp) != NULL) {
        LoopGetAngle(&lTemp, &aAtom1, &aAtom2, &aAtom3);

        saAngle.iAtom1 = iContainerTempInt(aAtom1);
        saAngle.iAtom2 = iContainerTempInt(aAtom2);
        saAngle.iAtom3 = iContainerTempInt(aAtom3);
        saAngle.iParmIndex = 0;
        saAngle.iPertParmIndex = 0;
        saAngle.fFlags = 0;

        strcpy(sAtom1, sAtomType(aAtom1));
        strcpy(sAtom2, sAtomType(aAtom2));
        strcpy(sAtom3, sAtomType(aAtom3));

/* TODO:Fix this UNGODLY HACK, the PARMSET type has to be modified to */
/* TODO:Allow the user to specify interactions that should be ignored */
/* TODO:Like HW-HW-OW angles for TIP3 waters */

        if (zbUnitIgnoreAngle(sAtom1, sAtom2, sAtom3))
            goto IGNORE1;

        iIndex = iParmSetFindAngle(uUnit->psParameters,
                                   sAtom1, sAtom2, sAtom3);
        if (iIndex == PARM_NOT_FOUND) {
            iTemp = PARM_NOT_FOUND;
            PARMLIB_LOOP(plLib, psTemp,
                         (iTemp = iParmSetFindAngle(psTemp, sAtom1,
                                                    sAtom2, sAtom3)));
            if (iTemp != PARM_NOT_FOUND) {
                ParmSetAngle(psTemp, iTemp,
                             sAtom1, sAtom2, sAtom3,
                             &dKt, &dT0, &dTkub, &dRkub, sDesc);
                iIndex = iParmSetAddAngle(uUnit->psParameters,
                                          sAtom1, sAtom2, sAtom3,
                                          dKt, dT0, dTkub, dRkub, sDesc);
            } else {
                bFailedGeneratingParameters = TRUE;
                iIndex = 0;
                VECTOR vPos1 = vAtomPosition(aAtom1);
                VECTOR vPos2 = vAtomPosition(aAtom2);
                VECTOR vPos3 = vAtomPosition(aAtom3);
                RESIDUE rRes1 = (RESIDUE)cContainerWithin(aAtom1);
                RESIDUE rRes2 = (RESIDUE)cContainerWithin(aAtom2);
                RESIDUE rRes3 = (RESIDUE)cContainerWithin(aAtom3);
                VPERROR("Could not find angle parameter for atom types: %s - %s - %s\n"
                        "        for atom %s.%s:%d.%s at position %8.3f,%8.3f,%8.3f,\n"
                        "            atom %s.%s:%d.%s at position %8.3f,%8.3f,%8.3f,\n"
                        "        and atom %s.%s:%d.%s at position %8.3f,%8.3f,%8.3f.\n",
                        sAtom1, sAtom2, sAtom3,
                        sContainerName(rRes1),sResidueChainId(rRes1),iResiduePdbSequence(rRes1),sContainerName(aAtom1),
                        vPos1.dX, vPos1.dY, vPos1.dZ,
                        sContainerName(rRes2),sResidueChainId(rRes2),iResiduePdbSequence(rRes2),sContainerName(aAtom2),
                        vPos2.dX, vPos2.dY, vPos2.dZ,
                        sContainerName(rRes3),sResidueChainId(rRes3),iResiduePdbSequence(rRes3),sContainerName(aAtom3),
                        vPos3.dX, vPos3.dY, vPos3.dZ );
            }
        }
        saAngle.iParmIndex = iIndex + 1;

      IGNORE1:
        ;

        if (bPert &&
            (bAtomFlagsSet(aAtom1, ATOMPERTURB) ||
             bAtomFlagsSet(aAtom2, ATOMPERTURB) ||
             bAtomFlagsSet(aAtom3, ATOMPERTURB))) {

            /* Note that the angle is perturbed and whether or */
            /* not it is on the boundary between perturbed and */
            /* non-perturbed */
            saAngle.fFlags |= PERTURBED;
            if (!(bAtomFlagsSet(aAtom1, ATOMPERTURB) &&
                  bAtomFlagsSet(aAtom2, ATOMPERTURB) &&
                  bAtomFlagsSet(aAtom3, ATOMPERTURB)))
                saAngle.fFlags |= BOUNDARY;
            if (bAtomFlagsSet(aAtom1, ATOMPERTURB)) {
                if (strlen(sAtomPertType(aAtom1)) != 0)
                    strcpy(sAtom1, sAtomPertType(aAtom1));
                else
                    strcpy(sAtom1, sAtomType(aAtom1));
            } else
                strcpy(sAtom1, sAtomType(aAtom1));
            if (bAtomFlagsSet(aAtom2, ATOMPERTURB)) {
                if (strlen(sAtomPertType(aAtom2)) != 0)
                    strcpy(sAtom2, sAtomPertType(aAtom2));
                else
                    strcpy(sAtom2, sAtomType(aAtom2));
            } else
                strcpy(sAtom2, sAtomType(aAtom2));
            if (bAtomFlagsSet(aAtom3, ATOMPERTURB)) {
                if (strlen(sAtomPertType(aAtom3)) != 0)
                    strcpy(sAtom3, sAtomPertType(aAtom3));
                else
                    strcpy(sAtom2, sAtomType(aAtom2));
            } else
                strcpy(sAtom3, sAtomType(aAtom3));

/* TODO:Fix this UNGODLY HACK, the PARMSET type has to be modified to */
/* TODO:Allow the user to specify interactions that should be ignored */
/* TODO:Like HW-HW-OW angles for TIP3 waters */

            if (zbUnitIgnoreAngle(sAtom1, sAtom2, sAtom3))
                goto IGNORE2;

            iIndex = iParmSetFindAngle(uUnit->psParameters,
                                       sAtom1, sAtom2, sAtom3);
            if (iIndex == PARM_NOT_FOUND) {
                iTemp = PARM_NOT_FOUND;
                PARMLIB_LOOP(plLib, psTemp,
                             (iTemp = iParmSetFindAngle(psTemp, sAtom1,
                                                        sAtom2, sAtom3)));
                if (iTemp != PARM_NOT_FOUND) {
                    ParmSetAngle(psTemp, iTemp,
                                 sAtom1, sAtom2, sAtom3,
                                 &dKt, &dT0, &dTkub, &dRkub, sDesc);
                    iIndex = iParmSetAddAngle(uUnit->psParameters,
                                              sAtom1, sAtom2, sAtom3,
                                              dKt, dT0, dTkub, dRkub,
                                              sDesc);
                } else {
                    bFailedGeneratingParameters = TRUE;
                    iIndex = 0;
                    VPERROR("Can't find angle parameter: %s - %s - %s\n",
                            sAtom1, sAtom2, sAtom3);
                }
            }
            saAngle.iPertParmIndex = iIndex + 1;
        }
      IGNORE2:

        if (saAngle.iPertParmIndex == 0 && saAngle.iParmIndex == 0)
            continue;

        /* Only add the angle interactions if there is a normal interaction */
        /* or a perturbed one */

        VarArrayAdd(uUnit->vaAngles, (GENP) & saAngle);
    }
    return (bFailedGeneratingParameters);
}


/*
 *  BoilTorsions() - reduce list of torsion params to
 *        unique numerical ones, updating param pointers
 *        in the topological list.
 *
 *        Bill Ross, May 1996
 */
static void
BoilTorsions(VARARRAY * vaPParms, int iParmOffset,
             VARARRAY vaTorsions, int iTorsionOffset)
{
    VARARRAY vaB;
    TORSIONPARMt *tpA, *tpB, tC;
    int iIndex, iA, iB, iParmCount, iTorsionCount;

    strcpy(tC.sDesc, "reduced params");
    strcpy(tC.sType1, "__");
    strcpy(tC.sType2, "__");
    strcpy(tC.sType3, "__");
    strcpy(tC.sType4, "__");
    strcpy(tC.sOrder, "___");

    iParmCount = iVarArrayElementCount(*vaPParms);
    iTorsionCount = iVarArrayElementCount(vaTorsions);

    iIndex = iParmOffset;
    vaB = vaVarArrayCreate(sizeof(TORSIONPARMt));
    tpA = PVAI(*vaPParms, TORSIONPARMt, 0);
    for (iA = 0; iA < iParmCount; iA++, tpA++) {
        int i, iOldIndex;
        SAVETORSIONt *stP;

        if (!strcmp(tpA->sType1, "__"))
            continue;
        /*
         *  torsion hasn't been marked as 'superfluous'
         *      so add to new array
         */
        tC.dKp = tpA->dKp;
        tC.iN = tpA->iN;
        tC.dP0 = tpA->dP0;
	tC.dScEE = tpA->dScEE;
	tC.dScNB = tpA->dScNB;
        VarArrayAdd(vaB, (GENP) & tC);
        iIndex++;
        iOldIndex = iParmOffset + iA + 1;

        /*
         *  update any affected torsions
         */
        if (iIndex != iOldIndex) {
            stP = PVAI(vaTorsions, SAVETORSIONt, iTorsionOffset);
            for (i = iTorsionOffset; i < iTorsionCount; i++, stP++) {
                if (stP->iParmIndex == iOldIndex)
                    stP->iParmIndex = iIndex;
                if (stP->iPertParmIndex == iOldIndex)
                    stP->iPertParmIndex = iIndex;
            }
        }

        /*
         *  mark any subsequent duplicates 'superfluous'
         *      and update indexes into the array
         */
        for (tpB = tpA + 1, iB = iA + 1; iB < iParmCount; iB++, tpB++) {
            if (!strcmp(tpB->sType1, "__"))
                continue;
            if (tpB->iN != tpA->iN)
                continue;
            if (tpB->dKp != tpA->dKp)
                continue;
            if (tpB->dP0 != tpA->dP0)
                continue;
	    if (tpB->dScEE != tpA->dScEE)
		continue;
	    if (tpB->dScNB != tpA->dScNB)
		continue;

            /*
             *  B is a duplicate of A
             */
            strcpy(tpB->sType1, "__");
            iOldIndex = iParmOffset + iB + 1;
            stP = PVAI(vaTorsions, SAVETORSIONt, iTorsionOffset);
            for (i = iTorsionOffset; i < iTorsionCount; i++, stP++) {
                if (stP->iParmIndex == iOldIndex)
                    stP->iParmIndex = iIndex;
                if (stP->iPertParmIndex == iOldIndex)
                    stP->iPertParmIndex = iIndex;
            }

        }
    }

    /*
     *  throw away the old parms array & put the new one in place
     */
    VarArrayDestroy(vaPParms);
    *vaPParms = vaB;
}




/*
 *        zbUnitIOIndexTorsionParameters
 *
 *        Author:        Christian Schafmeister (1991)
 *
 *        For all of the angles, search for their parameters
 *        within the UNITs PARMSET.  If they are
 *        found then set the index to the entry, otherwise
 *        add the angle parameter to the PARMSET and set the index.
 *        If any of the ATOMs in the angle  are to be perturbed then
 *        do the same with the perturbation parameter.
 *
 *        (bProper) is TRUE if generating parameters for PROPER torsions,
 *        otherwise IMPROPER torsions.
 *
 *        Return TRUE if there was a problem generating parameters.
 *
 *      Omitting the cacheing function of the UNIT's PARMSET, here
 *      is the intended algorithm:
 *
 *        foreach dihedral in molecule
 *          foreach parmset, working back from the most recently loaded
 *            foreach dihedral in parmset
 *              if an exact match is found
 *                throw away any previous wild card matches
 *                save term (and continue collecting terms in this parmset)
 *              endif
 *              if a wild card match is found and nothing else found yet
 *                  save it
 *              endif
 *            end/parmset_dihedrals
 *            *if an exact match found*
 *              don't search any more parmsets
 *            endif
 *          end/parmsets
 *          ..plug params in for molecule
 *        end/molecule_dihedrals
 */
extern int itest;
static BOOL
zbUnitIOIndexTorsionParameters(PARMLIB plLib, UNIT uUnit,
                               BOOL bProper, BOOL bPert )
{
    VPTRACEENTER("zbUnitIOIndexTorsionParameters" );
    LOOP lTemp;
    SAVETORSIONt stTorsion;
    ATOM aAtom1, aAtom2, aAtom3, aAtom4;
    STRING sAtom1, sAtom2, sAtom3, sAtom4;
    STRING sPert1, sPert2, sPert3, sPert4;
    STRING sOrigAtom1, sOrigAtom2, sOrigAtom3, sOrigAtom4;
    STRING sOrigPert1, sOrigPert2, sOrigPert3, sOrigPert4;
    TORSION tTorsion, tPertTorsion;
    BOOL bPerturbTorsion;
    PARMSET psTemp;
    int iTerm, iPertTerm;
    BOOL bDone, bUse, bUsePert, bCopy, bCopyPert, bEnd, bPertEnd;
    int iN, iPertIndex, iPertN, iLastN, iLastPertN;
    double dKp, dP0, dPertKp, dPertP0;
    double dScEE, dScNB, dPScEE, dPScNB;
    BOOL bCalc14, bCalcPert14;
#ifdef  DEBUG2
    STRING s1, s2, s3, s4;
    int iTParm, iTmp;
    double dTK, dTP;
    STRING sT1, sT2, sT3, sT4, sTemp;
#endif
        STRING sDesc;
    BOOL bFailedGeneratingParameters;
    int iTN, iTNPert = 0;
    int iaIndexes[4];
    char *cPaTypes[4];
    int iImproper = 0, iParmOffset = 0, iTorsionOffset = 0;
    int iCount = 0, i, iIndex;

#define                MAX_N                9999

    bFailedGeneratingParameters = FALSE;

    VP0("Building %s torsion parameters.\n",
         (bProper ? "proper" : "improper"));
    if (!bProper) {
        iParmOffset = iParmSetProperCount(uUnit->psParameters);
        iTorsionOffset = iVarArrayElementCount(uUnit->vaTorsions);
    }
    /*
     *  NOTE: In order for the 1-4 interactions to be calculated
     *  properly, add constraint bonds and angles AFTER
     *  the proper torsions are added.  This allows the 1-4
     *  interaction checker to use the bond lists and angle
     *  lists to check connectivity
     */
    if (bProper) {
        int iMax;

        /*
         *  set up 1-4 checking stuff by initializing index with bond
         *  and angle atom pairs
         *
         *  create index, non-duplicate keys, key length == 2 ints
         */
        create_index(&scr14_index, IX_DUPKEYREC, 2 * sizeof(int));

        /*
         *  set up convenience pointers for stuffing key
         *      with 2 ints
         */
        Pint1 = (int *) &e14.key;
        Pint2 = Pint1 + 1;

        /*
         *  Using the index pointer as a flag; 1==bond,angle, 0==torsion.
         *      The bonds should all be unique pairs. An angle could
         *      duplicate a bond if there is a 'triangle', and could
         *      duplicate an angle if there is a 'square'.
         */
        e14.recptr = (IX_RECPOS) 1;

        /* plop in bonded pairs if any */
        if ((iMax = iVarArrayElementCount(uUnit->vaBonds))) {
            SAVEBONDt *sbPBondT = PVAI(uUnit->vaBonds, SAVEBONDt, 0);

            for (i = 0; i < iMax; i++, sbPBondT++) {
                if (sbPBondT->iAtom1 < sbPBondT->iAtom2) {
                    *Pint1 = sbPBondT->iAtom1;
                    *Pint2 = sbPBondT->iAtom2;
                } else {
                    *Pint1 = sbPBondT->iAtom2;
                    *Pint2 = sbPBondT->iAtom1;
                }
                if (add_key(&e14, &scr14_index) != IX_OK)
                    DFATAL("1-4: cannot add bond %d %d\n%s", *Pint1, *Pint2,
                            "This may be caused by duplicate bond "
                            "specifications;\n"
                            "for example, explicit bond commands in addition "
                            "to PDB conect records.\n");
            }
        }
        /* plop in angle pairs */
        if ((iMax = iVarArrayElementCount(uUnit->vaAngles))) {
            SAVEANGLEt *saPAngleT = PVAI(uUnit->vaAngles, SAVEANGLEt, 0);
            for (i = 0; i < iMax; i++, saPAngleT++) {
                if (saPAngleT->iAtom1 < saPAngleT->iAtom3) {
                    *Pint1 = saPAngleT->iAtom1;
                    *Pint2 = saPAngleT->iAtom3;
                } else {
                    *Pint1 = saPAngleT->iAtom3;
                    *Pint2 = saPAngleT->iAtom1;
                }
                if (add_key(&e14, &scr14_index) != IX_OK){
                    VPNOTE("1-4: angle %d %d %s %s\n",
                         *Pint1, *Pint2,
                         "duplicates bond ('triangular' bond)",
                         "or angle ('square' bond)\n");
                }
            }
        }
    } else {
        /*
         *  Impropers - use an index to track instantiations of
         *      the pure types templated in the ff to actual
         *      residue & atom names. For this, we use string-based
         *      indexing to count duplicate instantiations, as opposed
         *      to the fixed-length (integer) records used for
         *      1-4 tracking in the proper torsions.
         *
         *  create index, 'count'-style duplicate keys, key length == string
         *  create record w/ key portion longer than the default
         */
        create_index(&improper_index, IX_DUPKEYREC, IX_LEN_CSTRING);
        MALLOC(ePImp, IX_REC *, sizeof(IX_REC) + 80);
        /*
         *  set recptr to avoid mem check complaint
         *      (since this index doesn't use its
         *      recptr)
         */
        ePImp->recptr = NULL;
    }

    lTemp = lLoop((OBJEKT) uUnit, (bProper ? PROPERS : IMPROPERS));

    while (oNext(&lTemp) != NULL) {
        if (bProper) {
            LoopGetTorsion(&lTemp, &aAtom1, &aAtom2, &aAtom3, &aAtom4);
        } else {
            LoopGetImproper(&lTemp, &aAtom1, &aAtom2, &aAtom3, &aAtom4);
        }

/*
        fprintf( stderr, "%s torsion %s %s %s %s\n",
                        (bProper ? "Proper" : "Improper"),
                        sAtomName(aAtom1), sAtomName(aAtom2),
                        sAtomName(aAtom3), sAtomName(aAtom4) );
*/

        MESSAGE("%s torsion %s %s %s %s\n",
                 (bProper ? "Proper" : "Improper"),
                 sAtomName(aAtom1), sAtomName(aAtom2),
                 sAtomName(aAtom3), sAtomName(aAtom4));

        stTorsion.bProper = bProper;
        stTorsion.bCalc14 = FALSE;
        stTorsion.bPertCalc14 = FALSE;


        /* iContainerTempInt contains atom indices into SAVEATOMt arrays */

        stTorsion.iAtom1 = iContainerTempInt(aAtom1);
        stTorsion.iAtom2 = iContainerTempInt(aAtom2);
        stTorsion.iAtom3 = iContainerTempInt(aAtom3);
        stTorsion.iAtom4 = iContainerTempInt(aAtom4);
        stTorsion.fFlags = 0;

        /* Define the names of the unperturbed ATOMs */

        strcpy(sAtom1, sAtomType(aAtom1));
        strcpy(sAtom2, sAtomType(aAtom2));
        strcpy(sAtom3, sAtomType(aAtom3));
        strcpy(sAtom4, sAtomType(aAtom4));
        strcpy(sOrigAtom1, sAtom1);
        strcpy(sOrigAtom2, sAtom2);
        strcpy(sOrigAtom3, sAtom3);
        strcpy(sOrigAtom4, sAtom4);

        /* skip if this torsion refers to an extra point:  */
        if( GDefaults.iDeleteExtraPointAngles ){
            if( strcmp( sAtom1, "EP" ) == 0 || strcmp( sAtom4, "EP" ) == 0 )
            continue;
        }

        /* Check if the torsion is to be perturbed, if it */
        /* is then set flags saying so, and create a TORSION for */
        /* the perturbation */
        bPerturbTorsion = FALSE;
        if (bPert &&
            (bAtomFlagsSet(aAtom1, ATOMPERTURB) ||
             bAtomFlagsSet(aAtom2, ATOMPERTURB) ||
             bAtomFlagsSet(aAtom3, ATOMPERTURB) ||
             bAtomFlagsSet(aAtom4, ATOMPERTURB))) {
            bPerturbTorsion = TRUE;


            /* Note that the torsion is perturbed and whether or */
            /* not it is on the boundary between perturbed and */
            /* non-perturbed */

            stTorsion.fFlags |= PERTURBED;
            if (!(bAtomFlagsSet(aAtom1, ATOMPERTURB) &&
                  bAtomFlagsSet(aAtom2, ATOMPERTURB) &&
                  bAtomFlagsSet(aAtom3, ATOMPERTURB) &&
                  bAtomFlagsSet(aAtom4, ATOMPERTURB))) {
                stTorsion.fFlags |= BOUNDARY;
                MESSAGE("Boundary torsion: %s-%s-%s-%s\n",
                         sAtom1, sAtom2, sAtom3, sAtom4);
            }

            /* Define the names of the perturbed atoms */
            if (bAtomFlagsSet(aAtom1, ATOMPERTURB) &&
                strlen(sAtomPertType(aAtom1)) != 0)
                strcpy(sPert1, sAtomPertType(aAtom1));
            else
                strcpy(sPert1, sAtom1);
            if (bAtomFlagsSet(aAtom2, ATOMPERTURB) &&
                strlen(sAtomPertType(aAtom2)) != 0)
                strcpy(sPert2, sAtomPertType(aAtom2));
            else
                strcpy(sPert2, sAtom2);
            if (bAtomFlagsSet(aAtom3, ATOMPERTURB) &&
                strlen(sAtomPertType(aAtom3)) != 0)
                strcpy(sPert3, sAtomPertType(aAtom3));
            else
                strcpy(sPert3, sAtom3);
            if (bAtomFlagsSet(aAtom4, ATOMPERTURB) &&
                strlen(sAtomPertType(aAtom4)) != 0)
                strcpy(sPert4, sAtomPertType(aAtom4));
            else
                strcpy(sPert4, sAtom4);
            strcpy(sOrigPert1, sPert1);
            strcpy(sOrigPert2, sPert2);
            strcpy(sOrigPert3, sPert3);
            strcpy(sOrigPert4, sPert4);
        }

        /* First search through the UNITs PARMSET for the */
        /* torsion parameters that we need. */

        tTorsion = tParmSetTORSIONCreate();

        if (bProper) {
            iTN = iParmSetFindProperTerms(uUnit->psParameters,
                                          tTorsion, TRUE,
                                          sAtom1, sAtom2, sAtom3, sAtom4);
        } else {
            iTN = iParmSetFindImproperTerms(uUnit->psParameters,
                                            tTorsion, TRUE,
                                            sAtom1, sAtom2,
                                            sAtom3, sAtom4);
        }

        if (iTN != PARM_NOT_FOUND) {
            MESSAGE("Found existing %s terms in UNIT PARMSET.\n",
                     (bProper ? "PROPER" : "IMPROPER"));
        }
        if (bPerturbTorsion) {
            tPertTorsion = tParmSetTORSIONCreate();
            if (bProper) {
                iTNPert = iParmSetFindProperTerms(uUnit->psParameters,
                                                  tPertTorsion, TRUE,
                                                  sPert1, sPert2,
                                                  sPert3, sPert4);
            } else {
                iTNPert = iParmSetFindImproperTerms(uUnit->psParameters,
                                                    tPertTorsion, TRUE,
                                                    sPert1, sPert2,
                                                    sPert3, sPert4);
            }
            if (iTNPert != PARM_NOT_FOUND) {
                MESSAGE("Found existing %s terms in UNIT PARMSET.\n",
                         (bProper ? "PROPER" : "IMPROPER"));
            }
        }

        if (iTN != PARM_FOUND_EXACT ||
            (bPerturbTorsion && iTNPert != PARM_FOUND_EXACT)) {

            /* search through the PARMLIBs */
            /* put 1st torsion into tTorsion since
               hopefully we're searching back from most
               recent parm loaded */

            PARMLIB_LOOP_ALL(plLib, psTemp) {
                if (iTN != PARM_FOUND_EXACT) {
                    if (bProper) {
                        iTN = iParmSetFindProperTerms(psTemp,
                                                      tTorsion, FALSE,
                                                      sAtom1, sAtom2,
                                                      sAtom3, sAtom4);
                    } else {
                        iTN = iParmSetFindImproperTerms(psTemp,
                                                        tTorsion, FALSE,
                                                        sAtom1, sAtom2,
                                                        sAtom3, sAtom4);
                    }
                }
                if (bPerturbTorsion && iTNPert != PARM_FOUND_EXACT) {
                    if (bProper) {
                        iTNPert = iParmSetFindProperTerms(psTemp,
                                                          tPertTorsion,
                                                          FALSE, sPert1,
                                                          sPert2, sPert3,
                                                          sPert4);
                    } else {
                        iTNPert = iParmSetFindImproperTerms(psTemp,
                                                            tPertTorsion,
                                                            FALSE, sPert1,
                                                            sPert2, sPert3,
                                                            sPert4);
                    }
                }
                if (iTN == PARM_FOUND_EXACT) {
                    if (bPerturbTorsion) {
                        if (iTNPert == PARM_FOUND_EXACT)
                            break;
                    } else
                        break;
                }
            }
        }
#ifdef        DEBUG2
        MESSAGE("%s %s-%s-%s-%s found %d terms\n",
                 (bProper ? "PROPER" : "IMPROPER"),
                 sAtom1, sAtom2, sAtom3, sAtom4,
                 iParmSetTORSIONTermCount(tTorsion));
        for (i = 0; i < iParmSetTORSIONTermCount(tTorsion); i++) {
            ParmSetTORSIONTerm(tTorsion, i,
                               &iTParm,
                               sT1, sT2, sT3, sT4,
                               &iTmp, &dTK, &dTP, sTemp);
            MESSAGE("Term %3d  %d %s-%s-%s-%s  %d  %lf  %lf\n",
                     i, iTParm, sT1, sT2, sT3, sT4, iTmp, dTK, dTP);
        }
        if (bPerturbTorsion) {
            MESSAGE("Pert%s %s-%s-%s-%s found %d terms\n",
                     (bProper ? "PROPER" : "IMPROPER"),
                     sPert1, sPert2, sPert3, sPert4,
                     iParmSetTORSIONTermCount(tPertTorsion));
            for (i = 0; i < iParmSetTORSIONTermCount(tPertTorsion); i++) {
                ParmSetTORSIONTerm(tPertTorsion, i,
                                   &iTParm,
                                   sT1, sT2, sT3, sT4,
                                   &iTmp, &dTK, &dTP, sTemp);
                MESSAGE("Term %3d  %d %s-%s-%s-%s  %d  %lf  %lf\n",
                         i, iTParm, sT1, sT2, sT3, sT4, iTmp, dTK, dTP);
            }
        }
#endif


        /* Now we have a complete TORSION in tTorsion, and */
        /* if the torsion is being perturbed we have another */
        /* complete TORSION in tPertTorsion */
        /* The (stTorsion) terms must now be created */

        iTerm = 0;
        iPertTerm = 0;
        bDone = FALSE;
        iIndex = PARM_NOT_FOUND;
        iPertIndex = PARM_NOT_FOUND;
        bEnd = FALSE;
        bPertEnd = FALSE;
        iN = MAX_N;
        iPertN = MAX_N;
        bCalc14 = TRUE;
        bCalcPert14 = TRUE;

        /*
         *  get 1st term
         */
        if (iParmSetTORSIONTermCount(tTorsion) != 0) {
            ParmSetTORSIONTerm(tTorsion, iTerm,
                               &iIndex,
                               sAtom1, sAtom2, sAtom3, sAtom4,
                               &iN, &dKp, &dP0, &dScEE, &dScNB,
			       sDesc);
            MESSAGE("First non-perturbed multiplicity: %d\n", iN);
        } else {
            if (bProper) {
                VPTRACE("sOrigAtom vs sAtom: %s, %s, %s, %s;\n"
                   "                           %s, %s, %s, %s.\n",
                        sOrigAtom1, sOrigAtom2, sOrigAtom3, sOrigAtom4,
                            sAtom1,     sAtom2,     sAtom3,     sAtom4 );
                /* My reading is that the parts of aAtom1, etc. used below */
                /* cannot change, so i use sAtom1 instead of sOrigAtom1, etc. */
                VECTOR vPos1 = vAtomPosition(aAtom1);
                VECTOR vPos2 = vAtomPosition(aAtom2);
                VECTOR vPos3 = vAtomPosition(aAtom3);
                VECTOR vPos4 = vAtomPosition(aAtom4);
                RESIDUE rRes1 = (RESIDUE)cContainerWithin(aAtom1);
                RESIDUE rRes2 = (RESIDUE)cContainerWithin(aAtom2);
                RESIDUE rRes3 = (RESIDUE)cContainerWithin(aAtom3);
                RESIDUE rRes4 = (RESIDUE)cContainerWithin(aAtom4);
                VPERROR(" ** No torsion terms for atom types: %-s-%-s-%-s-%-s\n"
                        "        for atom %s.%s:%d.%s at position %8.3f,%8.3f,%8.3f,\n"
                        "            atom %s.%s:%d.%s at position %8.3f,%8.3f,%8.3f,\n"
                        "            atom %s.%s:%d.%s at position %8.3f,%8.3f,%8.3f,\n"
                        "        and atom %s.%s:%d.%s at position %8.3f,%8.3f,%8.3f.\n",
                        sAtom1, sAtom2, sAtom3, sAtom4,
                        sContainerName(rRes1),sResidueChainId(rRes1),iResiduePdbSequence(rRes1),sContainerName(aAtom1),
                        vPos1.dX, vPos1.dY, vPos1.dZ,
                        sContainerName(rRes2),sResidueChainId(rRes2),iResiduePdbSequence(rRes2),sContainerName(aAtom2),
                        vPos2.dX, vPos2.dY, vPos2.dZ,
                        sContainerName(rRes3),sResidueChainId(rRes3),iResiduePdbSequence(rRes3),sContainerName(aAtom3),
                        vPos3.dX, vPos3.dY, vPos3.dZ,
                        sContainerName(rRes4),sResidueChainId(rRes4),iResiduePdbSequence(rRes4),sContainerName(aAtom4),
                        vPos4.dX, vPos4.dY, vPos4.dZ);
                bFailedGeneratingParameters = TRUE;
            } else if ( iAtomHybridization(aAtom3) == 2 ){
                RESIDUE rRes1 = (RESIDUE)cContainerWithin(aAtom1);
                RESIDUE rRes2 = (RESIDUE)cContainerWithin(aAtom2);
                RESIDUE rRes3 = (RESIDUE)cContainerWithin(aAtom3);
                RESIDUE rRes4 = (RESIDUE)cContainerWithin(aAtom4);
                VP1(" ** Warning: No sp2 improper torsion term for  %-s-%-s-%-s-%-s\n",
                     sOrigAtom1, sOrigAtom2, sOrigAtom3, sOrigAtom4);
                        VP1("        atoms are: %s.%s:%d.%s %s.%s:%d.%s %s.%s:%d.%s %s.%s:%d.%s\n",
                        sContainerName(rRes1),sResidueChainId(rRes1),iResiduePdbSequence(rRes1),sContainerName(aAtom1),
                        sContainerName(rRes2),sResidueChainId(rRes2),iResiduePdbSequence(rRes2),sContainerName(aAtom2),
                        sContainerName(rRes3),sResidueChainId(rRes3),iResiduePdbSequence(rRes3),sContainerName(aAtom3),
                        sContainerName(rRes4),sResidueChainId(rRes4),iResiduePdbSequence(rRes4),sContainerName(aAtom4));
            }
            bEnd = TRUE;
        }
        if (bPerturbTorsion) {
            if (iParmSetTORSIONTermCount(tPertTorsion) != 0) {
                ParmSetTORSIONTerm(tPertTorsion, iPertTerm,
                                   &iPertIndex,
                                   sPert1, sPert2, sPert3, sPert4,
                                   &iPertN, &dPertKp, &dPertP0, &dPScEE,
				   &dPScNB, sDesc);
                MESSAGE("First perturbed multiplicity: %d\n", iPertN);
            } else
                bPertEnd = TRUE;
        }
        if (!bPerturbTorsion)
            bDone = bEnd;
        else
            bDone = bEnd && bPertEnd;

        /* If the interaction is an improper then reorder */
        /* the atoms to get them in the same order as */
        /* was in the original parameter set */

        if (!bProper && iParmSetTORSIONTermCount(tTorsion) > 0) {
            cPaTypes[0] = sAtomType(aAtom1);
            cPaTypes[1] = sAtomType(aAtom2);
            cPaTypes[2] = sAtomType(aAtom3);
            cPaTypes[3] = sAtomType(aAtom4);
            iaIndexes[0] = stTorsion.iAtom1;
            iaIndexes[1] = stTorsion.iAtom2;
            iaIndexes[2] = stTorsion.iAtom3;
            iaIndexes[3] = stTorsion.iAtom4;
            MESSAGE("Old order: %d %d %d %d\n",
                     stTorsion.iAtom1,
                     stTorsion.iAtom2,
                     stTorsion.iAtom3, stTorsion.iAtom4);

            ParmSetImproperOrderAtoms(tTorsion, 0, cPaTypes, iaIndexes);

            stTorsion.iAtom1 = iaIndexes[0];
            stTorsion.iAtom2 = iaIndexes[1];
            stTorsion.iAtom3 = iaIndexes[2];
            stTorsion.iAtom4 = iaIndexes[3];
            MESSAGE("New order: %d %d %d %d\n",
                     stTorsion.iAtom1,
                     stTorsion.iAtom2,
                     stTorsion.iAtom3, stTorsion.iAtom4);
        }

        /* Loop over all of the terms */

        while (!bDone) {
            bUse = FALSE;
            bUsePert = FALSE;
            bCopy = FALSE;
            bCopyPert = FALSE;
            stTorsion.iParmIndex = 0;
            stTorsion.iPertParmIndex = 0;

            /* Get the next term, the one with the */
            /* lowest multiplicity, and advance to the */
            /* next multiplicity */

            if (!bPerturbTorsion) {
                /* Advance to the next term within the TORSION */
                bUse = TRUE;
                if (iIndex == PARM_NOT_FOUND)
                    bCopy = TRUE;
            } else {

                /* If the multiplicity is the same for */
                /* both the nonperturbed and perturbed */
                /* term then advance them both */
                if (iPertN == iN) {
                    bUse = TRUE;
                    bUsePert = TRUE;
                    if (iIndex == PARM_NOT_FOUND)
                        bCopy = TRUE;
                    if (iPertIndex == PARM_NOT_FOUND)
                        bCopyPert = TRUE;
                } else if (iN < iPertN) {
                    bUse = TRUE;
                    if (iIndex == PARM_NOT_FOUND)
                        bCopy = TRUE;
                } else {
                    bUsePert = TRUE;
                    if (iPertIndex == PARM_NOT_FOUND)
                        bCopyPert = TRUE;
                }
            }

            MESSAGE("Flags:  bUse:%d bCopy:%d  bUsePert:%d bCopyPert:%d\n",
                     bUse, bCopy, bUsePert, bCopyPert);

            /* Now save the terms into the UNITs PARMSET if */
            /* they are not already there */

            if (bCopy) {
                if (bProper)
                    iIndex = iParmSetAddProperTerm(uUnit->psParameters,
                                                   sAtom1, sAtom2, sAtom3,
                                                   sAtom4, iN, dKp, dP0,
						   dScEE, dScNB,
                                                   sDesc);
/*                else if ( !GDefaults.iCharmm )    ???---should I do this????     */
                else
                    iIndex = iParmSetAddImproperTerm(uUnit->psParameters,
                                                     sAtom1, sAtom2,
                                                     sAtom3, sAtom4, iN,
                                                     dKp, dP0, dScEE, dScNB,sDesc);
            }
            if (bCopyPert) {
                if (bProper) {
                    iPertIndex = iParmSetAddProperTerm(uUnit->psParameters,
                                                       sPert1, sPert2,
                                                       sPert3, sPert4,
                                                       iPertN, dPertKp,
                                                       dPertP0, dScEE,
						       dScNB, sDesc);
                } else {
                    iPertIndex =
                        iParmSetAddImproperTerm(uUnit->psParameters,
                                                sPert1, sPert2, sPert3,
                                                sPert4, iPertN, dPertKp,
                                                dPertP0,dScEE,  dScNB,  sDesc);
                }
                MESSAGE("iPertIndex = %d\n", iPertIndex);
            }

            /* Save the multiplicity */

            iLastN = iN;
            iLastPertN = iPertN;
            if (bUse) {
                stTorsion.iParmIndex = iParmOffset + iIndex + 1;
                iTerm++;
                if (iTerm >= iParmSetTORSIONTermCount(tTorsion)) {
                    bEnd = TRUE;
                    iN = MAX_N;
                } else {
                    ParmSetTORSIONTerm(tTorsion, iTerm,
                                       &iIndex,
                                       sAtom1, sAtom2, sAtom3, sAtom4,
                                       &iN, &dKp, &dP0, &dScEE, &dScNB,
				       sDesc);
                }
                MESSAGE("Advancing non-perturbed multiplicity to %d\n",
                         iN);
            }
            if (bUsePert) {
                stTorsion.iPertParmIndex = iParmOffset + iPertIndex + 1;
                iPertTerm++;
                if (iPertTerm >= iParmSetTORSIONTermCount(tPertTorsion)) {
                    bPertEnd = TRUE;
                    iPertN = MAX_N;
                } else {
                    ParmSetTORSIONTerm(tPertTorsion, iPertTerm,
                                       &iPertIndex,
                                       sPert1, sPert2, sPert3, sPert4,
                                       &iPertN, &dPertKp, &dPertP0, &dScEE,
				       &dScNB, sDesc);
                }
                MESSAGE("Advancing perturbed multiplicity to %d\n",
                         iPertN);
            }
            if (bProper)
                zUnitIOSetCalc14Flags(&stTorsion, &bCalc14, &bCalcPert14);
            /*
             *  If we are writing a perturbation topology file,
             *  for every unperturbed term for a given torsion,
             *  there must be a perturbed term (and vice-versa).
             *  This is because multiple torsional potentials
             *  may apply to a single torsion, and each is perturbed
             *  individually in gibbs.
             */

            /* At this point iLastN, iLastPertN should be the */
            /* multiplicity of the torsion term */

            if (bPerturbTorsion) {
                if (stTorsion.iParmIndex == 0 ||
                    stTorsion.iPertParmIndex == 0) {
                    bFailedGeneratingParameters = TRUE;
                    VPWARN("*** %s torsion parameters missing ***\n",
                         (bProper ? "Proper" : "Improper"));
                    VP0(" atom names: %-s-%-s-%-s-%-s\n",
                         sAtomName(aAtom1), sAtomName(aAtom2),
                         sAtomName(aAtom3), sAtomName(aAtom4));

                    VP0(" atom types: %-s-%-s-%-s-%-s  =pert=>  %-s-%-s-%-s-%-s\n",
                         sOrigAtom1, sOrigAtom2, sOrigAtom3, sOrigAtom4,
                         sOrigPert1, sOrigPert2, sOrigPert3, sOrigPert4);
                    /* note: missing multiplicity is for the _other_ state */
                    VP0("Please add a dummy parameter of multiplicity %d\n",
                         (stTorsion.iParmIndex !=
                          0 ? iLastN : iLastPertN));
                    VP0("for the %spert types to your parameter set.\n",
                         (stTorsion.iParmIndex == 0 ? "non-" : ""));
                    VP0(" - e.g. %-s-%-s-%-s-%-s  %s    0.0     0.       %d.\n",
                         (stTorsion.iParmIndex == 0 ? sOrigAtom1 : sOrigPert1),
                         (stTorsion.iParmIndex == 0 ? sOrigAtom2 : sOrigPert2),
                         (stTorsion.iParmIndex == 0 ? sOrigAtom3 : sOrigPert3),
                         (stTorsion.iParmIndex == 0 ? sOrigAtom4 : sOrigPert4),
                         (bProper ? "1" : ""),
                         (stTorsion.iParmIndex != 0 ? iLastN : iLastPertN));

                    VP0("%s %s\n%s %s\n", "(This is because multiple",
                         "torsional potentials may apply to a",
                         "single torsion, and each is perturbed",
                         "individually in gibbs.)");
                }
            }

            /* Add the term to the table */
            VarArrayAdd(uUnit->vaTorsions, (GENP) & stTorsion);
            iCount++;
            if (!bProper) {
                STRING sDesc1, sDesc2, sDesc3, sDesc4;

                AtomDescStr(aAtom1, FALSE, sDesc1);
                AtomDescStr(aAtom2, FALSE, sDesc2);
                AtomDescStr(aAtom3, FALSE, sDesc3);
                AtomDescStr(aAtom4, FALSE, sDesc4);
                sprintf(ePImp->key, "%s - %s - %s - %s",
                        sDesc1, sDesc2, sDesc3, sDesc4);
                if (add_key(ePImp, &improper_index) != IX_OK)
                    DFATAL("add_key() impropers\n");

                iImproper++;
            }
#ifdef        DEBUG2
            MESSAGE("Adding %s : %s - %s - %s - %s\n", bProper ? "PROPER" : "IMPROPER"),
                     sContainerFullDescriptor((CONTAINER) aAtom1, s1),
                     sContainerFullDescriptor((CONTAINER) aAtom2, s2),
                     sContainerFullDescriptor((CONTAINER) aAtom3, s3),
                     sContainerFullDescriptor((CONTAINER) aAtom4, s4);
            MESSAGE("Perturbed: %s\n", sBOOL(bPerturbTorsion));
#endif
            if (!bPerturbTorsion) {
                bDone = bEnd;
            } else {
                bDone = bEnd && bPertEnd;
            }
        }

        /* Release the memory used by the TORSIONs */

        ParmSetTORSIONDestroy(&tTorsion);

        if (bPerturbTorsion) {
            ParmSetTORSIONDestroy(&tPertTorsion);
        }
    }

    if (bProper) {

        /*
         *  throw away 1-4 scratch index
         */
        destroy_index(&scr14_index);

        /*
         *  'boil down' numerically redundant params
         */
        BoilTorsions(&uUnit->psParameters->vaTorsions, 0,        /* parm table */
                     uUnit->vaTorsions, 0);        /* atom lists */

    } else {
        RESIDUE rRes;
        IMPROPERt *IP;
        int iImproper2 = 0, iPrep = 0;

        BoilTorsions(&uUnit->psParameters->vaImpropers,        /* parm table */
                     iParmOffset, uUnit->vaTorsions,        /* atom lists */
                     iTorsionOffset);

        lTemp = lLoop((OBJEKT) uUnit, RESIDUES);
        while ((rRes = (RESIDUE) oNext(&lTemp))) {
            if (rRes->vaImpropers) {
                if (iPrep++ == 0)
                    VP0("old PREP-specified impropers:\n");
            }
            if ((iCount = iVarArrayElementCount(rRes->vaImpropers))) {
                sContainerDescriptor((CONTAINER) rRes, sDesc);
                IP = PVAI(rRes->vaImpropers, IMPROPERt, 0);
                for (i = 0; i < iCount; i++, IP++) {
                    iImproper2++;
                    VP0(" %s:  %.4s %.4s %.4s %.4s\n", sDesc + 1,
                         IP->sName1, IP->sName2, IP->sName3, IP->sName4);
                }
            }
        }
        if (iImproper && GiVerbosityLevel > 0) {

            VP0("--Impropers:\n");
            first_key(&improper_index);
            while (next_key(ePImp, &improper_index) == IX_OK)
                VP1("  %d\t%s\n", ePImp->count, ePImp->key);
        }
        VP0(" total %d improper torsion%s applied\n",
             iImproper, (iImproper != 1 ? "s" : ""));
        if (iPrep)
            VP0(" %d improper torsions in old prep form\n", iImproper2);
        destroy_index(&improper_index);
        FREE(ePImp);
    }


    VPTRACEEXIT("zbUnitIOIndexTorsionParameters" );
    return (bFailedGeneratingParameters);
}




/*
 *        zUnitDoResidue
 *
 *        Build a RESIDUE entry into the uUnit->vaResidues table.
 */
static void
zUnitDoResidue(UNIT uUnit, RESIDUE rRes, int *iPPos)
{
    SAVERESIDUEt *srPResidue;

    ContainerSetTempInt(rRes, (*iPPos) + 1);

    srPResidue = PVAI(uUnit->vaResidues, SAVERESIDUEt, *iPPos);
    srPResidue->rResidue = rRes;
    REF(rRes);
    strcpy(srPResidue->sName, sContainerName(rRes));
    srPResidue->sResidueType[0] = cResidueType(rRes);
    srPResidue->sResidueType[1] = '\0';
    srPResidue->iSequenceNumber = iContainerSequence(rRes);
    srPResidue->iNextChildSequence = iContainerNextChildsSequence(rRes);
    srPResidue->iPdbResSeq = iResiduePdbSequence(rRes);
    memcpy(srPResidue->sChainId,rRes->sChainId,sizeof(srPResidue->sChainId));
    if (rRes->cICode == ' ') srPResidue->sICode[0]=0;
    else { srPResidue->sICode[0]=rRes->cICode;srPResidue->sICode[1]=0; }
    (*iPPos)++;

}


/*
 *        zUnitDoAtoms
 *
 *        Loop over all of the ATOMs within the RESIDUE and add them
 *        to the UNITs tables.
 *
 *        If the RESIDUE is of type RESTYPESOLVENT then the first
 *        ATOM will be the imaging atom, this is to maintain compatibility
 *        with AMBER parm files.
 */
void
zUnitDoAtoms(UNIT uUnit, PARMLIB plParameters, RESIDUE rRes, int *iPPos,
             BOOL * bPFailed, BOOL bPert)
{
    SAVERESIDUEt *srPResidue;
    LOOP lTemp;
    ATOM aAtom, aIgnore;
    BOOL bFailed;

    srPResidue = PVAI(uUnit->vaResidues, SAVERESIDUEt,
                      iContainerTempInt(rRes) - 1);
    srPResidue->iAtomStartIndex = 0;

    /* Check if it is a solvent RESIDUE */
    /* If it is, put the imaging ATOM into the table first */

    aIgnore = NULL;
    if (cResidueType(rRes) == RESTYPESOLVENT) {
        if (aResidueImagingAtom(rRes) != NULL) {
            aIgnore = aResidueImagingAtom(rRes);
            zUnitIOTableAddAtom(uUnit, aIgnore, *iPPos,
                                plParameters, &bFailed, bPert);
            (*iPPos)++;
            *bPFailed |= bFailed;
            if (srPResidue->iAtomStartIndex == 0) {
                srPResidue->iAtomStartIndex = iContainerTempInt(aIgnore);
            }
        }
    }

    /* Now put the rest of the ATOMs */

    lTemp = lLoop((OBJEKT) rRes, DIRECTCONTENTSBYSEQNUM);
    while ((aAtom = (ATOM) oNext(&lTemp)) != NULL) {

        /* Ignore the imaging ATOM if there is one */

        if (aAtom == aIgnore)
            continue;

        zUnitIOTableAddAtom(uUnit, aAtom, *iPPos, plParameters, &bFailed,
                            bPert);
        (*iPPos)++;
        *bPFailed |= bFailed;
        if (srPResidue->iAtomStartIndex == 0) {
            srPResidue->iAtomStartIndex = iContainerTempInt(aAtom);
        }
    }
}

/*------------ Build SAVERESIDUEt array ------------*/

/* Put RESIDUEs in the following order: */
/* non-solvent residues */
/* solvent residues not in solvent cap */
/* solvent residues in solvent cap */
/* THIS IS ONLY DONE BECAUSE AMBER PARM FILES */
/* REQUIRE IT!!!!!!  PROGRAMS SHOULD ASSUME */
/* THAT RESIDUES HAVE ARBITRARY ORDERING WITHIN */
/* OBJECT FILE FORMAT FILES AND SORT */
/* THEM THEMSELVES TO PREVENT INCOMPATIBILITIES */
/* WITH FUTURE RELEASES!!!!!!!!!!!!!!!!!!!!!!!! */
// Also builds residue table, and marks RESIDUE TempInt with derived sequence order.

int
zUnitIOAmberOrderResidues( UNIT uUnit )
{
    LOOP    lResidues, lSpanning;
    RESIDUE rRes;
    int     i=0, iResidueCount=0, iMolecule=0;

    /* Clear the ATOMTOUCHED flag on all the ATOMs for molecuel labelling */
    ContainerResetAllAtomsFlags((CONTAINER) uUnit, ATOMTOUCHED);

    /* Loop through solvent RESIDUEs and:
     * 1) Count them
     *      Normally iResidueCount == iCollectionSize(cContainerContents(uUnit)))
     *      but leap allows container abuse, and it can contain anything.
     * 2) set the temp flag saying if they are in the cap or not
     * 3) label all molecule groups.
     *     ==> This is repeated in zUnitIOFindAndCountMolecules()
     *     but must be repeated after residue ordering
     *     TODO: zUnitIOFindAndCountMolecules can use labels instead
     *     of repeating the full spanning tree search
    */
    lResidues = lLoop((OBJEKT) uUnit, DIRECTCONTENTS);
    while ((rRes = (RESIDUE) oNext(&lResidues))) {
        iResidueCount++;
        if (cResidueType(rRes) == RESTYPESOLVENT) {
            if (bUnitCapContainsContainer(uUnit, (CONTAINER) rRes))
                ResidueSetFlags(rRes, RESIDUEINCAP);
            else
                ResidueResetFlags(rRes, RESIDUEINCAP);
        }
        /* Search for the next RESIDUE whose first ATOM has not */
        /* been touched */
        ATOM aAtom = (ATOM)oContainerFirstObject(rRes);
        if (!bAtomFlagsSet(aAtom, ATOMTOUCHED)) {
            /* Touch all of the ATOMs within the molecule that */
            /* contains the current RESIDUE */

            lSpanning = lLoop((OBJEKT) aAtom, SPANNINGTREE);
            FOREACH(aAtom, ATOM, lSpanning) {
                AtomSetFlags(aAtom, ATOMTOUCHED);
                CONTAINER cParent = cContainerWithin(aAtom);
                if (iObjectType(cParent) == RESIDUEid)
                    ((RESIDUE)cParent)->iTemp = iMolecule;
            }
            iMolecule++;
        }
    }

    if (iResidueCount == 0)
            return(0);

    /*
    **  allocate array for residues
    */
    uUnit->vaResidues = vaVarArrayCreate(sizeof(SAVERESIDUEt));
    VarArraySetSize((uUnit->vaResidues), iResidueCount);

    if ( !GDefaults.reorder_residues )
        printf("\"order_residues\" off: keep input residue order.\n");
    lResidues = lLoop((OBJEKT) uUnit, GDefaults.reorder_residues ?
            DIRECTCONTENTSPARMORDER : DIRECTCONTENTSBYSEQNUM);
    while ((rRes = (RESIDUE) oNext(&lResidues)) != NULL) {
       if (iObjectType(rRes) == RESIDUEid)
           zUnitDoResidue(uUnit, rRes, &i);
    }
    return(iResidueCount);
}

/*
 *      zUnitIOBuildTables
 *
 *        Author:        Christian Schafmeister (1991)
 *
 *      Build a table representation of the UNIT.
 *        plParameters specifices parameters (success returned in bPGeneratedParameters)
 *        If (bPert) then build tables for perturbation run.
 *        If (bCheck) then run UnitCheck()
 *
 *        Tables are built for both Object file (OFF) and for AmberParm output.
 *
 *        NOTE: Programmers should assume that the order of entries
 *        NOTE: within tables and OFF files is COMPLETELY ARBITRARY.
 *
 */
void
zUnitIOBuildTables(UNIT uUnit, PARMLIB plParameters,
            BOOL * bPGeneratedParameters, // set to TRUE if parameters generated
            BOOL bPert,
            BOOL bCheck)
{
    SAVECONNECTIVITYt *scPCon;
    SAVEMOLECULEt *smPMolecule;
    SAVERESTRAINTt *srPRestraint;
    SAVERESIDUEt *srPResidue;
    SAVEHIERARCHYt *shPHierarchy;
    SAVEGROUPSt sgAtomGroup;
    BAGLOOP blRestraint;
    RESTRAINT rRest;
    MOLECULE mMol;
    RESIDUE rRes;
    ATOM aAtom1, aAtom2, aAtom3, aAtom4;
    LOOP lMolecules;
    OBJEKT oAbove, oBelow;
    STRING sAtom1, sAtom2, sDesc;
    BOOL bGenerateParameters, bFailedGeneratingParameters;
    PARMSET psTemp;
    LOOP lTemp, lResidues;
    DICTLOOP dlGroups;
    LISTLOOP llAtoms;
    LIST lGroup;
    int i, j, iAtomCount, iMoleculeCount, iResidueCount;
    int iIndex, iCount, iGroup, iErrors = 0, iWarnings = 0;
    double dKx, dX0, dN, dA, dB, dEI, dEJ, dRI, dRJ;

    if (bCheck) {
        VP0("Checking Unit.\n");
        iFatal = 0;
        UnitCheck(uUnit, &iErrors, &iWarnings);
        if (iFatal) {
            VPFATALEXIT("Failed to generate parameters\n");
            *bPGeneratedParameters = FALSE;
            return;
        }
        if (iErrors || iWarnings) {
            /* just for fun, grammar */
            VPNOTE("Ignoring the %s%s%s%s%s from Unit Checking.\n\n",
                 (iErrors ? "error" : ""),
                 (iErrors > 1 ? "s" : ""),
                 ((iErrors && iWarnings) ? " and " : ""),
                 (iWarnings ? "warning" : ""),
                 (iWarnings > 1 ? "s" : ""));
        }
    }

    VP0("Building topology.\n");

    bFailedGeneratingParameters = FALSE;
    if (iUnitMode(uUnit) != UNITNORMAL) {
        DFATAL("The UNIT must be in NORMAL mode!");
    }

    UnitSetMode(uUnit, UNITTABLES);

    bGenerateParameters = FALSE;
    if (plParameters != NULL) {
        uUnit->psParameters = (PARMSET) oCreate(PARMSETid);
        bGenerateParameters = TRUE;
    }

    /* NOTE: molecules are not used, iMoleculeCount is always zero */
    /* Build the molecule information */

    iMoleculeCount = 0;
    lMolecules = lLoop((OBJEKT) uUnit, MOLECULES);
    while (oNext(&lMolecules) != NULL)
        iMoleculeCount++;

    if (iMoleculeCount) {
        uUnit->vaMolecules = vaVarArrayCreate(sizeof(SAVEMOLECULEt));
        VarArraySetSize((uUnit->vaMolecules), iMoleculeCount);

        lMolecules = lLoop((OBJEKT) uUnit, MOLECULES);
        smPMolecule = PVAI(uUnit->vaMolecules, SAVEMOLECULEt, 0);
        for (i = 0; (mMol = (MOLECULE) oNext(&lMolecules));
             i++, smPMolecule++) {
            ContainerSetTempInt(mMol, i + 1);
            smPMolecule->mMolecule = mMol;
            REF(mMol);
            strcpy(smPMolecule->sName, sContainerName(mMol));
            smPMolecule->iSequenceNumber = iContainerSequence(mMol);
            smPMolecule->iNextChildSequence =
                iContainerNextChildsSequence(mMol);
        }
    }


    /* Build the residue information */

    /* zUnitIOAmberOrderResidues generates the vaResidues array,
     * and puts residues in arbitrary amber order
     */
    iResidueCount = zUnitIOAmberOrderResidues( uUnit );

    if (iResidueCount) {

        /* Build the array for the atoms */
        // Now we repeat THE SAME LOGIC as in zUnitIOAmberOrderResidues()
        // or we are screwed!

        VP0("Building atom parameters.\n");

        iAtomCount = 0;
        lTemp = lLoop((OBJEKT) uUnit, ATOMS);
        while (oNext(&lTemp) != NULL)
            iAtomCount++;

        if (iAtomCount) {
            uUnit->vaAtoms = vaVarArrayCreate(sizeof(SAVEATOMt));
            VarArraySetSize((uUnit->vaAtoms), iAtomCount);
            i = 0;

            /* Put solvent ATOMs after other residues */
            /* THIS IS ONLY DONE BECAUSE AMBER PARM FILES */
            /* REQUIRE IT!!!!!!  PROGRAMS SHOULD ASSUME */
            /* THAT ATOMS HAVE ARBITRARY ORDERING AND SORT */
            /* THEM THEMSELVES TO PREVENT INCOMPATIBILITIES */
            /* WITH FUTURE RELEASES!!!!!!!!!!!!!!!!!!!!!!!! */

            lResidues = lLoop((OBJEKT) uUnit, GDefaults.reorder_residues ?
                      DIRECTCONTENTSPARMORDER : DIRECTCONTENTSBYSEQNUM);
            uUnit->iCapTempInt=0;
            while ((rRes = (RESIDUE) oNext(&lResidues)) != NULL) {
                /* Set the uUnit->iCapTempInt integer to point to */
                /* the last ATOM that is not CAP solvent */
                if ( !uUnit->iCapTempInt &&
                        cResidueType(rRes) == RESTYPESOLVENT &&
                        bResidueFlagsSet(rRes, RESIDUEINCAP))
                    uUnit->iCapTempInt = i;
                zUnitDoAtoms(uUnit, plParameters, rRes, &i,
                         &bFailedGeneratingParameters, bPert);
            }
        }
    }

/*
 *. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
 *
 *        Now the ATOMs have their indices stored in Container.iTempInt
 *
 */

    /* Generate the table storing the RESTRAINTs */

    if (iBagSize(uUnit->bRestraints)) {
        uUnit->vaRestraints = vaVarArrayCreate(sizeof(SAVERESTRAINTt));
        VarArraySetSize((uUnit->vaRestraints),
                        iBagSize(uUnit->bRestraints));
        srPRestraint = PVAI(uUnit->vaRestraints, SAVERESTRAINTt, 0);
        blRestraint = blBagLoop(uUnit->bRestraints);
        while ((rRest = (RESTRAINT) PBagNext(&blRestraint))) {
            srPRestraint->iType = iRestraintType(rRest);
            srPRestraint->fFlags = fRestraintFlags(rRest);
            switch (iRestraintType(rRest)) {
            case RESTRAINTBOND:
                RestraintBondGet(rRest, &aAtom1, &aAtom2, &dKx, &dX0);
                srPRestraint->iAtom1 = iContainerTempInt(aAtom1);
                srPRestraint->iAtom2 = iContainerTempInt(aAtom2);
                srPRestraint->iAtom3 = 0;
                srPRestraint->iAtom4 = 0;
                srPRestraint->dKx = dKx;
                srPRestraint->dX0 = dX0;
                srPRestraint->dN = 0.0;
                break;
            case RESTRAINTANGLE:
                RestraintAngleGet(rRest, &aAtom1, &aAtom2, &aAtom3,
                                  &dKx, &dX0);
                srPRestraint->iAtom1 = iContainerTempInt(aAtom1);
                srPRestraint->iAtom2 = iContainerTempInt(aAtom2);
                srPRestraint->iAtom3 = iContainerTempInt(aAtom3);
                srPRestraint->iAtom4 = 0;
                srPRestraint->dKx = dKx;
                srPRestraint->dX0 = dX0;
                srPRestraint->dN = 0.0;
                break;
            case RESTRAINTTORSION:
                RestraintTorsionGet(rRest, &aAtom1, &aAtom2,
                                    &aAtom3, &aAtom4, &dKx, &dX0, &dN);
                srPRestraint->iAtom1 = iContainerTempInt(aAtom1);
                srPRestraint->iAtom2 = iContainerTempInt(aAtom2);
                srPRestraint->iAtom3 = iContainerTempInt(aAtom3);
                srPRestraint->iAtom4 = iContainerTempInt(aAtom4);
                srPRestraint->dKx = dKx;
                srPRestraint->dX0 = dX0;
                srPRestraint->dN = dN;
                break;
            default:
                DFATAL("Invalid restraint type!");
            }
            srPRestraint++;
        }
    }


    /* Generate the indices for the RESIDUE connect atoms */
    /* and Imaging ATOMs */

    if (iVarArrayElementCount(uUnit->vaResidues)) {
        srPResidue = PVAI(uUnit->vaResidues, SAVERESIDUEt, 0);
        for (i = 0; i < iResidueCount; srPResidue++, i++) {
            srPResidue->iImagingAtomIndex = 0;
            if (aResidueImagingAtom(srPResidue->rResidue) != NULL) {
                srPResidue->iImagingAtomIndex =
                    iContainerTempInt(aResidueImagingAtom
                                      (srPResidue->rResidue));
            }
            for (j = 0; j < MAXCONNECT; j++) {
                if (bResidueConnectUsed(srPResidue->rResidue, j)) {
                    srPResidue->iaConnectIndex[j] =
                        iContainerTempInt(srPResidue->rResidue->
                                          aaConnect[j]);
                } else
                    srPResidue->iaConnectIndex[j] = 0;
            }
        }
    }

    /* Now generate the connectivity table */

    iCount = 0;
    lTemp = lLoop((OBJEKT) uUnit, BONDS);
    while (oNext(&lTemp) != NULL)
        iCount++;
    uUnit->vaConnectivity = vaVarArrayCreate(sizeof(SAVECONNECTIVITYt));
    VarArraySetSize((uUnit->vaConnectivity), iCount);
    if (iCount) {
        i = 0;
        lTemp = lLoop((OBJEKT) uUnit, BONDS);
        scPCon = PVAI(uUnit->vaConnectivity, SAVECONNECTIVITYt, 0);
        while (oNext(&lTemp) != NULL) {
            LoopGetBond(&lTemp, &aAtom1, &aAtom2);
            scPCon->fFlags = fAtomFindBondFlags(aAtom1, aAtom2);
            scPCon->iAtom1 = iContainerTempInt(aAtom1);
            scPCon->iAtom2 = iContainerTempInt(aAtom2);
            i++;
            scPCon++;
        }
    }


    /* Build a table for the Hierarchy information */

    iCount = 0;
    lTemp = lLoop((OBJEKT) uUnit, CONTAINERS);
    while (oNext(&lTemp) != NULL)
        iCount++;
    if (iCount) {
        uUnit->vaHierarchy = vaVarArrayCreate(sizeof(SAVEHIERARCHYt));
        VarArraySetSize((uUnit->vaHierarchy), iCount);
        lTemp = lLoop((OBJEKT) uUnit, CONTAINERS);
        i = 0;
        shPHierarchy = PVAI(uUnit->vaHierarchy, SAVEHIERARCHYt, 0);
        while ((oBelow = oNext(&lTemp)) != NULL) {
            oAbove = (OBJEKT) cContainerWithin(oBelow);
            shPHierarchy->sAboveType[0] = iObjectType(oAbove);
            shPHierarchy->sAboveType[1] = '\0';
            shPHierarchy->sBelowType[0] = iObjectType(oBelow);
            shPHierarchy->sBelowType[1] = '\0';
            if (iObjectType(oAbove) == UNITid) {
                shPHierarchy->iAboveIndex = 0;
            } else {
                shPHierarchy->iAboveIndex = iContainerTempInt(oAbove);
            }
            shPHierarchy->iBelowIndex = iContainerTempInt(oBelow);
            i++;
            shPHierarchy++;
        }
    }

    /* Build a table for the UNIT Connect information */

    uUnit->vaConnect = vaVarArrayCreate(sizeof(int));
    VarArraySetSize((uUnit->vaConnect), 2);
    if (bUnitHeadUsed(uUnit))
        iIndex = iContainerTempInt(aUnitHead(uUnit));
    else
        iIndex = 0;
    *PVAI(uUnit->vaConnect, int, 0) = iIndex;
    if (bUnitTailUsed(uUnit))
        iIndex = iContainerTempInt(aUnitTail(uUnit));
    else
        iIndex = 0;
    *PVAI(uUnit->vaConnect, int, 1) = iIndex;



    /* Build tables for the UNIT groups */

    if (iDictionaryElementCount(uUnit->dAtomGroups)) {
        uUnit->vaGroupNames = vaVarArrayCreate(sizeof(STRING));
        VarArraySetSize(uUnit->vaGroupNames,
                        iDictionaryElementCount(uUnit->dAtomGroups));

        /* assuming there are atoms; the array is grown on the fly */
        uUnit->vaGroupAtoms = vaVarArrayCreate(sizeof(SAVEGROUPSt));
        iGroup = 0;
        dlGroups = ydlDictionaryLoop(uUnit->dAtomGroups);
        while ((lGroup =
               (LIST) yPDictionaryNext(uUnit->dAtomGroups, &dlGroups))) {
            strcpy((char *) (PVAI(uUnit->vaGroupNames, STRING, iGroup)),
                   sDictLoopKey(dlGroups));
            llAtoms = llListLoop(lGroup);
            for (i = 0; i < iListSize(lGroup); i++) {
                sgAtomGroup.iGroupIndex = iGroup + 1;
                sgAtomGroup.iIndexAtom =
                    iContainerTempInt((ATOM) oListNext(&llAtoms));
                VarArrayAdd(uUnit->vaGroupAtoms, (GENP) & sgAtomGroup);
            }
            iGroup++;
        }
    }

/* The rest of the code should be executed ONLY if PARAMETERS are being */
/* generated */

    if (bFailedGeneratingParameters == FALSE && bGenerateParameters) {


        /* Now generate the BOND table */

        bFailedGeneratingParameters |=
            zbUnitIOIndexBondParameters(plParameters, uUnit, bPert);

	/* New now generate C4 table  */
	bFailedGeneratingParameters |=
            zbUnitIOIndexC4Pairwise(uUnit);

        /* Now generate the ANGLE table */

        bFailedGeneratingParameters |=
            zbUnitIOIndexAngleParameters(plParameters, uUnit, bPert);


        /* Now generate the TORSION and IMPROPER tables */
        /* Place them in the SAME VARARRAY !!!!! */

        if (uUnit->vaTorsions != NULL) {
            VP0("Regenerating proper and improper torsions.\n");
            VarArrayDestroy(&(uUnit->vaTorsions));
        }
        uUnit->vaTorsions = vaVarArrayCreate(sizeof(SAVETORSIONt));
        bFailedGeneratingParameters |=
            zbUnitIOIndexTorsionParameters(plParameters, uUnit, TRUE, bPert);
        bFailedGeneratingParameters |=
            zbUnitIOIndexTorsionParameters(plParameters, uUnit, FALSE, bPert);


        /* Generate the potential H-Bond parameters             */
        /* Do this by looping through ALL H-Bond parameters     */
        /* looking for those where both atoms are defined       */
        /* within this UNITs atom list                          */
        /* If there is one then add it to this UNITs PARMSET    */

        VP0("Building H-Bond parameters.\n");
        PARMLIB_LOOP_ALL(plParameters, psTemp) {
            for (i = 0; i < iParmSetTotalHBondParms(psTemp); i++) {
                ParmSetHBond(psTemp, i, sAtom1, sAtom2, &dA, &dB, sDesc);
                if ((iParmSetFindAtom(uUnit->psParameters, sAtom1)
                     != PARM_NOT_FOUND) &&
                    (iParmSetFindAtom(uUnit->psParameters, sAtom2)
                     != PARM_NOT_FOUND)) {
                    if (iParmSetFindHBond
                        (uUnit->psParameters, sAtom1,
                         sAtom2) == PARM_NOT_FOUND) {
                        iParmSetAddHBond(uUnit->psParameters, sAtom1,
                                         sAtom2, dA, dB, sDesc);
                    }
                }
            }
        }

	/* Carry LJ edits into the unit topology */

	VP0("Incorporating Non-Bonded adjustments.\n");
	PARMLIB_LOOP_ALL(plParameters, psTemp) {
	  for (i = 0; i < iParmSetTotalNBEdits(psTemp); i++) {
	    ParmSetNBEdit(psTemp, i, sAtom1, sAtom2, &dEI, &dEJ, &dRI, &dRJ,
			  sDesc);
	    if ((iParmSetFindAtom(uUnit->psParameters, sAtom1)
		 != PARM_NOT_FOUND) &&
		(iParmSetFindAtom(uUnit->psParameters, sAtom2)
		 != PARM_NOT_FOUND)) {
	      if (iParmSetFindNBEdit(uUnit->psParameters, sAtom1, sAtom2) ==
		  PARM_NOT_FOUND) {
		iParmSetAddNBEdit(uUnit->psParameters, sAtom1, sAtom2, dEI,
				  dEJ, dRI, dRJ, sDesc);
	      }
	    }
	  }
	}
    }
    if (bFailedGeneratingParameters == FALSE)
        *bPGeneratedParameters = TRUE;
    else
        *bPGeneratedParameters = FALSE;
}






/*
 *      zUnitIOBuildFromTables
 *
 *        Author:        Christian Schafmeister (1991)
 *
 *      Build a UNIT from its tables.
 */
void zUnitIOBuildFromTables(UNIT uUnit)
{
    SAVEATOMt *saPAtom;
    SAVECONNECTIVITYt *scPCon;
    SAVERESTRAINTt *srPRestraint;
    SAVERESIDUEt *srPResidue;
    SAVEMOLECULEt *smPMolecule;
    SAVEHIERARCHYt *shPHierarchy;
    SAVEGROUPSt *sgPGroupAtom;
    RESTRAINT rRest;
    ATOM aAtom, aAtom1, aAtom2, aAtom3, aAtom4;
    RESIDUE rRes;
    MOLECULE mMol;
    OBJEKT oObj1, oObj2;
    int i, j, iMax, iIndex, iChildSeq, iSeq, iGroup;
    STRING sGroup;
    LIST lGroup;

    /* Build the atoms */

    iMax = iVarArrayElementCount(uUnit->vaAtoms);
    if (!iMax)
        return;

    saPAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, 0);
    for (i = 0; i < iMax; i++, saPAtom++) {
        aAtom = (ATOM) oCreate(ATOMid);
        ContainerSetName(aAtom, saPAtom->sName);
        ContainerSetSequence(aAtom, saPAtom->iSequence);
        AtomSetElement(aAtom, saPAtom->iElement);
        AtomSetPertElement(aAtom, saPAtom->iPertElement);
        AtomSetPertName(aAtom, saPAtom->sPertName);
        AtomSetType(aAtom, saPAtom->sType);
        AtomSetPertType(aAtom, saPAtom->sPertType);
        AtomSetCharge(aAtom, saPAtom->dCharge);
        AtomSetPertCharge(aAtom, saPAtom->dPertCharge);
        AtomSetPosition(aAtom, saPAtom->vPos);
        AtomSetVelocity(aAtom, saPAtom->vVelocity);
        AtomDefineFlags(aAtom, saPAtom->fFlags);
        /*
         *  the ATOM objekt has 1 reference at this point..
         */
        saPAtom->aAtom = aAtom;
    }

    /* Build the Residues */

    srPResidue = PVAI(uUnit->vaResidues, SAVERESIDUEt, 0);
    iMax = iVarArrayElementCount(uUnit->vaResidues);
    rRes = NULL;
    for (i = 0; i < iMax; i++, srPResidue++) {
        rRes = (RESIDUE) oCreate(RESIDUEid);
        ContainerSetName(rRes, srPResidue->sName);
        ContainerSetSequence(rRes, srPResidue->iSequenceNumber);
        ContainerSetNextChildsSequence(rRes,
                                       srPResidue->iNextChildSequence);
        ResidueSetPdbSequence(rRes, srPResidue->iPdbResSeq);
        ResidueSetChainId(rRes, srPResidue->sChainId);
        ResidueSetType(rRes, srPResidue->sResidueType[0]);

        /* Define the imaging ATOM */

        iIndex = srPResidue->iImagingAtomIndex;
        if (iIndex != 0)
            aAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, iIndex - 1)->aAtom;
        else
            aAtom = NULL;
        ResidueSetImagingAtom(rRes, aAtom);

        for (j = 0; j < MAXCONNECT; j++) {
            iIndex = srPResidue->iaConnectIndex[j];
            if (iIndex != 0) {
                aAtom = (PVAI(uUnit->vaAtoms, SAVEATOMt,
                              iIndex - 1))->aAtom;
                ResidueSetConnectAtom(rRes, j, aAtom);
            }
        }
        /*
         *  the RESIDUE objekt has 1 reference at this point..
         */
        srPResidue->rResidue = rRes;
    }

    /* Build the Molecules */

    if ((iMax = iVarArrayElementCount(uUnit->vaMolecules))) {
        smPMolecule = PVAI(uUnit->vaMolecules, SAVEMOLECULEt, 0);
        for (i = 0; i < iMax; i++, smPMolecule++) {
            mMol = (MOLECULE) oCreate(MOLECULEid);
            ContainerSetName(mMol, smPMolecule->sName);
            ContainerSetSequence(mMol, smPMolecule->iSequenceNumber);
            // FIXME: added rRes=NULL initialization, may be used here; but is this dead code? --JMK
            ContainerSetNextChildsSequence(rRes,
                                           smPMolecule->
                                           iNextChildSequence);
            smPMolecule->mMolecule = mMol;
        }
    }

    /* Create bonds between atoms */

    if ((iMax = iVarArrayElementCount(uUnit->vaConnectivity))) {
        for (i = 0; i < iMax; i++) {
            scPCon = PVAI(uUnit->vaConnectivity, SAVECONNECTIVITYt, i);
            aAtom1 =
                PVAI(uUnit->vaAtoms, SAVEATOMt, scPCon->iAtom1 - 1)->aAtom;
            aAtom2 =
                PVAI(uUnit->vaAtoms, SAVEATOMt, scPCon->iAtom2 - 1)->aAtom;
            AtomBondToFlags(aAtom1, aAtom2, scPCon->fFlags);
        }
    }

    /* Build the hierarchy of Unit/Molecules/Residues/Atoms */

    if ((iMax = iVarArrayElementCount(uUnit->vaHierarchy))) {
        shPHierarchy = PVAI(uUnit->vaHierarchy, SAVEHIERARCHYt, 0);
        for (i = 0; i < iMax; i++, shPHierarchy++) {

            switch (shPHierarchy->sAboveType[0]) {
            case UNITid:
                oObj1 = (OBJEKT) uUnit;
                break;
            case MOLECULEid:
                oObj1 = (OBJEKT)
                    (PVAI(uUnit->vaMolecules, SAVEMOLECULEt,
                          shPHierarchy->iAboveIndex - 1))->mMolecule;
                break;
            case RESIDUEid:
                oObj1 = (OBJEKT)
                    (PVAI(uUnit->vaResidues, SAVERESIDUEt,
                          shPHierarchy->iAboveIndex - 1))->rResidue;
                break;
            case ATOMid:
                oObj1 = (OBJEKT)
                    (PVAI(uUnit->vaAtoms, SAVEATOMt,
                          shPHierarchy->iAboveIndex - 1))->aAtom;
                break;
            default:
                DFATAL("Coding error, object hierarchy");
            }
            switch (shPHierarchy->sBelowType[0]) {
            case UNITid:
                oObj2 = (OBJEKT) uUnit;
                break;
            case MOLECULEid:
                oObj2 = (OBJEKT)
                    (PVAI(uUnit->vaMolecules, SAVEMOLECULEt,
                          shPHierarchy->iBelowIndex - 1))->mMolecule;
                break;
            case RESIDUEid:
                oObj2 = (OBJEKT)
                    (PVAI(uUnit->vaResidues, SAVERESIDUEt,
                          shPHierarchy->iBelowIndex - 1))->rResidue;
                break;
            case ATOMid:
                oObj2 = (OBJEKT)
                    (PVAI(uUnit->vaAtoms, SAVEATOMt,
                          shPHierarchy->iBelowIndex - 1))->aAtom;
                break;
            default:
                DFATAL("Coding error, object hierarchy");
            }
            iChildSeq = iContainerNextChildsSequence(oObj1);
            iSeq = iContainerSequence(oObj2);
            /*
             *  ContainerAdd() increments Obj2's reference counter
             */
            ContainerAdd((CONTAINER) oObj1, oObj2);
            ContainerSetNextChildsSequence(oObj1, iChildSeq);
            ContainerSetSequence(oObj2, iSeq);
        }
    }

    /* Set up the RESTRAINTs */

    if ((iMax = iVarArrayElementCount(uUnit->vaRestraints))) {
        srPRestraint = PVAI(uUnit->vaRestraints, SAVERESTRAINTt, 0);
        for (i = 0; i < iMax; i++, srPRestraint++) {
            rRest = rRestraintCreate();
            RestraintDefineFlags(rRest, srPRestraint->fFlags);
            switch (srPRestraint->iType) {
            case RESTRAINTBOND:
                aAtom1 = PVAI(uUnit->vaAtoms, SAVEATOMt,
                              srPRestraint->iAtom1 - 1)->aAtom;
                aAtom2 = PVAI(uUnit->vaAtoms, SAVEATOMt,
                              srPRestraint->iAtom2 - 1)->aAtom;
                RestraintBondSet(rRest, aAtom1, aAtom2,
                                 srPRestraint->dKx, srPRestraint->dX0);
                break;
            case RESTRAINTANGLE:
                aAtom1 = PVAI(uUnit->vaAtoms, SAVEATOMt,
                              srPRestraint->iAtom1 - 1)->aAtom;
                aAtom2 = PVAI(uUnit->vaAtoms, SAVEATOMt,
                              srPRestraint->iAtom2 - 1)->aAtom;
                aAtom3 = PVAI(uUnit->vaAtoms, SAVEATOMt,
                              srPRestraint->iAtom3 - 1)->aAtom;
                RestraintAngleSet(rRest, aAtom1, aAtom2, aAtom3,
                                  srPRestraint->dKx, srPRestraint->dX0);
                break;
            case RESTRAINTTORSION:
                aAtom1 = PVAI(uUnit->vaAtoms, SAVEATOMt,
                              srPRestraint->iAtom1 - 1)->aAtom;
                aAtom2 = PVAI(uUnit->vaAtoms, SAVEATOMt,
                              srPRestraint->iAtom2 - 1)->aAtom;
                aAtom3 = PVAI(uUnit->vaAtoms, SAVEATOMt,
                              srPRestraint->iAtom3 - 1)->aAtom;
                aAtom4 = PVAI(uUnit->vaAtoms, SAVEATOMt,
                              srPRestraint->iAtom4 - 1)->aAtom;
                RestraintTorsionSet(rRest, aAtom1, aAtom2, aAtom3, aAtom4,
                                    srPRestraint->dKx, srPRestraint->dX0,
                                    srPRestraint->dN);
                break;
            default:
                DFATAL("Invalid RESTRAINT type loaded");
            }
            UnitAddRestraint(uUnit, rRest);
        }
    }

    /* Set up the UNIT Connect information */

    iIndex = *PVAI(uUnit->vaConnect, int, 0);
    aAtom = NULL;
    if (iIndex != 0)
        aAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, iIndex - 1)->aAtom;
    UnitSetHead(uUnit, aAtom);

    iIndex = *PVAI(uUnit->vaConnect, int, 1);
    aAtom = NULL;
    if (iIndex != 0)
        aAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, iIndex - 1)->aAtom;
    UnitSetTail(uUnit, aAtom);

    /* Set up the ATOM groups information */

    if ((iMax = iVarArrayElementCount(uUnit->vaGroupNames))) {
        for (i = 0; i < iMax; i++) {
            strcpy(sGroup,
                   (char *) (PVAI(uUnit->vaGroupNames, STRING, i)));
            bUnitGroupCreate(uUnit, sGroup);
        }
    }
    if ((iMax = iVarArrayElementCount(uUnit->vaGroupAtoms))) {
        sgPGroupAtom = PVAI(uUnit->vaGroupAtoms, SAVEGROUPSt, 0);
        for (i = 0; i < iMax; i++, sgPGroupAtom++) {
            iGroup = sgPGroupAtom->iGroupIndex - 1;
            iIndex = sgPGroupAtom->iIndexAtom;
            lGroup = lUnitGroup(uUnit,
                                (char *) PVAI(uUnit->vaGroupNames, STRING,
                                              iGroup));
            aAtom = PVAI(uUnit->vaAtoms, SAVEATOMt, iIndex - 1)->aAtom;
            ListAddUnique(lGroup, (GENP) aAtom);
        }
    }

}

/*
 *      zUnitIODestroyTables
 *
 *        Author:        Christian Schafmeister (1991)
 *
 *      Destroy all of the tables associated with the UNIT.
 */
void zUnitIODestroyTables(UNIT uUnit)
{
    int i, iCount;

    /*
     *  since the SAVEATOMt's have pointers to the ATOMs,
     *      need to deref the ATOMs so that they will
     *      be freed when the last pointer to them is
     *      deleted. Same goes for RESIDUEs..
     */
    if (uUnit->vaAtoms != NULL) {
        SAVEATOMt *aP = PVAI(uUnit->vaAtoms, SAVEATOMt, 0);
        iCount = iVarArrayElementCount(uUnit->vaAtoms);
        for (i = 0; i < iCount; i++, aP++)
            DEREF(aP->aAtom);
        VarArrayDestroy(&(uUnit->vaAtoms));
    }
    if (uUnit->vaResidues != NULL) {
        SAVERESIDUEt *rP = PVAI(uUnit->vaResidues, SAVERESIDUEt, 0);
        iCount = iVarArrayElementCount(uUnit->vaResidues);
        for (i = 0; i < iCount; i++, rP++)
            DEREF(rP->rResidue);
        VarArrayDestroy(&(uUnit->vaResidues));
    }

    if (uUnit->vaMolecules != NULL) {
        SAVEMOLECULEt *mP = PVAI(uUnit->vaMolecules, SAVEMOLECULEt, 0);
        iCount = iVarArrayElementCount(uUnit->vaMolecules);
        for (i = 0; i < iCount; i++, mP++)
            DEREF(mP->mMolecule);
        VarArrayDestroy(&(uUnit->vaMolecules));
    }
/* %%% */
    if (uUnit->vaHierarchy != NULL)
        VarArrayDestroy(&(uUnit->vaHierarchy));
    if (uUnit->vaBonds != NULL)
        VarArrayDestroy(&(uUnit->vaBonds));
    if (uUnit->vaC4Pairwise != NULL)
	    VarArrayDestroy(&(uUnit->vaC4Pairwise)); //New
    if (uUnit->vaAngles != NULL)
        VarArrayDestroy(&(uUnit->vaAngles));
    if (uUnit->vaTorsions != NULL)
        VarArrayDestroy(&(uUnit->vaTorsions));
    if (uUnit->vaRestraints != NULL)
        VarArrayDestroy(&(uUnit->vaRestraints));
    if (uUnit->vaConnect != NULL)
        VarArrayDestroy(&(uUnit->vaConnect));
    if (uUnit->vaConnectivity != NULL)
        VarArrayDestroy(&(uUnit->vaConnectivity));
    if (uUnit->vaGroupNames != NULL)
        VarArrayDestroy(&(uUnit->vaGroupNames));
    if (uUnit->vaGroupAtoms != NULL)
        VarArrayDestroy(&(uUnit->vaGroupAtoms));
    if (uUnit->vaAtomsPerMolecule != NULL)
        VarArrayDestroy(&(uUnit->vaAtomsPerMolecule));

    UnitSetMode(uUnit, UNITNORMAL);
}

static char *prepfmt = "   %-3d %-4s  %-4s  %c      %f  %f  %f    %f\n";

/*
 *  zPrintSideChain() - depth-first descent of side chain,
 *        printing atoms, noting loops
 */
static void
zPrintSideChain(FILE * fOut, ATOM aParentAtom, ATOM aAtom,
               int *iP, VARARRAY vaLoopAtoms)
{
    VECTOR vPos;
    int i, j, k, iChildTag = *iP;
    ATOM aChildAtom, aNbrs[MAXBONDS];

    /*
     *  figure tree type - count 'downstream',
     *      undesignated neighbors, i.e. omits
     *      parent and anything else already
     *      encountered. Mark all claimed atoms
     *      as belonging to this one for the
     *      benefit of recursion looping back
     *      around to this point before getting
     *      to the atom later in this routine
     */
/*
fprintf(stderr, "-- parent %s -> %s coord %d\n",
sAtomName(aParentAtom),sAtomName(aAtom), iAtomCoordination(aAtom));
*/
    j = 0;
    for (i = 0; i < iAtomCoordination(aAtom); i++) {
        aChildAtom = aAtomBondedNeighbor(aAtom, i);
        if (aChildAtom == aParentAtom)
            continue;
/*
fprintf(stderr, "   child %s type %c\n",
sAtomName(aChildAtom), (char)dAtomTemp( aChildAtom ));
*/
        if (dAtomTemp(aChildAtom) == (double) 'x' &&
            iAtomTempInt(aChildAtom) == -1) {
            AtomSetTempInt(aChildAtom, iChildTag);
            j++;
        }
    }
    zSetTreeType(aAtom, j);

    /*
     *  print
     */
    vPos = vAtomPosition(aAtom);
    fprintf(fOut, prepfmt, *iP,
            sAtomName(aAtom),
            sAtomType(aAtom),
            (char) dAtomTemp(aAtom),
            vPos.dX, vPos.dY, vPos.dZ, dAtomCharge(aAtom));
    AtomSetSeenId(aAtom, *iP);
    (*iP)++;

    /*
     *  put eligible children in order for readability -
     *      get obvious 'E' types 1st
     */
    k = 0;
    for (j = 0; j < iAtomCoordination(aAtom); j++) {
        aChildAtom = aAtomBondedNeighbor(aAtom, j);
        if (aChildAtom == aParentAtom)
            continue;
        if (dAtomTemp(aChildAtom) != (double) 'x')
            continue;
        if (iAtomCoordination(aChildAtom) == 1)
            aNbrs[k++] = aChildAtom;
    }
    /*
     *  bubble sort them into alphabetical order
     */
    while (1) {
        int iMore = 0;

        for (j = 0; j < k - 1; j++) {
            if (strcmp(sAtomName(aNbrs[j]), sAtomName(aNbrs[j + 1])) > 0) {
                ATOM aTmp;
                aTmp = aNbrs[j + 1];
                aNbrs[j + 1] = aNbrs[j];
                aNbrs[j] = aTmp;
                iMore++;
            }
        }
        if (!iMore)
            break;
    }
    /*
     *  print
     */
    for (j = 0; j < k; j++)
        zPrintSideChain(fOut, aAtom, aNbrs[j], iP, vaLoopAtoms);

    /*
     *  get remaining eligible children
     */
    k = 0;
    for (j = 0; j < iAtomCoordination(aAtom); j++) {
        aChildAtom = aAtomBondedNeighbor(aAtom, j);
        if (aChildAtom == aParentAtom)
            continue;
        if (iAtomCoordination(aChildAtom) == 1)        /* done above */
            continue;

        if (dAtomTemp(aChildAtom) == (double) 'x' &&
            iAtomTempInt(aChildAtom) == iChildTag) {
            aNbrs[k++] = aChildAtom;
        } else {
            STRING sTemp;

            /*
             *  atom seen already: may have found a new
             *      loop closing bond
             */
            if (strcmp(sAtomName(aAtom), sAtomName(aChildAtom)) < 0)
                sprintf(sTemp, "%-4s %-4s",
                        sAtomName(aAtom), sAtomName(aChildAtom));
            else
                sprintf(sTemp, "%-4s %-4s",
                        sAtomName(aChildAtom), sAtomName(aAtom));
            for (i = 0; i < iVarArrayElementCount(vaLoopAtoms); i++) {
                if (!strcmp(sTemp, PVAI(vaLoopAtoms, char, i)))
                     break;
            }
            if (i == iVarArrayElementCount(vaLoopAtoms))
                VarArrayAdd(vaLoopAtoms, (GENP) sTemp);
        }
    }
    /*
     *  bubble sort
     */
    while (1) {
        int iMore = 0;

        for (j = 0; j < k - 1; j++) {
            if (strcmp(sAtomName(aNbrs[j]), sAtomName(aNbrs[j + 1])) > 0) {
                ATOM aTmp;
                aTmp = aNbrs[j + 1];
                aNbrs[j + 1] = aNbrs[j];
                aNbrs[j] = aTmp;
                iMore++;
            }
        }
        if (!iMore)
            break;
    }
    /*
     *  print
     */
    for (j = 0; j < k; j++)
        zPrintSideChain(fOut, aAtom, aNbrs[j], iP, vaLoopAtoms);
}

/*
 *  zMarkSideChain() - depth-first descent of side chain
 */
static void
zMarkSideChain(ATOM aParentAtom, ATOM aAtom, int *iP)
{
    int i, j, k, iChildTag = *iP;
    ATOM aChildAtom, aNbrs[MAXBONDS];

    /*
     *  figure tree type - count 'downstream',
     *      undesignated neighbors, i.e. omits
     *      parent and anything else already
     *      encountered. Mark all claimed atoms
     *      as belonging to this one for the
     *      benefit of recursion looping back
     *      around to this point before getting
     *      to the atom later in this routine
     */
    j = 0;
    for (i = 0; i < iAtomCoordination(aAtom); i++) {
        aChildAtom = aAtomBondedNeighbor(aAtom, i);
        if (aChildAtom == aParentAtom)
            continue;
        if (dAtomTemp(aChildAtom) == (double) 'x' &&
            iAtomTempInt(aChildAtom) == -1) {
            AtomSetTempInt(aChildAtom, iChildTag);
            j++;
        }
    }
    zSetTreeType(aAtom, j);

    AtomSetSeenId(aAtom, *iP);
    (*iP)++;

    /*
     *  put eligible children in order for readability -
     *      get obvious 'E' types 1st
     */
    k = 0;
    for (j = 0; j < iAtomCoordination(aAtom); j++) {
        aChildAtom = aAtomBondedNeighbor(aAtom, j);
        if (aChildAtom == aParentAtom)
            continue;
        if (dAtomTemp(aChildAtom) != (double) 'x')
            continue;
        if (iAtomCoordination(aChildAtom) == 1)
            aNbrs[k++] = aChildAtom;
    }
    /*
     *  bubble sort them into alphabetical order
     */
    while (1) {
        int iMore = 0;

        for (j = 0; j < k - 1; j++) {
            if (strcmp(sAtomName(aNbrs[j]), sAtomName(aNbrs[j + 1])) > 0) {
                ATOM aTmp;
                aTmp = aNbrs[j + 1];
                aNbrs[j + 1] = aNbrs[j];
                aNbrs[j] = aTmp;
                iMore++;
            }
        }
        if (!iMore)
            break;
    }
    /*
     *  print
     */
    for (j = 0; j < k; j++)
        zMarkSideChain(aAtom, aNbrs[j], iP);

    /*
     *  get remaining eligible children
     */
    k = 0;
    for (j = 0; j < iAtomCoordination(aAtom); j++) {
        aChildAtom = aAtomBondedNeighbor(aAtom, j);
        if (aChildAtom == aParentAtom)
            continue;
        if (iAtomCoordination(aChildAtom) == 1)        /* done above */
            continue;

        if (dAtomTemp(aChildAtom) == (double) 'x' &&
            iAtomTempInt(aChildAtom) == iChildTag) {
            aNbrs[k++] = aChildAtom;
        }
    }
    /*
     *  bubble sort
     */
    while (1) {
        int iMore = 0;

        for (j = 0; j < k - 1; j++) {
            if (strcmp(sAtomName(aNbrs[j]), sAtomName(aNbrs[j + 1])) > 0) {
                ATOM aTmp;
                aTmp = aNbrs[j + 1];
                aNbrs[j + 1] = aNbrs[j];
                aNbrs[j] = aTmp;
                iMore++;
            }
        }
        if (!iMore)
            break;
    }
    /*
     *  print
     */
    for (j = 0; j < k; j++)
        zMarkSideChain(aAtom, aNbrs[j], iP);
}

int
iMarkMainChainAtoms(RESIDUE rRes, int complain)
{
    int iAtomCount, i, j, iLevel, iMin, iMax, iNext, ierr = 0;
    ATOM aAtom, aAtom0, aAtom1, aChildAtom;
    VARARRAY vaAtoms;
    LOOP lTemp;
    char *cPResName;
    cPResName = sContainerName(rRes);

    /*
     *  set up for breadth-1st search: count atoms,
     *      setting up array of atom pointers & initializing depth
     *  NOTE: 'x' will become 'BLA' in the PRMTOP, and will happen
     *  for any residue without CONNECT0 and CONNECT1
     */
    iAtomCount = 0;
    lTemp = lLoop((OBJEKT) rRes, ATOMS);
    while ((aAtom = (ATOM) oNext(&lTemp)) != NULL) {
        AtomSetTempInt(aAtom, -1);
        AtomSetTempDouble(aAtom, (double) 'x');
        iAtomCount++;
    }
    if (!iAtomCount) {
        VP0("  %s: no atoms\n", cPResName);
        return (0);
    }
    if (iAtomCount == 1) {
        /*
         *  call it a main chain whether connected or not
         */
        lTemp = lLoop((OBJEKT) rRes, ATOMS);
        aAtom = (ATOM) oNext(&lTemp);
        AtomSetTempDouble(aAtom, (double) 'M');
        return (1);
    }

    /*
     *  check connect atoms
     */
    aAtom0 = (ATOM) rRes->aaConnect[0];
    aAtom1 = (ATOM) rRes->aaConnect[1];
    if (aAtom0 == NULL) {
        if (complain)
            VP0("  %s:  connect0 not defined\n", cPResName);
        ierr++;
    }
    if (aAtom1 == NULL) {
        if (complain)
            VP0("  %s:  connect1 not defined\n", cPResName);
        ierr++;
    }
    if (ierr)
        return (-iAtomCount);
/*
fprintf(stderr," connect %s .. %s total %d\n",
sAtomName( aAtom0 ), sAtomName( aAtom1), iAtomCount);
*/
    /*
     *  prepare 'stack' for tracking depth-1st search
     */
    vaAtoms = vaVarArrayCreate(sizeof(ATOM));
    VarArraySetSize(vaAtoms, iAtomCount + 1);
    *PVAI(vaAtoms, ATOM, 0) = aAtom0;
    iLevel = 0;
    iMin = 0;
    iMax = 1;

    /*
     *  loop over each level in tree starting at aAtom0
     *      until target atom found, marking atoms w/ tree level
     */
    while (1) {
        /*
         *  prepare to accumulate next level down
         */
        iNext = iMax;
        /*
         *  loop over atoms in  current level
         */
/*
fprintf(stderr,"--- lev %d:  %d .. %d\n", iLevel, iMin, iMax);
*/
        for (i = iMin; i < iMax; i++) {
            /*
             *  1st time, i=0 and 0th elt of vaAtoms is aAtom0
             */
            aAtom = *PVAI(vaAtoms, ATOM, i);
/*
fprintf(stderr,"     down %s lev %d coord %d\n",
sAtomName(aAtom), iLevel, iAtomCoordination(aAtom));
*/
            /*
             *  mark level and see if this is target
             */
            AtomSetTempInt(aAtom, iLevel);
            if (aAtom == aAtom1)
                goto FOUND;

            /*
             *  put all unseen 'children' on
             *      next level stack
             */
            for (j = 0; j < iAtomCoordination(aAtom); j++) {
                aChildAtom = aAtomBondedNeighbor(aAtom, j);
/*
fprintf(stderr,"         child %s add %c\n",
sAtomName(aChildAtom),
(iAtomTempInt( aChildAtom ) == -1 ? 'Y' : 'N'));
*/
                /*
                 *  stay within residue
                 */
                if (cContainerWithin(aChildAtom) !=
                    cContainerWithin(aAtom)) continue;
                if (iAtomTempInt(aChildAtom) == -1) {
                    *PVAI(vaAtoms, ATOM, iNext++) = aChildAtom;
                    /*
                     *  re-flag atom to ensure
                     *      that it doesn't get included
                     *      again as a child of another
                     *      atom at this level
                     */
                    AtomSetTempInt(aChildAtom, -2);
                }
            }

        }
        iMin = iMax;
        iMax = iNext;
        iLevel++;
    }
  FOUND:
    /*
     *  Traverse marked tree from connect1 back up to
     *      connect0, marking main chain
     */
    aAtom = aAtom1;
    while (1) {
        ATOM aSelected = (ATOM) NULL;

        AtomSetTempDouble(aAtom, (double) 'M');
        if (aAtom == aAtom0)
            break;
        for (j = 0; j < iAtomCoordination(aAtom); j++) {
            aChildAtom = aAtomBondedNeighbor(aAtom, j);
            /*
             *  stay within residue
             */
            if (cContainerWithin(aChildAtom) != cContainerWithin(aAtom))
                continue;
            if (iAtomTempInt(aChildAtom) == iLevel - 1) {
                /*
                 *  choose lesser atom if >1, for readability
                 */
                if (aSelected == NULL) {
                    aSelected = aChildAtom;
                } else if (strcmp(sAtomName(aChildAtom),
                                  sAtomName(aSelected)) < 0) {
                    aSelected = aChildAtom;
                }
            }
        }
        aAtom = aSelected;
        iLevel--;
    }
    /*
     *  clean up
     */
    VarArrayDestroy(&vaAtoms);

    return (iAtomCount);
}

/*
 *  zMarkSideChains() - used for prmtop
 *
 *        NOTE: ziMarkMainChainAtoms() must be run 1st
 */
void
MarkSideChains(RESIDUE rRes)
{
    ATOM aAtom, aParentAtom, aAtom0, aAtom1;
    LOOP lTemp;
    int iCount;

    aAtom0 = (ATOM) rRes->aaConnect[0];
    aAtom1 = (ATOM) rRes->aaConnect[1];
    if (aAtom0 == NULL)                /* 1-atom residue */
        return;

    /*
     *  traverse tree from top, following main chain
     *
     *      reset atom temp ints to use in marking 'seen' atoms
     */
    lTemp = lLoop((OBJEKT) rRes, ATOMS);
    while ((aAtom = (ATOM) oNext(&lTemp)) != NULL) {
        AtomSetTempInt(aAtom, -1);
        AtomSetSeenId(aAtom, -1);
    }

    aAtom = aAtom0;
    aParentAtom = NULL;
    iCount = 0;
    while (1) {
        ATOM aNextMain;
        int j;

        AtomSetSeenId(aAtom, iCount);
        /*
         *  find next main chain down in this residue,
         *      marking side chains
         */
        aNextMain = NULL;
        for (j = 0; j < iAtomCoordination(aAtom); j++) {
            ATOM aChildAtom = aAtomBondedNeighbor(aAtom, j);
            if (aChildAtom == aParentAtom)
                continue;
            if (cContainerWithin(aChildAtom) != cContainerWithin(aAtom))
                continue;
            if (dAtomTemp(aChildAtom) == (double) 'M') {
                aNextMain = aChildAtom;
            } else if (dAtomTemp(aChildAtom) == (double) 'x') {
                /*
                 *  side chain not seen before
                 */
                zMarkSideChain(aAtom, aChildAtom, &iCount);
            }
        }
        if (aAtom == aAtom1)
            break;
        if (aNextMain == NULL)        /* 1 'M' atom in residue */
            break;
        aParentAtom = aAtom;
        aAtom = aNextMain;
    }
}




static void WritePrepRes(RESIDUE rRes, FILE * fOut)
{
    int i, j, iMax;
    LOOP lTemp;
    ATOM aAtom, aAtom0, aAtom1, aChildAtom, aParentAtom;
    VARARRAY vaLoopAtoms;
    char *cPResName;

    cPResName = sContainerName(rRes);
    VP0("  saving prep, residue %s\n", cPResName);

    /*
     *  mark main chain atoms; worry about side chains later
     */

    if (iMarkMainChainAtoms(rRes, 1) < 1)
        return;

    aAtom0 = (ATOM) rRes->aaConnect[0];
    aAtom1 = (ATOM) rRes->aaConnect[1];


    /*
     *  write residue header & dummy atoms
     */
    if (strlen(cPResName) > 3)
        cPResName += strlen(cPResName) - 3;
    fprintf(fOut, " leap-generated prep residue\n");
    fprintf(fOut, "%s.res\n", cPResName);
    fprintf(fOut, "%-3s  INT     0\n", cPResName);
    fprintf(fOut, "CHANGE   NOMIT DU   BEG \n");
    fprintf(fOut, "   0.00000\n");
    fprintf(fOut, "   %s\n   %s\n   %s\n",
            "1   DUMM  DU    M      0.000000  0.000000  0.000000  0.0",
            "2   DUMM  DU    M      1.000000  0.000000  0.000000  0.0",
            "3   DUMM  DU    M      1.000000  1.000000  0.000000  0.0");

    /*
     *  traverse tree from top, following main chain
     *      and ordering so that branches precede further
     *      main chain, printing atoms
     *
     *      reset atom temp ints to use in marking 'seen' atoms
     *      and SeenId for tracking order in file for ordering
     *      atoms in impropers; dAtomTemps were set in
     *      ziMarkMainChainAtoms()
     */
    lTemp = lLoop((OBJEKT) rRes, ATOMS);
    while ((aAtom = (ATOM) oNext(&lTemp)) != NULL) {
        AtomSetTempInt(aAtom, -1);
        AtomSetSeenId(aAtom, -1);
    }

    vaLoopAtoms = vaVarArrayCreate(2 * ATOMTYPELEN + 1);
    i = 4;
    aAtom = aAtom0;
    aParentAtom = NULL;
    while (1) {
        ATOM aNextMain = NULL;
        VECTOR vPos;

        vPos = vAtomPosition(aAtom);
        fprintf(fOut, prepfmt, i,
                sAtomName(aAtom),
                sAtomType(aAtom),
                (char) dAtomTemp(aAtom),
                vPos.dX, vPos.dY, vPos.dZ, dAtomCharge(aAtom));
        AtomSetSeenId(aAtom, i++);
        /*
         *  find next main chain down in this residue,
         *      marking side chains
         */
        for (j = 0; j < iAtomCoordination(aAtom); j++) {
            aChildAtom = aAtomBondedNeighbor(aAtom, j);
            if (aChildAtom == aParentAtom)
                continue;
            /*
             *  stay within residue
             */
            if (cContainerWithin(aChildAtom) != cContainerWithin(aAtom))
                continue;
            if (dAtomTemp(aChildAtom) == (double) 'M') {
                aNextMain = aChildAtom;
            } else if (dAtomTemp(aChildAtom) == (double) 'x') {
                /*
                 *  side chain not seen before
                 */
                zPrintSideChain(fOut, aAtom, aChildAtom, &i, vaLoopAtoms);
            }
        }
        if (aAtom == aAtom1)
            break;
        aParentAtom = aAtom;
        aAtom = aNextMain;
    }

    /*
     *  atoms written; write extra line to indicate no
     *      internal coords follow
     */
    fprintf(fOut, "\n");

    /*
     *  write any LOOP closing bonds
     */
    iMax = iVarArrayElementCount(vaLoopAtoms);
    if (iMax) {
        fprintf(fOut, "\nLOOP\n");
        for (i = 0; i < iMax; i++) {
            fprintf(fOut, " %s\n", PVAI(vaLoopAtoms, char, i));
        }
    }

    /*
     *  write list of potential IMPROPERs
     */
    lTemp = lLoop((OBJEKT) rRes, IMPROPERS);
    for (i = 0; oNext(&lTemp);) {
        ATOM aAtomA, aAtomB, aAtomC, aAtomD;
        int iSame;

        LoopGetImproper(&lTemp, &aAtomA, &aAtomB, &aAtomC, &aAtomD);

        if (aAtomC == aAtom0 || aAtomC == aAtom1)
            continue;
        if (iAtomCoordination(aAtomC) != 3)
            continue;

        if (!i++)
            fprintf(fOut, "\nIMPROPER\n");

        /*
         *  sort peripheral atoms by type
         */
        iSame = strcmp(sAtomType(aAtomA), sAtomType(aAtomB));
        if (iSame > 0 ||
            (!iSame && iAtomSeenId(aAtomA) > iAtomSeenId(aAtomB))) {
            aChildAtom = aAtomB;
            aAtomB = aAtomA;
            aAtomA = aChildAtom;
        }
        iSame = strcmp(sAtomType(aAtomB), sAtomType(aAtomD));
        if (iSame > 0 ||
            (!iSame && iAtomSeenId(aAtomB) > iAtomSeenId(aAtomD))) {
            aChildAtom = aAtomD;
            aAtomD = aAtomB;
            aAtomB = aChildAtom;
        }
        iSame = strcmp(sAtomType(aAtomA), sAtomType(aAtomB));
        if (iSame > 0 ||
            (!iSame && iAtomSeenId(aAtomA) > iAtomSeenId(aAtomB))) {
            aChildAtom = aAtomB;
            aAtomB = aAtomA;
            aAtomA = aChildAtom;
        }
        fprintf(fOut, " %-4s %-4s %-4s %-4s\n",
                sAtomName(aAtomA), sAtomName(aAtomB),
                sAtomName(aAtomC), sAtomName(aAtomD));
    }
    if (iAtomCoordination(aAtom0) == 2) {
        ATOM aAtomA, aAtomB;
        if (!i++)
            fprintf(fOut, "\nIMPROPER\n");
        aAtomA = aAtomBondedNeighbor(aAtom0, 0);
        aAtomB = aAtomBondedNeighbor(aAtom0, 1);
        if (strcmp(sAtomType(aAtomA), sAtomType(aAtomB)) > 0) {
            aChildAtom = aAtomB;
            aAtomB = aAtomA;
            aAtomA = aChildAtom;
        }
        fprintf(fOut, " %-4s %-4s %-4s %-4s\n",
                "-M", sAtomName(aAtomA), sAtomName(aAtom0),
                sAtomName(aAtomB));
    }
    if (iAtomCoordination(aAtom1) == 2) {
        ATOM aAtomA, aAtomB;
        if (!i++)
            fprintf(fOut, "\nIMPROPER\n");
        aAtomA = aAtomBondedNeighbor(aAtom1, 0);
        aAtomB = aAtomBondedNeighbor(aAtom1, 1);
        if (strcmp(sAtomType(aAtomA), sAtomType(aAtomB)) > 0) {
            aChildAtom = aAtomB;
            aAtomB = aAtomA;
            aAtomA = aChildAtom;
        }
        fprintf(fOut, " %-4s %-4s %-4s %-4s\n",
                sAtomName(aAtomA), "+M", sAtomName(aAtom1),
                sAtomName(aAtomB));
    }

    /*
     *  done
     */
    fprintf(fOut, "\nDONE\n");


    /*
     *  clean up
     */
    VarArrayDestroy(&vaLoopAtoms);

}

void UnitIOSaveAmberPrep(UNIT uUnit, FILE * fOut)
{
    LOOP lResidues;
    RESIDUE rRes;

    /*
     *  write a default beginning-of-prep
     */
    fprintf(fOut, " 0 0 0\n\n");

    /*
     *  put each residue in
     */
    lResidues = lLoop((OBJEKT) uUnit, RESIDUES);
    while ((rRes = (RESIDUE) oNext(&lResidues)))
        WritePrepRes(rRes, fOut);
    /*
     *  terminate file
     */
    fprintf(fOut, "STOP\n");

}


