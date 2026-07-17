#ifndef UNITIO_H
#define UNITIO_H
/*
 *      File:   unitio.h
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
 *      Description:
 *              Part of the UNIT object.
 *              All file input/output routines have been
 *              placed in the file 'unitio.c'
 *
 */

/*      Modifications induced by the implementation of the savemol2 command
*        Christine Cezard (2007) 
*        Universite de Picardie - Jules Verne, Amiens
*         http://q4md-forcefieldtools.org
*         zbUnitIOIndexBondParameters and zUnitDoAtoms are now "extern functions" 
*/ 
 

#define AMBERINDEX(i)   3*(i-1)

extern bool     bUnitIOLoadTables(UNIT uUnit, DATABASE db);
extern void     UnitIOSaveTables(UNIT uUnit, DATABASE db);

extern void     UnitIOBuildTables(UNIT uUnit, PARMLIB plParameters,
                        bool *bPGeneratedParameters, bool bPert, bool bCheck);
extern void     UnitIOBuildFromTables(UNIT uUnit);
extern void     UnitIODestroyTables(UNIT uUnit);
extern bool     bUnitIOIndexBondParameters(PARMLIB plLib, UNIT uUnit, bool bPert);

extern void     UnitDoAtoms(UNIT uUnit, PARMLIB plParameters, RESIDUE rRes, int *iPPos, bool * bPFailed, bool bPert);
extern void     UnitIOSaveAmberParmFormat(UNIT uUnit, char *prmtopName, bool bPolar, bool bPert, bool bNetcdf);

extern void     UnitIOSaveAmberPrep( UNIT uUnit, FILE *fOut );

extern int UnitIOAmberOrderResidues(UNIT);
extern int UnitLabelMolecules(UNIT);

extern int UnitIOSaveAmberParmNetcdf(const char *fname, UNITt *uUnit, bool bPert, bool bPolar,
                                      VARARRAY vaExcludedAtoms, VARARRAY vaExcludedCount, VARARRAY vaNBIndexMatrix,
                                      VARARRAY vaNBParameters, VARARRAY vaNBIndex);
extern void SaveAmberParmCMAP(UNIT uUnit, FILE * fOut);
extern int SaveAmberParmCMAPNetcdf(UNIT uUnit, int ncid);
extern void UnitIOSaveAmberCoordNetcdf(UNIT uUnit, char *crdName);
extern void UnitIOSaveAmberCoord(UNIT uUnit, char *crdName);
extern int BondAugmentationFound(UNIT uUnit);
void UnitIOFindAndCountMolecules(UNIT uUnit);

extern int iMarkMainChainAtoms(RESIDUE rRes, bool bComplain);
extern void MarkSideChainAtoms(RESIDUE rRes);

/*
 *        Private data types
 *
 */

#define PERTURBED       0x00000001
#define BOUNDARY        0x00000002

typedef struct {
    CONTAINERNAMEt sName;
    CONTAINERNAMEt sPertName;
    ATOMTYPEt sType;
    ATOMTYPEt sPertType;
    int iTypeIndex;
    int iPertTypeIndex;
    int iElement;
    int iPertElement;
    double dCharge;
    double dPertCharge;
    int iResidueIndex;
    VECTOR vPos;
    VECTOR vVelocity;
    int iSequence;
    FLAGS fFlags;
    ATOM aAtom;
} SAVEATOMt;

typedef struct {
    CONTAINERNAMEt sName; // was STRING
    int iSequenceNumber;
    int iaConnectIndex[MAXCONNECT];
    int iNextChildSequence;
    int iAtomStartIndex;
    int iImagingAtomIndex;
    int iPdbResSeq;
    char sChainId[3], sICode[2];
    char sResidueType[2];
    RESIDUE rResidue;
} SAVERESIDUEt;


typedef struct {
    int iType;
    FLAGS fFlags;
    int iAtom1;
    int iAtom2;
    int iAtom3;
    int iAtom4;
    double dKx;
    double dX0;
    double dN;
    int iParmIndex;  /* This is filled in when */
                     /* the parameters are written */
                     /* to the file */
} SAVERESTRAINTt;

typedef struct {
    int iAtom1;
    int iAtom2;
    int iParmIndex;
    int iPertParmIndex;
    FLAGS fFlags;
} SAVEBONDt;

typedef struct {
    int iAtom1;
    int iAtom2;
    int iAtom3;
    int iParmIndex;
    int iPertParmIndex;
    FLAGS fFlags;
} SAVEANGLEt;

typedef struct {
    int iAtom1;
    int iAtom2;
    int iAtom3;
    int iAtom4;
    int iParmIndex;
    bool bProper;
    bool bCalc14;
    int iPertParmIndex;
    bool bPertCalc14;
    FLAGS fFlags;
} SAVETORSIONt;


typedef struct {
    int iAtom1;
    int iAtom2;
    FLAGS fFlags;
} SAVECONNECTIVITYt;


typedef struct {
    STRING sName;
    int iSequenceNumber;
    int iNextChildSequence;
    MOLECULE mMolecule;
} SAVEMOLECULEt;

typedef struct {
    char sAboveType[2];
    int iAboveIndex;
    char sBelowType[2];
    int iBelowIndex;
} SAVEHIERARCHYt;

typedef struct {
    int iConnect;
} SAVECONNECTt;


typedef struct {
    double dE;
    double dR;
    double dE14;  // CHARMM ext
    double dR14;  // CHARMM ext
    bool bCapableOfHBonding;
    typeStr sType;
} NONBONDt;

typedef struct {
    double dA;
    double dB;
    double dA14;
    double dB14;
} NONBONDACt;

                        /* Used to save the bounding box info */
                        /* All of this is stored in one OFF entry */
                        /* dUseBox is greater than zero if the bounding */
                        /* box is used, and not greater than zero if it is */
                        /* not */

typedef struct {
    double dUseBox;
    double dBeta;
    double dXWidth;
    double dYWidth;
    double dZWidth;
} SAVEBOXt;

typedef struct {
    double dUseCap;
    double dX;
    double dY;
    double dZ;
    double dRadius;
} SAVECAPt;


typedef struct {
    int iGroupIndex;
    int iIndexAtom;
} SAVEGROUPSt;

extern double dGBRadiusForAtom(SAVEATOMt *sa, int iElement, double dMass, bool bLastAtom);
extern double dGBScreenForElement(int iElement);

#endif  /* UNITIO_H */
