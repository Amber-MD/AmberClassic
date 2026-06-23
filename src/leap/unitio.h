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
 
#ifndef UNITIO_H
#define UNITIO_H


extern BOOL     zbUnitIOLoadTables(UNIT uUnit, DATABASE db);
extern void     zUnitIOSaveTables(UNIT uUnit, DATABASE db);

extern void     zUnitIOBuildTables(UNIT uUnit, PARMLIB plParameters,
                        BOOL *bPGeneratedParameters, BOOL bPert, BOOL bCheck);
extern void     zUnitIOBuildFromTables(UNIT uUnit);
extern void     zUnitIODestroyTables(UNIT uUnit);
extern BOOL     zbUnitIOIndexBondParameters(PARMLIB plLib, UNIT uUnit, BOOL bPert);

extern BOOL     zbUNitIOIndexC4Pairwise(UNIT uUnit, double daC4Pairwise); //New
extern void     zUnitDoAtoms(UNIT uUnit, PARMLIB plParameters, RESIDUE rRes, int *iPPos, BOOL * bPFailed, BOOL bPert);
extern void     zUnitIOSaveAmberParmFormat(UNIT uUnit, FILE *fOut,
                        char *crdName, BOOL bPolar, BOOL bPert, BOOL bNetcdf, char sA[8][16], char sB[8][16], double daC4Type[16], int iC4count ); //NewT
extern void     zUnitIOSaveAmberParmFormat_old(UNIT uUnit, FILE *fOut,
                        char *crdName, BOOL bPolar, BOOL bPert, char sA[8][16], char sB[8][16], double daC4Type[16], int iC4count ); //NewT

extern void     UnitIOSaveAmberPrep( UNIT uUnit, FILE *fOut );

extern int      zUnitIOAmberOrderResidues(UNIT);

//extern void     UnitIOSaveC4Type( UNIT uUnit, char *sA, char *sB, double daC4Type ); //NewT


extern int iMarkMainChainAtoms(RESIDUE rRes, int complain);
extern void MarkSideChains(RESIDUE rRes);

/*
 *        Private data types
 *
 */


#define PERTURBED       0x00000001
#define BOUNDARY        0x00000002
#define JSBFAC          0.89089872 /* = 1/(2^(1/6)), needed for Jayaram et al. (M)GB radii */

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
    int iType;
    FLAGS fFlags;
    int iAtom1;
    int iAtom2;
    int iAtom3;
    int iAtom4;
    double dKx;
    double dX0;
    double dN;
    int iParmIndex;                /* This is filled in when */
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

// This is the New function added to help implementing
// atom-specific pariwise C4 interactions.
typedef struct {
    int iAtom1;
    int iAtom2;
    int iParmIndex;
    double daC4Pairwise;
} SAVEC4Pairwiset;


typedef struct {
    int iAtom1;
    int iAtom2;
    int iAtom3;
    int iParmIndex;
    int iPertParmIndex;
    FLAGS fFlags;
} SAVEANGLEt;

typedef struct {
    BOOL bProper;
    int iAtom1;
    int iAtom2;
    int iAtom3;
    int iAtom4;
    int iParmIndex;
    BOOL bCalc14;
    int iPertParmIndex;
    BOOL bPertCalc14;
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
    BOOL bCapableOfHBonding;
    double dE;
    double dR;
    double dE14;
    double dR14;
    typeStr sType;
} NONBONDt;

typedef struct {
  typeStr         sType1;
  typeStr         sType2;
  double          dEI;
  double          dEJ;
  double          dRI;
  double          dRJ;
  double          dA;
  double          dC;
  DESCRIPTION     sDesc;
} NBEDITt;

typedef struct {
    double dA;
    double dC;
    double d4; //NewT
    double dA14;
    double dC14;
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



#endif  /* UNITIO_H */
