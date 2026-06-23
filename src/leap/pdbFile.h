/*
 *      File:   pdbFile.h
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
 *              Read a PDB file, resolving problems of duplicate atoms
 *              and missing atoms.
 */
#include <stdint.h>
#include "avl.h"
#include "pdb.h"
#include "neighbors.h"
#include "varArray.h"
#include "classes.h"


extern void     PdbAppendToResMap( LIST lEntries );
extern void     PdbAppendToAtomMap( LIST lEntries );
extern void     PdbClearResMap(void);
extern void     PdbClearAtomMap(void);
extern void     PdbDisplayResMap(void);
extern void     PdbDisplayAtomMap(void);
extern int      fGetPdbResMapped(char *sResName);

extern void     PdbWrite(FILE *fOut, UNIT uUnit);
extern UNIT     uPdbRead(FILE *fPdb, VARARRAY vaUnits, BOOL bFormatCIF);


// Private/internal:

                /* PDBWRITEt is used to store data for writing PDB files */

typedef struct  {
        FILE    *fPdbFile;
        int     iRecordNumber;
        int     iResidueSeq;
        STRING  sResidueName;
        IX_DESC ixResCount; // Count unique residues
} PDBWRITEt;

typedef struct  {
        char    sName[8];
        int     iTerminator;
        char    sChainId[3], iCode;
        int     iPdbSequence;
        int     iFirstAtom;
        RESIDUE rResidue; // pointer to residue in growing UNIT, only for auto-link connections
} RESIDUENAMEt;

// ATOM info: name, xyz, element, atomSerial (only for CONECT), rest go to RESIDUENAMEt
typedef struct  {
        float   x,y,z;
        char    sName[6];
        uint8_t iElement;
        int     iAtomSerial;
        int     iResNameIndex; // parent index
} ATOMNAMEt;

typedef struct {
        char            type; // S=SSBND, C=covalent, M=metal, H=H-bond
        pdb_aname       name[2];
        char            alt_loc[2];
        pdb_residue     residues[2];
        struct {
            uint8_t op;
            int8_t dx,dy,dz;
        }               symop[2];
        double          distance;
} LINKt;

// Data layout for ixAtomIndex IX_RECORD with fixed key as struct members
typedef struct {
    int     resSeq;
    char    chainID[2];
    char    iCode; //, altLoc;
    char    name[5];
} ATOMKEYt;

typedef struct  {
        BOOL    bUsed;
        MATRIX  mTransform; // 4x4 transform (rotate+translate+scale+skew)
} PDBMATRIXt;

// PDBREADt is used to store data for reading PDB files
typedef struct  {
        FILE            *fPdbFile;
        UNIT            uUnit;         // UNIT being built
        RESIDUE         rRes;          // Most recent RESIDUE added to uUnit

        VARARRAY        vaResidues;    // RESIDUENAMEt
        VARARRAY        vaAtomRecs;    // ATOMNAMEt
        IX_DESC         ixAtomIndex;   // AVL atom identifier to ATOMNAMEt index (and unique check)
        VARARRAY        vaAtoms;       // ATOM (pointer into generated result UNIT members)
        VARARRAY        vaMatrices;    // PDBMATRIXt
        int             iMaxSerialNum; // Highest atomSerial, used to allocate vaAtoms

        VARARRAY        vaUnits;       // Template UNITs to be added to uUnit
        int             iNextUnit;     // index into vaUnits

        VARARRAY        vaPoints;        // Points for fast grid bond search
        unsigned int    *iPGroupStart;   // starting index of each point group (MALLOC)
        NeighborGrid    *ngBondGrid;     // Fast grid search data structure
        IX_DESC         ixResMap;        // AVL tree index auto-match ResMap
        //IX_DESC         ixPatchFilenames;// AVL tree to cache patch filenames found (not implemented yet)

        VARARRAY        vaConectRecs;    // struct pdb_conect
        VARARRAY        vaLinkRecs;      // struct pdb_link (OR converted pdb_ssbond)

        FILE            *fpPatchFileOut; // generated auto-link patch file

        int  iSubstRes1, iSubstRes2;
        char sSubstAtom1[32], sSubstAtom2[32];

} PDBREADt;

#define CHAINID_LIST_LEN (26+10+26)
//#define CHAINID_LIST_LEN (26+10+26+29)  // if we include special chars
extern const char *GsChainIdList;

extern void     CifReadFile( PDBREADt *prPRead);
