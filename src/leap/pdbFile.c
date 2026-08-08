/*
 *      File:   pdbFile.c
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

#include        <limits.h>

#include        "basics.h"

#include        "classes.h"
#include        "dictionary.h"
#include        "database.h"
#include        "pdb.h"
#include        "varArray.h"
#include        "avl.h"
#include        "sort.h"
#include        "matrix.h"
#include        "model.h"
#include        "minimizer.h"
#include        "leap.h"
#include        "defaults.h"
#include        "build.h"
#include        "atom.h"
#include        "elements.h"
#include        "hybrid36.h"

                /* PDBWRITEt is used to store data for writing PDB files */

typedef struct  {
        FILE    *fPdbFile;
        int     iRecordNumber;
        int     iResidueSeq;
        STRING  sResidueName;
} PDBWRITEt;

typedef struct  {
        char    sName[10];
        int     iTerminator;
        char    sChainId[3], iCode;
        int     iPdbSequence;
} RESIDUENAMEt;

typedef struct  {
        BOOL    bUsed;
        MATRIX  mTransform;
} PDBMATRIXt;

/* PDBREADt is used to store data for reading PDB files
 * (vaUnits) contains a VARARRAY of UNITs to take UNITs
 *      from for each residue read from the PDB file.
 *      The UNITs should all contain a single RESIDUE
 *      but don't need to, the UNITs are searched for
 *      ATOMs with names that match the names read from
 *      the PDB file.
 * (bRelax)(dRms) are used to determine whether
 *      the geometry of ATOMs not specified within the
 *      PDB file that are bonded to more than one other ATOM
 *      that is specified should be relaxed
 * (vaResidues) stores the residue names and termination
 *      status of residues read from the PDB file
 * (vaResidueSeq) stores the residue sequence number
 *      of residues read from the PDB file
 * (vaMatrices) stores the matrices that relate symmetry
 *      related monomers
 * (uUnit) will contain the UNIT read from the PDB file
 * (rRes) is the current RESIDUE to load atom records into
 * (iPdbSequence) is the sequence number read from the PDB
 *      file of the current residue.
 * (iNextUnit) is the index of the next UNIT to get from
 *      (vaUnits)
 */


typedef struct  {
        FILE            *fPdbFile;
        VARARRAY        vaUnits;
        VARARRAY        vaResidues;
        VARARRAY        vaResidueSeq;
        VARARRAY        vaAtoms;
        VARARRAY        vaMatrices;
        UNIT            uUnit;
        RESIDUE         rRes;
        int             iPdbSequence;
        int             iNextUnit;
        int             iMaxSerialNum;
} PDBREADt;



/*
----------------------------------------------------------------------

        Static variables

*/

static  DICTIONARY      SdResidueNameMap = NULL;
static  DICTIONARY      SdAtomNameMap = NULL;


extern int zUnitIOAmberOrderResidues( UNIT );


/*
 *----------------------------------------------------------------------
 *
 *        Private routines
 */


/*
 *      zcPPdbTerminationCode
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return a string that represents the termination
 *      type.  This is only for display purposes.
 */
static char *
zcPPdbTerminationCode( int iCode )
{
    switch ( iCode ) {
        case NOEND:
            return("Nonterminal");
            break;
        case FIRSTEND:
            return("Terminal/beginning");
            break;
        case LASTEND:
            return("Terminal/last");
            break;
        default:
            DFATAL(( "Illegal termination type" ));
            break;
    }
    return(NULL);       /* for lint */
}


/*
 *      zPdbToNameMapKey
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Build a Name Map DICTIONARY key from a termination type
 *      and name string.
 */
static void
zPdbToNameMapKey( int iPart, char *sPart, char *sKey )
{
    switch ( iPart ) {
        case NOEND:
            strcpy( sKey, sPart );
            break;
        case FIRSTEND:
            sprintf( sKey, "0%s", sPart );
            break;
        case LASTEND:
            sprintf( sKey, "1%s", sPart );
            break;
        default:
            DFATAL(( "Illegal termination type" ));
            break;
    }
}


/*
 *      zPdbFromNameMapKey
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Build a Name Map DICTIONARY key from a termination type
 *      and name string.
 */
static void
zPdbFromNameMapKey( char *sKey, int *iPPart, char *sPart )
{
    switch ( sKey[0] ) {
        case '0':
            *iPPart = FIRSTEND;
            strcpy( sPart, &(sKey[1]) );
            break;
        case '1':
            *iPPart = LASTEND;
            strcpy( sPart, &(sKey[1]) );
            break;
        default:
            *iPPart = NOEND;
            strcpy( sPart, sKey );
            break;
    }
}



/*
 *      zPdbNameMapAdd
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Add the name mapping to the Name Map DICTIONARY.
 */
static void
zPdbNameMapAdd( DICTIONARY *PSdNameMap, char *sKey, char *sData )
{
int             iLen;
char            *cPNew, *cPOld;

                /* If there is no DICTIONARY then create one */

    if ( !*PSdNameMap ) 
        *PSdNameMap = dDictionaryCreate();

                /* Allocate memory for the data */

    iLen = strlen(sData);
    MALLOC( cPNew, char*, iLen+1 );
    strcpy( cPNew, sData );

    MESSAGE(( "Adding %s:%s to name list.\n", sKey, sData ));
                /* Delete the old entry if it exists */
    cPOld = (char*)yPDictionaryDelete( *PSdNameMap, sKey );
    if ( cPOld != NULL ) {
        VP0(( "Substituting map %s -> %s  for  %s -> %s\n",
                sKey, cPNew, sKey, cPOld ));
        FREE(cPOld);
    }

                /* Create the new entry */
    DictionaryAdd( *PSdNameMap, sKey, (GENP)cPNew );

}



/*
 *      zPdbNameMapDestroy
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Empty the Name Map.
 */
static void
zPdbNameMapDestroy( DICTIONARY *PSdNameMap )
{
DICTLOOP        dlEntries;
char            *cPData;

    VP0(( "Clearing name map.\n" ));
    if ( *PSdNameMap == NULL ) return;
    dlEntries = ydlDictionaryLoop(*PSdNameMap);
    while ( (cPData = (char*)yPDictionaryNext(*PSdNameMap,&dlEntries))
                                != NULL ) {
        FREE(cPData);
    }
    DictionaryDestroy(PSdNameMap);
    *PSdNameMap = NULL;
}



/*
 *      zbPdbBuildResidueNameMapEntry
 *      zbPdbBuildAtomNameMapEntry
 *
 *      Author: Christian Schafmeister (1991)
 *      Author: Bill Ross (1993) - Added Atom name map
 *
 *      The NameMap is a DICTIONARY that either maps residue names 
 *      from PDB files to variable names of UNITs within LEaP
 *      or maps atom names to other possible atom names.
 *      This is used when the PDB file itself is used to
 *      construct the sequence of UNITs to concatenate together.
 *
 *      Parse the LIST and convert it into a key/data
 *      pair for entry into a DICTIONARY.
 *
 *      The LIST can have two or three entries.
 *
 *      Two entries:    { OSTRING_FROM OSTRING_TO }
 *
 *              OSTRING_FROM is the key.
 *              OSTRING_TO is the map.
 *
 *      For the residue map, this form designates a main chain residue type.
 *      For the atom map, this is the only allowed form.
 *
 *      Three entries:  { ODOUBLE_TERMINATE OSTRING_FROM OSTRING_TO }
 *
 *      For the residue map,
 *
 *              ODOUBLE_TERMINATE can be 0 or 1 for connect0 terminator
 *                      residues or connect1 terminator residues.
 *                      (First residue or last residue respectively).
 *              OSTRING_FROM is the residue name read from the PDB file.
 *              OSTRING_TO is the name of the variable within LEaP.
 *
 *      Print errors if the entries are not the right type.
 */
static BOOL
zbPdbBuildResidueNameMapEntry( LIST lEntry, int iMap, char *sKey, char *sData )
{
STRING          sTempKey;
ASSOC           aA, aB, aC, aD;
OBJEKT          oA, oB, oC;
LISTLOOP        llEntry;
int             iInt;

    if ( iObjectType(lEntry) != LISTid ) {
        VPWARN(( "Map entry %d is not a list. Ignored.\n", iMap ));
        return(FALSE);
    }
    aA = NULL;
    aB = NULL;
    aC = NULL;
    aD = NULL;
    llEntry = llListLoop(lEntry);
    aA = (ASSOC)oListNext(&llEntry);
    if ( aA != NULL ) {
        aB = (ASSOC)oListNext(&llEntry);
        if ( aB != NULL ) {
            aC = (ASSOC)oListNext(&llEntry);
            if ( aC != NULL ) {
                aD = (ASSOC)oListNext(&llEntry);
            }
        }
    }

    if ( aA == NULL || aB == NULL ) goto BADTYPE;
    oA = (OBJEKT)oAssocObject(aA);
    if ( iObjectType(oA) == ODOUBLEid ) {

        if ( aC == NULL ) goto BADTYPE;
        if ( aD != NULL ) goto BADTYPE;
        oB = (OBJEKT)oAssocObject(aB);
        oC = (OBJEKT)oAssocObject(aC);
        if ( iObjectType(oB) != OSTRINGid || 
                iObjectType(oC) != OSTRINGid ) goto BADTYPE;
        iInt = (int)dODouble(oA);
        strcpy( sTempKey, sOString(oB) );
        strcpy( sData, sOString(oC) );
        if ( (iInt != 0 && iInt != 1) ||
                strlen(sTempKey) == 0 ||
                strlen(sData) == 0 ) goto BADTYPE;
        if ( iInt == 0 ) iInt = FIRSTEND;
        else iInt = LASTEND;
        zPdbToNameMapKey( iInt, sTempKey, sKey );
    } else {
        if ( aC != NULL ) goto BADTYPE;
        oB = (OBJEKT)oAssocObject(aB);
        if ( iObjectType(oA) != OSTRINGid ||
                iObjectType(oB) != OSTRINGid ) goto BADTYPE;
        strcpy( sTempKey, sOString(oA) );
        strcpy( sData, sOString(oB) );
        if ( strlen(sTempKey) == 0 ||
                strlen(sData) == 0 ) goto BADTYPE;
        zPdbToNameMapKey( NOEND, sTempKey, sKey );
    }

    return(TRUE);

BADTYPE:
    VPWARN(( "Residue Map entry %d must have the form %s. Ignored\n",
                iMap,
                "{ [0 or 1] string string }" ));
    return(FALSE);
}

static BOOL
zbPdbBuildAtomNameMapEntry( LIST lEntry, int iMap, char *sKey, char *sData )
{
STRING          sTempKey;
ASSOC           aA, aB, aC;
OBJEKT          oA, oB;
LISTLOOP        llEntry;

    if ( iObjectType(lEntry) != LISTid ) {
        VPWARN(( "Map entry %d is not a list. Ignored.\n", iMap ));
        return(FALSE);
    }
    aA = NULL;
    aB = NULL;
    aC = NULL;
    llEntry = llListLoop(lEntry);
    aA = (ASSOC)oListNext(&llEntry);
    if ( aA != NULL ) {
        aB = (ASSOC)oListNext(&llEntry);
        if ( aB != NULL ) {
            aC = (ASSOC)oListNext(&llEntry);
            if ( aC != NULL ) {
                goto BADTYPE;
            }
        }
    }

    if ( aA == NULL || aB == NULL ) goto BADTYPE;
    oA = (OBJEKT)oAssocObject(aA);
    oB = (OBJEKT)oAssocObject(aB);
    if ( iObjectType(oA) != OSTRINGid || iObjectType(oB) != OSTRINGid ) 
        goto BADTYPE;
    strcpy( sTempKey, sOString(oA) );
    strcpy( sData, sOString(oB) );
    if ( strlen(sTempKey) == 0 || strlen(sData) == 0 ) 
        goto BADTYPE;

    zPdbToNameMapKey( NOEND, sTempKey, sKey );

    return(TRUE);

BADTYPE:
    VPWARN(( "Atom Map entry %d must have the form %s. Ignored\n",
                iMap, "{ string string }" ));
    return(FALSE);
}





/*
 *      zcPPdbMapName
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Lookup the variable name mapped to the residue or atom
 *      name, termination type.
 *      Return NULL if none was found.
 */
static char *
zcPPdbMapName( DICTIONARY SdNameMap, int iType, char *sName )
{
char            *cPData;
STRING          sKey;

    if ( SdNameMap == NULL ) 
        return(NULL);
    zPdbToNameMapKey( iType, sName, sKey );
    cPData = (char*)yPDictionaryFind( SdNameMap, sKey );
    return(cPData);
}


    


/*
 *      zPdbConvoluteAtomName
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Construct an ATOM name the way AMBER does it.
 *      PDB ATOM names are four characters long  
 *
 *              AABB
 *
 *              AA is the Atomic symbol RIGHT justified.
 *              BB is the identifier LEFT justified.
 *
 *
 *              'C1'   -> ' C1 '
 *              'C12'  -> ' C12'
 *              'CL2'  -> 'CL2 '
 *              'H12'  -> ' H12'
 *              'H123' -> 'H123'
 *              'HG11' -> 'HG11'  ** NOTE ** New atoms names can deviate from alignment
 */
static void
zPdbConvoluteAtomName( ATOM aAtom, STRING sName, STRING sResname, STRING sConvoluted )
{
STRING  sTemp, sElement, sDesc, sTemp2;
int     iElement;

    strcpy( sTemp, sName );

                /* If no name is defined then give it a default one */

    if ( strlen(sTemp) == 0 ) {
        strcpy( sTemp, "??" );
    }

                /* Shorten the name to four characters */

    if ( strlen(sTemp) > 4 ) {
        sTemp[4] = '\0';
    }

                /* Figure out which element the atom is, first try the ATOM */
                /* otherwise try to get it from the name */

    iElement = iAtomElement(aAtom);
    if ( bElementLegalNumber(iElement) == FALSE ) {
        iElement = iElementNumberFromAmber( sTemp );
    }

                /* Break the name into Element part and descriptor */

    if ( bElementLegalNumber(iElement) == FALSE ) {
        sElement[0] = sTemp[0];
        sElement[1] = '\0';
        if ( strlen(sTemp) <= 1 ) strcpy( sDesc, "" );
        else strcpy( sDesc, &(sTemp[1]) );
    } else {
        if ( strlen(sElementName(iElement,sTemp2)) == 2 ) {
            sElement[0] = sTemp[0];
            sElement[1] = sTemp[1];
            sElement[2] = '\0';
            if ( strlen(sTemp) <= 2 ) strcpy( sDesc, "" );
            else strcpy( sDesc, &(sTemp[2]) );
        } else {
            sElement[0] = sTemp[0];
            sElement[1] = '\0';
            if ( strlen(sTemp) <= 1 ) strcpy( sDesc, "" );
            else strcpy( sDesc, &(sTemp[1]) );
        }
    }

                /* Now combine the sElement and sDesc into sTemp */

        if ( strlen(sElement) == 2 ) {
                strcpy( sTemp, sElement );
        } else {
                strcpy( sTemp, " " );
                strcat( sTemp, sElement );
        }
        strcat( sTemp, sDesc );
        if ( strlen(sTemp) > 4 ) {
            sTemp[0] = sTemp[1];
            sTemp[1] = sTemp[2];
            sTemp[2] = sTemp[3];
            sTemp[3] = sTemp[4];
            sTemp[4] = '\0';
        }

                /* Make sure it's exactly four characters long */

    strcat( sTemp, "      " );
    sTemp[4] = '\0';
    strcpy( sConvoluted, sTemp );
}



/*
 *  zPdbDeconvoluteAtomName
 *
 *  Author: Christian Schafmeister (1991)
 *
 *  Deconvolute the ATOM name.
 *  See zPdbConvoluteAtomName.
 */
static void
zPdbDeconvoluteAtomName( STRING sConvoluted, STRING sName, int *iPElement )
{
STRING  sTemp, sElement;
int     iElement;

    iElement = iElementNumberFromAmber(sConvoluted);

    /* If element starts at char 2, element should be 1 character.
     * but in newer residues, single char element can start at column 1
     * e.g. there are cases where HG?? = H and not Hg, so element column is needed
     */
    /*
     *    if element is not legal, or its name begins with something other
     *     than a blank or a digit, just strip spaces and go on
     */
    if ( !bElementLegalNumber(iElement) ||
            (sConvoluted[0] != ' ' && !isdigit(sConvoluted[0]) )) {
        sRemoveSpaces( sConvoluted, sName );
    } else {
        /* Deconvolute PDB names into AMBER style */

        sElement[0] = sConvoluted[1];
        sElement[1] = '\0';
        iElement = iElementNumberFromAmber(sElement);
        if ( ( sConvoluted[3] == ' ' || sConvoluted[3] == '\0' ) 
                && sConvoluted[0] != ' ' ) {
            strcpy( sTemp, &(sConvoluted[1]) );
            sTemp[2] = sConvoluted[0];
            sTemp[3] = '\0';
        } else {
            strcpy( sTemp, &(sConvoluted[1]) );
            sTemp[3] = sConvoluted[0];
            sTemp[4] = '\0';
        }
        sRemoveSpaces( sTemp, sName );
    }
    *iPElement = iElement;
}





/*
 *      zPdbReadScan
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Scan through the PDB file and build a highlevel description
 *      of the contents that will be used later to read the PDB file.
 *      The things that are read from the PDB file are:
 *              Residue names are read into vaResidues, and
 *                      for each residue a value is kept that 
 *                      says whether the residue is a starting residue,
 *                      a stopping residue, or a mainchain residue.
 *                      This is all stored in the vaResidues VARARRAY.
 *              Transformation matrices are stored in the vaMatrices VARARRAY.
 */

static void
zPdbReadScan( PDBREADt *prPRead )
{
pdb_record      p;
int             iPdbSequence, iResNum;
BOOL            bLastReadPdbRecordWasTer = FALSE;
BOOL            bNewChain, bNewRes;
RESIDUENAMEt    rnName;
MATRIX          *mPMatrix;
int             iSerial, iStart, i, iRow;
int             iTerm, iLast, iSerialNum;
int             iResidueProblems = 0;
int             iSaidit;
IX_REC          atomentry;
IX_DESC         atomindex;
int             *iPRes = NULL;
char            *c2PChain = NULL;
char            c2Chain[2];
char            *cPAtom = NULL;
char            *cPInsertionCode;
int             iResIdMax, iAtomSerialMax;
int             iPdbResSeqOffset=0, iAtomSerialOffset=0;
int             iCurrentBioMT=0;

    VPTRACEENTER(( "zPdbReadScan" ));

    /*
     *  create index for checking uniqueness of atom names
     *          key = int, char[4], char[2], char
     */
    if ( sizeof(int) + 7 > sizeof(atomentry.key) )
        DFATAL(( "default key len too long\n" ));
    create_index( &atomindex, 2, sizeof(int) + 7 );
    iPRes = (int *) atomentry.key;
    c2PChain = (char *) iPRes + sizeof(int);
    cPInsertionCode = c2PChain + 2;
    cPAtom = cPInsertionCode + 1;

    iPdbSequence = iResNum = -9999;
    *cPInsertionCode = '-';  /* valid codes are alphabetic or blank */
    memcpy(c2PChain,"  ",2);
    bNewChain = TRUE;
    iSerialNum = 0; // Tracks *maximum* atomSerial
    iAtomSerialMax = HY36_WIDTH_5_MAX;
    iResIdMax = HY36_WIDTH_4_MAX;

    do {
        char    *resname;

        p = pdb_read_record(prPRead->fPdbFile);

                /* Process the records */
                
        switch ( p.record_type ) {
            case PDB_ATOM:
            case PDB_HETATM:
                // Old Leap format detected, adjust overflow wrapping limit
                if (p.pdb.atom.leap_expanded) {
                    iAtomSerialMax = 99999;
                    iResIdMax = 9999;
                }
                if (p.pdb.atom.serial_num == 0 && iSerialNum == iAtomSerialMax)
                    iAtomSerialOffset += iAtomSerialMax+1; // atomSerial overflow unwrapping
                p.pdb.atom.serial_num += iAtomSerialOffset;
                if ( p.pdb.atom.serial_num > iSerialNum )
                    iSerialNum = p.pdb.atom.serial_num;

                // Convert string to fixed width 2-char, right justified
                if (p.pdb.atom.residue.chain_id[0]==0)
                    memcpy(c2Chain,"  ",2);
                else if (p.pdb.atom.residue.chain_id[1]==0) {
                    c2Chain[0]=' '; c2Chain[1]=p.pdb.atom.residue.chain_id[0];
                }
                else memcpy(c2Chain,p.pdb.atom.residue.chain_id,2);

                VPTRACE(( "Read PDB_ATOM or PDB_HETATM record.\n" ));
                VPTRACE(( "bNewChain = %d; Chain = %.2s; iPdbSequence = %d; iResNum = %d; "
                          "iSerialNum = %d\n", bNewChain, c2Chain, iPdbSequence,
                          iResNum, iSerialNum ));
                /*
                 *  allow for right-shifted resnames XXX right justified is standard
                 */
                resname = p.pdb.atom.residue.name;
                while (*resname == ' ')
                        resname++;

                iTerm = NOEND;
                bNewChain = bNewChain || memcmp(c2PChain,c2Chain,2);
                if ( bNewChain ) {
/*
                    if (p.pdb.atom.residue.chain_id[0]!= ' ')
                        VP0(( " (starting new molecule for chain \"%.2s\")\n", p.pdb.atom.residue.chain_id ));
                    else
                        VP0(( " (starting new molecule for chain %c)\n", p.pdb.atom.residue.chain_id[1] ));
*/
                    iTerm = FIRSTEND;
                    if ( (iLast = iVarArrayElementCount( prPRead->vaResidues )) ) 
                        iLast--;        /* last element, maybe 0th */
                    PVAI( (prPRead->vaResidues), 
                                RESIDUENAMEt,iLast)->iTerminator = LASTEND;
                    memcpy(c2PChain,c2Chain,2);
                    iPdbResSeqOffset = 0;
                }
                bNewRes = (bNewChain || p.pdb.atom.residue.seq_num != iPdbSequence ||
                            p.pdb.atom.residue.insert_code != *cPInsertionCode ||
                            bLastReadPdbRecordWasTer );
                if (!bNewRes && strcmp(rnName.sName, resname) ) {
                    /* First detection of a residue with the same sequence
                     * number and insertion code but a different name.
                     */
                    VPWARN(( "Name change in pdb file residue %.2s %d%c;\n"
                        "this residue is split into %s and %s.\n",
                        c2PChain, iPdbSequence, *cPInsertionCode, rnName.sName, resname));
                    bNewRes = TRUE;
                    iResidueProblems++;
                }
                if (bNewRes) {
                    VPTRACE(( "Detected a new residue.\n" ));
                    if (p.pdb.atom.residue.seq_num == 0 && iPdbSequence == iResIdMax)
                         iPdbResSeqOffset += iResIdMax + 1;
                    if (GDefaults.bPdbKeepChainId) {
                        iResNum = p.pdb.atom.residue.seq_num + iPdbResSeqOffset; // ResSeq overflow unwrapping
                    } else if (iResNum < 0 )
                        iResNum = p.pdb.atom.residue.seq_num;
                    else
                        iResNum++;
                    rnName.iTerminator = iTerm;
                    rnName.iPdbSequence = iResNum;
                    if (p.pdb.atom.residue.chain_id[0]==' ') {
                         if (p.pdb.atom.residue.chain_id[1]==' ')
                              rnName.sChainId[0]=0;
                         else
                             strcpy( rnName.sChainId, p.pdb.atom.residue.chain_id+1);
                    }
                    else strcpy( rnName.sChainId, p.pdb.atom.residue.chain_id);
                    strcpy( rnName.sName, resname );

                    VarArrayAdd( prPRead->vaResidues, (GENP)&rnName );
                    bLastReadPdbRecordWasTer = FALSE;

                    MESSAGE(( "Reading residue: <%s>\n", rnName.sName ));
                    iPdbSequence = p.pdb.atom.residue.seq_num;
                    *cPInsertionCode = p.pdb.atom.residue.insert_code;

                    *iPRes = iResNum;
                }
                bNewChain = FALSE;

                /*
                 *  add to avl tree for uniqueness testing
                 *      (using arbitrary unique recptr to
                 *      overcome avl restriction on both
                 *      dup key and recptr)
                 */
                memset( cPAtom, 0, 5);
                strcpy( cPAtom, p.pdb.atom.name );
                if ( add_key( &atomentry,  &atomindex ) != IX_OK )
                        DFATAL(("Adding atom to index\n" ));
                break;

            case PDB_TER:
                VPTRACE(( "Read PDB_TER record.\n" ));
                bLastReadPdbRecordWasTer = TRUE;
                        /* If you read a TER card then make the */
                        /* last RESIDUE read a terminating RESIDUE */
                if ( (iLast = iVarArrayElementCount( prPRead->vaResidues )) ) {
                    iLast--;    /* last element, maybe 0th */
                    PVAI( (prPRead->vaResidues), 
                                RESIDUENAMEt,iLast)->iTerminator = LASTEND;
                    bNewChain = TRUE;
                }
                break;

            case PDB_REMARK:
                // Currently only REMARK 350 is parsed for BIOMT
                if (!GDefaults.iPdbReadBioMT || p.pdb.remark.num != 350) break;
                if (!strncmp(p.pdb.remark.text,"BIOMOLECULE:",12)) {
                    iCurrentBioMT = atoi(p.pdb.remark.text+12);
                    VP0(( "Found Biomolecule definition # %d.\n",iCurrentBioMT));
                    if (iCurrentBioMT != 1) {
                         VPWARN(("ReadBioMT only works correctly with a single biomolecule definition\n"));
                    }
                    break;
                }
                if (strncmp(p.pdb.remark.text,"  BIOMT",6)) break;
                if (iCurrentBioMT != GDefaults.iPdbReadBioMT) break;
/*
rrrrrr iii sss...  columns 12-70 = 59 chars                 -------->|
REMARK nnn text...
REMARK 350 BIOMOLECULE: 1
REMARK 350 APPLY THE FOLLOWING TO CHAINS: A, B, C, D, E, F, G, H, I, J,
REMARK 350                    AND CHAINS: K, L, M, N
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000
REMARK 350   BIOMT1   2  0.309017 -0.809017  0.500000      542.20813
REMARK 350   BIOMT2   2  0.809017  0.500000  0.309017     -335.10305
REMARK 350   BIOMT3   2 -0.500000  0.309017  0.809017      207.10519
REMARK 350   BIOMT1   3 -0.809017 -0.500000  0.309017     1084.41630
...
*/
                double m1, m2, m3, v;
                int n = sscanf(p.pdb.remark.text,
                       "  BIOMT%1d %3d %lf %lf %lf %lf",
                       &iRow, &iSerial, &m1, &m2, &m3, &v);
                if (iRow==1) VP2(( "Read PDB_REMARK 350 BIOMT1 record %d\n", iSerial ));
                if(n<6 || iRow<1 || iRow>3 || iSerial<=0) break; //error
                p.pdb.mtrix.serial_num = iSerial;
                p.pdb.mtrix.row_num = iRow;
                p.pdb.mtrix.m1 = m1;
                p.pdb.mtrix.m2 = m2;;
                p.pdb.mtrix.m3 = m3;
                p.pdb.mtrix.v = v;
                // First unit is identity, does not need replication
                p.pdb.mtrix.given = (iSerial == 1);
                // Now drop through to PDB_MTRIX processing of our BIOMT copied matix

            case PDB_MTRIX:
                if (p.record_type == PDB_MTRIX) {
                    if (GDefaults.iPdbReadBioMT) break;
                    VPTRACE(( "Read PDB_MTRIX record.\n" ));
                }
                iStart = iVarArrayElementCount(prPRead->vaMatrices);
                iSerial = p.pdb.mtrix.serial_num;
                if ( iSerial>iStart ) {
                    VarArraySetSize( (prPRead->vaMatrices), iSerial );
                    for ( i=iStart; i<iSerial; i++ ) {
                        PVAI(prPRead->vaMatrices,PDBMATRIXt,i)->bUsed = FALSE;
                        MatrixDiagonal(
                            PVAI(prPRead->vaMatrices,PDBMATRIXt,i)->mTransform,
                            0.0,
                            0.0,
                            0.0 );
                    }
                }
                PVAI(prPRead->vaMatrices,PDBMATRIXt,iSerial-1)->bUsed =
                    (p.pdb.mtrix.given != 1);
                iRow = p.pdb.mtrix.row_num-1;
                mPMatrix = (MATRIX*)
                    PVAI(prPRead->vaMatrices,PDBMATRIXt,iSerial-1)->mTransform;
                (*mPMatrix)[0][iRow] = p.pdb.mtrix.m1;
                (*mPMatrix)[1][iRow] = p.pdb.mtrix.m2;
                (*mPMatrix)[2][iRow] = p.pdb.mtrix.m3;
                (*mPMatrix)[3][iRow] = p.pdb.mtrix.v;
                break;
            default:
                break;
        }
    } while ( !feof(prPRead->fPdbFile) );

                /* Make the last residue a LASTEND */

    if ( (iLast = iVarArrayElementCount( prPRead->vaResidues )) ) {
        iLast--;        /* could be 0th element */
        PVAI( (prPRead->vaResidues), RESIDUENAMEt,iLast)->iTerminator = LASTEND;

        prPRead->iMaxSerialNum = iSerialNum+1;
    } else
        prPRead->iMaxSerialNum = 0;

        /* Rewind the PDB file */

    fseek( prPRead->fPdbFile, 0, 0 );

    if ( iResidueProblems ) {
        VPNOTE(( "%d residues had naming warnings.\n"
            "Thus, there are split residues;\n"
            "residue sequence numbers will not correspond to those in the pdb.\n",
            iResidueProblems ));
    }

    /*
     *  check uniqueness of atom names
     */
    iSaidit = 0;
    first_key( &atomindex );
    while ( next_key( &atomentry, &atomindex) == IX_OK ) {
        if ( atomentry.count > 1 ) {
            VP0(( "-- residue %d: duplicate [%s] atoms (total %d)\n",
                        *iPRes, cPAtom, atomentry.count ));
            iSaidit++;
        }
    }
    if ( iSaidit ) {
        VPWARN(( "Atom names in each residue should be unique.\n" ));
        VP0(( "     (Same-name atoms are handled by using the first\n" ));
        VP0(( "      occurrence and by ignoring the rest.\n" ));
        VP0(( "      Many instances of duplicate atom names usually come\n" ));
        VP0(( "      from alternate conformations in the PDB file.)\n\n" ));
    }

    destroy_index( &atomindex );
    VPTRACEEXIT (( "zPdbReadScan" ));
}





/*
 *      zPdbConvertNamesAndSequenceNumbers
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      If (vaUnits) is NULL then try to find UNITs within
 *      the parsers variable space that have names that can be
 *      mapped to by names from (vaResidues).
 *      If (vaUnits) is not NULL then 
 *      check that the (vaUnits) VARARRAY has the same number
 *      of entries as the (vaResidues) VARARRAY.  If it does then
 *      just return.  If it does not, then if there are too few entries
 *      try to find UNITs from the parser variable space the can be
 *      mapped to from residue names within (vaResidues).
 *
 *      Also renumber the sequence numbers so that
 *      Each TER card raises the sequence number by 10000
 *      if GDefaults.iPdbLoadSequential == 0
 */
void
zPdbConvertNamesAndSequenceNumbers( PDBREADt *prPPdb )
{
UNIT            uUnit;
int             i;
STRING          sName;
RESIDUENAMEt    *rnPName;

    VP1(( "Matching PDB residue names to LEaP variables.\n" ));

                /* If there is no vaUnits VARARRAY then create one */

    if ( !(prPPdb->vaUnits) ) {
        prPPdb->vaUnits = vaVarArrayCreate(sizeof(UNIT));
    }

                /* Now check to see if there are enough */
                /* UNITs in the vaUnits VARARRAY */

    if ( iVarArrayElementCount(prPPdb->vaUnits) == 
                iVarArrayElementCount(prPPdb->vaResidues) ) 
        return;

    if ( iVarArrayElementCount(prPPdb->vaUnits) <
                iVarArrayElementCount(prPPdb->vaResidues) ) {

        int iStart = iVarArrayElementCount(prPPdb->vaUnits);
        int iEnd = iVarArrayElementCount(prPPdb->vaResidues);

        VarArraySetSize( prPPdb->vaUnits, iEnd );

        rnPName = PVAI( prPPdb->vaResidues, RESIDUENAMEt, iStart );
        for ( i=iStart; i<iEnd; i++, rnPName++ ) {

            char *cPName = zcPPdbMapName( SdResidueNameMap, 
                        rnPName->iTerminator, rnPName->sName ); 

            if ( cPName ) {
                strcpy( sName, cPName );
                VP1(( "Mapped residue %s, term: %s, seq. number: %d to: %s.\n",
                        rnPName->sName,
                        zcPPdbTerminationCode( rnPName->iTerminator), i,sName));
            } else {
                strcpy( sName, rnPName->sName ); 
                MESSAGE(( "(Residue %d: %s, %s, was not found in name map.)\n",
                    i, sName, zcPPdbTerminationCode( rnPName->iTerminator)));
            }
            uUnit = (UNIT)oVariable(sName);
            if ( uUnit == NULL) {
                VPWARN(( "Unknown residue: %s   number: %d   type: %s\n", 
                        sName, i, zcPPdbTerminationCode(rnPName->iTerminator)));
                if ( rnPName->iTerminator != NOEND ) {
                    VP0(( "..relaxing end constraints to try for a dbase match\n"));
                    cPName = zcPPdbMapName( SdResidueNameMap, NOEND, 
                                                        rnPName->sName ); 
                    if ( cPName ) {
                        strcpy( sName, cPName );
                        uUnit = (UNIT)oVariable(sName);
                    }

                    if ( uUnit == NULL )
                        VPWARN(( "  -no luck\n" ));
                    else
                        VP0(( "  -matched to non-end type residue %s\n",
                                                        sName ));
                }
            } else if ( iObjectType(uUnit) != UNITid ) {
                VPWARN(( "Invalid unit: %s   sequence number: %d\n", sName, i ));
            }
            *PVAI(prPPdb->vaUnits,UNIT,i) = uUnit;
        }
    } else {
        VP0(( "There are %d units specified in the command,\n",
                iVarArrayElementCount(prPPdb->vaUnits) ));
        VP0(( "and only %d residues within the PDB file.\n",
                iVarArrayElementCount(prPPdb->vaResidues) ));
        VPWARN(( "The extra units will be ignored.\n" ));
    }

}



/*
 *      zuPdbGetNextUnit
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the next UNIT that will be connected into the UNIT
 *      being constructed from the PDB file.
 *      Make a copy of UNIT from the (vaUnits) VARARRAY unless 
 *      it is not a valid UNIT when we print a message and create
 *      a UNIT with a single RESIDUE.
 *      Return TRUE in *bPStart if this new RESIDUE starts
 *      a new chain.
 */
static UNIT
zuPdbGetNextUnit( PDBREADt *prPPdb, BOOL *bPStart, int iResNum )
{
int             iCurUnit, iTerm;
UNIT            uOrig, uUnit;
RESIDUE         rRes;

    iCurUnit = prPPdb->iNextUnit;
    prPPdb->iNextUnit++;
    uOrig = (UNIT)*PVAI( prPPdb->vaUnits, UNIT, iCurUnit );
    if ( iCurUnit == 0 ) {
        iTerm = LASTEND;
    } else {
        iTerm = 
            PVAI( prPPdb->vaResidues, RESIDUENAMEt, iCurUnit-1 )->iTerminator;
    }
    if ( iTerm == LASTEND )
        *bPStart = TRUE;
    else
        *bPStart = FALSE;
                        
    if ( iObjectType(uOrig) != UNITid ) {
        VP0(( "Creating new UNIT for residue: %s sequence: %d\n",
                PVAI(prPPdb->vaResidues,RESIDUENAMEt,iCurUnit)->sName, iResNum));
        uUnit = (UNIT)oCreate(UNITid);
        rRes = (RESIDUE)oCreate(RESIDUEid);
        ContainerSetName( (CONTAINER)rRes, 
                PVAI(prPPdb->vaResidues,RESIDUENAMEt,iCurUnit)->sName );
        ContainerAdd( (CONTAINER)uUnit, (OBJEKT)rRes );
    } else {
        MESSAGE(( "Getting UNIT: %s\n", sContainerName((CONTAINER)uOrig) ));
        uUnit = (UNIT)oCopy((OBJEKT)uOrig);
    }

    return(uUnit);
}






/*
 *      zPdbBuildCoordinatesForContainer
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Search through the CONTAINER for ATOMs that do not have
 *      ATOMPOSITIONKNOWN and build external coordinates for those.
 *      Also search through ATOMs that are bonded to ATOMs that have
 *      ATOMPOSITIONKNOWN and build those that do not have ATOMPOSITIONKNOWN
 *      This will guarantee that all ATOMs bonded to ATOMs within the 
 *      CONTAINER will be built.
 */
static void
zPdbBuildCoordinatesForContainer( CONTAINER cCont, int *iPAddH, 
        int *iPAddHeavy, int *iPAddUnk )
{
LOOP            lAtoms, lSpan;
ATOM            aAtom, aTemp, aSpan;
int             i;
#ifdef DEBUG
STRING          sSpan;
#endif
                /* Fix the internal coordinates */
    BuildFixInternals( (UNIT) cCont );

                /* Loop through all ATOMs looking for those that */
                /* do not have positions known and build externals */
                /* for them and neighbors that are bonded to them */
    lAtoms = lLoop( (OBJEKT)cCont, ATOMS );
    while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
        if ( !bAtomFlagsSet( aAtom, ATOMPOSITIONKNOWN ) ) {
            lSpan = lLoop( (OBJEKT)aAtom, SPANNINGTREE );
            LoopDefineInvisibleAtoms( &lSpan, ATOMPOSITIONKNOWN );

                        /* Look for a collision with an ATOM whose */
                        /* ATOMPOSITIONKNOWN flag is set */
            aSpan = NULL;
            while ( (aTemp = (ATOM)oNext(&lSpan)) ) {
                if ( aSpan == NULL ) aSpan = aTemp;
                if ( aLoopLastCollisionAtom(&lSpan) != NULL ) {
                    aSpan = aLoopLastCollisionAtom(&lSpan);
                    break;
                }
            }

            MESSAGE(( "Building externals from: %s\n",
                        sContainerFullDescriptor( (CONTAINER)aSpan, sSpan ) ));

                        /* Build external coordinates */
            lSpan = lLoop( (OBJEKT)aSpan, SPANNINGTREE );
            LoopDefineInvisibleAtoms( &lSpan, ATOMPOSITIONKNOWN );
            BuildExternalsUsingFlags( &lSpan,
                                        0, ATOMPOSITIONKNOWN, 
                                        ATOMPOSITIONKNOWN, 0,
                                        iPAddH, iPAddHeavy, iPAddUnk, TRUE );
        } else {
            for ( i=0; i<iAtomCoordination(aAtom); i++ ) {
                aSpan = aAtomBondedNeighbor(aAtom,i);
                if ( !bAtomFlagsSet( aSpan, ATOMPOSITIONKNOWN ) ) {
                    MESSAGE(( "Building externals from: %s\n",
                        sContainerFullDescriptor( (CONTAINER)aSpan, sSpan ) ));

                    lSpan = lLoop( (OBJEKT)aSpan, SPANNINGTREE );
                    LoopDefineInvisibleAtoms( &lSpan, ATOMPOSITIONKNOWN );
                    BuildExternalsUsingFlags( &lSpan,
                                                0, ATOMPOSITIONKNOWN,
                                                ATOMPOSITIONKNOWN, 0,
                                                iPAddH, iPAddHeavy, iPAddUnk, 
                                                TRUE );
                }
            }
        }
    }

                /* Destroy the INTERNALs in the CONTAINER */

    lAtoms = lLoop( (OBJEKT)cCont, ATOMS );
    BuildDestroyInternals( &lAtoms );

}








/*
 *      zPdbAddAtom
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Use the atom record to define the coordinate of an ATOM
 *      in the UNIT.
 */
static int
zPdbAddAtom( PDBREADt *prPPdb, pdb_record prRec, BOOL bNewRes, int iResNum, 
        int *iPAddH, int *iPAddHeavy, int *iPAddUnk )
{
UNIT            uNew;
int             iLowest, iElement, iCreate = 0;
LOOP            lResidues;
RESIDUE         rRes, rUse;
STRING          sName, sTemp;
ATOM            aAtom;
VECTOR          vPos;
BOOL            bStartChain;


                /* Check if the NEXT residue should be obtained */
                /* If it should be then make a copy of the next UNIT */
                /* in (vaUnits) and find the RESIDUE with the lowest */
                /* sequence number */

    if ( bNewRes == TRUE ) {
                        /* Build the geometry for the old one */

        if ( prPPdb->rRes ) {
            zPdbBuildCoordinatesForContainer( (CONTAINER)prPPdb->rRes,
                                iPAddH, iPAddHeavy, iPAddUnk );
        }
                        /* Get a copy of the next UNIT */
                        /* to read atom records into */

        uNew = zuPdbGetNextUnit( prPPdb, &bStartChain, iResNum );
        // JMK TODO:  if (iCollectionSize( cContainerContents((OBJEKT)uNew) ) == 1) { -- skip the loop
        iLowest = INT_MAX;
        rUse = NULL;
        lResidues = lLoop( (OBJEKT)uNew, RESIDUES );
        while ( (rRes = (RESIDUE)oNext(&lResidues)) ) {
            if ( iContainerSequence(rRes) < iLowest ) {
                iLowest = iContainerSequence(rRes);
                rUse = rRes;
            }
        }
        prPPdb->rRes = rUse;
        if ( !rUse ) {
            DFATAL(( "Could not find any RESIDUEs within the UNIT\n" ));
        }
                /* Build INTERNAL coordinates for the new part */

        BuildInternalsForContainer( &(uNew->cHeader), 0, ATOMPOSITIONKNOWN );

        ContainerSetNextChildsSequence( prPPdb->uUnit, iResNum );
        RESIDUENAMEt *rName = PVAI( prPPdb->vaResidues,RESIDUENAMEt,prPPdb->iNextUnit-1 );

                /* If the new UNIT is the first UNIT in a chain then */
                /* Don't sequence it to the last UNIT, just join them */

        if ( bStartChain ) {

            lResidues = lLoop( (OBJEKT)uNew, RESIDUES );
            rRes = (RESIDUE)oNext(&lResidues);
            UnitJoin( prPPdb->uUnit, uNew );
        } else {

                /* Build INTERNAL coordinates for the linkage between the */
                /* UNITs that is about to be formed */
                /* Sequence the new UNIT into the PDB UNIT */

            BuildInternalsBetweenUnitsUsingFlags( prPPdb->uUnit, uNew, 0, 0 );
            lResidues = lLoop( (OBJEKT)uNew, RESIDUES );
            rRes = (RESIDUE)oNext(&lResidues);
            UnitSequence( prPPdb->uUnit, uNew );
        }

            /* Define the PDB sequence number */

// TODO - what if rRes is NULL?
        rRes->iPdbResSeq = rName->iPdbSequence;
        strcpy( rRes->sChainId,rName->sChainId);
        rRes->cICode = rName->iCode;
    }

// TODO: use pdb element column to determine Element if not defined in Amber lib?
    zPdbDeconvoluteAtomName( prRec.pdb.atom.name, sName, &iElement );
    aAtom = (ATOM)cContainerFindName( (CONTAINER)prPPdb->rRes, ATOMid, sName );

    /*
     *  If atom name doesn't match, try aliases
     */
    if ( aAtom == NULL ) {
        char    *cPData;


        cPData = zcPPdbMapName( SdAtomNameMap, NOEND, sName );
        if (cPData != NULL )
                aAtom = (ATOM)cContainerFindName( (CONTAINER)prPPdb->rRes, ATOMid, cPData);
    }


                /* If the ATOM was not within the RESIDUE report */
                /* that you are going to create it */
    if ( aAtom == NULL ) {
        aAtom = (ATOM)oCreate(ATOMid);
        ContainerAdd( (CONTAINER)prPPdb->rRes, (OBJEKT)aAtom );
        ContainerSetName( (CONTAINER)aAtom, sName );
        AtomSetElement( aAtom, iElement );
        MESSAGE(( "Read atom: %s and adding it to: %s\n", sName,
                        sContainerName((CONTAINER)prPPdb->rRes) ));
        VP0(( "Created a new atom named: %s within residue: %s\n",
                sName, sContainerFullDescriptor((CONTAINER)prPPdb->rRes,sTemp) ));
        iCreate = 1;
    } else {
        MESSAGE(( "Read atom: %s and found it in the RESIDUE\n",
                        sName ));
    }
    if (prRec.pdb.atom.serial_num >= 0 ) { // Leap-built atoms get -1, ignore them
        *PVAI( prPPdb->vaAtoms, ATOM, prRec.pdb.atom.serial_num ) = aAtom;
    }

                        /* Define its position */
    VectorDef( &vPos, prRec.pdb.atom.x, prRec.pdb.atom.y, prRec.pdb.atom.z );
    AtomSetPosition( aAtom, vPos );
    return( iCreate );
}






/*
 *      zPdbReadAndCreateUnit
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Read the PDB file records and create a UNIT that has all ATOM
 *      coordinates defined.  For each new residue in the PDB file, get the
 *      next UNIT out of the (vaUnits) VARARRAY and pick the RESIDUE out of it
 *      with the lowest sequence number, this is the RESIDUE which is used to
 *      find ATOMs whose coordinates are read from the PDB file.
 *
 *      Added use of PDB TER records for new residue detection;
 *      Scott R. Brozell (2021).
 *
 */
static void
zPdbReadAndCreateUnit( PDBREADt *prPPdb )
{
pdb_record      p;
STRING          sResName;
int             i, iAdd, iAddH = 0, iAddHeavy = 0, iAddUnk = 0, iCreate = 0;
ATOM            aA, aB, *aPa;
int             iPdbSequence, iResNum, iAtoms = 0;
BOOL            bLastReadPdbRecordWasTer = FALSE;
BOOL            bNewRes;
char            cInsertionCode, c2ChainId[2];

    VPTRACEENTER(( "zPdbReadAndCreateUnit" ));
    iResNum = -9999; // sequential starting at 1, the rest are copied from PDB file
    iPdbSequence = -9999;
    cInsertionCode = '-';  /* valid codes are alphabetic or blank */
    c2ChainId[0]=0;

                /* Set up the PDBREADt structure to read a PDB file */

    prPPdb->rRes      = NULL;
    prPPdb->uUnit = (UNIT)oCreate(UNITid);
    ContainerSetName(prPPdb->uUnit, "default_name");  
    prPPdb->iNextUnit = 0;
    if ( prPPdb->iMaxSerialNum == 0 ) {
        VPFATALEXIT(( "\tNo atoms!\n" ));
        return;
    }
    // This array is an atom pointer arr indexed by atomSerial; only used for CONECT lookups
    VarArraySetSize( prPPdb->vaAtoms, prPPdb->iMaxSerialNum );
    aPa = PVAI( prPPdb->vaAtoms, ATOM, 0 );
    for ( i=0; i < prPPdb->iMaxSerialNum; i++, aPa++ ) {
        *aPa = NULL;
    }

    do {
        p = pdb_read_record(prPPdb->fPdbFile);

                /* Process the records */
                
        switch ( p.record_type ) {
        
            case PDB_ATOM:
            case PDB_HETATM:
                VPTRACE(( "Process PDB_ATOM or PDB_HETATM record.\n" ));
                if ( p.pdb.atom.residue.seq_num != iPdbSequence ||
                            strcmp( sResName, p.pdb.atom.residue.name ) ||
                            memcmp( c2ChainId, p.pdb.atom.residue.chain_id, 2) ||
                            p.pdb.atom.residue.insert_code != cInsertionCode ||
                            bLastReadPdbRecordWasTer ) {
                    VPTRACE(( "Detected a new residue.\n" ));
                    bNewRes = TRUE;
                    if ( iResNum < 0 ) // JMK: FIXME - always start at 1!!!
                        iResNum = p.pdb.atom.residue.seq_num;
                    else
                        iResNum++;
                    iPdbSequence = p.pdb.atom.residue.seq_num;
                    memcpy( c2ChainId, p.pdb.atom.residue.chain_id, 2);
                    strcpy( sResName, p.pdb.atom.residue.name );
                    cInsertionCode = p.pdb.atom.residue.insert_code;
                    bLastReadPdbRecordWasTer = FALSE;
                } else {
                    bNewRes = FALSE;
                }
                iCreate += zPdbAddAtom( prPPdb, p, bNewRes, iResNum, 
                                &iAddH, &iAddHeavy, &iAddUnk );
                iAtoms++;
                break;

            case PDB_TER:
                VPTRACE(( "Process PDB_TER record.\n" ));
                bLastReadPdbRecordWasTer = TRUE;
                break;

            case PDB_CONECT:
                VPTRACE(( "Process PDB_CONECT record for atom %i.\n",
                    p.pdb.conect.serial_num ));
                /*
                 * prPPdb->iMaxSerialNum is the size of prPPdb->vaAtoms and is
                 * 1 more than the last inputted atom's iSerialNum; thus, since
                 * serial numbers start at 1 and prPPdb->vaAtoms is zero based,
                 * it has wasted space at the start but none at the end.
                 * So guard against out of bounds high indexing.  srb 5-2022.
                 */
                if ( p.pdb.conect.serial_num >= prPPdb->iMaxSerialNum ) {
                    VPWARN(( "Ignoring CONECT record for atom serial number "
                        "(%i) that is greater than\nthe maximum serial"
                        " number inputted (%i) from the pdb file.\n",
                        p.pdb.conect.serial_num, prPPdb->iMaxSerialNum - 1 ));
                    break;
                }
                aA = *PVAI( prPPdb->vaAtoms, ATOM, p.pdb.conect.serial_num );
                if ( aA == NULL ) {
                    VPWARN(( "Illegal CONECT record in pdb file.\n" ));
                    break;
                }
                for ( i=0; i<4; i++ ) {
                    if ( p.pdb.conect.covalent[i] == 0 ) continue;
                    if ( p.pdb.conect.covalent[i] >= prPPdb->iMaxSerialNum ) {
                        VPWARN(( "In CONECT record for atom serial number (%i)\n"
                            "ignoring bonded atom serial number (%i) that is "
                            "greater than\nthe maximum serial"
                            " number inputted (%i) from the pdb file.\n",
                            p.pdb.conect.serial_num, p.pdb.conect.covalent[i],
                            prPPdb->iMaxSerialNum - 1 ));
                        break;
                    }
                    aB = *PVAI(prPPdb->vaAtoms,ATOM,p.pdb.conect.covalent[i]);
                    if ( aB == NULL ) {
                        VPWARN(( "Illegal CONECT record in pdb file.\n" ));
                        break;
                    }
                    if ( !bAtomBondedTo( aA, aB ) ) {
                        AtomBondTo( aA, aB );
                    }
                }
                break;

            case PDB_END:
                VPTRACE(( "Process PDB_END record.\n" ));
                break;
                
            default:
                VPTRACE(( "No Op PDB record.\n" ));
                break;
        }
        
        
    } while ( p.record_type != PDB_END );
    zPdbBuildCoordinatesForContainer( (CONTAINER)prPPdb->rRes,
                                &iAddH, &iAddHeavy, &iAddUnk );

                /* Build geometry for everything that has not */
                /* been build */
    zPdbBuildCoordinatesForContainer( (CONTAINER)prPPdb->uUnit,
                                &iAddH, &iAddHeavy, &iAddUnk );

    VP0(( "  total atoms in file: %d\n", iAtoms ));
    if ( (iAdd = iAddH + iAddHeavy + iAddUnk) ) {
        VP0(("  Leap added %d missing atom%s according to residue templates:\n",
                                iAdd, (iAdd > 1 ? "s" : "") ));
        if ( iAddHeavy )
                VP0(( "       %d Heavy\n", iAddHeavy ));
        if ( iAddH )
                VP0(( "       %d H / lone pairs\n", iAddH ));
        if ( iAddUnk )
                VP0(( "       %d unknown element\n", iAddUnk ));
    }
    if ( iCreate )
        VP0(( "  The file contained %d atoms not in residue templates\n", 
                                                        iCreate ));
    if ( iAdd  &&  iAdd == iCreate ) {
        VPWARN(( "Since the number of added atoms equals the number of missing "
                "atoms, it is likely\n" "that some atoms had incorrect names; "
                "you may want to use addPdbAtomMap to map\n"
                "these names, or change the names in the PDB file.\n\n" ));
    }
    VPTRACEEXIT (( "zPdbReadAndCreateUnit" ));
}       






/*
 *      zPdbCreateSymmetryRelatedPairs
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Use the transform matrices to build symmetry related
 *      monomers of the UNIT in (uUnit).
 *      Combine all of the monomers into a single UNIT.
 */
static void
zPdbCreateSymmetryRelatedMonomers( PDBREADt *prPPdb )
{
int             iTransforms, iLast, i;
VECTOR          vPos;
MATRIX          mTransform;
UNIT            uOrig, uCopy;
ATOM            aAtom;
LOOP            lAtoms;

    iTransforms = iVarArrayElementCount( prPPdb->vaMatrices );
    iLast = -1;
    for ( i=0; i<iTransforms; i++ ) {
        if ( PVAI(prPPdb->vaMatrices,PDBMATRIXt,i)->bUsed ) {
            iLast = i;
        }
    }

                /* If there are no transforms defined then */
                /* don't do any transforms */

    if ( iLast == -1 ) return;

                /* Get the original monomer */

    uOrig = (UNIT)oCopy((OBJEKT)prPPdb->uUnit);
    for ( i=0; i<iTransforms; i++ ) {
        if ( PVAI(prPPdb->vaMatrices,PDBMATRIXt,i)->bUsed ) {
            MatrixCopy( mTransform,
                        PVAI(prPPdb->vaMatrices,PDBMATRIXt,i)->mTransform );

                        /* If this is not the last transform then make */
                        /* copies of uOrig */

            if ( i != iLast ) {
                uCopy = (UNIT)oCopy((OBJEKT)uOrig);
            } else {
                uCopy = uOrig;
            }

                        /* Transform all of the ATOMs */

            lAtoms = lLoop( (OBJEKT)uCopy, ATOMS );
            while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
                MatrixTimesVector( vPos, mTransform, vAtomPosition(aAtom) );
                AtomSetPosition( aAtom, vPos );
            }

                /* Now combine the contents of the UNITs */

            VP1(( "Building symmetry related monomer %d.\n",i+1));
            UnitJoin( prPPdb->uUnit, uCopy );

        }
    }

}






/*
 *      zPdbFileBegin
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Initialize variables for writing a PDB file.
 */
static void
zPdbFileBegin( PDBWRITEt *pwPFile, FILE *fOut )
{
    pwPFile->fPdbFile = fOut;
    pwPFile->iRecordNumber= 1;
    pwPFile->iResidueSeq = 1;
}



/*
 *      zPdbFileWriteAtomRecord
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Write an ATOM record to the PDB file.
 *
 *      rRes argument added to copy ChainId (TODO: should also replace pwPFile->sResidueName)
 */
void
zPdbFileWriteAtomRecord( PDBWRITEt *pwPFile, ATOM aAtom, RESIDUE rRes)
{
pdb_record      p;
char            *sElement="",*sPName;
int             iElement;
STRING          sTemp;

    p.pdb.atom.serial_num = pwPFile->iRecordNumber++;
    strcpy( p.pdb.atom.residue.name, pwPFile->sResidueName );


    sTemp[0]=' ';
    sPName = &sTemp[1];
    strcpy( sPName, sContainerName((CONTAINER)aAtom) );
    iElement = iAtomElement(aAtom);
    if (!bElementLegalNumber(iElement)) {
        iElement = iElementNumberFromAmber(sPName);
        if (iElement == HELIUM || iElement == HOLMIUM || iElement == HAFNIUM ||  iElement == MERCURY) iElement = HYDROGEN;
    }
    if (iElement != NOELEMENT) sElement=GeaElements[iElement].sName;

    // Include leading space in atom name if it is a single letter element symbol *and* the name is not 4 characters
    if (strlen(sElement)<2 && strlen(sPName) <4) sPName--;
    strcpy( p.pdb.atom.name, sPName );
    MESSAGE(( "Element: |%s|   pdb_name=|%s|\n", sElement, sName ));

    if (rRes && GDefaults.bPdbKeepChainId) {
        strncpy(p.pdb.atom.residue.chain_id,rRes->sChainId,2);
        p.pdb.atom.residue.chain_id[2]=0;
        p.pdb.atom.residue.seq_num = rRes->iPdbResSeq;
        p.pdb.atom.residue.insert_code = rRes->cICode;
    } else {
        strcpy(p.pdb.atom.residue.chain_id,"  ");
        p.pdb.atom.residue.seq_num = pwPFile->iResidueSeq;
        p.pdb.atom.residue.insert_code=' ';
    }
    p.pdb.atom.alt_loc = ' ';
    p.pdb.atom.x = dVX(&vAtomPosition(aAtom));
    p.pdb.atom.y = dVY(&vAtomPosition(aAtom));
    p.pdb.atom.z = dVZ(&vAtomPosition(aAtom));
    p.pdb.atom.occupancy = 1.0;
    double chg = dAtomCharge(aAtom);
    if (GDefaults.pdbwritecharges)
        p.pdb.atom.temp_factor = chg;
    else
        p.pdb.atom.temp_factor = 0.0;
    strcpy(p.pdb.atom.element,sElement);
    if (floor(chg)==chg && chg != 0.0 && fabs(chg)<10.0) {
        p.pdb.atom.fcharge[0] = '0' + abs((int)chg);
        p.pdb.atom.fcharge[1] = chg > 0 ? '+' : '-';
        p.pdb.atom.fcharge[2] = '\0';
    } else
        strcpy(p.pdb.atom.fcharge,"  ");
    p.record_type = PDB_ATOM;
    pdb_write_record( pwPFile->fPdbFile, &p, NULL, 0 );

}



/*
 *      zPdbFileWriteTermination
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Write a termination record.
 */
static void
zPdbFileWriteTermination( PDBWRITEt *pwPFile )
{
pdb_record      p;

    p.pdb.ter.serial_num = pwPFile->iRecordNumber;

    strcpy( p.pdb.ter.residue.name, pwPFile->sResidueName );

    memcpy(p.pdb.ter.residue.chain_id,"  ",2);
    p.pdb.ter.residue.seq_num = pwPFile->iResidueSeq;
    p.pdb.ter.residue.insert_code=' ';
    p.record_type = PDB_TER;
    pdb_write_record( pwPFile->fPdbFile, &p, NULL, 0 );
}




/*
 *      zPdbFileEnd
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Write the END record.
 */
static void
zPdbFileEnd( PDBWRITEt *pwPFile )
{
pdb_record      p;

    p.record_type = PDB_END;
    pdb_write_record( pwPFile->fPdbFile, &p, NULL, 0 );
}




/*
 *      zPdbFileWriteContainer
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Write the ATOM records for the ATOMs within the
 *      CONTAINER.
 */
static void
zPdbFileWriteContainer( PDBWRITEt *pwPFile, CONTAINER cCont )
{
const char      C_TERMINAL_PREFIX = 'C';
const char      N_TERMINAL_PREFIX = 'N';
const int       RESIDUE_NAME_LENGTH = 3;
LOOP            lContents;
ATOM            aAtom;
char            *cPTemp;

    cPTemp = sContainerName(cCont);
    if ( strlen( cPTemp ) > 3 && isalpha( *(cPTemp + 3) ) &&
            isalpha( *(cPTemp + 2) ) && isalpha( *(cPTemp + 1) ) &&
            (*cPTemp == N_TERMINAL_PREFIX || *cPTemp == C_TERMINAL_PREFIX) ) {
        /* Probable N-terminal or C-terminal amino acid name */
        /* Do not copy the prefix. */
        strncpy( pwPFile->sResidueName, cPTemp + 1, RESIDUE_NAME_LENGTH );
        pwPFile->sResidueName[ RESIDUE_NAME_LENGTH ] = '\0';
        VPWARN(( " Converting %c-terminal residue name to PDB format: %s -> %s\n",
                    *cPTemp, sContainerName(cCont), pwPFile->sResidueName ));
    }
    else {
        strncpy( pwPFile->sResidueName, cPTemp, RESIDUE_NAME_LENGTH );
        /* The intentional side effect is to truncate long names. */
        pwPFile->sResidueName[ RESIDUE_NAME_LENGTH ] = '\0';
        if ( strlen(cPTemp) > RESIDUE_NAME_LENGTH ) {
            VPWARN(( " Truncating residue name for PDB format: %s -> %s\n",
                        sContainerName(cCont), pwPFile->sResidueName ));
        }
    }

    if ( iObjectType(cCont) == ATOMid ) {
        zPdbFileWriteAtomRecord( pwPFile, (ATOM)cCont, NULL );
        pwPFile->iResidueSeq++;
    } else if ( iObjectType(cCont) == RESIDUEid ) {
        RESIDUE rRes = (RESIDUE) cCont;
        if (GDefaults.bPdbKeepChainId) {
            lContents = lLoop( (OBJEKT)cCont, DIRECTCONTENTSBYSEQNUM );
            while ( (aAtom = (ATOM)oNext(&lContents)) ) {
                zPdbFileWriteAtomRecord( pwPFile, aAtom, rRes);
            }
        } else {
            /* Force left-justification for older format output */
            /* TODO right-justify regardless of ChainId option - JMK */
            if (strlen(pwPFile->sResidueName) < 3) {
                strcpy(&pwPFile->sResidueName[2]," ");
                if (pwPFile->sResidueName[1]==0) pwPFile->sResidueName[1]=' ';
            }
            pwPFile->iResidueSeq = rRes->iTemp;
            if ( pwPFile->iResidueSeq == 0 )
                    pwPFile->iResidueSeq = 1;
            lContents = lLoop( (OBJEKT)cCont, DIRECTCONTENTSBYSEQNUM );
            while ( (aAtom = (ATOM)oNext(&lContents)) ) {
                zPdbFileWriteAtomRecord( pwPFile, aAtom, NULL);
            }
        }
    } else if ( iObjectType(cCont) == MOLECULEid ) {
        lContents = lLoop( (OBJEKT)cCont, ATOMS );
        while ( (aAtom = (ATOM)oNext(&lContents)) ) {
            zPdbFileWriteAtomRecord( pwPFile, aAtom, NULL );
        }
        pwPFile->iResidueSeq++;
    }
}

/*
==========================================================================

        Public routines
*/



/*
 *      uPdbRead
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Read the PDB file using the UNITs in lUnits.
 *      For each successive residue within the PDB file, take the next
 *      UNIT out of the list of UNITs in (vaUnits).
 *      If (vaUnits) is NULL then an attempt is made to
 *      find UNITs within the parser variable space that
 *      are mapped to by residue names read from the PDB file.
 *
 */
UNIT
uPdbRead( FILE *fPdb, VARARRAY vaUnits )
{
PDBREADt        prPdb;
//LOOP            lResidues, lAtoms;
//RESIDUE         rRes;
//ATOM            aAtom;
int             i;

    VPTRACEENTER(( "uPdbRead" ));
    prPdb.fPdbFile = fPdb;
    prPdb.vaUnits = vaUnits;
    prPdb.vaResidues = vaVarArrayCreate( sizeof(RESIDUENAMEt) );
    prPdb.vaMatrices = vaVarArrayCreate( sizeof(PDBMATRIXt) );
    prPdb.vaAtoms    = vaVarArrayCreate( sizeof(ATOM) );

                /* Scan the PDB file looking for residue names, */
                /* residue termination status, and matrices */

    zPdbReadScan( &prPdb );

                /* Convert the residue names read in from the PDB file */
                /* to UNIT names unless the (vaUnits) array is defined */
                /* and complete */

    zPdbConvertNamesAndSequenceNumbers( &prPdb );

    if ( vaUnits != NULL ) {
        int     mismatch, total = 0;
        /*
         *  loadpdbusingseq: see if residues match
         */
        VP0(( "  matching pdb residues -> sequence template\n" ));
        VP0(( "\tres\tpdb\ttemplate\n" ));
        for (i=0; i<iVarArrayElementCount(vaUnits); i++) {
                if ( ! strstr( (*PVAI( vaUnits, UNIT, i ))->cHeader.sName,
                               PVAI(prPdb.vaResidues, RESIDUENAMEt, i)->sName) )
                        mismatch = 1;
                else
                        mismatch = 0;
                VP0(( "\t%d%s\t%s\t%s\n", i+1, (mismatch ? "*" : ""),
                        PVAI(prPdb.vaResidues, RESIDUENAMEt, i)->sName,
                        (*PVAI( vaUnits, UNIT, i ))->cHeader.sName ));
                total += mismatch;
        }
        if (total) {
                VP0(( "  * = possible mismatch; total %d\n", total ));
                VP0(( "      (i.e. pdb name not a substring of template)\n\n"));
        }
    }

                /* Read the PDB file and create the UNIT */

    zPdbReadAndCreateUnit( &prPdb );

                /* Build the symmetry related monomers */

    zPdbCreateSymmetryRelatedMonomers( &prPdb );

                /* Clean up */

    VarArrayDestroy( &(prPdb.vaResidues) );
    VarArrayDestroy( &(prPdb.vaMatrices) );
    VarArrayDestroy( &(prPdb.vaAtoms) );

    /*
     *  If (vaUnits) is not defined then we created
     *  (prPdb.vaUnits) and we should destroy it
     */

    if ( vaUnits == NULL )
        VarArrayDestroy( &(prPdb.vaUnits) );

    VPTRACEEXIT (( "uPdbRead" ));
    return ( prPdb.uUnit );
}

// Search for a bond to the next residue. Return 1 if none found.
static int
writeTER( RESIDUE rRes, int iRes )
{
        int     i, iNextRes = iRes + 1;
        ATOM    aAtom = (ATOM) rRes->aaConnect[CONNECT1];

        if (aAtom == NULL)
                return(1);

        for (i = 0; i < iAtomCoordination(aAtom); i++) {
                ATOM    aChildAtom = aAtomBondedNeighbor(aAtom,i);
                int     iCRes = aChildAtom->iSeenId;
                if (iCRes == iNextRes)
                        return(0);
        }
        return(1);
}

void
PdbWrite( FILE *fOut, UNIT uUnit )
{
        int             i, iResidueCount;
        LOOP            lContents;
        SAVERESIDUEt    *srPResidue;
        PDBWRITEt       pwFile;

        iResidueCount = zUnitIOAmberOrderResidues( uUnit );
        if ( iResidueCount == 0 ) {
                VP0((" no residues\n" ));
                return;
        }

        /*
        **  mark atoms w/ res numbers
        **  used for TER detection (atom->iSeenId), and for default resId (res->iTemp)
        */
        srPResidue = PVAI(uUnit->vaResidues, SAVERESIDUEt, 0);
        for (i = 0; i < iResidueCount; srPResidue++, i++) {
                ATOM    aAtom;
                RESIDUE rRes = srPResidue->rResidue;
                rRes->iTemp = i + 1;
                lContents = lLoop( (OBJEKT)rRes, DIRECTCONTENTSBYSEQNUM );
                while ( (aAtom = (ATOM)oNext(&lContents)) )
                        aAtom->iSeenId = i;
        }

        zPdbFileBegin( &pwFile, fOut );
        /*
         * If we use periodic boundary conditions, write the CRYST1 record
         */
        if (bUnitUseBox(uUnit)) {
            VP0(("   printing CRYST1 record to PDB file with box info\n" ));
            pdb_record p;
            double a, b, c;
            UnitGetBox(uUnit, &a, &b, &c);
            p.pdb.cryst1.a = a;
            p.pdb.cryst1.b = b;
            p.pdb.cryst1.c = c;
            p.pdb.cryst1.alpha = dUnitBeta(uUnit) / DEGTORAD;
            p.pdb.cryst1.beta = dUnitBeta(uUnit) / DEGTORAD;
            p.pdb.cryst1.gamma = dUnitBeta(uUnit) / DEGTORAD;
            p.pdb.cryst1.space_grp[0] = 'P';
            p.pdb.cryst1.space_grp[1] = ' ';
            p.pdb.cryst1.space_grp[2] = '1';
            p.pdb.cryst1.space_grp[3] = '\0';
            p.pdb.cryst1.z = 1;
            p.record_type = PDB_CRYST1;
            pdb_write_record(pwFile.fPdbFile, &p, NULL, 0);
        }

        if ( GDefaults.pdbwritecharges ) {
                pdb_record      p;

                p.record_type = PDB_REMARK;
                p.pdb.remark.num = 1;
                sprintf( p.pdb.remark.text,
                                "LEAP: TEMPERATURE FACTORS ARE CHARGES" );
                pdb_write_record( pwFile.fPdbFile, &p, NULL, 0 );
        }
        srPResidue = PVAI(uUnit->vaResidues, SAVERESIDUEt, 0);
        for (i = 0; i < iResidueCount; srPResidue++, i++) {
                RESIDUE rRes = srPResidue->rResidue;
                zPdbFileWriteContainer( &pwFile, (CONTAINER) rRes);
                if ( writeTER(rRes, i) )
                        zPdbFileWriteTermination( &pwFile );
        }
        zPdbFileEnd( &pwFile );
}


/*
 *      PdbAppendToNameMap
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Append the map entries to the Residue or Atom NameMap.
 *      (lEntries) is a LIST of LISTs where each
 *      sub-LIST contains two or three elements
 *      See zbPdbBuildNameMapEntry.
 *
 *      Print errors if the types are not correct.
 */
static void
PdbAppendToNameMap( DICTIONARY *PSdNameMap, LIST lEntries )
{
LIST            lEntry;
LISTLOOP        llEntries;
ASSOC           aEntry;
int             iMap;
STRING          sKey, sData;
BOOL            bOk;

    if ( iObjectType(lEntries) != LISTid ) {
        DFATAL(( "Need LIST" ));
    }
                /* Loop over all entries in the LIST */
                /* Parse the entries */
    iMap = 0;
    llEntries = llListLoop(lEntries);
    while ( (aEntry = (ASSOC)oListNext(&llEntries)) ) {
        iMap++;
        lEntry = (LIST)oAssocObject(aEntry);
        if ( PSdNameMap == &SdAtomNameMap )
            bOk = zbPdbBuildAtomNameMapEntry( lEntry, iMap, sKey, sData );
        else if ( PSdNameMap == &SdResidueNameMap ) {
            bOk = zbPdbBuildResidueNameMapEntry( lEntry, iMap, sKey, sData );
        } else {
            VPWARN(( "Programming error; skipping map %d %s %s\n" ,
                                        iMap, sKey, sData ));
            VP0(( " map is %d, Res %d Atom %d\n",
                                PSdNameMap, &SdResidueNameMap, &SdAtomNameMap));
            bOk = FALSE;
        }
        if ( bOk ) /* Add the name map entry */
                zPdbNameMapAdd( PSdNameMap, sKey, sData );
    }
}

void
PdbAppendToResMap( LIST lEntries )
{
        PdbAppendToNameMap( &SdResidueNameMap, lEntries );
}

void
PdbAppendToAtomMap( LIST lEntries )
{
        PdbAppendToNameMap( &SdAtomNameMap, lEntries );
}



/*
 *      PdbClearResMap
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Clear the residue name map.
 */
void
PdbClearResMap()
{
    zPdbNameMapDestroy( &SdResidueNameMap );
}

void
PdbClearAtomMap()
{
    zPdbNameMapDestroy( &SdAtomNameMap );
}



/*
 *      PdbDisplayResMap
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Display the Name Map.
 */
void
PdbDisplayResMap()
{
DICTLOOP        dlEntries;
char            *cPData;
char            *cPKey;
STRING          sKey;
int             iTerm;

    if ( SdResidueNameMap == NULL ) {
        VP0(( "The residue name map is empty.\n" ));
        return;
    }

    BasicsResetInterrupt();

                /* Print the main chain definitions */
    dlEntries = ydlDictionaryLoop( SdResidueNameMap );
    while ( (cPData = (char*)yPDictionaryNext( SdResidueNameMap, &dlEntries )) ) {
        cPKey = sDictLoopKey(dlEntries);
        zPdbFromNameMapKey( cPKey, &iTerm, sKey );
        if ( bBasicsInterrupt() ) goto CANCEL;
        if ( iTerm == NOEND ) {
            VP0(( "   %-8s --> %-8s\n", sKey, cPData ));
        }
    }
                /* Print the FIRSTEND chain definitions */
    dlEntries = ydlDictionaryLoop( SdResidueNameMap );
    while ( (cPData = (char*)yPDictionaryNext( SdResidueNameMap, &dlEntries )) ) {
        cPKey = sDictLoopKey(dlEntries);
        zPdbFromNameMapKey( cPKey, &iTerm, sKey );
        if ( bBasicsInterrupt() ) goto CANCEL;
        if ( iTerm == FIRSTEND ) {
            VP0(( " 0 %-8s --> %-8s\n", sKey, cPData ));
        }
    }
                /* Print the LASTEND chain definitions */
    dlEntries = ydlDictionaryLoop( SdResidueNameMap );
    while ( (cPData = (char*)yPDictionaryNext( SdResidueNameMap, &dlEntries )) ) {
        cPKey = sDictLoopKey(dlEntries);
        zPdbFromNameMapKey( cPKey, &iTerm, sKey );
        if ( bBasicsInterrupt() ) goto CANCEL;
        if ( iTerm == LASTEND ) {
            VP0(( " 1 %-8s --> %-8s\n", sKey, cPData ));
        }
    }
    return;
CANCEL:
    VP0(( "Interrupted.\n" ));
    return;
}


/*
 *      PdbDisplayAtomMap
 *
 *      Author: Bill Ross (1994)
 *
 *      Display the atom name map.
 */

void
PdbDisplayAtomMap()
{
DICTLOOP        dlEntries;
char            *cPData;
char            *cPKey;

    if ( SdAtomNameMap == NULL ) {
        VP0(( "The name map is empty.\n" ));
        return;
    }

    BasicsResetInterrupt();

    dlEntries = ydlDictionaryLoop( SdAtomNameMap );
    while ( (cPData = (char*)yPDictionaryNext( SdAtomNameMap, &dlEntries)) ) {
        cPKey = sDictLoopKey(dlEntries);
        if ( bBasicsInterrupt() ) goto CANCEL;
        VP0(( "   %-8s --> %-8s\n", cPKey, cPData ));
    }
    return;
CANCEL:
    VP0(( "Interrupted.\n" ));
    return;
}

