
/*
 *      File:   pdbFile2.c
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

#include        <assert.h>
#ifdef WIN32
#define strcasecmp _stricmp
#else
#include         <strings.h>
#endif

#include        "basics.h"
#include        "classes.h"
#include        "dictionary.h"
#include        "database.h"
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
#include        "pdbFile.h"
#include        "library.h"
#include        "parser.h"
#include        "symmetry.h"
#include        "tools.h" // we should put more into symmetry.c/h
#include        "unitio.h"
#include        "cifparse.h" // for ndb_close_cif()

_Static_assert(IX_DEFAULTKEYLEN >= CONTAINERNAMELEN, "IX_DEFAULTKEYLEN < CONTAINERNAMELEN");


/*
----------------------------------------------------------------------

        Static variables

*/

static DICTIONARY      SdResidueNameMap;
static  DICTIONARY      SdAtomNameMap;

#define MAX_BOND_LENGTH 2.21 // P-P bond, but Si-Si can be 2.33Å
static float dCovalentRadius[] = { 0.73,// use this for default radius
   0.00,                                                                                        0.00,
   0.00, 0.00,                                                    0.82, 0.77, 0.75, 0.73, 0.71, 0.00,
   0.00, 0.00,                                                    0.00, 1.11, 1.06, 1.02, 0.99, 0.00,
   0.00, 0.00, 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00, 0.00, 1.22, 1.19, 1.16, 1.14, 0.00,
   0.00, 0.00, 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00, 0.00, 0.00, 1.38, 1.35, 1.33, 0.00 };
static int iMaxCovalentElement = sizeof(dCovalentRadius) / sizeof(float);

// Possible chainID from best to worst.
// Note that AMBER Mask parsers don't support escape sequences
// with 2-char chainId, and upper+digits+lower we get (26+10+26)*(26+10+26) = 3844 chainId
// If we add special char at column 1 only, we get another 31*(26+10+26) = 1922
const char *GsChainIdList="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                          "012346789"
                          "abcdefghijklmnopqrstuvwxyz"
                          "!\"#$%&()*+,-./:;<=>?[\\]^_`{|}~";

/*
 *----------------------------------------------------------------------
 *
 *        Private routines
 */
static BOOL
bAtomsBondedDist(int iElem1, int iElem2, float dist2, float cutoff) {
    if (iElem1 <= 0 || iElem1 >= iMaxCovalentElement) iElem1 = 0;
    if (iElem2 <= 0 || iElem2 >= iMaxCovalentElement) iElem2 = 0;

    float r1 = dCovalentRadius[iElem1];
    float r2 = dCovalentRadius[iElem2];
    if (!r1 || !r2) return FALSE;

    float cd = (r1 + r2) * cutoff;
    return dist2 < cd * cd;
}

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
            DFATAL("Invalid termination type" );
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
 *
 *      Juno Krahn (2026)
 *      Added resName key for residue-specific atom renaming.
 *      And COMPID_HETID=2 type for COMP_ID to 3-letter HetID
 *      This maps input COMP_ID to AMBER/USER choice of 3-letter HetID code.
 *      HetID to COMP_ID is stored in UNIT.dHeterogens substring COMP_ID:<name>
 */
#define COMPID_HETID  0x10
#define HETID_COMPID  0x11
static void
zPdbToNameMapKey( int iPart, const char *sPart, char *sKey, const char *sResName )
{
    switch ( iPart ) {
        case NOEND:
            if (sResName)  sprintf( sKey, "%s:%s", sResName, sPart );
            else strcpy( sKey, sPart );
            break;
        case FIRSTEND:
            sprintf( sKey, "0:%s", sPart );
            break;
        case LASTEND:
            sprintf( sKey, "1:%s", sPart );
            break;
        case COMPID_HETID:
            sprintf( sKey, "H:%s", sPart );
            break;
        case HETID_COMPID:
            sprintf( sKey, "C:%s", sPart );
            break;
        default:
            DFATAL("Invalid termination type" );
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
    if (sKey[0]=='0' && sKey[1]==':') {
        strcpy( sPart, &(sKey[2]) );
        *iPPart = FIRSTEND;
    } else if (sKey[0]=='1' && sKey[1]==':') {
        strcpy( sPart, &(sKey[2]) );
        *iPPart = LASTEND;
    } else if (sKey[0]=='H' && sKey[1]==':') {
        strcpy( sPart, &(sKey[2]) );
        *iPPart = COMPID_HETID;
    } else if (sKey[0]=='C' && sKey[1]==':') {
        strcpy( sPart, &(sKey[2]) );
        *iPPart = HETID_COMPID;
    } else {
        strcpy( sPart, sKey );
        *iPPart = NOEND;
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

    MESSAGE("Adding %s:%s to name list.\n", sKey, sData );
                /* Delete the old entry if it exists */
    cPOld = (char*)yPDictionaryDelete( *PSdNameMap, sKey );
    if ( cPOld != NULL ) {
        VP0("Substituting map %s -> %s  for  %s -> %s\n",
                sKey, cPNew, sKey, cPOld );
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

    VP0("Clearing name map.\n" );
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
        VPWARN("Map entry %d is not a list. Ignored.\n", iMap );
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
                if ( aD != NULL ) goto BADTYPE;
            }
        }
    }

    if ( aA == NULL || aB == NULL ) goto BADTYPE;
    oA = (OBJEKT)oAssocObject(aA);
    if ( iObjectType(oA) == ODOUBLEid ) {

        if ( aC == NULL ) goto BADTYPE;
        oB = (OBJEKT)oAssocObject(aB);
        oC = (OBJEKT)oAssocObject(aC);
        if ( iObjectType(oB) != OSTRINGid ||
                iObjectType(oC) != OSTRINGid ) goto BADTYPE;
        iInt = (int)dODouble(oA);
        strcpy( sTempKey, sOString(oB) );
        strcpy( sData, sOString(oC) );
        if ( (iInt != 0 && iInt != 1 && iInt != 2) ||
                strlen(sTempKey) == 0 ||
                strlen(sData) == 0 ) goto BADTYPE;
        if ( iInt == 0 ) iInt = FIRSTEND;
        else if ( iInt == 1 ) iInt = LASTEND;
        else if ( iInt == 2 ) iInt = COMPID_HETID;
        else if ( iInt == 3 ) iInt = HETID_COMPID;
        zPdbToNameMapKey( iInt, sTempKey, sKey, NULL );
    } else {
        if ( aC != NULL ) goto BADTYPE;
        oB = (OBJEKT)oAssocObject(aB);
        if ( iObjectType(oA) != OSTRINGid ||
                iObjectType(oB) != OSTRINGid ) goto BADTYPE;
        strcpy( sTempKey, sOString(oA) );
        strcpy( sData, sOString(oB) );
        if ( strlen(sTempKey) == 0 ||
                strlen(sData) == 0 ) goto BADTYPE;
        zPdbToNameMapKey( NOEND, sTempKey, sKey, NULL );
    }

    return(TRUE);

BADTYPE:
    VPWARN("Residue Map entry %d must have the form %s. Ignored\n",
                iMap,
                "{ [0 or 1] string string }" );
    return(FALSE);
}

static BOOL
zbPdbBuildAtomNameMapEntry( LIST lEntry, int iMap, char *sKey, char *sData )
{
STRING          sTempKey;
ASSOC           aA, aB, aC, aD;
OBJEKT          oA, oB;
char            *sResName = NULL;
LISTLOOP        llEntry;

    if ( iObjectType(lEntry) != LISTid ) {
        VPWARN("Map entry %d is not a list. Ignored.\n", iMap );
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
                    aD = (ASSOC)oListNext(&llEntry);
                    if ( aD != NULL ) {
                    goto BADTYPE;
                }
            }
        }
    }

    if ( aA == NULL || aB == NULL ) goto BADTYPE;
    oA = (OBJEKT)oAssocObject(aA);
    oB = (OBJEKT)oAssocObject(aB);
    if ( aC != NULL ) {
        if ( iObjectType(oA) != OSTRINGid ) goto BADTYPE;
        sResName = sOString(oA);
        oA = oB;
        oB = (OBJEKT)oAssocObject(aC);
    }
    if ( iObjectType(oA) != OSTRINGid || iObjectType(oB) != OSTRINGid )
        goto BADTYPE;
    strcpy( sTempKey, sOString(oA) );
    strcpy( sData, sOString(oB) );
    if ( strlen(sTempKey) == 0 || strlen(sData) == 0 )
        goto BADTYPE;

    zPdbToNameMapKey( NOEND, sTempKey, sKey, sResName);

    return(TRUE);

BADTYPE:
    VPWARN("Atom Map entry %d must have the form %s. Ignored\n",
                iMap, "{ string string }" );
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
char *
zcPPdbMapName( DICTIONARY SdNameMap, int iType, const char *sName, const RESIDUE rRes )
{
char            *cPData;
STRING          sKey;

    if ( SdNameMap == NULL ) return(NULL);
    if (rRes) {
        zPdbToNameMapKey( iType, sName, sKey, sContainerName(rRes));
        cPData = (char*)yPDictionaryFind( SdNameMap, sKey );
        if (cPData) return cPData;
        zPdbToNameMapKey( iType, sName, sKey, sResidueTypeNameFromChar(cResidueType(rRes)));
        cPData = (char*)yPDictionaryFind( SdNameMap, sKey );
        if (cPData) return cPData;
    }
    zPdbToNameMapKey( iType, sName, sKey, NULL );
    cPData = (char*)yPDictionaryFind( SdNameMap, sKey );
    return(cPData);
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

    VP1("Matching PDB residue names to LEaP variables.\n" );

                /* If there is no vaUnits VARARRAY then create one */

    if ( !(prPPdb->vaUnits) ) {
        prPPdb->vaUnits = vaVarArrayCreate(sizeof(UNIT));
    }

                /* Now check to see if there are enough */
                /* UNITs in the vaUnits VARARRAY */

    else if ( iVarArrayElementCount(prPPdb->vaUnits) ==
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
                        rnPName->iTerminator, rnPName->sName, NULL );

            if ( cPName ) {
                strcpy( sName, cPName );
                VP1("Mapped residue %s, term: %s, seq. number: %d to: %s.\n",
                        rnPName->sName,
                        zcPPdbTerminationCode( rnPName->iTerminator), i,sName);
            } else {
                strcpy( sName, rnPName->sName );
                MESSAGE("(Residue %d: %s, %s, was not found in name map.)\n",
                    i, sName, zcPPdbTerminationCode( rnPName->iTerminator));
            }
            uUnit = (UNIT)oVariable(sName);
            if ( uUnit == NULL) {
                VPWARN("Unknown residue: %s   number: %d   type: %s\n",
                        sName, i, zcPPdbTerminationCode(rnPName->iTerminator));
                if ( rnPName->iTerminator != NOEND ) {
                    VP0("..relaxing end constraints to try for a dbase match\n");
                    cPName = zcPPdbMapName( SdResidueNameMap, NOEND,
                                                        rnPName->sName, NULL );
                    if ( cPName ) {
                        strcpy( sName, cPName );
                        uUnit = (UNIT)oVariable(sName);
                    }

                    if ( uUnit == NULL )
                        VPWARN("  -no luck\n" );
                    else
                        VP0("  -matched to non-end type residue %s\n",
                                                        sName );
                }
            } else if ( iObjectType(uUnit) != UNITid ) {
                VPWARN("Invalid unit: %s   sequence number: %d\n", sName, i );
            }
            *PVAI(prPPdb->vaUnits,UNIT,i) = uUnit;
        }
    } else {
        VP0("There are %d units specified in the command,\n",
                iVarArrayElementCount(prPPdb->vaUnits) );
        VP0("and only %d residues within the PDB file.\n",
                iVarArrayElementCount(prPPdb->vaResidues) );
        VPWARN("The extra units will be ignored.\n" );
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
        VP0("Creating new UNIT for residue: %s sequence: %d\n",
                PVAI(prPPdb->vaResidues,RESIDUENAMEt,iCurUnit)->sName, iResNum);
        uUnit = (UNIT)oCreate(UNITid);
        rRes = (RESIDUE)oCreate(RESIDUEid);
        ContainerSetName( (CONTAINER)rRes,
                PVAI(prPPdb->vaResidues,RESIDUENAMEt,iCurUnit)->sName );
        ContainerAdd( (CONTAINER)uUnit, (OBJEKT)rRes );
    } else {
        MESSAGE("Getting UNIT: %s\n", sContainerName((CONTAINER)uOrig) );
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
void
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

            MESSAGE("Building externals from: %s\n",
                        sContainerFullDescriptor( (CONTAINER)aSpan, sSpan ) );

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
                    MESSAGE("Building externals from: %s\n",
                        sContainerFullDescriptor( (CONTAINER)aSpan, sSpan ) );

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
 *      zPdbCreateSymmetryRelatedMonomers
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Use the transform matrices to build symmetry related
 *      monomers of the UNIT in (uUnit).
 *      Combine all of the monomers into a single UNIT.
 * -----------------------------------------------------------------------
 *              Juno Krahn (2026)
 * If crystal spacegroup symmetry expansion: For each symmetry transform,
 * compute the transformed center-of-mass, wrap it into the primary cell
 * (frac [0,1)) and apply the resulting Cartesian translation to all atoms
 * so the molecule center lands inside the unit cell.
 * --------------------------------------------------------------------- */
static void
zPdbCreateSymmetryRelatedMonomers( PDBREADt *prPPdb )
{
int             iTransforms, iLast, i;
VECTOR          vPos, vCen, vCenFrac, vCenWrap, vShift = {0,0,0};
MATRIX          mTransform, M, Mi;
UNIT            uOrig, uCopy;
ATOM            aAtom;
LOOP            lAtoms;

    iTransforms = iVarArrayElementCount( prPPdb->vaMatrices );
    iLast = -1;
    for ( i = 0; i < iTransforms; i++ ) {
        if ( PVAI(prPPdb->vaMatrices,PDBMATRIXt,i)->bUsed )
            iLast = i;
    }
    if ( iLast == -1 ) return;

    /* --- build fractional transforms from the unit cell --- */
    BuildFractionalTransforms( prPPdb->uUnit, M, Mi );

    /* --- get center of original monomer once, before any transforms --- */
    vCen = vContainerGeometricCenter( (CONTAINER)prPPdb->uUnit );

    uOrig = (UNIT)oCopy((OBJEKT)prPPdb->uUnit);

    for ( i = 0; i < iTransforms; i++ ) {
        if ( !PVAI(prPPdb->vaMatrices,PDBMATRIXt,i)->bUsed ) continue;

        MatrixCopy( mTransform,
                    PVAI(prPPdb->vaMatrices,PDBMATRIXt,i)->mTransform );

        if ( i != iLast )
            uCopy = (UNIT)oCopy((OBJEKT)uOrig);
        else
            uCopy = uOrig;

        if (GDefaults.bPdbExpandSymm) {
        /* --- transform center, wrap into primary cell, compute shift --- */
        MatrixTimesVector( vCen, mTransform, vCen );

        vCenFrac.dX = Mi[0][0]*vCen.dX + Mi[1][0]*vCen.dY + Mi[2][0]*vCen.dZ;
        vCenFrac.dY = Mi[0][1]*vCen.dX + Mi[1][1]*vCen.dY + Mi[2][1]*vCen.dZ;
        vCenFrac.dZ = Mi[0][2]*vCen.dX + Mi[1][2]*vCen.dY + Mi[2][2]*vCen.dZ;

        vCenWrap.dX = vCenFrac.dX - floor(vCenFrac.dX);
        vCenWrap.dY = vCenFrac.dY - floor(vCenFrac.dY);
        vCenWrap.dZ = vCenFrac.dZ - floor(vCenFrac.dZ);

        vShift.dX = M[0][0]*vCenWrap.dX + M[1][0]*vCenWrap.dY + M[2][0]*vCenWrap.dZ - vCen.dX;
        vShift.dY = M[0][1]*vCenWrap.dX + M[1][1]*vCenWrap.dY + M[2][1]*vCenWrap.dZ - vCen.dY;
        vShift.dZ = M[0][2]*vCenWrap.dX + M[1][2]*vCenWrap.dY + M[2][2]*vCenWrap.dZ - vCen.dZ;
        }

        /* --- transform all atoms and apply wrapping shift in one pass --- */
        lAtoms = lLoop( (OBJEKT)uCopy, ATOMS );
        while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
            MatrixTimesVector( vPos, mTransform, vAtomPosition(aAtom) );
            vPos.dX += vShift.dX;
            vPos.dY += vShift.dY;
            vPos.dZ += vShift.dZ;
            AtomSetPosition( aAtom, vPos );
        }

        VP1("Building symmetry related monomer %d.\n", i+1);
        UnitJoin( prPPdb->uUnit, uCopy );
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
    MESSAGE("Element: |%s|   pdb_name=|%s|\n", sElement, sPName );

    if (rRes) {
        strncpy(p.pdb.atom.residue.chain_id,rRes->sChainId,2);
        p.pdb.atom.residue.chain_id[2]=0;
        p.pdb.atom.residue.seq_num = rRes->iPdbResSeq;
        p.pdb.atom.residue.insert_code = rRes->cICode;
    } else {
        strcpy(p.pdb.atom.residue.chain_id,"");
        p.pdb.atom.residue.seq_num = pwPFile->iResidueSeq;
        p.pdb.atom.residue.insert_code=' ';
    }
    p.pdb.atom.alt_loc = ' ';
    p.pdb.atom.x = dVX(&vAtomPosition(aAtom));
    p.pdb.atom.y = dVY(&vAtomPosition(aAtom));
    p.pdb.atom.z = dVZ(&vAtomPosition(aAtom));
    p.pdb.atom.occupancy = bAtomFlagsSet(aAtom,ATOMBULKSOLVENT) ? 0.0 : 1.0;
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
    pdb_write_record( pwPFile->fPdbFile, &p);

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

    memcpy(p.pdb.ter.residue.chain_id,"  ",2); // FIXME do real values?
    p.pdb.ter.residue.seq_num = pwPFile->iResidueSeq;
    p.pdb.ter.residue.insert_code=' ';
    p.record_type = PDB_TER;
    pdb_write_record( pwPFile->fPdbFile, &p);
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
    pdb_write_record( pwPFile->fPdbFile, &p);
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
    char *cPHetID = zcPPdbMapName( SdResidueNameMap, COMPID_HETID, cPTemp, NULL );
    if (cPHetID) {
        VP2("Mapped Residue Comp_ID %s to HetID %s\n",cPTemp,cPHetID);
        cPTemp = cPHetID;
    }
    int iLen = strlen(cPTemp);
    if ( iLen <= RESIDUE_NAME_LENGTH ) {
        strcpy( pwPFile->sResidueName, cPTemp);
    } else {
        // Truncate long names by removing leading characters
        strcpy( pwPFile->sResidueName, cPTemp + iLen-RESIDUE_NAME_LENGTH);
        if ( isalpha( *(cPTemp + 3) ) &&
                isalpha( *(cPTemp + 2) ) && isalpha( *(cPTemp + 1) ) &&
                (*cPTemp == N_TERMINAL_PREFIX || *cPTemp == C_TERMINAL_PREFIX) ) {
            /* Probable N-terminal or C-terminal amino acid name */
            VPWARN(" Converting %c-terminal residue name to PDB format: %s -> %s\n",
                        *cPTemp, sContainerName(cCont), pwPFile->sResidueName );
        } else {
            // Left-truncate optimal for long COMP_ID name.
            // Also works for N-, C- prefix form of terminal naming scheme
            VPWARN(" Truncating residue name for PDB format: %s -> %s\n",
                        sContainerName(cCont), pwPFile->sResidueName );
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

/*
==========================================================================

        Public routines
*/


void
PdbWrite( FILE *fOut, UNIT uUnit )
{
        int             iResidueCount;
        LOOP            lContents;
        SAVERESIDUEt    *srPResidue;
        PDBWRITEt       pwFile;
        LOOP            lAtoms;
        ATOM            aAtom;
        ATOM            aNeighbor;
        RESIDUE         rNeighbor;
        int             iSSBOND;
        DICTLOOP        dlHeterogens;
        pdb_record      p;
        IX_REC          recResCount = { NULL }; // recptr is unused
        int             status;

        // Allocates unit->vaResidues
        iResidueCount = zUnitIOAmberOrderResidues( uUnit );
        if ( iResidueCount == 0 ) {
                VP0(" no residues\n" );
                if (uUnit->vaResidues) VarArrayDestroy(&uUnit->vaResidues);
                return;
        }

        /*
        **  mark residues and atoms atoms with sequential residue numbers
        **  used for TER detection (atom->iSeenId), and for default resSeq (res->iTemp)
        */

        create_index( &pwFile.ixResCount, IX_DUPKEY, IX_DEFAULTKEYLEN);

        srPResidue = PVAI(uUnit->vaResidues, SAVERESIDUEt, 0);
        for (int i = 0; i < iResidueCount; srPResidue++, i++) {
                RESIDUE rRes = srPResidue->rResidue;
                strcpy(recResCount.key,sContainerName(rRes));
                add_key(&recResCount, &pwFile.ixResCount);
                rRes->iTemp = i + 1;
                lContents = lLoop( (OBJEKT)rRes, DIRECTCONTENTSBYSEQNUM );
                while ( (aAtom = (ATOM)oNext(&lContents)) )
                        aAtom->iSeenId = i;
        }

        first_key(&pwFile.ixResCount);
        status = next_key(&recResCount, &pwFile.ixResCount);
        while (status == IX_OK ) {
            int ilen = strlen(recResCount.key);
            char *cPTemp = recResCount.key;
            if ( isalpha( *(cPTemp + 3) ) &&
                isalpha( *(cPTemp + 2) ) && isalpha( *(cPTemp + 1) ) &&
                (*cPTemp == 'N' || *cPTemp == 'C' ) ) {
                // pass
            } else {
                IX_REC recHetID = { NULL };
                char *cPHetID = zcPPdbMapName( SdResidueNameMap, COMPID_HETID, recResCount.key, NULL );
                if (cPHetID) {
                    strcpy(recHetID.key,cPHetID);
                    if (has_key(&recHetID, &pwFile.ixResCount)) {
                        VPFATAL("HetID=%s for CompID ResName=%s conflicts with other residue(s)\n",
                                "Define PDB ResName with: addPdbResMap { 2 \"%s\" \"<resName>\" }\n",
                                recHetID.key, recResCount.key, recHetID.key);
                    }
                } else if (ilen>3) {
                    // Default to left-truncated resName
                    strcpy(recHetID.key,recResCount.key+ilen-3);
                    if (has_key(&recHetID, &pwFile.ixResCount)) {
                        VPFATAL("Default HetID=%s for CompID=%s conflicts with existing residue\n"
                                "Define PDB ResName with: addPdbResMap { 2 \"%s\" \"<resName>\" }\n",
                                recHetID.key, recResCount.key, recHetID.key);
                    }
                }
            }
            status = next_key(&recResCount, &pwFile.ixResCount);
        }

        zPdbFileBegin( &pwFile, fOut );

        if (iDictionaryElementCount(uUnit->dHeterogens) != 0 ) {

            // FIXME: TODO: split descripton at first ';' into HETNAM+HETSYN
            dlHeterogens = ydlDictionaryLoop( uUnit->dHeterogens );
            while ( yPDictionaryNext( uUnit->dHeterogens, &dlHeterogens ) ) {
                p.pdb.hetnam.serial_num=0;
                STRING sDesc;
                strcpy(sDesc, (char*)PDictLoopData(dlHeterogens));
                char *cPHetID = sDictLoopKey(dlHeterogens);
                char *cPCompID = zcPPdbMapName( SdResidueNameMap, HETID_COMPID, cPHetID, NULL );
                //VP2("HETNAM: HETID=%s COMPID=%s\n",cPHetID,cPCompID);
                if (cPCompID && !strstr(sDesc,"COMP_ID:")) {
                    strcat(sDesc,"; COMP_ID:");
                    strcat(sDesc,cPCompID);
                }
                char *cPSep=strchr(sDesc,';');
                if (cPSep) *cPSep =0;
                if (cPCompID) strcpy(p.pdb.hetnam.extra,cPCompID);
                else p.pdb.hetnam.extra[0]=0;
                p.record_type = PDB_HETNAM;
                p.pdb.hetsyn.serial_num=0;
                assert(strlen(cPHetID)<=3); // FIXME: shouldn't happen?
                strcpy(p.pdb.hetnam.het_id,cPHetID);
                pdb_write_multiline(pwFile.fPdbFile, &p, &p.pdb.hetnam.serial_num,
                        p.pdb.hetnam.desc, sizeof(p.pdb.hetnam.desc), sDesc);
                if (cPSep) {
                    p.record_type = PDB_HETSYN;
                    p.pdb.hetsyn.serial_num=0;
                    if (*(++cPSep) == ' ') cPSep++;
                    pdb_write_multiline(pwFile.fPdbFile, &p, &p.pdb.hetsyn.serial_num,
                            p.pdb.hetsyn.desc, sizeof(p.pdb.hetsyn.desc), cPSep);
                }
            }
        }

        // Determine and write LINK, SSBOND records
        if (GDefaults.bPdbKeepChainId) {
            iSSBOND = 0;
            lAtoms = lLoop( (OBJEKT)uUnit, ATOMS );
            while ( (aAtom = (ATOM)oNext(&lAtoms)) ) {
                RESIDUE rResidue = (RESIDUE)cContainerWithin((CONTAINER) aAtom);
                int iResSeq = iContainerSequence((CONTAINER) rResidue);
                for ( int i=0; i<iAtomCoordination(aAtom); i++ ) {
                    aNeighbor = aAtomBondedNeighbor( aAtom, i );
                    rNeighbor = (RESIDUE)cContainerWithin((CONTAINER) aNeighbor);
                    if (iContainerSequence(cContainerWithin((CONTAINER) aNeighbor))-iResSeq > 1 ) {
                        // forward LINK or SSBOND
                        //
                        // Include leading space in atom name as for ATOM/HETAM
                        int iElement, pad1 = 1, pad2 = 1;
                        if (strlen(sContainerName(aAtom))>3) pad1=0;
                        else {
                            iElement = iAtomElement(aAtom);
                            if (iElement != NOELEMENT) pad1 = GeaElements[iElement].sName[1] == 0;
                        }
                        if (strlen(sContainerName(aNeighbor))>3) pad2=FALSE;
                        else {
                            iElement = iAtomElement(aNeighbor);
                            if (iElement != NOELEMENT) pad2 = GeaElements[iElement].sName[1] == 0;
                        }

                        strcpy(p.pdb.link.residues[0].name,sContainerName(rResidue));
                        strcpy(p.pdb.link.residues[0].chain_id,sResidueChainId(rResidue));
                        p.pdb.link.residues[0].seq_num=iResiduePdbSequence(rResidue);
                        p.pdb.link.residues[0].insert_code=' ';
                        strcpy(p.pdb.link.symop[0],"1555");

                        strcpy(p.pdb.link.residues[1].name,sContainerName(rNeighbor));
                        strcpy(p.pdb.link.residues[1].chain_id,sResidueChainId(rNeighbor));
                        p.pdb.link.residues[1].seq_num=iResiduePdbSequence(rNeighbor);
                        p.pdb.link.residues[1].insert_code=' ';
                        strcpy(p.pdb.link.symop[1],"1555");

                        double dX = aAtom->vPosition.dX - aNeighbor->vPosition.dX;
                        double dY = aAtom->vPosition.dY - aNeighbor->vPosition.dY;
                        double dZ = aAtom->vPosition.dZ - aNeighbor->vPosition.dZ;
                        p.pdb.link.distance = sqrt(dX*dX + dY*dY + dZ*dZ);

                        if (!strcmp(p.pdb.link.residues[0].name,"CYS") && !strcmp(p.pdb.link.residues[1].name,"CYS") &&
                            !strcmp(p.pdb.link.name[0],"SG") && !strcmp(p.pdb.link.name[1],"SG") ) {
                            p.pdb.ssbond.seq_num = ++iSSBOND;
                            p.record_type = PDB_SSBOND;
                        } else {
                            p.pdb.link.name[0][0]=' ';
                            strcpy(&p.pdb.link.name[0][pad1],sContainerName(aAtom));
                            p.pdb.link.alt_loc[0]=' ';
                            p.pdb.link.name[1][0]=' ';
                            strcpy(&p.pdb.link.name[1][pad2],sContainerName(aNeighbor));
                            p.pdb.link.alt_loc[1]=' ';
                            p.record_type = PDB_LINK;
                        }

                        pdb_write_record(pwFile.fPdbFile, &p);
                    }
                }
            }
        }
        /*
         * If we use periodic boundary conditions, write the CRYST1 record
         */
        if (bUnitUseBox(uUnit)) {
            VP0("   printing CRYST1 record to PDB file with box info\n" );
            double a, b, c;
            UnitGetBox(uUnit, &a, &b, &c);
            p.pdb.cryst1.a = a;
            p.pdb.cryst1.b = b;
            p.pdb.cryst1.c = c;
            p.pdb.cryst1.alpha = uUnit->dAlpha / DEGTORAD;
            p.pdb.cryst1.beta = dUnitBeta(uUnit) / DEGTORAD;
            p.pdb.cryst1.gamma = uUnit->dGamma / DEGTORAD;
            strcpy(p.pdb.cryst1.space_grp, "P 1");
            p.pdb.cryst1.z = 1;
            p.record_type = PDB_CRYST1;
            pdb_write_record(pwFile.fPdbFile, &p);
        }

        if ( GDefaults.pdbwritecharges ) {
                p.record_type = PDB_REMARK;
                p.pdb.remark.num = 1;
                sprintf( p.pdb.remark.text,
                                "LEAP: TEMPERATURE FACTORS ARE CHARGES" );
                pdb_write_record( pwFile.fPdbFile, &p);
        }
        srPResidue = PVAI(uUnit->vaResidues, SAVERESIDUEt, 0);
        for (int i = 0; i < iResidueCount; srPResidue++, i++) {
                RESIDUE rRes = srPResidue->rResidue;
                zPdbFileWriteContainer( &pwFile, (CONTAINER) rRes);
                if ( writeTER(rRes, i) )
                        zPdbFileWriteTermination( &pwFile );
        }
        zPdbFileEnd( &pwFile );
        VarArrayDestroy(&uUnit->vaResidues);
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
        DFATAL("Need LIST" );
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
            VPWARN("Programming error; skipping map %d %s %s\n" ,
                                        iMap, sKey, sData );
            VP0(" map is %d, Res %d Atom %d\n",
                                PSdNameMap, &SdResidueNameMap, &SdAtomNameMap);
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

// Given a PdbMap target name, get FIRST,LAST,NO END flags
int
fGetPdbResMapped(char *sResName) {
    int fResult=0;
    if (!SdResidueNameMap) return 0;
    DICTLOOP dlLoop = ydlDictionaryLoop( SdResidueNameMap );
    while ( yPDictionaryNext( SdResidueNameMap, &dlLoop )) {
        char *cPData = (char*)PDictLoopData(dlLoop);
        if ( strcmp(cPData,sResName) ) continue;
        UNIT uTemplate = (UNIT)oVariable(cPData);
        // Residue maps to undefined string
        if ( uTemplate == NULL) continue;
        char *cPKey = sDictLoopKey(dlLoop);
        int iTerm;
        STRING sKey;
        zPdbFromNameMapKey( cPKey, &iTerm, sKey );
        if (iTerm == FIRSTEND) fResult |= RESIDUEFIRSTEND;
        else if (iTerm == LASTEND) fResult |= RESIDUELASTEND;
        else if (iTerm != COMPID_HETID && iTerm != HETID_COMPID) fResult |= RESIDUENOEND;
    }
    return fResult;
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
        VP0("The residue name map is empty.\n" );
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
            VP0("   %-8s --> %-8s\n", sKey, cPData );
        }
    }
                /* Print the FIRSTEND chain definitions */
    dlEntries = ydlDictionaryLoop( SdResidueNameMap );
    while ( (cPData = (char*)yPDictionaryNext( SdResidueNameMap, &dlEntries )) ) {
        cPKey = sDictLoopKey(dlEntries);
        zPdbFromNameMapKey( cPKey, &iTerm, sKey );
        if ( bBasicsInterrupt() ) goto CANCEL;
        if ( iTerm == FIRSTEND ) {
            VP0(" 0 %-8s --> %-8s\n", sKey, cPData );
        }
    }
                /* Print the LASTEND chain definitions */
    dlEntries = ydlDictionaryLoop( SdResidueNameMap );
    while ( (cPData = (char*)yPDictionaryNext( SdResidueNameMap, &dlEntries )) ) {
        cPKey = sDictLoopKey(dlEntries);
        zPdbFromNameMapKey( cPKey, &iTerm, sKey );
        if ( bBasicsInterrupt() ) goto CANCEL;
        if ( iTerm == LASTEND ) {
            VP0(" 1 %-8s --> %-8s\n", sKey, cPData );
        }
    }
                /* Print the COMPID_HETID definitions */
    dlEntries = ydlDictionaryLoop( SdResidueNameMap );
    while ( (cPData = (char*)yPDictionaryNext( SdResidueNameMap, &dlEntries )) ) {
        cPKey = sDictLoopKey(dlEntries);
        zPdbFromNameMapKey( cPKey, &iTerm, sKey );
        if ( bBasicsInterrupt() ) goto CANCEL;
        if ( iTerm == COMPID_HETID ) {
            VP0(" COMPID %-8s --> HETID %-8s\n", sKey, cPData );
        }
    }
    dlEntries = ydlDictionaryLoop( SdResidueNameMap );
    while ( (cPData = (char*)yPDictionaryNext( SdResidueNameMap, &dlEntries )) ) {
        cPKey = sDictLoopKey(dlEntries);
        zPdbFromNameMapKey( cPKey, &iTerm, sKey );
        if ( bBasicsInterrupt() ) goto CANCEL;
        if ( iTerm == HETID_COMPID ) {
            VP0(" HETID %-8s --> COMPID %-8s\n", sKey, cPData );
        }
    }
    return;
CANCEL:
    VP0("Interrupted.\n" );
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
        VP0("The name map is empty.\n" );
        return;
    }

    BasicsResetInterrupt();

    dlEntries = ydlDictionaryLoop( SdAtomNameMap );
    while ( (cPData = (char*)yPDictionaryNext( SdAtomNameMap, &dlEntries)) ) {
        cPKey = sDictLoopKey(dlEntries);
        if ( bBasicsInterrupt() ) goto CANCEL;
        VP0("   %-8s --> %-8s\n", cPKey, cPData );
    }
    return;
CANCEL:
    VP0("Interrupted.\n" );
    return;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
 *      zPdbReadFile (from zPdbReadScan)
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
zPdbReadFile( PDBREADt *prPRead )
{
pdb_record      p;
int             iPdbSequence; // previous resSeq
BOOL            bLastReadPdbRecordWasTer = FALSE;
BOOL            bNewChain, bNewRes;
RESIDUENAMEt    rnName;
ATOMNAMEt       anAtom;
struct pdb_link pdbLink;
int             iSerial, iStart, iRow;
MATRIX          *mPMatrix;
int             iTerm, iLast, iSerialNum, iSerialNumMax;
char            c2CurrChain[2]={0}, cInsertionCode= ' ';
int             iMultipleResName = 0;
int             iResIdMax, iAtomSerialMax;
int             iPdbResSeqOffset=0, iAtomSerialOffset=0;
int             iCurrentModel=0;
int             iCurrentBioMT=0;
    VPTRACEENTER("zPdbReadFile" );

    bNewChain = TRUE;
    iSerialNum = 0, iSerialNumMax = 0; // Tracks previous and maximum atomSerial in FILE
    iAtomSerialMax = HY36_WIDTH_5_MAX; // MAX allowed in PDB format
    iResIdMax = HY36_WIDTH_4_MAX;

    do {
        p = pdb_read_record(prPRead->fPdbFile);

                /* Process the records */

        switch ( p.record_type ) {
            case PDB_CONECT:
                if (GDefaults.bPdbUseConect)
                    VarArrayAdd( prPRead->vaConectRecs, (GENP)&p.pdb.conect );
                break;

            case PDB_TITLE:
                strcpy(prPRead->uUnit->sDescription,p.pdb.title.title);
                StringTrim(prPRead->uUnit->sDescription);
                break;

            case PDB_HETNAM:
                {
                    StringTrim(p.pdb.hetnam.het_id);
                    STRING *desc = yPDictionaryFind(prPRead->uUnit->dHeterogens,p.pdb.hetnam.het_id);
                    if (!desc) {
                         MALLOC( desc, STRING *, sizeof(STRING) );
                         DictionaryAdd(prPRead->uUnit->dHeterogens, p.pdb.hetnam.het_id, desc);
                    }
                    if (!p.pdb.hetnam.serial_num) {
                        StringTrim(p.pdb.hetnam.desc);
                        strcpy(*desc,p.pdb.hetnam.desc);
                    } else {
                        StringRTrim(p.pdb.hetnam.desc);
                        strcat(*desc,p.pdb.hetnam.desc);
                    }
                }
                if (p.pdb.hetnam.extra[0] != 0 && p.pdb.hetnam.extra[0] != ' ') {
                    // Long to short name mapping: addPdbResMap { 2 COMP_ID HetID }
                    char sKey[32];
                    zPdbToNameMapKey(COMPID_HETID, p.pdb.hetnam.het_id, sKey, NULL);
                    zPdbNameMapAdd( &SdResidueNameMap, sKey, strdup(p.pdb.hetnam.extra));
                    zPdbToNameMapKey(HETID_COMPID, p.pdb.hetnam.het_id, sKey, NULL);
                    zPdbNameMapAdd( &SdResidueNameMap, sKey, strdup(p.pdb.hetnam.het_id));
                }
                break;

            case PDB_HETSYN:
                {
                    StringTrim(p.pdb.hetsyn.het_id);
                    STRING *desc = yPDictionaryFind(prPRead->uUnit->dHeterogens,p.pdb.hetsyn.het_id);
                    // There should ALWAYS be a HETNAM before HETSYN
                    if (!desc) {
                        VPWARN("HETSYN without matching HETNAME will be ignored: %s\n",p.pdb.hetsyn.het_id);
                        break;
                    }
                    // Add a ; seperator if it's the first HETSYN to append onto HETNAM
                    if (!p.pdb.hetsyn.serial_num) {
                        StringTrim(p.pdb.hetsyn.desc);
                        strcat(*desc,"; ");
                        strcat(*desc,p.pdb.hetsyn.desc);
                    } else {
                        StringRTrim(p.pdb.hetsyn.desc);
                        strcpy(*desc,p.pdb.hetsyn.desc);
                    }
                    char *pCompId = strstr(*desc,"COMP_ID:");
                    if (pCompId) {
                        char sCompId[32], sKey[32];
                        strncpy(sCompId,pCompId+strlen("COMP_ID:"),32);
                        pCompId = strtok(sCompId,"; \n\r");
                        VP0("Found COMP_ID=%s for HetID=%s\n",pCompId,p.pdb.hetsyn.het_id);
                        if ( strlen(pCompId) < 3 || strlen(pCompId) > 5 ) {
                            VPWARN("COMPID length unexpected size: \"%s\"\n",pCompId);
                        } else {
                            zPdbToNameMapKey(HETID_COMPID, p.pdb.hetsyn.het_id, sKey, NULL);
                            zPdbNameMapAdd( &SdResidueNameMap, sKey, strdup(pCompId));
                            zPdbToNameMapKey(COMPID_HETID, pCompId, sKey, NULL);
                            zPdbNameMapAdd( &SdResidueNameMap, sKey, strdup(p.pdb.hetsyn.het_id));
                        }
                    }
                }
                break;

            case PDB_ATOM:
            case PDB_HETATM:
                if (p.pdb.atom.alt_loc != ' ' && p.pdb.atom.alt_loc != GDefaults.cPdbAltLocSelect) break;

                if (iCurrentModel) {
                    if (!GDefaults.iPdbReadModel) GDefaults.iPdbReadModel = iCurrentModel;
                    if (GDefaults.iPdbReadModel > 0 && GDefaults.iPdbReadModel != iCurrentModel) break;
                    if (GDefaults.iPdbReadModel < 0 && p.pdb.atom.residue.chain_id[0]==' ') {
                        int i = iCurrentModel-1;
                        if (i < CHAINID_LIST_LEN) {
                            p.pdb.atom.residue.chain_id[0] = GsChainIdList[i];
                        } else VPFATAL("ChainID overflow for iReadModel < 0");
                    }
                }
                if (p.pdb.atom.residue.chain_id[0]==' ') {
                    p.pdb.atom.residue.chain_id[0] = p.pdb.atom.residue.chain_id[1];
                    p.pdb.atom.residue.chain_id[1] = 0;
                }
                if (p.pdb.atom.leap_expanded) {
                    iAtomSerialMax = 99999;
                    iResIdMax = 9999;
                }
                if (p.pdb.atom.serial_num == 0 && iSerialNum == iAtomSerialMax)
                    iAtomSerialOffset += iAtomSerialMax+1; // atomSerial overflow unwrapping
                iSerialNum = p.pdb.atom.serial_num;
                p.pdb.atom.serial_num += iAtomSerialOffset;
                if ( p.pdb.atom.serial_num > iSerialNumMax )
                    iSerialNumMax = p.pdb.atom.serial_num;

                VPTRACE("Read PDB_ATOM or PDB_HETATM record.\n" );
                VPTRACE("bNewChain = %d; Chain = %.2s; iPdbSequence = %d; "
                          "iSerialNum = %d\n", bNewChain, p.pdb.atom.residue.chain_id,
                          iPdbSequence, p.pdb.atom.serial_num );
                /*
                 *  allow for right-shifted resnames XXX right justified is standard
                 */
                char *cPResName = p.pdb.atom.residue.name;
                while (*cPResName == ' ') cPResName++;

                // Convert 3-letter HetID to global CompID if defined
                char *cPCompID = zcPPdbMapName( SdResidueNameMap, COMPID_HETID, cPResName, NULL );
                if (cPCompID) {
                    VP2("Mapped hetID %s to CompId %s\n", cPResName, cPCompID);
                    cPResName = cPCompID;
                }

                iTerm = NOEND;
                bNewChain = bNewChain || memcmp(c2CurrChain,p.pdb.atom.residue.chain_id,2);
                if ( bNewChain ) {
                    VP2(" (starting new molecule for chain \"%.2s\")\n", p.pdb.atom.residue.chain_id );
                    iTerm = FIRSTEND;
                    iLast = iVarArrayElementCount( prPRead->vaResidues );
                    if (iLast > 0) PVAI( (prPRead->vaResidues), RESIDUENAMEt,iLast-1)->iTerminator = LASTEND;
                    memcpy(c2CurrChain,p.pdb.atom.residue.chain_id,2);
                    if (!c2CurrChain[0]) c2CurrChain[1]=0; // double NUL for blank chain
                    iPdbResSeqOffset = 0;
                }
                bNewRes = (bNewChain || p.pdb.atom.residue.seq_num != iPdbSequence ||
                            p.pdb.atom.residue.insert_code != cInsertionCode ||
                            bLastReadPdbRecordWasTer );
                if (!bNewRes && strcmp(rnName.sName, cPResName) ) {
                    /* First detection of a residue with the same sequence
                     * number and insertion code but a different name.
                     */
                    VPWARN("Name change in pdb file residue %.2s %d%c;\n"
                        "this residue is split into %s and %s.\n",
                        c2CurrChain, iPdbSequence, cInsertionCode, rnName.sName, cPResName);
                    bNewRes = TRUE;
                    iMultipleResName++;
                }
                if (bNewRes) {
                    VPTRACE("Detected a new residue.\n" );
                    if (p.pdb.atom.residue.seq_num == 0 && iPdbSequence == iResIdMax)
                         iPdbResSeqOffset += iResIdMax + 1;
                    rnName.iPdbSequence = p.pdb.atom.residue.seq_num + iPdbResSeqOffset; // ResSeq overflow unwrapping
                    rnName.iTerminator = iTerm;
                    strcpy( rnName.sChainId, p.pdb.atom.residue.chain_id);
                    strcpy( rnName.sName, cPResName );
                    rnName.iCode = p.pdb.atom.residue.insert_code;

                    rnName.iFirstAtom = iVarArrayElementCount(prPRead->vaAtomRecs); // zero based array

                    anAtom.iResNameIndex = iVarArrayElementCount( prPRead->vaResidues );
                    VarArrayAdd( prPRead->vaResidues, (GENP)&rnName );
                    bLastReadPdbRecordWasTer = FALSE;

                    MESSAGE("Reading residue: <%s>\n", rnName.sName );
                    iPdbSequence = p.pdb.atom.residue.seq_num;

                }
                bNewChain = FALSE;

                if (p.pdb.atom.name[0]== ' ')
                    strcpy(anAtom.sName, p.pdb.atom.name+1);
                else
                    strcpy(anAtom.sName, p.pdb.atom.name);
                if (p.pdb.atom.element[1]!=0)
                    anAtom.iElement = iPdbElementNumber(p.pdb.atom.element);
                else {
                    int elem = iElementNumberFromAmber(anAtom.sName);
                    anAtom.iElement = (elem == NOELEMENT) ? 0 : elem;
                }
                anAtom.iAtomSerial = p.pdb.atom.serial_num;
                anAtom.x = p.pdb.atom.x;
                anAtom.y = p.pdb.atom.y;
                anAtom.z = p.pdb.atom.z;
                VarArrayAdd( prPRead->vaAtomRecs, (GENP)&anAtom );

                break;

            case PDB_TER:
                VPTRACE("Read PDB_TER record.\n" );
                bLastReadPdbRecordWasTer = TRUE;
                        /* If you read a TER card then make the */
                        /* last RESIDUE read a terminating RESIDUE */
                iLast = iVarArrayElementCount( prPRead->vaResidues );
                if (iLast > 0) PVAI( (prPRead->vaResidues), RESIDUENAMEt,iLast-1)->iTerminator = LASTEND;
                bNewChain = TRUE;
                break;

            case PDB_CRYST1:
                UnitSetBox(prPRead->uUnit,p.pdb.cryst1.a, p.pdb.cryst1.b, p.pdb.cryst1.c);
                prPRead->uUnit->dAlpha = p.pdb.cryst1.alpha * DEGTORAD;
                prPRead->uUnit->dBeta = p.pdb.cryst1.beta * DEGTORAD;
                prPRead->uUnit->dGamma = p.pdb.cryst1.gamma * DEGTORAD;
                UnitSetUseBox(prPRead->uUnit, TRUE );
                STRING sDesc;
                sprintf(sDesc,"SG:%s",p.pdb.cryst1.space_grp);
                if (sUnitDescription(prPRead->uUnit)[0])
                    strcat(sUnitDescription(prPRead->uUnit),";SG:");
                else
                    strcpy(sUnitDescription(prPRead->uUnit),"SG:");
                strcat(sUnitDescription(prPRead->uUnit),p.pdb.cryst1.space_grp);
                SPACEGROUPt sg;
                printf("CRYST1: sg=\"%s\"\n",p.pdb.cryst1.space_grp);

                if (!parse_spacegroup_file(-1,p.pdb.cryst1.space_grp, &sg)) {
                    VP0("Parsed spacegroup info for \"%s\"\n",p.pdb.cryst1.space_grp);
                    VP0("Number of symmetry ops = %d\n",sg.n_symops);
                    if (GDefaults.bPdbExpandSymm) {
                        VarArraySetSize( prPRead->vaMatrices, sg.n_symops );
                        MATRIX  M, Mi, mFrac, mTmp;
                        PDBMATRIXt  *maPSymops = PVAI(prPRead->vaMatrices,PDBMATRIXt,0);
                        BuildFractionalTransforms( prPRead->uUnit, M, Mi );
                        for (int i = 0; i < sg.n_symops; i++ ) {
                            /* --- build fractional symop as 4x4 ---
                             * rotation in upper-left 3x3, translation in column 3 */
                            for (int j = 0; j < 4; j++ )
                                mFrac[j][0] = mFrac[j][1] = mFrac[j][2] = mFrac[j][3] = 0.0;
                            mFrac[3][3] = 1.0;

                            for (int r = 0; r < 3; r++ )
                                for (int c = 0; c < 3; c++ )
                                    mFrac[c][r] = (double)sg.symops[i].rot[r][c];  /* col-major */

                            mFrac[3][0] = (double)sg.symops[i].trans[0] / 12.0;
                            mFrac[3][1] = (double)sg.symops[i].trans[1] / 12.0;
                            mFrac[3][2] = (double)sg.symops[i].trans[2] / 12.0;

                            /* Mcart = M * Mfrac * Mi */
                            MatrixMultiply( mTmp, M, mFrac );
                            MatrixMultiply( maPSymops[i].mTransform, mTmp, Mi);
                            maPSymops[i].bUsed = (i != 0);
                        }
                    }
                } else VP0("Failed to find spacegroup!\n");

                break;

            case PDB_REMARK:
                // Currently only REMARK 350 is parsed for BIOMT
                if (!GDefaults.iPdbReadBioMT || p.pdb.remark.num != 350) break;
                if (!strncmp(p.pdb.remark.text,"BIOMOLECULE:",12)) {
                    iCurrentBioMT = atoi(p.pdb.remark.text+12);
                    VP0("Found Biomolecule definition # %d.\n",iCurrentBioMT);
                    if (iCurrentBioMT != 1) {
                         VPWARN("ReadBioMT only works correctly with a single biomolecule definition\n");
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
                if (iRow==1) VP2("Read PDB_REMARK 350 BIOMT1 record %d\n", iSerial );
                if(n<6 || iRow<1 || iRow>3 || iSerial<=0) break; //error
                p.pdb.mtrix.serial_num = iSerial;
                p.pdb.mtrix.row_num = iRow;
                p.pdb.mtrix.m1 = m1;
                p.pdb.mtrix.m2 = m2;
                p.pdb.mtrix.m3 = m3;
                p.pdb.mtrix.v = v;
                // First unit is identity, does not need replication
                p.pdb.mtrix.given = (iSerial == 1);
                // Now drop through to PDB_MTRIX processing of our BIOMT copied matix
                /* fall through */

            case PDB_MTRIX:
                VPTRACE("Read PDB_MTRIX record.\n" );
                iStart = iVarArrayElementCount(prPRead->vaMatrices);
                iSerial = p.pdb.mtrix.serial_num;
                if ( iSerial>iStart ) {
                    VarArraySetSize( (prPRead->vaMatrices), iSerial );
                    for (int i=iStart; i<iSerial; i++ ) {
                        PVAI(prPRead->vaMatrices,PDBMATRIXt,i)->bUsed = FALSE;
                        MatrixDiagonal(
                            PVAI(prPRead->vaMatrices,PDBMATRIXt,i)->mTransform,
                            0.0, 0.0, 0.0 );
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

            case PDB_MODEL:
                iCurrentModel = p.pdb.model.num;
                break;

            case PDB_LINK:
                if (GDefaults.bPdbUseLinkRecords)
                    VarArrayAdd( prPRead->vaLinkRecs, (GENP)&p.pdb.link );
                break;

            case PDB_SSBOND:
                if (GDefaults.bPdbUseLinkRecords) {
                    strcpy(pdbLink.name[0],"SG");
                    strcpy(pdbLink.name[1],"SG");
                    pdbLink.alt_loc[0] = ' ';
                    pdbLink.alt_loc[1] = ' ';
                    pdbLink.residues[0] = p.pdb.link.residues[0];
                    pdbLink.residues[1] = p.pdb.link.residues[1];
                    strcpy(pdbLink.symop[0], p.pdb.link.symop[0]);
                    strcpy(pdbLink.symop[1], p.pdb.link.symop[1]);
                    pdbLink.distance = p.pdb.link.distance;
                    VarArrayAdd( prPRead->vaLinkRecs, (GENP)&pdbLink );
                }
                break;

            default:
                break;
        }
    } while ( !feof(prPRead->fPdbFile) );

                /* Make the last residue a LASTEND */

    if ( (iLast = iVarArrayElementCount( prPRead->vaResidues )) ) {
        iLast--;        /* could be 0th element */
        PVAI( (prPRead->vaResidues), RESIDUENAMEt,iLast)->iTerminator = LASTEND;
        prPRead->iMaxSerialNum = iSerialNumMax+1;
    }

    if ( iMultipleResName ) {
        // another place conflation of PDB ResId and residue Sequence #
        VPNOTE("%d residues had naming warnings.\n"
            "Thus, there are split residues;\n"
            "residue sequence numbers will not correspond to those in the pdb.\n",
            iMultipleResName );
    }

    VPTRACEEXIT(__func__);
}

extern OBJEKT oCmd_loadAmberParams( int iArgCount, ASSOC aaArgs[] );
extern OBJEKT oCmd_loadAmberPrep( int iArgCount, ASSOC aaArgs[] );
extern OBJEKT oCmd_loadMol2( int iArgCount, ASSOC aaArgs[] );
extern OBJEKT oCmd_loadOff( int iArgCount, ASSOC aaArgs[] );

static BOOL
zbFileReadable(char *sFilename)
{
    FILE *fp = FOPENNOCOMPLAIN(sFilename, "r");
    if (fp) { fclose(fp); return TRUE; }
    return FALSE;
}

static void
zDetectResidueType(UNIT uTemplate) {
    RESIDUE rRes = (RESIDUE)oContainerFirstObject(uTemplate);
    if (cResidueType(rRes) != RESTYPEUNDEFINED) return;
    ATOM aHead = aUnitHead(uTemplate);
    ATOM aTail = aUnitTail(uTemplate);
    if (!aHead && !aTail) return;
    char *cPHead = aHead ? sContainerName(aHead) : NULL;
    char *cPTail = aTail ? sContainerName(aTail) : NULL;
    char cResType = RESTYPEUNDEFINED;
    if ( (!cPHead || !strcmp(cPHead,"N")) && (!cPTail || !strcmp(cPTail,"C"))
            && cContainerFindName((CONTAINER)rRes,ATOMid,"CA")
            && cContainerFindName((CONTAINER)rRes,ATOMid,"O")) {
        cResType = RESTYPEPROTEIN;
    } else if ( (!cPHead || !strcmp(cPHead,"P")) && (!cPTail || !strcmp(cPTail,"O3'"))
            && cContainerFindName((CONTAINER)rRes,ATOMid,"O5'")
            && cContainerFindName((CONTAINER)rRes,ATOMid,"C3'")) {
        cResType = RESTYPENUCLEIC;
    } else if ( (!cPHead || !strcmp(cPHead,"C1"))
            && (!cPTail || cPTail[0]=='O' || !strcmp(cPTail,"C1"))
            && cContainerFindName((CONTAINER)rRes,ATOMid,"C2")
            && cContainerFindName((CONTAINER)rRes,ATOMid,"C3")
            && cContainerFindName((CONTAINER)rRes,ATOMid,"C4")) {
        cResType = RESTYPESACCHARIDE;
    }
    VP0("ResidueType for Template RESIDUE %s set to %s\n",
           sContainerName(uTemplate),sResidueTypeNameFromChar(cResType));
    ResidueSetType(rRes,cResType);
}

// Attempt to load a residue/ligand parameter definition.
// Note that we cannot load a leaprc/command file, because the parser
// state is waiting for use to return and is not re-entrant.
static UNIT
zLoadUnit(char *name) {
    STRING sFilename;
#define NFORMATS 4
    char *format_ext[NFORMATS] = {"lib","prepc","prepi","mol2"};
    int i;
    if (!GDefaults.bPdbAutoLoadRes) return NULL;
    for (i=0;i<NFORMATS;i++) {
        sprintf(sFilename,"%s.%s",name,format_ext[i]);
        VP2("Search for residue parameter file: %s\n", sFilename);
        if (zbFileReadable(sFilename)) break;
        sprintf(sFilename,"%c/%s.%s",cLower(name[0]),name,format_ext[i]);
        VP2("Search for residue parameter file: %s\n", sFilename);
        if (zbFileReadable(sFilename)) break;
    }
    if (i==NFORMATS) {
        VP2(" ... not found.\n", sFilename);
        return NULL;
    }

    ASSOC aAssoc = (ASSOC)oCreate(ASSOCid);
    AssocSetName( aAssoc, "" );
    OBJEKT oFilename = oCreate(OSTRINGid);
    OStringDefine( (OSTRING) oFilename, sFilename );
    AssocSetObject( aAssoc, oFilename );
    DEREF( oFilename );	/* keeps count = 1 */
    UNIT uTemplate;
    RESIDUE rRes;

    // loadOff and loadAmberPrep define variables for all UNITs read in
    // We must look up uTemplate for these. For Mol2, we already have uTemplate
    // and must set the matching variable.

    // HEAD, TAIL, ResidueType defined in the library
    if (i==0) {
        oCmd_loadOff( 1, &aAssoc );
        uTemplate = (UNIT)oVariable(name);

    // PREPC and PREPI automatically get HEAD = first main atom, TAIL = last main atom
    // and these can be the same atom (happens some terminal residues)
    } else if (i==1 || i==2 ) {
        oCmd_loadAmberPrep(1, &aAssoc );
        uTemplate = (UNIT)oVariable(name);

    } else { //  Mol2
        // Mol2 has no default head/tail.
        uTemplate = (UNIT)oCmd_loadMol2(1, &aAssoc );
        if (uTemplate) VariableSet( name, (OBJEKT)uTemplate ); /* adds 1 REF */
    }

    if (!(uTemplate &&
            iObjectType(uTemplate) == UNITid &&
            iContainerNumberOfChildren(uTemplate) == 1 &&
            (rRes = (RESIDUE)oContainerFirstObject(uTemplate)) &&
            iObjectType(rRes) == RESIDUEid )) {
        VPFATAL("Loaded object is not a Template Residue\n");
        return NULL;
    }
    if (i==1 || i==2 ) zDetectResidueType(uTemplate);
    else if (i==3) ResidueSetType(rRes,RESTYPELIGAND);

    sprintf(sFilename,"%s.frcmod",name);
    VP2("Search for FRCMOD file: %s\n", sFilename);
    if ( zbFileReadable(sFilename) ) {
        OStringDefine( (OSTRING) oFilename, sFilename );
        oCmd_loadAmberParams( 1, &aAssoc );
        VP2("Loading FRCMOD file: %s\n", GsBasicsFullName);
    } else VP2(" ... not found.\n", sFilename);
    Destroy((OBJEKT*)&aAssoc); // FIXME: does this destroy the contained string?
    return uTemplate;
}

// Candidate matches.
// Standard mechanism = only match by name, nothing else matters
// But if we use "ResMap { [01] A B }" converts to "UNIT=A, RES=B, HEAD/TAIL definitions"
// Then we map resName = B AND HEAD/TAIL = T,F or F,T and this should be only one match
// But, in the case of not FIRST or LAST, it can mean both HEAD and TAIL (residue) or neither (ligand)
// These can be differentiated by ResidueType, except that property is not consistently set!
// Generally, UNIT head tail are both set, even in non-polymer.
// polymer residue middle: will match name and head, tail
// non-polymer molecule: will match name and generally match HEAD+TAIL and get linked even though it should not
//
// here both map to LYS, HEAD/TAIL decides:
// UNIT NLYS RESIDUE LYS
// UNIT LYS RESIDUE LYS
//
// if we have ResMap { A B } in existing system.
// and both residue A and B already exist:
// All A get converted to B
// in my system we get:
//
//  UNIT A, RESIDUE B  <- renamed
//  UNIT A, RESIDUE A  <- this one must be hidden
//  UNIT B, RESIDUE B
//
// ResMap does not modify residue, just add lookup key by renamed value
// If we match, and key != RESIDUE name, this gets priority over head/tail match
//
// ResMap { 0 LYS NLYS }
// UNIT NLYS RESIDUE NLYS index=LYS   NLYS->NLYS  LYS=>NLYS if head/tail match
//
// or
//
// #no ResMap
// UNIT NLYS RESIDUE LYS index=NLYS+LYS  NLYS->LYS,  LYS->LYS
//
// if both, we do not double index, we note that NLYS->LYS is already in the system and skip that ResMap
// ResMap { 0 LYS NLYS }
// UNIT NLYS RESIDUE LYS index=NLYS+LYS -- NLYS->LYS as it is now



//START



#define MAXATOMS       512   // max atoms in a single residue (max in PDB = 439)
#define MAXCANDIDATES   40   // sanity cap
#define MAX_UNMATCHED_CONNECT 2  // spurious contact limit per residue

// -- Scoring (least tolerable -> most tolerable) ------------------------------
typedef struct {
    int mismatched_head_restype;
    int extra_connect;
    int missing_head_tail;
    int extra_head_tail;
    int missing_atoms;
    int missing_connect;
    int extra_atoms;   // raw count; patch logic checks == 1 for deprotonation
    int pdb_seq;
} MatchScore;

char *sPMatchResNames[MAXCANDIDATES];

typedef struct MatchCandidate {
    UNIT        uTemplate;
    const Pair *pairPMatched[MAXCONNECT];   // detected bonds matching CONNECT
    const Pair *pairPUnmatched[MAXCONNECT]; // detected bonds not matching CONNECT
    BOOL        atomMatched[MAXATOMS];      // TRUE if this atom matched the template
    char        sUnmatchedNames[80];        // List of PDB names not matched
    BOOL        bHeadIsCrossLink;           // Flag: HEAD bond to non-TAIL
    BOOL        bTailIsCrossLink;           // Flag: TAIL bonds to non-HEAD
    BOOL        bHasPendingForwardLinks;    // Flag: CONNECT links to forward atoms
                                            // (CONNECT status of forward atom is unknown)
    MatchScore  score;
} MatchCandidate;

// -- Comparator --------------------------------------------------------------

static int
ziCompareCandidates(const MatchCandidate *A, const MatchCandidate *B)
{
    const MatchScore *a = &A->score;
    const MatchScore *b = &B->score;
    #define CMP(field) if (a->field < b->field) return -1; if (a->field > b->field) return 1;
    CMP(mismatched_head_restype)
    CMP(extra_connect)
    CMP(missing_head_tail)
    CMP(extra_head_tail)
    CMP(missing_atoms)
    CMP(missing_connect)
    CMP(extra_atoms)
    CMP(pdb_seq)
    return 0;
    #undef CMP
}

// -- Running-best ------------------------------------------------------------

static void
zLogCandidateLoss(const MatchCandidate *winner, const MatchCandidate *loser)
{
    if (GiVerbosityLevel < 3) return;
    RESIDUE rW = (RESIDUE)oContainerFirstObject(winner->uTemplate);
    RESIDUE rL = (RESIDUE)oContainerFirstObject(loser->uTemplate);
    #define EXPLAIN(field) \
        if (winner->score.field < loser->score.field) { \
            VP2("  CANDIDATE DROP: %s.%s loses to %s.%s on " #field " (%d > %d)\n", \
                   sContainerName(loser->uTemplate), sContainerName(rL), \
                   sContainerName(winner->uTemplate), sContainerName(rW), \
                   loser->score.field, winner->score.field); \
            return; \
        }
    EXPLAIN(mismatched_head_restype)
    EXPLAIN(extra_connect)
    EXPLAIN(missing_head_tail)
    EXPLAIN(extra_head_tail)
    if (winner->score.missing_atoms < loser->score.missing_atoms) {
        VP2("  CANDIDATE DROP: %s.%s loses to %s.%s on missing_atoms (%d > %d)",
               sContainerName(loser->uTemplate), sContainerName(rL),
               sContainerName(winner->uTemplate), sContainerName(rW),
               loser->score.missing_atoms, winner->score.missing_atoms);
        if (loser->sUnmatchedNames[0])
            VP2(" missing: %s", loser->sUnmatchedNames);
        VP2("\n");
        return;
    }
    EXPLAIN(missing_connect)
    EXPLAIN(extra_atoms)
    EXPLAIN(pdb_seq)
    VP2("  CANDIDATE DROP: %s.%s loses to %s.%s (equal score, first wins)\n",
           sContainerName(loser->uTemplate), sContainerName(rL),
           sContainerName(winner->uTemplate), sContainerName(rW));
    #undef EXPLAIN
}

static void
zUpdateBestCandidate(MatchCandidate *best, const MatchCandidate *challenger)
{
    if (!best->uTemplate) { *best = *challenger; return; }
    if (ziCompareCandidates(challenger, best) < 0) {
        zLogCandidateLoss(challenger, best);
        *best = *challenger;
    } else {
        zLogCandidateLoss(best, challenger);
    }
}

// -- Helpers -----------------------------------------------------------------

static int
zCountTemplateConnects(ATOM *aPConnect)
{
    int i;
    for (i=2; i < MAXCONNECT; i++) { if (!aPConnect[i]) break; }
    return i-2;
}

static int qsort_strcmp(const void *a, const void *b)
    { return strcmp(*(const char *const *)a, *(const char *const *)b); }

static int
zBuildLinkPatchName(char *key_out, size_t key_max, BOOL bSort,
                    const char *res1, const char *atom1,
                    const char *res2, const char *atom2)
{
    char left[256], right[256];
    snprintf(left,  sizeof(left),  "%s@%s", res1, atom1);
    snprintf(right, sizeof(right), "%s@%s", res2, atom2);
    if (bSort && strcmp(left, right) > 0) {
        snprintf(key_out, key_max, "%s_link_%s", right, left); return 1;
    }
    snprintf(key_out, key_max, "%s_link_%s", left, right);
    return 0;
}

// -- Patch dispatch ----------------------------------------------------------

static BOOL
zProcessPatches(PDBREADt *prPPdb, STRING sPatchName[], int iNumPatchNames)
{
    STRING sFilename;
    int i;
    for (i = 0; i < iNumPatchNames; i++) {
        sprintf(sFilename, "%s.patch", sPatchName[i]);
        VP2("Search for patch file: %s\n", sFilename);
        if (zbFileReadable(sFilename)) break;
    }
    if (i == iNumPatchNames) { VP2(" ... not found\n"); return FALSE; }
    VP2("Found patch file: %s\n", GsBasicsFullName);
    VP2("# Patch %s\n",   sPatchName[i]);
    VP2("RESIDUE = %d\n", prPPdb->iSubstRes1);
    VP2("source %s\n",    GsBasicsFullName);
    if (prPPdb->fpPatchFileOut) {
        fprintf(prPPdb->fpPatchFileOut, "# Patch %s\n",   sPatchName[i]);
        fprintf(prPPdb->fpPatchFileOut, "RESIDUE = %d\n", prPPdb->iSubstRes1);
        fprintf(prPPdb->fpPatchFileOut, "source %s\n",    GsBasicsFullName);
    }
    return TRUE;
}

static BOOL
zProcessPatches2(PDBREADt *prPPdb,
                 STRING sPatchName[], BOOL bPatchReversed[], int iNumPatchNames)
{
    STRING sFilename;
    int i;
    for (i = 0; i < iNumPatchNames; i++) {
        sprintf(sFilename, "%s.patch", sPatchName[i]);
        VP2("Search for patch file: %s\n", sFilename);
        if (zbFileReadable(sFilename)) break;
    }
    if (i == iNumPatchNames) { VP2(" ... not found\n"); return FALSE; }
    VP2("Found patch file: %s\n", GsBasicsFullName);
    VP2("# Patch %s\n",    sPatchName[i]);
    VP2("RESIDUE1 = %d\n", bPatchReversed[i] ? prPPdb->iSubstRes2 : prPPdb->iSubstRes1);
    VP2("ATOM1 =    %s\n", bPatchReversed[i] ? prPPdb->sSubstAtom2 : prPPdb->sSubstAtom1);
    VP2("RESIDUE2 = %d\n", bPatchReversed[i] ? prPPdb->iSubstRes1 : prPPdb->iSubstRes2);
    VP2("ATOM2 =    %s\n", bPatchReversed[i] ? prPPdb->sSubstAtom1 : prPPdb->sSubstAtom2);
    VP2("source %s\n",     GsBasicsFullName);
    if (prPPdb->fpPatchFileOut) {
        fprintf(prPPdb->fpPatchFileOut, "# Patch %s\n", sPatchName[i]);
        if (bPatchReversed[i]) {
            fprintf(prPPdb->fpPatchFileOut, "RESIDUE1 = %d\n", prPPdb->iSubstRes2);
            fprintf(prPPdb->fpPatchFileOut, "ATOM1 =    %s\n", prPPdb->sSubstAtom2);
            fprintf(prPPdb->fpPatchFileOut, "RESIDUE2 = %d\n", prPPdb->iSubstRes1);
            fprintf(prPPdb->fpPatchFileOut, "ATOM2 =    %s\n", prPPdb->sSubstAtom1);
        } else {
            fprintf(prPPdb->fpPatchFileOut, "RESIDUE1 = %d\n", prPPdb->iSubstRes1);
            fprintf(prPPdb->fpPatchFileOut, "ATOM1 =    %s\n", prPPdb->sSubstAtom1);
            fprintf(prPPdb->fpPatchFileOut, "RESIDUE2 = %d\n", prPPdb->iSubstRes2);
            fprintf(prPPdb->fpPatchFileOut, "ATOM2 =    %s\n", prPPdb->sSubstAtom2);
        }
        fprintf(prPPdb->fpPatchFileOut, "source %s\n", GsBasicsFullName);
    }
    return TRUE;
}

// -- Bond emission + link patch search ---------------------------------------

static void
zMakeBond(PDBREADt *prPPdb, const Pair *p, BOOL bFromIsConnect,
          const char *cPResName, const char *cPResType, int iResGlobal)
{
    RESIDUENAMEt *rnPResidues  = PVAI(prPPdb->vaResidues, RESIDUENAMEt, 0);
    RESIDUENAMEt *rnTo         = &rnPResidues[p->to_group];
    ATOMNAMEt    *anFrom       = PVAI(prPPdb->vaAtomRecs, ATOMNAMEt, p->from_member);
    ATOMNAMEt    *anTo         = PVAI(prPPdb->vaAtomRecs, ATOMNAMEt, p->to_member);
    int           iNumPatchNames = 0;

    BOOL bToIsConnect = FALSE;
    if (rnTo->rResidue) {
        ATOM *aPToConnect = (ATOM*)rnTo->rResidue->aaConnect;
        for (int jc = 2; jc < MAXCONNECT && aPToConnect[jc]; jc++) {
            if (!strcmp(anTo->sName, sContainerName(aPToConnect[jc])))
                { bToIsConnect = TRUE; break; }
        }
        if (!bToIsConnect && aPToConnect[0] && !strcmp(anTo->sName, sContainerName(aPToConnect[0]))) bToIsConnect = TRUE;
        if (!bToIsConnect && aPToConnect[1] && !strcmp(anTo->sName, sContainerName(aPToConnect[0]))) bToIsConnect = TRUE;
    }

    ATOM    aFrom = (ATOM)cContainerFindName((CONTAINER)prPPdb->rRes,   ATOMid, anFrom->sName);
    ATOM    aTo   = (ATOM)cContainerFindName((CONTAINER)rnTo->rResidue, ATOMid, anTo->sName);
    RESIDUE rFrom = prPPdb->rRes;
    RESIDUE rTo   = rnTo->rResidue;

    int iNonConnectSides = !bFromIsConnect + !bToIsConnect;
    if ((iNonConnectSides >= 2 && GDefaults.iPdbIgnoreNonConnect >= 1) ||
        (iNonConnectSides >= 1 && GDefaults.iPdbIgnoreNonConnect >= 2)) {
        if (bToIsConnect)
            VPWARN("Filtered out bond from %s.%s:%d.%s to %s.%s:%d.%s\n"
                   "   with pdb_ignore_nonconnect = %d\n"
                   "   but earlier residue (from) match was based on this pending bond.\n",
                   sContainerName(rTo),   sResidueChainId(rTo),   iResiduePdbSequence(rTo),   sContainerName(aTo),
                   sContainerName(rFrom), sResidueChainId(rFrom), iResiduePdbSequence(rFrom), sContainerName(aFrom),
                   GDefaults.iPdbIgnoreNonConnect);
        return;
    }

    if (aFrom && aTo) {
        VP0("%s bond from %s.%s:%d.%s to %s.%s:%d.%s\n",
            GDefaults.bPdbAutoLink ? "Created" : "Detected",
            sContainerName(rFrom), sResidueChainId(rFrom), iResiduePdbSequence(rFrom), sContainerName(aFrom),
            sContainerName(rTo),   sResidueChainId(rTo),   iResiduePdbSequence(rTo),   sContainerName(aTo));
    } else {
        VP0("Warning: could not resolve atoms for bond %d.%s -> %d.%s\n",
            p->from_group, anFrom->sName, p->to_group, anTo->sName);
    }

    /* Both connect: patch optional -- autolink makes the bond if patch didn't.
     * Any non-connect: patch mandatory -- patch is responsible for bonding. */
    BOOL bPatchMandatory = !(bFromIsConnect && bToIsConnect);

    STRING sPatchName[4];
    BOOL   bPatchReversed[4];
    char   cResType2  = cResidueType(rnTo->rResidue);
    char  *cPRes2Type = strchr(RESTYPECLASSPOLYMER, cResType2)
                        ? sResidueTypeNameFromChar(cResType2) : NULL;

    prPPdb->iSubstRes1 = iResGlobal;
    strncpy(prPPdb->sSubstAtom1, anFrom->sName, sizeof(prPPdb->sSubstAtom1) - 1);
    prPPdb->iSubstRes2 = iContainerSequence(rnTo->rResidue);
    strncpy(prPPdb->sSubstAtom2, anTo->sName,   sizeof(prPPdb->sSubstAtom2) - 1);

    /* name/name sorted -- most specific */
    bPatchReversed[iNumPatchNames] = zBuildLinkPatchName(sPatchName[iNumPatchNames],
                            sizeof(STRING), TRUE, cPResName, anFrom->sName, rnTo->sName, anTo->sName);
    iNumPatchNames++;
    if (cPResType) {
        bPatchReversed[iNumPatchNames] = 0;
        zBuildLinkPatchName(sPatchName[iNumPatchNames], sizeof(STRING), FALSE,
                            cPResType, anFrom->sName, rnTo->sName, anTo->sName);
        iNumPatchNames++;
    }
    if (cPRes2Type) {
        bPatchReversed[iNumPatchNames] = 1;
        zBuildLinkPatchName(sPatchName[iNumPatchNames], sizeof(STRING), FALSE,
                            cPRes2Type, anTo->sName, cPResName, anFrom->sName);
        iNumPatchNames++;
    }
    if (cPResType && cPRes2Type) {
        bPatchReversed[iNumPatchNames] = zBuildLinkPatchName(sPatchName[iNumPatchNames],
                            sizeof(STRING), TRUE, cPResType, anFrom->sName, cPRes2Type, anTo->sName);
        iNumPatchNames++;
    }

    BOOL bPatchFound = zProcessPatches2(prPPdb, sPatchName, bPatchReversed, iNumPatchNames);

    if (!bPatchMandatory) {
        /* Optional (connect->connect): autolink makes the bond if patch didn't */
        if (aFrom && aTo && GDefaults.bPdbAutoLink && !bAtomBondedTo(aFrom, aTo))
            AtomBondTo(aFrom, aTo);
    } else if (!bPatchFound) {
        /* Mandatory (non-connect) and no patch: error, do not bond */
        VPERROR("Mandatory patch file was not found for non-CONNECT bond formation"
                " from %s.%s:%d.%s to %s.%s:%d.%s\n",
                sContainerName(rFrom), sResidueChainId(rFrom), iResiduePdbSequence(rFrom), sContainerName(aFrom),
                sContainerName(rTo),   sResidueChainId(rTo),   iResiduePdbSequence(rTo),   sContainerName(aTo));
    }
}

// -- Post-match modification / patch generation ------------------------------

static void
zProcessModifications(PDBREADt *prPPdb, MatchCandidate *match,
                      int iFirstAtom, int iLastAtom)
{
    if (!match->uTemplate) return;

    RESIDUE rResidue   = (RESIDUE)oContainerFirstObject(match->uTemplate);
    char   *cPResName  = sContainerName(rResidue);
    char    cResType   = cResidueType(rResidue);
    char   *cPResType  = strchr(RESTYPECLASSPOLYMER, cResType)
                         ? sResidueTypeNameFromChar(cResType) : NULL;
    int     iResGlobal = iContainerSequence(prPPdb->rRes);

    if (match->pairPMatched[0] && match->bHeadIsCrossLink &&
                match->pairPMatched[0]->to_group < match->pairPMatched[0]->from_group)
        zMakeBond(prPPdb, match->pairPMatched[0], TRUE, cPResName, cPResType, iResGlobal);
    if (match->pairPMatched[1] && match->bTailIsCrossLink &&
                match->pairPMatched[1]->to_group < match->pairPMatched[1]->from_group)
        zMakeBond(prPPdb, match->pairPMatched[1], TRUE, cPResName, cPResType, iResGlobal);
    for (int i = 2; i < MAXCONNECT; i++) {
        if (match->pairPMatched[i] &&
                match->pairPMatched[i]->to_group < match->pairPMatched[i]->from_group)
            zMakeBond(prPPdb, match->pairPMatched[i], TRUE, cPResName, cPResType, iResGlobal);
    }

    if (!(match->score.mismatched_head_restype ||
          match->score.extra_head_tail         ||
          match->score.missing_head_tail        ||
          match->score.extra_connect            ||
          match->score.missing_atoms            ||
          match->score.missing_connect          ||
          match->score.extra_atoms == 1)) return;

    if (match->bHasPendingForwardLinks)
        VP0("PATCH MAY BE REQUIRED (unresolved forward links pending):\n");
    else
        VP0("PATCH REQUIRED:\n");
    VP0(" UNIT=%s, RESIDUE=%s\n", sContainerName(match->uTemplate), cPResName);

    if (match->score.extra_connect)
        VPFATAL("Template residue with extra connects should never happen. A residue version\n"
                "without CONNECT atoms is expected to exist. Check your residue library.\n");

    if (match->score.mismatched_head_restype) {
        match->pairPUnmatched[match->score.missing_connect++] = match->pairPMatched[0];
        match->pairPMatched[0] = NULL;
    }

    BOOL bStructuralPatchEmitted = FALSE;
    STRING sPatchName[10];

    if (match->score.missing_head_tail) {
        if (match->pairPMatched[0] && !aUnitHead(match->uTemplate)) {
            ATOMNAMEt *an = PVAI(prPPdb->vaAtomRecs, ATOMNAMEt,
                                 match->pairPMatched[0]->from_member);
            prPPdb->iSubstRes1 = iResGlobal;
            sprintf(sPatchName[0], "%s@%s_head", cPResName, an->sName);
            if (cPResType) sprintf(sPatchName[1], "%s@%s_head", cPResType, an->sName);
            if (zProcessPatches(prPPdb, sPatchName, cPResType ? 2 : 1)) bStructuralPatchEmitted = TRUE;
        }
        if (match->pairPMatched[1] && !aUnitTail(match->uTemplate)) {
            ATOMNAMEt *an = PVAI(prPPdb->vaAtomRecs, ATOMNAMEt,
                                 match->pairPMatched[1]->from_member);
            prPPdb->iSubstRes1 = iResGlobal;
            sprintf(sPatchName[0], "%s@%s_tail", cPResName, an->sName);
            if (cPResType) sprintf(sPatchName[1], "%s@%s_tail", cPResType, an->sName);
            if (zProcessPatches(prPPdb, sPatchName, cPResType ? 2 : 1)) bStructuralPatchEmitted = TRUE;
        }
    }

    if (match->score.extra_head_tail) {
        ATOM aHead = aUnitHead(match->uTemplate);
        ATOM aTail = aUnitTail(match->uTemplate);
        if (!match->pairPMatched[0] && aHead) {
            prPPdb->iSubstRes1 = iResGlobal;
            sprintf(sPatchName[0], "%s@%s_cap", cPResName, sContainerName(aHead));
            if (cPResType) sprintf(sPatchName[1], "%s@%s_cap", cPResType, sContainerName(aHead));
            if (zProcessPatches(prPPdb, sPatchName, cPResType ? 2 : 1)) bStructuralPatchEmitted = TRUE;
        }
        if (!match->pairPMatched[1] && aTail) {
            prPPdb->iSubstRes1 = iResGlobal;
            sprintf(sPatchName[0], "%s@%s_cap", cPResName, sContainerName(aTail));
            if (cPResType) sprintf(sPatchName[1], "%s@%s_cap", cPResType, sContainerName(aTail));
            if (zProcessPatches(prPPdb, sPatchName, cPResType ? 2 : 1)) bStructuralPatchEmitted = TRUE;
        }
    }

    ATOM aExtraAtom = NULL;
    if (match->score.extra_atoms == 1) {
        LOOP lAtoms = lLoop((OBJEKT)rResidue, ATOMS);
        ATOM aCand;
        while ((aCand = (ATOM)oNext(&lAtoms)) != NULL) {
            BOOL found = FALSE;
            for (int i = iFirstAtom; i <= iLastAtom; i++) {
                if (!strcmp(sContainerName(aCand),
                            PVAI(prPPdb->vaAtomRecs, ATOMNAMEt, i)->sName))
                    { found = TRUE; break; }
            }
            if (!found) { aExtraAtom = aCand; break; }
        }
        if (!aExtraAtom)
            VPFATAL("Programming error looking for extra template atom in %s, %s:%d\n",
                     __func__, __FILE__, __LINE__);
    }

    if (match->score.extra_atoms == 1 && aExtraAtom && !bStructuralPatchEmitted) {
        sprintf(sPatchName[0], "%s_del_%s", cPResName, sContainerName(aExtraAtom));
        if (cPResType) sprintf(sPatchName[1], "%s_del_%s", cPResType, sContainerName(aExtraAtom));
        prPPdb->iSubstRes1 = iResGlobal;
        zProcessPatches(prPPdb, sPatchName, cPResType ? 2 : 1);
    }

    if (match->score.missing_atoms && !bStructuralPatchEmitted) {
        char *names[MAXATOMS];
        int   nNames = 0;
        for (int i = iFirstAtom, n = 0; i <= iLastAtom; i++, n++) {
            if (!match->atomMatched[n])
                names[nNames++] = PVAI(prPPdb->vaAtomRecs, ATOMNAMEt, i)->sName;
        }
        qsort(names, nNames, sizeof(char *), qsort_strcmp);
        sprintf(sPatchName[0], "%s_add_%s", cPResName, names[0]);
        if (cPResType) sprintf(sPatchName[1], "%s_add_%s", cPResType, names[0]);
        for (int i = 1; i < nNames; i++) {
            strcat(sPatchName[0], ","); strcat(sPatchName[0], names[i]);
            if (cPResType) { strcat(sPatchName[1], ","); strcat(sPatchName[1], names[i]); }
        }
        prPPdb->iSubstRes1 = iResGlobal;
        zProcessPatches(prPPdb, sPatchName, cPResType ? 2 : 1);
    }

    if (match->score.missing_connect) {
        for (int i = 0; i < match->score.missing_connect && i < MAXCONNECT; i++) {
            if (!match->pairPUnmatched[i]) continue;
            if (match->pairPUnmatched[i]->from_member <
                    match->pairPUnmatched[i]->to_member) continue;
            zMakeBond(prPPdb, match->pairPUnmatched[i], FALSE,
                      cPResName, cPResType, iResGlobal);
        }
    }
}


static void
zMatchResidueCandidate(PDBREADt *prPPdb, RESIDUENAMEt *rnPResidues, int iResIndex,
        ATOMNAMEt *anPAtoms, int iFirstAtom, int iLastAtom, BOOL bPDBHasHydrogens, // PDB info
        const Pair *pairPDetected[], int iNumDetected,           // bonding info
        char *key, UNIT uTemplate,                         // template info
        MatchCandidate *cnPCandidate)                      // result
{
     RESIDUE rTemplateRes = (RESIDUE)oContainerFirstObject(uTemplate);
     ATOM aHead = aUnitHead(uTemplate);
     ATOM aTail = aUnitTail(uTemplate);
     ATOM *aPTemplateConnect      = (ATOM *)rTemplateRes->aaConnect;
     int   numTemplateConnections = zCountTemplateConnects(aPTemplateConnect);

     if (GiVerbosityLevel > 6) {
         VP2("key=%s, unitName=%s, resName=%s, typ=%c, iNumDetected=%d, "
             "head=%s, tail=%s, conn=[%s,%s,%s]\n",
             key, sContainerName(uTemplate),
             sContainerName(rTemplateRes), cResidueType(rTemplateRes), iNumDetected,
             aHead ? sContainerName(aHead) : "None",
             aTail ? sContainerName(aTail) : "None",
             aPTemplateConnect[2] ? sContainerName(aPTemplateConnect[2]) : "None",
             aPTemplateConnect[3] ? sContainerName(aPTemplateConnect[3]) : "None",
             aPTemplateConnect[4] ? sContainerName(aPTemplateConnect[4]) : "None");
     }

     memset(cnPCandidate, 0, sizeof(*cnPCandidate));
     cnPCandidate->uTemplate = uTemplate;
     int matched_conn = 0, unmatched_conn = 0;

     // -- CONNECT Atom matching ---------------------------------------------------
     for (int id = 0; id < iNumDetected; id++) {
         const Pair *p         = pairPDetected[id];
         const char *sFromAtom = anPAtoms[p->from_member].sName;
         BOOL isBackward    = (p->to_group <  iResIndex);
         BOOL isAdjacentBack= (p->to_group == iResIndex - 1);
         BOOL isAdjacentFwd = (p->to_group == iResIndex + 1);
         BOOL isForward     = (p->to_group >  iResIndex);

         /* -- HEAD -- */
         if (aHead && !strcmp(sFromAtom, sContainerName(aHead))) {
             if (cnPCandidate->pairPMatched[0]) {
                 BOOL bCurrentIsAdj = cnPCandidate->pairPMatched[0] &&
                                      cnPCandidate->pairPMatched[0]->to_group == iResIndex - 1;
                 if (bCurrentIsAdj || !isAdjacentBack) continue;
             }
             cnPCandidate->pairPMatched[0] = p;
             // additional checks for TAIL(prev) -> HEAD(curr) because TAIL(prev) is mapped
             BOOL partnerIsTail = isAdjacentBack && prPPdb->rRes && aUnitTail(prPPdb->uUnit) &&
                              !strcmp(anPAtoms[p->to_member].sName, sContainerName(aUnitTail(prPPdb->uUnit)));
             cnPCandidate->bHeadIsCrossLink = !partnerIsTail;
             if (!cnPCandidate->bHeadIsCrossLink &&
                     cResidueType(prPPdb->rRes) != RESTYPEUNDEFINED &&
                     cResidueType(rTemplateRes) != RESTYPEUNDEFINED &&
                     cResidueType(prPPdb->rRes) != cResidueType(rTemplateRes)) {
                 VP0("HEAD/TAIL ResidueType mismatch (%s=%s,%s=%s)"
                     "   marking TAIL-HEAD connection as crosslink\n",
                     sContainerName(prPPdb->rRes), sResidueTypeNameFromChar(cResidueType(prPPdb->rRes)),
                     sContainerName(uTemplate), sResidueTypeNameFromChar(cResidueType(rTemplateRes)));
                 cnPCandidate->score.mismatched_head_restype = 1;
                 cnPCandidate->bHeadIsCrossLink = TRUE;
             }
             continue;
         }
         /* -- TAIL -- */
         if (aTail && !strcmp(sFromAtom, sContainerName(aTail))) {
             if (cnPCandidate->pairPMatched[1]) {
                 BOOL bCurrentIsAdj = cnPCandidate->pairPMatched[1] &&
                                   cnPCandidate->pairPMatched[1]->to_group == iResIndex + 1;
                 if (bCurrentIsAdj || !isAdjacentFwd) continue;
             }
             cnPCandidate->pairPMatched[1]  = p;
             cnPCandidate->bTailIsCrossLink = isBackward || (isForward && !isAdjacentFwd);
             continue;
         }
         /* -- CONNECT2 ... CONNECT5 -- */
         BOOL match=FALSE;
         for (int jc = 2; jc < MAXCONNECT && aPTemplateConnect[jc]; jc++) {
             if (!strcmp(sFromAtom, sContainerName(aPTemplateConnect[jc]))) {
                 cnPCandidate->pairPMatched[jc] = p;
                 matched_conn++;
                 match=TRUE;
                 break;
             }
         }
         if (match) continue;

         // Analyze non-CONNECT atom bonding contacts
         if (GDefaults.iPdbIgnoreNonConnect >= 2) continue;
         BOOL bToIsConnect = FALSE;
         BOOL bForwardLink = (pairPDetected[id]->to_group > iResIndex);
         if (bForwardLink) {
             bToIsConnect = TRUE; /* assume connect -- verified when partner is processed */
         } else {
             RESIDUENAMEt *rnTo = &rnPResidues[pairPDetected[id]->to_group];
             if (rnTo->rResidue) {
                 ATOMNAMEt *anTo        = PVAI(prPPdb->vaAtomRecs, ATOMNAMEt,
                                               pairPDetected[id]->to_member);
                 ATOM      *aPToConnect = (ATOM *)rnTo->rResidue->aaConnect;
                 for (int jc = 0; jc < MAXCONNECT && aPToConnect[jc]; jc++) {
                     if (!strcmp(anTo->sName, sContainerName(aPToConnect[jc])))
                         { bToIsConnect = TRUE; break; }
                 }
             }
         }
         if (bToIsConnect && GDefaults.iPdbIgnoreNonConnect >= 1) continue;
         if (bForwardLink) cnPCandidate->bHasPendingForwardLinks = TRUE;
         if (unmatched_conn < MAXCONNECT)
             cnPCandidate->pairPUnmatched[unmatched_conn++] = p;
     }

     cnPCandidate->score.extra_head_tail =
         ((!cnPCandidate->pairPMatched[0] && aHead) ? 1 : 0) +
         ((!cnPCandidate->pairPMatched[1] && aTail) ? 1 : 0);
     cnPCandidate->score.missing_head_tail =
         ((cnPCandidate->pairPMatched[0] && !aHead) ? 1 : 0) +
         ((cnPCandidate->pairPMatched[1] && !aTail) ? 1 : 0);

     cnPCandidate->score.missing_connect = unmatched_conn;
     cnPCandidate->score.extra_connect   = numTemplateConnections - matched_conn;

     // -- Atom matching ---------------------------------------------------

     int  num_unmatched = 0; // Number of atoms in residue not found in template
     cnPCandidate->sUnmatchedNames[0] = 0; // String list

     for (int i = iFirstAtom, n = 0; i <= iLastAtom; i++, n++) {
         cnPCandidate->atomMatched[n] = TRUE;
         const char *name    = anPAtoms[i].sName;
         BOOL found = (cContainerFindName((CONTAINER)rTemplateRes, ATOMid, name) != NULL);
         if (!found) {
             // Try atom name aliases
             char *cPAtomName = zcPPdbMapName(SdAtomNameMap, NOEND, name, rTemplateRes);
             if (cPAtomName)
                 found = (cContainerFindName((CONTAINER)rTemplateRes, ATOMid, cPAtomName) != NULL);
             if (!found) {
                 // Atom not in template, mark and append to string list
                 cnPCandidate->atomMatched[n] = FALSE;
                 if (strlen(cnPCandidate->sUnmatchedNames) <= sizeof(cnPCandidate->sUnmatchedNames)-10) {
                     if (strlen(cnPCandidate->sUnmatchedNames)+strlen(name) <= sizeof(cnPCandidate->sUnmatchedNames)-10) {
                         if (num_unmatched) strcat(cnPCandidate->sUnmatchedNames, ",");
                         strcat(cnPCandidate->sUnmatchedNames, name);
                     } else strcat(cnPCandidate->sUnmatchedNames, "...");
                 }
                 num_unmatched++;
             }
         }
     }

     // number of atoms in template, excluding hydrogen if not bPDBHasHydrogens
     int  na_template_counted = iContainerNumberOfChildren(rTemplateRes);
     if (bPDBHasHydrogens) {
         ATOM aAtom;
         LOOP lAtoms = lLoop((OBJEKT)rTemplateRes, ATOMS);
         while ((aAtom = (ATOM)oNext(&lAtoms)) != NULL)
             if (aAtom->iAtomicNumber == HYDROGEN) na_template_counted--;
     }

     int na_residue = iLastAtom - iFirstAtom + 1;
     cnPCandidate->score.missing_atoms = num_unmatched;
     cnPCandidate->score.extra_atoms   = na_template_counted - (na_residue - num_unmatched);
     cnPCandidate->score.pdb_seq       = iResiduePdbSequence(oContainerFirstObject(uTemplate));
}

// -- Primary residue template matching function ------------------------------

static UNIT
zPdbMatchResidueTemplate(PDBREADt *prPPdb, int iResIndex, MatchCandidate *cnPMatchReturn)
{
    IX_REC          ix_record;
    UNIT            uTemplate, uResult;
    int             status;
    int             iNumCandidates = 0;
    MatchCandidate  cnMatch, cnCandidate;
    memset(&cnMatch, 0, sizeof(cnMatch));
    RESIDUENAMEt   *rnPResidues   = PVAI(prPPdb->vaResidues, RESIDUENAMEt, 0);
    ATOMNAMEt      *anPAtoms      = PVAI(prPPdb->vaAtomRecs, ATOMNAMEt, 0);
    RESIDUENAMEt   *rnPName       = &rnPResidues[iResIndex];
    const int       iResidueCount = iVarArrayElementCount(prPPdb->vaResidues);
    const int       iFirstAtom    = rnPName->iFirstAtom;
    const int       iLastAtom     = (iResIndex < iResidueCount - 1)
                                    ? ((rnPName+1)->iFirstAtom - 1)
                                    : (iVarArrayElementCount(prPPdb->vaAtomRecs) - 1);
    const Pair     *ptPPairs      = NULL;
    unsigned int    iNumPairs     = 0;
    const Pair     *pairPDetected[MAXCONNECT] = { NULL };
    int             iNumDetected  = 0;

    VPTRACEENTER(__func__);
    VP2("------------------------------\n");

    // -- Link detection ------------------------------------------------------

    if (neighbor_grid_query_group(prPPdb->ngBondGrid, iResIndex,
                                  &ptPPairs, &iNumPairs) != 0) {
        VPFATAL("neighbor_grid_query_group(%d) failed\n", iResIndex);
        return NULL;
    }
    int iCloseContacts = 0;
    for (unsigned int i = 0; i < iNumPairs; i++) {
        const Pair *p = &ptPPairs[i];
        if ( p->d2 < (MINBONDLEN*MINBONDLEN) ) {
             iCloseContacts++;
             continue;
        }
        if (!bAtomsBondedDist(anPAtoms[p->from_member].iElement,
                              anPAtoms[p->to_member].iElement, p->d2,
                              (abs(p->to_group - p->from_group) <= 1)
                                  ? GDefaults.dPdbLinkCovalentCutoff
                                  : GDefaults.dPdbCrosslinkCovalentCutoff))
            continue;
        if (iNumDetected >= MAXCONNECT) {
            VPFATAL("Too many detected bonding contacts in residue %d\n", iResIndex);
            break;
        }
        pairPDetected[iNumDetected++] = p;
    }
    if (iCloseContacts>0) {
        VPWARN("%s close contacts (< %f A) encountered\n", iCloseContacts, MINBONDLEN);
    }

    // -- Template lookup -----------------------------------------------------

    // check HetID to CompID mapping -- for input PDB without HETNAM/HETSYN COMP_ID mapping
    char *cPResName = zcPPdbMapName( SdResidueNameMap,
                        HETID_COMPID, rnPName->sName, NULL );
    if (cPResName) {
        VP0("HetID %s -> COMP_ID %s\n",rnPName->sName,cPResName);
        strcpy(ix_record.key, cPResName);
        ix_record.recptr = NULL;
        // First lookup always fails with data=NULL (IX_DUPKEY)
        status = locate_key(&ix_record, &(prPPdb->ixResMap), 0);
        status = next_key(&ix_record, &(prPPdb->ixResMap));
        if (strcmp(cPResName, ix_record.key)) {
            VP0("Not found: %s\n", cPResName);
            uTemplate = zLoadUnit(cPResName);
            // If template loaded, add record to the index
            if (uTemplate) {
                strcpy(ix_record.key, cPResName);
                ix_record.recptr = uTemplate;
                add_key(&ix_record, &(prPPdb->ixResMap));
            } else cPResName = NULL;
        } else cPResName = NULL;
    }
    // Now try direct residue name lookup
    if (!cPResName) {
        cPResName = rnPName->sName;
        strcpy(ix_record.key, cPResName);
        ix_record.recptr = NULL;
        status = locate_key(&ix_record, &(prPPdb->ixResMap), 0);
        status = next_key(&ix_record, &(prPPdb->ixResMap));
        if (strcmp(cPResName, ix_record.key)) {
            VP0("Not found: %s\n", cPResName);
            uTemplate = zLoadUnit(cPResName);
            if (uTemplate) {
                strcpy(ix_record.key, cPResName);
                ix_record.recptr = uTemplate;
                add_key(&ix_record, &(prPPdb->ixResMap));
            }
        }
    }

    // -- Hydrogen check ------------------------------------------------------

    BOOL bPDBHasHydrogens = FALSE;
    for (int i = iFirstAtom; i <= iLastAtom; i++) {
        if (anPAtoms[i].iElement == HYDROGEN) bPDBHasHydrogens = TRUE;
        break;
    }

    // -- Candidate scoring loop ----------------------------------------------

    if (status != IX_END && !strcmp(cPResName, ix_record.key)) do {
        uTemplate = (UNIT)ix_record.recptr;

        zMatchResidueCandidate(prPPdb, rnPResidues, iResIndex,
                anPAtoms, iFirstAtom, iLastAtom, bPDBHasHydrogens,
                pairPDetected, iNumDetected,
                ix_record.key, uTemplate,
                &cnCandidate);

        zUpdateBestCandidate(&cnMatch, &cnCandidate);

        sPMatchResNames[iNumCandidates] = sContainerName(uTemplate);
        iNumCandidates++;
        if (iNumCandidates >= MAXCANDIDATES) break;
        status = next_key(&ix_record, &(prPPdb->ixResMap));

    } while (status != IX_END && !strcmp(cPResName, ix_record.key));

    // -- Report and return ---------------------------------------------------

    if (iNumCandidates > 0) {
        if (GiVerbosityLevel > 2) {
            VP2("%s.%s:%d Candidates: %s", cPResName,
                rnPName->sChainId, rnPName->iPdbSequence, sPMatchResNames[0]);
            for (int i = 1; i < iNumCandidates; i++) VP2(",%s", sPMatchResNames[i]);
            VP2(" Match: %s%s\n", sContainerName(cnMatch.uTemplate),
                bPDBHasHydrogens ? "" : " [non-H match]");
        }

        uTemplate = cnMatch.uTemplate;
        RESIDUE rTemplateRes = (RESIDUE)oContainerFirstObject(uTemplate);

        if (cnMatch.score.extra_connect        ||
                cnMatch.score.mismatched_head_restype ||
                cnMatch.score.missing_head_tail       ||
                cnMatch.score.extra_head_tail         ||
                cnMatch.score.missing_atoms           ||
                cnMatch.score.missing_connect         ||
                cnMatch.score.extra_atoms) {
            VP0("Imprecise match for residue: %3s %2.2s%4d%c\n",
                cPResName, rnPName->sChainId,
                rnPName->iPdbSequence, rnPName->iCode);
            VP0(" Template: UNIT=%s, RESIDUE=%s\n",
                sContainerName(cnMatch.uTemplate), sContainerName(rTemplateRes));
            if (cnMatch.score.mismatched_head_restype)
                VP0(" Previous-Head link ResidueType mismatch: %s,%s\n",
                    sResidueTypeNameFromChar(cResidueType(prPPdb->rRes)),
                    sResidueTypeNameFromChar(cResidueType(rTemplateRes)));
            if (cnMatch.score.extra_connect)
                VP0(" Unused connect=%d\n", cnMatch.score.extra_connect);
            if (cnMatch.score.missing_head_tail)
                VP0(" Head/Tail missing=%d\n", cnMatch.score.missing_head_tail);
            if (cnMatch.score.extra_head_tail)
                VP0(" Head/Tail unused=%d\n", cnMatch.score.extra_head_tail);
            if (cnMatch.score.missing_atoms) {
                VP0(" Missing atoms in template=%d", cnMatch.score.missing_atoms);
                if (cnMatch.sUnmatchedNames[0]) VP0(" (%s)", cnMatch.sUnmatchedNames);
                VP0("\n");
            }
            if (cnMatch.score.missing_connect) {
                VP0(" Missing connects=%d:", cnMatch.score.missing_connect);
                for (int i = 0; i < cnMatch.score.missing_connect && i < MAXCONNECT; i++) {
                    const Pair *p = cnMatch.pairPUnmatched[i];
                    if (!p) continue;
                    RESIDUENAMEt *rnTo = &rnPResidues[p->to_group];
                    VP0(" %s.%s:%d.%s->%s.%s:%d.%s",
                        cPResName, rnPName->sChainId, rnPName->iPdbSequence,
                        anPAtoms[p->from_member].sName,
                        rnTo->sName, rnTo->sChainId, rnTo->iPdbSequence,
                        anPAtoms[p->to_member].sName);
                }
                VP0("\n");
            }
            if (cnMatch.score.extra_atoms) {
                VP0(" Extra atoms in template=%d", cnMatch.score.extra_atoms);
                ATOM    aAtom;
                LOOP    lAtoms = lLoop((OBJEKT)rTemplateRes, ATOMS);
                BOOL    first  = TRUE;
                while ((aAtom = (ATOM)oNext(&lAtoms)) != NULL) {
                    if (!bPDBHasHydrogens && aAtom->iAtomicNumber == HYDROGEN) continue;
                    BOOL found = FALSE;
                    for (int i = iFirstAtom; i <= iLastAtom; i++) {
                        if (!strcmp(sContainerName(aAtom),
                                    PVAI(prPPdb->vaAtomRecs, ATOMNAMEt, i)->sName))
                            { found = TRUE; break; }
                    }
                    if (!found) {
                        VP0("%s%s", first ? " (" : ",", sContainerName(aAtom));
                        first = FALSE;
                    }
                }
                if (!first) VP0(")");
                VP0("\n");
            }
            if (cnMatch.score.missing_connect > MAX_UNMATCHED_CONNECT) {
                VP0("Warning: %s.%s:%d has %d unmatched cross-link contacts"
                    " -- excess likely spurious:\n",
                    cPResName, rnPName->sChainId, rnPName->iPdbSequence,
                    cnMatch.score.missing_connect);
                for (int i = 0; i < cnMatch.score.missing_connect && i < MAXCONNECT; i++) {
                    const Pair *p = cnMatch.pairPUnmatched[i];
                    if (!p) continue;
                    RESIDUENAMEt *rnTo = &rnPResidues[p->to_group];
                    VP0("  %s.%s:%d.%s -> %s.%s:%d.%s\n",
                        cPResName, rnPName->sChainId, rnPName->iPdbSequence,
                        anPAtoms[p->from_member].sName,
                        rnTo->sName, rnTo->sChainId, rnTo->iPdbSequence,
                        anPAtoms[p->to_member].sName);
                }
            }
        }

        if (cnMatch.bHeadIsCrossLink) {
            RESIDUENAMEt *rnTo = &rnPResidues[cnMatch.pairPMatched[0]->to_group];
            ATOMNAMEt    *anTo = PVAI(prPPdb->vaAtomRecs, ATOMNAMEt,
                                      cnMatch.pairPMatched[0]->to_member);
            VPNOTE("Head atom %s.%s:%d.%s is cross-linked to %s.%s:%d.%s%s\n",
                   cPResName, rnPName->sChainId, rnPName->iPdbSequence,
                   sContainerName(aUnitHead(cnMatch.uTemplate)),
                   rnTo->sName, rnTo->sChainId, rnTo->iPdbSequence, anTo->sName,
                   ((cnMatch.pairPMatched[0]->to_group >  iResIndex) ? " (forward link deferred)":""));
        }
        if (cnMatch.bTailIsCrossLink) {
            RESIDUENAMEt *rnTo = &rnPResidues[cnMatch.pairPMatched[1]->to_group];
            ATOMNAMEt    *anTo = PVAI(prPPdb->vaAtomRecs, ATOMNAMEt,
                                      cnMatch.pairPMatched[1]->to_member);
            VPNOTE("Tail atom %s.%s:%d.%s is cross-linked to %s.%s:%d.%s%s\n",
                   cPResName, rnPName->sChainId, rnPName->iPdbSequence,
                   sContainerName(aUnitTail(cnMatch.uTemplate)),
                   rnTo->sName, rnTo->sChainId, rnTo->iPdbSequence, anTo->sName,
                   ((cnMatch.pairPMatched[1]->to_group >  iResIndex) ? " (forward link deferred)":""));
        }
        MESSAGE("Getting UNIT: %s\n", sContainerName(uTemplate));
        uResult           = (UNIT)oCopy((OBJEKT)uTemplate);
        *cnPMatchReturn = cnMatch;

    } else {

        VP0("Creating new UNIT for residue: %3s %2.2s%4d%c\n",
            cPResName, rnPName->sChainId,
            rnPName->iPdbSequence, rnPName->iCode);
        uResult = (UNIT)oCreate(UNITid);
        RESIDUE rRes = (RESIDUE)oCreate(RESIDUEid);
        ContainerSetName(rRes, cPResName);
        ContainerAdd((CONTAINER)uResult, (OBJEKT)rRes);
        memset(cnPMatchReturn, 0, sizeof(*cnPMatchReturn));
        for (int i = 0; i < iNumDetected; i++)
            cnPMatchReturn->pairPMatched[i] = pairPDetected[i];
    }

    VPTRACEEXIT(__func__);
    return uResult;
}


// Technically we should increment reference counts, but this is a temporary dictionary.
// Each variable that is a template (UNIT containing  asingle RESIDUE) is indexed
// by the RESIDUE name (property matching variants) and by the VARIABLE name (direct
// traditional matching)
//
// There is also a LEAP command naming scheme for direct access of CONNECT attributes,
// which is currently not used. (Get pure property rank matching, then reconsider this.)
//
// Next, index all ResMap entries. These can be FIRSTEND or LASTEND which are validated
// for TAIL-only or HEAD-only. They can also be unlabelled, in which case they could be
// polymer or non-polymer aliases and therefore HEAD+TAIL correctness cannot be verified.
static void
zPdbIndexResidueTemplates(  PDBREADt *prPPdb )
{
DICTLOOP        dlLoop;
DICTIONARY      dVariables;
    VPTRACEENTER(__func__);

    create_index(&(prPPdb->ixResMap), IX_DUPKEY, IX_LEN_CSTRING);
    dVariables = dVariablesDictionary();
    dlLoop = ydlDictionaryLoop(dVariables);
    while ( yPDictionaryNext(dVariables, &dlLoop ) ) {
        UNIT uTemplate = (UNIT)PDictLoopData(dlLoop);
        if ( iObjectType(uTemplate) != UNITid ) continue;
        if ( iContainerNumberOfChildren(uTemplate) > 1 ) continue;
        RESIDUE rRes = (RESIDUE)oContainerFirstObject(uTemplate);
        if ( iObjectType(rRes) != RESIDUEid ) continue; // unlikely
        // VARIABLE is a UNIT containing a single RESIDUE
        // == Residue Template signature
        // Index it by VARIABLE name and RESIDUE name
        char *sVarName = sDictLoopKey(dlLoop);
        char *sResName = sContainerName(rRes);
        IX_REC ix_record;
        ix_record.recptr = (IX_RECPOS) uTemplate;
        char *key = ix_record.key;
        strcpy(key, sResName);
        // Add both RESIDUE and VARIABLE indices
        add_key(&ix_record,&(prPPdb->ixResMap));
        if (strcmp(sVarName,sResName)) {
            strcpy(key, sVarName);
            add_key(&ix_record,&(prPPdb->ixResMap));
        }
    }
    // Insert PdbResMap match entries
    if (SdResidueNameMap) {
        dlLoop = ydlDictionaryLoop( SdResidueNameMap );
        while ( yPDictionaryNext( SdResidueNameMap, &dlLoop )) {
            char *cPData = (char*)PDictLoopData(dlLoop);
            UNIT uTemplate = (UNIT)oVariable(cPData);
            // Residue maps to undefined string
            if ( uTemplate == NULL) continue;
            char *cPKey = sDictLoopKey(dlLoop);
            IX_REC ix_record;
            int iTerm;
            zPdbFromNameMapKey( cPKey, &iTerm, ix_record.key );
            if (iTerm == FIRSTEND && bUnitHeadUsed(uTemplate)) {
                VPNOTE("PdbResMap FIRSTEND template %s contains Head atom\n",cPData);
                UnitSetHead(uTemplate,NULL);
            } else if (iTerm == LASTEND && bUnitTailUsed(uTemplate)) {
                VPNOTE("PdbResMap LASTEND template %s contains Tail atom\n",cPData);
                UnitSetTail(uTemplate,NULL);
            } else if (iTerm == NOEND && bUnitHeadUsed(uTemplate)!=bUnitTailUsed(uTemplate)) {
                VPNOTE("PdbResMap NOEND template %s has mismatched Head/Tail linkage\n",cPData);
            }
            // TODO: enforce HEAD,TAIL presence matches iTerm flag
            ix_record.recptr = (IX_RECPOS) uTemplate;
            add_key(&ix_record,&(prPPdb->ixResMap));
        }
    }
    VPTRACEEXIT(__func__ );
}


/*
 *      zPdbCreateUnit (From zPdbReadAndCreateUnit)
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
zPdbCreateUnit( PDBREADt *prPPdb )
{
int             iAdd, iAddH = 0, iAddHeavy = 0, iAddUnk = 0, iCreate = 0, nAtoms;
ATOM            aA, aB, *aPa, aHead, aTail;
MatchCandidate  cnMatch;
UNIT            uNew;
BOOL            bStartChain;
int             iChainCount=0;

    VPTRACEENTER("zPdbCreateUnit" );

    nAtoms = iVarArrayElementCount(prPPdb->vaAtomRecs);
    if ( !nAtoms ) {
        VPFATALEXIT("\tNo atoms!\n" );
        return;
    }
    // This array is an atom pointer array indexed by atomSerial; only used for CONECT, LINK lookups
    if ( GDefaults.bPdbUseLinkRecords || GDefaults.bPdbUseConect) {
        VarArraySetSize( prPPdb->vaAtoms, prPPdb->iMaxSerialNum );
        aPa = PVAI( prPPdb->vaAtoms, ATOM, 0 );
        for (int i=0; i < prPPdb->iMaxSerialNum; i++, aPa++ ) *aPa = NULL;
    }

    // Main loop over residues and atoms:
    int iResidueCount = iVarArrayElementCount( prPPdb->vaResidues );
    //int iAtomCount = iVarArrayElementCount( prPPdb->vaAtomRecs );
    for (int iRes = 0; iRes < iResidueCount; iRes++) {
        RESIDUENAMEt *rnPName = PVAI( prPPdb->vaResidues, RESIDUENAMEt, iRes);

        if (GDefaults.bPdbAutoMatch) {
            uNew = zPdbMatchResidueTemplate(prPPdb, iRes, &cnMatch );
            bStartChain = !cnMatch.pairPMatched[0] || cnMatch.bHeadIsCrossLink;
        } else {
            uNew = zuPdbGetNextUnit( prPPdb, &bStartChain, rnPName->iPdbSequence );
        }

        RESIDUE rRes = (RESIDUE) oContainerFirstObject(uNew);
        prPPdb->rRes = rRes;
        rnPName->rResidue = rRes;

        BuildInternalsForContainer( &(uNew->cHeader), 0, ATOMPOSITIONKNOWN );
        //FIXME: don't need this?  ContainerSetNextChildsSequence( prPPdb->uUnit, iRes);

        // Save container name for AtomName mapping
        STRING sTemplateName;
        strcpy(sTemplateName,sContainerName(uNew));
                /* If the new UNIT is the first UNIT in a chain or non-polymer, */
                /* then don't sequence it to the last UNIT, just join them */
        if ( bStartChain ) {

            UnitJoin( prPPdb->uUnit, uNew );

        } else {

            // Save Head, Tail before join so we can check distances
            aTail = aUnitTail(prPPdb->uUnit);
            aHead = aUnitHead(uNew);

                /* Build INTERNAL coordinates for the linkage between the */
                /* UNITs that are about to be formed */
                /* Sequence the new UNIT into the PDB UNIT */
            // bUnitHeadUsed(uNew)==0 means the new residue does not have HEAD available!

            BuildInternalsBetweenUnitsUsingFlags( prPPdb->uUnit, uNew, 0, 0 );
            UnitSequence( prPPdb->uUnit, uNew );
        }
            /* Define the PDB chainID, resSeq, iCode */

        if ( bStartChain && GDefaults.bPdbResetChainID) {
            if (iChainCount < CHAINID_LIST_LEN ) {
                rRes->sChainId[0]=GsChainIdList[iChainCount];
                rRes->sChainId[1]=0;
            } else {
                int i1 = iChainCount/CHAINID_LIST_LEN;
                int i2 = iChainCount%CHAINID_LIST_LEN;
                if (i1==CHAINID_LIST_LEN && i2==0) VPWARN("2-char ChainID overflowing to special chars\n");
                else if (i1==CHAINID_LIST_LEN) {
                    VPFATAL("2-char ChainID overflow\n");
                    i1=0; iChainCount=0;
                }
                rRes->sChainId[0]=GsChainIdList[i1];
                rRes->sChainId[1]=GsChainIdList[i2];
                rRes->sChainId[2]=0;
            }
        } else strcpy( rRes->sChainId, rnPName->sChainId);
        rRes->iPdbResSeq = rnPName->iPdbSequence;
        rRes->cICode = rnPName->iCode;

                /* Apply ATOM information into the new residue */
        int iFirstAtom = rnPName->iFirstAtom;
        int iLastAtom = (iRes < iResidueCount - 1) ? ((rnPName + 1)->iFirstAtom - 1)
                                 : iVarArrayElementCount(prPPdb->vaAtomRecs) - 1;
        for (int iAtom = iFirstAtom; iAtom <= iLastAtom; iAtom++) {
            // ATOM processing part from AddAtom() // FIXME: use already matched index
            ATOMNAMEt *anAtom = PVAI( prPPdb->vaAtomRecs, ATOMNAMEt, iAtom);
            ATOM aAtom = (ATOM)cContainerFindName( (CONTAINER)rRes, ATOMid, anAtom->sName );

            /*
             *  If atom name doesn't match, try aliases
             */
            if ( aAtom == NULL ) {
                char *cPAtomName = zcPPdbMapName(SdAtomNameMap, NOEND, anAtom->sName, rRes);
                if (cPAtomName )
                        aAtom = (ATOM)cContainerFindName( (CONTAINER)rRes, ATOMid, cPAtomName);
            }

                        /* If the ATOM was not within the RESIDUE report */
                        /* that you are going to create it */
            if ( aAtom == NULL ) {
                aAtom = (ATOM)oCreate(ATOMid);
                ContainerAdd( (CONTAINER)rRes, (OBJEKT)aAtom );
                ContainerSetName( (CONTAINER)aAtom, anAtom->sName );
                AtomSetElement( aAtom, anAtom->iElement );
                MESSAGE("Read atom: %s and adding it to: %s\n", anAtom->sName,
                                sContainerName((CONTAINER)rRes) );
                STRING sTemp;
                VP0("Created a new atom named: %s within residue: %s\n",
                        anAtom->sName, sContainerFullDescriptor((CONTAINER)rRes,sTemp) );
                iCreate ++;
            } else {
                MESSAGE("Read atom: %s and found it in the RESIDUE\n", anAtom->sName );
                if (iAtomElement(aAtom)==NOELEMENT) AtomSetElement( aAtom, anAtom->iElement );
            }

            if ((GDefaults.bPdbUseLinkRecords || GDefaults.bPdbUseConect)
                    // Leap-built atoms get -1 in parmed, ignore them
                    && anAtom->iAtomSerial >= 0 ) {
                *PVAI( prPPdb->vaAtoms, ATOM, anAtom->iAtomSerial ) = aAtom;
            }

                                /* Define its position */
            VECTOR vPos;
            VectorDef( &vPos, anAtom->x, anAtom->y, anAtom->z );
            AtomSetPosition( aAtom, vPos );
        }

        if (!bStartChain && !GDefaults.bPdbAutoMatch &&
                    GDefaults.dPdbLinkCovalentCutoff > 0 && aHead && aTail) {
            float dX = aTail->vPosition.dX - aHead->vPosition.dX;
            float dY = aTail->vPosition.dY - aHead->vPosition.dY;
            float dZ = aTail->vPosition.dZ - aHead->vPosition.dZ;
            float d2 = dX*dX + dY*dY + dZ*dZ;
            if (!bAtomsBondedDist(iAtomElement(aTail),iAtomElement(aHead), d2,
                                  GDefaults.dPdbLinkCovalentCutoff )) {
                VP1("Starting new chain because dist=%g > cutoff=%g\n",
                                   sqrt(d2), GDefaults.dPdbLinkCovalentCutoff);
	        AtomRemoveBond( aTail, aHead );
            }
        }

        if (GDefaults.bPdbAutoMatch)
            zProcessModifications(prPPdb, &cnMatch, iFirstAtom, iLastAtom);

        zPdbBuildCoordinatesForContainer( (CONTAINER)rRes, &iAddH, &iAddHeavy, &iAddUnk );

    }

    if (GDefaults.bPdbUseConect) {
        for (int i=0; i<iVarArrayElementCount(prPPdb->vaConectRecs); i++) {
            struct pdb_conect *conect = PVAI( prPPdb->vaConectRecs, struct pdb_conect, i );
            VPTRACE("Process PDB_CONECT record for atom %i.\n",
                    conect->serial_num );
                /*
                 * prPPdb->iMaxSerialNum is the size of prPPdb->vaAtoms and is
                 * 1 more than the last inputted atom's iSerialNum; thus, since
                 * serial numbers start at 1 and prPPdb->vaAtoms is zero based,
                 * it has wasted space at the start but none at the end.
                 * So guard against out of bounds high indexing.  srb 5-2022.
                 */
            if (conect->serial_num >= prPPdb->iMaxSerialNum ) {
                VPWARN("Ignoring CONECT record for atom serial number "
                    "(%i) that is greater than\nthe maximum serial"
                    " number inputted (%i) from the pdb file.\n",
                    conect->serial_num, prPPdb->iMaxSerialNum - 1 );
                continue;
            }
            aA = *PVAI( prPPdb->vaAtoms, ATOM, conect->serial_num );
            if ( aA == NULL ) {
                VPWARN("Invalid CONECT record (atomSerial=%d) in pdb file.\n", conect->serial_num );
                continue;
            }
            for (int j=0; j<4; j++ ) {
                if ( conect->covalent[j] == 0 ) continue;
                if ( conect->covalent[j] >= prPPdb->iMaxSerialNum ) {
                    VPWARN("In CONECT record for atom serial number (%i)\n"
                        "ignoring bonded atom serial number (%i) that is "
                        "greater than\nthe maximum serial"
                        " number inputted (%i) from the pdb file.\n",
                        conect->serial_num, conect->covalent[j],
                        prPPdb->iMaxSerialNum - 1 );
                    continue;
                }
                aB = *PVAI(prPPdb->vaAtoms,ATOM, conect->covalent[j]);
                if ( aB == NULL ) {
                    VPWARN("Invalid CONECT record (atomSerial=%d) in pdb file.\n", conect->covalent[j] );
                    continue;
                }
                if ( !bAtomBondedTo( aA, aB ) ) {
                    if (GDefaults.bPdbLinkIons || bAtomsBondedDist(iAtomElement(aA),iAtomElement(aB),1.0,1.0) )
                        AtomBondTo( aA, aB );
                    else {
                        RESIDUE rA = (RESIDUE)cContainerWithin(aA);
                        RESIDUE rB = (RESIDUE)cContainerWithin(aB);
                        VP2("PDB CONECT %s %s %s%d to %s %s %s%d excluded (ionic)\n",
                               sContainerName(aA),sContainerName(rA),
                               sResidueChainId(rA),iResiduePdbSequence(rA),
                               sContainerName(aB),sContainerName(rB),
                               sResidueChainId(rB),iResiduePdbSequence(rB));
                    }
                }
            }
        }
    }
    if (GDefaults.bPdbUseLinkRecords) {
        VP0("Number LINK+SSBOND records parsed: %d\n",iVarArrayElementCount(prPPdb->vaLinkRecs));
        for (int i=0; i<iVarArrayElementCount(prPPdb->vaLinkRecs); i++) {
            struct pdb_link *link = PVAI( prPPdb->vaLinkRecs, struct pdb_link, i );
            VPTRACE("Process LINK/SSBOND record %i.\n", i );
            ATOM aAtom[2];
            for (int j=0;j<2;j++) {
                IX_REC rec = { NULL };
                ATOMKEYt *key = (ATOMKEYt *)rec.key;
                key->resSeq = link->residues[j].seq_num;
                if (link->residues[j].chain_id[0] != ' ') {
                    memcpy(key->chainID,link->residues[j].chain_id,2);
                } else {
                    key->chainID[1]=0;
                    if (link->residues[j].chain_id[1] != ' ')
                        key->chainID[0]=link->residues[j].chain_id[1];
                    else key->chainID[0]=0;
                }

                if (!key->chainID[0]) key->chainID[1]=0; // double NUL for blank chain
                key->iCode = link->residues[j].insert_code;
                memset(key->name,0,sizeof(key->name));
                if (link->name[j][0]==' ') strcpy(key->name,link->name[j]+1);
                else strcpy(key->name,link->name[j]);
                if ( locate_key( &rec,  &prPPdb->ixAtomIndex, 1 ) != IX_OK ) {
                    VPWARN("Unable to locate atom %d of link record %d: %s.%s:%d.%s\n",j,i,
                            link->residues[j].name,link->residues[j].chain_id,
                            link->residues[j].seq_num,link->name[j]);
                } else {
                    ATOMNAMEt *anAtom = (ATOMNAMEt *)rec.recptr;
                    aAtom[j] = *PVAI( prPPdb->vaAtoms, ATOM, anAtom->iAtomSerial );
                }
            }
            if ( aAtom[0] && aAtom[1] &&!bAtomBondedTo( aAtom[0], aAtom[1] ) ) {
                if (GDefaults.bPdbLinkIons || bAtomsBondedDist(iAtomElement(aAtom[0]),iAtomElement(aAtom[1]),1.0,1.0) ) {
                    AtomBondTo( aAtom[0], aAtom[1] );
                } else {
                    RESIDUE rA1 = (RESIDUE)cContainerWithin(aAtom[0]);
                    RESIDUE rA2 = (RESIDUE)cContainerWithin(aAtom[1]);
                    VP2("PDB LINK %s.%s:%d.%s to %s.%s:%d.%s excluded (ionic)\n",
                           sContainerName(rA1), sResidueChainId(rA1),
                           iResiduePdbSequence(rA1), sContainerName(aAtom[0]),
                           sContainerName(rA2), sResidueChainId(rA2),
                           iResiduePdbSequence(rA2), sContainerName(aAtom[1]) );
                }
            }
        }
    }

                /* Build geometry for everything that has not */
                /* been build */
    zPdbBuildCoordinatesForContainer( (CONTAINER)prPPdb->uUnit,
                                &iAddH, &iAddHeavy, &iAddUnk );

    VP0("  total atoms in file: %d\n", nAtoms );
    if ( (iAdd = iAddH + iAddHeavy + iAddUnk) ) {
        VP0("  Leap added %d missing atom%s according to residue templates:\n",
                                iAdd, (iAdd > 1 ? "s" : "") );
        if ( iAddHeavy )
                VP0("       %d Heavy\n", iAddHeavy );
        if ( iAddH )
                VP0("       %d H / lone pairs\n", iAddH );
        if ( iAddUnk )
                VP0("       %d unknown element\n", iAddUnk );
    }
    if ( iCreate )
        VP0("  The file contained %d atoms not in residue templates\n",
                                                        iCreate );
    if ( iAdd  &&  iAdd == iCreate ) {
        VPWARN("Since the number of added atoms equals the number of missing "
                "atoms, it is likely\n" "that some atoms had incorrect names; "
                "you may want to use addPdbAtomMap to map\n"
                "these names, or change the names in the PDB file.\n\n" );
    }
    VPTRACEEXIT("zPdbCreateUnit" );
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
uPdbRead( FILE *fPdb, VARARRAY vaUnits, BOOL bFormatCIF )
{
PDBREADt        prPdb = {0};
int             i;
    VPTRACEENTER("uPdbRead" );

    prPdb.fPdbFile = fPdb;
    prPdb.vaUnits = vaUnits;

    if (bFormatCIF) {
        VP0("CifRead:   Importing CIF format coordinate file\n");
        VP0("CifRead:   Prefer %s CIF columns\n",
            GDefaults.bCIFReadAuth ? "Author" : "Standardized");
    }
    if (GDefaults.bPdbAutoMatch) VP0("PdbRead:   AutoMatch enabled, crosslink cutoff=%.3f\n",
                   GDefaults.dPdbCrosslinkCovalentCutoff);
    else VP0("PdbRead:   AutoMatch disabled\n" );
    if (GDefaults.sPdbPatchFilename[0]) {
        VP0("PdbRead:   AutoMatch patch file = %s\n",GDefaults.sPdbPatchFilename);
        prPdb.fpPatchFileOut = fopen(GDefaults.sPdbPatchFilename,"w");
        if (!prPdb.fpPatchFileOut) {
            VPWARN("Cannot open autolink Patch file\n");
        }
    }
    VP0("PdbRead:   Tail-Head link cutoff=%.3f\n",GDefaults.dPdbLinkCovalentCutoff);
    if (GDefaults.bPdbAutoLoadRes)
        VP0("PdbRead:   AutoLoad enabled (attempt to load residue parameter files)\n" );
    if (!bFormatCIF)
        VP0("PdbRead:   CONECT processing %sabled\n", GDefaults.bPdbUseConect ? "en":"dis");
    VP0("PdbRead:   LINK,SSBOND processing %sabled\n", GDefaults.bPdbUseLinkRecords ? "en":"dis");
    if ( (!bFormatCIF && GDefaults.bPdbUseConect) || GDefaults.bPdbUseLinkRecords )
        VP0("PdbRead:   LINK/CONECT to metal ions will be %scluded\n", GDefaults.bPdbLinkIons ? "in":"ex");
    if ( GDefaults.iPdbReadModel == -2)
        VP0("PdbRead:   All MODELs retained with 2-char ChainId assignment\n");
    if ( GDefaults.iPdbReadModel == -1)
        VP0("PdbRead:   All MODELs retained with default ChainId assignment\n");
    else if ( GDefaults.iPdbReadModel == 0)
        VP0("PdbRead:   First MODEL will be retained, if MODEL is present\n");
    else
        VP0("PdbRead:   MODEL %d will be retained\n", GDefaults.iPdbReadModel);
    if (GDefaults.bPdbExpandBioMt)
        VP0("PdbRead:   REMARK BIOMT will be read and processed, MTRIXn will be ignored\n" );
    else if (GDefaults.bPdbExpandNCSMt)
        VP0("PdbRead:   MTRIXn records will be processed, if present\n" );
    else
        VP0("PdbRead:   MTRIXn records will not be processed\n" );

    prPdb.uUnit = (UNIT)oCreate(UNITid);
    ContainerSetName(prPdb.uUnit, "PDB_UNIT");
    prPdb.vaResidues = vaVarArrayCreate( sizeof(RESIDUENAMEt) );
    prPdb.vaAtomRecs  = vaVarArrayCreate( sizeof(ATOMNAMEt) );
    create_index( &prPdb.ixAtomIndex, IX_NODUPKEY, sizeof(ATOMKEYt));
    prPdb.vaMatrices = vaVarArrayCreate( sizeof(PDBMATRIXt) );
    if ( GDefaults.bPdbUseLinkRecords || GDefaults.bPdbUseConect)
        prPdb.vaAtoms = vaVarArrayCreate( sizeof(ATOM) );
    if (GDefaults.bPdbUseConect) prPdb.vaConectRecs = vaVarArrayCreate( sizeof(struct pdb_conect) );
    if (GDefaults.bPdbUseLinkRecords) prPdb.vaLinkRecs = vaVarArrayCreate( sizeof(struct pdb_link) );
    //create_index(&(prPdb.ixPatchFilenames), IX_NODUPKEY, IX_LEN_CSTRING);

                /* Scan the PDB file looking for residue names, */
                /* residue termination status, and matrices */

    if (bFormatCIF) CifReadFile( &prPdb );
    else zPdbReadFile( &prPdb );

    RESIDUENAMEt *rnRes = PVAI( prPdb.vaResidues, RESIDUENAMEt, 0);
    int iNumPdbResidues = iVarArrayElementCount(prPdb.vaResidues);
    if (GDefaults.bPdbAutoMatch) {
        prPdb.vaPoints = vaVarArrayCreate( sizeof(Point) );
        MALLOC( prPdb.iPGroupStart, unsigned int *, sizeof(int)*( iNumPdbResidues + 1) );
        prPdb.iPGroupStart[0] = 0;
    }
    int iDuplicates = 0;
    int iResIndex = -1;
    for (int iAtom=0; iAtom<iVarArrayElementCount(prPdb.vaAtomRecs); iAtom++) {
        ATOMNAMEt *anAtom = PVAI( prPdb.vaAtomRecs, ATOMNAMEt, iAtom);
        if (anAtom->iResNameIndex != iResIndex) {
            iResIndex = anAtom->iResNameIndex;
            rnRes = PVAI( prPdb.vaResidues, RESIDUENAMEt, anAtom->iResNameIndex);
            if (GDefaults.bPdbAutoMatch) {
                prPdb.iPGroupStart[iResIndex] = iVarArrayElementCount(prPdb.vaPoints);
            }
        }
        // Exclude ions (atom name = resName), water, and hydrogen in crosslink analysis
        if (GDefaults.bPdbAutoMatch
                && strcmp(anAtom->sName, rnRes->sName)
                && strcmp("HOH", rnRes->sName) // FIXME: better water matching, use SOLVENT template names!
                && strcmp("WAT", rnRes->sName) // FIXME: better water matching, use SOLVENT template names!
                && anAtom->iElement != HYDROGEN ) {
            // Not water or hydrogen, add to neighbor point lookup grid array:
            Point point;
            point.x = anAtom->x;
            point.y = anAtom->y;
            point.z = anAtom->z;
            point.group = iResIndex;
            point.member = iAtom;
            VarArrayAdd( prPdb.vaPoints, (GENP)&point );
        }

        IX_REC rec = { (IX_RECPOS)anAtom };
        ATOMKEYt *key = (ATOMKEYt *)rec.key;
        key->resSeq = rnRes->iPdbSequence;
        memcpy(key->chainID,rnRes->sChainId,2);
        key->chainID[1]=0; // double NUL for blank chain
        key->iCode = rnRes->iCode;
        memset(key->name,0,sizeof(key->name));
        strcpy(key->name,anAtom->sName);
        if ( add_key( &rec,  &prPdb.ixAtomIndex ) != IX_OK ) {
            if (iDuplicates++ < 50) {
                VP0("-- residue %d: duplicate [%4s %2.2s%4d%c] atom\n",
                        anAtom->iResNameIndex, key->name, key->chainID,
                             key->resSeq, key->iCode );
            }
        }

    }
    // add array sentinel for non-bond search grid points
    if (GDefaults.bPdbAutoMatch) {
        prPdb.iPGroupStart[iResIndex+1] = iVarArrayElementCount(prPdb.vaPoints);
    }

    if ( iDuplicates ) {
        VPWARN("Atom names in each residue should be unique. %d duplicate atoms found.\n",iDuplicates );
        VP0("     (Same-name atoms are handled by using the first\n" );
        VP0("      occurrence and by ignoring the rest.\n" );
        VP0("      Many instances of duplicate atom names usually come\n" );
        VP0("      from alternate conformations in the PDB file.)\n\n" );
    }
    if (GDefaults.bPdbAutoMatch) {
        prPdb.ngBondGrid = neighbor_grid_setup( PVAI(prPdb.vaPoints, Point, 0),
                iVarArrayElementCount(prPdb.vaPoints), iNumPdbResidues, prPdb.iPGroupStart,
                MAX_BOND_LENGTH * GDefaults.dPdbLinkCovalentCutoff ); // max expected bond length * cutoff multiplier

        // Generate an index into all defined Residue Templates
        zPdbIndexResidueTemplates( &prPdb );

    } else {

        // Use vaUnits, resName and PdbMap to create RESIDUE array
        zPdbConvertNamesAndSequenceNumbers( &prPdb );

        if (vaUnits) {
            int     mismatch, total = 0;
            /*
             *  loadPdbUsingSeq: see if residues match
             */
            VP0("  matching pdb residues -> sequence template\n" );
            VP0("\tres\tpdb\ttemplate\n" );
            for (i=0; i<iVarArrayElementCount(vaUnits); i++) {
                if ( ! strstr( (*PVAI( vaUnits, UNIT, i ))->cHeader.sName,
                               PVAI(prPdb.vaResidues, RESIDUENAMEt, i)->sName) )
                        mismatch = 1;
                else
                        mismatch = 0;
                VP0("\t%d%s\t%s\t%s\n", i+1, (mismatch ? "*" : ""),
                        PVAI(prPdb.vaResidues, RESIDUENAMEt, i)->sName,
                        (*PVAI( vaUnits, UNIT, i ))->cHeader.sName );
                total += mismatch;
            }
            if (total) {
                    VP0("  * = possible mismatch; total %d\n", total );
                    VP0("      (i.e. pdb name not a substring of template)\n\n");
            }
        }
    }

    zPdbCreateUnit(&prPdb);
    if (prPdb.fpPatchFileOut) fclose(prPdb.fpPatchFileOut);

                /* Build the symmetry related monomers */

    if (GDefaults.bPdbExpandBioMt || 
            GDefaults.bPdbExpandNCSMt ||
            GDefaults.bPdbExpandSymm)
        zPdbCreateSymmetryRelatedMonomers( &prPdb );

                /* Clean up */

    VarArrayDestroy(&prPdb.vaResidues);
    VarArrayDestroy(&prPdb.vaMatrices);
    VarArrayDestroy(&prPdb.vaAtomRecs);
    destroy_index(&prPdb.ixAtomIndex );
    if (prPdb.vaAtoms) VarArrayDestroy(&prPdb.vaAtoms);
    if (prPdb.vaConectRecs) VarArrayDestroy(&prPdb.vaConectRecs);
    if (prPdb.vaPoints) VarArrayDestroy(&prPdb.vaPoints);
    if (prPdb.ngBondGrid) neighbor_grid_free(prPdb.ngBondGrid);
    if (prPdb.ixResMap.root) destroy_index(&prPdb.ixResMap);
    if (prPdb.iPGroupStart) FREE( prPdb.iPGroupStart );
    if (bFormatCIF) ndb_cif_close();

    VariableSet( "PDB_UNIT", (OBJEKT)prPdb.uUnit );
    ContainerSetName(prPdb.uUnit, "default_name");
    VPTRACEEXIT("uPdbRead" );
    return ( prPdb.uUnit );
}

