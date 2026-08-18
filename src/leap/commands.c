/*
 *      File:   commands.c
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
 *             David A. Rivkin                                          *
 *                                                                      *
 *     Principal Investigator: Peter A. Kollman                         *
 *                                                                      *
 ************************************************************************
 *
 *      All commands typed into the parser are
 *      executed here.
 *
 *      Each function name is prefixed with Cmd
 *      and has two arguments, an integer with
 *      the number of arguments, and an array
 *      of objects which are the arguments.
 *
 */

/*
 *TODO: Add the following commands to the documentation.
 *      New command line options.
 *      addPath
 *      measure
 *      select   - change this one
 *      deSelect
 *
 *
 * This help text pattern:
 *
 *        hTemp = hHelp( "addions" );
 *        if ( hTemp == NULL ) {
 *              VPFATALDELAYEDEXIT("No help available on addIons\n" );
 *        } else {
 *              VPFATALDELAYEDEXIT("\n%s\n", sHelpText(hTemp) );
 *        }
 *
 * Should be replaced by:
 *        VPFATALDELAYEDEXIT("%s",sHelpText("addIons") );
 *
 *  where 1) sHelptext() searche tolower(arg) and the return includes newline character
 *        2) if nor found, return a pointer to buffer with sprintf(buffer,"No help available on %s\n",arg)
 */

/*        Modifications
 *        Christine Cezard  (2007)
 *        Universite de Picardie - Jules Verne, Amiens
 *        http://q4md-forcefieldtools.org
 *
 *        Added SaveMol2
 *              Flip
 *           &  Relax
 *              commands
 *        Added case i (integer) in function bCmdMatchTypes
 *
 *
 *        Robin Betz (2011, 2013)
 *        UC San Diego
 *
 *        Added addIonsRand
 *        Added saveAmberParmNetCDF
 *
*/

/*        Modifications
 *        Mason Louchart  (2011)
 *        Universite de Picardie - Jules Verne, Amiens
 *        http://q4md-forcefieldtools.org
 *
 *        Added SaveMol3
 *              LoadMol3
 */

#include        <stdio.h>
#include        <stdlib.h>
#include        <float.h>
#include        <string.h>

#include        "basics.h"
#include        "vector.h"
#include        "matrix.h"
#include        "classes.h"
#include        "dictionary.h"
#include        "database.h"
#include        "library.h"
#include        "parmLib.h"
#include        "pdbFile.h"
#include        "help.h"
#include        "parser.h"
#include        "tools.h"
#include        "amber.h"
#include        "commands.h"
#include        "defaults.h"
#include        "leap.h"
#include        "octree.h"
#include        "tripos.h"
#include        "build.h"
#include        "zMatrix.h"
#include        "unitio.h"
#include        "mol2File.h"
#include        "mol3File.h"
#include        "minimizer.h"
#include        "model.h"
#include        "elements.h"
#include        "select_mask.h"

int     iMemDebug = 0;

//extern DICTIONARY       GdVariables;

#define ATOMSINBOND     2
#define ATOMSINANGLE    3
#define ATOMSINTORSION  4

void VarArrayDeleteMore( VARARRAY header, int pos, int num);
void SelectRelaxInFramework(UNIT uUnit, MINIMIZER mMinimizer);

/*
 *  COMMANDt    cCommands[]
 *
 *      table mapping command oCmd_ routines to strings
 *      - moved to end of this file.
 */




/*
 *---------------------------------------------------------------------
 *
 *        Private routines
 *
 */


/*
 *      bCmdMatchTypes
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return true if the type of the OBJEKT matches one
 *      of the types in the string up to the first space or
 *      end of string.
 *
 *      (*sPNeedType) will return the name of the type required.
 *
 *      See 'bCmdGoodArguments' for a description of the type characters.
 */
static  bool
bCmdMatchTypes( OBJEKT oObj, char **sPTypes, char *sNeedType )
{
char    *sTemp;

    strcpy( sNeedType, "" );

    sTemp = *sPTypes;
    while ( (**sPTypes!=' ') && (**sPTypes!='\0') ) {
        (*sPTypes)++;
    }
    if ( **sPTypes==' ' )
        (*sPTypes)++;
    do {
        switch ( *sTemp ) {
            case '*':
                return true;
            case 'u':
                if ( iObjectType(oObj) == UNITid )
                        return true;
                strcat( sNeedType, "unit" );
                break;
            case 'm':
                if ( iObjectType(oObj) == MOLECULEid )
                        return true;
                strcat( sNeedType, "molecule" );
                break;
            case 'r':
                if ( iObjectType(oObj) == RESIDUEid )
                        return true;
                strcat( sNeedType, "residue" );
                break;
            case 'a':
                if ( iObjectType(oObj) == ATOMid )
                        return true;
                strcat( sNeedType, "atom" );
                break;
            case 'l':
                if ( iObjectType(oObj) == LISTid )
                        return true;
                strcat( sNeedType, "list" );
                break;
            case 'n':
                if ( iObjectType(oObj) == ODOUBLEid )
                        return true;
                strcat( sNeedType, "number" );
                break;
            case 'i':
                if ( iObjectType(oObj) == OINTEGERid )
                        return true;
                strcat( sNeedType, "integer" );
                break;
            case 's':
                if ( iObjectType(oObj) == OSTRINGid )
                        return true;
                strcat( sNeedType, "string" );
                break;
            case 'p':
                if ( iObjectType(oObj) == PARMSETid )
                        return true;
                strcat( sNeedType, "parameter_set" );
                break;
            case 'z':
                if ( oObj == NULL )
                        return true;
                strcat( sNeedType, "null" );
                break;
            case '\0':
            case ' ':
                return false;
            default:
                DFATAL("ILLEGAL type character in bCmdMatchTypes" );
        }
        sTemp++;
        if ( *sTemp != ' ' && *sTemp ) {
            strcat( sNeedType, " " );
        }
    } while ( *sTemp != '\0' );
    return false;
}


/*
 *      bCmdGoodArguments
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return true if iArgCount equals iMixNeeded...iMaxNeeded, otherwise return
 *      false and print an error message.
 *      Also check the types of the arguments against the type string
 *      A string of characters separated by spaces.
 *      The characters corresponding to:
 *
 *      *       - Matches all types.
 *      u       - Matches a UNIT.
 *      m       - Matches a MOLECULE.
 *      r       - Matches a RESIDUE.
 *      a       - Matches an ATOM.
 *      l       - Matches a LIST.
 *      n       - Matches a NUMBER (ODOUBLE).
 *      i       - Matches an INTEGER
 *      s       - Matches an OSTRING.
 *      p       - Matches a PARMSET.
 *      z       - Matches a NULL. (literal string "null")
 *      ,       - Option list allowed to end here
 *
 *      Eg:  "* umra l s, s" will match anything in [0],
 *              a UNIT/MOLECULE/RESIDUE/ATOM in [1],
 *              a LIST in [2], and an OSTRING in [3],
 *              and optionally OSTRING in [4].
 *
 *      The arguments in the oaArgs array are either OBJEKTs or
 *      they are ASSOCs.  If they are ASSOCs then the OBJEKT within
 *      the ASSOC is tested.
 */
static  bool
bCmdGoodArguments( char *sCmd, int iArgCount, ASSOC aaArgs[], char *sTypes )
{
int             i;
OBJEKT          oObj;
int             iMinNeeded, iMaxNeeded;
STRING          sNeed;

        /* Count how many arguments are required */

    if ( strlen(sTypes) == 0 ) {
        iMinNeeded = 0;
        iMaxNeeded = 0;
    } else {
        iMinNeeded = -1;
        iMaxNeeded = 1;
        for ( i=0; i<strlen(sTypes); i++ ) {
            if ( sTypes[i] == ' ' )
                iMaxNeeded++;
            if ( sTypes[i] == ',' )
                iMinNeeded = iMaxNeeded;
        }
        if (iMinNeeded<0) iMinNeeded = iMaxNeeded;
    }

                /* Check the number of arguments */

    if ( iArgCount < iMinNeeded || iArgCount > iMaxNeeded ) {
        VPFATAL("%s: Improper number of arguments!\n", sCmd );
        return false;
    }

                /* Check the types of arguments */
                /* If the OBJEKT is an ASSOC then get the OBJEKT */
                /* within the ASSOC */

    for ( i=0; i<iArgCount; i++ ) {
        oObj = OBJEKT_from(aaArgs[i]);
        if ( iObjectType(oObj) == ASSOCid )
                oObj = oAssocObject(oObj);
        if ( !bCmdMatchTypes( oObj, &sTypes, sNeed ) ) {
            VPFATAL("%s: Argument #%d is of type %s must be of type: [%s]\n"
                      "    Here are some suggestions for correcting this error:\n"
                      "    Verify the type of an argument with the desc command.\n"
                      "    Check for alternate argument names with the list command.\n",
                sCmd, i+1, sObjectType(oObj), sNeed );
            return false;
        }
    }
    return true;
}


#if 0
/*
 *      vaCmdListToVarArray
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Copy the OBJEKTs from the LIST into a VARARRAY and
 *      return the VARARRAY.
 */
static VARARRAY
vaCmdListToVarArray( LIST lList )
{
VARARRAY        vaList;
LISTLOOP        llElements;
int             i, iSize;
OBJEKT          oObj;

    iSize = iListSize(lList);
    vaList = vaVarArrayCreate( sizeofOBJEKT_from( ));
    VarArraySetSize( vaList, iSize );

    i = 0;
    llElements = llListLoop(lList);
    while ( oObj = oListNext(&llElements) ) {
        *PVAI( vaList, OBJEKT, i ) = oObj;
    }
    return vaList;
}
#endif






/*
 *=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 *
 *      Actual commands
 *
 */



/*
 *      oCmd_quit
 *
 *      Author: Christian Schafmeister (1991)
 *
 */
OBJEKT
oCmd_quit( int iArgCount, ASSOC aaArgs[] )
{
    if ( !bCmdGoodArguments( "quit", iArgCount, aaArgs, (char *)"" ) ) {
        VPFATALDELAYEDEXIT("usage:  quit\n" );
        return NULL;
    }

    GrMainResult.iCommand = CMD_QUIT;
    VP0("\tQuit\n" );
    VPEPILOG( );
    return NULL;
}



/*
 *      oCmd_describe
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Describe the first OBJEKT in the array.
 */
OBJEKT
oCmd_describe( int iArgCount, ASSOC aaArgs[] )
{
    if ( !bCmdGoodArguments( "desc", iArgCount, aaArgs, "*" ) ) {
        VPFATALDELAYEDEXIT("usage:  desc <variable>\n" );
        return NULL;
    }

    Describe( oAssocObject(aaArgs[0]) );

    return NULL;
}



/*
 *      oCmd_debugOn
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Turn on debugging.
 */
OBJEKT
oCmd_debugOn( int iArgCount, ASSOC aaArgs[] )
{
    if ( !bCmdGoodArguments( "debugOn", iArgCount, aaArgs, "s" ) ) {
        VPFATALDELAYEDEXIT("usage:  debugOn <filename>\n" );
        return NULL;
    }

    MessageAddFile( sOString(oAssocObject(aaArgs[0])) );
    VP0("Messages will be displayed from the files:\n" );
    MessageFileList();
    return NULL;
}


/*
 *      oCmd_debugOff
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Turn off debugging.
 */
OBJEKT
oCmd_debugOff( int iArgCount, ASSOC aaArgs[] )
{
    if ( !bCmdGoodArguments( "debugOff", iArgCount, aaArgs, "s" ) ) {
        VPFATALDELAYEDEXIT("usage:  debugOff <filename>\n" );
        return NULL;
    }

    MessageRemoveFile( sOString(oAssocObject(aaArgs[0])) );
    VP0("Messages will be displayed from the files (if any):\n" );
    MessageFileList();
    return NULL;
}



/*
 *      oCmd_debugStatus
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Display status of LEaP, like memory usage etc.
 *
 */
OBJEKT
oCmd_debugStatus( int iArgCount, ASSOC aaArgs[] )
{
char            *sCmd = "debugStatus";

    if ( !bCmdGoodArguments( sCmd,  iArgCount, aaArgs, "" ) ) {
        VPFATALDELAYEDEXIT( "usage:  debugStatus\n" );
        return NULL;
    }
#ifdef DEBUG
    VP0("DEBUG build\n");
    VP0("Messages will be displayed from the files (if any):\n" );
    MessageFileList();
#else
    VP0("Release build, no debugging\n");
#endif
    PrintMemoryStats();
    return NULL;
}


/*
 *      oCmd_help
 *
 *      Author: Christian Schafmeister (1991)
 *      Revised by:  David A. Rivkin (1992)
 *
 *      Offer help.
 */
OBJEKT
oCmd_help( int iArgCount, ASSOC aaArgs[] )
{
int             i = 1;
int             iColumns;
HELP            hTemp;

    iColumns = 4;

    if ( iArgCount == 0 ) {
        VP0("Help is available on the following subjects:\n\n" );
        HelpLoop();
        while ( (hTemp = hHelpNext()) ) {
            if ( i % iColumns ) {
                /* columns 0-(end-1) */
                VP0("%-22s", sHelpUpSubject(hTemp) );
            } else {
                /* last column in line */
                VP0("%s\n", sHelpUpSubject(hTemp) );
            }
            i++;
        }
        if ( (i-1) % iColumns )
            VP0("\n" );      /* close final part line */

        i = 0;
        VP0("\nFor a list of the current aliases, type \"alias\".\n" );
    } else if ( !bCmdGoodArguments( "help", iArgCount, aaArgs, "s" )) {
            VP0("usage:  help <command>\n" );
    } else {
        hTemp = hHelp( sOString(oAssocObject(aaArgs[0])));
        if ( hTemp == NULL ) {
                VP0("No help available.\n" );
        } else {
                VP0("\n" );
                VP0("%s\n", sHelpText(hTemp) );
        }
    }
    return NULL;
}

#ifdef DEBUG

OBJEKT
oCmd_lex( int iArgCount, ASSOC aaArgs[] )
{
    VP0("num args = %d\n",iArgCount);
    for (int i=0;i<iArgCount;i++) {
       VP0("Arg %d name=%s, desc=",i,sAssocName(aaArgs[i]));
        Describe( oAssocObject(aaArgs[i]) );
    }
    return NULL;
}

#endif

/*
 *      oCmd_list
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      List all variables declared.
 */
OBJEKT
oCmd_list( int iArgCount, ASSOC aaArgs[] )
{
    if ( !bCmdGoodArguments( "list", iArgCount, aaArgs, "" ) )
        VPFATALDELAYEDEXIT("usage:  list\n" );
    else
        VariablesList();
    return NULL;
}





/*
 *      oCmd_loadOff
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Load an OFF file and define all of the variables.
 */
OBJEKT
oCmd_loadOff( int iArgCount, ASSOC aaArgs[] )
{
STRING          sFilename;
char            *sNext;
OBJEKT          oObj;
LIBRARY         ul;

    if ( !bCmdGoodArguments( "loadOff", iArgCount, aaArgs, "s" ) ) {
        VPFATALDELAYEDEXIT("usage:  loadOff <filename>\n" );
        return NULL;
    }
    strcpy( sFilename, sOString(oAssocObject(aaArgs[0])) );

    ul = lLibraryOpen( sFilename, OPENREADONLY );
    if ( ul == NULL )
        return NULL;

    VP0("Loading library: %s\n", GsBasicsFullName );

    LibraryLoop( ul );
    do {
        sNext = sLibraryNext( ul );
        if ( sNext != NULL ) {
            VP1("Loading: %s\n", sNext );
            oObj = oLibraryLoad( ul, sNext );   /* comes w/ 1 REF */
            VariableSet( sNext, oObj );         /* adds 1 REF */
            DEREF( oObj );                      /* balance REF */

                        /* If the object loaded is a PARMSET then */
                        /* Add it to the PARMLIBRARY */
                        /* And make it the default PARMLIBRARY */

            if ( iObjectType(oObj)==PARMSETid ) {
                PARMSET psTemp = PARMSET_from(oObj);

                strcpy( sParmName(psTemp), sFilename );
                ParmLibAddParmSet( GplAllParameters, PARMSET_from(oObj) );
                ParmLibDefineDefault( GplAllParameters );
            }
        }
    } while ( sNext != NULL );

    LibraryClose( &ul );

    return NULL;
}





/*
 *      oCmd_sequence
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Create a molecule from a sequence of residues.
 *
 *      Arguments:
 *              [0]     A LIST of units.
 */
OBJEKT
oCmd_sequence( int iArgCount, ASSOC aaArgs[] )
{
LISTLOOP        llElements;
UNIT            uFirst, uSecond, uUnit;
ASSOC           aAss;
LOOP            lTemp, lAtoms;
ATOM            aConnect;
char            *sCmd = "sequence";
RESIDUE         rRes;
int             iDum;

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "l" ) ) {
        VPFATALDELAYEDEXIT("usage:  sequence <LIST>\n" );
        return NULL;
    }
    llElements = llListLoop( LIST_from(oAssocObject(aaArgs[0])) );

                /* Get the first element from the list */

    aAss = ASSOC_from(oListNext(&llElements));
    MESSAGE("Copying the first UNIT\n" );
    if ( iObjectType(oAssocObject(aAss)) != UNITid ) {
        VPFATALEXIT("%s: Invalid UNIT at position #1.\n", sCmd );
        return NULL;
    }

                /* Get the first UNIT and build INTERNALs */

    uFirst = uCopyUnit(oAssocObject(aAss));
    VP1("Sequence: %s\n", sContainerName(CONTAINER_from( uFirst)) );
    while ( (aAss = ASSOC_from(oListNext(&llElements))) ) {
        uUnit = UNIT_from(oAssocObject(aAss));
        if ( uUnit == NULL ) {
                VP0("Unknown UNIT: %s\n", sAssocName(aAss) );
        } else {

            MESSAGE("Copying a subsequent UNIT\n" );

                /* If the object is not a UNIT then destroy what we have */
                /* up to this point and return */

            if ( iObjectType(uUnit) != UNITid ) {
                Destroy((OBJEKT *)&uFirst);
                VPFATALEXIT("%s: Invalid UNIT named: %s\n", sCmd,
                        sAssocName(aAss) );
                return NULL;
            }

                        /* Copy the next UNIT */

            uSecond = uCopyUnit(uUnit);
            VP1("Sequence: %s\n", sContainerName(CONTAINER_from( uSecond)) );

                        /* Build INTERNALs for the next UNIT */

            MESSAGE("Building internals for subsequent UNIT\n" );
            BuildInternalsForContainer( CONTAINER_from( uSecond),
                        ATOMNEEDSBUILD, ATOMPOSITIONKNOWN );

            aConnect = aUnitHead( uSecond );

            if ( aConnect != NULL ) {
                BuildInternalsBetweenUnitsUsingFlags( uFirst, uSecond,
                                        ATOMNEEDSBUILD,
                                        0 );
            }
            MESSAGE("Joining two UNITS deleting the second\n" );
            UnitSequence( uFirst, uSecond );

            if ( aConnect != NULL ) {
                lAtoms = lLoop( OBJEKT_from(aConnect), SPANNINGTREE );
                iDum = 0;       /* for purify */
                BuildExternalsUsingFlags( &lAtoms, ATOMNEEDSBUILD, 0,
                                        ATOMPOSITIONKNOWN|ATOMPOSITIONBUILT,
                                        ATOMNEEDSBUILD|ATOMPOSITIONFIXED,
                                        &iDum, &iDum, &iDum, false );
            }
        }
    }

                /* Destroy INTERNALs to clean up */

    lAtoms = lLoop( OBJEKT_from(uFirst), ATOMS );
    BuildDestroyInternals( &lAtoms );

                /* Define PDB sequence */

    lTemp = lLoop( OBJEKT_from(uFirst), RESIDUES );
    while ( (rRes = RESIDUE_from(oNext(&lTemp))) ) {
        ResidueSetPdbSequence( rRes, iContainerSequence(CONTAINER_from( rRes)) );
    }

    return OBJEKT_from(uFirst);
}



/*
 *      oCmd_loadMol2
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Load a UNIT from a Mol2 file.
 *
 *      Arguments:
 *              [0]     - OSTRING, filename.
 */
OBJEKT
oCmd_loadMol2( int iArgCount, ASSOC aaArgs[] )
{
FILE            *fMol2;
UNIT            uUnit;

    if ( !bCmdGoodArguments( "loadMol2", iArgCount, aaArgs, "s" ) ) {
        VPFATALDELAYEDEXIT("usage:  <variable> = loadMol2 <filename>\n" );
        return NULL;
    }

    fMol2 = FOPENCOMPLAIN( sOString(oAssocObject(aaArgs[0])), "r" );
    if ( fMol2 == NULL ) return NULL;

    VP0("Loading Mol2 file: %s\n", GsBasicsFullName );
    uUnit = uTriposReadUnit( fMol2 );
    fclose(fMol2);

    return OBJEKT_from(uUnit);
}


/*___ oCmd_loadMol3 ___________________________________________________________
|                                                                              |
|       Author: Mason Louchart (2011)                                          |
|       http://q4md-forcefieldtools.org                                        |
|       Universite de Picardie - Jules Verne, Amiens                           |
|                                                                              |
|       Tutorial available at                                                  |
|       http://q4md-forcefieldtools.org/Tutorial/leap-mol3.php                 |
|                                                                              |
|       Load a UNIT from a Mol3 file.                                          |
|                                                                              |
|       Arguments:                                                             |
|               [0]     - OSTRING, filename.                                   |
|_____________________________________________________________________________*/
OBJEKT
oCmd_loadMol3( int iArgCount, ASSOC aaArgs[] )
{
FILE            *fMol3;
UNIT            uUnit;

    if ( !bCmdGoodArguments( "loadMol3", iArgCount, aaArgs, "s" ) ) {
        VPFATALDELAYEDEXIT("usage:  <variable> = loadMol3 <filename>\n" );
        return NULL;
    }

    fMol3 = FOPENCOMPLAIN( sOString(oAssocObject(aaArgs[0])), "r" );
    if ( fMol3 == NULL ) return NULL;

    VP0("Loading Mol3 file: %s\n", GsBasicsFullName );
    uUnit = uTriposReadUnit( fMol3 );

    fclose(fMol3);
    return OBJEKT_from(uUnit);
}


/*
 *      oCmd_loadCif
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Load a UNIT from a PDB file.
 *
 *      Arguments:
 *              [0]     - OSTRING, filename.
 */
OBJEKT
oCmd_loadCif( int iArgCount, ASSOC aaArgs[] )
{
FILE            *fCif;
UNIT            uUnit;
    if ( !bCmdGoodArguments( "loadCif", iArgCount, aaArgs, "s" ) ) {
        VPFATALDELAYEDEXIT("usage:  <variable> = loadCif <filename>\n" );
        return NULL;
    }

    fCif = FOPENCOMPLAIN( sOString(oAssocObject(aaArgs[0])), "r" );
    if ( fCif == NULL ) return NULL;

    VP0("Loading CIF file: %s\n", GsBasicsFullName );
    uUnit = uPdbRead( fCif, NULL, true );
    fclose(fCif);

    return OBJEKT_from(uUnit);
}

/*
 *      oCmd_loadPdb
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Load a UNIT from a PDB file.
 *
 *      Arguments:
 *              [0]     - OSTRING, filename.
 */
OBJEKT
oCmd_loadPdb( int iArgCount, ASSOC aaArgs[] )
{
FILE            *fPdb;
UNIT            uUnit;
char            *sFileName;
    if ( !bCmdGoodArguments( "loadPdb", iArgCount, aaArgs, "s" ) ) {
        VPFATALDELAYEDEXIT("usage:  <variable> = loadPdb <filename>\n" );
        return NULL;
    }

    sFileName = sOString(oAssocObject(aaArgs[0]));
    fPdb = FOPENCOMPLAIN( sFileName, "r" );
    if ( fPdb == NULL ) return NULL;

    VP0("Loading PDB file: %s\n", GsBasicsFullName );
    uUnit = uPdbRead( fPdb, NULL, false );
    fclose(fPdb);

    if (!GDefaults.bCompatible) {
    // Name UNIT based on filename
    size_t len = strlen(sFileName);
    if (len > 4 && !strcasecmp(&sFileName[len-4],".pdb")) sFileName[len-4]=0;
    char *p = strrchr(sFileName,'/');
    if (!p) p = sFileName;
    strcpy(sContainerName(uUnit),p);
    }

    return OBJEKT_from(uUnit);
}



/*
 *      oCmd_loadPdbUsingSeq
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Load a UNIT from a PDB file.
 *
 *      Arguments:
 *              [0]     - OSTRING, filename.
 *              [1]     - LIST, list of UNITs to use.
 */
OBJEKT
oCmd_loadPdbUsingSeq( int iArgCount, ASSOC aaArgs[] )
{
FILE            *fPdb;
UNIT            uUnit;
VARARRAY        vaUnits;
LIST            lUnits;
int             i, iErr;
ASSOC           aObj;
OBJEKT          oObj;
LISTLOOP        llLoop;
char            *sCmd = "loadPdbUsingSeq";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "s l" ) ) {
        VPFATALDELAYEDEXIT("usage:  <variable> = loadPdbUsingSeq <filename> <unitLIST>\n" );
        return NULL;
    }

    fPdb = FOPENCOMPLAIN( sOString(oAssocObject(aaArgs[0])), "r" );
    if ( fPdb == NULL ) return NULL;

                /* Copy the list of UNITs into a VARARRAY */

    lUnits = LIST_from(oAssocObject(aaArgs[1]));
    vaUnits = vaVarArrayCreate( sizeof(UNIT) );
    VarArraySetSize( vaUnits, iCollectionSize(lUnits) );
    llLoop = llListLoop(lUnits);
    i = 0;
    iErr = 0;
    while ( (aObj = ASSOC_from(oListNext(&llLoop))) ) {
        oObj = oAssocObject(aObj);
        if ( iObjectType(oObj) != UNITid ) {
            VP0("%s: %s is not a unit!\n", sCmd, sAssocName(aObj) );
            iErr++;
        }
        *PVAI( vaUnits, UNIT, i ) = UNIT_from(oObj);
        i++;
    }
    if ( iErr ) {
        VarArrayDestroy( &vaUnits );
        fclose(fPdb);
        VP0("Not loaded\n" );
        return NULL;
    }


    VP0("Loading PDB file: %s using sequence %s\n",
                GsBasicsFullName, sAssocName(aaArgs[1]) );
    uUnit = uPdbRead( fPdb, vaUnits, false);

    VarArrayDestroy( &vaUnits );

    fclose(fPdb);

    return OBJEKT_from(uUnit);
}






/*
 *      oCmd_saveOff
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Save a UNIT/PARMSET or a LIST of UNITs/PARMSETs to a UNITLIBRARY.
 *
 *      Arguments:
 *              [0]     - UNIT or LIST of UNITs to save.
 *              [1]     - OSTRING filename.
 */
OBJEKT
oCmd_saveOff( int iArgCount, ASSOC aaArgs[] )
{
LIBRARY         ul;
LISTLOOP        llUnits;
OBJEKT          oObj;
ASSOC           aObj;
char            *sCmd = "saveOff";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "upl s" ) ) {
        VPFATALDELAYEDEXIT("usage:  saveOff <object> <filename>\n" );
        return NULL;
    }

    ul = lLibraryOpen( sOString(oAssocObject(aaArgs[1])), OPENREADWRITE );
    if ( ul==NULL ) return NULL;

    DisplayerAccumulateUpdates();
    if ( iObjectType( oAssocObject(aaArgs[0]) ) == UNITid ||
         iObjectType( oAssocObject(aaArgs[0]) ) == PARMSETid ) {
        VP1("Saving %s.\n", sAssocName(aaArgs[0]) );
        LibrarySave( ul, sAssocName(aaArgs[0]),
                        oAssocObject(aaArgs[0]), NULL );
    } else {
        oObj = oAssocObject(aaArgs[0]);
        llUnits = llListLoop( LIST_from(oObj));
        while ( (aObj = ASSOC_from(oListNext(&llUnits)) ) != NULL ) {
            oObj = oAssocObject(aObj);
            if ( iObjectType(oObj) != UNITid &&
                 iObjectType(oObj) != PARMSETid ) {
                VP0("%s: Cannot save %s - type %s (ignoring).\n",
                        sCmd, sAssocName(aObj), sObjectType(oObj) );
            } else {
                VP1("Saving %s.\n", sAssocName(aObj) );
                LibrarySave( ul, sAssocName(aObj),
                                oAssocObject(aObj), NULL );
            }
        }
    }
    DisplayerReleaseUpdates();

    LibraryClose( &ul );
    return NULL;
}


/*
 *      Davids Changes
 *
 */

/*
 *      oCmd_createParmset
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return a newly created PARMSET.
 *
 *      Arguments:
 *              [0]     - OSTRING PARMSET name.
 */
OBJEKT
oCmd_createParmset( int iArgCount, ASSOC aaArgs[] )
{
PARMSET         psParmSet;

    if ( !bCmdGoodArguments( "createParmset",  iArgCount, aaArgs, "s" ) ) {
        VPFATALDELAYEDEXIT("usage:  <variable> = createParmset <name>\n" );
        return NULL;
    }

    psParmSet = PARMSET_from(oCreate(PARMSETid));
    AssocSetObject( aaArgs[0], OBJEKT_from(psParmSet ));
    return OBJEKT_from(psParmSet);
}




/*
 *      oCmd_createUnit
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return a newly created UNIT.
 *
 *      Arguments:
 *              [0]     - OSTRING UNIT name.
 */
OBJEKT
oCmd_createUnit( int iArgCount, ASSOC aaArgs[] )
{
UNIT            uUnit;

    if ( !bCmdGoodArguments( "createUnit", iArgCount, aaArgs, "s" ) ) {
        VPFATALDELAYEDEXIT("usage:  <variable> = createUnit <name>\n" );
        return NULL;
    }

    uUnit = UNIT_from(oCreate(UNITid));
    ContainerSetName( CONTAINER_from( uUnit), sOString(oAssocObject(aaArgs[0])) );
    return OBJEKT_from(uUnit);
}


/*
 *      oCmd_createResidue
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return a newly created RESIDUE.
 *
 *      Arguments:
 *              [0]     - OSTRING RESIDUE name.
 */
OBJEKT
oCmd_createResidue( int iArgCount, ASSOC aaArgs[] )
{
RESIDUE         rRes;

    if ( !bCmdGoodArguments( "createResidue", iArgCount, aaArgs, "s" ) ) {
        VPFATALDELAYEDEXIT("usage:  <variable> = createResidue <name>\n" );
        return NULL;
    }

    rRes = RESIDUE_from(oCreate(RESIDUEid));
    ContainerSetName( CONTAINER_from( rRes), sOString(oAssocObject(aaArgs[0])) );
    return OBJEKT_from(rRes);
}






/*
 *      oCmd_createAtom
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return a newly created ATOM.
 *
 *      Arguments:
 *              [0]     - OSTRING ATOM name.
 *              [1]     - OSTRING ATOM type.
 *              [2]     - ODOUBLE ATOM charge.
 */
OBJEKT
oCmd_createAtom( int iArgCount, ASSOC aaArgs[] )
{
ATOM            aAtom;
char            *sCmd = "createAtom";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "s s n" ) ) {
        VPFATALDELAYEDEXIT("usage:  <variable> = createAtom <name> <type> <charge>\n" );
        return NULL;
    }

    aAtom = ATOM_from(oCreate(ATOMid));
    ContainerSetName( CONTAINER_from( aAtom), sOString(oAssocObject(aaArgs[0])) );
    AtomSetType( aAtom, sOString(oAssocObject(aaArgs[1])) );
    AtomSetCharge( aAtom, dODouble(oAssocObject(aaArgs[2])) );

    return OBJEKT_from(aAtom);
}





/*
 *      oCmd_add
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Add the oB to oA.
 *      This can only be done if oB is not contained by anything else.
 *
 *      Arguments:
 *              [0]     - OBJEKT oA
 *              [1]     - OBJEKT oB
 */
OBJEKT
oCmd_add( int iArgCount, ASSOC aaArgs[] )
{
OBJEKT          oA, oB;
int             iA, iB;
STRING          sTemp;
char            *sCmd = "add";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "umr mra" ) ) {
        VPFATALDELAYEDEXIT("usage:  add <unit/residue/atom> <unit/residue/atom>\n" );
        return NULL;
    }

    oA = oAssocObject(aaArgs[0]);
    oB = oAssocObject(aaArgs[1]);
    iA = iObjectType(oA);
    iB = iObjectType(oB);
    if ( (iA==UNITid) &&
        !((iB==MOLECULEid)||(iB==RESIDUEid)||(iB==ATOMid)) ) {
        VP0("%s: UNITs cannot contain UNITs.\n", sCmd );
        return NULL;
    }
    if ( (iA==MOLECULEid) &&
        !((iB==RESIDUEid)||(iB==ATOMid)) ) {
        VP0("%s: MOLECULEs can only contain RESIDUEs and ATOMs.\n", sCmd );
        return NULL;
    }
    if ( (iA==RESIDUEid) && (iB!=ATOMid) ) {
        VP0("%s: Residues can only contain ATOMs.\n", sCmd );
        return NULL;
    }
    if ( cContainerWithin(CONTAINER_from( oB)) != NULL ) {
        VP0("%s: The object %s is already contained within %s\n",
                sCmd, sAssocName(aaArgs[1]),
                sContainerFullDescriptor( cContainerWithin(CONTAINER_from( oB)), sTemp ) );
        return NULL;
    }
    ContainerAdd( CONTAINER_from(oA), oB );
    return NULL;
}




/*
 *      oCmd_remove
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Remove oB from oA if it is contained by oA.
 *
 *      Arguments:
 *              [0]     - OBJEKT oA
 *              [1]     - OBJEKT oB
 */
OBJEKT
oCmd_remove( int iArgCount, ASSOC aaArgs[] )
{
OBJEKT          oA, oB;
char            *sCmd = "remove";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "umr mra" ) ) {
        VPFATALDELAYEDEXIT("usage:  remove <unit/residue/atom> <unit/residue/atom>\n" );
        return NULL;
    }

    oA = oAssocObject(aaArgs[0]);
    oB = oAssocObject(aaArgs[1]);

    REF( oB );  /* bContainerRemove() needs this */
    if ( !bContainerRemove( CONTAINER_from(oA), oB ) ) {
        VP0("%s: Could not find %s within %s.\n",
                sCmd, sAssocName(aaArgs[1]), sAssocName(aaArgs[0]) );
    }
    DEREF( oB ); /* reset count after bContainerRemove */

    return NULL;
}



/*
 *      oCmd_bond
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Create a bond between two atoms of the appropriate order.
 *      If the bond order is not given then use a single bond.
 *
 *      Arguments:
 *              [0]     - ATOM aA
 *              [1]     - ATOM aB
 *      option  [2]     - OSTRING, bond order.
 */
OBJEKT
oCmd_bond( int iArgCount, ASSOC aaArgs[] )
{
ATOM            aA, aB;
int             iOrder;
char            cOrder;
STRING          sTemp;
char            *sCmd = "bond";

   if ( iArgCount == 2 ) {
        cOrder = 'S';
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "a a" ) ) {
            VPFATALDELAYEDEXIT("usage:  bond <atom1> <atom2> [order]\n" );
            return NULL;
        }
    } else {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "a a s" ) ) {
            VPFATALDELAYEDEXIT("usage:  bond <atom1> <atom2> [order]\n" );
            return NULL;
        }
        cOrder = sOString(oAssocObject(aaArgs[2]))[0];
    }

    DisplayerAccumulateUpdates();

    aA = ATOM_from(oAssocObject(aaArgs[0]));
    aB = ATOM_from(oAssocObject(aaArgs[1]));

    switch ( cOrder ) {
        case 'S':
            iOrder = BONDSINGLE;
            break;
        case 'D':
            iOrder = BONDDOUBLE;
            break;
        case 'T':
            iOrder = BONDTRIPLE;
            break;
        case 'A':
            iOrder = BONDAROMATIC;
            break;
        default:
            VPFATALEXIT("%s: Unknown bond order (%c), no bond made.\n"
                    "Valid bond orders are: S, D, T, A.\n", sCmd, cOrder );
            goto DONE;
    }
    extern bool bIsConnectAtom(ATOM aAtom);
    if (iVerbosity()>2) {
        if (!bIsConnectAtom(aA))
            VPWARN("Bond: %s is not a CONNECT atom\n",
                        sContainerFullDescriptor(CONTAINER_from(aA), sTemp) );
        if (!bIsConnectAtom(aB))
            VPWARN("Bond: %s is not a CONNECT atom\n",
                        sContainerFullDescriptor(CONTAINER_from(aB), sTemp ) );
    }
    AtomBondToOrder( aA, aB, iOrder );
DONE:
    DisplayerReleaseUpdates();
    return NULL;
}


/*
 *      oCmd_deleteBond
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Remove a bond between two atoms.
 *      Print a message if the bond does not exist.
 *
 *      Arguments:
 *              [0]     - ATOM aA
 *              [1]     - ATOM aB
 */
OBJEKT
oCmd_deleteBond( int iArgCount, ASSOC aaArgs[] )
{
ATOM            aA, aB;
char            *sCmd = "deleteBond";

                /* Check the arguments */

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "a a" ) ) {
        VPFATALDELAYEDEXIT("usage:  deleteBond <atom1> <atom2>\n" );
        return NULL;
    }

    DisplayerAccumulateUpdates();

    aA = ATOM_from(oAssocObject(aaArgs[0]));
    aB = ATOM_from(oAssocObject(aaArgs[1]));

    if ( bAtomBondedTo( aA, aB ) ) {
        AtomRemoveBond( aA, aB );
    } else
        VP0("%s: Atoms are not bonded.\n", sCmd );

    DisplayerReleaseUpdates();

    return NULL;
}



/*
 *      oCmd_zMatrix
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Construct the external coordinates for the atoms.
 *
 *      Arguments:
 *              [0]     - Container that contains the atoms.
 *              [1]     - A list of atoms and internal coordinates.
 *
 *      The entries in the list of atoms and internal coordinates can
 *      look like:
 *
 *      a1 a2 b12
 *      a1 a2 a3 b12 t123
 *      a1 a2 a3 a4 b12 t123 p1234
 *      a1 a2 a3 a4 b12 t123 t124 orientation
 *
 *      Where a1,a2,a3,a4 can be an atom or an atom name which exists
 *      in the container.
 */
OBJEKT
oCmd_zMatrix( int iArgCount, ASSOC aaArgs[] )
{
LISTLOOP        llLines, llElements;
LIST            lLine;
int             iCount;
OBJEKT          oaElements[10];
VECTOR          vNew, vAtom2, vAtom3, vAtom4;
CONTAINER       cCont;
OBJEKT          oObj;
ASSOC           aAss, aAss2;
char            *sCmd = "zMatrix";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "umr l" ) ) {
        VPFATALDELAYEDEXIT("usage:  zMatrix <unit/residue> <LIST>\n" );
        return NULL;
    }

    DisplayerAccumulateUpdates();

    cCont = CONTAINER_from(oAssocObject(aaArgs[0]));
    llLines = llListLoop( LIST_from(oAssocObject(aaArgs[1])) );
    while ( (aAss = ASSOC_from(oListNext(&llLines)) ) != NULL ) {
        lLine = LIST_from(oAssocObject(aAss));
        if ( iObjectType(lLine) != LISTid ) {
            VPWARN("%s: Invalid object in zMatrix list was ignored.\n", sCmd );
            continue;
        }
        llElements = llListLoop( lLine );
        iCount = 0;
        while ( ( aAss2 = ASSOC_from(oListNext(&llElements))) != NULL ) {
            oObj = oAssocObject(aAss2);
            if ( iObjectType(oObj) == OSTRINGid ) {
                oObj = OBJEKT_from(cContainerFindName( cCont,
                                ATOMid, sOString(oObj)) );
            }
            oaElements[iCount] = oObj;
            iCount++;
        }

        switch ( iCount ) {
            case 3:
                if ( !bCmdGoodArguments( sCmd, iCount, (ASSOC *)oaElements,
                                                        "a a n" ) )
                    goto BAD;
                ZMatrixNothing( &vAtom2 );
                AtomSetPosition( oaElements[1], vAtom2 );
                ZMatrixBond( &vNew, &vAtom2, dODouble(oaElements[2]) );
                AtomSetPosition( oaElements[0], vNew );
                break;
            case 5:
                if ( !bCmdGoodArguments( sCmd, iCount, (ASSOC *)oaElements,
                                                        "a a a n n"))
                    goto BAD;
                vAtom2 = vAtomPosition( oaElements[1] );
                vAtom3 = vAtomPosition( oaElements[2] );
                ZMatrixBondAngle( &vNew, &vAtom2, &vAtom3,
                                dODouble(oaElements[3]),
                                dODouble(oaElements[4])*DEGTORAD );
                AtomSetPosition( oaElements[0], vNew );
                break;
            case 7:
                if ( !bCmdGoodArguments( sCmd, iCount, (ASSOC *)oaElements,
                                                "a a a a n n n" ) )
                    goto BAD;
                vAtom2 = vAtomPosition( oaElements[1] );
                vAtom3 = vAtomPosition( oaElements[2] );
                vAtom4 = vAtomPosition( oaElements[3] );
                ZMatrixBondAngleTorsion( &vNew, &vAtom2, &vAtom3, &vAtom4,
                                dODouble(oaElements[4]),
                                dODouble(oaElements[5])*DEGTORAD,
                                dODouble(oaElements[6])*DEGTORAD );
                AtomSetPosition( oaElements[0], vNew );
                break;
            case 8:
                if (!bCmdGoodArguments( sCmd, iCount, (ASSOC *)oaElements,
                                                "a a a a n n n n" ))
                    goto BAD;
                vAtom2 = vAtomPosition( oaElements[1] );
                vAtom3 = vAtomPosition( oaElements[2] );
                vAtom4 = vAtomPosition( oaElements[3] );
                ZMatrixBondTwoAnglesOrientation( &vNew,
                                &vAtom2, &vAtom3, &vAtom4,
                                dODouble(oaElements[4]),
                                dODouble(oaElements[5])*DEGTORAD,
                                dODouble(oaElements[6])*DEGTORAD,
                                dODouble(oaElements[7]) );
                AtomSetPosition( oaElements[0], vNew );
                break;
            default:
                VPWARN("%s: Invalid zMatrix entry was ignored.\n", sCmd );
                goto DONE;
        }
        continue;
BAD:
        VPWARN("%s: Invalid object in zMatrix entry.  Entry was ignored\n",
                                                sCmd );
        continue;
    }

DONE:
    DisplayerReleaseUpdates();
    return NULL;
}



/*
 *      oCmd_saveOffParm
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Save a UNIT or a LIST of UNITs to a UNITLIBRARY.
 *      Save it WITH parameters!!!!!!!!!!
 *
 *      Arguments:
 *              [0]     - UNIT or LIST of UNITs to save.
 *              [1]     - OSTRING filename.
 */
OBJEKT
oCmd_saveOffParm( int iArgCount, ASSOC aaArgs[] )
{
LIBRARY         ul;
LISTLOOP        llUnits;
UNIT            uUnit;
char            *sCmd = "saveOffParm";

    VP0("saveOffParm: command deactivated\n" );
    return NULL;
    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "ul s" ) ) {
        VPFATALDELAYEDEXIT("usage:  saveOffParm <object> <filename>\n" );
        return NULL;
    }

    if ( iParmLibSize(GplAllParameters) == 0 ) {
        VPWARN("%s: There are no parameter sets loaded\n", sCmd );
        return NULL;
    }

    ul = lLibraryOpen( sOString( oAssocObject(aaArgs[1]) ), OPENREADWRITE );

    if ( iObjectType( oAssocObject(aaArgs[0]) ) == UNITid ) {
        LibrarySave( ul, sContainerName( oAssocObject(aaArgs[0])),
                        oAssocObject(aaArgs[0]),
                        GplAllParameters );
    } else {
        llUnits = llListLoop( LIST_from(oAssocObject(aaArgs[0])) );
        while ( (uUnit = UNIT_from(oListNext(&llUnits)) ) != NULL ) {
            if ( iObjectType(uUnit) != UNITid ) {
                VPWARN("%s: Invalid UNIT in list was ignored.\n", sCmd );
                continue;
            }
            LibrarySave( ul, sContainerName(uUnit), OBJEKT_from( uUnit), GplAllParameters );
        }
    }

    LibraryClose( &ul );
    return NULL;
}




/*
 *      oCmd_loadAmberPrep
 *
 *      Author: Christian Schafmeister (1991)
 *      Modified: David A. Rivkin ( 1-15-93 ) Fixed the Dummy atom reading
 *
 *      Load an AMBER PREP file, add all of the UNITs into
 *      the systems Variable list.
 *
 *      Arguments:
 *              [0] -   OSTRING, filename.
 *      option  [1] -   OSTRING, string to prefix to each unit name.
 *
 *      The option [1] is provided to change the names of Old AMBER
 *      types to distinguish All Atom residues from United Atom residues
 *      and to distinguish Terminating residues from main chain
 *      residues.
 *
 *      FIXME: This leaks a dictionary with every invokation
 */
OBJEKT
oCmd_loadAmberPrep( int iArgCount, ASSOC aaArgs[] )
{
DICTIONARY      dUnits;
DICTLOOP        dlLoop;
UNIT            uUnit;

STRING          sName, sPrefix;
LOOP            lResidues;
RESIDUE         rRes;
char            *sCmd = "loadAmberPrep";

    if ( iArgCount == 1 ) {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "s" ) ) {
            VPFATALDELAYEDEXIT("usage:  loadAmberPrep <filename> [prefix]\n" );
            return NULL;
        }
        strcpy( sPrefix, "" );
    } else {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "s s" ) ) {
            VPFATALDELAYEDEXIT("usage:  loadAmberPrep <filename> [prefix]\n" );
            return NULL;
        }
        strcpy( sPrefix, sOString(oAssocObject(aaArgs[1])) );
    }
    dUnits = dAmberReadPrepFile( sOString(oAssocObject(aaArgs[0])) );
    if ( dUnits != NULL ) {
        dlLoop = ydlDictionaryLoop( dUnits );
        while ( (uUnit=UNIT_from(yPDictionaryNext( dUnits, &dlLoop ))) != NULL ) {
            strcpy( sName, sPrefix );
            strcat( sName, sContainerName(CONTAINER_from( uUnit)) );
            VP1("Loaded UNIT: %s\n", sName );

                /* Set the name for the UNIT */

            ContainerSetName( CONTAINER_from( uUnit), sName );

                /* Set the name for the only residue in the unit */

            lResidues = lLoop( OBJEKT_from(uUnit), RESIDUES );
            rRes = RESIDUE_from(oNext(&lResidues));
            ContainerSetName( CONTAINER_from( rRes), sName );
            VariableSet( sName, OBJEKT_from( uUnit ));       /* adds 1 REF */
            DEREF(uUnit);
        }
        DictionaryDestroy(&dUnits);
    }
    return NULL;
}

/*
 *      oCmd_loadAmberParams
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Load an AMBER parameter file, return the PARMSET
 *      and add the PARMSET to the systems PARMLIBRARY.
 *
 *      Arguments:
 *              [0] -   OSTRING, filename.
 *
 *      TODO: fix the ugly hack of checking a subset of already loaded libs
 *      Better: disable double load of all files? If so, by name or full path?
 */
OBJEKT
oCmd_loadAmberParams( int iArgCount, ASSOC aaArgs[] )
{
PARMSET         psParms;
STRING          sFile;
FILE            *fIn;
char            *sUsage =
                  "usage:  <variable> = loadAmberParams <filename> \n";


/*
 *   nasty, special kludge to prevent loading parm10 more than once;
 *      This is not generalizable....
 */
static          int parm10_loaded = 0;
static          int parm99_loaded = 0;
static          int parm15_loaded = 0;
static          char parm10[] = "parm10.dat";
static          char parm99[] = "parm99.dat";
static          char parm15[] = "parm15";

    if ( !bCmdGoodArguments( "loadAmberParams", iArgCount, aaArgs, "s" ) ) {
        VPFATALDELAYEDEXIT( sUsage );
        return NULL;
    }
    strcpy( sFile, sOString(oAssocObject(aaArgs[0])) );

    if( strstr( sFile, parm10 ) ) {
        parm10_loaded += 1;
        if( parm10_loaded > 1 ){
            VPNOTE( "Skipping %s: already loaded\n", sFile );
            parm10_loaded -= 1;
            return NULL;
        }
    }
    if( strstr( sFile, parm99 ) ) {
        parm99_loaded += 1;
        if( parm99_loaded > 1 ){
            VPNOTE( "Skipping %s: already loaded\n", sFile );
            parm99_loaded -= 1;
            return NULL;
        }
    }
    if( strstr( sFile, parm15 ) ) {
        parm15_loaded += 1;
        if( parm15_loaded > 1 ){
            VPNOTE( "Skipping %s: already loaded\n", sFile );
            parm15_loaded -= 1;
            return NULL;
        }
    }

    fIn = FOPENCOMPLAIN( sFile, "r" );
    if ( fIn == NULL )
        return NULL;

    if( parm99_loaded + parm15_loaded + parm10_loaded > 1 ){
        VPFATALEXIT( "Cannot load more than one of parm99/10/15.dat\n"
                "If you are running interactively then you should save your work,"
                "\nquit LEaP, and retry being careful to load only one parm of"
                " above.\n" );
    }

    VP0("Loading parameters: %s\n", GsBasicsFullName );
    psParms = psAmberReadParmSet( fIn, sFile );

    if ( psParms != NULL ) {
        ParmLibAddParmSet( GplAllParameters, psParms );
        ParmLibDefineDefault( GplAllParameters );
    } else
        VPNOTE( "-- no parameters loaded" );

    return OBJEKT_from(psParms);
}


/*
 *      oCmd_savePdb
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Save the UNIT as a PDB file.
 *      Create a duplicate of the UNIT and impose a hierarchy on
 *      it.
 *
 *      Arguments:
 *              [0] -   Unit to save.
 *              [1] -   Name of file to save to.
 */
OBJEKT
oCmd_savePdb( int iArgCount, ASSOC aaArgs[] )
{
char            *sString;
FILE            *fOut;

    if ( !bCmdGoodArguments( "savePdb", iArgCount, aaArgs, "u s" ) ) {
        VPFATALDELAYEDEXIT("usage:  savePdb <object> <filename>\n" );
        return NULL;
    }
    sString = sOString( oAssocObject(aaArgs[1]) );
    fOut = FOPENCOMPLAIN( sString, "w" );
    if ( fOut == NULL )
        return NULL;
    VP0("Writing pdb file: %s\n", sString);
    PdbWrite( fOut, UNIT_from(oAssocObject(aaArgs[0]) ));
    fclose( fOut );
    return NULL;
}


/*
 *      oCmd_saveMol2
 *      Based on savepdb
 *      Author: Christine Cezard (2007)
 *      Universite de Picardie - Jules Verne, Amiens
 *      http://q4md-forcefieldtools.org
 *
 *      Tutorial available at
 *      http://q4md-forcefieldtools.org/Tutorial/leap.php
 *
 *      Save the UNIT as a Mol2 file.
 *      Create a duplicate of the UNIT and impose a hierarchy on it.
 *
 *      Arguments:
 *              [0] -   Unit to save.
 *              [1] -   Name of file to save to.
 *              [2] -   Option for column 6 (0 = Default, 1 = Amber Atom type)
 */
OBJEKT
oCmd_saveMol2( int iArgCount, ASSOC aaArgs[] )
{
char            *sString;
FILE            *fOut;
double          choice;

    if ( !bCmdGoodArguments( "saveMol2", iArgCount, aaArgs, "u s n" ) ) {
        VP0("usage:  saveMol2 <object> <filename> <option>\n" );
        VP0("<option> = 0 for Default atom types \n" );
        VPFATALDELAYEDEXIT("<option> = 1 for AMBER atom types \n" );
        return NULL;
    }
    sString = sOString( oAssocObject(aaArgs[1]) );
    choice =  (int) dODouble( oAssocObject(aaArgs[2]) ) ;
    fOut = FOPENCOMPLAIN( sString, "w" );
    if ( fOut == NULL )
        return NULL;
    VP0("Writing mol2 file: %s\n", sString);
    Mol2Write( fOut, UNIT_from(oAssocObject(aaArgs[0])), choice);
    fclose( fOut );
    return NULL;
}


/*_____ oCmd_saveMol3 _________________________________________________________
|                                                                              |
|       Based on saveMol2                                                      |
|       Author: Mason Louchart (2011)                                          |
|       http://q4md-forcefieldtools.org                                        |
|       Universite de Picardie - Jules Verne, Amiens                           |
|                                                                              |
|       Tutorial available at                                                  |
|       http://q4md-forcefieldtools.org/Tutorial/leap-mol3.php                 |
|                                                                              |
|       Save the UNIT as a Mol3 file.                                          |
|                                                                              |
|       Arguments:                                                             |
|               [0] -  Unit to save.                                           |
|               [1] -  Name of file to save to.                                |
|               [2] -  Option for column 6 (0 = Default, 1 = Amber Atom type)  |
|_____________________________________________________________________________*/
OBJEKT
oCmd_saveMol3( int iArgCount, ASSOC aaArgs[] )
{
char            *sString;
FILE            *fOut;
double          choice;

    if ( !bCmdGoodArguments( "saveMol3", iArgCount, aaArgs, "u s n" ) ) {
        VP0("usage:  saveMol3 <object> <filename> <option>\n" );
        VP0("<option> = 0 for Default atom types \n" );
        VPFATALDELAYEDEXIT("<option> = 1 for AMBER atom types \n" );
        return NULL;
    }
    sString = sOString( oAssocObject(aaArgs[1]) );
    choice =  (int) dODouble( oAssocObject(aaArgs[2]) );
    fOut = FOPENCOMPLAIN( sString, "w" );
    if ( fOut == NULL ) return NULL;

    VP0("Writing mol3 file: %s\n", sString);

    Mol3Write( fOut, UNIT_from(oAssocObject(aaArgs[0])), choice);
    fclose( fOut );
    return NULL;
}


/*
 *      oCmd_solvateBox
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Solvate the UNIT [0] within copies of a box of SOLVENT.
 *
 *      Arguments:
 *              [0] -   Unit to add solvent to.
 *              [1] -   Unit with solvent.
 *              [2] -   LIST with ( x, y, z ) distances or
 *                      xyz, the closest the wall of the solvent box
 *                      can come to the smallest box which contains
 *                      the entire unit that is centered on the
 *                      origin.  This is called the BufferZone.
 *                      Or a single ODOUBLE when x,y,z distances are the
 *                      same.
 *      Option  [3] -   ODOUBLE with closeness parameter.
 */
OBJEKT
oCmd_solvateBox( int iArgCount, ASSOC aaArgs[] )
{
UNIT            uSolvent, uSolute;
double          dCloseness = 1.0, daBuffer[3];
ASSOC           aAss;
LISTLOOP        llNumbers;
int             i, iInitialSize, iFinalSize;
bool            bIsotropic = false;
char            *sCmd = "solvateBox";
char            *sUsage =
                 "usage:  solvateBox <solute> <solvent> <buffer> [iso] [closeness]\n";

    /*
     *  check args - always need 2 units & cutoff
     */
    if ( !bCmdGoodArguments( sCmd, 3, aaArgs, "u u ln" ) ) {
        VPFATALDELAYEDEXIT("%s",sUsage );
        return NULL;
    }
    if ( iArgCount > 5 ) {
        VPFATALDELAYEDEXIT("%s",sUsage );
        return NULL;
    }
    if ( iArgCount != 3 ) {
        OBJEKT  oObj;

        /*
         *  handle possible 'iso' &/or closeness
         */
        oObj = OBJEKT_from(aaArgs[3]);
        if ( iObjectType(oObj) != ASSOCid )
                DFATAL("unexpected objtype %s\n", sObjectType(oObj) );

        if ( strcmp("iso", sAssocName(oObj) ) ) {
                oObj = oAssocObject(oObj);
        } else {
                bIsotropic = true;
                if ( iArgCount == 5 ) {
                        oObj = oAssocObject(aaArgs[4]);
                } else
                        oObj = NULL;
        }
        if ( oObj != NULL ) {
                if ( iObjectType(oObj) != ODOUBLEid ) {
                        VPFATALEXIT("%s",sUsage );
                        return NULL;
                }
                dCloseness = dODouble( oObj );
                if ( dCloseness <= 0.0 ) {
                        VPFATALEXIT("%s",sUsage );
                        return NULL;
                }
        }
    }

        /* Make sure there is a parameter set loaded */

    if ( iParmLibSize(GplAllParameters) == 0 ) {
        VPWARN("%s: There are no parameter sets loaded\n", sCmd );
        return NULL;
    }

    uSolute = UNIT_from(oAssocObject(aaArgs[0]));
    uSolvent = UNIT_from(oAssocObject(aaArgs[1]));

    if ( iObjectType( oAssocObject(aaArgs[2]) ) == LISTid ) {
        /*
         *  x,y,z box clearances
         */
        if ( iListSize(oAssocObject(aaArgs[2])) != 3 ) {
            VPFATALEXIT("%s: Argument #3 must have three values: "
                    "{ x y z } or one.\n%s", sCmd, sUsage );
            return NULL;
        }
        llNumbers = llListLoop(LIST_from(oAssocObject(aaArgs[2])));
        for ( i=0; i<3; i++ ) {
            aAss = ASSOC_from(oListNext(&llNumbers));
            if ( iObjectType(oAssocObject(aAss)) != ODOUBLEid ) {
                VPFATALEXIT("%s: Bad value #%d in the third argument.\n%s",
                        sCmd, i, sUsage );
                return NULL;
            }
            daBuffer[i] = dODouble(oAssocObject(aAss));
            if ( daBuffer[i] < 0.0 ) {
                VPFATALEXIT("%s: Bad value #%d in the third argument.\n%s",
                        sCmd, i+1, sUsage );
                return NULL;
            }
        }
    } else {
        daBuffer[0] = dODouble(oAssocObject(aaArgs[2]));
        if ( daBuffer[0] < 0.0 ) {
                VPFATALEXIT("%s: Bad box clearance.\n%s", sCmd, sUsage );
                return NULL;
        }
        daBuffer[2] = daBuffer[1] = daBuffer[0];
    }

    /*
     *  make copy of solvent w/ box & solv residue types set
     */
    uSolvent = zToolSetupSolvent( uSolvent );

    /*
     *  set up solute - centered, & if iso, w/ principal axes aligned
     */
    TurnOffDisplayerUpdates();

    if ( bIsotropic )   /* orient principal axes */
        ToolCenterUnitByRadii( uSolute, true );
    else
        ToolCenterUnitByRadii( uSolute, false );

    TurnOnDisplayerUpdates();
    ContainerDisplayerUpdate( CONTAINER_from( uSolute ));

    /*
     *  do the solvation
     */
    TurnOffDisplayerUpdates();

    iInitialSize = iContainerNumberOfChildren( CONTAINER_from( uSolute ));

    zToolSolvateAndShell( uSolute, uSolvent,
                daBuffer[0], daBuffer[1], daBuffer[2], dCloseness,
                NOSHELL, false, true, false, bIsotropic );

    /*
     *  Get rid of solvent copy
     */
    ContainerDestroy( (CONTAINER *) &uSolvent );
    iFinalSize = iContainerNumberOfChildren( CONTAINER_from( uSolute ));

    TurnOnDisplayerUpdates();
    ContainerDisplayerUpdate(CONTAINER_from( uSolute));


    VP0("  Added %d residues.\n", iFinalSize - iInitialSize );

    return NULL;
}


/*
 *      oCmd_solvateBox
 *
 *      Author: Juno Krahn (2026)
 *
 *      Solvate the UNIT [0] within copies of a box of SOLVENT.
 *
 *      Arguments:
 *              [0] -   Unit to add solvent to.
 *              [1] -   Unit with solvent.
 *      Option  [2] -   ODOUBLE with closeness parameter.
 */
OBJEKT
oCmd_solvateCell( int iArgCount, ASSOC aaArgs[] )
{
UNIT            uSolvent, uSolute;
double          dCloseness = 1.0;
int             iInitialSize, iFinalSize;
bool            bUsage;
char            *sCmd = "solvateCell";
char            *sUsage =
                 "usage:  \n";

    if (iArgCount == 2)
        bUsage = !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u" );
    else
        bUsage = !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u n" );
    if (bUsage) {
        VPFATALDELAYEDEXIT("%s",sUsage );
        return NULL;
    }

        /* Make sure there is a parameter set loaded */

    if ( iParmLibSize(GplAllParameters) == 0 ) {
        VPWARN("%s: There are no parameter sets loaded\n", sCmd );
        return NULL;
    }

    uSolute = UNIT_from(oAssocObject(aaArgs[0]));
    uSolvent = UNIT_from(oAssocObject(aaArgs[1]));
    if ( iArgCount == 3 )
        dCloseness = dODouble(oAssocObject(aaArgs[2]));

    /*
     *  make copy of solvent w/ box & solv residue types set
     */
    uSolvent = zToolSetupSolvent( uSolvent );

    ContainerDisplayerUpdate( CONTAINER_from( uSolute ));

    /*
     *  do the solvation
     */
    TurnOffDisplayerUpdates();

    iInitialSize = iContainerNumberOfChildren( CONTAINER_from( uSolute ));

    ToolSolvateCell( uSolute, uSolvent, dCloseness);

    /*
     *  Get rid of solvent copy
     */
    ContainerDestroy( (CONTAINER *) &uSolvent );
    iFinalSize = iContainerNumberOfChildren( CONTAINER_from( uSolute ));

    TurnOnDisplayerUpdates();
    ContainerDisplayerUpdate(CONTAINER_from( uSolute));


    VP0("  Added %d residues.\n", iFinalSize - iInitialSize );

    return NULL;
}



/*
 *      oCmd_solvateOct
 *
 *      Author: Bill Ross (1998)
 *
 *      Solvate the UNIT with copies of a box of SOLVENT in a
 *      box with shaved corners.
 *
 *      Arguments:
 *              [0] -   Unit to add solvent to.
 *              [1] -   Unit with solvent.
 *              [2] -   LIST with ( x, y, z, d ) distances or
 *                      Or a single ODOUBLE when x,y,z,d distances are the
 *                      same.
 *      Option  [3] -   ODOUBLE with closeness parameter.
 */
OBJEKT
oCmd_solvateOct( int iArgCount, ASSOC aaArgs[] )
{
OBJEKT          oObj;
UNIT            uSolvent, uSolute;
double          dCloseness = 1.0, daBuffer[4];
ASSOC           aAss;
LISTLOOP        llNumbers;
int             i, iInitialSize, iFinalSize;
bool            bIsotropic = true;
char            *sCmd = "solvateOct";
char            *sUsage =
                 "usage:  solvateOct <solute> <solvent> <buffer> [aniso] [closeness]\n";

        /*
         *  check args - always need 2 units & cutoff
         */
        if ( !bCmdGoodArguments( sCmd, 3, aaArgs, "u u ln" ) ) {
                VPFATALDELAYEDEXIT("%s",sUsage );
                return NULL;
        }
        if ( iArgCount > 5 ) {
                VPFATALDELAYEDEXIT("%s",sUsage );
                return NULL;
        }
        if ( iArgCount != 3 ) {

                /*
                *  handle possible 'aniso' &/or closeness
                */

                oObj = OBJEKT_from(aaArgs[3]);
                if ( iObjectType(oObj) != ASSOCid )
                DFATAL("unexpected objtype %s\n", sObjectType(oObj) );

                if ( strcmp("aniso", sAssocName(oObj) ) ) {

                        /* aniso not found;
                           check for old-fashioned iso keyword, and ignore it:  */
                        if ( strcmp("iso", sAssocName(oObj) ) ) {
                                oObj = oAssocObject(oObj);
                        } else {
                                if ( iArgCount == 5 ) {
                                        oObj = oAssocObject(aaArgs[4]);
                                } else {
                                        oObj = NULL;
                                }
                        }

                } else {

                        /* found aniso keyword: */
                        bIsotropic = false;
                        if ( iArgCount == 5 ) {
                                oObj = oAssocObject(aaArgs[4]);
                        } else {
                                oObj = NULL;
                        }
                }

                if ( oObj != NULL ) {
                        if ( iObjectType(oObj) != ODOUBLEid ) {
                                VPFATALEXIT("%s",sUsage );
                                return NULL;
                        }
                        dCloseness = dODouble( oObj );
                        if ( dCloseness <= 0.0 ) {
                                VPFATALEXIT("%s",sUsage );
                                return NULL;
                        }
                }
        }

        /* Make sure there is a parameter set loaded */

    if ( iParmLibSize(GplAllParameters) == 0 ) {
        VPWARN("%s: There are no parameter sets loaded\n", sCmd );
        return NULL;
    }

    uSolute = UNIT_from(oAssocObject(aaArgs[0]));
    uSolvent = UNIT_from(oAssocObject(aaArgs[1]));

    if ( iObjectType( oAssocObject(aaArgs[2]) ) == LISTid ) {

        if ( bIsotropic ) {
                VPFATALEXIT("%s: 'iso' requires a single clearance value\n", sCmd );
                return NULL;
        }

        /*
         *  x,y,z,d box clearances
         */
        if ( iListSize(oAssocObject(aaArgs[2])) != 4 ) {
            VPFATALEXIT("%s: Argument #3 must have four values: "
                    "{ x y z d } or one.\n%s", sCmd, sUsage );
            return NULL;
        }
        llNumbers = llListLoop(LIST_from(oAssocObject(aaArgs[2])));
        for ( i=0; i<4; i++ ) {
            aAss = ASSOC_from(oListNext(&llNumbers));
            if ( iObjectType(oAssocObject(aAss)) != ODOUBLEid ) {
                VPFATALEXIT("%s: Bad value #%d in the third argument.\n%s",
                        sCmd, i, sUsage );
                return NULL;
            }
            daBuffer[i] = dODouble(oAssocObject(aAss));
            if ( daBuffer[i] < 0.0 ) {
                VPFATALEXIT("%s: Bad value #%d in the third argument.\n%s",
                        sCmd, i+1, sUsage );
                return NULL;
            }
        }
    } else {
        daBuffer[0] = dODouble(oAssocObject(aaArgs[2]));
        if ( daBuffer[0] < 0.0 ) {
                VPFATALEXIT("%s: Bad box clearance.\n%s", sCmd, sUsage );
                return NULL;
        }
        daBuffer[3] = daBuffer[2] = daBuffer[1] = daBuffer[0];
    }

    /*
     *  make copy of solvent w/ box & solv residue types set
     */
    uSolvent = zToolSetupSolvent( uSolvent );

    /*
     *  center solute by vdw
     */
    TurnOffDisplayerUpdates();
    if ( bIsotropic )   /* orient principal axes */
        ToolCenterUnitByRadii( uSolute, true );
    else
        ToolCenterUnitByRadii( uSolute, false );
    TurnOnDisplayerUpdates();
    ContainerDisplayerUpdate( CONTAINER_from( uSolute ));

                /* adjust box for diagonal clearance if necc */

    ToolOctBoxCheck( uSolute, daBuffer, true, bIsotropic );

                /* Solvate the solute */

    TurnOffDisplayerUpdates();

    iInitialSize = iContainerNumberOfChildren( CONTAINER_from( uSolute ));
    zToolSolvateAndShell( uSolute, uSolvent,
                daBuffer[0], daBuffer[1], daBuffer[2],
                dCloseness, NOSHELL, false, true, true, bIsotropic );

    /*
     *  Get rid of solvent copy
     */
    ContainerDestroy( (CONTAINER *) &uSolvent );
    iFinalSize = iContainerNumberOfChildren( CONTAINER_from( uSolute ));

    TurnOnDisplayerUpdates();
    ContainerDisplayerUpdate(CONTAINER_from( uSolute));


    VP0("  Added %d residues.\n", iFinalSize - iInitialSize );

    return NULL;
}


/*
 *      oCmd_solvateDontClip
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Solvate the UNIT [0] within copies of a box of SOLVENT.
 *      Do not clip the large SOLVENT box that is created around
 *      the solute by duplicating the SOLVENT box.
 *
 *      Arguments:
 *              [0] -   Unit to add solvent to.
 *              [1] -   Unit with solvent.
 *              [2] -   LIST with ( x, y, z ) distances or
 *                      xyz, the closest the wall of the solvent box
 *                      can come to the smallest box which contains
 *                      the entire unit that is centered on the
 *                      origin.  This is called the BufferZone.
 *                      Or a single ODOUBLE when x,y,z distances are the
 *                      same.
 *      Option  [3] -   ODOUBLE with closeness parameter.
 */
OBJEKT
oCmd_solvateDontClip( int iArgCount, ASSOC aaArgs[] )
{
UNIT            uSolvent, uSolute;
double          dCloseness = 1.0, daBuffer[3];
ASSOC           aAss;
LISTLOOP        llNumbers;
int             i, iInitialSize, iFinalSize;
bool            bIsotropic = false;
char            *sCmd = "solvateDontClip";
char            *sUsage =
                 "usage:  solvateDontClip <solute> <solvent> <buffer> [closeness]\n";

    if ( iArgCount == 3 ) {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u ln" ) ) {
            VPFATALDELAYEDEXIT("%s",sUsage );
            return NULL;
        }
        dCloseness = 1.0;
    } else {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u ln n" ) ) {
            VPFATALDELAYEDEXIT("%s",sUsage );
            return NULL;
        }
        dCloseness = dODouble( oAssocObject(aaArgs[3]) );
    }

        /* Make sure there is a parameter set loaded */

    if ( iParmLibSize(GplAllParameters) == 0 ) {
        VPWARN("%s: There are no parameter sets loaded\n", sCmd );
        return NULL;
    }

    uSolute = UNIT_from(oAssocObject(aaArgs[0]));
    uSolvent = UNIT_from(oAssocObject(aaArgs[1]));

    if ( iObjectType( oAssocObject(aaArgs[2]) ) == LISTid ) {
        if ( iListSize(oAssocObject(aaArgs[2])) != 3 ) {
            VPFATALEXIT("%s: Argument #3 must have three values: "
                    "{ x y z } or one.\n%s", sCmd, sUsage );
            return NULL;
        }
        llNumbers = llListLoop(LIST_from(oAssocObject(aaArgs[2])));
        for ( i=0; i<3; i++ ) {
            aAss = ASSOC_from(oListNext(&llNumbers));
            if ( iObjectType(oAssocObject(aAss)) != ODOUBLEid ) {
                VPFATALEXIT("%s: Bad value #%d in the third argument.\n%s",
                        sCmd, i, sUsage );
                return NULL;
            } else daBuffer[i] = dODouble(oAssocObject(aAss));
        }
    } else {
        daBuffer[0] = dODouble(oAssocObject(aaArgs[2]));
        daBuffer[2] = daBuffer[1] = daBuffer[0];
    }

    /*
     *  make copy of solvent w/ box & solv residue types set
     */
    uSolvent = zToolSetupSolvent( uSolvent );

        /* Orient the principle axis of the solute along the */
        /* coordinate axis, setting a box  */

    TurnOffDisplayerUpdates();
    ToolCenterUnitByRadii( uSolute, false );
    TurnOnDisplayerUpdates();

    ContainerDisplayerUpdate(CONTAINER_from( uSolute));

    TurnOffDisplayerUpdates();
    iInitialSize = iContainerNumberOfChildren( CONTAINER_from( uSolute ));
    zToolSolvateAndShell( uSolute, uSolvent,
                daBuffer[0], daBuffer[1], daBuffer[2],
                dCloseness, NOSHELL, false, false, false, bIsotropic );

    /*
     *  Get rid of solvent copy
     */
    ContainerDestroy( (CONTAINER *) &uSolvent );
    iFinalSize = iContainerNumberOfChildren( CONTAINER_from( uSolute ));

    TurnOnDisplayerUpdates();
    ContainerDisplayerUpdate(CONTAINER_from( uSolute));

    VP0("  Added %d residues.\n", iFinalSize - iInitialSize );

    return NULL;
}





/*
 *      oCmd_solvateShell
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Solvate the UNIT [0] within copies of a box of SOLVENT.
 *
 *      Arguments:
 *              [0] -   Unit to add solvent to.
 *              [1] -   Unit with solvent.
 *              [2] -   ODOUBLE with thickness parameter.
 *      Option  [3] -   ODOUBLE with closeness parameter.
 */
OBJEKT
oCmd_solvateShell( int iArgCount, ASSOC aaArgs[] )
{
UNIT            uSolvent, uSolute;
double          dCloseness = 1.0, dThickness;
int             iInitialSize, iFinalSize;
char            *sCmd = "solvateShell";
char            *sUsage =
        "usage:  solvateShell <solute> <solvent> <buffer> [closeness]\n";

    if ( iArgCount == 3 ) {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u n" ) ) {
            VPFATALDELAYEDEXIT("%s",sUsage );
            return NULL;
        }
        dCloseness = 1.0;
    } else {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u n n" ) ) {
            VPFATALDELAYEDEXIT("%s",sUsage );
            return NULL;
        }
        dCloseness = dODouble( oAssocObject(aaArgs[3]) );
    }

        /* Make sure there is a parameter set loaded */

    if ( iParmLibSize(GplAllParameters) == 0 ) {
        VPWARN("%s: There are no parameter sets loaded\n", sCmd );
        return NULL;
    }

    uSolute = UNIT_from(oAssocObject(aaArgs[0]));
    uSolvent = UNIT_from(oAssocObject(aaArgs[1]));

    /*
     *  make copy of solvent w/ box & solv residue types set
     */
    uSolvent = zToolSetupSolvent( uSolvent );

        /* Orient the principle axis of the solute along the */
        /* coordinate axis, setting a box */

    TurnOffDisplayerUpdates();
    ToolCenterUnitByRadii( uSolute, false );
    TurnOnDisplayerUpdates();
    ContainerDisplayerUpdate(CONTAINER_from( uSolute));

    dThickness = dODouble(oAssocObject(aaArgs[2]));
    if ( dThickness <= 0.0 ) {
        VPFATALEXIT("%s: bad thickness\n%s", sCmd, sUsage );
        return  NULL ;
    }

    TurnOffDisplayerUpdates();
    iInitialSize = iContainerNumberOfChildren( CONTAINER_from( uSolute ));
    zToolSolvateAndShell( uSolute, uSolvent,
                dThickness, dThickness, dThickness,
                dCloseness, dThickness, true, true, false, false );

    /*
     *  Get rid of solvent copy
     */
    ContainerDestroy( (CONTAINER *) &uSolvent );
    iFinalSize = iContainerNumberOfChildren( CONTAINER_from( uSolute ));

    TurnOnDisplayerUpdates();
    ContainerDisplayerUpdate(CONTAINER_from( uSolute));

    VP0(" Added %d residues.\n", iFinalSize - iInitialSize );


    return NULL;
}








/*
 *      oCmd_verbosity
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Set the verbosity to the value passed.
 *
 *      Arguments:
 *              [0] -   Verbosity level 0-2.
 */
OBJEKT
oCmd_verbosity( int iArgCount, ASSOC aaArgs[] )
{
int     iVerb;

    if ( !bCmdGoodArguments( "verbosity", iArgCount, aaArgs, "n" ) ) {
        VPFATALDELAYEDEXIT("usage:  verbosity <level>\n" );
        return NULL;
    }

    iVerb = (int)dODouble(oAssocObject(aaArgs[0]));
    VerbositySet(iVerb);
    GrMainResult.iCommand     = CMD_VERBOSITY;
    GrMainResult.sVariable[0] = (char)iVerb;

    VP2("Verbosity level: %d\n", iVerb );
    return NULL;
}


//#include <stdio.h>
//#include <string.h>
#include <sys/stat.h>

#define MAX_LOG_SIZE   (10 * 1024 * 1024)   /* trigger threshold */
#define KEEP_LOG_SIZE  (5  * 1024 * 1024)   /* how much tail to keep */
// use getenv AMBERLEAP_LOGMAX
void
CheckAndTruncateLog( const char *sLogPath )
{
struct stat     statBuf;
FILE            *fpIn, *fpOut;
char            sTempPath[1024];
long            lSkip;
char            sBuf[65536];
size_t          nRead;

    if ( stat( sLogPath, &statBuf ) != 0 ) return;      /* doesn't exist yet */
    if ( statBuf.st_size <= MAX_LOG_SIZE ) return;       /* under threshold */

    lSkip = statBuf.st_size - KEEP_LOG_SIZE;

    fpIn = fopen( sLogPath, "r" );
    if ( fpIn == NULL ) return;

    if ( fseek( fpIn, lSkip, SEEK_SET ) != 0 ) {
        fclose(fpIn);
        return;
    }

		/* Optional: skip forward to the next newline so we */
		/* don't keep a truncated partial line at the top   */
    {
    int c;
        while ( (c = fgetc(fpIn)) != EOF && c != '\n' ) { }
    }

    snprintf( sTempPath, sizeof(sTempPath), "%s.tmp", sLogPath );
    fpOut = fopen( sTempPath, "w" );
    if ( fpOut == NULL ) {
        fclose(fpIn);
        return;
    }

    fprintf( fpOut, "[earlier log entries truncated]\n" );
    while ( (nRead = fread(sBuf, 1, sizeof(sBuf), fpIn)) > 0 ) {
        fwrite( sBuf, 1, nRead, fpOut );
    }

    fclose(fpIn);
    fclose(fpOut);

    rename( sTempPath, sLogPath );   /* atomic on POSIX, same filesystem */
}

/*
 *      oCmd_logFile
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Define the LogFile.
 *
 *      Arguments:
 *              [0] -   OSTRING with Log File name.
 */
OBJEKT
oCmd_logFile( int iArgCount, ASSOC aaArgs[] )
{
STRING          sFilename;
FILE            *fLog;

    if ( !bCmdGoodArguments( "logFile", iArgCount, aaArgs, "s" ) ) {
        VPFATALDELAYEDEXIT("usage:  logFile <filename>\n" );
        return NULL;
    }

    strcpy( sFilename, sOString(oAssocObject(aaArgs[0])) );

                /* Close the old log file if it exists */

    if ( fVerbosityLogFile() != NULL ) fclose( fVerbosityLogFile() );
    VerbositySetLogFile( NULL );

                /* Open the new log file */

    CheckAndTruncateLog( sFilename );
    fLog = FOPENCOMPLAIN( sFilename, "a" );
    if ( fLog != NULL ) {
        time_t  t = time( (time_t *)0 );

        fprintf(fLog, "log started: %s\n", ctime(&t) );
        VerbositySetLogFile( fLog );
        VP0("Log file: %s\n", GsBasicsFullName );
    }
    return NULL;
}




/*
 *      oCmd_combine
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Combine the contents of several UNITs.
 *      Return a UNIT that is the combined contents.
 *
 *      Arguments:
 *              [0]     A LIST of units.
 */
OBJEKT
oCmd_combine( int iArgCount, ASSOC aaArgs[] )
{
LISTLOOP        llElements;
UNIT            uCombined, uCurrent;
ASSOC           aAss;
LOOP            lTemp;
RESIDUE         rRes;
OBJEKT          oObj;
char            *sCmd = "combine";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "l" ) ) {
        VPFATALDELAYEDEXIT("usage:  <variable> = combine <LIST>\n" );
        return NULL;
    }
    llElements = llListLoop( LIST_from(oAssocObject(aaArgs[0])) );

                /* Get the first element from the list */

    uCombined = NULL;
    while ( (aAss = ASSOC_from(oListNext(&llElements))) ) {
        oObj = oAssocObject(aAss);
        if ( iObjectType( oObj ) != UNITid ) {
                VP0("%s: %s is type %s\n", sCmd,
                                sAssocName(aAss), sObjectType(oObj) );
                VP0("  expected UNIT\n" );
                continue;
        }
        /*
        **      objekt is a unit, so make a copy
        */
        uCurrent = uCopyUnit(oObj);
        VP1("  Sequence: %s\n", sContainerName(CONTAINER_from( uCurrent)) );
        if ( uCombined == NULL ) {
                MESSAGE("Copying the first UNIT\n" );
                uCombined = uCurrent;
        } else {
                MESSAGE("Copied a subsequent UNIT\n" );
                MESSAGE("Joining two UNITS deleting the second\n" );
                UnitJoin( uCombined, uCurrent );
        }
    }
    if ( uCombined == NULL ) {
        VP0("No UNITS, so no combine performed\n" );
        return NULL;
    }

                /* Define PDB sequence */

    lTemp = lLoop( OBJEKT_from(uCombined), RESIDUES );
    while ( (rRes = RESIDUE_from(oNext(&lTemp))) ) {
        ResidueSetPdbSequence( rRes, iContainerSequence(CONTAINER_from( rRes)) );
    }

    return OBJEKT_from(uCombined);
}




/*
 *      oCmd_source
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Execute a stream of commands from a file.
 *
 *      NOTE: This function MODIFIES the global variable GfCurrentInput.
 *
 *      Arguments:
 *              [0]     OSTRING containing the file.
 */
OBJEKT
oCmd_source( int iArgCount, ASSOC aaArgs[] )
{
FILE            *fCmds;

    if ( !bCmdGoodArguments( "source", iArgCount, aaArgs, "s" ) ) {
        VPFATALDELAYEDEXIT("usage:  source <filename>\n" );
        return NULL;
    }

    fCmds = FOPENCOMPLAIN( sOString(oAssocObject(aaArgs[0])), "r" );

    if ( fCmds != NULL ) {
        if ( bINPUTMAXDEPTH() ) {
            VP0("Source commands are nested too deep!\n" );

        } else {
                /* Push the file onto the input file stack.  The main */
                /* parsing routine will continue parsing until there are */
                /* no more files on the input file stack. */
                VP0("----- Source: %s\n", GsBasicsFullName );
                INPUTPUSHFILE( fCmds );
                VP0("----- Source of %s done\n", GsBasicsFullName );
        }
    }

    return NULL;
}




/*
 *      oCmd_check
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Check the OBJEKT to see if it is ready to have calculations
 *      run on it.  Various tests are run on the OBJEKT like,
 *      do the atoms have types defined, are all the RESIDUES known, etc.
 *
 *      Arguments:
 *              [0]     UNIT, MOLECULE, RESIDUE, ATOM to check.
 *      Davids Changes
 *      Option  [1]     optional PARMSET to place missing parms into.
 *
 */
OBJEKT
oCmd_check( int iArgCount, ASSOC aaArgs[] )
{
int             iErrors, iWarnings;
UNIT            uUnit;
PARMSET         psParmSet;
char            *sCmd = "check";
double          dCloseness = 1.5;
bool            bAbsoluteDist = true;
    switch ( iArgCount ) {
        case 1:
            if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "umra" )) {
                VPFATALDELAYEDEXIT("usage:  check <unit> [parmset] [dCloseness]\n" );
                return NULL;
            }
            break;
        case 2:
            if ( bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u n" )) {
                dCloseness = dODouble(oAssocObject( aaArgs[1] ));
                iArgCount = 1;
            } else if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u ps" )) {
                VPFATALDELAYEDEXIT("usage:  check <unit> [parmset] [dCloseness]\n" );
                return NULL;
            }
            break;
        case 3:
            if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u ps n" )) {
                VPFATALDELAYEDEXIT("usage:  check <unit> [parmset] [dCloseness]\n" );
                return NULL;
            }
            dCloseness = dODouble(oAssocObject( aaArgs[2] ));
            iArgCount = 2;
            break;
        default:
            VP0("%s: A maximum of two parameters are acceptable:\n", sCmd );
            VP0("     a UNIT and an (optional) PARMSET.\n" );
            return  NULL ;
            break;
    }

    if (dCloseness < 0)  {
         bAbsoluteDist = false;
         dCloseness *= -1.0;
    }

                /* Run the checks on the objects */

    iErrors = 0;
    iWarnings = 0;
    VP0("Checking '%s'....\n", sAssocName( aaArgs[0] ));
    ContainerCheck( CONTAINER_from( oAssocObject( aaArgs[0] )), &iErrors, &iWarnings );

                        /* Look for close contacts */
    iWarnings += iToolDistanceSearch( CONTAINER_from(oAssocObject(aaArgs[0])), dCloseness,
                                        bAbsoluteDist, DISTANCE_SEARCH_PRINT_WARNINGS );

    if ( iObjectType( oAssocObject( aaArgs[0] )) == UNITid ) {
        VP0("Checking parameters for unit '%s'.\n", sAssocName( aaArgs[0] ));
        uUnit = UNIT_from(oAssocObject( aaArgs[0] ));
        if ( iArgCount == 2 ) {
            if ( iObjectType(oAssocObject( aaArgs[1] )) == OSTRINGid ) {
                VP0("Creating empty parmset %s\n", sAssocName( aaArgs[1] ) );
                psParmSet = PARMSET_from(oCreate(PARMSETid));
            } else {
                psParmSet = PARMSET_from(oAssocObject( aaArgs[1]));
            }
            UnitCheckForParms( uUnit, GplAllParameters, psParmSet );
            if ( iObjectType(oAssocObject( aaArgs[1] )) == OSTRINGid )
                Destroy((OBJEKT *) &psParmSet );
        } else {
            UnitCheckForParms( uUnit, GplAllParameters, NULL );
        }
    }
    if ( iErrors || iWarnings )
        VP0("%s:  ", sCmd );
    if ( iErrors )
        VP0("Errors:  %d   ", iErrors );
    if ( iWarnings )
        VP0("Warnings: %d\n", iWarnings );

    if ( iErrors == 0 )
        VP0("%s is OK.\n", sObjectType( oAssocObject( aaArgs[0] )));

    return NULL;
}

static void
getUsage()
{
        VPFATALDELAYEDEXIT("usage: variable = get <container> <parameter> [sum|ave|elem]\n");
}

OBJEKT
oCmd_get( int iArgCount, ASSOC aaArgs[] )
{
char    *sCmd = "get";
OBJEKT  oResult;
int     iMode=0;
    if ( iArgCount == 2 && !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "lumra s" ) ) {
        getUsage();
        return NULL;
    } else if ( iArgCount == 3 ) {
        if (!bCmdGoodArguments( sCmd, iArgCount, aaArgs, "lumra s s" ) ) {
            getUsage();
            return NULL;
        }
        char *mode = sOString(oAssocObject(aaArgs[2]));
        if (!strcasecmp(mode,"ave")) iMode=1;
        else if (!strcasecmp(mode,"elem")) iMode=2;
        else if (strcasecmp(mode,"sum")) {
            getUsage();
            return NULL;
        }
    }
    OBJEKT oQuery = oAssocObject(aaArgs[0]);
    char *sItem = sOString(oAssocObject(aaArgs[1]));
    if (!strcasecmp(sItem,"count") && iObjectType(oQuery)==LISTid) {
        int iCount=0;
        OBJEKT oElement;
        LISTLOOP llElements = llListLoop(LIST_from(oQuery));
        while ( (oElement=oListNext(&llElements)) ) iCount++;
        oResult = oCreate(ODOUBLEid);
        ODoubleSet(ODOUBLE_from(oResult), (double)iCount);
    } else if (iObjectType(oQuery)==LISTid) {
        OBJEKT oElement;
        oResult = NULL;
        LISTLOOP llElements = llListLoop(LIST_from(oQuery));
        int iCount=0;
        while ( (oElement=oListNext(&llElements)) ) {
            ASSOC aAssoc = NULL;
            if (iObjectType(oElement)==ASSOCid) {
                aAssoc = ASSOC_from(oElement);
                oElement = oAssocObject(aAssoc);
                MESSAGE("Getting attribute for: %s\n", sAssocName(aAssoc) );
            }
            if ( !bObjectInClass( oElement, CONTAINERid ) ) {
                VPFATALEXIT("%s: Cannot get attribute for %s - not a 'container'\n"
                        "\ttype %s\n", sCmd,
                         (iObjectType(oElement)==ASSOCid ? sAssocName(aAssoc) : "object"),
                         sObjectType(oElement) );
                continue;
            }
            iCount++;
            CONTAINER cElem = CONTAINER_from(oElement);
            OBJEKT oElemResult = oContainerGetAttribute(cElem, sItem);
            STRING sTemp;
            sContainerFullDescriptor(cElem,sTemp);
            if (iObjectType(oElemResult)==OSTRINGid)
                VP2("%s: \"%s\"\n",sTemp,sOString(oElemResult));
            else if (iObjectType(oElemResult)==ODOUBLEid)
                VP2("%s: %g\n",cElem,dODouble(oElemResult));
            else
                VP2("%s: object type %c\n",sTemp,iObjectType(oElemResult));
            if (!oResult) { oResult = oElemResult; }
            else if (iObjectType(oElemResult)==ODOUBLEid)  {
               ODoubleSet(oResult, dODouble(oResult)+dODouble(oElemResult));
            }
            if (iMode==2) break;
            // TODO: how to combine strings?? add vectors, vector GET
        }
        if (iMode==1 && iCount > 0) ODoubleSet(oResult, dODouble(oResult) / (double)iCount);
    } else oResult = oContainerGetAttribute(CONTAINER_from(oQuery), sItem);
    if (iObjectType(oResult)==OSTRINGid)
        VP2("Result: \"%s\"\n",sOString(oResult));
    else if (iObjectType(oResult)==ODOUBLEid)
        VP2("Result: %g\n",dODouble(oResult));
    else
        VP2("Result: object type %c\n",iObjectType(oResult));
    return oResult;

}

/*
 *      oCmd_set
 *
 *      Author: Christian Schafmeister (1991)
 *      Modified:  David A. Rivkin (1-22-93)
 *              Added support for the "Defaults" table.
 *
 *      Set attributes of an object or set a default value.
 *
 *      Arguments:
 *              [0] -  OBJEKT, LIST or STRING whose attribute is to be modified.
 *              [1] -  OSTRING with attribute to modify.
 *              [2] -  OBJEKT new value of attribute.
 */
static void
setUsage()
{
        VPFATALDELAYEDEXIT("usage:  set <container> <parameter> <object>\n"
              "   or:  set [default] <parameter> <value>\n" );
}

OBJEKT
oCmd_set( int iArgCount, ASSOC aaArgs[] )
{
LISTLOOP        llElements;
OBJEKT          oObject;
char            *sString;
char            *sCmd = "set";

    if ( (iArgCount == 2 && !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "s ns")) ||
         (iArgCount != 2 && !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "lumras s *" )) ){
        setUsage();
        return NULL;
    }

    /************ process   set [default] <setting> <value> ************/
    if ( iArgCount == 2 ) {
        return SetDefault(sOString(oAssocObject(aaArgs[0])), oAssocObject(aaArgs[1]));
    } else if (iObjectType(oAssocObject(aaArgs[0])) == OSTRINGid ) {
        /*
         *  setting an 'environmental' default;
         *      catch user attempt to reference non-existent unit
         */
        sString = sOString(oAssocObject(aaArgs[0]));
        if (strchr(sString, '.')) {
                VPFATAL("%s: not a container (e.g. residue)\n", sString );
                setUsage();
                return NULL;
        }
        if ( strcasecmp( sString, "default" )) {
                VPFATAL("%s: expected 'default'\n", sString );
                setUsage();
                return NULL;
        }
        /*
        **  handle the default
        */
        return SetDefault(sOString(oAssocObject(aaArgs[1])), oAssocObject(aaArgs[2]));
    }

    /************ process   set <object> <value> ************/
    DisplayerAccumulateUpdates();
    // object is List of objects
    char *sParam = sOString(oAssocObject(aaArgs[1]));
    OBJEKT oValue = oAssocObject(aaArgs[2]);
    if ( iObjectType(oAssocObject(aaArgs[0])) == LISTid ) {

        llElements = llListLoop(LIST_from(oAssocObject(aaArgs[0])));
        while ( (oObject = oListNext(&llElements)) ) {
            ASSOC aAssoc=NULL;
            if (iObjectType(oObject)==ASSOCid) {
                aAssoc = ASSOC_from(oObject);
                oObject = oAssocObject(aAssoc);
                MESSAGE("Setting attribute for: %s\n", sAssocName(aAssoc) );
            }
            if ( bObjectInClass( oObject, CONTAINERid ) ) {
                ContainerSetAttribute( CONTAINER_from(oObject),sParam,oValue);
            } else {
                VPFATALEXIT("%s: Cannot set attribute for %s - not a 'container'\n"
                        "\ttype %s\n", sCmd,
                        (iObjectType(oObject)==ASSOCid ?  sAssocName(aAssoc) : "object"),
                         sObjectType(oObject) );
            }
        }
    }
    // Scalar/single object
    else {
        MESSAGE("Setting attribute for: %s\n", sAssocName(oAssocObject(aaArgs[0])) );
        if ( bObjectInClass( oAssocObject(aaArgs[0]), CONTAINERid ) ) {
            ContainerSetAttribute( CONTAINER_from( oAssocObject(aaArgs[0])),
                    sOString(oAssocObject(aaArgs[1])),
                    oAssocObject(aaArgs[2]) );
        } else {
            VPFATALEXIT("%s: Cannot set attribute for %s - not a 'container'\n"
                    "\ttype %s\n", sCmd, sAssocName(aaArgs[0]),
                     sObjectType(oAssocObject(aaArgs[0])) );
        }
    }
    DisplayerReleaseUpdates();
    return NULL;
}

/*
 *      oCmd_setBox
 *
 *      Author: Bill Ross (1996)
 *
 *      Arguments:
 *              [0] -   OBJEKT to center/set vdw box.
 *              [1] -   'vdw' or 'centers'
 *              [2] -   optional offset (1 number or list of 3)
 */
OBJEKT
oCmd_setBox( int iArgCount, ASSOC aaArgs[] )
{
UNIT    uUnit;
double  dX, dY, dZ, daBuffer[3];
int     i;
bool    bUsage;
char    *sCmd = "setBox";
char    *sUsage = "usage:  setBox <unit> vdw|centers [ clearance | <clearance_xyz_list> ]\n";
char    *sOpt;
bool    bVdw;

    if ( iArgCount == 2 )
        bUsage = !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u s" );
    else
        bUsage = !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u s ln" );

    if ( bUsage ) {
        VPFATALDELAYEDEXIT("%s",sUsage );
        return NULL;
    }

    uUnit = UNIT_from(oAssocObject(aaArgs[0]));
    sOpt = (char *)sOString(oAssocObject(aaArgs[1]));
    if ( !strcmp(sOpt, "vdw")  ) {
        bVdw = true;
    } else if ( !strcmp(sOpt, "centers")  ) {
        bVdw = false;
    } else {
        VPFATAL("%s: Expected 'vdw' or 'centers' for second argument.\n", sCmd );
        VPFATALDELAYEDEXIT("%s",sUsage );
        return NULL;
    }

    if ( bUnitBoxOct(uUnit) ) {
        VP0(" removing previous octbox..\n" );
        UnitResetFlags(uUnit, UNITBOXOCT);
        UnitResetFlags(uUnit, UNITUSEBOUNDINGBOX);
    } else if ( bUnitUseBox( uUnit ) ) {
        VP0(" removing previous box..\n" );
        UnitResetFlags(uUnit, UNITUSEBOUNDINGBOX);
    }

    /*
     *  get any box offset
     */
    if ( iArgCount == 3 ) {
        if ( iObjectType( oAssocObject(aaArgs[2]) ) == LISTid ) {
            LISTLOOP        llNumbers;

            /*
             *  x,y,z box clearances
             */
            if ( iListSize(oAssocObject(aaArgs[2])) != 3 ) {
                VPFATAL("%s: Expected 3 non-negative floating point numbers"
                        " { x y z } for third argument.\n",
                        sCmd );
                VPFATALDELAYEDEXIT("%s",sUsage );
                return NULL;
            }
            llNumbers = llListLoop(LIST_from(oAssocObject(aaArgs[2])));
            for ( i=0; i<3; i++ ) {
                ASSOC   aAss = ASSOC_from(oListNext(&llNumbers));
                if ( iObjectType(oAssocObject(aAss)) != ODOUBLEid ) {
                    VPFATAL("%s: Bad value #%d in third argument.\n", sCmd, i );
                    VP0("%s: Expected 3 non-negative floating point numbers"
                            " { x y z } for third argument.\n",
                            sCmd );
                    VPFATALDELAYEDEXIT("%s",sUsage );
                    return NULL;
                }
                daBuffer[i] = 2.0 * dODouble(oAssocObject(aAss));
                if ( daBuffer[i] < 0.0 ) {
                    VPFATAL("%s: Bad value #%d in third argument.\n", sCmd, i );
                    VP0("%s: Expected 3 non-negative floating point numbers"
                            " { x y z } for third argument.\n",
                            sCmd );
                    VPFATALDELAYEDEXIT("%s",sUsage );
                    return NULL;
                }
            }
        } else {
            daBuffer[0] = 2.0 * dODouble(oAssocObject(aaArgs[2]));
            if ( daBuffer[0] < 0.0 ) {
                VPFATAL("%s: Bad box clearance.\n", sCmd );
                VP0("%s: Expected a non-negative floating point number"
                        " for third argument.\n",
                        sCmd );
                VPFATALDELAYEDEXIT("%s",sUsage );
                return NULL;
            }
            daBuffer[2] = daBuffer[1] = daBuffer[0];
        }
    }

    DisplayerAccumulateUpdates();
    if( GDefaults.bNoCenter == 0 ){
        if (bVdw) {
            ToolCenterUnitByRadii( uUnit, false );
        } else {
//XXX This actually centers the box, which is then shifted back at unitio.c
//    This does not support arbitrary box, it just uses 0,1 interval instead of -1/2,1/2
            VECTOR vOrigin;
            VectorDef( &vOrigin, 0.0, 0.0, 0.0 );
            ContainerCenterAt( CONTAINER_from(uUnit), vOrigin );
            ToolSetUnitBoxByCenters( uUnit );
        }
    }
    DisplayerReleaseUpdates();

    UnitSetUseBox( uUnit, true );
    UnitSetBeta( uUnit, 90.0*DEGTORAD );
    UnitGetBox( uUnit, &dX, &dY, &dZ );
    if ( iArgCount == 3 ) {
        dX += daBuffer[0];
        dY += daBuffer[1];
        dZ += daBuffer[2];
        UnitSetBox( uUnit, dX, dY, dZ );
    }
    VP0("Box dimensions:  %f %f %f\n", dX, dY, dZ );
    return NULL;
}

/*
 *      oCmd_setCell
 *
 *      Author: Juno Krahn (2026)
 *
 *      Arguments:
 *              [0] -   OBJEKT to center/set vdw box.
 *              [1] - [6] -- full unit cell parameters
 */
OBJEKT
oCmd_setCell( int iArgCount, ASSOC aaArgs[] )
{
UNIT    uUnit;
bool    bUsage;
char    *sCmd = "setCell";
char    *sUsage = "usage:  setCell <unit> <a> <b> <c> [ <alpha> <beta> <gamma> | <beta> ]\n";

    if ( iArgCount == 4 )
        bUsage = !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u n n n" );
    else if ( iArgCount == 5 )
        bUsage = !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u n n n n" );
    else
        bUsage = !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u n n n n n n" );

    if ( bUsage ) {
        VPFATALDELAYEDEXIT("%s",sUsage );
        return NULL;
    }

    uUnit = UNIT_from(oAssocObject(aaArgs[0]));
    double a = dODouble(oAssocObject(aaArgs[1]));
    double b = dODouble(oAssocObject(aaArgs[2]));
    double c = dODouble(oAssocObject(aaArgs[3]));
    double alpha=90.0, beta=90.0, gamma=90.0;
    if (iArgCount == 7) {
       alpha = dODouble(oAssocObject(aaArgs[4]));
       beta = dODouble(oAssocObject(aaArgs[5]));
       gamma = dODouble(oAssocObject(aaArgs[6]));
    } else if (iArgCount == 5)
       beta = dODouble(oAssocObject(aaArgs[4]));
    UnitSetCell(uUnit,a,b,c,alpha,beta,gamma);
    UnitSetUseBox( uUnit, true );
    VP0("Cell dimensions:  %f %f %f %f %f %f\n", a,b,c,alpha,beta,gamma );
    return NULL;
}

/*
 * Show default variables
 *      Author: Yong Duan
 *      (adapted from set default)
 */
OBJEKT
oCmd_showDefault( int iArgCount, ASSOC aaArgs[] )
{
char    *sCmd = "showDefault";
    if (!bCmdGoodArguments( sCmd, iArgCount, aaArgs, (iArgCount == 1) ? "s" : "" )) {
        VPFATALDELAYEDEXIT("usage:  showDefault [<parameter_name>]\n" );
        return NULL;
    }
    if ( iArgCount == 1 ) PrintDefaultSettings(sOString(oAssocObject(aaArgs[0])));
    else PrintDefaultSettings(NULL);
    return NULL;
}

/*
 *      oCmd_charge
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Arguments:
 *              [0] -   uNIT whose charge should be calculated.
 *
 */

typedef struct {
    double dFrac;
    double dCharge;
    RESIDUE rRes;
} CHARGERESt;

int qsort_rescharge(const void *a, const void *b) {
    double dA = ((const CHARGERESt*)a)->dFrac;
    double dB = ((const CHARGERESt*)b)->dFrac;
    if (dA < dB) return 1;
    if (dA > dB) return -1;
    return 0;
}

OBJEKT
oCmd_charge( int iArgCount, ASSOC aaArgs[] )
{
CONTAINER               cCont;
double                  dCharge, dPertCharge;

    if ( !bCmdGoodArguments( "charge", iArgCount, aaArgs, "umral" ) ) {
        VPFATALDELAYEDEXIT("usage:  charge <object>\n" );
        return NULL;
    }

    cCont = CONTAINER_from(oAssocObject(aaArgs[0]));

    ContainerTotalCharge( cCont, &dCharge, &dPertCharge );
    VP0("Total unperturbed charge: %10.6lf\n", dCharge );
    /*dPertCharge is now the delta of the charge and not the actual perturbed charge*/
    VP0("Total perturbed charge:   %10.6lf\n", (dCharge+dPertCharge) );

    int iCount=0, iNumRes = iContainerNumberOfChildren(cCont);
    CHARGERESt *crResidues = malloc(sizeof(CHARGERESt)*iNumRes);
    RESIDUE rRes;
    LOOP lRes = lLoop( OBJEKT_from(cCont), RESIDUES );
    while ( ( rRes = RESIDUE_from(oNext(&lRes)) ) != NULL ) {
            double charge=0;
            ATOM aAtom;
            FOR_ATOMS_IN_RESIDUE(aAtom,rRes)
                charge += dAtomCharge(aAtom);
            crResidues[iCount].dFrac = fabs(charge - round(charge));
            crResidues[iCount].dCharge = charge;
            crResidues[iCount++].rRes = rRes;
    }
    qsort((void *)crResidues, iNumRes, sizeof(CHARGERESt), qsort_rescharge);
    if (crResidues[0].dFrac > 0.00005) {
        VP0("Residues with the largest deviation from unit charge:\n");
        for (int i=0;i<iNumRes && i<20 && crResidues[i].dFrac > 0.00005;i++) {
            STRING sDesc;
            VP1("%15s  Q=%7.3f  dq=%8.4f\n",sContainerDescriptor(CONTAINER_from(crResidues[i].rRes),sDesc),
                     crResidues[i].dCharge,crResidues[i].dFrac);
        }
    }
    free(crResidues);
    return NULL;
}



/*
 *      oCmd_saveAmberParm
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Arguments:
 *              [0] -   UNIT which should be saved.
 *              [1] -   OSTRING Filename to which parmfile should be written.
 *              [2] -   OSTRING Filename to which coordfile should be written.
 *
 */
OBJEKT
oCmd_saveAmberParm( int iArgCount, ASSOC aaArgs[] )
{
UNIT    uUnit;
char    *sCmd = "saveAmberParm";

    VPTRACEENTER("oCmd_saveAmberParm" );
    VPTRACEMULTIPLEEXIT("oCmd_saveAmberParm" );
    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u s s" ) ) {
        VPFATALDELAYEDEXIT("usage:  saveAmberParm <unit> <topologyfile> <coordfile> \n" );
        return NULL;
    }
    if ( iParmLibSize(GplAllParameters) == 0 ) {
        VPWARN("%s: There are no parameter sets loaded\n", sCmd );
        return NULL;
    }
    uUnit = UNIT_from(oAssocObject(aaArgs[0]));
    TurnOffDisplayerUpdates();
    UnitSaveAmberParmFile(uUnit, sOString(oAssocObject(aaArgs[1])), sOString(oAssocObject(aaArgs[2])),
                GplAllParameters, false, false, false);
    TurnOnDisplayerUpdates();

    return NULL;
}

/*
 *      oCmd_saveAmberParmNetCDF
 *
 *      Author: Robin Betz (2013)
 *
 *      Arguments:
 *              [0] -   UNIT which should be saved.
 *              [1] -   OSTRING Filename to which parmfile should be written.
 *              [2] -   OSTRING Filename to which coordfile should be written.
 *
 */
OBJEKT
oCmd_saveAmberParmNetCDF( int iArgCount, ASSOC aaArgs[] )
{
UNIT    uUnit;
char    *sCmd = "saveAmberParmNetCDF";

    VPTRACEENTER("oCmd_saveAmberParmNetCDF" );
    VPTRACEMULTIPLEEXIT("oCmd_saveAmberParmNetCDF" );
    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u s s" ) ) {
        VPFATALDELAYEDEXIT("usage: saveAmberParmNetCDF <unit> <topologyfile>"
                " <coordfile> \n" );
        return NULL;
    }

    if ( iParmLibSize(GplAllParameters) == 0 ) {
        VPWARN("%s: There are no parameter sets loaded\n", sCmd );
        return NULL;
    }
    uUnit = UNIT_from(oAssocObject(aaArgs[0]));
    TurnOffDisplayerUpdates();
    UnitSaveAmberParmFile( uUnit, sOString(oAssocObject(aaArgs[1])), sOString(oAssocObject(aaArgs[2])),
                GplAllParameters, false, false, true);
    TurnOnDisplayerUpdates();

    return NULL;
}


/*
 *      oCmd_saveAmberParmPol
 *
 *      Arguments:
 *              [0] -   UNIT which should be saved.
 *              [1] -   OSTRING Filename to which parmfile should be written.
 *              [2] -   OSTRING Filename to which coordfile should be written.
 *
 */
OBJEKT
oCmd_saveAmberParmPol( int iArgCount, ASSOC aaArgs[] )
{
UNIT    uUnit;
char    *sCmd = "saveAmberParmPol";

    VPTRACEENTER("oCmd_saveAmberParmPol" );
    VPTRACEMULTIPLEEXIT("oCmd_saveAmberParmPol" );
    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u s s" ) ) {
        VPFATALDELAYEDEXIT("usage:  %s <unit> <topologyfile> <coordfile>\n",
                sCmd );
        return NULL;
    }

    if ( iParmLibSize(GplAllParameters) == 0 ) {
        VPWARN("%s: There are no parameter sets loaded\n", sCmd );
        return NULL;
    }
    uUnit = UNIT_from(oAssocObject(aaArgs[0]));
    TurnOffDisplayerUpdates();
    UnitSaveAmberParmFile( uUnit, sOString(oAssocObject(aaArgs[1])), sOString(oAssocObject(aaArgs[2])),
                GplAllParameters, true, false, false);
    TurnOnDisplayerUpdates();
    return NULL;
}


/*
 *      oCmd_saveAmberParmPert
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Arguments:
 *              [0] -   UNIT which should be saved.
 *              [1] -   OSTRING Filename to which parmfile should be written.
 *              [2] -   OSTRING Filename to which coordfile should be written.
 *
 */
OBJEKT
oCmd_saveAmberParmPert( int iArgCount, ASSOC aaArgs[] )
{
UNIT    uUnit;
char    *sCmd = "saveAmberParmPert";

    VPTRACEENTER("oCmd_saveAmberParmPert" );
    VPTRACEMULTIPLEEXIT("oCmd_saveAmberParmPert" );
    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u s s" ) ) {
        VPFATALDELAYEDEXIT("usage:  %s <unit> <topologyfile> <coordfile>\n",
                sCmd );
        return NULL;
    }

    if ( iParmLibSize(GplAllParameters) == 0 ) {
        VPWARN("%s: There are no parameter sets loaded\n", sCmd );
        return NULL;
    }
    uUnit = UNIT_from(oAssocObject(aaArgs[0]));
    TurnOffDisplayerUpdates();
    UnitSaveAmberParmFile( uUnit, sOString(oAssocObject(aaArgs[1])), sOString(oAssocObject(aaArgs[1])),
                GplAllParameters, false, true, false);
    TurnOnDisplayerUpdates();

    return NULL;
}


/*
 *      oCmd_saveAmberParmPolPert
 *
 *      Arguments:
 *              [0] -   UNIT which should be saved.
 *              [1] -   OSTRING Filename to which parmfile should be written.
 *              [2] -   OSTRING Filename to which coordfile should be written.
 *
 */
OBJEKT
oCmd_saveAmberParmPolPert( int iArgCount, ASSOC aaArgs[] )
{
UNIT    uUnit;
char    *sCmd = "saveAmberParmPert";

    VPTRACEENTER("oCmd_saveAmberParmPolPert" );
    VPTRACEMULTIPLEEXIT("oCmd_saveAmberParmPolPert" );
    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u s s" ) ) {
        VPFATALDELAYEDEXIT("usage:  %s <unit> <topologyfile> <coordfile>\n",
                sCmd );
        return NULL;
    }

    if ( iParmLibSize(GplAllParameters) == 0 ) {
        VPWARN("%s: There are no parameter sets loaded\n", sCmd );
        return NULL;
    }
    uUnit = UNIT_from(oAssocObject(aaArgs[0]));
    TurnOffDisplayerUpdates();
    UnitSaveAmberParmFile( uUnit, sOString(oAssocObject(aaArgs[1])), sOString(oAssocObject(aaArgs[2])),
                GplAllParameters, true, true, false);
    TurnOnDisplayerUpdates();
    return NULL;
}



OBJEKT
oCmd_saveAmberPrep( int iArgCount, ASSOC aaArgs[] )
{
UNIT            uUnit;
FILE            *fOut;
char            *sCmd = "saveAmberPrep";

    VPTRACEENTER("oCmd_saveAmberPrep" );
    VPTRACEMULTIPLEEXIT("oCmd_saveAmberPrep" );
    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u s" ) ) {
        VPFATALDELAYEDEXIT("usage:  saveAmberPrep <unit> <file>\n" );
        return NULL;
    }

    uUnit = UNIT_from(oAssocObject(aaArgs[0]));

    fOut = FOPENCOMPLAIN( sOString(oAssocObject(aaArgs[1])), "w" );
    if ( fOut == NULL ) {
        VPFATALEXIT("%s: Could not open file: %s\n",
                sCmd, sOString(oAssocObject(aaArgs[1])));
        return NULL;
    }

    TurnOffDisplayerUpdates();
    UnitIOSaveAmberPrep( uUnit, fOut );
    TurnOnDisplayerUpdates();

    fclose( fOut );
    VP0("  -- Remember to delete unwanted IMPROPER terms!\n" );
    return NULL;
}


/*
 *      oCmd_clearVariables
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Erase all variables.
 */
OBJEKT
oCmd_clearVariables( int iArgCount, ASSOC aaArgs[] )
{
LISTLOOP        llVariables;
ASSOC           aAssoc;

    if ( iArgCount == 0 ) {
        VP0("Clearing all variables\n" );
        VariablesDestroy();
        VariablesInit();
    } else {
        if ( !bCmdGoodArguments( "clearVariables", iArgCount, aaArgs, "l" ) ) {
            VPFATALDELAYEDEXIT("usage:  clearVariables [LIST]\n" );
            return NULL;
        }
        llVariables = llListLoop( LIST_from(oAssocObject(aaArgs[0])) );
        while ( (aAssoc = ASSOC_from(oListNext(&llVariables))) ) {
            VariableRemove( sAssocName(aAssoc) );
        }
    }
    return NULL;
}






/*
 *      oCmd_matchVariables
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return a list of all variables that match the name.
 *
 *      Arguments:
 *              [0] -   OSTRING, name with '*' and '?' wildcards to match.
 *
 */
OBJEKT
oCmd_matchVariables( int iArgCount, ASSOC aaArgs[] )
{
LIST            lVars;
DICTIONARY      dVariables;
DICTLOOP        dlEntries;
ASSOC           aAssoc;
char            *sCmd = "matchVariables";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "s" ) ) {
        VPFATALDELAYEDEXIT("usage:  <variable> = matchVariables <string>\n" );
        return NULL;
    }
    dVariables = dVariablesDictionary();
    lVars = LIST_from(oCreate(LISTid));

    dlEntries = ydlDictionaryLoop(dVariables);
    while ( yPDictionaryNext( dVariables, &dlEntries ) ) {
        if ( bStringMatchPattern( sDictLoopKey(dlEntries),
                                  sOString(oAssocObject(aaArgs[0])) ) ) {
            aAssoc = ASSOC_from(oCreate(ASSOCid));
            AssocSetName( aAssoc, sDictLoopKey(dlEntries) );
            AssocSetObject( aAssoc, PDictLoopData(dlEntries) );
            ListAddToEnd( lVars, OBJEKT_from(aAssoc ));
        }
    }

    return OBJEKT_from(lVars);
}




/*
 *      oCmd_edit
 *
 *      Author: Christian Schafmeister (1991)
 *      Modified: David A. Rivkin (1992)
 *
 *      If in a Graphical environment setup the GrMainResult structure
 *      to return the UNIT that the USER wishes to edit.  Otherwise
 *      print a disclaimer that 'edit' does not work in a non-graphical
 *      environment.
 *
 *      This command causes a COPY of the OBJEKT to be edited.
 *      If the OBJEKT that the user wants to edit is NULL then
 *      this command assumes that they want to edit a new
 *      UNIT with a single RESIDUE within it.
 *
 *      Arguments:
 *              [0] -   OBJEKT to edit.
 *      Davids Changes - changed "uz" to "ups" in bCmdGoodArguments call.
 *              Moved and copied the "GrMainResult.oObject = uUnit" to the
 *              first two if statements and added a check and
 *              GrMainResult.oObject = psParmSet added PARMSET psParmSet
 *              to the local variables.
 */
OBJEKT
oCmd_edit( int iArgCount, ASSOC aaArgs[] )
{
UNIT            uUnit;
PARMSET         psParmSet;
RESIDUE         rRes;

    if ( !GbGraphicalEnvironment ) {
        VPWARN("The edit command only works in a graphical environment.\n" );
        GrMainResult.iCommand = CMD_NONE;
        return NULL;
    }
    if ( !bCmdGoodArguments( "edit", iArgCount, aaArgs, "ups" ) ) {
        VPFATALDELAYEDEXIT("usage:  edit <unit/parmset>\n" );
        return NULL;
    }

    GrMainResult.iCommand = CMD_EDIT;
    strcpy( GrMainResult.sVariable, sAssocName(aaArgs[0]) );

                /* Check if the UNIT exists, if it doesn't then */
                /* create a new one with a single RESIDUE in it */

    if ( oAssocObject(aaArgs[0]) == NULL  ||
                        iObjectType(oAssocObject(aaArgs[0])) == OSTRINGid ) {
            VP0("Creating a new, empty UNIT \"%s\"\n",
                        sOString( oAssocObject(aaArgs[0])));
            uUnit = UNIT_from(oCreate(UNITid));
            ContainerSetName( CONTAINER_from( uUnit),
                                sOString( oAssocObject(aaArgs[0])));
            rRes  = RESIDUE_from(oCreate(RESIDUEid));
            ContainerSetName( CONTAINER_from( rRes),
                                sOString(oAssocObject(aaArgs[0])));
            ContainerAdd( CONTAINER_from(uUnit), OBJEKT_from(rRes ));
            VariableSet( sOString( oAssocObject(aaArgs[0])), OBJEKT_from(uUnit ));   /* adds 1 REF */
            GrMainResult.oObject = OBJEKT_from( uUnit);

    } else if ( iObjectType(oAssocObject(aaArgs[0])) == UNITid ) {
            uUnit = UNIT_from(oAssocObject(aaArgs[0]));
            GrMainResult.oObject = OBJEKT_from( uUnit);

    } else if ( iObjectType(oAssocObject(aaArgs[0])) == PARMSETid) {
            psParmSet = PARMSET_from(oAssocObject(aaArgs[0]));
            GrMainResult.oObject = OBJEKT_from( psParmSet);

    } else {
            VPWARN("Can't edit type %s\n",
                                sObjectType(oAssocObject(aaArgs[0])) );
            GrMainResult.iCommand = CMD_NONE;
    }

    return NULL;
}




/*
 *      oCmd_alignAxes
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Align the principle axis of the molecule along the
 *      coordinate axis.
 *
 *      Arguments:
 *              [0] -   UNIT to align.
 *
 */
OBJEKT
oCmd_alignAxes( int iArgCount, ASSOC aaArgs[] )
{

    if ( !bCmdGoodArguments( "alignAxes", iArgCount, aaArgs, "u" ) ) {
        VPFATALDELAYEDEXIT("usage:  alignAxes <unit>\n" );
        return NULL;
    }

    DisplayerAccumulateUpdates();

    ToolOrientPrincipleAxisAlongCoordinateAxis( UNIT_from(oAssocObject(aaArgs[0]) ));

                /* Instruct the graphics system to reset the viewing */
                /* matrices of the UNIT */

    DisplayerReleaseUpdates();

    return NULL;
}

OBJEKT
oCmd_resequence( int iArgCount, ASSOC aaArgs[] )
{
OBJEKT oUnit, oResidue, oAtom;
LOOP lResidues, lAtoms;
int iResidue, iAtom, iNonResidue=0, iNonAtom=0, iAllAtoms=0;
STRING sResidue, sAtom;
char *sUnit;

    if ( !bCmdGoodArguments( "resequence", iArgCount, aaArgs, "u" ) ) {
        VPFATALDELAYEDEXIT("usage:  resequence <unit>\n" );
        return NULL;
    }
    oUnit = oAssocObject(aaArgs[0]);
    sUnit = sAssocName(aaArgs[0]);

    iResidue=0;
    LOOPOVERALL(oUnit,DIRECTCONTENTSBYSEQNUM,oResidue,OBJEKT,lResidues) {
        ContainerSetSequence(oResidue,++iResidue);
        if (iObjectType(oResidue)!=RESIDUEid) {
            VPWARN("UNIT %s contains a non-RESIDUE member %s\n",
                    sUnit, sContainerFullDescriptor((CONTAINER)oResidue, sResidue));
            continue;
        }
        iAtom=0;
        LOOPOVERALL(oResidue,DIRECTCONTENTSBYSEQNUM,oAtom,OBJEKT,lAtoms) {
            ContainerSetSequence(oAtom,++iAtom);
            if (iObjectType(oAtom)!=ATOMid) {
                VPWARN("RESIDUE %s contains a non-ATOM member %s\n",
                        sContainerFullDescriptor((CONTAINER)oResidue, sResidue),
                        sContainerFullDescriptor((CONTAINER)oAtom, sAtom));
            }
        }
        ContainerSetNextChildsSequence(oResidue,iAtom+1);
        iAllAtoms += iAtom;
    }
    ContainerSetNextChildsSequence(oUnit,iResidue+1);
    VP0("Unit %s contains %d residues and %d atoms\n", sUnit, iResidue, iAllAtoms);
    if (iNonResidue || iNonAtom)
        VP0("Unit contains %d non-residues, with residues containing %d non-atoms\n"
            "Contents of unexpected objects were not renumbered.\n",
            sUnit, iNonResidue, iNonAtom);
    return NULL;
}



/*
 *      oCmd_select
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      This command sets the ATOMSELECTED flag of atoms
 *      within the OBJEKT
 *
 *      Arguments:
 *              [0]     - CONTAINER/LIST.
 */

OBJEKT
oCmd_select( int iArgCount, ASSOC aaArgs[] )
{
OBJEKT          oOver;
LOOP            lAtoms;
ATOM            aAtom;

    if ( !bCmdGoodArguments( "select", iArgCount, aaArgs, "umral" ) ) {
        VPFATALDELAYEDEXIT("usage:  select <object>\n" );
        return NULL;
    }

    DisplayerAccumulateUpdates();

    oOver = oAssocObject(aaArgs[0]);
    lAtoms = lLoop( oOver, ATOMS );
    while ( (aAtom = ATOM_from(oNext(&lAtoms)))) {
        AtomSetFlags( aAtom, ATOMSELECTED );
    }
    DisplayerReleaseUpdates();
    return NULL;
}

OBJEKT
oCmd_selectMask( int iArgCount, ASSOC aaArgs[] )
{
    if ( !bCmdGoodArguments( "selectMask", iArgCount, aaArgs, "u s" ) ) {
        VPFATALDELAYEDEXIT("usage:  selectMask <unit> <string>\n" );
        return NULL;
    }
    VARARRAY vaAtoms = vaAtomMaskSelect(
                             UNIT_from(oAssocObject(aaArgs[0])),
                             sOString(oAssocObject(aaArgs[1])));
    ATOM *aPAtom = PVAI( vaAtoms, ATOM, 0 );
    int iNumAtoms = iVarArrayElementCount(vaAtoms);
    int iSet = 0;
    for (int i=0; i<iNumAtoms; i++, aPAtom++) {
        if (!bAtomFlagsSet( *aPAtom, ATOMSELECTED )) {
           iSet++;
           AtomSetFlags( *aPAtom, ATOMSELECTED );
        }
        if (iVerbosity()>2) {
            STRING sTemp;
            VP0("%d) %s\n",i, sContainerFullDescriptor(&(*aPAtom)->cHeader,sTemp));
        }
    }
    VP0("%d atoms match Atom Mask; %d atoms added to selected status\n",iNumAtoms, iSet);
    VarArrayDestroy(&vaAtoms);
    return NULL;
}

OBJEKT
oCmd_deSelectMask( int iArgCount, ASSOC aaArgs[] )
{
    if ( !bCmdGoodArguments( "deSelectMask", iArgCount, aaArgs, "u s" ) ) {
        VPFATALDELAYEDEXIT("usage:  deSelectMask <unit> <string>\n" );
        return NULL;
    }
    UNIT uUnit = UNIT_from(oAssocObject(aaArgs[0]));
    SELNODE selNode = selParseAtomMask(sOString(oAssocObject(aaArgs[1])));
    VARARRAY vaAtoms = vaUnitEvalSelection(selNode, uUnit);
    SelFree(selNode);
    ATOM *aPAtom = PVAI( vaAtoms, ATOM, 0 );
    int iNumAtoms = iVarArrayElementCount(vaAtoms);
    int iReset = 0;
    for (int i=0; i<iNumAtoms; i++, aPAtom++) {
        if (bAtomFlagsSet( *aPAtom, ATOMSELECTED )) {
           iReset++;
           AtomResetFlags( *aPAtom, ATOMSELECTED );
        }
        if (iVerbosity()>2) {
            STRING sTemp;
            VP0("%d) %s\n",i, sContainerFullDescriptor(&(*aPAtom)->cHeader,sTemp));
        }
    }
    VP0("%d atoms match Atom Mask; %d atoms removed from selected status\n",iNumAtoms, iReset);
    VarArrayDestroy(&vaAtoms);
    return NULL;
}

/*
 *      oCmd_scaleCharges
 *
 *      scale charges; useful for polar setup
 */
OBJEKT
oCmd_scaleCharges( int iArgCount, ASSOC aaArgs[] )
{
OBJEKT          oOver;
LOOP            lAtoms;
ATOM            aAtom;
double          dScale;

    if ( !bCmdGoodArguments( "scaleCharges", iArgCount, aaArgs, "umral n" ) ) {
        VPFATALDELAYEDEXIT("usage:  scaleCharges <object> <scale_factor>\n" );
        return NULL;
    }

    dScale = dODouble(oAssocObject(aaArgs[1]));
    if ( dScale <= 0.0 ) {
        VPFATAL("scaleCharges: scale_factor must be > 0\n" );
        VPFATALDELAYEDEXIT("usage:  scaleCharges <object> <scale_factor>\n" );
        return NULL;
    }

    oOver = oAssocObject(aaArgs[0]);
    lAtoms = lLoop( oOver, ATOMS );
    while ( (aAtom = ATOM_from(oNext(&lAtoms)))) {
        /* HACK - resetting charge without marking for update */
        dAtomCharge( aAtom ) *= dScale;
        dAtomPertCharge( aAtom ) *= dScale;
    }
    return NULL;
}

/*
 *      oCmd_scaleCharges
 *
 *      scale charges; useful for polar setup
 */
OBJEKT
oCmd_scaleCoor( int iArgCount, ASSOC aaArgs[] )
{
OBJEKT          oOver;
LOOP            lAtoms;
ATOM            aAtom;
double          dScale;

    if ( !bCmdGoodArguments( "scaleCoor", iArgCount, aaArgs, "umral n" ) ) {
        VPFATALDELAYEDEXIT("usage:  scaleObject <object> <scale_factor>\n" );
        return NULL;
    }

    oOver = oAssocObject(aaArgs[0]);
    dScale = dODouble(oAssocObject(aaArgs[1]));

    if ( iObjectType(oOver) == UNITid ) {
        UNIT uUnit = UNIT_from(oOver);
        uUnit->dXWidth *= dScale;
        uUnit->dYWidth *= dScale;
        uUnit->dZWidth *= dScale;
    }

    lAtoms = lLoop( oOver, ATOMS );
    while ( (aAtom = ATOM_from(oNext(&lAtoms)))) {
        VECTOR vNewPosition = {
            aAtom->vPosition.dX * dScale,
            aAtom->vPosition.dY * dScale,
            aAtom->vPosition.dZ * dScale
        };
        AtomSetPosition(aAtom, vNewPosition);
    }
    return NULL;
}


/*
 *      oCmd_deSelect
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      This command resets the ATOMSELECTED flag of atoms
 *      within the OBJEKT
 *
 *      Arguments:
 *              [0]     - CONTAINER/LIST.
 */

OBJEKT
oCmd_deSelect( int iArgCount, ASSOC aaArgs[] )
{
OBJEKT          oOver;
LOOP            lAtoms;
ATOM            aAtom;

    if ( !bCmdGoodArguments( "deSelect", iArgCount, aaArgs, "umral" ) ) {
        VPFATALDELAYEDEXIT("usage:  deSelect <object>\n" );
        return NULL;
    }

    DisplayerAccumulateUpdates();

    oOver = oAssocObject(aaArgs[0]);
    lAtoms = lLoop( oOver, ATOMS );
    while ( (aAtom = ATOM_from(oNext(&lAtoms)))) {
        AtomResetFlags( aAtom, ATOMSELECTED );
    }

    DisplayerReleaseUpdates();

    return NULL;
}




/*
 *      oCmd_restrainBond
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Add a RESTRAINT bond to the UNIT, make sure that the atoms
 *      are contained within the UNIT.
 *
 *      Arguments:
 *              [0]     - The UNIT to add the RESTRAINT bond to.
 *              [1]     - The first ATOM of the RESTRAINT bond.
 *              [2]     - The second ATOM of the RESTRAINT bond.
 *              [3]     - ODOUBLE: The value of Kr.
 *              [4]     - ODOUBLE: The value of R0.
 */
OBJEKT
oCmd_restrainBond( int iArgCount, ASSOC aaArgs[] )
{
int             i;
UNIT            uUnit;
ATOM            aaAtoms[ATOMSINBOND];
double          dKr, dR0;
RESTRAINT       rRest;
char            *sCmd = "restrainBond";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u a a n n" ) ) {
        VPFATALDELAYEDEXIT("usage:  restrainBond <unit> <a> <b> <force> <length>\n" );
        return NULL;
    }

    uUnit  = UNIT_from(oAssocObject(aaArgs[0]));
    for ( i=0; i<ATOMSINBOND; i++ )
        aaAtoms[i] = ATOM_from(oAssocObject(aaArgs[1+i]));
    dKr    = dODouble(oAssocObject(aaArgs[3]));
    dR0    = dODouble(oAssocObject(aaArgs[4]));

    for ( i=0; i<ATOMSINBOND; i++ ) {
        if ( !bContainerContainedBy( CONTAINER_from(aaAtoms[i]), CONTAINER_from(uUnit )) ) {
            VPFATALEXIT("%s: Atom#%d is not part of the UNIT\n", sCmd, i+1 );
            return NULL;
        }
    }

                /* Create the RESTRAINT and add it to the UNIT */

    rRest = rRestraintCreate();
    RestraintBondSet( rRest, aaAtoms[0], aaAtoms[1], dKr, dR0 );
    UnitAddRestraint( uUnit, rRest );

    return NULL;
}





/*
 *      oCmd_restrainAngle
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Add a RESTRAINT angle to the UNIT, make sure that the atoms
 *      are contained within the UNIT.
 *
 *      Arguments:
 *              [0]     - The UNIT to add the RESTRAINT to.
 *              [1]     - The first ATOM of the RESTRAINT.
 *              [2]     - The second ATOM of the RESTRAINT.
 *              [3]     - The third ATOM of the RESTRAINT.
 *              [4]     - ODOUBLE: The value of Kt.
 *              [5]     - ODOUBLE: The value of T0.
 */
OBJEKT
oCmd_restrainAngle( int iArgCount, ASSOC aaArgs[] )
{
int             i;
UNIT            uUnit;
ATOM            aaAtoms[ATOMSINANGLE];
double          dKt, dT0;
RESTRAINT       rRest;
char            *sCmd = "restrainAngle";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u a a a n n" ) ) {
        VPFATALDELAYEDEXIT("usage:  restrainAngle <unit> <a> <b> <c> <force> <length>\n" );
        return NULL;
    }

    uUnit  = UNIT_from(oAssocObject(aaArgs[0]));
    for ( i=0; i<ATOMSINANGLE; i++ )
        aaAtoms[i] = ATOM_from(oAssocObject(aaArgs[1+i]));
    dKt    = dODouble(oAssocObject(aaArgs[4]));
    dT0    = DEGTORAD * dODouble(oAssocObject(aaArgs[5]));

    for ( i=0; i<ATOMSINANGLE; i++ ) {
        if ( !bContainerContainedBy( CONTAINER_from(aaAtoms[i]), CONTAINER_from(uUnit )) ) {
            VPFATALEXIT("%s: Atom#%d is not part of the UNIT\n", sCmd, i+1 );
            return NULL;
        }
    }

                /* Create the RESTRAINT and add it to the UNIT */

    rRest = rRestraintCreate();
    RestraintAngleSet( rRest, aaAtoms[0], aaAtoms[1], aaAtoms[2], dKt, dT0 );
    UnitAddRestraint( uUnit, rRest );

    return NULL;
}





/*
 *      oCmd_restrainTorsion
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Add a RESTRAINT torsion to the UNIT, make sure that the atoms
 *      are contained within the UNIT.
 *
 *      Arguments:
 *              [0]     - The UNIT to add the RESTRAINT bond to.
 *              [1]     - The first ATOM of the RESTRAINT bond.
 *              [2]     - The second ATOM of the RESTRAINT bond.
 *              [3]     - The third ATOM of the RESTRAINT bond.
 *              [4]     - The fourth ATOM of the RESTRAINT bond.
 *              [5]     - ODOUBLE: The value of Kp.
 *              [6]     - ODOUBLE: The value of P0.
 *              [7]     - ODOUBLE: The value of N.
 */
OBJEKT
oCmd_restrainTorsion( int iArgCount, ASSOC aaArgs[] )
{
int             i;
UNIT            uUnit;
ATOM            aaAtoms[ATOMSINTORSION];
double          dKp, dP0, dN;
RESTRAINT       rRest;
char            *sCmd = "restrainTorsion";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u a a a a n n n" ) ) {
        VPFATALDELAYEDEXIT("usage:  restrainTorsion <unit> <a> <b> <c> <d> <force> <length>\n" );
        return NULL;
    }

    uUnit  = UNIT_from(oAssocObject(aaArgs[0]));
    for ( i=0; i<ATOMSINTORSION; i++ )
        aaAtoms[i] = ATOM_from(oAssocObject(aaArgs[1+i]));
    dKp    = dODouble(oAssocObject(aaArgs[5]));
    dP0    = DEGTORAD * dODouble(oAssocObject(aaArgs[6]));
    dN     = dODouble(oAssocObject(aaArgs[7]));

    for ( i=0; i<ATOMSINTORSION; i++ ) {
        if ( !bContainerContainedBy( CONTAINER_from(aaAtoms[i]), CONTAINER_from(uUnit )) ) {
            VPFATALEXIT("%s: Atom#%d is not part of the UNIT\n", sCmd, i+1 );
            return NULL;
        }
    }

                /* Create the RESTRAINT and add it to the UNIT */

    rRest = rRestraintCreate();
    RestraintTorsionSet( rRest, aaAtoms[0], aaAtoms[1], aaAtoms[2], aaAtoms[3],
                         dKp, dP0, dN );
    UnitAddRestraint( uUnit, rRest );

    return NULL;
}



/*
 *      oCmd_delete
 *
 *      Author: Juno Krahn (2026)
 *
 *      Delete an object
 */
static void
DeleteObject(OBJEKT oObject) {
OBJEKT oParent, oElement;
LISTLOOP llElements;

    switch (iObjectType(oObject)) {
    case RESIDUEid:
    case ATOMid:
        oParent = OBJEKT_from(cContainerWithin(oObject));
        REF( oParent );  /* hold reference until transferred */
        if ( !bContainerRemove( CONTAINER_from(oParent), oObject ) ) {
            VPFATALEXIT("Could not remove %s\n", sContainerName(oObject) );
        }
        DEREF( oParent ); /* release our refernece */
        break;
    case LISTid:
        llElements = llListLoop(LIST_from(oObject));
        while ( (oElement=oListNext(&llElements)) ) {
            ASSOC aAssoc = NULL;
            if (iObjectType(oElement)==ASSOCid) {
                aAssoc = ASSOC_from(oElement);
                oElement = oAssocObject(aAssoc);
            }
            if ( !bObjectInClass( oElement, CONTAINERid ) ) {
                VPFATALEXIT("Cannot delete %s - not a 'container'\n"
                        "\ttype %s\n", 
                         (iObjectType(oElement)==ASSOCid ? sAssocName(aAssoc) : "object"),
                         sObjectType(oElement) );
                continue;
            }
            oParent = OBJEKT_from(cContainerWithin(oElement));
            REF( oParent );  /* hold reference until transferred */
            if ( !bContainerRemove( CONTAINER_from(oParent), oElement ) ) {
                VPFATALEXIT("Could not remove %s\n", sContainerName(oObject) );
            }
            DEREF( oParent ); /* release our refernece */
        }
        break;
    default:
        VPFATAL("Cannot remove object of type '%c'\n",iObjectType(oObject));
    }
}

OBJEKT
oCmd_delete( int iArgCount, ASSOC aaArgs[] )
{

    if ( !bCmdGoodArguments( "delete", iArgCount, aaArgs, "ural" ) ) {
      VPFATALDELAYEDEXIT("usage:  delete <unit/residue/atom/list>\n" );
      return NULL;
    }
    OBJEKT oObject = oAssocObject(aaArgs[0]);
    if (iObjectType(oObject)==UNITid)
        VariableRemove( sAssocName(aaArgs[0]) );
    else 
        DeleteObject(oObject);

    return NULL;
}




/*
 *      oCmd_deleteRestraint
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Remove a RESTRAINT bond, angle, or torsion depending
 *      on the number of atoms passed by the caller.
 */
OBJEKT
oCmd_deleteRestraint( int iArgCount, ASSOC aaArgs[] )
{
int             i, iAtomCount;
UNIT            uUnit;
ATOM            aaAtoms[4];
RESTRAINT       rRest;
bool            bGotOne;
bool            bFound;
char            *sCmd = "deleteRestraint";

    if ( iArgCount == 3 ) {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u a a" ) ) {
          VPFATALDELAYEDEXIT("usage:  deleteRestraint <unit> <a> <b> [<c> <d>]\n" );
          return NULL;
        }
        iAtomCount = 2;
    } else if ( iArgCount == 4 ) {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u a a a" ) ) {
          VPFATALDELAYEDEXIT("usage:  deleteRestraint <unit> <a> <b> [<c> <d>]\n" );
          return NULL;
        }
        iAtomCount = 3;
    } else {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u a a a a" ) ) {
          VPFATALDELAYEDEXIT("usage:  deleteRestraint <unit> <a> <b> [<c> <d>]\n" );
          return NULL;
        }
        iAtomCount = 4;
    }

        /* Get the UNIT and the ATOMs that form the RESTRAINT */
    uUnit  = UNIT_from(oAssocObject(aaArgs[0]));
    for ( i=0; i<iArgCount-1; i++ )
        aaAtoms[i] = ATOM_from(oAssocObject(aaArgs[1+i]));

                /* Now loop through all the RESTRAINTs and find */
                /* the one the user specified */
    bGotOne = false;
    UnitLoopRestraints( uUnit );
    while ( (rRest = rUnitNextRestraint(uUnit)) ) {
        switch ( iRestraintType(rRest) ) {
            case RESTRAINTBOND:
                if ( iAtomCount == 2 ) {
                    if ( bRestraintBondMatchAtoms( rRest, aaAtoms[0],
                                                   aaAtoms[1] ) ) {
                        bGotOne = true;
                        goto DONE;
                    }
                }
                break;
            case RESTRAINTANGLE:
                if ( iAtomCount == 3 ) {
                    if ( bRestraintAngleMatchAtoms( rRest, aaAtoms[0],
                                                    aaAtoms[1], aaAtoms[2] )){
                        bGotOne = true;
                        goto DONE;
                    }
                }
                break;
            case RESTRAINTTORSION:
                if ( iAtomCount == 4 ) {
                    if ( bRestraintTorsionMatchAtoms( rRest, aaAtoms[0],
                                                aaAtoms[1], aaAtoms[2],
                                                aaAtoms[3] ) ) {
                        bGotOne = true;
                        goto DONE;
                    }
                }
                break;
            default:
                DFATAL("%s: Invalid RESTRAINT!\n", sCmd );
        }
    }

DONE:
    if ( !bGotOne ) {
        VPFATALEXIT("%s: No such restraint could be found.\n", sCmd );
    } else {
        VP1("Removing restraint.\n" );
        bFound = bUnitRemoveRestraint( uUnit, rRest );
        VP1("Restraint was removed = %i \n", bFound );
    }
    return NULL;
}





/*
 *      oCmd_impose
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Impose the internal coordinates the user specified onto
 *      the RESIDUEs of the UNIT.
 *
 *      Arguments:
 *              [0]     - The UNIT to add change INTERNALs for.
 *              [1]     - The list of RESIDUEs to change.
 *              [2]     - The second ATOM of the RESTRAINT bond.
 *              [3]     - The third ATOM of the RESTRAINT bond.
 *              [4]     - The fourth ATOM of the RESTRAINT bond.
 *              [5]     - ODOUBLE: The value of Kp.
 *              [6]     - ODOUBLE: The value of P0.
 *              [7]     - ODOUBLE: The value of N.
 */
OBJEKT
oCmd_impose( int iArgCount, ASSOC aaArgs[] )
{
LIST            lResidues, lInternals;
OBJEKT          oOne;
LISTLOOP        llResidues, llOne;
RESIDUE         rRes;
#define         MAXOBJEKTSININTERNALLIST 101
OBJEKT          oaIntObjekts[MAXOBJEKTSININTERNALLIST];
int             iObjekts, iStart, iDum;
RESIDUE         rResModify;
LOOP            lAtoms, lSpanning;
UNIT            uUnit;
ASSOC           aInternal, aOne;
LISTLOOP        llInternals;
bool            bSkipSubList;
LIST            lOne;
char            *sCmd = "impose";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u l l" ) ) {
        VPFATALDELAYEDEXIT("usage:  %s <unit> <residueseqlist> <internalslistlist>\n",
                        sCmd );
        return NULL;
    }

    DisplayerAccumulateUpdates();

        /* Build the INTERNALs */

    uUnit = UNIT_from(oAssocObject( aaArgs[0] ));
    lAtoms = lLoop( OBJEKT_from(uUnit), ATOMS );
    BuildInternalsUsingFlags( &lAtoms, 0, 0, 0, ATOMPOSITIONKNOWN );

        /* First get a list of RESIDUEs using oCmd_Residues */

    lResidues = lToolListOfResidues( uUnit, LIST_from(oAssocObject(aaArgs[1])) );

        /* Now apply the INTERNALs to each of the RESIDUEs in turn */

    llResidues = llListLoop(lResidues);
    while ( (rRes = RESIDUE_from(oListNext(&llResidues))) ) {

                /* Loop over all the sub-lists */

        lInternals = LIST_from(oAssocObject(aaArgs[2]));
        llInternals = llListLoop(lInternals);
        while ( (aInternal = ASSOC_from(oListNext(&llInternals))) ) {

            if ( iObjectType(oAssocObject(aInternal)) != LISTid ) {
                VPFATALEXIT("%s: Invalid internal list in internalslistlist !\n"
                      "        Note that argument #3 is a list of lists.\n"
                      "        Here is an example:\n"
                      "        impose x {167 168} { { \"CH3\" $C $N $CA 180.0"
                      " } }\n",
                      sCmd );
                continue;
            }

                        /* Copy what is in the sub-list into an array */

            bSkipSubList = false;
            iObjekts = 0;
            lOne = LIST_from(oAssocObject(aInternal));
            llOne = llListLoop(lOne);
            while( (aOne = ASSOC_from(oListNext(&llOne))) ) {
                oOne = OBJEKT_from(oAssocObject(aOne));
                oaIntObjekts[iObjekts++] = oOne;
                if ( iObjekts > MAXOBJEKTSININTERNALLIST ) {
                    VPFATALEXIT("%s: Too many lists in argument #3\n", sCmd );
                    bSkipSubList = true;
                    break;
                }
            }

            if ( bSkipSubList ) continue;

                        /* Now apply the INTERNAL onto the RESIDUE */
            iStart = 0;
            if ( iObjectType(oaIntObjekts[0]) == ODOUBLEid ) iStart = 1;

            switch ( iObjekts - iStart ) {
                        /* Do a bond INTERNAL */
                case 3:
                    MESSAGE("Imposing a bond INTERNAL\n" );
                    if ( iObjectType(oaIntObjekts[iStart+0]) != OSTRINGid ||
                         iObjectType(oaIntObjekts[iStart+1]) != OSTRINGid ||
                         iObjectType(oaIntObjekts[iStart+2]) != ODOUBLEid ) {
                        VPFATALEXIT("%s: Invalid bond internal definition\n",
                                sCmd );
                    } else {
                        if ( iStart == 1 ) {
                           rResModify = rResidueConnected( rRes,
                                        (int)dODouble(oaIntObjekts[0]) );
                        } else rResModify = rRes;
                        MESSAGE("Imposing on Residue: %s:%d\n",
                                sContainerName(CONTAINER_from( rResModify)),
                                iContainerSequence(CONTAINER_from( rResModify)) );
                        bBuildChangeInternalBond( CONTAINER_from( rResModify),
                                sOString(oaIntObjekts[iStart+0]),
                                sOString(oaIntObjekts[iStart+1]),
                                dODouble(oaIntObjekts[iStart+2]) );
                    }
                    break;

                        /* Do an angle INTERNAL */
                case 4:
                    MESSAGE("Imposing an angle INTERNAL\n" );
                    if ( iObjectType(oaIntObjekts[iStart+0]) != OSTRINGid ||
                         iObjectType(oaIntObjekts[iStart+1]) != OSTRINGid ||
                         iObjectType(oaIntObjekts[iStart+2]) != OSTRINGid ||
                         iObjectType(oaIntObjekts[iStart+3]) != ODOUBLEid ) {
                        VPFATALEXIT("%s: Invalid angle internal definition\n",
                                sCmd );
                    } else {
                        if ( iStart == 1 ) {
                           rResModify = rResidueConnected( rRes,
                                        (int)dODouble(oaIntObjekts[0]) );
                        } else rResModify = rRes;
                        MESSAGE("Imposing on Residue: %s:%d\n",
                                sContainerName(CONTAINER_from( rResModify)),
                                iContainerSequence(CONTAINER_from( rResModify)) );
                        bBuildChangeInternalAngle( CONTAINER_from( rResModify),
                                sOString(oaIntObjekts[iStart+0]),
                                sOString(oaIntObjekts[iStart+1]),
                                sOString(oaIntObjekts[iStart+2]),
                                dODouble(oaIntObjekts[iStart+3])*DEGTORAD );
                    }
                    break;

                        /* Do a torsion INTERNAL */
                case 5:
                    MESSAGE("Imposing a torsion INTERNAL\n" );
                    if ( iObjectType(oaIntObjekts[iStart+0]) != OSTRINGid ||
                         iObjectType(oaIntObjekts[iStart+1]) != OSTRINGid ||
                         iObjectType(oaIntObjekts[iStart+2]) != OSTRINGid ||
                         iObjectType(oaIntObjekts[iStart+3]) != OSTRINGid ||
                         iObjectType(oaIntObjekts[iStart+4]) != ODOUBLEid ) {
                        VPFATALEXIT("%s: Invalid angle internal definition\n",
                                sCmd );
                    } else {
                        if ( iStart == 1 ) {
                           rResModify = rResidueConnected( rRes,
                                        (int)dODouble(oaIntObjekts[0]) );
                        } else rResModify = rRes;
                        MESSAGE("Imposing on Residue: %s:%d\n",
                                sContainerName(CONTAINER_from( rResModify)),
                                iContainerSequence(CONTAINER_from( rResModify)) );
                        bBuildChangeInternalTorsion( CONTAINER_from( rResModify),
                                sOString(oaIntObjekts[iStart+0]),
                                sOString(oaIntObjekts[iStart+1]),
                                sOString(oaIntObjekts[iStart+2]),
                                sOString(oaIntObjekts[iStart+3]),
                                dODouble(oaIntObjekts[iStart+4])*DEGTORAD );
                    }
                    break;
                default:
                    VPFATALEXIT("%s: Improper internal definition\n", sCmd );
                    break;
            }
        }
    }

                /* Now build the externals again */

    lAtoms = lLoop( OBJEKT_from(uUnit), ATOMS );
    lSpanning = lLoop( oNext(&lAtoms), SPANNINGTREE );
    iDum = 0;   /* for purify */
    BuildExternalsUsingFlags( &lSpanning, 0, 0,
                                ATOMPOSITIONKNOWN|ATOMPOSITIONBUILT,
                                0,
                                &iDum, &iDum, &iDum, false );

                /* Destroy all of the INTERNALs */

    lAtoms = lLoop( OBJEKT_from(uUnit), ATOMS );
    BuildDestroyInternals( &lAtoms );

    DisplayerReleaseUpdates();

    return NULL;
}




/*
 *      oCmd_translate
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Translate all of the ATOMs within the CONTAINER by
 *      the vector in [1].
 *
 *      Arguments:
 *              [0]     - The CONTAINER whose atoms to translate.
 *              [1]     - A list with the three coordinates of the vector.
 */
OBJEKT
oCmd_translate( int iArgCount, ASSOC aaArgs[] )
{
LIST            lVector;
LISTLOOP        llElements;
CONTAINER       cCont;
double          daVector[3];
int             i;
VECTOR          vOffset;
ASSOC           aAssoc;
char            *sCmd = "translate";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "umral l" ) ) {
        VPFATALDELAYEDEXIT("usage:  %s <unit/residue/atom> <directionlist>\n",
                                sCmd );
        return NULL;
    }

    DisplayerAccumulateUpdates();

                /* Get the CONTAINER to translate */

    cCont = CONTAINER_from(oAssocObject( aaArgs[0] ));

    lVector = LIST_from(oAssocObject( aaArgs[1] ));
    i = 0;
    llElements = llListLoop(lVector);
    while ( (aAssoc = ASSOC_from(oListNext(&llElements))) ) {
        if ( iObjectType(oAssocObject(aAssoc)) == ODOUBLEid ) {
            if ( i<3 ) {
                daVector[i] = dODouble(oAssocObject(aAssoc));
                i++;
            } else {
                VPFATALEXIT("%s: Invalid vector\n", sCmd );
                break;
            }
        } else {
            VPFATALEXIT("%s: Invalid vector\n", sCmd );
            break;
        }
    }

    if ( i==3 ) {
        VectorDef( &vOffset, daVector[0], daVector[1], daVector[2] );
        ContainerTranslateBy( cCont, vOffset );
    }

    DisplayerReleaseUpdates();

    return NULL;
}






/*
 *      oCmd_center
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Display the geometric center of the ATOMs in the CONTAINER.
 *
 *      Arguments:
 *              [0]     - The CONTAINER to use to find the geometric
 *                        center of ATOMs.
 */
OBJEKT
oCmd_center( int iArgCount, ASSOC aaArgs[] )
{
CONTAINER       cCont;
VECTOR          vOffset;

    if ( !bCmdGoodArguments( "center", iArgCount, aaArgs, "umral" ) ) {
          VPFATALDELAYEDEXIT("usage:  center <unit/residue/atom>\n" );
          return NULL;
    }

                /* Get the CONTAINER to translate */

    cCont = CONTAINER_from(oAssocObject( aaArgs[0] ));

    vOffset = vContainerGeometricCenter(cCont);

    VP0("The center is at: %4.2lf, %4.2lf, %4.2lf\n",
                dVX(&vOffset), dVY(&vOffset), dVZ(&vOffset) );

    return NULL;
}



/*
 *      oCmd_solvateCap
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Arguments:
 *              [0] -   UNIT to add cap of solvent to.
 *              [1] -   UNIT of solvent.
 *              [2] -   RESIDUE, MOLECULE, ATOM, LIST of atoms, or
 *                              LIST of 3 doubles to use to calculate
 *                              the point around which to place the cap.
 *              [3] -   ODOUBLE, radius of cap.
 *      Option  [4] -   ODOUBLE with closeness parameter.
 *
 */
OBJEKT
oCmd_solvateCap( int iArgCount, ASSOC aaArgs[] )
{
UNIT            uSolute;
UNIT            uSolvent;
OBJEKT          oPosition;
VECTOR          vPos;
double          dRadius, dCloseness = 1.0;
int             iInitialSize, iFinalSize;
char            *sCmd = "solvateCap";

    if ( iArgCount == 4 ) {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u mral n" ) ) {
            VPFATALDELAYEDEXIT("usage:  solvateCap <solute> <solvent> <position>"
                    " <radius> <closeness>\n" );
            return NULL;
        }
        dCloseness = 1.0;
    } else {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u mral n n" ) ) {
            VPFATALDELAYEDEXIT("usage:  solvateCap <solute> <solvent> <position>"
                    " <radius> <closeness>\n" );
            return NULL;
        }
        dCloseness = dODouble(oAssocObject(aaArgs[4]));
    }

    uSolute = UNIT_from(oAssocObject(aaArgs[0]));
    uSolvent = UNIT_from(oAssocObject(aaArgs[1]));
    /*
     *  make copy of solvent w/ box & solv residue types set
     */
    uSolvent = zToolSetupSolvent( uSolvent );
    oPosition = oAssocObject(aaArgs[2]);

                /* Check if the oPosition argument is an array of three */
                /* doubles, if it is then treat it like a vector */

    if ( !bToolGeometricCenter( oPosition, &vPos ) ) {
        return NULL;
    }


    dRadius = dODouble(oAssocObject(aaArgs[3]));
    if ( dRadius < 0.0 ) {
        VPFATALEXIT("radius (%f) must be > 0\n", dRadius );
        return NULL;
    }

    TurnOffDisplayerUpdates();
    iInitialSize = iContainerNumberOfChildren( CONTAINER_from( uSolute ));
    zToolSolvateInSphere( uSolute, uSolvent, &vPos, dRadius, dCloseness );

    iFinalSize = iContainerNumberOfChildren( CONTAINER_from( uSolute ));

    VP0("Added %d residues.\n", iFinalSize - iInitialSize );

                /* Define the solvent cap within the UNIT */

    UnitSetSolventCap( uSolute, dVX(&vPos), dVY(&vPos), dVZ(&vPos), dRadius );
    UnitSetUseSolventCap( uSolute, true );

    TurnOnDisplayerUpdates();
    ContainerDisplayerUpdate( CONTAINER_from( uSolute ));

    return NULL;
}





/*
 *      oCmd_addPath
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Arguments:
 *              [0] -   OSTRING, directory to add to file search
 *                      path list.
 *
 */
OBJEKT
oCmd_addPath( int iArgCount, ASSOC aaArgs[] )
{
    if ( !bCmdGoodArguments( "addPath", iArgCount, aaArgs, "s" ) ) {
        VPFATALDELAYEDEXIT("usage:  addPath <path>\n" );
        return NULL;
    }

    if ( BasicsAddDirectory( sOString(oAssocObject(aaArgs[0])), 0 ) )
        VP0("%s added to file search path.\n",
                        sOString(oAssocObject(aaArgs[0])) );

    return NULL;
}




/*
 *      oCmd_crossLink
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Create a crosslink between two RESIDUEs
 *      according to the named connect ATOMs.
 *
 *      Arguments:
 *              [0] -   RESIDUE
 *              [1] -   OSTRING, name of connect ATOM on [0]
 *              [2] -   RESIDUE
 *              [3] -   OSTRING, name of connect ATOM on [2]
 *      option  [4] -   OSTRING, bond order, see BOND command.
 */
OBJEKT
oCmd_crossLink( int iArgCount, ASSOC aaArgs[] )
{
int             iConnectA, iConnectB;
RESIDUE         rA, rB;
OBJEKT          oObj;
int             iOrder;
char            *sCmd = "crossLink";
char            *usage =
        "usage:  crossLink <res1> <connect> <res2> <connect> [bondorder]\n";

    if ( iArgCount == 5 ) {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "r s r s s" ) ) {
          VP0(usage );
          VP0("(For multiple molecules, note that residue numbering jumps\n");
          VPFATALDELAYEDEXIT(" by 1001 for each new molecule)\n");
          return NULL;
        }
        oObj = oAssocObject(aaArgs[4]);
        iOrder = iAtomBondOrderFromName(sOString(oObj));
        if ( iOrder == BONDNONE ) {
            VPFATALEXIT("%s: Invalid bond order\n", sCmd );
            return NULL;
        }
    } else {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "r s r s" ) ) {
          VPFATALDELAYEDEXIT("usage:  crossLink <res1> <connect> <res2> "
                                                 "<connect> [bondorder]\n" );
          return NULL;
        }
        iOrder = BONDSINGLE;
    }

    oObj = oAssocObject(aaArgs[1]);
    iConnectA = iResidueConnectFromName(sOString(oObj));
    if ( iConnectA == NOEND ) {
        VPFATALEXIT("%s: Invalid connect atom: %s\n", sCmd, sOString(oObj) );
        return NULL;
    }
    oObj = oAssocObject(aaArgs[3]);
    iConnectB = iResidueConnectFromName(sOString(oObj));
    if ( iConnectB == NOEND ) {
        VPFATALEXIT("%s: Invalid connect atom: %s\n", sCmd, sOString(oObj) );
        return NULL;
    }

    rA = RESIDUE_from(oAssocObject(aaArgs[0]));
    rB = RESIDUE_from(oAssocObject(aaArgs[2]));

    DisplayerAccumulateUpdates();

    if ( !bResidueCrossLink( rA, iConnectA, rB, iConnectB, iOrder ) ) {
        VPFATALEXIT("%s: Could not form cross link, invalid connection atom\n",
                                sCmd );
    }

    DisplayerReleaseUpdates();

    return NULL;
}


/*
 *      oCmd_addPdbResMap
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Append.
 *
 *      Arguments:
 *              [0] -   LIST of entries to add.
 *
 *      The Name Map is used to map residue names read from
 *      PDB files to variable names within LEaP.
 *      The LIST is a LIST of LISTs each sub-list containing
 *      two or three entries.  Each entry has the form:
 *
 *      { ODOUBLE OSTRING OSTRING }
 *      { OSTRING OSTRING }
 *
 *      ODOUBLE can be 0 or 1
 *      The first OSTRING is the name within the PDB file.
 *      The second OSTRING is the variable name to map to.
 */
OBJEKT
oCmd_addPdbResMap( int iArgCount, ASSOC aaArgs[] )
{
    if ( !bCmdGoodArguments( "addPdbResMap", iArgCount, aaArgs, "l" ) ) {
         VPFATALDELAYEDEXIT("usage:  addPdbResMap <list_of_lists>\n" );
         return NULL;
    }
    PdbAppendToResMap( LIST_from(oAssocObject(aaArgs[0])) );
    return NULL;
}

/*
 *      oCmd_addPdbAtomMap
 *
 *      Author: Bill Ross
 */
OBJEKT
oCmd_addPdbAtomMap( int iArgCount, ASSOC aaArgs[] )
{

    if ( !bCmdGoodArguments( "addPdbMap", iArgCount, aaArgs, "l" ) ) {
         VPFATALDELAYEDEXIT("usage:  addPdbAtomMap <list_of_lists>\n" );
         return NULL;
    }

    PdbAppendToAtomMap( LIST_from(oAssocObject(aaArgs[0])) );

    return NULL;
}

/*
 *      oCmd_addAtomTypes
 *
 *      Author: Bill Ross (1996)
 */
OBJEKT
oCmd_addAtomTypes( int iArgCount, ASSOC aaArgs[] )
{
LIST            lList;
        if ( !bCmdGoodArguments( "addAtomTypes", iArgCount, aaArgs, "l" ) ) {
                VPFATALDELAYEDEXIT("usage:  addAtomTypes <list_of_lists>\n" );
                return NULL;
        }

        lList =  LIST_from(oAssocObject(aaArgs[0]));

        AmberAddAtomTypes( lList );
        return NULL;
}



/*
 *      oCmd_displayPdbResMap
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Display the Name Map.
 */
OBJEKT
oCmd_displayPdbResMap( int iArgCount, ASSOC aaArgs[] )
{
    if ( !bCmdGoodArguments( "displayPdbResMap", iArgCount, aaArgs, "" ) ) {
         VPFATALDELAYEDEXIT("usage:  displayPdbResMap\n" );
         return NULL;
    }
    PdbDisplayResMap();
    return NULL;
}

OBJEKT
oCmd_displayPdbAtomMap( int iArgCount, ASSOC aaArgs[] )
{
    if ( !bCmdGoodArguments( "displayPdbAtomMap", iArgCount, aaArgs, "" ) ) {
         VPFATALDELAYEDEXIT("usage:  displayPdbAtomMap\n" );
         return NULL;
    }
    PdbDisplayAtomMap();
    return NULL;
}





/*
 *      oCmd_clearPdbResMap
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Clear the Name Map.
 *
 *
 */
OBJEKT
oCmd_clearPdbResMap( int iArgCount, ASSOC aaArgs[] )
{
    if ( !bCmdGoodArguments( "clearPdbResMap", iArgCount, aaArgs, "" ) ) {
         VPFATALDELAYEDEXIT("usage:  clearPdbResMap\n" );
         return NULL;
    }
    PdbClearResMap();
    return NULL;
}

OBJEKT
oCmd_clearPdbAtomMap( int iArgCount, ASSOC aaArgs[] )
{

    if ( !bCmdGoodArguments( "clearPdbAtomMap", iArgCount, aaArgs, "" ) ) {
         VPFATALDELAYEDEXIT("usage:  clearPdbAtomMap\n" );
         return NULL;
    }
    PdbClearAtomMap();
    return NULL;
}





/*
 *      oCmd_measureGeom
 *
 *      Author: Christian Schafmeister (1991)
 *
 *              Measure the distance, angle, torsion between
 *              two, three, or four ATOMs.
 *
 *      Arguments:
 *              [0] -   ATOM
 *              [1] -   ATOM
 *      option  [2] -   ATOM
 *      option  [3] -   ATOM
 *
 *      Return:
 *              ODOUBLE, return the internal coordinates value.
 */
OBJEKT
oCmd_measureGeom( int iArgCount, ASSOC aaArgs[] )
{
ATOM            aA, aB, aC, aD;
ODOUBLE         odVal;
double          dVal;
char            *sCmd = "measureGeom";

    if ( iArgCount == 4 ) {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "a a a a" ) ) {
          VPFATALDELAYEDEXIT("usage:  measureGeom <atom> <atom> [atom [atom]]\n" );
          return NULL;
        }
        aA = ATOM_from(oAssocObject(aaArgs[0]));
        aB = ATOM_from(oAssocObject(aaArgs[1]));
        aC = ATOM_from(oAssocObject(aaArgs[2]));
        aD = ATOM_from(oAssocObject(aaArgs[3]));
        dVal = dVectorAtomTorsion(
                &vAtomPosition(aA), &vAtomPosition(aB),
                &vAtomPosition(aC), &vAtomPosition(aD) )*RADTODEG;
        VP0("Torsion angle: %4.2lf degrees\n", dVal );
    } else if ( iArgCount == 3 ) {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "a a a" ) ) {
          VPFATALDELAYEDEXIT("usage:  measureGeom <atom> <atom> [atom [atom]]\n" );
          return NULL;
        }
        aA = ATOM_from(oAssocObject(aaArgs[0]));
        aB = ATOM_from(oAssocObject(aaArgs[1]));
        aC = ATOM_from(oAssocObject(aaArgs[2]));
        dVal = dVectorAtomAngle(
                &vAtomPosition(aA), &vAtomPosition(aB),
                &vAtomPosition(aC) )*RADTODEG;
        VP0("Angle: %4.2lf degrees\n", dVal );
    } else {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "a a" ) ) {
          VPFATALDELAYEDEXIT("usage:  measureGeom <atom> <atom> [atom [atom]]\n" );
          return NULL;
        }
        aA = ATOM_from(oAssocObject(aaArgs[0]));
        aB = ATOM_from(oAssocObject(aaArgs[1]));
        dVal = dVectorAtomLength(
                &vAtomPosition(aA), &vAtomPosition(aB) );
        VP0("Distance: %4.2lf angstroms\n", dVal );
    }
    odVal = ODOUBLE_from(oCreate(ODOUBLEid));
    ODoubleSet( odVal, dVal );
    return OBJEKT_from(odVal);
}


OBJEKT
oCmd_memDebug( int iArgCount, ASSOC aaArgs[] )
{
        iMemDebug = 1;
        return NULL;
}

/*
 *      oCmd_bondByDistance
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Perform an N^2 search of all ATOMs within the
 *      container and create single bonds between all ATOMs
 *      that are within a certain distance of each other.
 *
 *      Arguments:
 *              [0] -   UNIT/MOLECULE/RESIDUE/LIST
 *      option  [1] -   ODOUBLE  maximum bonding distance
 */
OBJEKT
oCmd_bondByDistance( int iArgCount, ASSOC aaArgs[] )
{
double          dDist;
int             iBonds;
char            *sCmd = "bondByDistance";
char            *sUsage = "usage:  bondByDistance <unit> [maxdistance]\n";

    if ( iArgCount == 1 ) {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "umrl" ) ) {
            VPFATALDELAYEDEXIT("%s",sUsage );
            return NULL;
        }
        if (!(dDist = GDefaults.dDSearchDistance))
            dDist = DEFAULT_DISTANCE_SEARCH;
    } else {
        if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "umrl n" ) ) {
            VPFATALDELAYEDEXIT("%s",sUsage );
            return NULL;
        }
        dDist = dODouble(oAssocObject(aaArgs[1]));
    }

    TurnOffDisplayerUpdates();

    iBonds = iToolDistanceSearch( CONTAINER_from(oAssocObject(aaArgs[0])), dDist,
                                        true, DISTANCE_SEARCH_CREATE_BONDS );

    VP0("Created %d bonds.\n", iBonds );

    TurnOnDisplayerUpdates();
    ContainerDisplayerUpdate( CONTAINER_from( oAssocObject(aaArgs[0])));

    return NULL;
}



/*
 *      oCmd_groupSelectedAtoms
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Create a group within the UNIT using all
 *      of the ATOMs that have ATOMSELECTED flag set.
 *
 *              [0] -   UNIT
 *              [1] -   OSTRING  name of group
 */
OBJEKT
oCmd_groupSelectedAtoms( int iArgCount, ASSOC aaArgs[] )
{
char            *sGroup;
UNIT            uUnit;
LOOP            lAtoms;
ATOM            aAtom;
int             iAtoms;
char            *sCmd = "groupSelectedAtoms";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u s" ) ) {
          VPFATALDELAYEDEXIT("usage:  groupSelectedAtoms <unit> <groupname>\n" );
          return NULL;
    }

    DisplayerAccumulateUpdates();

    uUnit = UNIT_from(oAssocObject(aaArgs[0]));
    sGroup = sOString(oAssocObject(aaArgs[1]));

    if ( lUnitGroup( uUnit, sGroup ) ) {
        bUnitGroupDestroy( uUnit, sGroup );
    }

    bUnitGroupCreate( uUnit, sGroup );

    iAtoms = 0;
    lAtoms = lLoop( OBJEKT_from(uUnit), ATOMS );
    while ( (aAtom = ATOM_from(oNext(&lAtoms))) ) {
        if ( bAtomFlagsSet( aAtom, ATOMSELECTED ) ) {
            bUnitGroupAddAtom( uUnit, sGroup, aAtom );
            iAtoms++;
        }
    }
    VP0("Group %s contains %d atoms.\n", sGroup, iAtoms );

    DisplayerReleaseUpdates();

    return NULL;
}







/*
 *      oCmd_transform
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Translate all of the ATOMs within the CONTAINER by
 *      the vector in [1].
 *
 *      Arguments:
 *              [0]     - The CONTAINER whose atoms to transform.
 *              [1]     - A list of lists representing a 4x4 matrix.
 */
OBJEKT
oCmd_transform( int iArgCount, ASSOC aaArgs[] )
{
LIST            lVectorX;
LIST            lVectorY;
LISTLOOP        llElementsX;
LISTLOOP        llElementsY;
CONTAINER       cCont;
MATRIX          mTransform;
int             iX, iY;
ASSOC           aAssocX;
ASSOC           aAssocY;
char            *sCmd = "transform";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "umral l" ) ) {
          VPFATALDELAYEDEXIT("usage:  transform <atoms> <matrixlist>\n" );
          return NULL;
    }

                /* Get the CONTAINER to translate */

    cCont = CONTAINER_from(oAssocObject( aaArgs[0] ));

    lVectorY = LIST_from(oAssocObject( aaArgs[1] ));
    iX = 0;
    iY = 0;
    MatrixIdentity( mTransform );
    llElementsY = llListLoop(lVectorY);
    while ( (aAssocY = ASSOC_from(oListNext(&llElementsY))) ) {
        if ( iObjectType(oAssocObject(aAssocY)) == LISTid ) {
            if ( iY<4 ) {
                lVectorX = LIST_from(oAssocObject(aAssocY));
                llElementsX = llListLoop(lVectorX);
                while ( (aAssocX = ASSOC_from(oListNext(&llElementsX))) ) {
                    if ( iObjectType(oAssocObject(aAssocX)) == ODOUBLEid ) {
                        if ( iX<4 ) {
                            mTransform[iX][iY] = dODouble(oAssocObject(aAssocX));
                            iX++;
                        } else { goto ERROR; }
                    } else { goto ERROR; }
                }
                iX = 0;
                iY++;
            } else { goto ERROR; }
        } else { goto ERROR; }
    }

    DisplayerAccumulateUpdates();

    ContainerTransformBy( cCont, mTransform );

    DisplayerReleaseUpdates();

    return NULL;

ERROR:

    VPFATALEXIT("%s: Invalid matrix\n", sCmd );
    return NULL;
}



/*
 *      oCmd_copy
 *
 *      Author: Christian Schafmeister (1991)
 *
 */
OBJEKT
oCmd_copy( int iArgCount, ASSOC aaArgs[] )
{
OBJEKT          oNew;

    if ( !bCmdGoodArguments( "copy", iArgCount, aaArgs, "umranp" ) ) {
          VPFATALDELAYEDEXIT("usage:  <newvariable> = copy <variable>\n" );
          return NULL;
    }

    oNew = oCopy(oAssocObject(aaArgs[0]));

    return oNew;
}


/*
 *      oCmd_listParmSets
 *
 *      Author: Juno Krahn
 *
 *      List all ParmSets
 */
OBJEKT
oCmd_listParmSets( int iArgCount, ASSOC aaArgs[] )
{
PARMSET         psSet;

    if ( !bCmdGoodArguments( "listParmSets", iArgCount, aaArgs, "" ) ) {
          VPFATALDELAYEDEXIT("usage:  listParmSets\n" );
          return NULL;
    }

    VPTRACEENTER(__func__ );
    int sumAtoms=0, sumBonds=0, sumAngles=0;
    int sumTorsions=0, sumImpropers=0, sumHBonds=0, sumNBEdits=0;
    int nParmSets=0;
    // XXX: what is the differnece between GplAllParameters and GplDefaultParmLib ?
    ParmLibParmSetLoop(GplAllParameters);
    while ( bParmLibNextParmSet(GplAllParameters, &psSet) ) {
        int nAtoms =     iVarArrayElementCount(psSet->vaAtoms);
        int nBonds =     iVarArrayElementCount(psSet->vaBonds);
        int nAngles =    iVarArrayElementCount(psSet->vaAngles);
        int nTorsions =  iVarArrayElementCount(psSet->vaTorsions);
        int nImpropers = iVarArrayElementCount(psSet->vaImpropers);
        int nHBonds =    iVarArrayElementCount(psSet->vaHBonds);
        int nNBEdits =   iVarArrayElementCount(psSet->vaNBEdits);
        nParmSets++;
        VP0("\nFilename: %s\n",psSet->sFname);
        VP0("Title: %s\n",psSet->sTitle);
        VP0("Atoms:%d Bonds:%d Angles:%d Tors:%d Impr:%d HBnds:%d NBFix:%d\n",
                nAtoms, nBonds, nAngles, nTorsions, nImpropers, nHBonds, nNBEdits);
        sumAtoms += nAtoms;
        sumBonds += nBonds;
        sumAngles += nAngles;
        sumTorsions += nTorsions;
        sumImpropers += nImpropers;
        sumHBonds += nHBonds;
        sumNBEdits += nNBEdits;
    }
    VP0("\nTotal ParmSets: %d\n"
         "Atoms:%d Bonds:%d Angles:%d Tors:%d Impr:%d HBnds:%d NBFix:%d\n",
                nParmSets, sumAtoms, sumBonds, sumAngles,
                sumTorsions, sumImpropers, sumHBonds, sumNBEdits);
    return NULL;
}

/*
 *      oCmd_listResidues
 *
 *      Author: Juno Krahn
 *
 *      List all template residues in the variable table
 */
OBJEKT
oCmd_listResidues( int iArgCount, ASSOC aaArgs[] )
{
DICTLOOP        dlLoop;
DICTIONARY      dVariables;
STRING          sErrors;
int             iCount=0, iErrorCount=0;

    if ( !bCmdGoodArguments( "listResidues", iArgCount, aaArgs, "" ) ) {
          VPFATALDELAYEDEXIT("usage:  listResidues\n" );
          return NULL;
    }

    VPTRACEENTER(__func__ );
    dVariables = dVariablesDictionary();
    dlLoop = ydlDictionaryLoop(dVariables);
    int iUpdatedElements = 0;
    while ( yPDictionaryNext(dVariables, &dlLoop ) ) {
        OBJEKT oObj = PDictLoopData(dlLoop);
        if ( iObjectType(oObj) != UNITid ) continue;
        UNIT uUnit = UNIT_from(oObj);
        if ( iContainerNumberOfChildren(uUnit) > 1 ) continue;
        RESIDUE rRes = RESIDUE_from(oContainerFirstObject(uUnit));
        if ( iObjectType(rRes) != RESIDUEid ) continue; // unlikely
        // Variable is a UNIT containing a single RESIDUE
        if (iCount % 40 == 0) {
            VP0("--------+-------+---+--------+---------+-------+---+----+----+----+----+----+----+---------\n");
            VP0("   NAME   UNIT   End HEAD(el) TAIL(el)  RESIDUE Typ  C0   C1   C2   C3   C4   C5   Flags\n");
            VP0("--------+-------+---+--------+---------+-------+---+----+----+----+----+----+----+---------\n");
        }
        iCount++;
        char *cPVarName = sDictLoopKey(dlLoop);
        char *cPUnitName = sContainerName(uUnit);
        char *cPResName = sContainerName(rRes);
        int fEndFlag = fGetPdbResMapped(cPUnitName);
        char sEndFlags[4]="";
        if (fEndFlag & RESIDUEFIRSTEND) strcat(sEndFlags,"F");
        if (fEndFlag & RESIDUENOEND) strcat(sEndFlags,"N");
        if (fEndFlag & RESIDUELASTEND) strcat(sEndFlags,"L");
        ATOM aHead = ATOM_from(uUnit->aHead);
        ATOM aTail = ATOM_from(uUnit->aTail);
        ATOM aAtom;
        bool bNonAtom=false, bLongAtomName=false;
        bool bMissingElement=false, bUpdatedElement=false;
        LOOP lContents = lLoop( OBJEKT_from(rRes), DIRECTCONTENTSBYSEQNUM );
        while ( (aAtom = ATOM_from(oNext(&lContents))) ) {
            char *cPAtomName = sContainerName(aAtom);
            if ( iObjectType(aAtom) != ATOMid ) {
                 bNonAtom = true;
                 break;
            }
            if ( iAtomElement(aAtom) < 0 ) {
                 iUpdatedElements++;
                 AtomSetElement(aAtom,iElementNumberFromAmber(cPAtomName));
                 if ( iAtomElement(aAtom) < 0 ) bMissingElement=true;
                 else bUpdatedElement=true;
            }
            if ( strlen(sContainerName(aAtom))>4) bLongAtomName=true;
        }
        if (bNonAtom) sprintf(sErrors,"Contains non-Atoms(%c), ",iObjectType(uUnit));
        else sErrors[0]=0;
        bool bMissingConnect01 = iContainerNumberOfChildren(rRes) > 1 &&
                   !(bResidueConnectUsed(rRes,0) && bResidueConnectUsed(rRes,1));
        if (bMissingConnect01) strcat(sErrors,"Missing Connect01, ");
        if (strlen(sEndFlags)>1) strcat(sErrors,"Mixed END flags, ");
        if ( (fEndFlag & RESIDUEFIRSTEND && (aHead || !aTail)) ||
                 (fEndFlag & RESIDUELASTEND && (!aHead || aTail)) ||
                 (fEndFlag & RESIDUENOEND && (!aHead) != (!aTail) ) ) {
             strcat(sErrors,"PdbMap mismatch, ");
        }
        if (bUpdatedElement) strcat(sErrors,"Updated elem, ");
        if (bMissingElement) strcat(sErrors,"Undefined elem, ");
        if (bLongAtomName) strcat(sErrors,"name>4, ");
        if (strlen(cPResName)>4 || (strlen(cPResName)==4 &&
                (cPResName[0] != 'N' &&cPResName[0] != 'C'))) strcat(sErrors,"ResName>3, ");
        if (cResidueType(rRes) == '?' && (aHead || aTail)) strcat(sErrors,"Undef. ResType, ");
        // design note: connect0, connect1 always defined for tree purposes, (is it really used??)
        //  but only unit head/tail used to link polymers
        if (aHead && aHead != aResidueConnectAtom(rRes,0)) strcat(sErrors,"Head mismatch, ");
        if (aTail && aTail != aResidueConnectAtom(rRes,1)) strcat(sErrors,"Tail mismatch, ");
        for (int i=0;i<MAXCONNECT;i++) {
            if ( aResidueConnectAtom(rRes,i) && iAtomElement(aResidueConnectAtom(rRes,i)) == HELIUM) {
                strcat(sErrors,"Links to Hydrogen, ");
                break;
            }
        }
        if (strlen(sErrors)>2) {
            sErrors[strlen(sErrors)-2]=0;
            iErrorCount++;
        }
        VP0("%8s=<U %4s;%3s,%4s(%2d),%4s(%2d)>,<R %4s;%3.3s;%4s,%4s,%4s,%4s,%4s,%4s> "
#ifdef DEBUG
               "(%d,%d) "
#endif
               "%s\n",
               cPVarName,cPUnitName,sEndFlags,
               (aHead?sContainerName(aHead):""),
               (aHead?iAtomElement(aHead):-1),
               (aTail?sContainerName(aTail):""),
               (aTail?iAtomElement(aTail):-1),
               cPResName,
               sResidueTypeNameFromChar(cResidueType(rRes)),
               bResidueConnectUsed(rRes,0) ? sContainerName(aResidueConnectAtom(rRes,0)):"",
               bResidueConnectUsed(rRes,1) ? sContainerName(aResidueConnectAtom(rRes,1)):"",
               bResidueConnectUsed(rRes,2) ? sContainerName(aResidueConnectAtom(rRes,2)):"",
               bResidueConnectUsed(rRes,3) ? sContainerName(aResidueConnectAtom(rRes,3)):"",
               bResidueConnectUsed(rRes,4) ? sContainerName(aResidueConnectAtom(rRes,4)):"",
               bResidueConnectUsed(rRes,5) ? sContainerName(aResidueConnectAtom(rRes,5)):"",
#ifdef DEBUG // reference counting is hopelelssly corrupted so don't even try
               OBJEKT_from(uUnit)->iReferences,
               OBJEKT_from(rRes)->iReferences,
#endif
               sErrors );
    }
    VP0("--------+-------+--------+---------+-------+---+----+----+----+----+----+----+---------\n");
    VP0("Found %d Residue Template variables and %d with errors or warnings\n",iCount,iErrorCount);
    if (iUpdatedElements) VP0("Updated %d missing Atom element defintions\n",iUpdatedElements);
    VPTRACEEXIT(__func__ );
    return NULL;
}

/*
 *      oCmd_listOff
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      List the contents of a library.
 */
OBJEKT
oCmd_listOff( int iArgCount, ASSOC aaArgs[] )
{
STRING          sFilename;
LIBRARY         lLib;
char            *cPNext;

    if ( !bCmdGoodArguments( "listOff", iArgCount, aaArgs, "s" ) ) {
          VPFATALDELAYEDEXIT("usage:  listOff <filename>\n" );
          return NULL;
    }
    strcpy( sFilename, sOString(oAssocObject(aaArgs[0])) );

    lLib = lLibraryOpen( sFilename, OPENREADONLY );
    if ( lLib == NULL ) return NULL;

    VP0("Index of library: %s\n", sFilename );

    LibraryLoop( lLib );
    while ( (cPNext = sLibraryNext(lLib)) ) {
        VP0("%s\n", cPNext );
    }

    LibraryClose( &lLib );

    return NULL;
}




/*
 *      oCmd_deleteOffLibEntry
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Remove an entry from a library.
 */
OBJEKT
oCmd_deleteOffLibEntry( int iArgCount, ASSOC aaArgs[] )
{
STRING          sFilename;
LIBRARY         lLib;
bool            bRemoved;
STRING          sEntry;
char            *sCmd = "deleteOffLibEntry";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "s s" ) ) {
          VPFATALDELAYEDEXIT("usage:  deleteOffLibEntry <filename> <entry>\n" );
          return NULL;
    }
    strcpy( sFilename, sOString(oAssocObject(aaArgs[0])) );
    strcpy( sEntry, sOString(oAssocObject(aaArgs[1])) );

    lLib = lLibraryOpen( sFilename, OPENREADWRITE );
    if ( lLib == NULL ) return NULL;

    bRemoved = bLibraryRemove( lLib, sEntry );

    if ( !bRemoved ) {
        VPFATALEXIT("%s: %s was not found.\n", sCmd, sEntry );
    } else {
        VP0("%s was removed.\n", sEntry );
    }

    LibraryClose(&lLib);

    return NULL;
}




/*
 *      oCmd_mutate
 *
 *      Author: Christian Schafmeister (1991)
 *
 */
OBJEKT
oCmd_mutate( int iArgCount, ASSOC aaArgs[] )
{
UNIT            uUnit;
ODOUBLE         odSeqNum;
RESIDUE         rNew;
RESIDUE         rOld;
RESIDUE         rCopy;
int             iSeqNum;
LOOP            lAtoms;
ATOM            aAtom;
ATOM            aNeighbor;
int             i, iNext;
char            *sCmd = "mutate";

    if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u n r" ) ) {
          VPFATALDELAYEDEXIT("usage:  mutate <unit> <number> <residue>\n" );
          return NULL;
    }
    uUnit = UNIT_from(oAssocObject( aaArgs[0] ));
    odSeqNum = ODOUBLE_from(oAssocObject( aaArgs[1] ));
    rNew = RESIDUE_from(oAssocObject( aaArgs[2] ));

    DisplayerAccumulateUpdates();

    iSeqNum = (int)dODouble(odSeqNum);

    rOld = RESIDUE_from(cContainerFindSequence( CONTAINER_from(uUnit),
                                                RESIDUEid, iSeqNum ));
    if ( rOld == NULL ) {
        VPFATAL("%s: Could not find residue with sequence number: %d\n",
                                        sCmd, iSeqNum );
        goto FAIL;
    }

                /* Make sure there are no bonds out of the new residue */

    lAtoms = lLoop( OBJEKT_from(rNew), ATOMS );
    while ( (aAtom = ATOM_from(oNext(&lAtoms))) ) {
        for ( i=0; i<iAtomCoordination(aAtom); i++ ) {
            aNeighbor = aAtomBondedNeighbor( aAtom, i );
            if ( rNew != RESIDUE_from(cContainerWithin(CONTAINER_from( aNeighbor)))) {
                VPFATAL("%s: The mutant residue cannot be bonded to anything.\n",
                                                sCmd );
                goto FAIL;
            }
        }
    }

    rCopy = rCopyResidue(rNew);

                /* Perform the mutation */

    ResidueMutate( rCopy, rOld );

                /* Remove the old RESIDUE from the UNIT */
                /* And add the new RESIDUE, without incrementing */
                /* the UNITs next child sequence number */

    REF( rOld );  /* hold temporary refenece */
    bContainerRemove( CONTAINER_from(uUnit), OBJEKT_from(rOld ));
    iNext = iContainerNextChildsSequence( CONTAINER_from( uUnit ));
    ContainerAdd( CONTAINER_from(uUnit), OBJEKT_from(rCopy ));
    ContainerSetSequence( CONTAINER_from( rCopy),
                        iContainerSequence(CONTAINER_from( rOld)) );
    ContainerSetNextChildsSequence( CONTAINER_from( uUnit), iNext );

    DEREF( rOld ); /* release temorary reference */
    goto RET;

FAIL:
    VPFATALDELAYEDEXIT("Mutation failed.\n" );

RET:
    DisplayerReleaseUpdates();
    return NULL;
}


/*
 *      oCmd_addIons
 *
 *      Author: Bill Ross (1993)
 *
 *      addIons unit, ion1, #ion1, [ion2, #ion2]
 *
 *              [0] -   UNIT
 *              [1] -   UNIT    ion to add
 *              [2] -   NUMBER  number of ions to add (0 = neutralize)
 *      Optional[3] -   UNIT    second ion to add
 *      Optional[4] -   NUMBER  number of the second ion to add
 *
 *      Adds Counter Ions.
 *      A Coulomb's Law ESP grid is used to calculate
 *      the appropriate location for the ion.
 *
 *      How it works:
 *
 *      1.  If [2] = 0, [0] must be charged and [1] must be opposite in charge;
 *              [0] is neutralized with [1].
 *      2.  If [2] != 0, [1] (and optionally [3]) are added (in alternation).
 *
 */

/*
 *  vaSolventResidues() - make array of residue pointers
 *      TODO - maybe restrict it to dShellExtent?
 */
VARARRAY        vaSolventResidues( uUnit )
UNIT            uUnit;
{
VARARRAY        vaSolvent;
LOOP            lRes;
RESIDUE         rRes;

        vaSolvent = vaVarArrayCreate( sizeof(RESIDUE) );
        lRes = lLoop( OBJEKT_from(uUnit), RESIDUES );
        while ( (rRes = RESIDUE_from(oNext( &lRes ))) ) {
            if ( cResidueType( rRes ) != RESTYPESOLVENT )
                continue;
            if (GDefaults.bKeepInputSolvent && !bResidueFlagsSet(rRes,RESIDUEBULKSOLVENT))
                continue;
            VarArrayAdd( vaSolvent, (GENP)&rRes );
        }
        if ( !iVarArrayElementCount( vaSolvent ) )
            VarArrayDestroy( &vaSolvent );
        return  vaSolvent ;
}

// Find closest vaSolvent atom to PvIon coordinate. If d < 3A,
// delete it from both uUnit (removed) and vaSolvent (destroy -> NULL, fast)
static void
CheckSolvent( UNIT uUnit, VARARRAY vaSolvent, UNIT uIon, VECTOR *PvIon, NeighborGrid *ngSolvent )
{
RESIDUE         *PrRes=NULL, *PrClosest=NULL;
VECTOR          vClosest;
double          dmin2;
        /*
         *  PrRes (pointer to residue pointer) is used so that
         *      the residue pointer can be set null if the residue
         *      is deleted
         */
        dmin2 = FLT_MAX;
        const Pair *Pairs;
        size_t count;
        neighbor_grid_query_point(ngSolvent,PvIon->dX,PvIon->dY,PvIon->dZ,1,1,&Pairs,&count);
        for (int i=0;i<count;i++) {
            if (Pairs[i].d2 < dmin2) {
                PrRes = PVAI( vaSolvent, RESIDUE, Pairs[i].to_group );
                if ( *PrRes ) {
                    dmin2 = Pairs[i].d2;
                    PrClosest = PrRes;
                    vClosest = vAtomPosition( (ATOM)Pairs[i].to_p );
                }
            }
        }
        if ( dmin2 < 9 ) {  /* HACK test 3A */
                VP0("(Replacing solvent molecule)\n");
                REF( *PrClosest );  /* hold temporary refernece */
                if ( bContainerRemove( CONTAINER_from(uUnit), OBJEKT_from(*PrClosest ))) {
                        ContainerDestroy( (CONTAINER *) PrClosest );
                        *PrClosest = NULL;
                        *PvIon = vClosest;
                } else
                        VPFATALEXIT("solvent removal failed\n" );
                DEREF( *PrClosest );  /* release our temporary reference */
        } else
                VP0("(No solvent overlap)\n");
        return;
}

OBJEKT
addIons( int iArgCount, ASSOC aaArgs[], bool bIncludeSolvent, char *sCmd)
{
UNIT            uUnit=NULL, uIon1=NULL, uIon2=NULL, uPlace=NULL;
int             iIon1=0, iIon2=0;
double          dCharge, dPertCharge, dICharge1, dICharge2;
double          dIonSize1, dIonSize2, dMinSize, dMinSeparation = 0.0;
int             i, iUnknown, ierr, iSystemCharge;
VECTOR          vNewPoint, vMaxPot, vMinPot;
HELP            hTemp;
LOOP            lAtoms;
ATOM            aAtom;
OCTREE          octTreeSolute = NULL;
VARARRAY        vaSolvent = NULL;
NeighborGrid    *ngSolventAtoms = NULL;
unsigned int    *iGroupStart = NULL;
VARARRAY        vaSolventPoints = NULL;
char            *sName1, *sName2=NULL, *sUnitName;

    VPTRACEENTER("addIons" );
    VPTRACEMULTIPLEEXIT("addIons" );
    BasicsResetInterrupt();

    /*
     *  Test args
     */
    ierr = 0;
    switch( iArgCount ) {
        case 3:  // One ion and desired number
          if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u n" )) ierr++;
          break;
        case 5:  // Two ions and desired number of each of them
          if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u n u n" ))
          { ierr++; break; }
          // Get the arguments for the second ion
          uIon2 = UNIT_from(oAssocObject( aaArgs[3] ));
          sName2 = sAssocName( aaArgs[3] );
          iIon2 = (int)dODouble( oAssocObject( aaArgs[4] ));
          break;
        default:
          ierr++;
          break;
    } /* end of switch */

    if ( uIon2  &&  iIon2 == 0 && iIon1 != 0)
    {
      VPFATAL("%s: '0' is not allowed as the value for the second ion.\n", sCmd );
      ++ierr;
    }
    if (dMinSeparation < 0.0)
    {
      VPFATAL("%s: %d is not a valid minimum distance between ions.\n",
            sCmd, dMinSeparation );
      ++ierr;
    }

    if ( ierr ) {
          hTemp = hHelp( sCmd );
          if ( hTemp == NULL ) {
                VPFATALDELAYEDEXIT("No help available on %s\n", sCmd );
          } else {
                VPFATALDELAYEDEXIT("\n%s\n", sHelpText(hTemp) );
          }
          return NULL;
    }

    /*
     *  Translate args common to both cases
     */
    uUnit = UNIT_from(oAssocObject( aaArgs[0] ));
    uIon1 = UNIT_from(oAssocObject( aaArgs[1] ));
    iIon1 = (int)dODouble( oAssocObject( aaArgs[2] ));

    sUnitName = sAssocName( aaArgs[0] );
    sName1 = sAssocName( aaArgs[1] );

    /*
     *  Consider target unit's charge
     */
    ContainerTotalCharge( CONTAINER_from( uUnit), &dCharge, &dPertCharge );
    iSystemCharge = (int)round( dCharge );
    if ( !iSystemCharge ) {
        VP0("%s nominal net charge is zero, with a residual charge of %g\n", sUnitName, dCharge);
        if ( iIon1 == 0 ) {
                VP0("%s: Can't neutralize.\n", sCmd );
                return NULL;
        }
        //VP0("Adding the ions anyway.\n");
    } else
        MESSAGE("dCharge:  %4.2lf\n", dCharge );

    /*
     *  Consider ion(s) charge
     */
    ContainerTotalCharge(CONTAINER_from(uIon1), &dICharge1, &dPertCharge );
    if ( !dICharge1 ) {
        VPFATALEXIT("%s: %s is not an ion and is not appropriate for placement.\n",
                sCmd, sName1);
        return NULL;
    }
    if ( uIon2 ) {
        ContainerTotalCharge(CONTAINER_from(uIon2), &dICharge2, &dPertCharge );
        if ( !dICharge2 ) {
           VPFATALEXIT("%s: %s is not an ion and is not appropriate for placement.\n",
                sCmd, sName2);
           return NULL;
        }
        if (dICharge1 * dICharge2 >0) {
           VPFATALEXIT("%s: Ion1 and Ion2 charge must have opposite sign.\n", sCmd);
           return NULL;
        }
    }

    /*
     *  Consider neutralization
     */
    if ( iIon1 == 0 ) {
        if ( dICharge1 * dCharge > 0) {
            if (!uIon2) {
                VPWARN( "%s: 1st Ion & target unit have charges of the same "
                        "sign:\n" "     unit charge = %g; ion1 charge = %g;\n"
                        "     can't neutralize.\n" , sCmd,
                              dCharge, dICharge1 );
                return(NULL);
            }
            VPWARN("%s: 1st Ion & target unit have charges of the same sign:\n"
                   "     unit charge = %g; ion1 charge = %g; ion2 charge = %g\n"
                   "     Swapping Ion1 & Ion2.\n" , sCmd,
                         dCharge, dICharge1, dICharge2 );
            uIon1 = uIon2;
            dICharge1 = dICharge2;
            char *s = sName1; sName1 = sName2; sName2 = s;
        }
        /*
         *  Get the nearest integer number of ions that
         *      we need to add to get as close as possible
         *      to neutral
         */
        iIon1 = (int)lrint( fabs(dCharge) / fabs(dICharge1) );
        if ( iIon1 == 0 )
            VP0(" %f %d %d %d\n", fabs( dCharge),
                (int)fabs( dCharge), (int)fabs( dICharge1 ),
                (int)(fabs( dCharge) / fabs( dICharge1 )) );
        if ( iIon2 ) {
                VP0("%s: Neutralization - can't specify 2nd ion count.\n", sCmd );
                return NULL;
        }
        VP0("%d %s ion%s required to neutralize.\n", iIon1, sName1,(iIon1 > 1 ? "s" : "") );
    }

    /*
     *  Consider ion sizes and positions.
     */
    dIonSize1 = 0.0;
    iUnknown = 0;
    lAtoms = lLoop( OBJEKT_from(uIon1), ATOMS );
    for(i=0; (aAtom = ATOM_from(oNext(&lAtoms))); i++) {
        if ( bAtomSetTmpRadius( aAtom ) )
                VP0("Using default radius %5.2f for ion %s\n",
                        ATOM_DEFAULT_RADIUS, sName1);
        dIonSize1 = MAX( dIonSize1, dAtomTemp( aAtom ) );
        if ( !bAtomFlagsSet( aAtom, ATOMPOSITIONKNOWN ) )
                iUnknown++;
    }
    if ( i > 1 ) {
        if ( iUnknown ) {
            VPFATALEXIT("Ion %s is polyatomic and has %d atoms w/ no position\n",
                sName1, iUnknown );
            return  NULL ;
        }
        VP0("Ion %s is polyatomic; multiplying max radius %5.2f by # atoms\n",
                sName1, dIonSize1 );
        dIonSize1 *= (double) i;
    } else if ( iUnknown ) {
        VECTOR          vPos;

        VectorDef( &vPos, 0.0, 0.0, 0.0 );
        lAtoms = lLoop( OBJEKT_from(uIon1), ATOMS );
        aAtom = ATOM_from(oNext(&lAtoms));
        AtomSetPosition( aAtom, vPos );
        AtomSetFlags( aAtom, ATOMPOSITIONKNOWN );
    }
    dMinSize = dIonSize1;
    dIonSize2 = 0.0;
    if ( uIon2 ) {
        iUnknown = 0;
        lAtoms = lLoop( OBJEKT_from(uIon2), ATOMS );
        for(i=0; (aAtom = ATOM_from(oNext(&lAtoms))); i++) {
                if ( bAtomSetTmpRadius( aAtom ) )
                        VP0("Using default radius %5.2f for ion %s\n",
                                ATOM_DEFAULT_RADIUS, sName2 );
                dIonSize2 = MAX( dIonSize2, dAtomTemp( aAtom ) );
                if ( !bAtomFlagsSet( aAtom, ATOMPOSITIONKNOWN ) )
                    iUnknown++;
        }
        if ( i > 1 ) {
            if ( iUnknown ) {
                VPFATALEXIT("Ion %s is polyatomic and has %d atoms w/ no position\n",
                    sName1, iUnknown );
                return  NULL ;
            }
            VP0("Ion %s is polyatomic; multiplying max radius %5.2f by # atoms",
                sName2, dIonSize2 );
            dIonSize2 *= (double) i;
        } else if ( iUnknown ) {
            VECTOR              vPos;

            VectorDef( &vPos, 0.0, 0.0, 0.0 );
            lAtoms = lLoop( OBJEKT_from(uIon1), ATOMS );
            aAtom = ATOM_from(oNext(&lAtoms));
            AtomSetPosition( aAtom, vPos );
            AtomSetFlags( aAtom, ATOMPOSITIONKNOWN );
        }
        dMinSize = MIN( dIonSize1, dIonSize2 );
    }

    VP0("Adding %d counter ions to \"%s\" using %gÅ grid, shell extent %gÅ\n",
                iIon1 + iIon2, sUnitName,
                GDefaults.dGridSpace,GDefaults.dShellExtent);

    int iTotalIons = iIon1 + iIon2;
    if ( !iTotalIons ) return NULL;

    if (!bIncludeSolvent) {// IncludeSolvent option groups solvent as part of solute
        vaSolvent = vaSolventResidues( uUnit );

        if ( vaSolvent ) {
            VP0("Solvent present: replacing closest with ion\n" );
            VP0("\t when steric overlaps occur ( < 3Å )\n" );
            int count = iVarArrayElementCount(vaSolvent);
            vaSolventPoints = vaVarArrayCreate(sizeof(Point));
            iGroupStart = MALLOC(sizeof(*iGroupStart)*(count+1));
            RESIDUE *rPRes = PVAI(vaSolvent,RESIDUE,0);
            for ( i=0; i<count; i++) {
                iGroupStart[i] = i;
                aAtom=ATOM_from(oContainerFirstObject(rPRes[i]));
                Point p= {
                    vAtomPosition(aAtom).dX,
                    vAtomPosition(aAtom).dY,
                    vAtomPosition(aAtom).dZ,
                    i, { .p = (void*)aAtom } };
                VarArrayAdd(vaSolventPoints, &p);
            }
            iGroupStart[count] = iVarArrayElementCount(vaSolventPoints);
            ngSolventAtoms = neighbor_grid_setup(
                   PVAI(vaSolventPoints, Point, 0),
                  iVarArrayElementCount(vaSolventPoints),count,iGroupStart,3.0);
        } else
            VP1(" (no solvent present)\n" );
    }
    TurnOffDisplayerUpdates();
    double dGridSize = GDefaults.dGridSpace;
    double dShellExtent = GDefaults.dShellExtent;
    if (iIon1 + iIon2 > 5) {
        const double dIonCount = (double) (iIon1 + iIon2);
        const double dDist = exp(log(dIonCount + 1.0)/3.0); // cubed root of total ion count + 1
        double dMaxSize = (dIonSize1 > dIonSize2 ? dIonSize1 : dIonSize2);
        double dMinShellExtent = dMaxSize * (dDist > 1.0 ? dDist : 1.0); // minimum size for saturated close packing
        if (dShellExtent < dMinShellExtent) {
            VPWARN("Inflating ion packing shell extent from %g to %g\n",
                    dShellExtent, dMinShellExtent);
            dShellExtent = dMinShellExtent;
        }
    }
    /*
     *  Build grid and calc potential on it.
     */
    START(tOct1);
    octTreeSolute = octOctTreeCreate( uUnit, OCT_SHELL, dGridSize, dMinSize, dShellExtent,
                 bIncludeSolvent );
    STOP(tOct1,"build OCTTREE");
    if ( !octTreeSolute ) {
        VP0("%s: No %s to add ions to\n", sCmd, bIncludeSolvent?"atoms":"solute" );
        return NULL;
    }
    double dExclusionRadius = dIonSize1 + dIonSize2;
    double dPointsPerPlacement = (4.0/3.0) * M_PI * pow(dExclusionRadius, 3) / pow(dGridSize, 3);
    if ( dPointsPerPlacement < 1.0 ) dPointsPerPlacement = 1.0;  /* floor at 1 - can't consume less than the landing point itself */
    double dSafetyFactor = 2.0;  /* arbitrary but conservative headroom - tune from real test data */
    double dRequiredPoints = (double)(iIon1 + iIon2) * dPointsPerPlacement;
    if ( (double)octTreeSolute->iTreePoints < dRequiredPoints * dSafetyFactor) {
        DFATAL("%s: requested %d ions may not fit - grid has %d candidate points, "
                 "estimated capacity needed ~%.0f (exclusion radius %.2f vs grid spacing %.2f). "
                 "Consider a finer grid, larger shell extent, or fewer ions.\n",
                 sCmd, iIon1+iIon2, octTreeSolute->iTreePoints, dRequiredPoints,
                 dExclusionRadius, dGridSize );
        return NULL;
    }

    START(tOct2);
    OctTreeInitCharges( octTreeSolute, GDefaults.iDielectricFlag,
                        dIonSize1, &vMinPot, &vMaxPot );
    STOP(tOct2,"OCTTREE init charges");
/*
OctTreePrintGrid( octTreeSolute, "Charge", COLOR_RANGE );
*/

    while ( iIon1 || iIon2 ) {
        if ( bBasicsInterrupt() ) goto CANCEL;
        if ( iIon1 ) {
                if ( dICharge1 < 0 )
                        vNewPoint = vMaxPot;
                else
                        vNewPoint = vMinPot;
                if ( vaSolvent )
                        CheckSolvent( uUnit, vaSolvent, uIon1, &vNewPoint, ngSolventAtoms );

                /*
                 *  Make a copy of ion unit and give it new point.
                 */
                uPlace = uCopyUnit(uIon1);
                ContainerCenterAt( CONTAINER_from( uPlace), vNewPoint );

                /*
                 *  Add ion to solute.
                 */
                UnitJoin( uUnit, uPlace );
                VP0("Placed %s in %s at (%7.2lf,%7.2lf,%7.2lf).\n",
                        sName1, sUnitName,
                        dVX(&(vNewPoint)),
                        dVY(&(vNewPoint)),
                        dVZ(&(vNewPoint)));
                /*
                 *  Delete ion from grid (allowing clearance to most likely
                 *      future adjacent ion) and update esp.
                 */
                OctTreeDeleteSphere( octTreeSolute, &vNewPoint,
                                dIonSize1 + (iIon2 ? dIonSize2 : dIonSize1) );
                OctTreeUpdateCharge( octTreeSolute, &vNewPoint,
                        (float)dICharge1, (iIon2 ? dIonSize2 : dIonSize1),
                        &vMaxPot, &vMinPot );
                iIon1--;
        }
        if ( iIon2 ) {
                if ( dICharge2 < 0 )
                        vNewPoint = vMaxPot;
                else
                        vNewPoint = vMinPot;
                if ( vaSolvent )
                        CheckSolvent( uUnit, vaSolvent, uIon2, &vNewPoint, ngSolventAtoms );

                /*
                 *  Make a copy of ion unit and give it new point.
                 */
                uPlace = uCopyUnit(uIon2);
                ContainerCenterAt( CONTAINER_from( uPlace), vNewPoint );

                /*
                 *  Add ion to solute.
                 */
                UnitJoin( uUnit, uPlace );
                VP0("Placed %s in %s at (%7.2lf,%7.2lf,%7.2lf).\n",
                        sName2, sUnitName,
                        dVX(&(vNewPoint)),
                        dVY(&(vNewPoint)),
                        dVZ(&(vNewPoint)));
                /*
                 *  Delete ion from grid (allowing clearance to most likely
                 *      future adjacent ion) and update esp.
                 */
                OctTreeDeleteSphere( octTreeSolute, &vNewPoint,
                                dIonSize2 + (iIon1 ? dIonSize1 : dIonSize2) );
                OctTreeUpdateCharge( octTreeSolute, &vNewPoint,
                        (float)dICharge2, (iIon1 ? dIonSize1 : dIonSize2),
                        &vMaxPot, &vMinPot );
                iIon2--;
        }
    }
/*
OctTreePrintGrid( octTreeSolute, "Charge2", COLOR_RANGE );
*/
    VP0("\nDone adding %d ions.\n", iTotalIons );
DONE:
    if (octTreeSolute)
        OctTreeDestroy( &octTreeSolute );
    if ( vaSolvent )
        VarArrayDestroy( &vaSolvent );
    if (vaSolventPoints)
        VarArrayDestroy( &vaSolventPoints );
    if (ngSolventAtoms)
        neighbor_grid_free(ngSolventAtoms);
    if (iGroupStart)
        FREE(iGroupStart);
    TurnOnDisplayerUpdates();
    ContainerDisplayerUpdate( CONTAINER_from( uUnit ));
    return NULL;

CANCEL:
    VP0("\n%s: Interrupted.\n", sCmd );
    BasicsResetInterrupt();
    goto DONE;
}

OBJEKT
oCmd_addIons( int iArgCount, ASSOC aaArgs[] )
{
    return addIons(iArgCount, aaArgs, false, "addIons" );
}

OBJEKT
oCmd_addIons2( int iArgCount, ASSOC aaArgs[] )
{
    return addIons(iArgCount, aaArgs, true, "addIons2" );
}

/*
 *     oCmd_addIonsRand
 *
 *     Robin Betz (2011)
 */
OBJEKT
oCmd_addIonsRand( int iArgCount, ASSOC aaArgs[] )
{
  UNIT            uUnit=NULL, uIon1=NULL, uIon2=NULL, uPlace=NULL;
  int             iIon1=0, iIon2=0;
  double          dCharge, dPertCharge, dICharge1, dICharge2;
  double          dMinSeparation = 0.0, dMolarity = 0.0;
  int             i, iUnknown, ierr, iRandom;
  HELP            hTemp;
  RESIDUE         *rPRes;
  ATOM            aAtom;
  LOOP            lAtoms;
  VARARRAY        vaSolvent = NULL;
  VECTOR          *vIonCenters;
  bool            bPlaceIon;
  int             iIonCount = 0;
  int             iFailCounter = 0;
  char            *sName1=NULL, *sName2=NULL, *sUnitName;
  char            *sCmd = "addIonsRand";

  // Setup
  BasicsResetInterrupt();

  // Test arguments
  ierr = 0;
  switch( iArgCount )
  {
    case 3:  // One ion and desired number / charge,
      if (bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u n" ))
          iIon1 = (int)dODouble( oAssocObject( aaArgs[2] ));
            // or two ions (auto neutralize, same as: uSolute uIon1 0 uIon2 0)
      else if (bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u u" ))
          uIon2 = UNIT_from(oAssocObject( aaArgs[2] ));
      else
          ++ierr;
      break;

    case 4:  // One ion and desired number / charge and minimum separation
      if (iObjectType(oAssocObject(aaArgs[2]))==ODOUBLEid) {
          if (bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u n n" )) {
              // Get the minimum separation
              dMinSeparation = dODouble( oAssocObject( aaArgs[3] ));
              break;
          }
      // Two ions and desired molarity
      } else if (bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u u n" )) {
          // Get the arguments for the second ion
          uIon2 = UNIT_from(oAssocObject( aaArgs[2] ));
          sName2 = sAssocName( aaArgs[2] );
          dMolarity = dODouble( oAssocObject( aaArgs[3] ));
          break;
      }
      ++ierr;
      break;

    case 5: // Two ions and number of each of them
      if (iObjectType(oAssocObject(aaArgs[2]))==ODOUBLEid) {
          if ( bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u n u n" )) {
              // Get the arguments for the second ion
              uIon2 = UNIT_from(oAssocObject( aaArgs[3] ));
              sName2 = sAssocName( aaArgs[3] );
              iIon2 = (int)dODouble( oAssocObject( aaArgs[4] ));
              break;
          }
      // Two ions and desired molarity and min sep
      } else if ( bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u u n n" )) {
          uIon2 = UNIT_from(oAssocObject( aaArgs[2] ));
          sName2 = sAssocName( aaArgs[2] );
          dMolarity = dODouble( oAssocObject( aaArgs[3] ));
          dMinSeparation = dODouble( oAssocObject( aaArgs[4] ));
          break;
      }
      ++ierr;
      break;

    case 6: // Two ions and number of each of them and minimum separation
      if (bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u n u n n" )) {
          // Get the arguments for the second ion
          uIon2 = UNIT_from(oAssocObject( aaArgs[3] ));
          sName2 = sAssocName( aaArgs[3] );
          iIon2 = (int)dODouble( oAssocObject( aaArgs[4] ));
          // Get the minimum separation
          dMinSeparation = dODouble( oAssocObject( aaArgs[5] ));
      } else ++ierr;
      break;

    default:
      ++ierr;
      break;
  }

  if (dMolarity < 0.0) {
    VPFATAL("%s: %g is not a valid molarity.\n",
          sCmd, dMolarity );
    ++ierr;
  }
  if (dMinSeparation < 0.0) {
    VPFATAL("%s: %g is not a valid minimum distance between ions.\n",
          sCmd, dMinSeparation );
    ++ierr;
  }

  // Display help if command is malformed
  if ( ierr ) {
    hTemp = hHelp( "addIonsRand" );
    if ( hTemp == NULL )
      VPFATALDELAYEDEXIT("No help available on addIonsRand\n" );
    else
      VPFATALDELAYEDEXIT("\n%s\n", sHelpText(hTemp) );
    return NULL;
  }

  // Translate unit, ion, and charge arguments
  uUnit = UNIT_from(oAssocObject( aaArgs[0] ));
  uIon1 = UNIT_from(oAssocObject( aaArgs[1] ));
  sUnitName = sAssocName( aaArgs[0] );
  sName1 = sAssocName( aaArgs[1] );

  // Check the unit's validity
  ContainerTotalCharge( CONTAINER_from( uUnit), &dCharge, &dPertCharge );
  int iSystemCharge = (int)round( dCharge );
  if ( !iSystemCharge && !dMolarity && !iIon1 ) {
    VP0("%s nominal net charge is zero, with a residual charge of %g\n", sUnitName, dCharge);
    VP0("%s: Can't neutralize.\n", sCmd );
    return NULL;
  }
  MESSAGE("dCharge:  %4.2lf\n", dCharge );

  // Make sure the ions are actually ions
  ContainerTotalCharge(CONTAINER_from(uIon1), &dICharge1, &dPertCharge );
  if ( !dICharge1 )
  {
    VPFATALEXIT("%s: %s is not an ion and is not appropriate for placement.\n",
          sCmd, sAssocName( aaArgs[1] ));
    return NULL;
  }
  if ( uIon2 )
  {
    ContainerTotalCharge(CONTAINER_from(uIon2), &dICharge2, &dPertCharge );
    if ( !dICharge2 )
    {
      VPFATALEXIT("%s: %s is not an ion and is not appropriate for placement.\n",
            sCmd, sAssocName( aaArgs[3] ));
      return NULL;
    }
    if (dICharge1 * dICharge2 >0) {
       VPFATALEXIT("%s: Ion1 and Ion2 charge must have opposite sign.\n", sCmd);
       return NULL;
    }
  }

  // Check validity of neutralization
  if ( !dMolarity && iIon1 == 0 && uIon2) {
    if ( dICharge1 * dCharge > 0) {
        if (!uIon2) {
            VPWARN( "%s: 1st Ion & target unit have charges of the same "
                     "sign:\n" "     unit charge = %g; ion1 charge = %g;\n"
                        "     can't neutralize.\n" , sCmd,
                     dCharge, dICharge1 );
            return(NULL);
        }
        VPWARN("%s: 1st Ion & target unit have charges of the same sign:\n"
               "     unit charge = %g; ion1 charge = %g; ion2 charge = %g\n"
               "     Swapping Ion1 & Ion2.\n" , sCmd,
                         dCharge, dICharge1, dICharge2);
        uIon1 = uIon2;
        dICharge1 = dICharge2;
        char *s = sName1; sName1 = sName2; sName2 = s;
    }
    /*
     *  Get the nearest integer number of ions that
     *      we need to add to get as close as possible
     *      to neutral
     */
    iIon1 = (int)lrint( fabs(dCharge) / fabs(dICharge1) );
    VP0("%d %s ion%s required to neutralize.\n", iIon1,
          sAssocName( aaArgs[1] ), (iIon1 > 1 ? "s" : "") );
  }

  // Check ion size and position
  iUnknown = 0;
  lAtoms = lLoop( OBJEKT_from(uIon1), ATOMS );
  for(i=0; (aAtom = ATOM_from(oNext(&lAtoms))); i++) {
    if ( bAtomSetTmpRadius( aAtom ) )
      VP0("Using default radius %5.2f for ion %s\n",
            ATOM_DEFAULT_RADIUS, sAssocName( aaArgs[1] ) );
    if ( !bAtomFlagsSet( aAtom, ATOMPOSITIONKNOWN ) )
      iUnknown++;
  }
  if ( i > 1 ) {
    if ( iUnknown ) {
      VPFATALEXIT("Ion %s is polyatomic and has %d atoms w/ no position\n",
            sAssocName( aaArgs[1] ), iUnknown );
      return  NULL ;
    }
  } else if ( iUnknown ) {
    VECTOR          vPos;

    VectorDef( &vPos, 0.0, 0.0, 0.0 );
    lAtoms = lLoop( OBJEKT_from(uIon1), ATOMS );
    aAtom = ATOM_from(oNext(&lAtoms));
    AtomSetPosition( aAtom, vPos );
    AtomSetFlags( aAtom, ATOMPOSITIONKNOWN );
  }
  if ( uIon2 ) {
    iUnknown = 0;
    lAtoms = lLoop( OBJEKT_from(uIon2), ATOMS );
    for(i=0; (aAtom = ATOM_from(oNext(&lAtoms))); i++) {
      if ( bAtomSetTmpRadius( aAtom ) )
        VP0("Using default radius %5.2f for ion %s\n",
              ATOM_DEFAULT_RADIUS, sAssocName( aaArgs[3] ) );
      if ( !bAtomFlagsSet( aAtom, ATOMPOSITIONKNOWN ) )
        iUnknown++;
    }
    if ( i > 1 ) {
      if ( iUnknown ) {
        VPFATALEXIT("Ion %s is polyatomic and has %d atoms w/ no position\n",
              sAssocName( aaArgs[1] ), iUnknown );
        return  NULL ;
      }
    } else if ( iUnknown ) {
      VECTOR              vPos;
      VectorDef( &vPos, 0.0, 0.0, 0.0 );
      lAtoms = lLoop( OBJEKT_from(uIon1), ATOMS );
      aAtom = ATOM_from(oNext(&lAtoms));
      AtomSetPosition( aAtom, vPos );
      AtomSetFlags( aAtom, ATOMPOSITIONKNOWN );
    }
  }

  vaSolvent = vaSolventResidues( uUnit );
  if ( !vaSolvent )
  {
    VPFATALEXIT("No solvent present. Add solvent first.\n");
    return NULL;
  }
  //  { 'B', "keep_input_solvent", "Keep_Input_Solvent", &GDefaults.bKeepInputSolvent, .defval.integer=0 },
  int iNumSolvent = iVarArrayElementCount(vaSolvent);
  if ( dMolarity ) {
    iIon1 = iIon2 = round((double)iNumSolvent * dMolarity / 55.51);
    VP0("Adding %d ion pairs to %d solvent residues to approximate %gM ionic strength\n",
        iIon1, iNumSolvent, dMolarity );
  }
  if ( iIon1 + iIon2 == 0 ) {
    VP0("No ions to add.\n");
    return NULL;
  }
  if ( iNumSolvent-iIon1-iIon2 <= 0)
  {
    VPFATALEXIT("Too few solvent molecules to add ions.\n" );
    return NULL;
  }
  VP0("Adding %d counter ions to \"%s\". %d solvent molecules will remain.\n",
        iIon1 + iIon2, sAssocName( aaArgs[0] ), iVarArrayElementCount(vaSolvent)-iIon1-iIon2);

  TurnOffDisplayerUpdates();
  if ( dMinSeparation )
    vIonCenters = (VECTOR*)MALLOC((iIon1+iIon2)*sizeof(VECTOR));

  // now actually add the ions
  while ( iIon1 || iIon2 )
  {
    if ( bBasicsInterrupt() ) goto CANCEL;
    if ( iIon1 ) {
      // Pick random solvent molecule to replace
      iRandom = iUniformRandom(iVarArrayElementCount( vaSolvent ));

      // Get position of solvent residue atom
      rPRes = (RESIDUE*)PVarArrayIndex ( vaSolvent, iRandom );
      aAtom = ATOM_from(oContainerFirstObject(*rPRes));
      VECTOR vNewPoint = vAtomPosition( aAtom );

      // Check that new point isn't too close to other ions
      bPlaceIon = true;
      if (dMinSeparation) {
        for (i=0; i<iIonCount; ++i) {
          if ( dVectorDistance(&vIonCenters[i], &vNewPoint) < dMinSeparation) {
            ++iFailCounter;
            bPlaceIon = false;
            break;
          }
        }
      } // end min sep
      if ( bPlaceIon ) {
        if (dMinSeparation) VP2("%d: ", iIonCount);
        VP0("Placed %s in %s at (%7.2lf,%7.2lf,%7.2lf).\n",
              sName1, sUnitName,
              dVX(&(vNewPoint)),
              dVY(&(vNewPoint)),
              dVZ(&(vNewPoint)));

        // Save this ion's position if desired
        if ( dMinSeparation ) vIonCenters[iIonCount++] = vNewPoint;

        // Copy ion unit, position, and add it to the unit
        uPlace = uCopyUnit(uIon1);
        ContainerCenterAt(CONTAINER_from( uPlace), vNewPoint );
        UnitJoin( uUnit, uPlace );

        // Delete the solvent residue that was replaced
        REF( *rPRes );  /* bContainerRemove() needs this */
        ResidueYouAreBeingRemoved( *rPRes );
        if ( bContainerRemove( CONTAINER_from(uUnit), OBJEKT_from(*rPRes )) == false)
          DFATAL("rmv solv %d failed\n", iRandom );
        ContainerDestroy((CONTAINER *) rPRes );
        rPRes = NULL;
        VarArrayDelete(vaSolvent, iRandom);

        --iIon1;
      }
    } // end ion1
    if ( iIon2 ) {
      // Pick random solvent molecule to replace
      iRandom = iUniformRandom(iVarArrayElementCount( vaSolvent ));

      // Get position of solvent residue atom
      rPRes = (RESIDUE*)PVarArrayIndex ( vaSolvent, iRandom );
      aAtom = ATOM_from(oContainerFirstObject(*rPRes));
      VECTOR vNewPoint = vAtomPosition( aAtom );

      // Check that new point isn't too close to other ions
      bPlaceIon = true;
      if (dMinSeparation) {
        for (i=0; i<iIonCount; ++i) {
          if ( dVectorDistance(&vIonCenters[i], &vNewPoint) < dMinSeparation) {
            ++iFailCounter;
            bPlaceIon = false;
            break;
          }
        }
      }
      if ( bPlaceIon ) {
        if (dMinSeparation) VP2("%d: ", iIonCount);
        VP0("Placed %s in %s at (%7.2lf,%7.2lf,%7.2lf).\n",
              sName2, sUnitName,
              dVX(&(vNewPoint)),
              dVY(&(vNewPoint)),
              dVZ(&(vNewPoint)));

        // Save this ion's position if desired
        if ( dMinSeparation ) vIonCenters[iIonCount++] = vNewPoint;

        // Copy ion unit, position, and add it to the unit
        uPlace = uCopyUnit(uIon2);
        ContainerCenterAt(CONTAINER_from( uPlace), vNewPoint );
        UnitJoin( uUnit, uPlace );

        // Delete the solvent residue that was replaced
        REF( *rPRes );  /* bContainerRemove() needs this */
        ResidueYouAreBeingRemoved( *rPRes );
        if ( bContainerRemove( CONTAINER_from(uUnit), OBJEKT_from(*rPRes )) == false)
          DFATAL("rmv solv %d failed\n", iRandom );
        ContainerDestroy((CONTAINER *) rPRes );
        rPRes = NULL;
        VarArrayDelete(vaSolvent, iRandom);
        --iIon2;
      }
    }
    if ( iFailCounter > 100 ) {
      DFATAL("Impossible to place %d ions with minimum separation of %f.\n",
               iIon1 + iIon2, dMinSeparation );
      break;
    }
  }
DONE:
  VarArrayDestroy(&vaSolvent);
  return NULL;
CANCEL:
  VP0("\n%s: Interrupted.\n", sCmd );
  BasicsResetInterrupt();
  goto DONE;
}

OBJEKT
oCmd_addIonSolv( int iArgCount, ASSOC aaArgs[] )
{
UNIT            uUnit=NULL, uIon1=NULL, uIon2=NULL;
int             iIon1=0, iIon2=0;
double          dCharge, dPertCharge, dICharge1, dICharge2;
double          dIonSize1, dIonSize2;
int             i, iUnknown, ierr, iMinPotRes, iMaxPotRes, iReplace;
HELP            hTemp;
LOOP            lAtoms;
ATOM            aAtom;
VARARRAY        vaSolvent;
char            *sCmd = "addIonSolv";

    BasicsResetInterrupt();

    /*
     *  Test args
     */
    ierr = 0;
    switch( iArgCount ) {
        case 3:
          if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u n" ))
                ierr++;
          break;
        case 5:
          if ( !bCmdGoodArguments( sCmd, iArgCount, aaArgs, "u u n u n" )) {
                ierr++;
                break;
          }
          /*
           *  Translate the 2 extra args
           */
          uIon2 = UNIT_from(oAssocObject( aaArgs[3] ));
          iIon2 = (int)dODouble( oAssocObject( aaArgs[4] ));
          if ( uIon2  &&  iIon2 == 0 ) {
              VPFATAL("%s: '0' is not allowed as the value for the second ion.\n",
                                                sCmd );
              ierr++;
          }
          break;
        default:
          ierr++;
          break;
    } /* end of switch */

    if ( ierr ) {
          hTemp = hHelp( "addionsolv" );
          if ( hTemp == NULL ) {
                VPFATALDELAYEDEXIT("No help available on addIonSolv\n" );
          } else {
                VPFATALDELAYEDEXIT("\n%s\n", sHelpText(hTemp) );
          }
          return NULL;
    }

    /*
     *  Translate args common to both cases
     */
    uUnit = UNIT_from(oAssocObject( aaArgs[0] ));
    uIon1 = UNIT_from(oAssocObject( aaArgs[1] ));
    iIon1 = (int)dODouble( oAssocObject( aaArgs[2] ));

    /*
     *  make sure unit has solvent
     */
    vaSolvent = vaSolventResidues( uUnit );
    if ( vaSolvent == NULL ) {
        VPFATALEXIT("No solvent present: solvate 1st or use addIons\n" );
        return NULL;
    }

    /*
     *  Consider target unit's charge
     */
    ContainerTotalCharge( CONTAINER_from( uUnit), &dCharge, &dPertCharge );
    if ( !dCharge ) {
        VP0("%s has a charge of 0.\n", sAssocName( aaArgs[1] ));
        if ( iIon1 == 0 ) {
                VP0("%s: Can't neutralize.\n", sCmd );
                return NULL;
        }
        VP0("Adding the ions anyway.\n");
    } else
        MESSAGE("dCharge:  %4.2lf\n", dCharge );

    /*
     *  Consider ion(s) charge
     */
    ContainerTotalCharge(CONTAINER_from(uIon1), &dICharge1, &dPertCharge );
    if ( !dICharge1 ) {
        VPFATALEXIT("%s: %s is not an ion and is not appropriate for placement.\n",
                sCmd, sAssocName( aaArgs[1] ));
        return NULL;
    }
    if ( uIon2 ) {
        ContainerTotalCharge(CONTAINER_from(uIon2), &dICharge2, &dPertCharge );
        if ( !dICharge2 ) {
           VPFATALEXIT("%s: %s is not an ion and is not appropriate for placement.\n",
                sCmd, sAssocName( aaArgs[3] ));
           return NULL;
        }
    }

    /*
     *  Consider neutralization
     */
    if ( iIon1 == 0 ) {
        if ( (dICharge1 < 0  &&  dCharge < 0) ||
             (dICharge1 > 0  &&  dCharge > 0)) {
                VPWARN("%s: 1st Ion & target unit have charges of the same "
                         "sign:\n" "     unit charge = %g; ion1 charge = %g;\n"
                            "     can't neutralize.\n" , sCmd,
                         dCharge, dICharge1 );
                return NULL;
        }
        /*
         *  Get the nearest integer number of ions that
         *      we need to add to get as close as possible
         *      to neutral
         */
        iIon1 = (int)lrint( fabs(dCharge) / fabs(dICharge1) );
        if ( iIon1 == 0 )
            VP0(" %f %d %d %d\n", fabs( dCharge),
                (int)fabs( dCharge), (int)fabs( dICharge1 ),
                (int)(fabs( dCharge) / fabs( dICharge1 )) );
        if ( uIon2 ) {
                VP0("%s: Neutralization - can't do 2nd ion.\n", sCmd );
                return NULL;
        }
        VP0("%d %s ion%s required to neutralize.\n", iIon1,
                sAssocName( aaArgs[1] ), (iIon1 > 1 ? "s" : "") );
    }

    /*
     *  Consider ion sizes and positions.
     */
    dIonSize1 = 0.0;
    iUnknown = 0;
    lAtoms = lLoop( OBJEKT_from(uIon1), ATOMS );
    for(i=0; (aAtom = ATOM_from(oNext(&lAtoms))); i++) {
        if ( bAtomSetTmpRadius( aAtom ) )
                VP0("Using default radius %5.2f for ion %s\n",
                        ATOM_DEFAULT_RADIUS, sAssocName( aaArgs[1] ) );
        dIonSize1 = MAX( dIonSize1, dAtomTemp( aAtom ) );
        if ( !bAtomFlagsSet( aAtom, ATOMPOSITIONKNOWN ) )
                iUnknown++;
    }
    if ( i > 1 ) {
        if ( iUnknown ) {
            VPFATALEXIT("Ion %s is polyatomic and has %d atoms w/ no position\n",
                sAssocName( aaArgs[1] ), iUnknown );
            return  NULL ;
        }
        VP0("Ion %s is polyatomic; multiplying max radius %5.2f by # atoms\n",
                sAssocName( aaArgs[1] ), dIonSize1 );
        dIonSize1 *= (double) i;
    } else if ( iUnknown ) {
        VECTOR          vPos;

        VectorDef( &vPos, 0.0, 0.0, 0.0 );
        lAtoms = lLoop( OBJEKT_from(uIon1), ATOMS );
        aAtom = ATOM_from(oNext(&lAtoms));
        AtomSetPosition( aAtom, vPos );
        AtomSetFlags( aAtom, ATOMPOSITIONKNOWN );
    }
    if ( uIon2 ) {
        dIonSize2 = 0.0;
        iUnknown = 0;
        lAtoms = lLoop( OBJEKT_from(uIon2), ATOMS );
        for(i=0; (aAtom = ATOM_from(oNext(&lAtoms))); i++) {
                if ( bAtomSetTmpRadius( aAtom ) )
                        VP0("Using default radius %5.2f for ion %s\n",
                                ATOM_DEFAULT_RADIUS, sAssocName( aaArgs[3] ) );
                dIonSize2 = MAX( dIonSize2, dAtomTemp( aAtom ) );
                if ( !bAtomFlagsSet( aAtom, ATOMPOSITIONKNOWN ) )
                    iUnknown++;
        }
        if ( i > 1 ) {
            if ( iUnknown ) {
                VPFATALEXIT("Ion %s is polyatomic and has %d atoms w/ no position\n",
                    sAssocName( aaArgs[1] ), iUnknown );
                return  NULL ;
            }
            VP0("Ion %s is polyatomic; multiplying max radius %5.2f by # atoms",
                sAssocName( aaArgs[3] ), dIonSize2 );
            dIonSize2 *= (double) i;
        } else if ( iUnknown ) {
            VECTOR              vPos;

            VectorDef( &vPos, 0.0, 0.0, 0.0 );
            lAtoms = lLoop( OBJEKT_from(uIon1), ATOMS );
            aAtom = ATOM_from(oNext(&lAtoms));
            AtomSetPosition( aAtom, vPos );
            AtomSetFlags( aAtom, ATOMPOSITIONKNOWN );
        }
    }

    VP0("Adding %d counter ions to \"%s\", substituting solvent\n",
                        iIon1 + iIon2, sAssocName( aaArgs[0] ));

    if ( iIon1 + iIon2 == 0 ) {
        VarArrayDestroy( &vaSolvent );
        return NULL;
    }

    if ( iIon1 + iIon2 > iVarArrayElementCount( vaSolvent ) ) {
        VarArrayDestroy( &vaSolvent );
        VPFATALEXIT("Can't do it - more ions than solvent\n" );
        return NULL;
    }

    /*
     *  calc potential on solvent centers.
     */
    VP0("calculating initial potential at 1st atom in each solvent res..\n");
    ToolInitSolventPotential( uUnit, vaSolvent, &iMinPotRes, &iMaxPotRes );

    VP0("placing ions..\n" );
    TurnOffDisplayerUpdates();
    while ( iIon1 || iIon2 ) {
        if ( bBasicsInterrupt() ) goto CANCEL;
        if ( iIon1 ) {
                if ( dICharge1 < 0 )
                        iReplace = iMaxPotRes;
                else
                        iReplace = iMinPotRes;
                ToolReplaceSolvent( uUnit, vaSolvent,
                                        iReplace, uIon1, dICharge1,
                                        &iMinPotRes, &iMaxPotRes );
                VP0("Placed %s in %s.\n",
                        sAssocName( aaArgs[1] ), sAssocName( aaArgs[0] ) );

                iIon1--;
        }
        if ( iIon2 ) {
                if ( dICharge2 < 0 )
                        iReplace = iMaxPotRes;
                else
                        iReplace = iMinPotRes;

                ToolReplaceSolvent( uUnit, vaSolvent,
                                        iReplace, uIon2, dICharge2,
                                        &iMinPotRes, &iMaxPotRes );
                VP0("Placed %s in %s.\n",
                        sAssocName( aaArgs[3] ), sAssocName( aaArgs[0] ) );
                iIon2--;
        }
    }
    VP0("\nDone adding ions.\n" );
    goto DONE;
CANCEL:
    VP0("\n%s: Interrupted.\n", sCmd );
DONE:
    VarArrayDestroy( &vaSolvent );
    TurnOnDisplayerUpdates();
    ContainerDisplayerUpdate( CONTAINER_from( uUnit ));
    return NULL;
}

/*
 *      oCmd_addIonsNear
 *
 */
OBJEKT
oCmd_addIonsNear( int iArgCount, ASSOC aaArgs[] )
{
        VP0("Not implemented\n");
        return NULL;
}

/*
 *      oCmd_alias
 *
 *      Author:  David A. Rivkin (1992)
 *
 *      Add or remove an entry to or list entries in the Alias table.
 *      If both strings are present, then add the alias to the table.
 *      If only the first one string is there then remove the alias.
 *      If no arguments are given, then list the current aliases.
 *
 *      alias [alias[, string]]
 *
 *      Arguments:
 *      optional[0] - OSTRING alias is the alias to use
 *      optional[1] - OSTRING string is the original command word
 *
 */

VARARRAY GvaAlias;


OBJEKT
oCmd_alias( int iArgCount, ASSOC aaArgs[] )
{
STRING  sCommand;
STRING  sAlias;
ALIASt  aAlias, *PaAlias;
bool    bOK;
int     iAliases, i;
HELP    hTemp;
char    *sCmd = "alias";

    switch ( iArgCount) {
        case 2 : /* Adds an alias */
                MESSAGE("Add alias:  %p, %s\n",
                        oAssocObject(aaArgs[0]),
                        sObjectType( oAssocObject(aaArgs[1])));
                if ( bCmdGoodArguments( sCmd, iArgCount, aaArgs, "s s" )) {
                    strcpy( sAlias, sOString( oAssocObject( aaArgs[0] )));
                    strcpy( sCommand, sOString( oAssocObject( aaArgs[1] )));
                } else {
                    hTemp = hHelp( "alias" );
                    if ( hTemp == NULL ) {
                        VPFATALDELAYEDEXIT("No help available on \"alias\".\n" );
                    } else {
                        VPFATALDELAYEDEXIT("\n%s\n", sHelpText(hTemp) );
                    }
                    return NULL;
                }
                StringLower( sAlias );
                StringLower( sCommand );
                bOK = false;

                /* Make sure that sCommand is a command */
                for ( i=0; cCommands[i].fCallback; i++ ) {
                    if ( strcasecmp( sCommand, cCommands[i].sName ) == 0 ) {
                        bOK = true;
                        break;
                    }
                }
                if ( !bOK ) {
                    VPFATALEXIT("%s: '%s' is not a command.\n"
                        "Please check the spelling and try again.\n", sCmd, sCommand );
                    return  NULL ;
                }

                /* Make sure that the alias is not an existing command */
                for ( i=0; cCommands[i].sName[0] != 0; i++ ) {
                    if ( strcasecmp( sAlias, cCommands[i].sName ) == 0 ) {
                        bOK = false;
                        break;
                    }
                }
                if ( !bOK ) {
                    VPFATALEXIT( "%s: '%s' is already one of the commands.\n"
                        "Please try something different.\n", sCmd, sAlias );
                    return( NULL );
                }

                if ( GvaAlias == 0 ) {
                    MESSAGE("Creating a new Global Alias Structure.\n" );
                    GvaAlias = vaVarArrayCreate( sizeof( ALIASt ));
                }
                /*
                 *  Make sure that the alias does not already exist,
                 *      if so, replace it
                 */
                bOK = false;
                iAliases = iVarArrayElementCount( GvaAlias );
                if ( iAliases ) {
                    PaAlias = PVAI( GvaAlias, ALIASt, 0 );
                    for ( i = 0; i < iAliases; i++, PaAlias++ ) {
                        if ( strcmp( sAlias, PaAlias->sName ) == 0 ) {
                            bOK = true;
                            strcpy( PaAlias->sCommand, sCommand );
                            break;
                        }
                    }
                }
                if ( !bOK ) {
                    memset(&aAlias, 0, sizeof(aAlias)); /* for Purify */
                    strcpy( aAlias.sName, sAlias );
                    strcpy( aAlias.sCommand, sCommand );
                    VarArrayAdd( GvaAlias, (GENP)&aAlias );
                }
                break;

        case 1 : /* Remove an alias from the list */
                if ( bCmdGoodArguments( sCmd, iArgCount, aaArgs, "s" )) {
                    strcpy( sAlias, sOString( oAssocObject( aaArgs[0] )));
                } else if ( bCmdGoodArguments( sCmd, iArgCount, aaArgs, "z" )) {
                    strcpy( sAlias, sAssocName( aaArgs[0] ));
                } else {
                    hTemp = hHelp( "alias" );
                    if ( hTemp == NULL ) {
                        VPFATALDELAYEDEXIT("No help available on \"alias\".\n" );
                    } else {
                        VPFATALDELAYEDEXIT("\n%s\n", sHelpText(hTemp) );
                    }
                    return NULL;
                }
                StringLower( sAlias );
                iAliases = iVarArrayElementCount( GvaAlias );
                if ( !iAliases ) {
                    VPWARN("%s: There are no aliases loaded.\n", sCmd );
                    return NULL;
                }
                PaAlias = PVAI( GvaAlias, ALIASt, 0 );
                for ( i = 0; i < iAliases; i++, PaAlias++ ) {
                    if ( strcmp( sAlias, PaAlias->sName ) == 0 ) {
                        VarArrayDelete( GvaAlias, i );
                        break;
                    }
                }
                break;
        case 0 : /*  List all the aliases */
                iAliases = iVarArrayElementCount( GvaAlias );
                if ( !iAliases ) {
                    VPWARN("There are no aliases loaded.\n" );
                    return NULL;
                }
                VP0("Current Aliases  [alias....command]\n" );
                for ( i = 0; i < iAliases; i++ ) {
                    aAlias = *PVAI( GvaAlias, ALIASt, i );
                    if (( i % 2 ) == 0 ) { /* An odd entry */
                        VP0("%-10s..%-24s",
                                aAlias.sName, aAlias.sCommand );
                    } else {
                        VP0("%-10s..%-10s\n",
                                aAlias.sName, aAlias.sCommand );
                    }
                }
                if (( i % 2 ) == 1 ) {
                    /* Left over entry (odd number of entries) */
                    VP0("\n" );
                }
                break;

        default: hTemp = hHelp( "alias" );
                if ( hTemp == NULL ) {
                    VPFATALEXIT("No help available on alias.\n" );
                } else {
                    VPFATALEXIT("\n%s\n", sHelpText(hTemp) );
                }
                break;

   } /* end of switch */
   return NULL;
}


#if 0
/*
 *      oCmd_update
 *
 *      Author: David A. Rivkin (1993)
 *
 *      Scans a UNIT for any residue that is different from
 *      the residues in the UNIT and changes them to that residue.
 *      OR
 *      Checks a particular residue in a unit to see if the template
 *      residues have changed and modifies it as needed.
 *
 *      update unit/residue
 *              ARG[0] - UNIT/RESIDUE
 *
 */
OBJEKT
oCmd_update( int iArgCount, ASSOC aaArgs[] )
{
UNIT            uUnit;
RESIDUE         rNew, rOld, rCopy, rTemp;
int             iResidues;
LOOP            lRes, lResidue;
int             i, iNext;
HELP            hTemp;
STRING          sResName, sTemp;
DICTLOOP        dlLoop;
char            *sCmd = "update";

    if ( iArgCount != 1 ) {
        VPFATAL("Invalid number of arguments.\n" );
        hTemp = hHelp( "update" );
        if ( hTemp == NULL ) {
            VPFATALDELAYEDEXIT("No help available on update.\n" );
        } else {
            VPFATALDELAYEDEXIT("%s\n", sHelpText(hTemp) );
        }
        return NULL;
    }
    if (! bCmdGoodArguments( sCmd, iArgCount, aaArgs, "ru" )) {
        hTemp = hHelp( "update" );
        if ( hTemp == NULL ) {
            VPFATALDELAYEDEXIT("No help available on update.\n" );
        } else {
            VPFATALDELAYEDEXIT("%s\n", sHelpText(hTemp) );
        }
        return NULL;
    }
    switch ( iObjectType(oAssocObject( aaArgs[0] ))) {
        case    UNITid:
            DisplayerAccumulateUpdates();
            uUnit = UNIT_from(oAssocObject( aaArgs[0] ));
            lResidue = lLoop( OBJEKT_from(uUnit), RESIDUES );
            while ( (rOld = RESIDUE_from(oNext( &lResidue ))) ) {
                strcpy( sTemp, sContainerName(CONTAINER_from(rOld)));
                for ( i = 0; i < strlen(sTemp); i++){
                    if ( sTemp[i] == ' ') break;
                    sResName[i] = sTemp[i];
                }
                sResName[i] = '\0';
                // FIXME, bad code, do this:
                //  rNew =  RESIDUE_from(yPDictionaryFind( GdVariables, sResName))
                dlLoop = ydlDictionaryLoop(GdVariables);
                while ( yPDictionaryNext( GdVariables, &dlLoop )) {
                    if (!strcmp(sResName, sDictLoopKey(dlLoop))) {
                        iResidues = 0;
                        lRes = lLoop( OBJEKT_from(PDictLoopData)( dlLoop ), RESIDUES);
                        rNew = rCopy = NULL;
                        while ((rNew = RESIDUE_from(oNext( &lRes )))) {
                            iResidues++;
                                /* Make sure that there are not more than
                                        1 residue in the unit */
                            if ( iResidues > 1 ) {
                                MESSAGE("Unit name matches but has more than 1 residue.\n"
                                );
                                goto NEXTRES1;
                            }
                            rCopy = (RESIDUE)oCopy(OBJEKT_from(rNew ));
                        }
                        ResidueMutate( rCopy, rOld );

                        /* Remove the old RESIDUE from the UNIT */
                        /* And add the new RESIDUE, without incrementing */
                        /* the UNITs next child sequence number */

                        REF( rOld );  /* bContainerRemove() needs this */
                        bContainerRemove( CONTAINER_from( uUnit), OBJEKT_from(rOld ));
                        iNext = iContainerNextChildsSequence( CONTAINER_from( uUnit ));
                        ContainerAdd( CONTAINER_from(uUnit), OBJEKT_from(rCopy ));
                        ContainerSetSequence( CONTAINER_from( rCopy), iContainerSequence(CONTAINER_from( rOld)) );
                        ContainerSetNextChildsSequence( CONTAINER_from( uUnit), iNext );

                        DEREF( rOld );

                        VP0("Updating residue %s.\n", sResName );
                        goto NEXTSEQ;
                    }
NEXTRES1:       ;
                }
NEXTSEQ:        ;
            }
            DisplayerReleaseUpdates();
            break;

        case    RESIDUEid:
            DisplayerAccumulateUpdates();
            rOld = RESIDUE_from(oAssocObject( aaArgs[0] ));
            uUnit = (UNIT)cContainerWithin(CONTAINER_from(rOld ));
            strcpy( sTemp, sContainerName(CONTAINER_from(rOld )));
            for ( i = 0; i < strlen(sTemp); i++){
                if ( sTemp[i] == ' ') break;
                    sResName[i] = sTemp[i];
            }
            sResName[i] = '\0';

            dlLoop = ydlDictionaryLoop(GdVariables);
            // FIXME, bad code, do this:
            //  rNew =  RESIDUE_from(yPDictionaryFind( GdVariables, sResName))
            while ( yPDictionaryNext( GdVariables, &dlLoop )) {
                if (!strcmp(sResName, sDictLoopKey(dlLoop))) {
                    iResidues = 0;
                    lRes = lLoop( OBJEKT_from(PDictLoopData)( dlLoop ), RESIDUES);
                    rNew = NULL;
                    while ( (rTemp = RESIDUE_from(oNext( &lRes))) ) {
                        rNew = rTemp;
                        iResidues++;
                    }
                        /* Make sure that there are not more than
                                1 residue in the unit */
                    if ( iResidues > 1 ) {
                        MESSAGE("Unit name matches but has more than 1 residue.\n" );
                        goto NEXTRES2;
                    }

                    rCopy = rCopyResidue(rNew);
                    ResidueMutate( rCopy, rOld );

                        /* Remove the old RESIDUE from the UNIT */
                        /* And add the new RESIDUE, without incrementing */
                        /* the UNITs next child sequence number */

                    REF( rOld );  /* bContainerRemove() needs this */
                    bContainerRemove( CONTAINER_from( uUnit), OBJEKT_from(rOld ));
                    iNext = iContainerNextChildsSequence( CONTAINER_from( uUnit ));
                    ContainerAdd( CONTAINER_from(uUnit), OBJEKT_from(rCopy ));
                    ContainerSetSequence( CONTAINER_from( rCopy), iContainerSequence(CONTAINER_from( rOld)) );
                    ContainerSetNextChildsSequence( CONTAINER_from( uUnit), iNext );

                    DEREF( rOld );

                    VP0("Updating residue %s.\n", sResName );
                    goto NEXTRES2;
                }
NEXTRES2:       ;
            }
            DisplayerReleaseUpdates();
            break;
        default:
            DFATAL("Impossible control flow.\n" );
            break;
    }
    return NULL;
}
#endif

/*
 *      oCmd_flip
 *      Based on XAUESelectedAtomsFlipChirality in xaUnitEditor.c (xleap source code)
 *      Author: Christine cezard (2007)
 *      Universite de Picardie - Jules Verne, Amiens
 *      http://q4md-forcefieldtools.org
 *
 */
OBJEKT
oCmd_flip(int iArgCount, ASSOC aaArgs[])
{
UNIT        uUnit;
LOOP        lAtoms;
ATOM        aAtom;

    uUnit = UNIT_from(oAssocObject(aaArgs[0]));
    if ( uUnit == NULL ) return NULL;

    lAtoms = lLoop(OBJEKT_from( uUnit), ATOMS );
    while ( (aAtom = ATOM_from(oNext(&lAtoms)))) {
        if ( bAtomFlagsSet( aAtom, ATOMSELECTED ) ) {
            bBuildFlipChiralityFor(CONTAINER_from( uUnit), aAtom );
        }
    }


    return NULL;
}

/*
 *      oCmd_relax
 *      Based on Based on XAUERelaxSelectionInFramework in xaUnitEditor.c (xleap source code)
 *      Author: Christine Cezard (2007)
 *      Universite de Picardie - Jules Verne, Amiens
 *      http://q4md-forcefieldtools.org
 *
 */
OBJEKT
oCmd_relax(int iArgCount, ASSOC aaArgs[])
{
MINIMIZER       mStrain;
UNIT            uUnit;

    uUnit = UNIT_from(oAssocObject(aaArgs[0]));
    if ( uUnit == NULL ) return NULL;

        /* Setup a MINIMIZER and give it a callback to use */
        /* to update the display every step of the minimization */

    mStrain = mMinimizerCreate();
        /* Set up the MINIMIZER to use, and turn off any */
        /* control-c that may have been hit before */

    BasicsResetInterrupt();
    SelectRelaxInFramework( uUnit, mStrain );

    return NULL;
}

/*
 *      oCmd_addH
 *      Based on Based on XAUEAddHydrogensBuildExternals in xaUnitEditor.c (xleap source code)
 *      Author: D. Roe (2011)
 *      Rutgers University
 *
 */
OBJEKT
oCmd_addH(int iArgCount, ASSOC aaArgs[])
{
LOOP            lAtom, lSpan;
ATOM            aAtom, aStart;
UNIT            uUnit;
int             iDum;

    DisplayerAccumulateUpdates();
    uUnit = UNIT_from(oAssocObject(aaArgs[0]));
    if ( uUnit == NULL ) return NULL;

        /* Add hydrogens */

    ModelAddHydrogens( uUnit );

        /* Try to build geometries for simple rings */

    BuildInternalsForSimpleRings( CONTAINER_from(uUnit ));

        /* Assign internal coordinates for all internals that */
        /* include atoms that need building */

    lAtom = lLoop(OBJEKT_from( uUnit), ATOMS );
    BuildInternalsUsingFlags( &lAtom, ATOMPOSITIONDRAWN, 0,
                                ATOMNEEDSBUILD,
                                ATOMPOSITIONDRAWN );

        /* Build spanning trees for all atoms that need building */
        /* and build external coordinates for those atoms */

    lAtom = lLoop(OBJEKT_from( uUnit), ATOMS );
    while ( (aAtom = ATOM_from(oNext(&lAtom)))) {
        if ( bAtomFlagsSet( aAtom, ATOMNEEDSBUILD ) ) {

                        /* Look for a collision with an ATOM that has */
                        /* already been built */
                        /* Then start building from there */

            lSpan = lLoop(OBJEKT_from( aAtom), SPANNINGTREE );
            LoopDefineVisibleAtoms( &lSpan, ATOMNEEDSBUILD );
            while ( oNext(&lSpan) );
            if ( iLoopInvisibleCollisionCount(&lSpan) > 0 ) {
                aStart = aLoopLastCollisionAtom(&lSpan);
                lSpan = lLoop(OBJEKT_from( aStart), SPANNINGTREE );
            } else {
                lSpan = lLoop(OBJEKT_from( aAtom), SPANNINGTREE );
            }
            LoopDefineVisibleAtoms( &lSpan, ATOMNEEDSBUILD );
            iDum = 0;   /* for purify */
            BuildExternalsUsingFlags( &lSpan, ATOMNEEDSBUILD, 0,
                                        ATOMPOSITIONKNOWN|ATOMPOSITIONBUILT,
                                        ATOMNEEDSBUILD,
                                        &iDum, &iDum, &iDum, true );
        }
    }

                /* Destroy all of the INTERNALs */

    lAtom = lLoop(OBJEKT_from( uUnit), ATOMS );
    BuildDestroyInternals( &lAtom );

    DisplayerReleaseUpdates();

    return NULL;
}


COMMANDt        cCommands[] = {

        { "add",                oCmd_add },
        { "addH",               oCmd_addH },
        { "addIons",            oCmd_addIons },
        { "addIons2",           oCmd_addIons2 },
        { "addIonSolv",         oCmd_addIonSolv },
        { "addIonsRand",        oCmd_addIonsRand },
        { "addPath",            oCmd_addPath },
        { "addPdbAtomMap",      oCmd_addPdbAtomMap },
        { "addPdbResMap",       oCmd_addPdbResMap },
        { "addAtomTypes",       oCmd_addAtomTypes },
        { "alias",              oCmd_alias },
        { "alignAxes",          oCmd_alignAxes },
        { "bond",               oCmd_bond },
        { "bondByDistance",     oCmd_bondByDistance },
        { "center",             oCmd_center },
        { "charge",             oCmd_charge },
        { "check",              oCmd_check },
        { "clearPdbAtomMap",    oCmd_clearPdbAtomMap },
        { "clearPdbResMap",     oCmd_clearPdbResMap },
        { "clearVariables",     oCmd_clearVariables },
        { "combine",            oCmd_combine },
        { "copy",               oCmd_copy },
        { "createAtom",         oCmd_createAtom },
        { "createParmset",      oCmd_createParmset },
        { "createResidue",      oCmd_createResidue },
        { "createUnit",         oCmd_createUnit },
        { "crossLink",          oCmd_crossLink },
        { "debugOff",           oCmd_debugOff },
        { "debugOn",            oCmd_debugOn },
        { "debugStatus",        oCmd_debugStatus },
        { "delete",             oCmd_delete },
        { "deleteBond",         oCmd_deleteBond },
        { "deleteOffLibEntry",  oCmd_deleteOffLibEntry },
        { "deleteRestraint",    oCmd_deleteRestraint },
        { "desc",               oCmd_describe },
        { "deSelect",           oCmd_deSelect },
        { "deSelectMask",       oCmd_deSelectMask },
        { "displayPdbAtomMap",  oCmd_displayPdbAtomMap },
        { "displayPdbResMap",   oCmd_displayPdbResMap },
        { "edit",               oCmd_edit },
        { "flip",               oCmd_flip },
        { "get",                oCmd_get },
        { "groupSelectedAtoms", oCmd_groupSelectedAtoms },
        { "help",               oCmd_help },
        { "impose",             oCmd_impose },
#ifdef DEBUG
        { "lex",                oCmd_lex },
#endif
        { "list",               oCmd_list },
        { "listOff",            oCmd_listOff },
        { "listParmSets",       oCmd_listParmSets },
        { "listResidues",       oCmd_listResidues },
        { "loadAmberParams",    oCmd_loadAmberParams },
        { "loadAmberPrep",      oCmd_loadAmberPrep },
        { "loadOff",            oCmd_loadOff },
        { "loadMol2",           oCmd_loadMol2 },
        { "loadMol3",           oCmd_loadMol3 },
        { "loadCif",            oCmd_loadCif },
        { "loadPdb",            oCmd_loadPdb },
        { "loadPdbUsingSeq",    oCmd_loadPdbUsingSeq },
        { "logFile",            oCmd_logFile },
        { "matchVariables",     oCmd_matchVariables },
        { "measureGeom",        oCmd_measureGeom },
        { "memDebug",           oCmd_memDebug},
        { "quit",               oCmd_quit },
        { "relax",              oCmd_relax },
        { "remove",             oCmd_remove },
        { "restrainAngle",      oCmd_restrainAngle },
        { "restrainBond",       oCmd_restrainBond },
        { "restrainTorsion",    oCmd_restrainTorsion },
        { "resequence",         oCmd_resequence },
        { "saveAmberParm",      oCmd_saveAmberParm },
        { "saveAmberParmNetCDF",oCmd_saveAmberParmNetCDF },
        { "saveAmberParmPert",  oCmd_saveAmberParmPert },
        { "saveAmberParmPol",   oCmd_saveAmberParmPol },
        { "saveAmberParmPolPert", oCmd_saveAmberParmPolPert },
        { "saveAmberPrep",      oCmd_saveAmberPrep },
        { "saveMol2",           oCmd_saveMol2 },
        { "saveMol3",           oCmd_saveMol3 },
        { "saveOff",            oCmd_saveOff},
        { "saveOffParm",        oCmd_saveOffParm },
        { "savePdb",            oCmd_savePdb },
        { "scaleCharges",       oCmd_scaleCharges },
        { "scaleCoor",          oCmd_scaleCoor },
        { "select",             oCmd_select },
        { "selectMask",         oCmd_selectMask },
        { "sequence",           oCmd_sequence },
        { "set",                oCmd_set },
        { "setBox",             oCmd_setBox},
        { "setCell",            oCmd_setCell},
        { "showdefault",        oCmd_showDefault},
        { "solvateBox",         oCmd_solvateBox },
        { "solvateCell",        oCmd_solvateCell },
        { "solvateCap",         oCmd_solvateCap },
        { "solvateDontClip",    oCmd_solvateDontClip },
        { "solvateOct",         oCmd_solvateOct },
        { "solvateShell",       oCmd_solvateShell },
        { "source",             oCmd_source },
        { "translate",          oCmd_translate },
        { "transform",          oCmd_transform },
/*      { "update",             oCmd_update },          */
        { "verbosity",          oCmd_verbosity },
        { "zMatrix",            oCmd_zMatrix },
/* The last command must be blank */
        { "", NULL }
};
