/*
 *      File:   fortran.c
 *
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
 *
 *              This file contains routines for performing
 *              FORTRAN style formatted output.
 *              It contains a routine that allows the caller
 *              to specify the number of entries that should
 *              appear on each line and the format of each
 *              entry on the line.
 *              Then there is a routine that the user uses
 *              to write each subsequent piece of data.
 *
 *		Debugging assistance is given with 
 *			FortranDebugOn()
 *			FortranDebug(str)
 *			GsFortranDebug
 *
 *		When FortranDebugOn is called, subsequent
 *		calls to FortranDebugStr write strings into
 *		the output file that are always ignored by
 *		these fortran input routines, but can
 *		be used by the programmer to get their
 *		bearings within these huge ugly files.
 *		Whenever a FortranDebug string is found
 *		by the reading routines, it is copied into
 *		the global variable GsFortranDebug.
 *
 */



#include	"basics.h"
#include	"fortran.h"

#include        "stringExtra.h"


/*
 *--------------------------------------------------------------------
 *
 *      Static variables 
 *
 */


#define	FORTRAN_DEBUG_COMMENT_CHAR	')'


static  STRING  SsFormat;
static  FILE*   SfFile;
static  int     SiPerLine;
static  int     SiOnLine = 0;   /* Stores the # of entries already on line */
static  bool    SbWroteNothing;
static  int     SiIntWidth;

static  STRING  SsInput;
static  bool    SbNeedInput = true;

static	bool	SbFortranDebug = false;
STRING	GsFortranDebug;



/*
 *----------------------------------------------------------
 *
 *	Private routines
 *
 */

/* 
 *	zFortranGetInputLine
 *
 *	Read the next line into the STATIC variable SsInput.
 */
static void
zFortranGetInputLine()
{
   while ( !feof(SfFile) ) {
	SsInput[0] = '\0';
	fgets( SsInput, sizeof(SsInput), SfFile );
	if ( SsInput[0] == FORTRAN_DEBUG_COMMENT_CHAR ) {
	    strcpy( GsFortranDebug, &(SsInput[1]) );
	    MESSAGE("FORTRAN DEBUG:%s", GsFortranDebug );
	    continue;
	}
	break;
   }
}



/*
 *==========================================================
 *
 *	Public routines
 *
 */


/*
 *      FortranFile
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Save the file where output is to be placed.
 */
void
FortranFile( FILE *fOut )
{
    SbNeedInput = true;
    SfFile = fOut;
}

static const char *ScPSection = NULL;

void FortranSection(const char *label) {
  char line[80];
  ScPSection = label;
  sprintf(line,"%%FLAG %s",label);
  FortranFormat(1, "%-80s");
  FortranWriteString(line);
}

void FortranIntFormat(int maxval) {
  char fmt[32],line[80];
  int width, count;
  width = snprintf(NULL, 0, "%d", maxval);
  snprintf(fmt, sizeof(fmt), "%%%dd", width);  /* e.g. "%4d" */
  count = 80 / width;                 /* ints per line */
  if (GDefaults.dPrmtopFormat>1.0) {
    sprintf(line,"%%COMMENT CFORMAT:%s;", fmt);
    FortranWriteString(line);
  }
  sprintf(line,"%%FORMAT(%dI%d)",count,width);
  FortranWriteString(line);
  FortranFormat(count, fmt);
  SiIntWidth = width;
}



/*
 *      FortranFormat
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Save the format and number of entries per line.
 */
void
FortranFormat( int iPerLine, char *sFormat )
{
    SiPerLine = iPerLine;
    strcpy( SsFormat, sFormat );
    SiOnLine = 0;
    SbWroteNothing = true;
}



/*
 *      FortranWriteInt
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Write an integer out to the file, increment the number of entries on 
 *      the line, if it is equal to SiPerLine, then move to the next line.
 */
void
FortranWriteInt( int iVal )
{
    static const char *sWarnSection=NULL;
    int n = fprintf( SfFile, SsFormat, iVal );
    if (n > SiIntWidth) {
        // warn on ove
        if (sWarnSection != ScPSection)
            VPFATAL("Integer overflow in %s, writing value %d with format width %d\n", ScPSection, iVal, SiIntWidth);
        sWarnSection = ScPSection;
    }
    SiOnLine++;
    SbWroteNothing = false;
    if ( SiOnLine >= SiPerLine ) {
        fprintf( SfFile, "\n" );
        SiOnLine = 0;
    }
}

  



/*
 *      FortranWriteDouble
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Write a double out to the file, increment the number of entries on 
 *      the line, if it is equal to SiPerLine, then move to the next line.
 */
void
FortranWriteDouble( double dVal )
{
    fprintf( SfFile, SsFormat, dVal );
    SiOnLine++;
    SbWroteNothing = false;
    if ( SiOnLine >= SiPerLine ) {
        fprintf( SfFile, "\n" );
        SiOnLine = 0;
    }
}




/*
 *      FortranWriteString
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Write a String out to the file, increment the number of entries on 
 *      the line, if it is equal to SiPerLine, then move to the next line.
 */
void
FortranWriteString( char *sVal )
{
    fprintf( SfFile, SsFormat, sVal );
    SiOnLine++;
    SbWroteNothing = false;
    if ( SiOnLine >= SiPerLine ) {
        fprintf( SfFile, "\n" );
        SiOnLine = 0;
    }
}




/*
 *      FortranEndLine
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      If the number of objects on the line is not zero then print
 *      and end of line.
 */
void
FortranEndLine()
{
    if ( SbWroteNothing || SiOnLine != 0 ) 
	fprintf( SfFile, "\n" );
    SbWroteNothing = true;
    SiOnLine = 0;
}


/*
 *      sFortranReadString
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Read the next string from the file
 */
char *
sFortranReadString( char *sString )
{
bool            bGotOne;
int             iLen;

    bGotOne = false;
    strcpy( sString, "" );
    while ( !bGotOne ) {
        if ( SbNeedInput ) {
            if ( feof(SfFile) ) return NULL;
            strcpy( SsInput, "" );
	    zFortranGetInputLine();
            if ( (iLen=strlen(SsInput)) > 0 ) SsInput[iLen-1] = '\0';
        }
        SbNeedInput = false;
        sRemoveLeadingSpaces( SsInput );
        sRemoveFirstString( SsInput, sString );
        if ( strlen(sString)!=0 ) bGotOne = true;
        SbNeedInput = ( strlen(SsInput)==0 );
   }
   return sString;
}




/*
 *      sFortranReadLabel
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Read the next label from the file
 */
char *
sFortranReadLabel( char *sString )
{
bool            bGotOne;
int             iLen;

    bGotOne = false;
    strcpy( sString, "" );
    while ( !bGotOne ) {
        if ( SbNeedInput ) {
            if ( feof(SfFile) ) return NULL;
	    zFortranGetInputLine();
            if ( (iLen=strlen(SsInput)) > 0 ) SsInput[iLen-1] = '\0';
        }
        SbNeedInput = false;
        memcpy( sString, SsInput, 4 );
        sString[4] = '\0';
        strcpy( SsInput, SsInput+4 );
        if ( strlen(sString)!=0 ) bGotOne = true;
        SbNeedInput = ( strlen(SsInput)==0 );
   }
   return sString;
}



/*
 *      iFortranReadInt
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Read the next string and return it as an integer.
 */
int
iFortranReadInt( )
{
STRING          sInt;
int             iVal;

    sFortranReadString( sInt );
    sscanf( sInt, "%d", &iVal );
    return iVal;
}




/*
 *      dFortranReadDouble
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Read the next string and return it as a double.
 */
double  
dFortranReadDouble( )
{
STRING          sDouble;
double          dVal;

    sFortranReadString( sDouble );
    sscanf( sDouble, "%lf", &dVal );
    return dVal; 
}



/*
 *      FortranSkipLine
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Skip a line of input
 */
void
FortranSkipLine()
{
    if ( feof(SfFile) ) 
	return;
    zFortranGetInputLine();
}




/*
 *	FortranDebugOn
 *
 *	Turn fortran debugging ON.  Default is OFF.
 *
 *	When on, FortranDebug calls write single strings
 *	Into the fortran format file which will be ignored
 *	by fortran reading routines, but can be
 *	used by the programmer to get their bearings within
 *	the nasty looking files.
 */
void
FortranDebugOn()
{
    SbFortranDebug = true;
}

/*
 *	FortranDebug
 *
 *	If SbFortranDebug is true then write the
 *	string to the output file with a comment character
 *	prefixed to it.
 */
void
FortranDebug( char *sStr )
{
    if ( SbFortranDebug ) {
        fprintf( SfFile, "%c%s\n", FORTRAN_DEBUG_COMMENT_CHAR, sStr );
    }
}
