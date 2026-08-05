/*
 *      File: basics.c
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
 *              Declare global variables, and some routines which
 *              are always required.
 */



#include	"basics.h"
#include	"defaults.h"
#include	"zMatrix.h"
#include        "cmap.h"

#include	<signal.h>
#include	<stdarg.h>
#include	<stdlib.h>


STRING	GsBasicsFullName;

bool	GbInterrupt = FALSE;

int	GiUnitEditors = 0;

int     GiVPerrorCount = 0;   /* Count the number of Fatal errors */
int     GiVPwarningCount = 0; /* Count the number of Warnings */
int     GiVPnoteCount = 0;    /* Count the number of Notes */

/*
 * ------------------------------------------------------------------
 *
 *      typedefs
 *
 *	For Fail-safe memory management.
 *
 */

#define FILENAMELEN     40
#define TRAILERLEN      100
                        /* MUST not be longer than TRAILERLEN+1 */
#define CHECKSTR        \
"usedUSEDThis is a buffer at the end of the object to check for overwrites"
                        /* MUST not be longer than TRAILERLEN+1 */
#define FREESTR         \
"FREEfreeNow this is free also and will be free until it is malloced again" 


/*
 *      myAcos
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the acos of a number, but don't crap
 *      out on domain errors.
 */
double  
myAcos( double d )
{
	if ( d >= 1.0 ) 
		return(0.0);
	if ( d <= -1.0 )
		return(PI);
	return(acos(d));
}



/*
 *	myPow
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return the power of a number, but don't crap out
 *	on possible underflow errors.
 */
double	
myPow( double x, double y )
{
	if ( fabs(x) < VERYSMALL ) 
		return(0.0);
	return(pow(x,y));
}




/*
 *      iDoubleCompare
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Compare two double precision values.
 *
 *      The values are considered to be the same if
 *      the ranges    [dA-TOLERANCE,dA+TOLERANCE] contains dB.
 *
 *      Return 0 if they are the same, a value <0 if dA<dB and
 *      a value >0 if dA>dB.  Just like strcmp.
 */
int     
iDoubleCompare( double dA, double dB )
{
	if ( (dA-TOLERANCE<dB) && (dA+TOLERANCE>dB) )
		return(0);
	if ( dA < dB )
		return(-1);
	return(1);
}



/*
 *	bStringToDouble
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return TRUE if the string passed to this routine is
 *	could be completely converted into a double precision
 *	value.  Only change the value of *dPData if the string
 *	completely represents a double precision value.
 */
bool	
bStringToDouble( char *cPData, double *dPData )
{
	char	*cPEnd;
	double	dValue;

	dValue = (double)strtod( cPData, &cPEnd );
	if ( cPEnd - cPData == (long int)strlen(cPData) ) {
		*dPData = dValue;
		return(TRUE);
	}
	return(FALSE);
}

/*
 *	bStringToInt - By David Rivkin
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return TRUE if the string passed to this routine is
 *	could be completely converted into a integer
 *	value.  Only change the value of *dPData if the string
 *	completely represents a integer value.
 */
bool	
bStringToInt( char *cPData, int *iPData )
{
	char	*cp;

	for (cp=cPData; *cp; cp++) {
		if (*cp == '-') {
			if (cp != cPData) {
				VP0("'-' embedded in number %s\n", cPData );
				return(FALSE);
			}
			continue;
		}
		if (!isdigit((unsigned char)*cp)) {
			VP0("non-digit in %s\n", cPData );
			return(FALSE);
		}
		
	}
	*iPData = atoi( cPData );
	return(TRUE);
}



/*
 *      StringLower
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Convert the string to lower case.
 *      Only attempt to modify letters that are uppercase
 *      This allows const char* sStr if already lowercase.
 */
void
StringLower( char *sStr )
{
    while ( *sStr != '\0' ) {
        char c = cLower(*sStr);
        if (*sStr != c) *sStr = c;
        sStr++;
    }
}

/*
 *      StringUpper
 *
 *	Author:	Juno Krahn (2026)
 *
 *      Uppercase version of StringLower
 */
void
StringUpper( char *sStr )
{
    while ( *sStr != '\0' ) {
        char c = cUpper(*sStr);
        if (*sStr != c) *sStr = c;
        sStr++;
    }
}

/*
 *      StringTrim
 *
 *	Author:	Juno Krahn (2026)
 *
 *      Trim whitespace at beginning and end.
 */
void
StringTrim(char *sStr)
{
    char *start = sStr;
    while (*start && isspace((unsigned char)*start)) start++;
    if (start != sStr) memmove(sStr, start, strlen(start) + 1);
    StringRTrim(sStr);
}


/*
 *      StringRTrim
 *
 *	Author:	Juno Krahn (2026)
 *
 *      Trim whitespace at end.
 */
void
StringRTrim(char *sStr)
{
    char *end = sStr + strlen(sStr);   /* points at NUL */
    while (end > sStr && isspace((unsigned char)end[-1])) end--;
    *end = '\0';
}


//---------------------------------------------------------------
//      Memory status via rusage, GLIBC or libasan
#ifdef __GLIBC__
#  include <malloc.h>
#endif
#if defined(__SANITIZE_ADDRESS__)
#  ifdef __has_include
#    if __has_include(<sanitizer/allocator_interface.h>)
#      include <sanitizer/allocator_interface.h>
#      define USE_MODERN_SAN_HEADERS 1
#    endif
#  endif
#  ifndef USE_MODERN_SAN_HEADERS
     extern size_t __sanitizer_get_current_allocated_bytes(void);
     extern size_t __sanitizer_get_heap_size(void);
     extern size_t __sanitizer_get_free_bytes(void);
#  endif
#endif
#include <sys/resource.h>

void
PrintMemoryStats(void) {
#ifdef __GLIBC__
    struct mallinfo mi = mallinfo();
    VP0("Current memory usage:%9ld KB\n", mi.uordblks / 1024 );
#endif
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        long kb = usage.ru_maxrss;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
        kb = kb / 1024;
#endif
        // On Linux, ru_maxrss is in Kilobytes.
        VP0("Peak RAM Usage (HWM):%9ld KB\n", kb);
    }
#if defined(__SANITIZE_ADDRESS__)
    VP0("Built with Address Sanitization library\n");
    size_t current_heap = __sanitizer_get_current_allocated_bytes();
    size_t free_bytes = __sanitizer_get_free_bytes();
    size_t heap_size = __sanitizer_get_heap_size();    // Current total bytes allocated by the application via ASan's wrapper
    VP0("--- ASan Memory Snapshot ---\n");
    VP0("Actively Used Heap:  %9zu KB\n", current_heap / 1024);
    VP0("Total Managed Pool:  %9zu KB\n", heap_size / 1024);
    VP0("Cached Free List:    %9zu KB\n", free_bytes / 1024);
#endif
}


/*
 *---------------------------------------------------------------
 *
 *      Message printing.
 *
 *      Maintain and look up source files that currently display
 *      MESSAGE statements.
 */
 
static  int     SiMessageFiles = 0;
static  STRING  SsaMessageFiles[MAXMESSAGEFILES];


 

/*
 *      bMessageCheck
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Check if a message should be printed.
 */
bool    
bMessageCheck( char *sFile )
{
int     i;
    if ( SiMessageFiles == 0 ) 
	return FALSE;

    for ( i=0; i<SiMessageFiles; i++ ) {
        if ( SsaMessageFiles[i][0] == '*' ) 
		return TRUE;
        if ( strncmp( sFile, SsaMessageFiles[i], 
                        strlen(SsaMessageFiles[i]) ) == 0 ) 
		return TRUE;
    }
    return FALSE;
}


/*
 *      MessageAddFile
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add a file to the list of files which should produce
 *      messages.
 */
void    
MessageAddFile( char *sFile )
{

    if ( SiMessageFiles >= MAXMESSAGEFILES ) 
	return;
    strcpy( SsaMessageFiles[SiMessageFiles], sFile );
    SiMessageFiles++;

    PRINTF("TURNING ON MESSAGES FROM: <%s>\n", sFile );
}


/*
 *      MessageRemoveFile
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Remove a file from the list of files that should produce
 *      messages.
 */
void    
MessageRemoveFile( char *sFile )
{
int             i, j;

                /* Find the file */
    for ( i=0; i<SiMessageFiles; i++ ) {
        if ( strcmp( SsaMessageFiles[i], sFile )== 0 ) {
            for ( j=i;j<SiMessageFiles-1;j++ ) {
                strcpy( SsaMessageFiles[j], SsaMessageFiles[j+1] );
            }
            SiMessageFiles--;
            break;
        }
    }
    PRINTF("TURNING OFF MESSAGES FROM: <%s>\n", sFile );
}


/*
 *      MessageFileList
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Print a list of the files which will generate messages.
 */
void    
MessageFileList(void)
{
	int	i;


    for ( i=0; i<SiMessageFiles; i++ ) {
        myPrintf(0, "%s\n", SsaMessageFiles[i] );
    }
    myPrintf(0, "------\n" );
}    


#ifdef DEBUG
/*
 *	MessageInitialize
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Initialize message printing from 
 *	those files mentioned in the MESSAGEON environment variable.
 *
 *	The MESSAGEON variable contains a sequence of names
 *	separated by spaces.
 */
static void	
MessageInitialize(void)
{
STRING	sName;
char	*cPGet, *cPPut;

    PRINTF("Reading MESSAGEON environment variable.\n" );
    cPGet = (char*)getenv("MESSAGEON");

    if ( cPGet == NULL ) 
	return;

    while ( *cPGet ) {
        cPPut = sName;
	while ( *cPGet && *cPGet != ' ' ) {
	    *cPPut = *cPGet;
	    cPPut++;
	    cPGet++;
	}
	*cPPut = '\0';
	if ( *cPGet ) cPGet++;
	MESSAGEON(sName);
    }
}

void
myPrintTrace( const char *prefix, int depth, const char *fmt, ... )
{
    char    sBuf[MAXCHARSPERPRINTF];
    char    sMsg[MAXCHARSPERPRINTF];
    va_list args;

    va_start( args, fmt );
    vsnprintf( sBuf, sizeof(sBuf), fmt, args );
    va_end( args );

    snprintf( sMsg, sizeof(sMsg),
              "Trace: %s %.*s from call depth %2d.\n",
              prefix, (int)(MAXSTRINGLENGTH - 32), sBuf, depth );

    myPrintString( sMsg, TRACE_VERBOSITY );
}
#endif



/*
 *=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
 *
 *	File handling routines.
 *
 */
#define	MAXDIRECTORIES	10

static STRING		SsaDirectories[MAXDIRECTORIES];
static int		SiCurDirectory = 0;

static int	
iExpandDir( char *sExpanded, char *sOriginal )
{
#if (defined WIN32)
    strcpy( sExpanded, sOriginal );
    return(0);
#else
char		user[100];
int		i;
struct passwd 	*pw, *getpwnam();

/* 
TODO: add string size protection
*/
    if ( sOriginal[0] == '$' ) {
        for (i=0; i<99; i++) {
	    if ( sOriginal[i+1] == '\0' || sOriginal[i+1] == '/' ) {
		user[i] = '\0';
		break;
	    }
	    user[i] = sOriginal[i+1];
        }
        i++;
        char *var = getenv(user);
        if ( var == NULL ) {
	    VP0("Could not get environment value for: %s\n", user);
            return 1;
        }
        strcpy( sExpanded, var );
        if ( sOriginal[i] == '/' )
	    strcat( sExpanded, &sOriginal[i] );
        return(0);
    }
    if ( sOriginal[0] != '~' ) {
	/*
	**	nothing to expand
	*/
    	strcpy( sExpanded, sOriginal );
    	return(0);
    }

    /*
     *  there is a leading ~ which needs to be expanded
     */

    if ( sOriginal[1] == '\0'  ||  sOriginal[1] == '/' ) {

	/*
	 *  relative to user's home dir
	 */

	char	*home;

	home = (char*)getenv("HOME");
	if ( home == (char*)NULL ) {
		VP0("~: Could not get HOME from environment\n" );
		return(1);
	}
	strcpy( sExpanded, home);
	if ( sOriginal[1] == '/' )
		strcat( sExpanded, &sOriginal[1] );
	return(0);
    }

    /*
     *  relative to a specified user's home dir - get user name &
     *	look it up in pwd file
     */

    for (i=0; i<99; i++) {
	if ( sOriginal[i+1] == '\0' || sOriginal[i+1] == '/' ) {
		user[i] = '\0';
		break;
	}
	user[i] = sOriginal[i+1];
    }
    i++;
    pw = getpwnam( user );
    if ( pw == NULL ) {
	VP0("%s: Could not get from password file: %s\n", user,
						strerror(errno) );
	return(1);
    }

    strcpy( sExpanded, pw->pw_dir );
    if ( sOriginal[i] == '/' )
	strcat( sExpanded, &sOriginal[i] );
	
    return(0);
#endif
}

/*
 *	BasicsAddDirectory
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Add a directory to the list of directories to search
 *	for a file when it cannot be opened in the current directory.
 */
int		
BasicsAddDirectory( STRING sDirectory, int bomb )
{
STRING		sTmpDir;
FILESTATUSt	fsStatus;

    if ( SiCurDirectory + 1 >= MAXDIRECTORIES ) {
	VP0("Limit reached on include directories - can't add dir\n" );
	if ( bomb )
		exit(1);
	return(0);
    }
    if ( iExpandDir( sTmpDir, sDirectory ) ) {
	if ( bomb )
		exit(1);
	return(0);
    }
    fsStatus = fsSysdependFileStatus(sTmpDir);
    if ( !(fsStatus.fMode & FILEDIRECTORY) ) {
	if (errno)
		VP0("%s: %s\n", sTmpDir, strerror(errno) );
	else
		VP0("%s: not a directory %s\n", sTmpDir );
	if ( bomb )
		exit(1);
	return(0);
    }
    strcpy( SsaDirectories[SiCurDirectory], sTmpDir );
    SiCurDirectory++;
    return(1);
}







/*
 *	fBasicsMyFopen
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	This routine is a wrapper to the library fopen routine.
 *	If the path is absolute (starts w/ root or '.'), try it. 
 *	If not, try the current directory then search through the 
 *	list of directories in SvaDirectories.
 */
FILE *
fBasicsMyFopen( char *sFilename, char *sAttributes, bool bComplain )
{
FILE		*fFile;
int		i, iExistErr;
FILESTATUSt	fsStatus;
STRING		sExpanded;
int absolutePath;

    VPTRACEENTER("fBasicsMyFopen" );
    VPTRACEMULTIPLEEXIT("fBasicsMyFopen" );
    if ( strlen(sFilename) == 0 ) 
	return(NULL);

    if ( strlen(sFilename) == 1 && sFilename[0] == '-' ) {
	VP2("Reading from standard input.\n" );
	return(stdin);
    }

    if ( iExpandDir( sExpanded, sFilename ) )
	return(NULL);
    fsStatus = fsSysdependFileStatus(sExpanded);

    if ( fsStatus.fMode & FILEDIRECTORY ) {
	VP0("%s is a directory\n", sExpanded );
	return(NULL);
    }

    #if (defined WIN32)
    absolutePath = sExpanded[0] == 'C' && sExpanded[1] == ':' && sExpanded[2] == '\\';
    #else
    absolutePath = sExpanded[0] == '/' || sExpanded[0] == '.';
    #endif

    if (absolutePath) {
	/* (this is where an sExpanded ~ dir would fall) */
        strcpy ( GsBasicsFullName, sExpanded );
        fFile = fopen( sExpanded, sAttributes );
	if ( fFile == NULL  &&  bComplain ) {
		VPFATALEXIT("Could not open file %s: %s\n", 
						sExpanded, strerror(errno) );
	}
	return(fFile);
    }
    int len = snprintf ( GsBasicsFullName, sizeof(GsBasicsFullName), "./%s", sExpanded );
    if (len >= (int)sizeof(GsBasicsFullName) ) fFile = NULL;
    else fFile = fopen( sExpanded, sAttributes );
    if ( fFile != NULL ) 
	return(fFile);
    iExistErr = 1;
    for ( i=0; i<SiCurDirectory; i++ ) {
	strcpy( GsBasicsFullName, SsaDirectories[i] );
	strcat( GsBasicsFullName, "/" );
	strcat( GsBasicsFullName, sExpanded );
	fFile = fopen( GsBasicsFullName, sAttributes );
	if ( fFile != NULL ) 
	    return(fFile);
	if ( errno != ENOENT ) {
	    VP0("Opening %s: %s\n", GsBasicsFullName, strerror(errno) );
	    iExistErr = 0;
	}
    }
    if ( bComplain ) {
	VPFATALEXIT("Could not open file %s: %s\n", sExpanded,
			( iExistErr ? "not found" : "system error" )  );
    }
    return(NULL);
}




/*
 *------------------------------------------------------------------
 *
 *	Handle multiple output sinks.
 *
 *	The current output sink is the one where all output from
 *	the VPx commands is sent.
 *
 */


typedef	struct	{
	bool		bSinkUsed;
	bool		bPrintPrefix;
	STRING		sPrefix;
	VFUNCTION	fCallback;
	GENP		PData;
} SINKINFOt;

static	int		SiNumberOfSinks = 0;
static	SINKINFOt	*SsiPSinks = NULL;

#define	MAX_SINK_STACK		30

static int		SiaSinkStack[MAX_SINK_STACK];
static int		SiNextSink = 0;
static int		iCurrentPrintSink = -1;

/*
 *	iCreatePrintSink
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Create a new print sink and return an integer
 *	handle for the sink that can be used by subsequent
 *	calls to DefineCurrentPrintSink or DestroyPrintSink.
 *
 *	The PData field can be used to pass additional info to
 *	the output callback, like a Widget.
 *
 *	The Callback function has the following form.
 *
 *	callback( cPString, PData )
 *
 */
int
iCreatePrintSink( VFUNCTION fOutputCallback, char *sPrefix, GENP PData )
{
int		i;
bool		bFoundOne;

	/* First check if there isn't a free sink */

    bFoundOne = FALSE;
    for ( i=0; i<SiNumberOfSinks; i++ ) {
	if ( !SsiPSinks[i].bSinkUsed ) {
	    bFoundOne = TRUE;
	    break;
	}
    }

		/* If one is not found then allocate space for one */

    if ( !bFoundOne ) {
	if ( SsiPSinks == NULL ) {
	    i = 0;
	    SsiPSinks = (SINKINFOt*)MALLOC(sizeof(SINKINFOt) );
	    SiNumberOfSinks = 1;
	} else {
	    i = SiNumberOfSinks;
	    SsiPSinks = (SINKINFOt*)REALLOC(SsiPSinks, sizeof(SINKINFOt)*(SiNumberOfSinks+1) );
	    SiNumberOfSinks++;
	    /*
	     *  realloc can mean that globals pointing to the
	     *	current print sink now point to freed memory
	     */
	    if (iCurrentPrintSink == -1)
		DFATAL("iCurrentPrintSink == -1\n");
	    GcPPrefix = SsiPSinks[iCurrentPrintSink].sPrefix;
	    GbPrintPrefix = SsiPSinks[iCurrentPrintSink].bPrintPrefix;
	    GfPrintStringCallback = SsiPSinks[iCurrentPrintSink].fCallback;
	    GPData = SsiPSinks[iCurrentPrintSink].PData;
	}
    }

		/* Define the information */

    SsiPSinks[i].bSinkUsed = TRUE;
    SsiPSinks[i].bPrintPrefix = TRUE;
    SsiPSinks[i].fCallback = fOutputCallback;
    strcpy( SsiPSinks[i].sPrefix, sPrefix );
    SsiPSinks[i].PData = PData;
    MESSAGE("Created print sink: %d\n", i );

    return(i);
}



GENP	GPData;

/*
 *	zDefineCurrentPrintSink
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Using the handle passed in iHandle, define the
 *	current print sink.
 */
static void
zDefineCurrentPrintSink( int iHandle )
{
    iCurrentPrintSink = iHandle;

    if ( iHandle >= SiNumberOfSinks || iHandle < 0 ) {
	DFATAL("Invalid print sink handle: %d, must be in [0..%d]\n",
			iHandle, SiNumberOfSinks-1 );
    }

    MESSAGE("Selecting print sink: %d\n", iHandle );

    if ( SsiPSinks[iHandle].bSinkUsed ) {
	GcPPrefix = SsiPSinks[iHandle].sPrefix;
	GbPrintPrefix = SsiPSinks[iHandle].bPrintPrefix;
	GfPrintStringCallback = SsiPSinks[iHandle].fCallback;
	GPData = SsiPSinks[iHandle].PData;
    } else {
	DFATAL("Unused print sink: %d\n", iHandle );
    }
}


/*
 *	DestroyPrintSink
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Destroy the print sink indicated by iHandle.
 */
void	
DestroyPrintSink( int iHandle )
{
    if ( iHandle >= SiNumberOfSinks || iHandle < 0 ) {
	DFATAL("Invalid print sink handle: %d, must be in [0..%d]\n",
			iHandle, SiNumberOfSinks-1 );
    }

    MESSAGE("Destroying print sink: %d\n", iHandle );

    if ( SsiPSinks[iHandle].bSinkUsed ) {
	SsiPSinks[iHandle].bSinkUsed = FALSE;
    } else {
	DFATAL("Unused print sink: %d\n", iHandle );
    }
}



/*
 *	PushCurrentPrintSink
 *
 *	Push the sink in (iSink) onto the stack as the
 *	current print sink.
 */
void
PushCurrentPrintSink( int iSink )
{
    if ( SiNextSink >= MAX_SINK_STACK ) {
	DFATAL("Exhausted Print Sink Stack\n" );
    }
    SiaSinkStack[SiNextSink] = iSink;
    SiNextSink++;
    zDefineCurrentPrintSink(iSink);
}



/*
 *	PopCurrentPrintSink
 *
 *	Pop the top print sink from the stack.
 */
void
PopCurrentPrintSink(void)
{
	if ( SiNextSink <= 0 ) {
		DFATAL("Underflow in Print Sink Stack\n" );
	}
	SiNextSink--;
	if ( SiNextSink == 0 ) {
	    GfPrintStringCallback = NULL;
	} else {
	    zDefineCurrentPrintSink(SiaSinkStack[SiNextSink-1]);
	}
}





/*
 *------------------------------------------------------------------
 *
 *      My own printf which checks verbosity levels
 *      and writes output to an optional log file.
 */

	/* Define variables used by myPrintf which maintain what */
	/* function to call to print a string on the display */
	/* and a buffer that stores the output to be printed */
	/* By default myPrintf will write its non-log output to */
	/* stdout using 'puts'.  But this can be re-directed using */
	/* DefinePrintStringCallback(func(char*)). */
	/* This is useful for output to Widgets */


FILE	*GfLog = NULL;
int     GiTraceIndentationLevel = 0;
int     GiVerbosityLevel = 0;
int     GiVerbosity;                    /* This changes for every P# */
bool	GbPrintPrefix = TRUE;
char	*GcPPrefix = NULL;


/*
 *	myPuts
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Put a string to the stdout without appending a newline to the end.
 *	PData is not used.
 */
static void	
myPuts( char *sLine, GENP PData )
{
    printf( "%s", sLine );
}

void	(*GfPrintStringCallback)() = (VFUNCTION)myPuts;


/*
 *  myPrintString
 *
 *  Author: Christian Schafmeister (1991) — updated (2025)
 *
 *  Print cPString to appropriate devices.
 *  Breaks on \n characters; prefixes lines with GcPPrefix when enabled.
 *  iVerbosity replaces the former GiVerbosity global.
 */
void
myPrintString(const char *cPString, int iVerbosity)
{
    const char  *cPStart, *cPStop;
    char         sTempBuf[MAXCHARSPERPRINTF];
    char        *cPPrint;
    bool         bPrintPrefix;
    size_t       len;

    cPStart = cPString;
    cPStop  = cPString;

    while (*cPStop) {
        /* Advance cPStop to next \n or end of string */
        while (*cPStop && *cPStop != '\n')
            cPStop++;

        if (*cPStop == '\n') {
            /* Copy segment including the \n, null-terminate */
            len = (size_t)(cPStop - cPStart) + 1;   /* include \n */
            if (len >= sizeof(sTempBuf))
                len = sizeof(sTempBuf) - 1;          /* truncate safely */
            memcpy(sTempBuf, cPStart, len);
            sTempBuf[len] = '\0';
            cPPrint      = sTempBuf;
            bPrintPrefix = TRUE;
            cPStart      = ++cPStop;                 /* next segment */
        } else {
            /* End of string — print remainder as-is (no copy needed) */
            cPPrint      = (char *)cPStart;          /* cast: read-only use */
            bPrintPrefix = FALSE;
        }

        /* Emit prefix if required */
        if (GbPrintPrefix && GcPPrefix != NULL) {
            if (GfPrintStringCallback == NULL) {
                DFATAL("Invalid print string callback!");
            }
            if (GiVerbosityLevel >= iVerbosity)
                GfPrintStringCallback(GcPPrefix, GPData);
            if (GfLog != NULL && iVerbosity >= 0)
                fputs(GcPPrefix, GfLog);
        }

        GbPrintPrefix = bPrintPrefix;

        if (GiVerbosityLevel >= iVerbosity)
            GfPrintStringCallback(cPPrint, GPData);

        if (GfLog != NULL && iVerbosity >= 0)
            fputs(cPPrint, GfLog);
    }
}

static struct timespec StsElapsedStart;
static clock_t SctCPUStart;
/*
 *	BasicsInitialize
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Initialize basic routines.
 */
void	
BasicsInitialize(void)
{
        clock_gettime(CLOCK_MONOTONIC, &StsElapsedStart);
        SctCPUStart=clock();
	(void)zMatrixInit();
        InitializeDefaults();
/*
    signal( SIGSEGV, zBasicsTrapSEGV );
    signal( SIGBUS, zBasicsTrapBUS );
*/
#ifdef	DEBUG
    MessageInitialize();
#endif

} 

void	
BasicsFinalize(void)
{
    if (iVerbosity()>1) PrintMemoryStats();
    if (!GDefaults.bTiming) return;
    struct timespec tsElapsedEnd;
    clock_t ctCPUEnd = clock();
    clock_gettime(CLOCK_MONOTONIC, &tsElapsedEnd);
    double dt = (double)(tsElapsedEnd.tv_sec - StsElapsedStart.tv_sec)
              + (double)(tsElapsedEnd.tv_nsec - StsElapsedStart.tv_nsec) / 1e9;
    double dCPU = ((double) (ctCPUEnd - SctCPUStart)) / (double)CLOCKS_PER_SEC;
    int elapsed_m   = (int)(dt / 60.0);
    double elapsed_s = fmod(dt, 60.0);
    int cpu_m   = (int)(dCPU / 60.0);
    double cpu_s = fmod(dCPU, 60.0);
    VP0("Elapsed time: %d:%05.2fs   CPU time: %d:%05.2fs\n",
           elapsed_m, elapsed_s, cpu_m, cpu_s);
}


/*
 *	BasicsKlassMisMatchPanic
 *
 *	A Klass mismatch has occured so stop.
 */
void	
BasicsKlassMismatchPanic( GENP PObj, char *cPFile, int iLine )
{
    printf( "ERROR: Klass mismatch in file: %s  line: %d\n",
		cPFile, iLine );
    if ( !PObj ) {
	printf( "Object is NULL!\n" );
    } else {
	printf( "Object is not NULL but of unknown type.\n" );
    }

    abort();
}
