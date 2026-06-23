/*
 *      File:   basics.h
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
 *              Define all basic symbols, which must always be
 *              present in my code.
 *
 *              There are two macros that have a profound effect
 *              on what happens in this include file.
 *
 *              DEBUG:  If DEBUG is defined then all debugging macros
 *                      are defined, including memory debugging.
 *                      If NOT defined then all debugging macros
 *                      are defined as nothing.
 *
 *              MEMORY_DEBUG:   Can have one of the following five values.
 *                              If MEMORY_DEBUG is anything other than 0
 *                              then DEBUG must be defined.
 *              0 -     Dont do any MEMORY debugging.
 *              1 -     Use my memory usage tracking.
 *                              This adds 4 bytes to every MALLOC and stores
 *                              the size of the allocation in an int
 *                              at the front of the memory block.
 *                              This is useful for finding memory leaks.
 *              2 -     Use my memory debugging.
 *                              Every MALLOC,FREE,REALLOC first tests
 *                              memory for overwrites etc.
 *              3 -     Use SUN memory debugging, requires linking with
 *                              /usr/lib/debug/malloc.o
 *              4 -     Use my memory debugging, but test memory EVERY TIME
 *                              a routine is entered or exited.
 */


#ifndef BASICS_H
#define BASICS_H

#include        <stdio.h>
#include        <math.h>
#include        <ctype.h>
#include        <time.h>
#include        <errno.h>
#include        <string.h>
#include        <stdlib.h>

// Windows compatibility
#ifdef WIN32

// the windows headers have some other ATOM and ERROR types that have to, as we say in the business, go for a swim in concrete boots

# define ATOM WIN32_BUILTIN_ATOM
# define ERROR WIN32_BUILTIN_ERROR

# include <Windows.h>
# undef ATOM
# undef ERROR

# define strcasecmp _stricmp
#else
# include <unistd.h>
# include <strings.h>
#endif


#ifdef CMAKE
// use new config.h-based logic

#include		<config.h>

// headers
// ---------------------------------------------------------------------

#ifdef HAVE_STRING_H
# include <string.h>
#endif

#ifdef HAVE_SYS_PARAM_H
# include <sys/param.h>
#endif

#ifdef HAVE_STRINGS_H
# include <strings.h>
#endif

#ifdef HAVE_DIRENT_H
# include <dirent.h>
#endif

#ifdef HAVE_SYS_DIR_H
# include <sys/dir.h>
#endif

#ifdef HAVE_PWD_H
# include <pwd.h>
#endif

#ifdef HAVE_STDBOOL_H
# include <stdbool.h>
#endif

#ifdef HAVE_SYS_TIME_H
# include <sys/time.h>
#endif

// types
// ---------------------------------------------------------------------

#ifdef BOOL_IN_X11_XMD
# include <X11/Xmd.h>
#else
# if HAVE_BOOL || defined(BOOL)
// do nothing
# elif defined(HAVE__BOOL)
typedef _Bool BOOL;
# else
typedef unsigned char   BOOL;
# endif
#endif

#if defined(WIN32) && !defined(MAXPATHLEN)
# define MAXPATHLEN MAX_PATH
#endif


#else /*CMAKE*/
// use creaky platform-name-based logic

// headers
// ---------------------------------------------------------------------

#if (!defined WIN32)
# include <pwd.h>
#endif

#if (defined(__unix__) || defined(unix)) && !defined(USG)
# include <sys/param.h>
#endif

#if (defined WIN32)
# include <dirent.h>
# include <strings.h>
#elif (defined SYSV || defined __i386__)
# include        <dirent.h>
# include        <string.h>
#else
# include		<sys/dir.h>
# include        <strings.h>
#endif

#include	<sys/time.h>

#ifdef __APPLE__
  // for PATH_MAX
# include <sys/syslimits.h> 
#endif

// types
// ---------------------------------------------------------------------

        /* BOOL is defined on NeXT appkit/appkit.h */
#if !defined(APPKIT_H) && !defined(XtSpecificationRelease)
typedef unsigned char   BOOL;
#endif

#ifndef MAXPATHLEN
# ifdef MAX_PATH
#  define MAXPATHLEN MAX_PATH
# elif defined(PATH_MAX)
#  define MAXPATHLEN PATH_MAX
# endif
#endif


#endif /*CMAKE*/

#ifdef XtSpecificationRelease
#include        <X11/Xmd.h>
#endif

#define MAXSTRINGLENGTH 1069
/* The first four-digit quasiall-even-digits non-quasi-repdigit emirp ! */

typedef char            STRING[MAXSTRINGLENGTH];
typedef char *          STRINGref;

/*
 *
 *      Define the FILESTATUS type that is returned by a sysdepend.c
 *      function.
 */

typedef int             FLAGS;

typedef struct  {
        FLAGS           fMode;
        int             iSize;
} FILESTATUSt;

#define FILEDOESNOTEXIST        0x00000000
#define FILENORMAL              0x00000001
#define FILEDIRECTORY           0x00000002


/*
 *---------------------------------------------------------------------
 *
 *      Define types that are CONSTANTLY used.
 *
 */


                        /* Initialize basics routines */

#define sBOOL(b)        ( b ? "TRUE" : "FALSE" )

#ifndef TRUE
# define        TRUE    1
#endif
#ifndef FALSE
# define        FALSE   0
#endif
#ifndef NULL
# define        NULL    0
#endif
#ifndef MAX
# define        MAX( x, y )     ( (x) < (y) ? (y) : (x) )
#endif
#ifndef MIN
# define        MIN( x, y )     ( (x) > (y) ? (y) : (x) )
#endif
#ifndef SWAP
# define        SWAP(a,b,temp)  {temp=a;a=b;b=temp;}
#endif
#ifndef SWAP_STRINGS
# define        SWAP_STRINGS(a,b,temp)  { \
        strcpy(temp,a);strcpy(a,b);strcpy(b,temp);}
#endif

/*
 *      Define some scientific programming constants and macros.
 */

#define MAXBONDLEN      3.0             /* 3.0 angstroms max covalent bond */
#define MINBONDLEN      0.5             /* 0.5 angstroms min covalent bond */

/* Use M_PI from math.h */
#ifndef PI
# define PI             3.1415926535897932384626433832795
#endif
/* 10^{-12} */
#define VERYSMALL       0.000000000001
#define TOLERANCE       0.001           /* Used for comparing doubles */
#define DEGTORAD        0.01745329251994329576 // Was 0.0174533


extern  double  zzdS,zzdS6;

#define pow2(x)         ( zzdS = x, zzdS*zzdS )
#define pow6(x)         ( zzdS6 = pow2(x), zzdS6*zzdS6*zzdS6 )




typedef void    GEN;
typedef GEN*    GENP;



/*
 *-------------------------------------------------------------------
 *
 *      Object verification
 *
 *      Object verification is performed by comparing the objects
 *      class structure pointer to a caller provided class structure.
 *      Failure of the pointers to match causes a panic.
 *
 *      To use this, add a GENP field to each object called (PKlass) and
 *      assign it to point to an arbitrary class structure at object
 *      creation time.  Then place CK(object,&KlassStructure) in every
 *      object function and macro that requires a valid object of the
 *      correct type.
 *
 *      I use the name KLASS rather than CLASS because other
 *      systems like X-windows use CLASS.
 */

typedef struct  KLASSs {
                char            *cPName;
                struct  KLASSs  *kSuperKlass;
                } KLASSt;

typedef KLASSt  *KLASS;


#ifdef  DEBUG

#define CK(x,c) ((!(x)||(GENP)((x)->PKlass)!=(GENP)(c))?BasicsKlassMismatchPanic((x),__FILE__,__LINE__):0)

#else

#define CK(x,c) 1

#endif




/*
 *-------------------------------------------------------------------
 *
 *      Interrupt routines are used to check if CONTROL-C has been
 *      hit.  The interrupt is set if an 'INT' signal is caught
 *      or in the graphical system, a CONTROL-C is received by
 *      the command window.
 *
 *      The way to use these routines is:
 *
 *              When checking for an interrupt:
 *                      1) Clear the interrupt using BasicsResetInterrupt
 *                      2) Keep checking the interrupt using bBasicsInterrupt
 *                      3) Clear the interrupt using BasicsResetInterrupt
 *
 *              When setting the interrupt:
 *                      1) Set the interrupt using BasicsInterrupt.
 */

                        /* bBasicsInterrupt() can be used to check */
                        /*      if CONTROL-C has been hit */
                        /* BasicsResetInterrupt clears the interrupt */
                        /* BasicsInterrupt will set the interrupt */

extern  BOOL    GbInterrupt;
#define bBasicsInterrupt()      ( GbInterrupt )
#define BasicsResetInterrupt()  ( GbInterrupt = FALSE )
#define BasicsSetInterrupt()    ( GbInterrupt = TRUE )


/*
 *      Some basic functions found in basics.c
 */


                /* Convert cases */
#define cUpper(c) (((c) >= 'a' && (c) <= 'z') ? ((c) & ~0x20) : (c))
#define cLower(c) (((c) >= 'A' && (c) <= 'Z') ? ((c) | 0x20) : (c))



/*
        General useful macros
*/


#ifndef isdigit
#define isdigit( c )    ( (c>='0')&&(c<='9'))
#endif

// for avoiding compiler warnings on compile-condition ignored variables (e.g. DEBUG)
#define UNUSED(x) ((void)(x))



/*
 *      Verbose printing.
 *
 *      This is a general facility for printing information
 *      to the user.
 *      The global variable GiVerbosityLevel contains a
 *      value that determines whether or not to print a message
 *      to the stdout.
 *      Also, if the file variable GhAllOut is not -1 then
 *      ALL output is written to the stream GhAllOut regardless
 *      of the verbosity level.
 *
 *      The effect is that VP0 messages are always printed, VP2 messages
 *      are only printed when the user requests maximum verbosity.
 *
 *      The code for myPrintf is system dependent and is found in
 *      'sysdepend.c'
 *
 *      The user can define many destinations for output.
 *      These are defined by calling iCreatePrintSink, which
 *      must be passed a pointer to a function that will print
 *      character strings, and a string that will be prefixed to
 *      the start of each output line.  iCreatePrintSink returns
 *      an integer that is used as a handle to redirect output
 *      by subsequent calls to DefineCurrentPrintSink, or
 *      to destroy the print sink using DestroyPrintSink.
 *
 *      By default there is only one print sink that directs output
 *      to the TTY screen using 'printf'.
 */

#define MAXCHARSPERPRINTF       4096    /* 4096 characters max */

typedef void            (*VFUNCTION)();

extern BOOL             GbPrintPrefix;
extern VFUNCTION        GfPrintStringCallback;
extern char             *GcPPrefix;
extern GENP             GPData;

extern int              GiTraceIndentationLevel;
extern int              GiVerbosityLevel;
extern FILE             *GfLog;

#define VerbositySet( x )        ( GiVerbosityLevel = x )
#define iVerbosity()             ( GiVerbosityLevel )
#define VerbositySetLogFile( x ) ( GfLog = x )
#define fVerbosityLogFile()      ( GfLog )

/*
 *      Verbose print macros.
 *
 *      Message level is the first argument; compared against GiVerbosityLevel
 *      inside myPrintf.  No double-parens, no global side effects.
 *
 *      VP(level, fmt, ...)  — generic, explicit level
 *      VP0..VP4             — convenience wrappers for common levels
 *      VPDISPLAY            — always prints (level -1, never filtered)
 *      VPLOG                — log-only / trace-level (level 99999)
 */
#define VPDISPLAY(...)      myPrintf(-1,     __VA_ARGS__)
#define VP0(...)            myPrintf(0,      __VA_ARGS__)
#define VP1(...)            myPrintf(1,      __VA_ARGS__)
#define VP2(...)            myPrintf(2,      __VA_ARGS__)
#define VP3(...)            myPrintf(3,      __VA_ARGS__)
#define VP4(...)            myPrintf(4,      __VA_ARGS__)
#define VPLOG(...)          myPrintf(99999,  __VA_ARGS__)
#define VP(level, ...)      myPrintf((level),__VA_ARGS__)


/*
 *      VPERROR, VPFATAL, VPFATALDELAYEDEXIT, VPFATALEXIT, VPWARN, and VPNOTE
 *      are a set of error handling macros based on the original VP* macros.
 *      They are the preferred mechanism for user notifications to ensure
 *      consistency in the user interface and from a software methodology
 *      perspective to follow the principle of information hiding.
 *
 *      Their behavior depends on whether LEaP is getting input from
 *      standard input (note that input to LEaP cannot be redirected)
 *      or from an input file.  The VPFATAL*EXIT macros check for stdin:
 *      they emit "Error" and exit if not stdin else they just emit and
 *      keep going.  The intention is to maintain LEaP's very forgiving
 *      behavior when used interactively, but to implement an immediate
 *      reporting and failure behavior for scripted use.
 */

extern STRING  GsProgramName;  /* in parser.h, but not easy to include that */
extern FILE    *fINPUTFILE(void);  /* in parser.h, but not easy to include that */
extern int     GiVPerrorCount;   /* Count the number of Fatal errors */
extern int     GiVPwarningCount; /* Count the number of Warnings */
extern int     GiVPnoteCount;    /* Count the number of Notes */

/*
 *      VPEPILOG emits the numbers of errors, warnings, and notes.
 */
#define VPEPILOG() \
    myPrintf(0, "\nExiting LEaP: Errors = %i; Warnings = %i; Notes = %i.\n", \
             GiVPerrorCount, GiVPwarningCount, GiVPnoteCount)

/*
 *      VPERROR is for the uncommon case of error reporting where immediate
 *      aborting is not desired, usually in order to check for more issues.
 */
#define VPERROR( fmt, ... ) do { \
        ++GiVPerrorCount; \
        if ( fINPUTFILE() == NULL ) { \
            myPrintf(0, "\nError: " fmt, ##__VA_ARGS__); \
        } else { \
            myPrintf(0, "\n%s: Error!\n" fmt, GsProgramName, ##__VA_ARGS__); \
        } \
    } while (0)

/*
 *      VPFATAL and VPFATALDELAYEDEXIT are for error handling that spans
 *      the function call graph.  A typical case is where initial error
 *      detection and reporting use VPFATAL and subsequent termination
 *      occurs in a parent function via VPFATALDELAYEDEXIT.
 */
#define VPFATAL( fmt, ... ) do { \
        ++GiVPerrorCount; \
        if ( fINPUTFILE() == NULL ) { \
            myPrintf(0, "\nError: " fmt, ##__VA_ARGS__); \
        } else { \
            myPrintf(0, "\n%s: Fatal Error!\n" fmt, GsProgramName, ##__VA_ARGS__); \
        } \
    } while (0)

#define VPFATALDELAYEDEXIT( fmt, ... ) do { \
        myPrintf(0, fmt, ##__VA_ARGS__); \
        if ( fINPUTFILE() != NULL ) { \
            myPrintf(0, "\nExiting LEaP: Errors = %i; Warnings = %i; Notes = %i.\n", \
                     GiVPerrorCount, GiVPwarningCount, GiVPnoteCount); \
            exit(21); \
        } \
    } while (0)

/*
 *      VPFATALEXIT is for the common case of error reporting and immediate
 *      aborting (in the non interactive mode of usage of LEaP) at the
 *      point of error detection.
 */
#define VPFATALEXIT( fmt, ... ) do { \
        ++GiVPerrorCount; \
        if ( fINPUTFILE() == NULL ) { \
            myPrintf(0, "\nError: " fmt, ##__VA_ARGS__); \
        } else { \
            myPrintf(0, "\n%s: Fatal Error!\n" fmt, GsProgramName, ##__VA_ARGS__); \
            myPrintf(0, "\nExiting LEaP: Errors = %i; Warnings = %i; Notes = %i.\n", \
                     GiVPerrorCount, GiVPwarningCount, GiVPnoteCount); \
            exit(31); \
        } \
    } while (0)

/*
 *      VPWARN does not abort.  VPWARN indicates that a message is
 *      a warning, and the intention is that users should investigate.
 */
#define VPWARN( fmt, ... ) do { \
        ++GiVPwarningCount; \
        if ( fINPUTFILE() == NULL ) { \
            myPrintf(0, "\nWarning: " fmt, ##__VA_ARGS__); \
        } else { \
            myPrintf(0, "\n%s: Warning!\n" fmt, GsProgramName, ##__VA_ARGS__); \
        } \
    } while (0)

/*
 *      VPNOTE does not abort.  VPNOTE indicates that a message is
 *      informational, and the intention is that users should note it.
 */
#define VPNOTE( fmt, ... ) do { \
        ++GiVPnoteCount; \
        if ( fINPUTFILE() == NULL ) { \
            myPrintf(0, "\nNote: " fmt, ##__VA_ARGS__); \
        } else { \
            myPrintf(0, "\n%s: Note.\n" fmt, GsProgramName, ##__VA_ARGS__); \
        } \
    } while (0)


#ifdef  DEBUG
/*
 *      Tracing macros based on the original VP* macros.
 *      But these are guarded with DEBUG.
 *
 *      myPrintTrace() builds "Trace: <prefix> <caller-msg> from call depth NN.\n"
 *      in a single call, avoiding the old sequence of 3 myPrintf calls that
 *      relied on GiVerbosity being set as a global side effect.
 */
#define TRACE_VERBOSITY 1069

extern void myPrintTrace(const char *prefix, int depth, const char *fmt, ...);

#define VPTRACEENTER( ... ) \
    myPrintTrace("Enter", GiTraceIndentationLevel++, __VA_ARGS__)
#define VPTRACEEXIT( ... ) \
    myPrintTrace("Exit ", --GiTraceIndentationLevel, __VA_ARGS__)
#define VPTRACEMULTIPLEEXIT( ... ) \
    myPrintTrace("Exit (premature due to multiple returns)", \
                 --GiTraceIndentationLevel, __VA_ARGS__)
#define VPTRACE( ... ) \
    myPrintf(TRACE_VERBOSITY, "Trace: " __VA_ARGS__)

#else
#define VPTRACEENTER( ... )
#define VPTRACEEXIT( ... )
#define VPTRACEMULTIPLEEXIT( ... )
#define VPTRACE( ... )
#endif  /* DEBUG */


                /* PRINTF is only to be used to dump messages */
                /* to stdout, which might be picked up by xMessageFilter */

#define PRINTF(fmt, ...)         printf( "+" fmt, ##__VA_ARGS__ )
#define PRINTF_NO_PREFIX(...)    printf( __VA_ARGS__ )

/*
 *-------------------------------------------------------------------
 *
 *      DEBUGGING messages.  It is better to use MESSAGE!!!!
 */
#ifdef  DEBUG
#define DDEBUG( ... )           printf( __VA_ARGS__ )
#else
#define DDEBUG( ... )
#endif



/*
 *--------------------------------------------------------------------
 *
 *      Messaging
 *
 *      By calling the macros MESSAGEON and MESSAGEOFF the
 *      caller can selectively turn on and turn off messages
 *      from different source files.
 *
 *      The flag DEBUG determines whether or not MESSAGE macros
 *      are even expanded.  If DEBUG is defined then MESSAGE macros
 *      are expanded, and whether or not the MESSAGE code is
 *      executed depends on bMessageCheck(__FILE__).
 */

                /* The following are characters prefixed to different */
                /* stdout messages */


#define MESSAGEON( s )          ( MessageAddFile(s) )
#define MESSAGEOFF( s )         ( MessageRemoveFile(s) )
#define MESSAGELISTFILES()      ( MessageFileList() )

#define PRINT_MESSAGE   '@'
#define PRINT_TRACE     '#'
#define PRINT_ALWAYS    '!'

#define MAXMESSAGEFILES 10

#ifdef  DEBUG
#define MESSAGE( ... ) do { \
        if (bMessageCheck(__FILE__)) { \
            printf( "%c%s %d|", PRINT_MESSAGE, __FILE__, __LINE__ ); \
            printf( __VA_ARGS__ ); \
            fflush(stdout); \
        } \
    } while(0)

#define MESSAGEEXECUTE( c )     if (bMessageCheck(__FILE__)) c;

#define FILEMESSAGE(f, ... ) do { \
        if ( bMessageCheck(f) ) { \
            printf( "%c%s %d|", PRINT_MESSAGE, __FILE__, __LINE__ ); \
            printf( __VA_ARGS__ ); \
            fflush(stdout); \
        } \
    } while(0)

#else
#define MESSAGE( ... )
#define MESSAGEEXECUTE( c )
#define FILEMESSAGE(f, ... )
#endif

#define ALWAYS( ... ) do { \
        printf( "%c%s %d|", PRINT_MESSAGE, __FILE__, __LINE__ ); \
        printf( __VA_ARGS__ ); \
        fflush(stdout); \
    } while(0)



/*
 *-----------------------------------------------------------------------
 *
 *      File handling function wrappers.
 *
 */



#define FOPENCOMPLAIN(a,b)      fBasicsMyFopen(a,b,TRUE)
#define FOPENNOCOMPLAIN(a,b)    fBasicsMyFopen(a,b,FALSE)

extern STRING   GsBasicsFullName;

#define fsBasicsFileStatus(fn)  fsSysdependFileStatus(fn)
#define BasicsCurrentWorkingDirectory(path) \
                SysdependCurrentWorkingDirectory(path)


/*
 *-----------------------------------------------------------------------
 *
 *      Fail safe memory allocation
 *
 *      Only use the first set of memory allocation tools when
 *      debugging MEMORY PROBLEMS!!!  They are very slow because
 *      they constantly check malloc'd memory for overruns and
 *      other errors.
 */

/*
 *      MEMORY_DEBUG can be 0, 1, 2, 3, 4
 *
 *              0 -     Dont do any MEMORY debugging.
 *              1 -     Use my memory usage tracking.
 *                              This adds 4 bytes to every MALLOC and stores
 *                              the size of the allocation in an int
 *                              at the front of the memory block.
 *              2 -     Use my memory debugging.
 *              3 -     Use SUN memory debugging, requires linking with
 *                              /usr/lib/debug/malloc.o
 *              4 -     Use my memory debugging, but test memory every time
 *                              a routine is entered or exited.
 *
 *      NOTE: Whenever it is changed, EVERYTHING must be recompiled!!!!!!!
 */

#define MEM_NO_DEBUG                    0
#define MEM_TRACKING                    1
#define MEM_DEBUG                       2
#define MEM_SUN_DEBUG                   3
#define MEM_DEBUG_ULTRA_PARANOID        4


                                /* By default don't use memory debugging */
#ifndef MEMORY_DEBUG
# define MEMORY_DEBUG     MEM_NO_DEBUG
#endif


extern int      GiMemoryAllocated;

extern STRING   GsLastMemOpFile;
extern int      GiLastMemOpLine;
extern BOOL     GbTestMemory;


                        /* Turn my memory testing on/off */

#define TEST_MEMORY_ON(b)       ( GbTestMemory = b )
#define bTEST_MEMORY_ON()       ( GbTestMemory )

#define REGISTERMEMOP() ( strcpy( GsLastMemOpFile, __FILE__ ),\
                          GiLastMemOpLine = __LINE__ )

#if MEMORY_DEBUG == MEM_NO_DEBUG
# define        INITMEMORYDEBUG()
# define        TESTMEMORY()
# define        MALLOC( dest, type, size ) { \
    (dest) = (type)malloc(size); \
    if ( (dest)==(type)NULL ) { \
        DFATAL( "Malloc: %s", strerror(errno) ); \
        }}

# define        FREE(x) { \
                        *(char*)(x)='?';free(x); x=NULL; }

# define        REALLOC( dest, type, src, size ) { \
        (dest) = (type) realloc( src, size ); \
        if ( (dest)==(type)NULL ) { \
                DFATAL( "Realloc: %s", strerror(errno) ); \
        }}

# define        LISTUNFREEDMEMORYTOLOGFILE()

#endif


#if MEMORY_DEBUG == MEM_TRACKING
# define        TRUEPOS(x)      ((int*)(x)-1)
# define        INITMEMORYDEBUG()
# define        TESTMEMORY()
# define        MALLOC( dest, type, size ) { \
    (dest) = (type)((char*)malloc(size+sizeof(int))+sizeof(int)); \
    FILEMESSAGE( "MEMORY", ("Malloc %d bytes at: %lX\n", size, (dest) ) ); \
    *(int*)(TRUEPOS(dest)) = size; \
    GiMemoryAllocated += size; \
    if ( (dest)==(type)NULL ) { \
        printf( "BAD MALLOC file: %s line: %d\n", __FILE__, __LINE__); \
        exit(1);}}

# define        FREE(x) {  \
        *(char*)(x)='?'; \
        FILEMESSAGE( "MEMORY", ("Free %d bytes at: %lX\n", *(int*)TRUEPOS(x), x ) ); \
        GiMemoryAllocated -= *(int*)TRUEPOS(x); \
        free(TRUEPOS(x)); \
        x=NULL; }

# define        REALLOC( dest, type, src, size ) { \
        FILEMESSAGE( "MEMORY", ("<<<Realloc %d bytes at: %lX\n", size, src ) ); \
        GiMemoryAllocated -= *TRUEPOS(src); \
        (dest) = (type)((int*)realloc( TRUEPOS(src), size+sizeof(int) )+1); \
        FILEMESSAGE( "MEMORY", \
                        (">>>Realloc %d bytes at: %lX\n", size, (dest) ) ); \
        GiMemoryAllocated += size; \
        *(int*)(TRUEPOS(dest)) = size; \
}

# define        LISTUNFREEDMEMORYTOLOGFILE()

#endif

#if (MEMORY_DEBUG==MEM_DEBUG)||(MEMORY_DEBUG==MEM_DEBUG_ULTRA_PARANOID)

# define        INITMEMORYDEBUG()

#define TESTMEMORY()    { \
        DebugMemoryTest( __FILE__, __LINE__, FALSE ); \
        REGISTERMEMOP(void); \
}

#define MALLOC( dest, type, size ) { \
    (dest) = (type)DebugMalloc((long)(size), __FILE__, __LINE__ ); \
    if ( (dest)==(type)NULL ) { \
        DFATAL( "Bad malloc" ); \
    } REGISTERMEMOP(void); \
}

#define FREE(x)        { \
        *(char*)(x)='?'; \
        DebugFree(x,__FILE__,__LINE__);x=NULL; \
        REGISTERMEMOP(void); \
}

#define REALLOC( dest, type, src, size ) { \
        (dest)=(type)DebugRealloc((src),(long)(size),__FILE__,__LINE__); \
                REGISTERMEMOP(void); \
}

#define LISTUNFREEDMEMORYTOLOGFILE()    { \
        DebugMemoryTest( __FILE__,__LINE__, TRUE ); \
        REGISTERMEMOP(void); \
}
#endif

#if (MEMORY_DEBUG == MEM_DEBUG_ULTRA_PARANOID)
# define        MEMORY_TEST_MACRO()     TESTMEMORY()
#else
# define        MEMORY_TEST_MACRO()
#endif


#if MEMORY_DEBUG == MEM_SUN_DEBUG

# define        INITMEMORYDEBUG() malloc_debug(2);
# define        TESTMEMORY()    (malloc_verify(),REGISTERMEMOP())

# define        MALLOC( dest, type, size ) { \
    (dest) = (type)malloc(size); \
    if ( (dest)==(type)NULL ) { \
        printf( "BAD MALLOC file: %s line: %d\n", __FILE__, __LINE__); \
        exit(1);}REGISTERMEMOP(void); \
}

# define        FREE(x)    { \
                        *(char*)(x)='?',free(x); \
                        x=NULL;REGISTERMEMOP(void); \
}

# define        REALLOC( dest, type, src, size ) { \
                        (dest = (type) realloc( src, size ),REGISTERMEMOP(void);) }
# define        LISTUNFREEDMEMORYTOLOGFILE()
#endif


#define DFATAL(fmt, ...) do { \
    printf( "%cFATAL ERROR----------------------------------------\n", PRINT_ALWAYS ); \
    printf( "%cFATAL:    In file [%s], line %d\n", PRINT_ALWAYS, __FILE__, __LINE__ ); \
    printf( "%cFATAL:    Message: " fmt "\n", PRINT_ALWAYS, ##__VA_ARGS__ ); \
    printf( "%cABORTING.\n", PRINT_ALWAYS ); \
    exit(1); \
} while(0)



/*
#include        "objekt.h"
#include        "collection.h"
#include        "dictionary.h"
#include        "parmSet.h"
#include        "list.h"
#include        "assoc.h"
#include        "vector.h"
#include        "container.h"
#include        "atom.h"
#include        "bag.h"
#include        "block.h"
#include        "loop.h"
#include        "unit.h"
#include        "minimizer.h"
#include        "byteArray.h"
*/
extern void            IMem(void);
extern void            TMem(void);

/*  sysdepend.c  */

extern void             myPrintf(int iVerbosity, const char *fmt, ...);
extern void             myPrintString(const char *cPString, int iVerbosity);
extern void             myPrintTrace(const char *prefix, int depth,
                                     const char *fmt, ...);
extern void             SysdependDirectoryList( char *cPPath,
                                STRING *saPNames[], int *iPNumber);
extern FILESTATUSt      fsSysdependFileStatus(char *cPName);
extern void             SysdependCurrentWorkingDirectory(STRING sPath);

/*  basics.c  */

extern double           myAcos( double d );
extern double           myPow( double x, double y );
extern int              iDoubleCompare( double dA, double dB );
extern BOOL             bStringToDouble( char *cPData, double *dPData );
extern BOOL             bStringToInt( char *cPData, int *iPData );
extern void             StringLower( char *sStr );
extern void             StringUpper( char *sStr );
extern void             StringTrim( char *sStr );
extern void             StringRTrim( char *sStr );
extern char             *DebugMalloc( long lSize, char *sFile, int iLine );
extern char             *DebugRealloc( char *cPBlock, long lSize, char *sFile,
                                int iLine );
extern void             DebugFree( char *cPBlock, char *sFile, int iLine );
extern BOOL             bMessageCheck( char *sFile );
extern void             MessageAddFile( char *sFile );
extern void             MessageRemoveFile( char *sFile );
extern void             MessageFileList(void);
extern int              BasicsAddDirectory( STRING sDirectory, int bomb );
extern FILE             *fBasicsMyFopen( char *sFilename, char *sAttributes,
                                BOOL bComplain );
extern int              iCreatePrintSink( VFUNCTION fOutputCallback,
                                char *sPrefix, GENP PData );
extern void             DestroyPrintSink( int iHandle );
extern void             PushCurrentPrintSink( int iSink );
extern void             PopCurrentPrintSink(void);
extern void             BasicsInitialize(void);
extern void             BasicsKlassMismatchPanic( GENP PObj, char *cPFile,
                                int iLine );


#endif  /* #ifndef BASICS_H */
