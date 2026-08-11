/*
 *      File:   parser.y
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
 *      A YACC program for parsing commands to LEaP.
 *
 *      The syntax is described in this file.
 *
 *      The syntax has the form:        [] delimit necessary syntax elements
 *
 *      [command] [arg1], [arg2], [arg3], ... ;
 *      [variable] = [expression];
 *      [variable] = ;                           Clears a variable 
 *      
 *
 *      The parser requires all commands to be ended with a semicolon ';' or newline
 *      Comments are started with a '#' and terminated with a '\n'.
 *      Strings can be delimited with double quotes '"', or
 *      started with a single '$' character and ended with a comma, space,
 *      curly bracket '}', semicolon, equals sign, etc.
 *
        
 */


%token  LVARIABLE
%token  LSTRING
%token  LNUMBER
%token  LASSIGN 
%token  LENDOFCOMMAND
%token  LOPENLIST
%token  LASSOC
%token  LCLOSELIST
%token  LOPENPAREN
%token  LCLOSEPAREN
%token  LQUIT
%token  LCOMMA
%token  LDOT
%token  LCOMMAND
%token  LDUMMY
%token  LNULL
%token  LNOTSINGLECHAR

%token  LEVAL
%token  LINT LFRAC LSTR
%token  LPLUS LMINUS LMUL LDIV LPOW
%left LPLUS LMINUS
%left LMUL LDIV
%right UMINUS LPOW
%type <aVal> rawexp scalar evalcmd evalexpr evalterm evalfactor

%{
#include	<unistd.h>
#include	<math.h>
#include        "basics.h"

#include        "classes.h"

#include        "dictionary.h"
#include        "parmLib.h"

#include        "commands.h"
#include	"block.h"
#include	"parser.h"

#include        "leap.h"
#include        "block.h"
#include        "select_mask.h"

#include        "help.h"
#include        "defaults.h"

#define         MESSAGEFILTER   MESSPARSER




#define         NULLSTR         "null"

/*
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

        GLOBAL VARIABLES

*/

                /* The global variable GfCurrentInput defines the file */
                /* from which input is currently read.  When the file */
                /* is empty, then switch back to stdin */

#define	MAXINPUT	1000
#define	MAXINPUTFILES	10		/* Maximum 10 input files */
					/* can be open at once */


char            GsInputLine[MAXINPUT] = "";
bool		GbLastLine = false;
int             GiInputPos = 0;
int             GiInputStartPos = 0;
PARMLIB		GplAllParameters;
RESULTt		GrMainResult;
BLOCK		GbCommand = NULL;
BLOCK		GbExecute = NULL;
int		GiClipPrompts = 0;
bool		GbGraphicalEnvironment;
STRING		GsProgramName;

extern int	iMemDebug;

static	STRING	*SbFirstSourceFiles = NULL;
static  bool    SbLastTokenWasCommand = false;
static	int	iFirstSource = 0;
static	bool	SbUseStartup = true;



/*
 *-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
 *
 *	The following is used by the parser to maintain a stack of
 *	files where input is received from.  If the file is NULL
 *	then input is received from the main program in the form
 *	of BLOCKS.  The main program is then responsible for reading
 *	the stdin ( in the command line interface ) or for gathering
 *	keypress events from X-Windows ( in the graphical interface ).
 *
 */


int             GiInputFileStackPos = 0;
FILE*           GfaInputFileStack[MAXINPUTFILES];



/*
 *----------------------------------------------------------------
 *
 *	Not quite GLOBAL variables used by the parser.
 */

                                /* Arguments to routines are passed through */
                                /* an array */
#define MAXARGS         10
#define MAXLISTNEXT     10


ATOM            aDummy;
ASSOC           aaArgs[MAXARGS];
int             iArgCount;

                                /* List stuff is used for input of nested */
                                /* lists */
#define MAXLISTNEST     10
ASSOC           aaLists[MAXLISTNEST];
int             iCurrentList = -1;
#define PUSHLIST()      iCurrentList++
#define POPLIST()       iCurrentList--
#define CURRENTLIST     aaLists[iCurrentList]


OBJEKT          o0;
double          dTemp;
ASSOC           aAssoc;
STRING          sTemp;
bool		bQuit = false;
bool		bCommandFound = false;

                /* There seems to be a problem with YACC not properly */
                /* declaring yylval and yyval */
typedef struct  {
	ASSOC		aVal;
	double		dVal;
	STRING		sVal;
	CMD_FUNCTION	fCallback;
} YYSTYPEt;

#define YYSTYPE YYSTYPEt


static char * zsScanMaskLiteral(void);
static  OBJEKT zoGetObject( char *sName );
static  int    yyerror( const char *sStr );
static  int    yylex();
static  int    yyparse();
static ASSOC zaEvalBinary(char operator, ASSOC arg1, ASSOC arg2 );
%}


%start  input





%%
/*------------------------------------------------------------

        RULES
*/

/*
 *	Bogus function for xaUtilMessageFilter
 *
void
yyparse()
{
*/

input   :       line
		{
			return 0;
		}
		;

line    :       LENDOFCOMMAND
	|	instruct
                        {
                        bCommandFound = false;
                        }
        |       error LENDOFCOMMAND
                        {
                            VP0("\n" );
                            yyerrok;
                        }
        ;

instruct:       assign LENDOFCOMMAND
        |       function LENDOFCOMMAND { DEREF($<aVal>1); }
        ;

assign  :       LVARIABLE LASSIGN express 
                        {
                                /* Set the value of the variable */
                            aAssoc = $<aVal>3;
			    if ( aAssoc != NULL ) {
                                MESSAGE("Assigning a value to %s\n", 
                                                        $<sVal>1 );
                                VariableSet( $<sVal>1, oAssocObject(aAssoc) );
				MESSAGE("DEREF (assign) - %s\n",
							sAssocName(aAssoc) );
			    } else {
				MESSAGE("Not assigning value to %s - rmving\n",
							$<sVal>1 );
				VariableRemove( $<sVal>1 );
			    }
                            DEREF( aAssoc );

                        }
        |       LVARIABLE LASSIGN 
                        {
                            MESSAGE("Removing variable %s\n", $<sVal>1 );
                            VariableRemove( $<sVal>1 );
                        }
        ;

express :       rawexp
        |       function
        |       evalcmd
        ;

rawexp  :       LOPENLIST 
                        {
                            PUSHLIST();
                                /* Create an ASSOC for the list */
                            aAssoc = (ASSOC)oCreate(ASSOCid);
                            AssocSetName( aAssoc, "" );
                            OBJEKT oObject = oCreate(LISTid);
                            LIST_from(oObject)->bFreeChildren = true;
                            AssocTakeObject( aAssoc, oObject );
                            CURRENTLIST = aAssoc;
                        }
                    elements LCLOSELIST 
                        {
                            $<aVal>$ = $<aVal>3;
                            POPLIST();
                        }
        |       LASSOC 
                        {
                            $<aVal>$ = $<aVal>1;
                        }
        |       scalar
        |       LDUMMY     /*  '*' wildcard   */
                        {
                            aAssoc = (ASSOC)oCreate(ASSOCid);
                            AssocSetName( aAssoc, "" );
                            AssocSetObject( aAssoc, aDummy );
                            $<aVal>$ = aAssoc;
                        }
        |       LNULL
                        {
                            MESSAGE("Parsed a null\n" );
                            aAssoc = (ASSOC)oCreate(ASSOCid);
                            AssocSetName( aAssoc, "" );
                            AssocTakeObject( aAssoc, NULL );
                            $<aVal>$ = aAssoc;
                        }
        ;

scalar :       LNUMBER
                        {
                            aAssoc = (ASSOC)oCreate(ASSOCid);
                            AssocSetName( aAssoc, "" );
                            o0 = oCreate(ODOUBLEid);
                            ODoubleSet( o0, $<dVal>1 );
                            AssocTakeObject( aAssoc, o0 );
                            $<aVal>$ = aAssoc;
                        }
        |       LSTRING
                        {
                            aAssoc = (ASSOC)oCreate(ASSOCid);
                            AssocSetName( aAssoc, "" );
                            o0 = oCreate(OSTRINGid);
                            OStringDefine( (OSTRING) o0, $<sVal>1 );
                            AssocTakeObject( aAssoc, o0 );
                            $<aVal>$ = aAssoc;
                        }
        |       LVARIABLE
                        {
			    OBJEKT o = zoGetObject($<sVal>1); /* returns with coutned REF */
                            aAssoc = (ASSOC)oCreate(ASSOCid);
                            AssocSetName( aAssoc, $<sVal>1 );
                            AssocTakeObject( aAssoc, o );
                            $<aVal>$ = aAssoc;
                        }
         ;

function:       cmdname 
                        { iArgCount = 0;
                        } 
                    args  
                        {
                            struct timespec tsStart, tsEnd;
                            clock_t ctStart, ctEnd;
                            clock_gettime(CLOCK_MONOTONIC, &tsStart);
                            ctStart = clock();
                            /* Execute the command */
			    MESSAGE("executing function\n");
                            /* Command execute;returns with counted REF */
                            o0 = $<fCallback>1( iArgCount, aaArgs );
			    if ( o0 != NULL ) {
                            	aAssoc = (ASSOC)oCreate(ASSOCid);
                            	AssocTakeObject( aAssoc, o0 );
                            	$<aVal>$ = aAssoc;
			    } else {
				MESSAGE("func == NULL---\n");
				$<aVal>$ = NULL;
			    }
                             /* DEREF each of the arguments */
                            for (int i=0; i<iArgCount; i++ ) {
				MESSAGE("DEREF (function) - %s\n", sAssocName(aaArgs[i]));
                                DEREF( aaArgs[i] );
                            }
                            if(GDefaults.bTiming ){
                                clock_gettime(CLOCK_MONOTONIC, &tsEnd);
                                ctEnd = clock();
                                double dt = (double)(tsEnd.tv_sec - tsStart.tv_sec)
                                          + (double)(tsEnd.tv_nsec - tsStart.tv_nsec) / 1e9;
                                double dCPU = ((double) (ctEnd - ctStart)) / (double)CLOCKS_PER_SEC;
                                VP0("Command: %12s Elapsed time: %8g sec   CPU time: %8g sec\n",$<sVal>1,dt,dCPU);
                            }
                        }
        ;

evalcmd :       LEVAL evalexpr
                        {
                            $<aVal>$ = $<aVal>2;
                            if ($<aVal>$ != NULL && iVerbosity() > 2 ) {
                                OBJEKT oVar = oAssocObject($<aVal>$);
                                if (iObjectType(oVar) == OSTRINGid)
                                    VP2("Expression evaluated to string: \"%s\"\n", sOString(oVar));
                                else
                                    VP2("Expression evaluated to number: %g\n", dODouble(oVar));
                            }
                        }
        ;

evalexpr: evalexpr LPLUS evalterm
                        {
                            $<aVal>$ = zaEvalBinary('+', $<aVal>1, $<aVal>3);
                            DEREF($<aVal>1);
                            DEREF($<aVal>3);
                        }
                        | evalexpr LMINUS evalterm
                        {
                            $<aVal>$ = zaEvalBinary('-', $<aVal>1, $<aVal>3);
                            DEREF($<aVal>1);
                            DEREF($<aVal>3);
                        }
                        | evalterm
                        {
                            $<aVal>$ = $<aVal>1;
                        }
                        ;

evalterm : evalterm LDUMMY evalpower  /* LDUMMY == '*' */
                        {
                            $<aVal>$ = zaEvalBinary('*', $<aVal>1, $<aVal>3);
                            DEREF($<aVal>1);
                            DEREF($<aVal>3);
                        }
                        | evalterm LDIV evalpower
                        {
                            $<aVal>$ = zaEvalBinary('/', $<aVal>1, $<aVal>3);
                            DEREF($<aVal>1);
                            DEREF($<aVal>3);
                        }
                        | evalpower
                        {
                            $<aVal>$ = $<aVal>1;
                        }
         ;

evalpower: evalfactor LPOW evalpower
                        {
                            $<aVal>$ = zaEvalBinary('^', $<aVal>1, $<aVal>3);
                            DEREF($<aVal>1);
                            DEREF($<aVal>3);
                        }
                        | evalfactor
                        {
                            $<aVal>$ = $<aVal>1;
                        }
          ;

evalfactor: scalar
                        | LMINUS evalfactor
                        {
                            // This is the only unary operator
                            $<aVal>$ = $<aVal>2; // transfer ref
                            dODouble(oAssocObject($<aVal>$)) *= -1.0;
                        }
                        | LINT LOPENPAREN evalexpr LCLOSEPAREN
                        {
                            OBJEKT doVal = oAssocObject($<aVal>3);
                            double int_part;
                            modf(dODouble(doVal), &int_part); 
                            ODoubleSet(doVal, int_part);
                            $<aVal>$ = $<aVal>3; //transfer ref
                        }
                        | LOPENPAREN evalexpr LCLOSEPAREN
                        {
                            $<aVal>$ = $<aVal>2; //transfer ref
                        }
          ;

elements:       /* nothing */
        |       elements rawexp 
                        {
                                /* Get the element and add it to the list */
                            MESSAGE("Adding to list:\n" );
                            ListAddToEnd( (LIST)oAssocObject(CURRENTLIST), 
							(OBJEKT) $<aVal>2 ); // adds REF
			    DEREF( $<aVal>2 );
                            $<aVal>$ = CURRENTLIST;
                        } 
        ;

cmdname :	LCOMMAND
        ;

args    :       /* Nothing */
        |       arg
        |       args arg
        ;

arg     :       rawexp  
                        {
                            aaArgs[iArgCount++] = $<aVal>1;
                        }
        ;



/*
 *
 *	Bogus stuff for xaUtilMessageFilter
 *
}
 */



%%
/*------------------------------------------------------------

        ROUTINES

*/

static  bool    SbGotUngetc = false;
static  char    ScUngetc;


/*
 *      yyerror
 *
 *      Respond to errors.
 */
int
yyerror( const char *sStr )
{
    VP0("\n%s%*s%.*s",GsInputLine,GiInputStartPos-1,"",
              GiInputPos-GiInputStartPos+1,
              "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^...");
    VPFATALEXIT("\nError from the parser: \n       %s.\n"
            "       Check for typos, misspellings, etc.\n"
            "       Hint: previous token was%s a command\n"
            "       Try help on the command name and desc on the command arguments.\n",
            sStr, SbLastTokenWasCommand? "":" not" );
    return 1;
}

FILE *
fINPUTFILE()            
{
	if ( GiInputFileStackPos < 0 )
		return(NULL);
	return( GfaInputFileStack[GiInputFileStackPos] );
}

/*
 *	zbGetLine
 *
 *	Get the next line of input from the current input source.
 *	This may be GbCommand or GbExecute depending if
 *	the input is coming from the user or from an execute file.
 *	Return FALSE if there is no more lines to be received from 
 *	the BLOCK.
 */
bool
zbGetLine( char *sLine, bool *bPFromExecute )
{
bool		bGotBlock;
char		c;

    /* VPTRACEENTER("zbGetLine" ); */

    if ( fINPUTFILE() == NULL ) {
	*bPFromExecute = false;
	/* VPTRACEEXIT("zbGetLine" ); */
	return(bBlockReadLine( GbCommand, sLine ));
    }

    *bPFromExecute = true;

		/* If there is no line in the execute BLOCK */
		/* then read in another block from the current file */

    if ( bBlockEndOfRead(GbExecute) ) {
	    BlockEmpty( GbExecute );
	    bGotBlock = false;
	    while ( !feof(fINPUTFILE()) && !bGotBlock ) {
		c = fgetc(fINPUTFILE());
		if ( feof(fINPUTFILE()) ) break;
		bGotBlock = bBlockAddChar( GbExecute, c );
	    }

		/* If a complete BLOCK was not read then */
		/* append a final '\n' character, if that doesn't */
		/* make this a complete BLOCK then there is an error */
		/* in the input, also we are at the end of the file */
		/* so pop the file off the stack and say that we have */
		/* a complete BLOCK */

	    if ( !bGotBlock ) {
	    	bGotBlock = bBlockAddChar( GbExecute, '\n' );
	    	if ( fINPUTFILE() != NULL )
	    		fclose( fINPUTFILE() );
	    	INPUTPOPFILE();
	    	VPTRACE("After pop in zbGetLine GiInputFileStackPos = %d\n",
	    	          GiInputFileStackPos );
	    }
    }
    /* VPTRACEEXIT("zbGetLine" ); */
    return(bBlockReadLine( GbExecute, sLine ));
}

void zStripComment(char *s)
{
    bool bInQuote = false;
    for (; *s; s++) {
        if (*s == '"')
            bInQuote = !bInQuote;
        if (*s == '#' && !bInQuote) {
            s[0] = '\n';
            s[1] = 0;
            break;
        }
    }
}


/*
 *      cGetChar
 *
 *      Get the next character from the line buffer.
 *	If there are no more characters in the line buffer and
 *	GbLastLine is TRUE then return '\0', otherwise if the
 *	line buffer is empty, fill it and return the next character
 *	in it.
 */
char
cGetChar(void)
{
char            c;
bool		bFromExecute;

                /* If there is a pushed character then return it */

    if ( SbGotUngetc ) {
        SbGotUngetc = false;
        c = ScUngetc;
        goto DONE;
    }

		/* Now if the input line is empty then fill it */
    if ( GsInputLine[GiInputPos] == '\0' ) {
	if ( GbLastLine ) {
	    c = '\0';
	    GbLastLine = false;
	    goto DONE;
	}
	GbLastLine = zbGetLine( GsInputLine, &bFromExecute );
	GiInputPos = 0;
	if ( bFromExecute ) {
	    for (int i=GiClipPrompts; i<=iINPUTSTACKDEPTH(); i++ ) VP2(">" );
	    VP2(" " );
	    VP2("%s", GsInputLine );
	} else {
	    VPLOG("> %s", GsInputLine );
	}
        zStripComment(GsInputLine);
    }

    c = GsInputLine[GiInputPos++];
    
DONE:
    return(c);

}



/*
 *      zUngetc
 *
 *      Push one character back to the file.
 */
void
zUngetc( char c )
{
    SbGotUngetc = true;
    ScUngetc = c;
}



int
iOneCharToken( int c )
{
	switch(c) {
		case ';':
        		return(LENDOFCOMMAND);
		case '=':
        		return(LASSIGN);
		case '*':
        		return(LDUMMY);
		case '(':
        		return(LOPENPAREN);
		case ')':
        		return(LCLOSEPAREN);
		case '{':
        		return(LOPENLIST);
		case '}':
        		return(LCLOSELIST);
		default:
        		return(LNOTSINGLECHAR);
	}
}

/////////////////////////////////////////
//#include <string.h>
//#include <stdio.h>

/*
 * Assumed external declarations (adjust to match your actual headers):
 *
 *   typedef struct { ... } OBJEKT;
 *
 *   OBJEKT *oVariable(const char *sVarName);   // returns NULL if not found
 *   char   *sOString(OBJEKT *oVar);            // returns string form of oVar
 */

#define MAX_VARNAME_LEN 128

/*
 * ReplaceVarPattern
 * -----------------
 * Scans sInput for occurrences of "{name}" style placeholders
 * (a '[', an identifier, a ']'). For each one found:
 *
 *   - Extracts the identifier text between the braces -> sVarName
 *   - Calls oVariable(sVarName) to look up the object
 *   - If the lookup is NOT NULL, replaces "{sVarName}" in the output
 *     with the string returned by sOString(oVar)
 *   - If the lookup IS NULL, the placeholder is left untouched
 *     (including its braces) in the output
 *
 * Parameters:
 *   sInput      - the source string to scan (unmodified)
 *   sOutput     - caller-supplied buffer to receive the result
 *   nOutputSize - size of sOutput buffer, in bytes (including NUL)
 *
 * Returns:
 *   0 on success
 *  -1 if the output buffer was too small (result is truncated but
 *     still NUL-terminated)
 */
int ReplaceVarPattern(const char *sInput, char *sOutput, size_t nOutputSize)
{
    size_t nOutPos = 0;
    const char *p = sInput;
    int bTruncated = 0;

    if (!sInput || !sOutput || nOutputSize == 0)
        return -1;

    #define APPEND_CHAR(c)                                  \
        do {                                                \
            if (nOutPos + 1 < nOutputSize) {                \
                sOutput[nOutPos++] = (c);                    \
            } else {                                         \
                bTruncated = 1;                              \
            }                                                \
        } while (0)

    #define APPEND_STR(s)                                   \
        do {                                                \
            const char *_s = (s);                            \
            while (*_s) {                                    \
                APPEND_CHAR(*_s);                             \
                _s++;                                          \
            }                                                \
        } while (0)
    const char open='[', close=']';

    while (*p != '\0')
    {
        if (*p == open)
        {
            const char *pOpen  = p;
            const char *pScan  = p + 1;
            char sVarName[MAX_VARNAME_LEN];
            size_t nLen = 0;
            int bFoundClose = 0;

            /* Try to read an identifier up to a matching '}' */
            while (*pScan != '\0' && *pScan != open  )
            {
                if (*pScan == close)
                {
                    bFoundClose = 1;
                    break;
                }
                if (nLen + 1 < sizeof(sVarName))
                {
                    sVarName[nLen++] = *pScan;
                }
                else
                {
                    /* identifier too long to be a valid placeholder;
                       bail out and treat '{' as a literal char */
                    break;
                }
                pScan++;
            }
            sVarName[nLen] = '\0';

            if (bFoundClose && nLen > 0)
            {
                OBJEKT oVar = oVariable(sVarName);

                if (oVar != NULL && iObjectType(oVar) == OSTRINGid)
                {
                    APPEND_STR(sOString(oVar));
                    /* advance past the closing '}' */
                    p = pScan + 1;
                    continue;
                }
                else
                {
                    /* oVariable() returned NULL: leave "{name}" as-is */
                    size_t nSpan = (size_t)(pScan - pOpen) + 1; /* includes braces */
                    size_t i;
                    for (i = 0; i < nSpan; i++)
                    {
                        APPEND_CHAR(pOpen[i]);
                    }
                    p = pScan + 1;
                    continue;
                }
            }
            /* Not a valid {identifier} pattern; copy '{' literally */
            APPEND_CHAR(*p);
            p++;
        }
        else
        {
            APPEND_CHAR(*p);
            p++;
        }
    }

    sOutput[nOutPos < nOutputSize ? nOutPos : nOutputSize - 1] = '\0';

    #undef APPEND_CHAR
    #undef APPEND_STR

    return bTruncated ? -1 : 0;
}


/*
 *      yylex
 *
 *      Lexical analyzer for the parser.
 *      Read characters from stdin and return the token types
 *      read and place the value read into the global UNION
 *      yylval.
 *
 *      The things that it recognizes are:
 *              LDOUBLE         [-]###.###E## or ###
 *              LSTRING         "xx xx xxx" or '$'everything up to ' ' ',' ';'
 *              commands        xxxxxxx
 *              LVARIABLE       xxxxxxx which are not commands
 *              LTERMINATOR     ;
 *              LASSIGN         =
 *
 *	Modified 17 November 1992 - David A. Rivkin
 *		Added checking the alias table for command matches.
 *	Total rewrite October 1993 - Bill Ross
 *
 */
int
yylex(void)
{
STRING          sStr;
int             iMax, tok;
bool            bGotExp, bGotDot;
char            c;
STRING		sPossibleCmd;

    SbLastTokenWasCommand = false;
    /*
     * Skip whitespace: blanks, tabs, end of lines, carriage returns, and commas.
     */
    while ( (c=cGetChar())==' ' || c=='\t' || c=='\n' || c=='\r' || c == ',' );

    if ( c == '\0' ) 
	return(LENDOFCOMMAND);

    GiInputStartPos = GiInputPos;
    /*
     *  Check the 1-character possibilities: , ; = * ( ) { }
     */
    tok = iOneCharToken( c );
    if ( tok != LNOTSINGLECHAR ) {
        MESSAGE("Parsed /%c/\n", c );
	return(tok);
    }

    /*
     *  it isn't a 1-char thing; read in the rest 
     *	and push back the terminating char
     */
    sStr[0] = c;
    for (int j=1;;j++) {

	if ( j >= sizeof(STRING) )
	    DFATAL("string too long");

	c = cGetChar();
	/*
	 *  NULL terminates anything (?)
	 */
	if ( c == '\0' ) {
	    	sStr[j] = '\0';
	    	break;
	}

        // Backslash gets \ plus next char, regardless of anything else
        if ( c == '\\' ) {
	    sStr[j++] = '\\';
            c = cGetChar();
            sStr[j] = c;
            if ( c == '\0' ) break;
            continue;
        }
             
	/*
	 *  allow anything inside quotes; chop closing quote
	 */
	if ( sStr[0] == '"' ) {
		if ( c == '"' ) {
			sStr[j] = '\0';
			break;
		}
		sStr[j] = c;
		continue;
	}
	/*
	 *  whitespace is a delimiter outside of quotes
	 */
	if ( c == ' ' || c == '\t' || c == '\n' || c=='\r' || c == ',' ) {
	    	sStr[j] = '\0';
	    	if ( c=='\r' ) {
	    	    VPNOTE("A carriage return character has been read.\n"
	    	            "This is an indicator of DOS line endings.\n"
	    	            "Although LEaP treats it as whitespace, other"
	    	            " programs may not.\n"
	    	            "One can convert line endings from DOS to UNIX with"
	    	            " various tools including:\n    dos2unix\n" );
	    	}
	    	break;
	}
	/*
	 *  special case for $-type strings: allow embedded single-char
	 *	tokens, except ';'
	 */
	if ( sStr[0] == '$' ) {
		if ( c == ';' ) {
			zUngetc( c );
	    		sStr[j] = '\0';
	    		break;
		}
		sStr[j] = c;
		continue;
	}


	/*
	 *  '(' is a single-char token UNLESS the identifier collected
	 *  so far names a UNIT, in which case this opens an inline
	 *  AtomMask selection expression: UNIT(mask)
	 */
	if ( c == '(' ) {
	    sStr[j] = '\0';   /* terminate identifier collected so far */
	    OBJEKT oVar = oVariable( sStr );
	    if ( oVar != NULL && iObjectType(oVar) == UNITid ) {
                OBJEKT oResult;
		char *sMask = zsScanMaskLiteral();   /* consumes up to matching ')' */
                bool bSelectResidues = false;
                if (*sMask=='R') {
                    bSelectResidues = true;
                    sMask++;
                }
		LIST lAtoms = lAtomMaskSelect((UNIT)oVar, sMask );
                lAtoms->bFreeChildren = true;
		aAssoc = (ASSOC)oCreate(ASSOCid);
                int iSize = iCollectionSize(lAtoms);
                VP2("AtomMask selected %d atoms\n", iSize );
                //{
                //    OBJEKT oObject;
                //    LISTLOOP llAtoms = llListLoop(lAtoms);
                //    while ( (oObject = oListNext(&llAtoms)) ) {
                //        printf("Type = %c\n",iObjectType(oObject));
                //    }
                //}
                if (bSelectResidues) {
                    RESIDUE rRes = NULL;
                    OBJEKT oObject;
                    LIST lResidues = LIST_from(oCreate(LISTid));
                    LISTLOOP llAtoms = llListLoop(lAtoms);
                    while ( (oObject = oListNext(&llAtoms)) ) {
                        ATOM aAtom = ATOM_from(oObject);
                        RESIDUE rAtomRes = RESIDUE_from(cContainerWithin(aAtom));
                        if (rAtomRes != rRes) {
                            rRes = rAtomRes;
                            ListAddToEnd(lResidues, OBJEKT_from(rAtomRes));
                        }
                    }
                    VP2("AtomMask selected %d residues\n",
                          iCollectionSize(lResidues));
                    lResidues->bFreeChildren = true;
                    ListDestroy(&lAtoms);
                    oResult = OBJEKT_from(lResidues);
                } else {
                    oResult = OBJEKT_from(lAtoms);
                    VP2("AtomMask selected %d atoms\n",
                          iCollectionSize(lAtoms));
                }
                AssocSetName( aAssoc, "atomMask" );
                AssocTakeObject( aAssoc, oResult );
		yylval.aVal = aAssoc;
		printf /*MESSAGE*/("Parsed inline mask on UNIT %s: %s\n", sStr, sMask );
		return(LASSOC);
	    }
	    /* not a UNIT: '(' really is just the single-char token */
	    zUngetc( c );
	    break;
	}


	tok = iOneCharToken( c );
	if ( tok != LNOTSINGLECHAR ) {
		zUngetc( c );
	    	sStr[j] = '\0';
	    	break;
	}
	sStr[j] = c;
    }

    // Unary minus or subtraction. Not in single char set because it
    // may be a leading sign in a negative number or be in filenames.
    if ( sStr[1] == 0 ) {
        if ( sStr[0] == '-' ) return LMINUS;
        if ( sStr[0] == '^' ) return(LPOW);
        if ( sStr[0] == '+' ) return(LPLUS);
        if ( sStr[0] == '/' ) return(LDIV);
    }

    /*
     *  see if it's a number
     */
    bGotExp = false;
    bGotDot = false;
    if ( isdigit((unsigned char)sStr[0]) || sStr[0] == '-' || sStr[0] == '+' || 
				( sStr[0] == '.' && isdigit((unsigned char)sStr[1]) ) ) {

        for (int j=0; j<sizeof(STRING); j++ ) {
            MESSAGE("Thinking NUMBER got: %c\n", sStr[j] );
	    switch ( sStr[j] ) {
		case '\0':
        	    if ( sscanf( sStr, "%lf", &yylval.dVal ) != 1 ) {
			VPWARN(" Couldn't scan NUMBER from (%s)\n", sStr );
			return(LNULL);
		    }
        	    MESSAGE("Parsed a number: %lf\n", yylval.dVal );
        	    return(LNUMBER);
		case '.':
		    if ( bGotDot ) {
			VPERROR("Multiple decimal points in NUMBER-like token"
					" (%s).\n", sStr );
			goto notnum;
		    }
		    if ( bGotExp ) {
			VPERROR("Decimal point follows exponent in NUMBER-like"
					" token (%s).\n", sStr );
			goto notnum;
		    }
        	    bGotDot = true;
		    break;
		case 'e':
		case 'E':
            	    if ( bGotExp ) {
			VPERROR("Multiple exponent indicators ('e') in "
					"NUMBER-like token (%s).\n", sStr );
			goto notnum;
		    }
                    bGotExp = true;
                    break;
		case '+':
		case '-':
                    break;
		default:
            	    if ( !isdigit((unsigned char)sStr[j]) ) {
			goto notnum;
		    }
            	    break;
	    }
	}
    }
notnum:
    /* 
     *  see if it's a string in quotes
     */
    if ( sStr[0] == '"' ) {
        strcpy( yylval.sVal, &sStr[1] );
        MESSAGE("Parsed a STRING: %s\n", sStr );
        return(LSTRING);
    }

    /* 
     *  see if it's a string prefixed w/ '$'
     */
    if ( sStr[0] == '$' ) {
        strcpy( yylval.sVal, &sStr[1] );
        MESSAGE("Parsed a STRING: %s\n", sStr );
        return(LSTRING);
    }

    /*
     *  see if it's a variable reference prefixed w/ '@'
     */
    if ( sStr[0] == '@' ) {
        strcpy( yylval.sVal, &sStr[1] );
        MESSAGE("Parsed a @variable: %s\n", sStr );
        return(LVARIABLE);
    }

    /*
     *  see if it's an eval function)
     */
    if ( !strcmp(sStr,"int") ) {
        strcpy( yylval.sVal, sStr );
        return(LINT);
    }
    if ( !strcmp(sStr,"frac") ) {
        strcpy( yylval.sVal, sStr );
        return(LFRAC);
    }
    if ( !strcmp(sStr,"str") ) {
        strcpy( yylval.sVal, sStr );
        return(LSTR);
    }

    /* Check for parser based eval command first */
    if ( strcasecmp( sStr, "eval" ) == 0 ) return(LEVAL);

                /* LASTLY!!!!!!!! */
    /* 
     *  see if it's a variable/command 
     */
    strcpy( yylval.sVal, sStr );
    strcpy( sPossibleCmd, sStr );
    StringLower( sPossibleCmd );
    
    /* Check if there is an alias that is an exact match */
    if ( (iMax = iVarArrayElementCount( GvaAlias )) ) {
	ALIAS		aAlias;
	aAlias = PVAI( GvaAlias, ALIASt, 0 );
	for (int i=0; i<iMax; i++, aAlias++ ) {
	    if ( strcasecmp( aAlias->sName, sPossibleCmd ) == 0 ) {
	    	strcpy( sPossibleCmd, aAlias->sCommand );
	    }
        }
    }

    /* Check if there is an exact match of the command */
    /* If a command has already been found for this input
    	line, then do not consider the string a command
    	but rather as a STRING variable */
                	
    if ( !bCommandFound ) {
	for (int j=0; cCommands[j].fCallback; j++ ) {
            if ( strcasecmp( sPossibleCmd, cCommands[j].sName ) == 0 ) {
		yylval.fCallback = cCommands[j].fCallback;
		MESSAGE("Parsed a command: %s\n", sStr );
		bCommandFound = true;
                SbLastTokenWasCommand = true;
		return(LCOMMAND);
	    }
        }
    }


                /* If the variable name is null then return LNULL */

    if ( strcmp( sStr, NULLSTR ) == 0 ) 
	return(LNULL);
    
                /* Return the variable name */

    strcpy( yylval.sVal, sStr );
    MESSAGE("Parsed a variable: %s\n", sStr );
    return(LVARIABLE);

}

#define NOBODYSMASKWILLOVERRUNME 4096
/*
 *      sScanMaskLiteral
 *
 *      Called immediately after consuming the opening '(' of an
 *      inline AtomMask expression: UNIT(mask).
 *
 *      Reads raw characters up to the matching ')', treating the
 *      mask text as a foreign sub-grammar:
 *        - whitespace and commas are passed through verbatim
 *          (AtomMask is whitespace-sensitive; LEaP's usual
 *          comma-as-whitespace rule does NOT apply inside here)
 *        - '(' / ')' are depth-counted so AtomMask's own grouping
 *          parens don't terminate the scan early
 *        - '[' ... ']' is treated as opaque: '(' or ')' characters
 *          inside a bracket class (fnmatch-style patterns) do NOT
 *          affect paren depth and are not otherwise interpreted
 *
 *      Returns the mask text with the closing ')' already consumed.
 */
static char *
zsScanMaskLiteral()
{
static char sBuf[NOBODYSMASKWILLOVERRUNME];
int     j = 0;
int     iDepth = 1;       /* opening '(' already consumed by caller */
bool    bInBracket = false;
char    c;

    for ( ;; ) {

	if ( j >= NOBODYSMASKWILLOVERRUNME - 1 )
	    DFATAL("AtomMask expression too long");

	c = cGetChar();

	if ( c == '\0' )
	    DFATAL("unterminated AtomMask expression, missing ')'");

	if ( bInBracket ) {
	    sBuf[j++] = c;
	    if ( c == ']' )
		bInBracket = false;
	    continue;
	}

	if ( c == '[' ) {
	    bInBracket = true;
	    sBuf[j++] = c;
	    continue;
	}

	if ( c == '(' ) {
	    iDepth++;
	    sBuf[j++] = c;
	    continue;
	}

	if ( c == ')' ) {
	    iDepth--;
	    if ( iDepth == 0 ) {
		sBuf[j] = '\0';
		return(sBuf);
	    }
	    sBuf[j++] = c;
	    continue;
	}

	sBuf[j++] = c;
    }
}

//----------- Molecular object hierarchy parser ---------------
// main function is zoFindObject() which uses the following
// internal tokenizer helpers. If not resolved, the full token
// is returned as a string.
//
// Parser syntax:
// Initial object resolved as variable, which can be indirect:
//      $var1 resolves to varible named by var1
// input stream requires 2 
// <COLLECTION>.<integer> -- single depth list/vector member
// <UNIT>.<RESIDUE>.<ATOM>

typedef enum {
    CTOK_ERROR    = -1,
    CTOK_INT      =  0,
    CTOK_STRING   =  1,
    CTOK_CHAINRES =  2,
} eContTokenType;

static int
zsNextObjectToken(char *string, char **next_return)
{
char    *next;
    for (next = string; *next; next++) {
        if (strchr(".@%", *next)) {
            char sep = *next;
            *next = 0;
            *next_return = next + 1;
            return sep;
        }
    }
    return 0;
}

/* Returns TRUE if s is a non-empty string of digits */
static bool
zbResidNumeric(const char *s)
{
    if (!*s) return false;
    for (const char *p = s; *p; p++)
        if (!isdigit((unsigned char)*p)) return false;
    return true;
}

/*
 * Resolve a token (which may have a leading '$' and/or a ':' separator)
 * into a typed value.
 *
 * Literals are untyped: digit-only → CTOK_INT, otherwise → CTOK_STRING,
 *   except that a ':' anywhere means CTOK_CHAINRES (left side may be empty).
 * Variables are typed by their object type: ODOUBLEid → CTOK_INT,
 *   OSTRINGid → CTOK_STRING or CTOK_CHAINRES if value contains ':'.
 *   A numeric-looking OSTRINGid is NOT converted to CTOK_INT.
 *
 * For CTOK_CHAINRES: sOut receives chain (may be ""), *piValue receives resid.
 * For CTOK_STRING:   sOut receives the name.
 * For CTOK_INT:      *piValue receives the index.
 * Emits VPWARN and returns CTOK_ERROR on any problem.
 */
static eContTokenType
eResolveToken(const char *cPPos, char *sOut, int *piValue)
{
OBJEKT      oVar;
const char  *sVal;
char        sBuf[MAXSTRINGLENGTH];
char        *cPColon;
bool        bFromVar = false;

    /* --- expand leading '$' if present --- */
    if (*cPPos == '$') {
        cPPos++;
        oVar = oVariable(cPPos);
        if (!oVar) {
            VPWARN("Variable %s not found\n", cPPos);
            return CTOK_ERROR;
        }
        if (iObjectType(oVar) == ODOUBLEid) {
            *piValue = (int)dODouble(oVar);
            return CTOK_INT;
        }
        if (iObjectType(oVar) != OSTRINGid) {
            VPWARN("Variable %s is not a string or numeric\n", cPPos);
            return CTOK_ERROR;
        }
        sVal = sOString(oVar);
        bFromVar = true;
    } else {
        sVal = cPPos;
    }

    /* --- check for 'chain:resid' pattern (chain may be empty) --- */
    strcpy(sBuf, sVal);
    cPColon = strchr(sBuf, ':');
    if (cPColon) {
        const char *sResid;
        *cPColon = '\0';
        strcpy(sOut, sBuf);     /* chain -- preserve as-is, may be "" or "042" etc */
        sResid = cPColon + 1;

        /* right side: may itself have a leading '$' */
        if (*sResid == '$') {
            sResid++;
            oVar = oVariable(sResid);
            if (!oVar) {
                VPWARN("Variable %s not found\n", sResid);
                return CTOK_ERROR;
            }
            if (iObjectType(oVar) == ODOUBLEid) {
                *piValue = (int)dODouble(oVar);
            } else if (iObjectType(oVar) == OSTRINGid) {
                const char *sR = sOString(oVar);
                if (!zbResidNumeric(sR)) {
                    VPWARN("Resid variable %s value '%s' is not numeric\n", sResid, sR);
                    return CTOK_ERROR;
                }
                *piValue = atoi(sR);
            } else {
                VPWARN("Resid variable %s is not a string or numeric\n", sResid);
                return CTOK_ERROR;
            }
        } else {
            if (!zbResidNumeric(sResid)) {
                VPWARN("Resid '%s' is not numeric\n", sResid);
                return CTOK_ERROR;
            }
            *piValue = atoi(sResid);
        }
        return CTOK_CHAINRES;
    }

    /* --- no colon ---
     * vars: OSTRINGid is always a name, even if digits only
     * literals: digit-only → int, otherwise → string             */
    if (!bFromVar && zbResidNumeric(sVal)) {
        *piValue = atoi(sVal);
        return CTOK_INT;
    }
    strcpy(sOut, sVal);
    return CTOK_STRING;
}


OBJEKT
zoFindObject(char *sIdentifier)
{
CONTAINER       cCont;
int             iType, iValue;
OBJEKT          oObj;
STRING          sLine;
STRING          sOut;
char            cSep, cSepNext;
char            *cPPos, *cPNext;
LOOP            lRes;
RESIDUE         rRes;
eContTokenType  eTok;

    strcpy(sLine, sIdentifier);

    cPPos    = sLine;
    cSepNext = zsNextObjectToken(cPPos, &cPNext);

    if (*cPPos == '$') {
        cPPos++;
        OBJEKT oVar = oVariable(cPPos);
        if (!oVar) {
            VPWARN("Variable %s not found\n", cPPos);
            return NULL;
        }
        if (iObjectType(oVar) != OSTRINGid) {
            VPWARN("Variable %s is not a string\n", cPPos);
            return NULL;
        }
        cPPos = sOString(oVar);
    }
    oObj = oVariable(cPPos);
    if (!cSepNext) return oObj;

    cPPos    = cPNext;
    cSep     = cSepNext;
    cSepNext = zsNextObjectToken(cPPos, &cPNext);

    if ( bObjectInClass(oObj, COLLECTIONid ) && isdigit((unsigned char)*cPPos))  {
        int index = atoi(cPPos);
        ASSOC aAssocElem;
        int i=0;
        LISTLOOP ll = llListLoop((LIST)oObj);
        while ( (aAssocElem = (ASSOC)oListNext(&ll)) ) {
            i=i+1;
            if (i == index) return oAssocObject(aAssocElem);
        }
        VPWARN("Index %d out of bounds trying to index object\n", index );
        return NULL;
    }

    while (cSep) {
        if (!bObjectInClass(oObj, CONTAINERid)) return NULL;
        iType = iObjectType(oObj);
        cCont = NULL;

        switch (cSep) {
        case '.':
            eTok = eResolveToken(cPPos, sOut, &iValue);
            switch (eTok) {
            case CTOK_ERROR:
                return NULL;

            case CTOK_CHAINRES:
                if (iObjectType(oObj) != UNITid) {
                    VPWARN("ChainId:Resid selection only valid within UNIT\n");
                    return NULL;
                }
                lRes = lLoop(oObj, RESIDUES);
                while ((rRes = (RESIDUE)oNext(&lRes))) {
                    if (iResiduePdbSequence(rRes) == iValue &&
                            !strcmp(sOut, rRes->sChainId)) {
                        cCont = (CONTAINER)rRes;
                        break;
                    }
                }
                if (!cCont) {
                    VPWARN("ChainId:Resid %s:%d not found in %s %s\n",
                            sOut, iValue, sObjectIndexType(iType), sContainerName(oObj));
                    return NULL;
                }
                break;

            case CTOK_INT:
                /* try flat index first; fall through to name on failure */
                cCont = cContainerFindSequence((CONTAINER)oObj, DIRECTCONTENTS, iValue);
                if (cCont) break;
                /* not found by index -- try as name using the original token */
                strcpy(sOut, cPPos);
                /* FALLTHROUGH */
            case CTOK_STRING:
                cCont = cContainerFindName((CONTAINER)oObj, DIRECTCONTENTS, sOut);
                if (!cCont) {
                    VPWARN("Name/index %s not found in %s %s\n",
                            sOut, sObjectIndexType(iType), sContainerName(oObj));
                    return NULL;
                }
                break;
            }
            break;

        case '@':
            /* '@' always means name -- no numeric interpretation.
             * Only '$' expansion allowed; var must be OSTRINGid.        */
            if (*cPPos == '$') {
                const char *sVarName = cPPos + 1;
                OBJEKT oVar = oVariable(sVarName);
                if (!oVar) {
                    VPWARN("Variable %s not found\n", sVarName);
                    return NULL;
                }
                if (iObjectType(oVar) != OSTRINGid) {
                    VPWARN("Atom name variable %s is not a string\n", sVarName);
                    return NULL;
                }
                strcpy(sOut, sOString(oVar));
            } else {
                strcpy(sOut, cPPos);
            }
            /* check for group notation if parent is UNIT */
            if (iObjectType(oObj) == UNITid) {
                OBJEKT oGroup = (OBJEKT)lUnitGroup((UNIT)oObj, sOut);
                if (oGroup) return oGroup;
            }
            cCont = cContainerFindName((CONTAINER)oObj, DIRECTCONTENTS, sOut);
            if (!cCont) {
                VPWARN("Name %s not found in %s %s\n",
                        sOut, sObjectIndexType(iType), sContainerName(oObj));
                return NULL;
            }
            break;

        case '%':
            /* Legacy: resid without chainid. Must be UNIT, must be int. */
            if (iObjectType(oObj) != UNITid) {
                VPWARN("Selection delimiter %% is only valid within UNIT\n");
                return NULL;
            }
            eTok = eResolveToken(cPPos, sOut, &iValue);
            if (eTok == CTOK_ERROR) return NULL;
            if (eTok != CTOK_INT) {
                VPWARN("Resid selection requires an integer\n");
                return NULL;
            }
            lRes = lLoop(oObj, RESIDUES);
            while ((rRes = (RESIDUE)oNext(&lRes))) {
                if (iResiduePdbSequence(rRes) == iValue) {
                    cCont = (CONTAINER)rRes;
                    break;
                }
            }
            if (!cCont) {
                VPWARN("PDB resId %d not found in %s %s\n",
                        iValue, sObjectIndexType(iType), sContainerName(oObj));
                return NULL;
            }
            break;

        default:
            VPFATAL("Programming error (cSep %c) %s %s\n", cSep, __FILE__, __func__);
            return NULL;
        }

        oObj     = (OBJEKT)cCont;
        cPPos    = cPNext;
        cSep     = cSepNext;
        if (cSep) cSepNext = zsNextObjectToken(cPPos, &cPNext);
    }

    return oObj;

}
 
static OBJEKT
zoGetObject(char *sName)
{
    printf("get object %s\n",sName);
    OBJEKT oObj;
    OSTRING osString;

    oObj = zoFindObject(sName);
    if (oObj != NULL) return REF(oObj);

    /*
     * Not an object variable or resolvable object reference.
     * Return the whole string as an OSTRING.
     */
    osString = (OSTRING)oCreate(OSTRINGid);
    STRING sResult;
    ReplaceVarPattern(sName, sResult, sizeof(sResult));
    OStringDefine(osString, sResult);
    return (OBJEKT)osString;
}

static ASSOC
zaEvalBinary(char operator, ASSOC arg1, ASSOC arg2 ) {
OBJEKT oResult = NULL;
int iType1 = iObjectType(oAssocObject(arg1));
int iType2 = iObjectType(oAssocObject(arg2));
OBJEKT oArg1 = oAssocObject(arg1);
OBJEKT oArg2 = oAssocObject(arg2);
    if (iType1 == OSTRINGid && iType2 == OSTRINGid) {
        if (operator == '+') {
            oResult = oCreate(OSTRINGid);
            OStringCopy((OSTRING)oResult, (OSTRING)oArg1);
            OStringConcat((OSTRING)oResult, (OSTRING)oArg2);
        } else oResult = NULL;
    }
    else if (iType1 == OSTRINGid && iType2 == ODOUBLEid) {
        if (operator == '+') {
            STRING sLine;
            double d2 = dODouble(oAssocObject(arg2));
            sprintf(sLine,"%s%lg",sOString(oArg1), d2);
            oResult = oCreate(OSTRINGid);
            OStringDefine((OSTRING)oResult, sLine);
        } else oResult = NULL;
    }
    else if (iType1 == ODOUBLEid && iType2 == ODOUBLEid) {
        oResult = oCreate(ODOUBLEid);
        double d1 = dODouble(oArg1);
        double d2 = dODouble(oArg2);
        switch (operator) {
            case '*': ODoubleSet(oResult,d1*d2); break;
            case '/': ODoubleSet(oResult,d1/d2); break; // check for divide by zero? Or just let NaN take over?
            case '+': ODoubleSet(oResult,d1+d2); break;
            case '-': ODoubleSet(oResult,d1-d2); break;
            case '^': ODoubleSet(oResult,pow(d1,d2)); break;
            default: Destroy(&oResult); break;
        }
    }

    if (!oResult) {
        VPFATAL("Invalid operator/type: %s %c %s\n",sObjectIndexType(iType1),
                         operator, sObjectIndexType(iType2));
        return NULL;
    }
    ASSOC aResult = (ASSOC)oCreate(ASSOCid);
    AssocTakeObject( aResult, oResult);
    return aResult;
}


    
/*
 *================================================================
 *
 *	Public routines
 */




/*
 *	ParseArguments
 *
 *	Parse the arguments that the user passes to LEaP from
 *	the command line arguments.
 */
void
ParseArguments( int argc, char *argv[] )
{
char		c;
extern	char	*optarg;

    while ( (c = getopt( argc, argv, "hsI:f:" )) != (char)(EOF) ) {
	switch (c) {
	    case 'h':
		printf( "Usage: %s [options]\n", argv[0] );
		printf( "Options:\n" );
		printf( " -h         Generate this message.\n" );
		printf( " -s         Ignore %s startup file.\n", LEAPRC );
		printf( " -I {dir}   Add {dir} to search path.\n" );
		printf( " -f {file}  Source {file}.\n" );
		exit(0);
	    case 's':
		printf( "-s: Ignoring all %s startup files.\n", LEAPRC );
		SbUseStartup = false;
		break;
	    case 'I':
		printf( "-I: Adding %s to search path.\n", optarg );
		BasicsAddDirectory( optarg, 1 );
		break;
	    case 'f':
		printf( "-f: Source %s.\n", optarg );
		if ( iFirstSource == 0 ) {
			SbFirstSourceFiles = (STRING *)MALLOC(sizeof(STRING));
			iFirstSource = 1;
		} else {
			iFirstSource++;
			SbFirstSourceFiles = (STRING *)REALLOC(SbFirstSourceFiles,
					iFirstSource * sizeof(STRING));
		}
		strcpy( SbFirstSourceFiles[iFirstSource-1], optarg );
		break;
	}
    }
}






/*
 *	ParseInit
 *
 *	Initialize the parser.
 *	If SbStartup is TRUE the execute the LEAPRC script.
 */
void
ParseInit( RESULTt *rPResult )
{
FILE	*fStartup;
int	iFile;

    VPTRACEENTER("ParseInit" );
    VP0("\nWelcome to LEaP!\n" );

#ifdef  DEBUG
    VP0("LEaP is running in DEBUG mode!\n" );
#endif

    HelpInitialize();

                /* Initialize the first file in the stack to be stdin */
                
    GfaInputFileStack[0] = NULL;
    GbExecute = bBlockCreate();

                /* Create a few OBJEKTs that will be used by the parser */

    aDummy = (ATOM)oCreate(ATOMid);
    ContainerSetName( aDummy, "DUMMY" );
    GplAllParameters = plParmLibCreate();

    VariablesInit();
    GrMainResult.iCommand = CNONE;
    rPResult->iCommand = CNONE;
    
                /* Parse the LEAPRC file if bUseStartup is TRUE */

    if ( SbUseStartup ) {
        fStartup = FOPENNOCOMPLAIN( LEAPRC, "r" );
        if ( fStartup == NULL ) {
	    VP0("(no %s in search path)\n", LEAPRC );
	} else {
	    /*
	     *  source the leaprc
	     */
	    VP0("Sourcing %s: %s\n", LEAPRC, GsBasicsFullName );
	    INPUTPUSHFILE( fStartup );
	    GiClipPrompts = 1;
	    while ( fINPUTFILE() != NULL ) {
	        yyparse();
		if ( GrMainResult.iCommand == CQUIT ) {
	    	    if ( fINPUTFILE() != NULL )
	    	    	fclose( fINPUTFILE() );
	    	    INPUTPOPFILE();
		}
	    }
	    *rPResult = GrMainResult;
	    GiClipPrompts = 0;
	}
    }

		/* Parse the source files specified */
		/* on the command line using the -f options */

    for (iFile=0; iFile<iFirstSource; iFile++) {
	fStartup = FOPENCOMPLAIN( SbFirstSourceFiles[iFile], "r" );
	if ( fStartup != NULL ) {
	    VP0("Sourcing: %s\n", GsBasicsFullName );
	    INPUTPUSHFILE( fStartup );
	    VPTRACE("After push GiInputFileStackPos = %d\n",
	               GiInputFileStackPos );
	    GiClipPrompts = 1;
	    while ( fINPUTFILE() != NULL ) {
		yyparse();
		if ( GrMainResult.iCommand == CQUIT ) {
	    	    if ( fINPUTFILE() != NULL )
	    	    	fclose( fINPUTFILE() );
	    	    INPUTPOPFILE();
	    	    VPTRACE("After pop GiInputFileStackPos = %d\n",
	    	              GiInputFileStackPos );
		}
	    }
	    *rPResult = GrMainResult;
	    GiClipPrompts = 0;
	} else {
	    exit(41);
	}
    }
    if ( SbFirstSourceFiles )
	FREE( SbFirstSourceFiles );
    VPTRACEEXIT("ParseInit" );
}




/*
 *	ParseBlock
 *
 *	Parse a BLOCK containing one complete command.
 *	Return in rResult the result of the command.
 */
void
ParseBlock( BLOCK bBlock, RESULTt *rPResult )
{
		/* Set up the BLOCK from which to read the command */

    VPTRACEENTER("ParseBlock" );
    MESSAGE("Parsing block: %s\n", sBlockText(bBlock) );
    GbCommand = bBlock;
    BlockResetRead( GbCommand );

    GrMainResult.iCommand = CNONE;

		/* Parse the BLOCK */
		/* Keep parsing as long as the 'execute' command */
		/* keeps setting the GLOBAL variable GbContinueParsing */

    do {
	yyparse();
	if ( GrMainResult.iCommand == CQUIT ) {
	    if ( fINPUTFILE() != NULL )
	    	fclose( fINPUTFILE() );
	    INPUTPOPFILE();
	    VPTRACE("After pop GiInputFileStackPos = %d\n",
	              GiInputFileStackPos );
	}
    } while ( fINPUTFILE() != NULL );

	/* Reset the bCommandFound variable as a new command may be
	   available in a new input */
    bCommandFound = false;
    
	/* Copy the RESULT from the global result variable */

    *rPResult = GrMainResult;

    VPTRACEEXIT("ParseBlock" );
}




extern void CMAPDestroy(void);
/*
 *	ParseShutdown
 *
 *	Shutdown the parser, release all variables setup in ParseInit
 *	Only needed if debugging memory mgt, since prog mem is all
 *	freed when process exits anyway.
 */
void
ParseShutdown()
{
	if ( !iMemDebug )

	Destroy( (OBJEKT *)&aDummy );
	VariablesDestroy();
	ParmLibDestroy( &GplAllParameters );

	BlockDestroy( &GbExecute );

	HelpShutdown();
}
