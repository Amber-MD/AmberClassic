/*
 *	File:	message.c
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
 *	Description:
 *		Store messages in a long list.
 *		Keep track of what function they come from
 *		using 'function.h'.
 */



#include	"basics.h"

#include	"function.h"

#include	"mess.h"



MESS		GmMesss = NULL;
MESS		GmLastMess = NULL;


/*
 *	mMessAdd
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Add a message to the list, determine which function it is
 *	in and set the messages Print flag to that of its function.
 *
 *	The first part of the message must look like:
 *		@(\w+) (\w+)|(text)
 *		Where $1 == Filename.
 *		      $2 == line number.
 *		      $3 == Mess text.
 *	Or:
 *		#(\w+) (\w+)
 *		Where $1 == Filename.
 *		      $2 == line number.
 *	Or:
 *		!(text)
 *		      $1 == Text.
 *
 *	The first example is a normal message, the second example
 *	comes from a stack trace. Stack trace messages are ALWAYS
 *	printed.
 *
 */
MESS
mMessAdd( char *sText )
{
STRING		sFile;
int		iLine;
bool		bPrint, bTrace;
MESS		mMess;
char		*cPText;
char		cType;
int		iScanned, iFunction;

    iScanned = sscanf( sText, "%c%s %d", &cType, sFile, &iLine );

    cType = sText[0];
    switch ( cType ) {
	case PRINT_MESSAGE:
	case PRINT_TRACE:
	    iScanned = sscanf( sText, "%c%s %d", &cType, sFile, &iLine );
	    if ( iScanned != 3 ) {
		DFATAL("Invalid (file:line number) combination in message: |%s|\n",
			sText );
	    }
	    iFunction = iFunctionFindWithFilenameLine( sFile, iLine );
	    break;
	case PRINT_ALWAYS:
	    iFunction = NO_FUNCTION;
	    break;
	default:
	    cType = PRINT_ALWAYS;
	    iFunction = NO_FUNCTION;
	    break;
    }

    mMess = (MESS)MALLOC(sizeof(MESSt) ); 
    cPText = (char*)MALLOC(strlen(sText)+1 );

    mMess->cType = cType;
    mMess->iFunction = iFunction;
    strcpy( cPText, sText );
    mMess->cPText = cPText;
    switch ( cType ) {
	case PRINT_MESSAGE:
	    mMess->bPrint = bFunctionPrint(iFunction);
	    break;
	case PRINT_TRACE:
	case PRINT_ALWAYS:
	    mMess->bPrint = true;
	    break;
    }

    mMess->mNext = NULL;
    if ( GmMesss == NULL ) {
	GmMesss = mMess;
	GmLastMess = mMess;
    } else {
	GmLastMess->mNext = mMess;
	GmLastMess = mMess;
    }

    return mMess;
}



/*
 *	mMessLoop
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return a MESS that can be passed to MessNext
 *	to return subsequent MESSs.
 */
MESS
mMessLoop()
{
    return((MESS)NULL);
}



/*
 *	mMessNext
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Change (mPMess) to point to the next MESS.
 */
MESS
mMessNext( MESS *mPMess )
{
    if ( *mPMess == NULL ) {
	*mPMess = GmMesss;
    } else {
	*mPMess = (*mPMess)->mNext;
    }

    return *mPMess;
}



