/*
 *      File:   parmLib.c
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
 *              A PARMLIB maintains an ordered list of PARMSETs
 *              which are searched sequentially for parameters.
 */





#include	"basics.h"

#include        "classes.h"

#include        "parmLib.h"



PARMLIB	GplDefaultParmLib = NULL;



/*
 *================================================================
 *
 *	Private routines
 *
 */

/*
 *	zbParmLibOK
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Check if the variable passed is a valid PARMLIB.
 */
static bool
zbParmLibOk( PARMLIB pl )
{
bool		bOk;

    bOk = ( pl != NULL );
#if 0
    if ( !bOk ) {
	VP0("There are no PARMSETs loaded!\n" );
    }
#endif
    return bOk;
}




/*
 *---------------------------------------------------------------
 *
 *	Public routines.
 */



/*
 *      plParmLibCreate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return a PARMLIB
 */
PARMLIB
plParmLibCreate()
{
PARMLIB     plTemp;

    plTemp = (PARMLIB)MALLOC(sizeof(PARMLIBt) );

    plTemp->lParmSets = (LIST)oCreate(LISTid);
    return plTemp;
}


/*
 *      ParmLibDestroy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Destroy a PARMLIB
 */
void
ParmLibDestroy( PARMLIB *plPLib )
{
PARMSET		psSet;

    ParmLibParmSetLoop( *plPLib );
    while( bParmLibNextParmSet( *plPLib, &psSet ) ) {
	ParmSetDestroy( &psSet );
    }
    Destroy( (OBJEKT *)&((*plPLib)->lParmSets) );
    FREE( *plPLib );
    *plPLib = NULL;
}


/*
 *      ParmLibAddParmSet
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add a PARMSET to the PARMLIB
 */
void    
ParmLibAddParmSet( PARMLIB plLib, PARMSET psSet )
{
    if ( !zbParmLibOk(plLib) ) 
	return;
    ListAdd( plLib->lParmSets, (OBJEKT)psSet );
}




/*
 *      ParmLibParmSetLoop
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Create a loop over all PARMSETs in the PARMLIB.
 */
void    
ParmLibParmSetLoop( PARMLIB plLib )
{
    if ( !zbParmLibOk(plLib) ) 
	return;

    plLib->llParmSetLoop = llListLoop( plLib->lParmSets );
}





/*
 *      bParmLibNextParmSet
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return each PARMSET one at a time.
 *      Return false if there are no more PARMSETs.
 */
bool
bParmLibNextParmSet( PARMLIB plLib, PARMSET *psPSet )
{
    if ( !zbParmLibOk(plLib) ) 
	return false;
    *psPSet = (PARMSET)oListNext(&(plLib->llParmSetLoop));
    if ( *psPSet == NULL ) 
	return false;
    return true;
}


