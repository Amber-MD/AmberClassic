/*
 *      File: assoc.c
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
 *      Class: 
 *              ASSOC
 *      Superclass: 
 *              OBJEKT
 *
 *      Description:
 *
 *              An ASSOC (short for association) is an object that
 *              stores a short name and an OBJEKT.
 *              It associates the name with the OBJEKT.
 *              This is used by the command line interface to
 *              store objects with the names that are bound to them.
 *
 */



#include	"basics.h"

#include        "classes.h"



/*
-------------------------------------------------------------------

        Define static variables here.
*/




/*
===================================================================

        Define private routines here.
*/


/*
*******************************************************************

        Define public message routines here.
*/




/*
 *      aAssocCreate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Create an ASSOC.
 */
ASSOC   
aAssocCreate()
{
ASSOC   aNew;

    aNew = (ASSOC)MALLOC(sizeof(ASSOCt) );
    AssocSetName( aNew, "" );
    AssocSetObject( aNew, NULL );
    return aNew;
}



/*
 *      AssocDestroy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Destroy an ASSOC.
 *      DEREF the OBJEKT within it.
 */
void
AssocDestroy( ASSOC *aPAssoc )
{
    DEREF( oAssocObject(*aPAssoc) );
    FREE( *aPAssoc );
    *aPAssoc = NULL;
}



/*
 *      AssocDescribe
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Describe the ASSOC.
 */
void
AssocDescribe( ASSOC aAssoc )
{
	if ( strlen(sAssocName(aAssoc)) == 0)
		VP0(" <assoc>: " );
	else
		VP0("%s: ", sAssocName(aAssoc) );
	Describe( oAssocObject(aAssoc) );
}


// Create a vector as LIST of ASCOC of DOUBLE from VECTOR 
OBJEKT
CreateAssocVector(VECTOR vVector) {
    LIST lResult = (LIST)oCreate(LISTid);

    OBJEKT oDouble = oCreate(ODOUBLEid);
    ODoubleSet( oDouble,  vVector.dX );
    ASSOC aAssoc = (ASSOC)oCreate(ASSOCid);
    AssocSetName( aAssoc, "x" );
    AssocSetObject( aAssoc, oDouble);
    ListAddToEnd( lResult, (OBJEKT)aAssoc );

    oDouble = oCreate(ODOUBLEid);
    ODoubleSet( oDouble,  vVector.dY );
    aAssoc = (ASSOC)oCreate(ASSOCid);
    AssocSetName( aAssoc, "y" );
    AssocSetObject( aAssoc, oDouble);
    ListAddToEnd( lResult, (OBJEKT)aAssoc );

    oDouble = oCreate(ODOUBLEid);
    ODoubleSet( oDouble,  vVector.dZ );
    aAssoc = (ASSOC)oCreate(ASSOCid);
    AssocSetName( aAssoc, "z" );
    AssocSetObject( aAssoc, oDouble);
    ListAddToEnd( lResult, (OBJEKT)aAssoc );

    return (OBJEKT)lResult;
}
