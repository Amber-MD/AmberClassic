/*
 *      File: byteArray.c
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
 *              BYTEARRAY
 *      Superclass: 
 *              OBJEKT
 *
 *      Description:
 *
 *              A BYTEARRAY is a basic object that contains no other
 *              objects.  It can coorespond to any primative C type.
 *              The point of using a BYTEARRAY is that often objects
 *              contain other objects.  In order to allow
 *              objects to contain integers, strings, floats etc.
 *              there must either be an object type defined for each
 *              of these or one general object which can encompass them 
 *              all.  The basic property that the BYTEARRAY
 *              has is that it knows how to 'Destroy' itself.
 *              This allows more comples objects to destroy 
 *              this object just by sending it the 'Destroy' message.
 *
 */





#include	"basics.h"

#include	"classes.h"




/*
*******************************************************************

        Define public message routines here.
*/



/*
 *      baByteArrayCreate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Create a new BYTEARRAY of the size iSize
 */
 
BYTEARRAY
baByteArrayCreate( int iSize )
{
int             iTotalSize;
BYTEARRAY       baPtr;

    iTotalSize = iSize + sizeof(BYTEARRAYHEADERt);
    baPtr = (BYTEARRAY)MALLOC(iTotalSize );
    baPtr->bahHeader.iDataSize = iSize;
    return baPtr;
}




/*
 *      ByteArrayDestroy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Destroy the byte array, release the memory.
 */
void
ByteArrayDestroy( BYTEARRAY *baPByteArray )
{
    VERIFYOBJEKT( *baPByteArray, BYTEARRAYid );
    FREE(*baPByteArray);
    *baPByteArray = NULL;
}


/*
 *      ByteArrayDescribe
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return a short string that describes the
 *      byte array.
 */
void
ByteArrayDescribe( BYTEARRAY baByteArray )
{
    VERIFYOBJEKT( baByteArray, BYTEARRAYid );
    VP0("BYTEARRAY size: %d\n", iByteArraySize(baByteArray) );
    VP0("BYTEARRAY= |%s|\n", PByteArray(baByteArray) );
}




/*
 *      baByteArrayDuplicate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return a pointer to a duplicate of the BYTEARRAY
 *      OBJEKTS contained within the BYTEARRAY will not
 *      be copied, BYTEARRAYs should NEVER have internal
 *      OBJEKTS!
 */
BYTEARRAY
baByteArrayDuplicate( BYTEARRAY baOld )
{
OBJEKT          baNew;
int             iSize;

    iSize = iByteArraySize(baOld);
    baNew = oCreateSize( BYTEARRAYid, iSize );
    memcpy( PByteArray(baNew), PByteArray(baOld), iSize );
    return((BYTEARRAY)baNew);
}

