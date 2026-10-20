/*
 *      File: collection.h
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
 *              COLLECTION
 *      Superclass: 
 *              OBJEKT
 *
 */
 
#ifndef COLLECTION_H
#define COLLECTION_H

#include "objekt.h"

/*
-----------------------------------------------------------------------

        Define object typedefs here.
        
        Object typedef MUST include the superclass object as
        its first structure element.
*/


typedef struct  {
	OBJEKTt         oObjectHeader;
	int             iElementCount;
} COLLECTIONt;

typedef COLLECTIONt     *COLLECTION;    



/*
======================================================================

        Define object messages here.
        
        There must be at least a Create, Destroy, and Describe message.
        Hook into the messages of the superclasses so that
        when the message is sent to the most primative superclass
        of this class that it will eventually make it into these routines.
*/



#define iCollectionSize( Cx ) ( ((COLLECTION)(Cx))->iElementCount )
#define CollectionSetSize( Cx, S ) ( ((COLLECTION)Cx)->iElementCount = S )

extern COLLECTION		cCollectionCreate( int iType );
extern void			CollectionDestroy( COLLECTION *cPCollect );
extern void			CollectionDescribe( COLLECTION cCollect );
extern GENP			PCollectionLoop( COLLECTION cCollect );
extern OBJEKT			oCollectionNext( COLLECTION cCollect, 
					GENP *PPNode );
extern COLLECTION		cCollectionDuplicate( COLLECTION cCol );

#endif /* COLLECTION_H */
