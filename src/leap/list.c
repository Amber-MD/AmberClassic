/*
 *	File:	list.c
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
 *	Description:
 *		List is a doubly linked list of OBJEKTS.
 *
 *      **NOTE**: The original design did not free object references
 *              in ListDestroy() (I don't know why). This means that
 *              DEREF(list) will cause memory leaks esp in parser
 *              list arguments where a general object can be a list.
 *              Set bFreeChildren=true for automatic content DEREF()
 *
 *	Note:
 *	        ListAdd() inserts at the head node, reversing
 *		order. Use ListAddToEnd() to maintain insertion order.
 *
 *	WARNING:
 *		NEVER NEVER NEVER remove OBJEKTs from
 *		a LIST that you are currently looping over.
 *		This can lead to inconsistencies.
 *		The solution is to change the structure 
 *		LISTLOOP so that it contains a pointer
 *		to the next element in the list instead
 *		of just to the current one which is to
 *		be Removed from the list.
 *
 *TODO:	Change LISTLOOP to contain the pointer to the next
 *TODO:	element in the list.
 */




#include	"basics.h"

#include	"classes.h"
#include        "defaults.h"


/*
 *-------------------------------------------------------------------------
 *
 *	Private routines
 */



/*
 *	searchList
 *
 *	Find the node that points to the object pointed to by the
 *	argument.  Return a pointer to the pointer to the node
 *	which points to the object.
 */
static NODEP
searchList( NODEP nPPList, OBJEKT oObj )
{
        for (NODEP nPCur = nPPList; nPCur; nPCur = nPCur->nPNextNode)
            if ( nPCur->PObject == oObj ) return nPCur;
        return NULL;
}
        


/*
 *-----------------------------------------------------------------
 *
 *        Public routines
 *        
 */


/*
 *	lListCreate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return a new list.
 */
LIST
lListCreate()
{
LIST    lList;

        lList = (LIST)CALLOC(sizeof(LISTt) );
        return lList;
}


/*
 *      ListDestroy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Destroy the list.  The caller is responsible for destroying the
 *      contents of the list. UNLESS bFreeChildren has been set.
 */
void
ListDestroy( LIST *lPList )
{
NODEP   nPNext, nPFree;
LIST    l = *lPList;

        /* Destroy the LIST nodes */

    nPNext = l->nPFirstNode;
    while ( nPNext != NULL ) {
        if (l->bFreeChildren && nPNext->PObject) DEREF(nPNext->PObject);
        nPFree = nPNext;
        nPNext = nPNext->nPNextNode;
        
                /* Free the node memory. */
        FREE(nPFree);
    }

        /* Now destroy the LIST itself */
    FREE( l );
    *lPList = NULL;

}



/*
 *        ListAdd
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *        Add an element to the head of a list.
 *        Simply return if the user tries to add NULL since
 *        NULLs take up no space.
 */
void
ListAdd( LIST lList, OBJEKT oObj )
{
NODEP   nPTempList;

    if ( oObj == NULL ) 
	return;

                /* Insert the new node at the head of the list */

    nPTempList = lList->nPFirstNode;
    lList->nPFirstNode = (NODEP)MALLOC(sizeof(NODE) );

                /* Keep up to date the last node of the list */

    if ( nPTempList == NULL ) {
        lList->nPLastNode = lList->nPFirstNode;
    }
    else nPTempList->nPPrevNode = lList->nPFirstNode;

    (lList->nPFirstNode)->PObject = oObj;
    (lList->nPFirstNode)->nPNextNode = nPTempList;
    (lList->nPFirstNode)->nPPrevNode = NULL;
    CollectionSetSize( lList, iCollectionSize(lList)+1 );

                /* Reference the object because it is being added */
    REF( oObj );
}





/*
 *      ListAddToEnd
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Add an element to the tail of a list.
 *      Simply return if the user tries to add NULL since
 *      NULLs take up no space.
 */
void
ListAddToEnd( LIST lList, OBJEKT oObj )
{
NODEP   nPTempList;

    if ( oObj == NULL ) 
	return;

                /* Insert the new node at the tail of the list */

    nPTempList = lList->nPLastNode;
    lList->nPLastNode = (NODEP)MALLOC(sizeof(NODE) );

                /* Keep up to date the last node of the list */

    if ( nPTempList == NULL ) {
        lList->nPFirstNode = lList->nPLastNode;
    } else {
        nPTempList->nPNextNode = lList->nPLastNode;
    }

    (lList->nPLastNode)->PObject = oObj;
    (lList->nPLastNode)->nPNextNode = NULL;
    (lList->nPLastNode)->nPPrevNode = nPTempList;
    CollectionSetSize( lList, iCollectionSize(lList)+1 );

                /* Reference the object because it is being added */
    REF( oObj );
}




/*
 *	ListAddUnique
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Add the object to the LIST only if it isn't already in
 *	there.
 */
void
ListAddUnique( LIST lList, OBJEKT oObject )
{
    if ( bListContains( lList, oObject ) ) 
	return;
    ListAddToEnd( lList, oObject);
}





/*
 *	ListConcat
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Add the contents of lList2 to lList1.
 *	On returning, lList2 will be empty and should be destroyed.
 */
void
ListConcat( LIST lList1, LIST lList2 )
{
    VERIFYOBJEKT( lList1, LISTid );
    VERIFYOBJEKT( lList2, LISTid );

    MESSAGE("Before concat List1 size = %d,    List2 size = %d\n",
		iCollectionSize(lList1),
		iCollectionSize(lList2) );

    if ( iListSize(lList2) == 0 ) return;

	/* Make the LIST larger */

    CollectionSetSize( lList1, iCollectionSize(lList1)+
				iCollectionSize(lList2) );

	/* Point the last node of lList1 to the first node of lList2 */

    if ( lList1->nPLastNode != NULL ) {
	lList1->nPLastNode->nPNextNode = lList2->nPFirstNode;
	lList2->nPFirstNode->nPPrevNode = lList1->nPLastNode;
    } else {
	lList1->nPFirstNode = lList2->nPFirstNode;
    }

    lList1->nPLastNode = lList2->nPLastNode;

	/* Make List2 empty */

    lList2->nPFirstNode = NULL;
    lList2->nPLastNode = NULL;
    CollectionSetSize( lList2, 0 );

    MESSAGE("After concat List1 size = %d,    List2 size = %d\n",
		iCollectionSize(lList1),
		iCollectionSize(lList2) );

    
}




/*
 *      ListDescribe
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Fill the string with a description of the list.
 *      It is the callers responsibility to ensure that
 *      the string is long enough to store the description.
 */
void
ListDescribe( LIST lList )
{
OBJEKT          oObj;
int             iSize;
LISTLOOP        llL;

    iSize = iCollectionSize(lList);
    VP0("List size=%d\n", iSize );
    llL = (LISTLOOP)PCollectionLoop( (COLLECTION)lList );
    if ( llL == NULL ) 
	return;
    VP0("--list contents:\n" );
    while ( (oObj = oCollectionNext( (COLLECTION)lList, (GENP *)&llL )) 
								!= NULL ) {
        Describe( oObj );
    }
    VP0("--End of list\n" );
}




        
/*
 *      bListRemove
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Remove an object from the list.
 *      DEREF the object, possible Destroying it if there are
 *      no more references to it.
 *      
 *      Return:
 *              false if the element is not in the list.
 */
bool
bListRemove( LIST lList, OBJEKT oObject )
{
NODEP   nPNode;

                /* If the list is empty then return */

    if ( lList->nPFirstNode==NULL ) return false;
    if (bObjectInClass( oObject, CONTAINERid ) && CONTAINER_from(oObject)->nPListNode ) {
        nPNode = CONTAINER_from(oObject)->nPListNode;
    } else {
                /* Search the list for the object */
        nPNode = searchList( lList->nPFirstNode, oObject);
        if ( !nPNode ) return false;
    }

                /* Now remove the object from the list */
    if ( nPNode == lList->nPFirstNode ) {
        lList->nPFirstNode = nPNode->nPNextNode;
        if (nPNode->nPNextNode) nPNode->nPNextNode->nPPrevNode = NULL;
    } else nPNode->nPPrevNode->nPNextNode = nPNode->nPNextNode;
    if ( nPNode == lList->nPLastNode ) {
        lList->nPLastNode = nPNode->nPPrevNode;
        if (nPNode->nPPrevNode) nPNode->nPPrevNode->nPNextNode = NULL;
    } else nPNode->nPNextNode->nPPrevNode = nPNode->nPPrevNode;

                /* Destroy the node itself */
    FREE(nPNode);

    CollectionSetSize( lList, iCollectionSize(lList) - 1 );

                /* DEREF the object */
    DEREF( oObject );

    return true;
}





        
/*
 *      bListContains
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Find out if a list contains a particular object.
 *      
 *      Return:
 *              false if the element is not in the list.
 */
bool
bListContains( LIST lList, OBJEKT oObject )
{
                /* If the list is empty then return */
    if ( lList->nPFirstNode==NULL ) return false;

                /* Search the list for the object */
    return searchList( lList->nPFirstNode, oObject ) != NULL;
}




/*
 *      llListLoop
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Prepare the list to be looped over.
 *      Calls to ListNext will return the next object in the list.
 */
LISTLOOP
llListLoop( LIST lList )
{
    if ( lList == NULL )
	DFATAL("llListLoop called with NULL list\n" );
    return lList->nPFirstNode;
}




/*
 *      oListNext
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the next element in the list.
 *
 */
OBJEKT
oListNext( LISTLOOP *llPListLoop )
{
NODEP   nPNode;

    if ( *llPListLoop == NULL ) 
	return NULL;
    nPNode = *llPListLoop;
    *llPListLoop = (*llPListLoop)->nPNextNode;
    return nPNode->PObject;
}







/*
 *      lListDuplicate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return a duplicate of the list.
 *	Should only be called by the collectionduplicate
 *	routine (TODO - delete collection), which is
 *	called by oObjectDuplicate(),
 *	which sets objekt attributes.
 */
LIST
lListDuplicate( LIST lOld )
{
OBJEKT          lNew, oObj, oNew;
LISTLOOP        llLoop;

    lNew = oCreate(LISTid);
    llLoop = llListLoop(lOld);
    while ( ( oObj = oListNext(&llLoop) ) != NULL ) {
        oNew = oObjectDuplicate(oObj);
        /* Call ListAddTeEnd() to maintain the same order! JMK 2026 */
        if (GDefaults.bReverseLists) {
            ListAdd( (LIST)lNew, oNew );
            CONTAINER_from(oNew)->nPListNode = LIST_from(lNew)->nPFirstNode;
        } else {
            ListAddToEnd( (LIST)lNew, oNew );
            CONTAINER_from(oNew)->nPListNode = LIST_from(lNew)->nPLastNode;
        }
	oNew->iReferences = 1;	/* since ListAdd() increments */
    }
    return((LIST)lNew);
}
    
