/*
 *      File: list.h
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
 *      Define LIST class.
 *
 */

#ifndef LIST_H
#define LIST_H

#include "basics.h"
#include "collection.h"
#include "objekt.h"

typedef struct  NodeStruct      {
	struct NodeStruct       *nPNextNode;
	struct NodeStruct       *nPPrevNode;
	OBJEKT                  PObject;
} NODE;

typedef NODE    *NODEP;
typedef NODE	*LISTLOOP;

typedef struct LISTSTRUCT {
	COLLECTIONt     oSuper;
	NODEP           nPFirstNode;
	NODEP           nPLastNode;
        bool            bFreeChildren;
} LISTt;
                
typedef LISTt	*LIST;


#ifdef DEBUG
static inline LIST list_from_genp(GENP p)
  { assert(bIsObjekt(p)); assert(iObjectType((OBJEKT)p)==LISTid); return (LIST)p; }
static inline LIST list_from_objekt(OBJEKT o)
  { assert ( iObjectType(o)==LISTid ); return (LIST)o; }
static inline OBJEKT objekt_from_list(LIST l) { return &(l->oSuper.oObjectHeader); }
static inline LIST list_from_list(LIST l) { return l; }

#define LIST_from(x) _Generic((x), \
    OBJEKT: list_from_objekt, \
    GENP: list_from_genp, \
    LIST: list_from_list \
)(x)
#else
#define LIST_from(x) ((LIST)(x))
#endif

/*
--------------------------------------------------------------------

        LIST messages
        
*/

extern LIST		lListCreate(void);
extern void		ListDestroy(LIST *lPList);
extern void		ListDescribe(LIST lList);

extern void		ListAdd(LIST lList, OBJEKT oObj);
extern void		ListAddToEnd(LIST lList, OBJEKT oObj);
extern void		ListAddUnique(LIST lList, OBJEKT oObj);
extern void		ListConcat(LIST lList1, LIST lList2);
extern bool		bListRemove(LIST lList, OBJEKT oObj);
extern bool		bListContains(LIST lList, OBJEKT oObj);
extern LISTLOOP		llListLoop(LIST lList);
extern OBJEKT		oListNext(LISTLOOP *llPListLoop);

extern LIST		lListDuplicate(LIST lOld);

#define         iListSize(l)            iCollectionSize(l)

#endif /* LIST_H */
