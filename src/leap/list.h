/*
 *      File: list.h
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
	OBJEKT                  PObject;
} NODE;

typedef NODE    *NODEP;
typedef NODE	*LISTLOOP;

typedef struct  {
	COLLECTIONt     oSuper;
	NODEP           nPFirstNode;
	NODEP           nPLastNode;
} LISTt;
                
typedef LISTt	*LIST;


#ifdef DEBUG
static inline LIST list_from_genp(GENP p)
  { assert(bIsObjekt(p)); assert(iObjectType((OBJEKT)p)==LISTid); return (LIST)p; }
static inline LIST list_from_objekt(OBJEKT o)
  { assert ( iObjectType(o)==LISTid ); return (LIST)o; }
static inline OBJEKT objekt_from_list(LIST l) { return &(l->oSuper.oObjectHeader); }

#define LIST_from(x) _Generic((x), \
    OBJEKT: list_from_objekt, \
    GENP: list_from_genp \
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
extern void		ListAddUnique(LIST lList, GENP PData);
extern void		ListConcat(LIST lList1, LIST lList2);
extern bool		bListRemove(LIST lList, GENP PPtr);
extern bool		bListContains(LIST lList, GENP PPtr);
extern LISTLOOP		llListLoop(LIST lList);
extern OBJEKT		oListNext(LISTLOOP *llPListLoop);

extern LIST		lListDuplicate(LIST lOld);

#define         iListSize(l)            iCollectionSize(l)

#endif /* LIST_H */
