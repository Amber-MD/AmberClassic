/*
 *      File: bag.h
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
 *      Define BAG class.
 *	A BAG maintains a collection of objects in a linked
 *	list.  The objects are arbitrary pointers.
 *	This code was taken from list.c/list.h and
 *	modified to handle general objects.
 *
 */

#ifndef	BAG_H
#define BAG_H 

typedef	struct	BNodeStruct	{
	struct BNodeStruct	*bnPNextNode;
	GENP			PObject;
} BNODE;

typedef BNODE    *BNODEP;

typedef	struct	{
	BNODE	*bnPCur;
	BNODE	*bnPNext;
} BAGLOOP;

typedef struct  {
	int		iElementCount;
	BNODEP		bnPFirstNode;
	BNODEP		bnPLastNode;
} BAGt;
                
typedef BAGt	*BAG;


#define		iBagSize(b)            ( b->iElementCount )


/*  bag.c  */

extern BAG		bBagCreate(void);
extern void 		BagDestroy( BAG *bPBag );
extern void		BagAdd( BAG bBag, GENP PObj );
extern void		BagAddToEnd( BAG bBag, GENP PObj );
extern void		BagDescribe( BAG bBag );
extern bool		bBagRemove( BAG bBag, GENP PPtr );
extern bool		bBagContains( BAG bBag, GENP PPtr );
extern void		BagUniqueAdd( BAG bBag, GENP PObj );
extern BAGLOOP		blBagLoop( BAG bBag );
extern GENP		PBagNext( BAGLOOP *blPBagLoop );


#endif	/* BAG_H */

