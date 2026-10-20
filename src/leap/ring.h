/*
 *	File:	ring.h
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
 *		Maintain a RING structure.  A doubly linked list of
 *		nodes where insertions can be made relative to other
 *		elements within the RING.  Each node of the RING
 *		has a pointer to an arbitrary data type.
 */


#ifndef	RING_H
#define	RING_H



typedef	struct	RINGNODEs {
	GENP			PData;
	struct RINGNODEs	*rnPPrev;
	struct RINGNODEs	*rnPNext;
} RINGNODEt;

typedef	struct	{
	RINGNODEt		*rnPFirst;
	int			iElements;
} RINGt;

typedef	struct	{
	RINGNODEt	*rnPNext;
	RINGNODEt	*rnPCur;
	int		iCount;
} RINGLOOP;

typedef	RINGt *		RING;



extern RING	rRingCreate(void);
extern void	RingDestroy(RING *rPRing);

#define	bRingEmpty(r)	(((RING)(r))->rnPFirst == NULL )
#define	iRingSize(r)	(((RING)(r))->iElements)

extern void	RingAfterAdd(RING rRing, GENP PAfter, GENP PData);
extern void	RingBeforeAdd(RING rRing, GENP PBefore, GENP PData);

extern void	rlRingLoop(RING rRing, RINGLOOP *rlP);
extern GENP	PRingNext(RINGLOOP *rlPLoop);
extern bool	bRingRemove(RING rRing, GENP PData);

#define	rnPRingNode(rl)		(rl.rnPCur)




#endif	/* RING_H */














