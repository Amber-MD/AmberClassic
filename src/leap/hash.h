/*
 *	File:	hash.h
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
 *		Handle hash tables.
 *
 */

#ifndef HASH_H
#define	HASH_H

typedef	struct	{
	char	*cPKey;
	GENP	PData;
} HASH_ENTRYt;


typedef	struct HASH_LISTs {
	HASH_ENTRYt		heEntry;
	struct HASH_LISTs	*hlPNext;
} HASH_LISTt;


typedef	struct	{
	HASH_LISTt	**hlPaTable;
	int		iTableSize;
	int		iElementCount;
} HASH_TABLEt;


typedef	HASH_TABLEt	*HASH_TABLE;


extern HASH_TABLE	htHTCreate(int iTableSize);
extern void		HTDestroy(HASH_TABLE *htPTable);
extern void		HTDescribe(HASH_TABLE htTable);

extern void		HTAdd(HASH_TABLE htTable, char *sKey, GENP PData); // sKey becomes owned by the hash!!
extern GENP		PHTFind(HASH_TABLE htTable, const char *sKey);
extern bool		bHTDelete(HASH_TABLE htTable, const char *sKey, 
				char **cPPKey, GENP *PPData);

extern void		HTWalk(HASH_TABLE htTable, GENP PCallback, 
				VFUNCTION VFunc);


#define	iHTTableSize(h)		((h)->iTableSize)
#define	iHTElementCount(h)	((h)->iElementCount)


#endif /* HASH_H */
