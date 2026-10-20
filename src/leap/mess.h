/*
 *	File:	message.h
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
 *		Store messages in a long list.
 *		Keep track of what function they come from
 *		using 'function.h'.
 */


#ifndef	BASICS_H
#include	"basics.h"
#endif

typedef	struct MESS_STRUCT	{
	struct MESS_STRUCT	*mNext;
	char			cType;
	int			iFunction;
	char			*cPText;
	bool			bPrint;
} MESSt;

typedef	MESSt	*MESS;


extern MESS	GmMesss;
extern MESS	GmLastMess;

extern MESS	mMessAdd(char *sText);
extern MESS	mMessLoop(void);
extern MESS	mMessNext(MESS *mPMess);

#define	sMessText(m)		(m->cPText)
#define	bMessPrint(m)		(m->bPrint)
#define	iMessFunction(m)	(m->iFunction)
#define	MessSetPrint(m,b)	(m->bPrint = b)

