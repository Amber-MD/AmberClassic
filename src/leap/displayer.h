/*
 *	File:	displayer.h
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
 *		Manage an object that calls callback routines
 *		to update objects that are being displayed in
 *		a graphical user interface.
 */

#ifndef DISPLAYER_H
#define DISPLAYER_H

#define	D_MODIFIED	0x01
#define	D_SENSITIVE	0x02

typedef	struct	DISPLAYERNODEs {
	struct DISPLAYERNODEs	*dnNext;
	VFUNCTION		fCallback;
	GENP			PData;
} DISPLAYERNODEt;

typedef	DISPLAYERNODEt	*DISPLAYERNODE;

typedef	struct	DISPLAYERs {
	char			cStatus;
	DISPLAYERNODE		dnFirst;
	GENP			PObject;
	struct DISPLAYERs	*dNext;
} DISPLAYERt;

typedef	DISPLAYERt	*DISPLAYER;

extern	int		GiDisplayerAccumulateUpdates;


extern DISPLAYER	dDisplayerCreate(GENP PObject);
extern void		DisplayerDestroy(DISPLAYER *dPOld);

extern void		DisplayerAdd(DISPLAYER dDisp, VFUNCTION vFunc, 
				GENP PData );
extern bool		bDisplayerRemove(DISPLAYER dDisp, VFUNCTION vFunc, 
				GENP PData);
extern void		DisplayerUpdate(DISPLAYER dDisp);

#define	DisplayerAccumulateUpdates()	(GiDisplayerAccumulateUpdates++)
extern void		DisplayerReleaseUpdates(void);
extern void		TurnOffDisplayerUpdates(void);
extern void		TurnOnDisplayerUpdates(void);

#define	DisplayerSetSensitive(d,s) \
(d->cStatus = (s ? (d->cStatus|D_SENSITIVE) : (d->cStatus&(~D_SENSITIVE))) )

#define	bDisplayerSensitive(d)	(d->cStatus&D_SENSITIVE)

#define	DisplayerSetModified(d)	(d->cStatus|=D_MODIFIED)
#define	bDisplayerModified(d)	(d->cStatus&D_MODIFIED)


#endif /* DISPLAYER_H */
