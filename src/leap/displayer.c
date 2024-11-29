/*
 *	File:	displayer.c
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
 *		Manage an object that calls callback routines
 *		to update objects that are being displayed in
 *		a graphical user interface.
 *
 *		The way to use these things is within any object
 *		that can be displayed within a window, add a DISPLAYER
 *		as one of its instance variables.
 *
 *		When you open a window that contains the object
 *		call DisplayerAdd and add the callback
 *		routine and an argument that will identify the
 *		object to the calback; the callback is one
 *		that will update the window.  Call DisplayerRemove
 *		when the window is closed.
 *
 *		Then whenever the object has been modified
 *		call DisplayerUpdate. 
 *	
 *		A performance enhancement that can be used to 
 *		buffer DISPLAYER updates until a large number
 *		of updates has been done is the concept of
 *		sensitivity.  DISPLAYERs can be made inSENSITIVE
 *		which will cause them to flag that they should
 *		be updated when they become sensitive again.
 */


#include	"basics.h"
#include	"displayer.h"

extern BOOL	GbGraphicalEnvironment;

extern int	GiUnitEditors;

/*
 *----------------------------------------------------------
 *
 *	Private variables.
 */
 
int		GiDisplayerAccumulateUpdates = 0;
static	DISPLAYER	SdDisplayerAccumulateList = NULL;



/*
 *	dDisplayerCreate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Create a DISPLAYER for the object in PObject.
 */
DISPLAYER
dDisplayerCreate( GENP PObject )
{
DISPLAYER	dNew;

	return(NULL);

}




/*
 *	DisplayerDestroy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Destroy the DISPLAYER and the list of nodes.
 */
void	
DisplayerDestroy( DISPLAYER *dPOld )
{
DISPLAYERNODE	dnCur, dnNext;
DISPLAYER	dPrev, dCur;

	return;

}




/*
 *	DisplayerAdd
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Add a callback routine to the DISPLAYER.
 */
void
DisplayerAdd( DISPLAYER dDisp, VFUNCTION vFunc, GENP PData )
{
}



/*
 *	bDisplayerRemove
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Remove a callback function from the DISPLAYER list.
 *	Return TRUE if it was found.
 */
BOOL
bDisplayerRemove( DISPLAYER dDisp, VFUNCTION vFunc, GENP PData )
{
   return(TRUE);
}

/*
 *	DisplayerUpdate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Loop through the list of callbacks, calling
 *	each one in turn with the object and the argument.
 *	Then clear the D_MODIFIED property of the DISPLAYER.
 *	
 *	If the DISPLAYER is insensitive then set the MODIFIED
 *	flag of the DISPLAYER and return.
 */
void
DisplayerUpdate( DISPLAYER dDisp )
{

	return;

}


/*
 *	DisplayerReleaseUpdates
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Decrement the GiDisplayerAccumulateUpdates variable, if it is
 *	zero then:
 *	Update all DISPLAYERs in the SdDisplayerAccumulateList.
 */
void
DisplayerReleaseUpdates()
{
}

/*
 *  TurnOffDisplayerUpdates(), TurnOnDisplayerUpdates()
 *
 *	Recommended that these only be called from commands.c
 *	and that ContainerDisplayerUpdate( <unit> ) be called
 *	right after. 
 */
static int	iOffed = 0;

void
TurnOffDisplayerUpdates()
{
}

void
TurnOnDisplayerUpdates()
{
}
