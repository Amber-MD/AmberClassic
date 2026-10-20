/*
 *      File:   leap.h
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
 *      Contains all definitions that are used in LEaP
 */
 

#ifndef	LEAP_H
#define	LEAP_H

#ifndef	PARMLIB_H
# include        "parmLib.h"
#endif

#ifndef	AMBER_H
# include        "amber.h"
#endif




/* 
 *----------------------------------------------------
 *
 *      Global variables for LEaP
 */

#define LEAPRC  "leaprc"                /* Startup filename */

extern	bool	GbGraphicalEnvironment;	/* TRUE if in a graphical environment*/







/*
 *      variables.c
 */
extern void	VariablesInit(void);
extern void	VariablesList(void);
extern void	VariableSet(const char *sName, OBJEKT oObj);
extern void	VariableRemove(const char *sName );
extern OBJEKT	oVariable(const char *sName);
extern void	VariablesDestroy(void);
extern DICTIONARY dVariablesDictionary(void);



/*
 *      set.c
 */
extern  void    SetAttribute(void);         /* ( OBJEKT, STRING, OBJEKT ) */


#endif
