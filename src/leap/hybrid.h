/*
 *	File:		hybrid.h
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
 *	Author:		David Rivkin
 *	Date Created:	17 August 1992
 *	Dates Changed:	
 *	Description:
 *		Header file for hybrid.c functions.
 *		These functions translate between integer and string 
 *		representations of the Hybridization methods.
 *
 */
 
 
#define HYBRIDNAMELEN	10	/* "undefined" is the longest possible string */

extern int	iHybridNumber(char *sName);
extern char	*sHybridName(int iHybrid, char *sName);
