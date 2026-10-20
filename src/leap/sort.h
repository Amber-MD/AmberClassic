/*
 *      File: sort.h
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
 *      Description:
 *              This is an interface to the general sorting routine
 *              'Sort'.  It allows for sorting by string, integer, and
 *              double keys.
 */

#ifndef	SORT_H
#define	SORT_H

#define	SORT_ASCENDING		TRUE
#define	SORT_DESCENDING		FALSE

extern void	SortByInteger(GENP PStart, int iElements, int iSize, 
			GENP PFirst, bool bAscending);
extern void	SortByDouble( GENP PStart, int iElements, int iSize, 
			GENP PFirst, bool bAscending);
extern void	SortByString(GENP PStart, int iElements, int iSize,
			GENP PFirst, bool bAscending);

typedef bool	(*SIFTFUNCTION)(void*); 

extern void	Sift(GENP PData, int iElementSize, int iElements, 
			SIFTFUNCTION bFCriteria, int *iPFirstFalse);
#endif
