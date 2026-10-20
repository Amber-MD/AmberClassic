/*
 *	File:	graphUtil.h
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
 *		This file contains functions that perform
 *		various operations on graphs.  Some of the
 *		functions are finding the shortest paths between
 *		two vertices, finding all of the smallest cycles,
 *		so on and so on.
 */





/*
 *	Functions
 */

extern bool	bGraphUtilFindShortestPath( UNIT uUnit, ATOM aStart, ATOM aStop,
                	void (*fPCallback)());
extern void	zGraphUtilDescribeRing( INTERNAL iRing );
extern bool	zbGraphUtilSeparate( INTERNAL iBig, int iBigIndex, 
			INTERNAL iSmall, int iSmallIndex );
extern int	ziGraphUtilJoinRingGroups( VARARRAY vaRingGroup, int iGroup1, 				int iGroup2 );
extern void	GraphUtilFindAllSmallestRingsAndRingGroups( UNIT uUnit, 
			VARARRAY *vaPRingGroups);
extern void	GraphUtilDestroyRingGroupVarArray( VARARRAY *vaPRingGroup );
extern void	GraphUtilFindAllSmallestRings( UNIT uUnit );



