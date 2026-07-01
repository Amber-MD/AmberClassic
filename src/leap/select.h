/*
 *	File:	select.h
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
 *		Select different parts of a UNIT.
 */


#include "minimizer.h"

extern void	SelectAtom(ATOM aAtom, bool bOn);
extern bool	bSelectRingWithAtom(UNIT uUnit, ATOM aAtom, bool bOn);
extern void	SelectResidueWithAtom(UNIT uUnit, ATOM aAtom, bool bOn);
extern void	SelectMoleculeWithAtom(UNIT uUnit, ATOM aAtom, bool bOn);
extern void	SelectEverything(UNIT uUnit, bool bOn);
extern void	SelectChainBetween(UNIT uUnit, ATOM aA, ATOM aB, bool bOn);
extern void	SelectRelaxInFramework(UNIT uUnit, MINIMIZER mMinimizer);
extern bool     bSelectChainBetween( UNIT uUnit, ATOM aA, ATOM aB, bool bOn );

#define	SelectSetFlags(u,f)	ContainerWithFlagsSetAtomFlags(u,ATOMSELECTED,f)
#define	SelectResetFlags(u,f)\
			ContainerWithFlagsResetAtomFlags(u,ATOMSELECTED,f)



