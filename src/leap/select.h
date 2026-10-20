/*
 *	File:	select.h
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



