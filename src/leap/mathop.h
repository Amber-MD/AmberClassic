/*
 *      File:   mathop.h
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
 *              This file tries to contain many of the mathematical
 *              operations that are used in molecular modeling.
 *              Eg: converting nonbond parameters for two individual
 *                      atoms into a single A,C pair. etc.
 */



#ifndef	MATHOP_H
#define	MATHOP_H

extern void	MathOpConvertNonBondToAC(double dE1, double dR1, 
			double dE2, double dR2, double *dPA, double *dPC);
extern void	MathOpDiagonalize( MATRIX mA, VECTOR *vEigen, MATRIX mV, 
			int *iPNrot );
extern void	MathOpMomentOfInertia( UNIT uUnit, MATRIX mMoment );


#endif
