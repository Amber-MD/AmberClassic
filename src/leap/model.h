/*
 *      File:   model.h
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
 *              The MODEL object is an object with a single instance
 *              that interacts with the BUILD object and generates
 *              model bond lengths, angles, torsions etc which
 *              are determined from atom types.
 */


#ifndef	MODEL_H
#define	MODEL_H

extern double	dModelBondLength(ATOM aAtom1, ATOM aAtom2);
extern double	dModelBondAngle(ATOM aAtom1, ATOM aAtom2, ATOM aAtom3 );
extern void	ModelAddHydrogens( UNIT uUnit );
extern void	ModelProperTerms(ATOM a1, ATOM a2, ATOM a3, ATOM a4, TORSION tTorsion);
extern void	ModelAngleParm(ATOM a1, ATOM a2, ATOM a3, double *dPK, double *dPE);
extern void	ModelBondParm(ATOM a1, ATOM a2, double *dPK, double *dPE);
extern void	ModelAssignTorsionsAround(ATOM aX, ATOM aY, ATOM aOnly);

#endif
