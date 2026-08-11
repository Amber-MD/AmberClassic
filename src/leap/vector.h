/*
 *	File:	vector.h
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
 *             David A. Rivkin                                             *
 *                                                                      *
 *     Principal Investigator: Peter A. Kollman                         *
 *                                                                      *
 ************************************************************************
 *
 *	Description:
 *		Manage three dimensional double precision vectors.
 */

#ifndef VECTOR_H
#define VECTOR_H

typedef struct  {
	double  dX, dY, dZ;
} VECTOR;


/*
 *      vector.c messages
 */

#define VectorDef( vP, x, y, z )	( (vP)->dX=x, (vP)->dY=y, (vP)->dZ=z )


#define dVX(v)                 ((v)->dX)
#define dVY(v)                 ((v)->dY)
#define dVZ(v)                 ((v)->dZ)

#define VectorSetX(v,d)       ((v)->dX = d)
#define VectorSetY(v,d)       ((v)->dY = d)
#define VectorSetZ(v,d)       ((v)->dZ = d)


/*  vector.c  */

extern VECTOR		vVectorAdd( VECTOR *vPX, VECTOR *vPY );
extern VECTOR		vVectorSub( VECTOR *vPX, VECTOR *vPY );
extern double		dVectorLen( VECTOR *vPX );
extern double		dVectorAngle( VECTOR *vPX, VECTOR *vPY );
extern double		dVectorAbsAngle( VECTOR *vPX, VECTOR *vPY, 
				VECTOR *vPRef );
extern VECTOR		vVectorTimesScalar( VECTOR *vPX, double dS );
extern double		dVectorDot( VECTOR *vPX, VECTOR *vPY );
extern VECTOR		vVectorCross( VECTOR *vPX, VECTOR *vPY );
extern double		dVectorAtomChirality( VECTOR *vPCenter, 
				VECTOR *vPA, VECTOR *vPB, VECTOR *vPC );
extern double		dVectorAtomNormalizedChirality( VECTOR *vPCenter, 
				int iCenterCoord,
				VECTOR *vPA, bool bA, VECTOR *vPB, bool bB,
				VECTOR *vPC, bool bC, VECTOR *vPD, bool bD );
extern double		dVectorAtomLength( VECTOR *vPA, VECTOR *vPB );
extern double		dVectorAtomAngle( VECTOR *vPA, VECTOR *vPB, 
				VECTOR *vPC );
extern double		dVectorAtomTorsion( VECTOR *vPA, VECTOR *vPB, 
				VECTOR *vPC, VECTOR *vPD );
extern double           dVectorDistance( VECTOR *Pv1, VECTOR *Pv2 );
extern double           dVectorDistanceSq( VECTOR *Pv1, VECTOR *Pv2 );

#endif          /* ifndef VECTOR_H */
