/*
 *      File:   nVector.h
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
 *              Routines to maintain N dimensional vectors
 */

#ifndef NVECTOR_H

#define NVECTOR_H

typedef struct  {
	int     iSize;
	double  daValues[1];
} NVECTORt;

typedef NVECTORt	*NVECTOR;
                
                
/*
        nVector.c messages
*/

extern NVECTOR	nvNVectorCreate(int iElements);
extern void	NVectorDestroy(NVECTOR *nvPX);
extern void	NVectorDescribe(NVECTOR nvVector);

extern NVECTOR	nvNVectorCopy(NVECTOR nvVector);
extern void	NVectorCopy(NVECTOR nvDest, NVECTOR nvSource);
extern void	NVectorZero(NVECTOR nvVector);
extern void	NVectorFill(NVECTOR nvVector, double dValue);
extern void	NVectorAdd(NVECTOR nvA, NVECTOR nvX, NVECTOR nvY);
extern void	NVectorSub(NVECTOR nvA, NVECTOR nvX, NVECTOR nvY);
extern void	NVectorTimesScalar(NVECTOR nvA, NVECTOR nvX, double dS);
extern double	dNVectorDot(NVECTOR nvX, NVECTOR nvY);
extern double	dNVectorLen(NVECTOR nvX);
#define iNVectorSize(N)                 (N->iSize)
#define dNVectorElement(N,I)            (((N->daValues))[I] )
#define NVectorSetElement(N,I,E)        (((N->daValues))[I] = E)


#endif          /* ifndef NVECTOR_H */
