#ifndef INTERNAL_H
#define INTERNAL_H
/*
 *      File: internal.h
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
 *      Class: 
 *              INTERNAL
 *      Superclass: 
 *              OBJEKT
 *
 *      Description:
 *
 *              An INTERNAL is an object representing an internal coordinate.
 *              Either a chirality, a bond length, a bond angle, 
 *              a torsion angle, or an improper torsion angle
 *              (not implemented yet).
 *              INTERNALs are contained within ATOMs.  Each atom
 *              that is referenced by the INTERNAL also contains the
 *              INTERNAL.  ATOMs contain internals as part of their
 *              CONTAINER contents.
 *		A special type of INTERNAL is the RING.  This allows
 *		the system to calculate and keep track of where RINGs
 *		are.
 *		
 *
 *
 */
 

#include	"ring.h"


/*
 *-----------------------------------------------------------------------
 *
 *      Define object typedefs here.
 *        
 *      Object typedef MUST include the superclass object as
 *      its first structure element.
*/


#define INTERNALUNDEFINED       '?'
#define INTERNALCHIRALITY       'C'
#define INTERNALBOND            'B'
#define INTERNALANGLE           'A'
#define INTERNALTORSION         'T'
#define INTERNALIMPROPER        'I'
#define	INTERNALRING		'R'

typedef	struct	{
	ATOM		aAtom1, aAtom2, aAtom3, aAtom4;
	double		dValue;
} NORMALINTERNALt;


/*
 *	New type used to store ring INTERNALs
 */


typedef	struct	{
	RING		rAtoms;		/* Store RING atoms here */
	RINGLOOP	rlLoop;		/* Used to loop over ATOMs */
	int		iRingTempInt;	/* Used to separate rings */
} RINGINTERNALt;

typedef struct INTERNALSTRUCT {
	OBJEKTt         oSuper;
	char            cInternalType;
	union	{
		NORMALINTERNALt		niNormal;
		RINGINTERNALt		riRing;
	} iType;
} INTERNALt;
                
typedef INTERNALt	*INTERNAL;


#define         MAXTORSIONSAROUNDBOND   12      /* Actually never more than 9 */


#ifdef DEBUG
static inline OBJEKT objekt_from_internal(INTERNAL i) { return &(i->oSuper); }
static inline INTERNAL internal_from_internal(INTERNAL i) { return i; }
static inline INTERNAL internal_from_objekt(OBJEKT o)
  { return o ? (assert(iObjectType(o)==INTERNALid),(INTERNAL)o) : NULL; }
#define INTERNAL_from(x) _Generic((x), \
    OBJEKT: internal_from_objekt, \
    INTERNAL: internal_from_internal \
)(x)
#else
#define INTERNAL_from(x) ((INTERNAL)(x))
#endif

/*
======================================================================

        Define object messages here.
        
        There must be at least a Create, Destroy, and Describe message.
        Hook into the messages of the superclasses so that
        when the message is sent to the most primative superclass
        of this class that it will eventually make it into these routines.
*/


/*      Define Create, Destroy, Describe methods */

extern INTERNAL		iInternalCreate(void);
extern void		InternalDestroy(INTERNAL *iPInt);
extern void		InternalDescribe(INTERNAL iInt);

extern INTERNAL		iInternalChirality(ATOM aAtom1, double dValue);
extern INTERNAL		iInternalBond(ATOM aAtom1, ATOM aAtom2, double dValue);
extern INTERNAL		iInternalAngle(ATOM aAtom1, ATOM aAtom2, ATOM aAtom3, 
				double dValue);
extern INTERNAL		iInternalTorsion(ATOM aAtom1, ATOM aAtom2, 
				ATOM aAtom3, ATOM aAtom4, double dValue);
extern INTERNAL		iInternalRing(void);

extern void		InternalRingAddAtomAfter(INTERNAL iRing, ATOM aAtom, 
				ATOM aPrev);
extern void		InternalRingAddAtomBefore(INTERNAL iInt, ATOM aAtom, 
				ATOM aBefore);
#define	iInternalRingSize(II) iRingSize((INTERNAL_from(II))->iType.riRing.rAtoms)
extern void		InternalRingLoopAtoms(INTERNAL iInt);
extern ATOM		aInternalRingNextAtom(INTERNAL iInt);
extern bool		bInternalRingRemoveAtom(INTERNAL iInt, ATOM aAtom);
#define	InternalRingSetTempInt(i,in)	\
		(INTERNAL_from(i)->iType.riRing.iTempInt = in )
#define	iInternalRingTempInt(i)	(INTERNAL_from(i)->iType.riRing.iTempInt)


extern INTERNAL		iInternalFindChirality(ATOM aAtom1);
extern INTERNAL		iInternalFindTorsion(ATOM aAtom1, ATOM aAtom2, 
				ATOM aAtom3, ATOM aAtom4);
extern INTERNAL		iInternalFindAngle(ATOM aAtom1, ATOM aAtom2, 
				ATOM aAtom3);
extern INTERNAL		iInternalFindBond(ATOM aAtom1, ATOM aAtom2);

extern bool		bInternalGoodTorsion(INTERNAL iInt, ATOM aAtom, 
				INTERNAL *iPTorsion, INTERNAL *iPAngle, 
				INTERNAL *iPBond);
extern bool		bInternalGoodAngle(INTERNAL iInt, ATOM aAtom, 
				INTERNAL *iPAngle, INTERNAL *iPBond);
extern bool		bInternalGoodPairOfAngles( INTERNAL iInt, ATOM aAtom,
				INTERNAL *iPAngle1, INTERNAL *iPAngle2, 
				INTERNAL *iPBond); 
extern bool		bInternalGoodBond(INTERNAL iInt, ATOM aAtom, 
				INTERNAL *iPBond);


#define InternalSetType(II,TT)  (INTERNAL_from(II)->cInternalType=TT)
#define iInternalType(II)       ((int)(INTERNAL_from(II)->cInternalType))
#define aInternalAtom1(II)      (INTERNAL_from(II)->iType.niNormal.aAtom1)
#define aInternalAtom2(II)      (INTERNAL_from(II)->iType.niNormal.aAtom2)
#define aInternalAtom3(II)      (INTERNAL_from(II)->iType.niNormal.aAtom3)
#define aInternalAtom4(II)      (INTERNAL_from(II)->iType.niNormal.aAtom4)
#define dInternalValue(II)      (INTERNAL_from(II)->iType.niNormal.dValue)
#define InternalSetValue(II,VV) (INTERNAL_from(II)->iType.niNormal.dValue=VV)


extern int	iInternalFindAllTorsionInternalsAround(ATOM aAtom2, 
			ATOM aAtom3, INTERNAL iaTorsions[]);

#endif/* INTERNAL_H */
