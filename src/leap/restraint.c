/*
 *	File:	restraint.c
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
 *		Maintain RESTRAINTs.
 */



#include	"basics.h"

#include	"classes.h"
#include	"restraint.h"



/*
 *	rRestraintCreate
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return a pointer to a RESTRAINT of undefined type.
 */
RESTRAINT
rRestraintCreate()
{
RESTRAINT	rNew;

    rNew = (RESTRAINT)MALLOC(sizeof(RESTRAINTt) );
    rNew->iType = RESTRAINTNONE;
    RestraintDefineFlags(rNew, RESTRAINTALWAYS );
#ifdef DEBUG
    rNew->iMagic = RESTRAINT_MAGIC;
#endif
    return rNew;
}




/*
 *	RestraintDestroy
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	FREE the memory allocated to the RESTRAINT.
 */
void
RestraintDestroy( RESTRAINT *rPRes )
{
    FREE( *rPRes );
    *rPRes = NULL;
}





/*
 *	bRestraintBondMatchAtoms
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return true if the RESTRAINT bond has the
 *	same atoms as the user passed in.
 */
bool
bRestraintBondMatchAtoms( RESTRAINT rRes, ATOM aAtom1, ATOM aAtom2 )
{
    if ( (rRes->rType.rbBond.aAtom1 == aAtom1 &&
	  rRes->rType.rbBond.aAtom2 == aAtom2 ) ||
	 (rRes->rType.rbBond.aAtom1 == aAtom2 &&
	  rRes->rType.rbBond.aAtom2 == aAtom1 ) ) return true;
    return false;
}

	
/*
 *	bRestraintAngleMatchAtoms
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return true if the RESTRAINT angle has the
 *	same atoms as the user passed in.
 */
bool
bRestraintAngleMatchAtoms( RESTRAINT rRes, 
	ATOM aAtom1, ATOM aAtom2, ATOM aAtom3 )
{

    if ( (rRes->rType.raAngle.aAtom1 == aAtom1 &&
          rRes->rType.raAngle.aAtom2 == aAtom2 &&
	  rRes->rType.raAngle.aAtom3 == aAtom3 ) ||
	 (rRes->rType.raAngle.aAtom1 == aAtom3 &&
	  rRes->rType.raAngle.aAtom2 == aAtom2 &&
	  rRes->rType.raAngle.aAtom3 == aAtom1 ) ) return true;
    return false;
}



/*	
 *	bRestraintTorsionMatchAtoms
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Return true if the RESTRAINT torsion has the
 *	same atoms as the user passed in.
 */
bool
bRestraintTorsionMatchAtoms( RESTRAINT rRes, 
	ATOM aAtom1, ATOM aAtom2, ATOM aAtom3, ATOM aAtom4 )
{
    if ( (rRes->rType.rtTorsion.aAtom1 == aAtom1 &&
          rRes->rType.rtTorsion.aAtom2 == aAtom2 &&
	  rRes->rType.rtTorsion.aAtom3 == aAtom3 &&
	  rRes->rType.rtTorsion.aAtom4 == aAtom4 ) ||
	 (rRes->rType.rtTorsion.aAtom1 == aAtom4 &&
	  rRes->rType.rtTorsion.aAtom2 == aAtom3 &&
	  rRes->rType.rtTorsion.aAtom3 == aAtom2 &&
	  rRes->rType.rtTorsion.aAtom4 == aAtom1 ) ) return true;
    return false;
}



/*
 *	bRestraintContainsAtom
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *	Test to see if the atom is one of the atoms of
 *	the RESTRAINT.
 *	If it is return true.
 */
bool
bRestraintContainsAtom( RESTRAINT rRes, ATOM aAtom )
{
    switch ( rRes->iType ) {
	case RESTRAINTTORSION:
	    if ( rRes->rType.rtTorsion.aAtom1 == aAtom ||
		 rRes->rType.rtTorsion.aAtom2 == aAtom ||
		 rRes->rType.rtTorsion.aAtom3 == aAtom ||
		 rRes->rType.rtTorsion.aAtom4 == aAtom ) return true;
	    return false;
	case RESTRAINTANGLE:
	    if ( rRes->rType.raAngle.aAtom1 == aAtom ||
		 rRes->rType.raAngle.aAtom2 == aAtom ||
		 rRes->rType.raAngle.aAtom3 == aAtom ) return true;
	    return false;
	case RESTRAINTBOND:
	    if ( rRes->rType.rbBond.aAtom1 == aAtom ||
		 rRes->rType.rbBond.aAtom2 == aAtom ) return true;
	    return false;
	default:
	    DFATAL("Invalid RESTRAINT type!" );
    }
    return false;	/* for lint */
}



