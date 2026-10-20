/*
 *      File: oInteger.h
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
 *              OINTEGER        
 *      Superclass: 
 *              OBJEKT
 *
 *      Description:
 *
 *              OINTEGER is an object that has the same properties
 *              as normal integers.  It is only declared to
 *              make passing arguments to functions as homogenous as
 *              possible.  It should not be used for any other purpose. 
 *
 *
 */
 
#ifndef	OINTEGER_H
#define OINTEGER_H


/*
-----------------------------------------------------------------------

        Define object typedefs here.
        
        Object typedef MUST include the superclass object as
        its first structure element.
*/


typedef struct  {
	OBJEKTt         oSuper;
	int             iValue;
} OINTEGERt;

typedef OINTEGERt	*OINTEGER;



/*
======================================================================

        Define object messages here.
        
        There must be at least a Create, Destroy, and Describe message.
        Hook into the messages of the superclasses so that
        when the message is sent to the most primative superclass
        of this class that it will eventually make it into these routines.
*/


/*      Define Create, Destroy, Describe methods */


extern OINTEGER		oiOIntegerCreate(void);
extern void		OIntegerDestroy(OINTEGER *oiPInt);
extern void		OIntegerDescribe(OINTEGER oiInt);
extern OINTEGER		oiOIntegerDuplicate( OINTEGER oi );

#define iOInteger( OI ) ( ((OINTEGER)OI)->iValue )
#define OIntegerSet( OI, V ) ( ((OINTEGER)OI)->iValue = V )



#endif  /* OINTEGER_H */
