/*
 *      File: object.h
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
 *      Class: 
 *              OBJEKT
 *      Superclass:
 *              none.   (Class heirarchy is defined in object.c)
 *
 *      Description:
 *              Defines ALL classes, the OBJEKT class and
 *              OBJEKT messages.
 */

#ifndef OBJEKT_H
#define OBJEKT_H
#include <stdint.h> 
/*
#####################################################################

        Object class types

        Every Object subclass of OBJEKT must have a define here.
        So that the OBJEKT messages can recognize where
        to call to handle subclass messages.
        Suffix the macro with "id".
*/

#define NULLid          'N'
#define OBJEKTid        'O'
#define BYTEARRAYid     'B'
#define COLLECTIONid    'c'
#define LISTid          'L'
#define CONTAINERid     'C'
#define UNITid          'U'
#define MOLECULEid      'M'
#define RESIDUEid       'R'
#define ATOMid          'A'
#define INTERNALid      'I'
#define PARMSETid       'P'
#define OINTEGERid      'i'
#define ODOUBLEid       'd'
#define OSTRINGid       's'
#define ASSOCid         'a'

#define	NO_OBJEKTid	'?'


/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        If you are adding new CLASSES . . .
        DO NOT MODIFY BELOW HERE! ----------------------------------
        DO NOT MODIFY BELOW HERE! ----------------------------------
        DO NOT MODIFY BELOW HERE! ----------------------------------
        DO NOT MODIFY BELOW HERE! ----------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

/*
=====================================================================

        OBJEKT instance variables.
        
        Either this, or the instance variables of a subclass of
        OBJEKT MUST!!!! be included as the first thing in a subclass'
        instance variable typedef.
        
        Suffix each type with "t".

        OBJEKTt contains a field that describes the object type and
        a field that is used for counting the references to the object
        for garbage collection.  There is also a linked list used by
        the garbage collector.
        OBJEKTs that are garbage collected have an iReferences field
        of 0 or greator.  An iReferences value of 0 causes the OBJEKT
        to be Destroyed and storage reclaimed.

        Reference counting is done using two macros:
                REF(x)  which increments the reference count.
                DEREF(x) which decrements the reference count and
                        Destroys the OBJEKT if necessary.
        REF and DEREF need the OBJEKTs, NOT POINTERS to OBJEKTs
        It is the PROGRAMMERS responsibility to make sure that
        OBJEKTs references are properly maintained.
        The action of Creating an OBJEKT sets its reference count
        to 1.

        Garbage collection is mainly done in the command line interface.

*/

// gcc object poiner tracker! TODO
//#if defined(__GNUC__) && defined(__analyzer__)
//    // Explicitly pair these functions together for the analyzer
//    __attribute__((malloc(decref, 1))) Object* create_object(int initial_value);
//    __attribute__((malloc(decref, 1))) Object* addref(Object* obj);
//#endif

typedef struct OBJEKTSTRUCT {
	char            cObjType;
#ifdef DEBUG
	GENP            PKlass;
#endif
	int             iReferences; 
} OBJEKTt;

typedef OBJEKTt *OBJEKT;

extern	KLASSt	GkTObjekt;

#ifdef DEBUG
#define VERIFYOBJEKT_loc( O, Oid, file, line ) {\
    if ( iObjectType(O) != Oid ) {\
        fprintf(stderr, "ERROR, Object of type: %c should be type: %c\n",\
                iObjectType(O), Oid );\
        fprintf(stderr, "In file: %s   line: %d\n",file,line );\
        abort(); \
    }}
#define CK_loc(x,c,file,line) ((!(x)||(GENP)((x)->PKlass)!=(GENP)(c))?BasicsKlassMismatchPanic((x),file,line):0)
#define assert_loc(cond, file, line) \
    ((cond) ? (void)0 : (fprintf(stderr, "%s:%d: assertion failed: %s\n", (file), (line), #cond), abort()))
static inline bool bIsObjekt(void *p) {
    OBJEKT o = (OBJEKT)p;
    return (o->PKlass == (GENP)&GkTObjekt);
}
static inline uint32_t iGetObjektType( void *o ) {
    if (!o) return NULLid;
    assert(bIsObjekt(o));
    return ((OBJEKT)o)->cObjType;
}
static inline OBJEKT objekt_from_objekt(OBJEKT o) { return o; }
static inline OBJEKT objekt_from_genp(GENP p) { return p ? (assert(bIsObjekt(p)), (OBJEKT)p) : NULL; }

#define OBJEKT_from(x) _Generic((x), \
    OBJEKT: objekt_from_objekt, \
    CONTAINER: objekt_from_container, \
    UNIT: objekt_from_unit, \
    MOLECULE: objekt_from_molecule, \
    RESIDUE: objekt_from_residue, \
    ATOM: objekt_from_atom, \
    INTERNAL: objekt_from_internal, \
    GENP: objekt_from_genp, \
    LIST: objekt_from_list, \
    ASSOC: objekt_from_assoc, \
    PARMSET: objekt_from_parmset, \
    OSTRING: objekt_from_string, \
    ODOUBLE: objekt_from_double \
)(x)

#define iObjectType( o ) iGetObjektType((void*)o)

#else

#define OBJEKT_from(x) ((OBJEKT)(x))
#define iObjectType( o ) \
        ( (o) ? (int)(OBJEKT_from(o)->cObjType) : NULLid )

#endif


/*
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        OBJEKT messages.
        
        New classes being defined will most likely want to to
        modify the code for the following routines in 'object.c'
        There are two Create methods, one is for normal objects,
        the other is for objects with indices, eg: ARRAYs.
        
*/

extern OBJEKT	oCreateSize(int iType, int iSize);
extern OBJEKT	oCreate(int iType);
extern OBJEKT	oCopy(OBJEKT oCur);
extern void	Destroy(OBJEKT  *oPObject);
extern void	Describe(OBJEKT oObject);

extern bool	bObjectInClass(OBJEKT oObject, int iClass);
extern char	*sObjectIndexType(int iType);
extern char	*sObjectType(OBJEKT oObj);
extern bool	bObjektWarnType(OBJEKT oObj, int iType);


	/* DO NOT USE oObjectDuplicate from applications */
	/* INSTEAD USE oCopy */
	/* oCopy calls oObjectDuplicate and then cleans up internal */
	/* references between objects */

extern OBJEKT	oObjectDuplicate(OBJEKT oOld);

#define	INVALIDATE_OBJEKT(o)		(OBJEKT_from(o)->cObjType = NO_OBJEKTid )



#define REF( o ) (!(o)?NULL:(OBJEKT_from(o)->iReferences++,(o)))

#define DEREF( o )      {\
	if ( (OBJEKT)(o)!=NULL ) {\
        	OBJEKT_from(o)->iReferences--;\
        	if ( OBJEKT_from(o)->iReferences<=0 ) \
			Destroy((OBJEKT *)&(o));\
	}\
}


/*
####################################################################

        Common routines for debuggin objects.
*/

#ifdef  DEBUG
#define VERIFYOBJEKT( O, Oid ) {\
    if ( iObjectType(O) != Oid ) {\
        MESSAGE( "ERROR, Object of type: %c should be type: %c\n",\
                iObjectType(O), Oid );\
        MESSAGE( "In file: %s   line: %d\n", __FILE__, __LINE__ );\
    }}
#define BADOBJEKT( O ) {\
    if ( iObjectType(O) != Oid ) {\
        MESSAGE( "ERROR, Bad object of type: %c\n",\
                iObjectType(O) );\
        MESSAGE( "In file: %s   line: %d\n", __FILE__, __LINE__ );\
    }}
#else
#define VERIFYOBJEKT( O, Oid )
#define BADOBJEKT( O )
#endif


        


#endif  /* #ifdef OBJEKT_H */
