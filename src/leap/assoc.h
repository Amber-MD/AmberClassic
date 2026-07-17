/*
 *      File:   assoc.h
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
 *              ASSOC
 *      Superclass: 
 *              OBJEKT
 *
 *      Description:
 *
 *              An ASSOC (short for association) is an object that
 *              stores a short name and an OBJEKT.
 *              It associates the name with the OBJEKT.
 *              This is used by the command line interface to
 *              store objects with the names that are bound to them.
 *
 */
 

#ifndef ASSOC_H
#define ASSOC_H


/*
-----------------------------------------------------------------------

        Define object typedefs here.
        
        Object typedef MUST include the superclass object as
        its first structure element.
*/

typedef struct  ASSOCSTRUCT {
        OBJEKTt         oSuper;
        STRING          sName;
        OBJEKT          oObj;
} ASSOCt;

typedef ASSOCt  *ASSOC;

#ifdef DEBUG
static inline ASSOC assoc_from_assoc(ASSOC a,const char *file,int line) { return a; }
static inline ASSOC assoc_from_objekt(OBJEKT o,const char *file,int line) {
    if (!o) return NULL;
    VERIFYOBJEKT_loc(o,ASSOCid,file,line);
    return (ASSOC)o;
}
static inline OBJEKT objekt_from_assoc(ASSOC a) { return &(a->oSuper); }

#define ASSOC_from(x) _Generic((x), \
    ASSOC: assoc_from_assoc, \
    OBJEKT: assoc_from_objekt \
)(x,__FILE__,__LINE__)
#else
#define ASSOC_from(x) ((ASSOC)(x))
#endif

/*
======================================================================

        Define object messages here.
        
        There must be at least a Create, Destroy, and Describe message.
        Hook into the messages of the superclasses so that
        when the message is sent to the most primitive superclass
        of this class that it will eventually make it into these routines.
*/


#define AssocSetName(a,n)       (strcpy(ASSOC_from(a)->sName,n))
#define sAssocName(a)           (ASSOC_from(a)->sName)
#define AssocSetObject(a,o)     { OBJEKT _o = OBJEKT_from(o); \
                                ASSOC_from(a)->oObj = _o; \
                                REF(_o); }
#define AssocTakeObject(a,o)    (ASSOC_from(a)->oObj = OBJEKT_from(o))
#define oAssocObject(a)         ( ASSOC_from(a)->oObj )


/*  assoc.c  */

extern ASSOC            aAssocCreate(void);
extern void             AssocDestroy( ASSOC *aPAssoc );
extern void             AssocDescribe( ASSOC aAssoc );

// helper function:
extern OBJEKT           CreateAssocVector(VECTOR vVector);


#endif /* ASSOC_H */
