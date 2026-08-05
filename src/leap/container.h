/*
 *      File: container.h
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
 *              CONTAINER
 *      Superclass: 
 *              OBJEKT, LOOP
 *      SubClasses:
 *              MOLSYSTEM, MOLECULE, RESIDUE, ATOM
 *
 *      Description:
 *
 *              CONTAINER is a superclass for all objects used to
 *              define a system of molecules.  The properties that
 *              it defines are:
 *              Every CONTAINER can be contained by another.
 *              Every CONTAINER can have contents.
 *              Every CONTAINER has a name.
 *              Looping over CONTAINER contents is handled by
 *              the LOOP class.
 */
 
#ifndef CONTAINER_H
#define CONTAINER_H


# include       "objekt.h"
# include       "displayer.h"
# include       "matrix.h"
# include       "stringExtra.h"

/*
 *-----------------------------------------------------------------------
 *
 *        Define object typedefs here.
 *        
 *        Object typedef MUST include the superclass object as
 *        its first structure element.
 */



#define CONTAINERNAMELEN        32
typedef char    CONTAINERNAMEt[CONTAINERNAMELEN];

typedef struct  CONTAINERSTRUCT {
        OBJEKTt                 oHeader;
        struct CONTAINERSTRUCT  *cCopy;
        CONTAINERNAMEt          sName;
        struct CONTAINERSTRUCT  *cContainedBy;
        int                     iNextChildsSequence;
        DISPLAYER               dDisp;
        int                     iTempInt;
        int                     iSequence;
        struct CONTAINERSTRUCT  *cLoopNext;     /* Used for LOOPing over DIRECTCONTENTSBYSEQNUM */
        OBJEKT                  lContents;   // Should be  struct LISTSTRUCT *lContents;
        struct NodeStruct       *nPListNode;  // tracks member node in lContents
} CONTAINERt;

typedef CONTAINERt      *CONTAINER;

#ifdef DEBUG
static inline CONTAINER container_from_objekt(OBJEKT o)
  { assert ( bObjectInClass( o, CONTAINERid ) ); return (CONTAINER)o; }
static inline CONTAINER container_from_genp(void *p)
  { assert(bIsObjekt(p)); assert(bObjectInClass((OBJEKT)p,CONTAINERid)); return (CONTAINER)p; }
static inline CONTAINER container_from_container(CONTAINER c) { return c; }
static inline OBJEKT objekt_from_container(CONTAINER c) { return c ? &(c->oHeader) : NULL; }

#define CONTAINER_from(x) _Generic((x), \
    OBJEKT: container_from_objekt, \
    CONTAINER: container_from_container, \
    MOLECULE: container_from_molecule, \
    UNIT: container_from_unit, \
    RESIDUE: container_from_residue, \
    ATOM: container_from_atom, \
    GENP: container_from_genp \
)(x)
#else
#define CONTAINER_from(x) ((CONTAINER)(x))
#endif

/*
======================================================================

        Define object messages here.
        
        There must be at least a Create, Destroy, and Describe message.
        Hook into the messages of the superclasses so that
        when the message is sent to the most primitive superclass
        of this class that it will eventually make it into these routines.
*/

#ifdef NOGUI
#define CDU(c) ((void)0)
#else
#define CDU(c)  ContainerDisplayerUpdate(CONTAINER_from(c))
#endif

#define sContainerName( c )     ( CONTAINER_from(c)->sName )
#define ContainerSetName( c, n ) ( StringCopyMax( CONTAINER_from(c)->sName, n,\
                                        sizeof(CONTAINERNAMEt) ),CDU(c) )
#define cContainerWithin( c )   ( CONTAINER_from(c)->cContainedBy )
#define ContainerSetWithin( c, w ) (CONTAINER_from(c)->cContainedBy = w,CDU(c) )
#define cContainerContents( c )    (CONTAINER_from(c)->lContents)

#define iContainerNumberOfChildren(c) (iListSize(cContainerContents(c)))
// Convenience speed-up to grab the first object in a single-element container // JMK
#define oContainerFirstObject( c ) (LIST_from(cContainerContents(c))->nPFirstNode->PObject)
#define oContainerLastObject( c ) (LIST_from(cContainerContents(c))->nPLastNode->PObject)


                /* If the copy is NULL then return the original */
                /* otherwise return the copy                    */
#define cContainerCopyPointer( c ) \
( c == NULL ? NULL : ( ((CONTAINER_from(c)->cCopy)!=NULL) ? \
(CONTAINER)(CONTAINER_from(c)->cCopy) : CONTAINER_from(c) ))

#define ContainerSetCopyPointer( c, n )   (CONTAINER_from(c)->cCopy=n)

#define ContainerSetTempInt(c,i)        (CONTAINER_from(c)->iTempInt = i )
#define iContainerTempInt(c)            (CONTAINER_from(c)->iTempInt)
#define iContainerSequence( c )         (CONTAINER_from(c)->iSequence)
#define ContainerSetSequence( c,xn )    (CONTAINER_from(c)->iSequence = xn,CDU(c))
#define iContainerNextChildsSequenceInc( c ) \
        (CONTAINER_from(c)->iNextChildsSequence++)
#define iContainerNextChildsSequence( c ) \
        (CONTAINER_from(c)->iNextChildsSequence)
#define ContainerSetNextChildsSequence( c, i ) \
        (CONTAINER_from(c)->iNextChildsSequence = i,CDU(c))


#define ContainerSetLoopNext( c, n )    ( CONTAINER_from(c)->cLoopNext = n )
#define cContainerLoopNext(c)           ( CONTAINER_from(c)->cLoopNext )

#define dContainerDisplayer(c)  (CONTAINER_from(c)->dDisp)

#define sFullDescriptor(x,s) sContainerFullDescriptor(CONTAINER_from(x),s)

extern CONTAINER        cContainerCreate( int iType );
extern void             ContainerDestroy( CONTAINER *cPContainer );
extern void             ContainerDescribe( CONTAINER cContainer  ) ;
extern void             ContainerAdd( CONTAINER cContainer, OBJEKT oObject );
extern bool             bContainerRemove( CONTAINER cContainer, 
                                OBJEKT oObject );
extern CONTAINER        cContainerDuplicate( CONTAINER cOld );
extern void             ContainerResetPointers( CONTAINER cContainer );
extern CONTAINER        cContainerFindSequence( CONTAINER cCont, 
                                int iContainerType, int iSeq );
extern CONTAINER        cContainerFindName( CONTAINER cCont, 
                                int iContainerType, const char *sName );
extern VECTOR           vContainerGeometricCenter( CONTAINER cCont );
extern void             ContainerBoundingBox( CONTAINER cCont, 
                                VECTOR *vPMin, VECTOR *vPMax );
extern void             ContainerCenterAt( CONTAINER cCont, VECTOR vCenter );
extern void             ContainerTransformBy( CONTAINER cCont, 
                                MATRIX mTransform );
extern void             ContainerTranslateBy( CONTAINER cCont, VECTOR vOffset );
extern void             ContainerSetAllAtomsFlags( CONTAINER cCont, 
                                FLAGS fFlags );
extern void             ContainerResetAllAtomsFlags( CONTAINER cCont, 
                                FLAGS fFlags );
extern void             ContainerWithFlagsSetAtomFlags( CONTAINER cCont, 
                                FLAGS fNeed, FLAGS fFlags );
extern void             ContainerWithFlagsResetAtomFlags( CONTAINER cCont, 
                                FLAGS fNeed, FLAGS fFlags );
extern void             ContainerWithoutFlagsSetAtomFlags( CONTAINER cCont, 
                                FLAGS fNeed, FLAGS fFlags );
extern void             ContainerWithoutFlagsResetAtomFlags( CONTAINER cCont,
                                FLAGS fNeed, FLAGS fFlags );
extern char             *sContainerDescriptor( CONTAINER cCont, char *sDesc );
extern char             *sContainerFullDescriptor( CONTAINER cCont, 
                                char *sFullDesc );
extern void             ContainerCheck( CONTAINER cCont, int *iPErrors, 
                                int *iPWarnings );
extern bool             bContainerContainedBy( CONTAINER cIn, CONTAINER cOut );
extern void             ContainerYouAreBeingRemoved( CONTAINER cCont ) ;
extern void             ContainerIAmBeingRemoved( CONTAINER cCont, 
                                CONTAINER cRemoved );
extern void             ContainerSetAttribute( CONTAINER cCont, 
                                STRINGref sAttribute, OBJEKT oValue );
extern OBJEKT           oContainerGetAttribute( CONTAINER cCont,
                                STRINGref sAttribute );
extern void             ContainerTotalCharge( CONTAINER cCont, 
                                double *dPCharge, double *dPPertCharge );
extern void             ContainerDisplayerUpdate( CONTAINER cCont );
extern void             ContainerTreeMakeInsensitive( CONTAINER cCont );
extern void             ContainerTreeMakeSensitive( CONTAINER cCont );
extern void             ContainerResetAllCopyPointers( CONTAINER cTop );
extern bool             bContainerSpaceConflict( CONTAINER cCont1, 
                                CONTAINER cCont2 );


#endif /* CONTAINER_H */
