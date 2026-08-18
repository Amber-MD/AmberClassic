#ifndef ATOM_H
#define ATOM_H

/*
 *      File: atom.h
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
 *             David A. Rivkin                                          *
 *                                                                      *
 *     Principal Investigator: Peter A. Kollman                         *
 *                                                                      *
 ************************************************************************
 *
 *      Class: 
 *              ATOM
 *      Superclass: 
 *              CONTAINER
 *
 *      Description:
 *
 *              Atoms can contain anything.
 */
 
#include        "basics.h"
#include        "container.h"
#include        "vector.h"

#define ATOM_DEFAULT_RADIUS     1.5

                        /* 5 characters to represent the type */
#define ATOMTYPELEN     5
                        /* 10 characters to represent the name */
#define ATOMNAMELEN     10
                        /* maximum 8 bonds out of each atom */
#define MAXBONDS        8


/*
-----------------------------------------------------------------------

        Define object typedefs here.
        
        Object typedef MUST include the superclass object as
        its first structure element.
*/

#include        "elements.h"
//#include        "unitio.h" //New

                /* Hybridization type */

#define HUNDEFINED      -1
#define HUNKNOWN        0
#define HSP1            1
#define HSP2            2
#define HSP3            3


                /* Bond flags */
#define BONDORDERONLY   0x0000000F
#define BONDNONE        0x00000000
#define BONDSINGLE      0x00000001
#define BONDDOUBLE      0x00000002
#define BONDTRIPLE      0x00000003
#define BONDAROMATIC    0x00000004

#define BONDTEMPORARY   0x00000010


typedef char            BONDt;
typedef char            ATOMTYPEt[ATOMTYPELEN];
typedef char            atomNameStr[ATOMNAMELEN];

typedef struct  ATOMSTRUCT {
        CONTAINERt              cHeader;
        int                     iUniqueId;
        int                     iAtomicNumber;
        int                     iPertAtomicNumber;
        CONTAINERNAMEt          sPertName;
        ATOMTYPEt               sType;
        char                    iCode, altLoc;
        ATOMTYPEt               sPertType;
        double                  dCharge;
        double                  dPertCharge;
        double                  dPolar;
        double                  dPertPolar;
        double                  dScreenF;
        int                     iIndex;
        VECTOR                  vPosition;
        VECTOR                  vVelocity;
        FLAGS                   fFlags;         /* Atom flags */
        int                     iCoordination;
        struct ATOMSTRUCT       *aaBonds[MAXBONDS];
        FLAGS                   faBondFlags[MAXBONDS];
        double                  dTemp;
        GENP                    PTemp;

        /* Spanning tree stuff */
        /* Used for LOOPing over SPANNINGTREE */
        struct ATOMSTRUCT       *aNextSpan;
        struct ATOMSTRUCT       *aBackSpan;
        int                     iSeenId;
        int                     iBackCount;
                /* Graphics stuff */
        GENP                    PGraphicsData;
} ATOMt;

typedef ATOMt   *ATOM;


                /* Atom flags */

#define ATOMALLFLAGS            0xFFFFFFFF

#define ATOMPERMANENTFLAGS      0x0000FFFF

                                /* Differentiate between the atoms position */
                                /* being FIXED or being BUILT               */
                                /* Depending on whether the atoms position  */
                                /* was read in from PDB or built using the  */
                                /* builder, Also define whether the position*/
                                /* is known or not.  This makes visibility */
                                /* testing in spanning tree loops easy      */
#define ATOMSELECTED            0x00000001
                                /* This flag must be set if the ATOM is to */
                                /* be perturbed */

#define ATOMPERTURB             0x00000002
#define ATOMNOTDISPLAYED        0x00000004
#define RESIDUEIMAGEATOM        0x00000010      /* not used (so far) */
#define ATOMBULKSOLVENT         0x00000020

#define ATOMTEMPORARYFLAGS      0xFFFF0000
        
#define ATOMTOUCHED             0x00010000
#define ATOMPOSITIONKNOWN       0x00020000      /* ATOMPOSITIONEXTERNAL */
#define ATOMPOSITIONINTERNAL    0x00040000
#define ATOMNEEDSMINIMIZER      0x00080000
#define ATOMNEEDSBUILD          0x00100000
#define ATOMTOUCHED2            0x00200000
#define ATOMPOSITIONFIXED       0x00400000
#define ATOMPOSITIONBUILT       0x00800000
#define ATOMPOSITIONDRAWN       0x01000000
#define ATOMMOCKPOSITIONKNOWN   0x02000000
#define ATOMFLAGUNUSED0         0x04000000
#define ATOMFLAGUNUSED1         0x08000000
#define ATOMFLAG0               0x10000000
#define ATOMFLAG1               0x20000000
#define ATOMFLAG2               0x40000000
#define ATOMFLAG3               0x80000000


#ifdef DEBUG
static inline CONTAINER container_from_atom(ATOM a,const char *f, int l) { return a ? &(a->cHeader) : NULL; }
static inline OBJEKT objekt_from_atom(ATOM a) { return a ? &(a->cHeader.oHeader) : NULL; }
static inline ATOM atom_from_objekt(OBJEKT o, const char *file, int line)
  { return o ? (assert_loc(iObjectType(o)==ATOMid,file,line),(ATOM)o) : NULL; }
static inline ATOM atom_from_container(CONTAINER c, const char *file, int line)
  { return c ? (assert_loc(iObjectType(&(c->oHeader))==ATOMid,file,line), (ATOM)c) : NULL; }
static inline ATOM atom_from_atom(ATOM a,const char *file,int line) { return a; }

#define ATOM_from(x) _Generic((x), \
    OBJEKT: atom_from_objekt, \
    CONTAINER: atom_from_container, \
    ATOM: atom_from_atom \
)(x,__FILE__,__LINE__)
#else
#define ATOM_from(x) ((ATOM)(x))
#endif

/*
======================================================================

        Define object messages here.
        
        There must be at least a Create, Destroy, and Describe message.
        Hook into the messages of the superclasses so that
        when the message is sent to the most primitive superclass
        of this class that it will eventually make it into these routines.
*/
#define AtomDefineBondFlags(a,i,f)      ((ATOM_from(a))->faBondFlags[i] = f,CDU(a) )
#define AtomSetBondFlags(a,i,f)         ((ATOM_from(a))->faBondFlags[i] |= f,CDU(a) )
#define AtomResetBondFlags(a,i,f)       \
                ((ATOM_from(a))->faBondFlags[i] &= ~(f),CDU(a) )
#define AtomSetBondOrder(a,i,f) ( AtomResetBondFlags(a,i,BONDORDERONLY),\
                                  AtomSetBondFlags(a,i,f),CDU(a) )
#define iAtomBondOrder(a,i)     (((ATOM_from(a))->faBondFlags[i])&(BONDORDERONLY))
#define fAtomBondFlags(a,i)     ((ATOM_from(a))->faBondFlags[i])
#define bAtomBondFlagsSet(a,i,f)        (((ATOM_from(a))->faBondFlags & f)!=0)
#define aAtomBondedNeighbor(a,i)        ((ATOM_from(a))->aaBonds[i])
#define AtomSetTempPtr( a, p )  ((ATOM_from(a))->PTemp = (GENP)(p))
#define PAtomTempPtr(a)                 ((ATOM_from(a))->PTemp)
#define AtomSetElement( a, n )          ((ATOM_from(a))->iAtomicNumber = n,CDU(a))
#define iAtomElement(a)                 ((ATOM_from(a))->iAtomicNumber)
#define AtomSetPertElement( a, n )      ((ATOM_from(a))->iPertAtomicNumber = n,CDU(a))
#define iAtomPertElement(a)             ((ATOM_from(a))->iPertAtomicNumber)
#define AtomDupPosition(a,vv)   {(ATOM_from(a))->fFlags|=ATOMPOSITIONKNOWN;\
                                        VectorCopy(vv, vAtomPosition(ATOM_from(a))); CDU(a);}
#define AtomSetPosition(a,vv)   ((ATOM_from(a))->fFlags|=ATOMPOSITIONKNOWN,\
                                        vAtomPosition(ATOM_from(a)) = vv, CDU(a))
#define AtomSetPositionNoFlags(a,vv)    {VectorCopy(vv, vAtomPosition(ATOM_from(a))); CDU(a);}
#define vAtomPosition(a)                ((ATOM_from(a))->vPosition)
#define AtomDupVelocity(a,vv)           {VectorCopy(vv, vAtomVelocity(ATOM_from(a))); CDU(a);}
#define AtomSetVelocity(a,vv)           (vAtomVelocity(ATOM_from(a)) = vv, CDU(a))
#define vAtomVelocity(a)                ((ATOM_from(a))->vVelocity)
#define iAtomCoordination(a)            ((ATOM_from(a))->iCoordination)
#define iAtomId(a)                      ((ATOM_from(a))->iUniqueId)
#define AtomSetType( a, x )             (strcpy((ATOM_from(a))->sType,x),CDU(a))
#define sAtomType( a )                  ((ATOM_from(a))->sType )
#define AtomSetName( a, s )             ContainerSetName( a, s )
#define sAtomName( a )                  sContainerName( a )
#define AtomSetPertName( a, s )         (strcpy((ATOM_from(a))->sPertName, s ),CDU(a))
#define sAtomPertName( a )              ((ATOM_from(a))->sPertName)
#define AtomSetPertType( a, x )         (strcpy((ATOM_from(a))->sPertType,x),CDU(a))
#define sAtomPertType( a )              ((ATOM_from(a))->sPertType )
#define bAtomPerturbed( a )             (bAtomFlagsSet((ATOM_from(a)),ATOMPERTURB))
#define AtomDefineFlags(a,f)            ((ATOM_from(a))->fFlags = f,CDU(a) )
#define AtomSetFlags(a,f)               ((ATOM_from(a))->fFlags |= f,CDU(a) )
#define AtomResetFlags(a,f)             ((ATOM_from(a))->fFlags &= (~f),CDU(a) )
#define fAtomFlags(a)                   ((ATOM_from(a))->fFlags)
#define bAtomFlagsSet(a,f)              (((ATOM_from(a))->fFlags&(f))==(f))
#define bAtomFlagsReset(a,f)     ((~((ATOM_from(a))->fFlags)&(f))==(f))
#define aAtomNextSpan(AA)               (ATOM)(((ATOM)(AA))->aNextSpan)
#define AtomSetNextSpan(AA,SS)          (((ATOM)(AA))->aNextSpan = SS )
#define aAtomBackSpan(AA)               (ATOM)(((ATOM)(AA))->aBackSpan) 
#define AtomSetBackSpan(AA,SS)          (((ATOM)(AA))->aBackSpan = SS )
#define iAtomBackCount(AA)              (((ATOM)(AA))->iBackCount)
#define AtomSetBackCount(AA,I)          (((ATOM)(AA))->iBackCount=(I))
#define iAtomSeenId(AA)                 (((ATOM)(AA))->iSeenId)
#define AtomSetSeenId(AA,II)            (((ATOM)(AA))->iSeenId =(II))
#define AtomSetTempInt(AA,II)           (ContainerSetTempInt(AA,II))
#define iAtomTempInt(AA)                (iContainerTempInt(AA))
#define AtomSetCharge( a, x )           ((ATOM_from(a))->dCharge = x,CDU(a))
#define AtomSetPolar( a, x )            ((ATOM_from(a))->dPolar = x,CDU(a))
#define dAtomCharge( a )                ((ATOM_from(a))->dCharge)
#define dAtomPolar( a )                 ((ATOM_from(a))->dPolar)
#define AtomSetPertCharge( a, x )       ((ATOM_from(a))->dPertCharge = x,CDU(a))
#define AtomSetPertPolar( a, x )        ((ATOM_from(a))->dPertPolar = x,CDU(a))
#define dAtomPertCharge( a )            ((ATOM_from(a))->dPertCharge)
#define dAtomPertPolar( a )             ((ATOM_from(a))->dPertPolar)
#define iAtomIndex( a )                 ((ATOM_from(a))->iIndex)
#define AtomSetIndex( a, i )            ((ATOM_from(a))->iIndex = i)
#define AtomSetGraphicsPointer(a,p) ((ATOM_from(a))->PGraphicsData=((GENP)p))
#define PAtomGraphicsPointer(a)  ((ATOM_from(a))->PGraphicsData)

#define AtomSetTempDouble(a,v)          ((ATOM_from(a))->dTemp = v)
#define AtomTempDoubleIncrement(a,v)    ((ATOM_from(a))->dTemp += v)
#define AtomTempDoubleSquare(a)         ((ATOM_from(a))->dTemp *= \
                                                (ATOM_from(a))->dTemp)
#define AtomTempDoubleSquareRoot(a)     ((ATOM_from(a))->dTemp = \
                                                sqrt((ATOM_from(a))->dTemp) )
#define dAtomTemp(a)                    ((ATOM_from(a))->dTemp)

#define bAtomVisible(a) ((bAtomFlagsSet(a,ATOMPOSITIONKNOWN)||\
                          bAtomFlagsSet(a,ATOMPOSITIONDRAWN))&&\
                          !bAtomFlagsSet(a,ATOMNOTDISPLAYED))
#define bAtomHasPosition(a)     (bAtomFlagsSet(a,ATOMPOSITIONKNOWN)||\
                                bAtomFlagsSet(a,ATOMPOSITIONDRAWN))

#define cAtomAltLoc(a)           ((a)->cAltLoc)
#define AtomSetAltLoc(a,c)       ((a)->cAltLoc=(c==' ')?c:0)
#define cAtomInsertCode(a)       ((a)->cInsertCode)
#define AtomSetInsertCode(a,c)   ((a)->cInsertCode=(c==' ')?c:0)


/*  atom.c  */

extern ATOM             aAtomCreate(void);
extern void             AtomDestroy( ATOM *aPAtom );
extern void             AtomDescribe( ATOM aAtom );
extern void             AtomDescStr( ATOM aA, bool bResNum, char *cPDesc );
extern bool             bAtomCoordinationSaturated( ATOM aAtom );
extern bool             AtomTmpBondTo( ATOM aAtom1, ATOM aAtom2 );
extern void             AtomBondToOrder( ATOM aAtom1, ATOM aAtom2, int iOrder );
extern void             AtomBondToFlags( ATOM aAtom1, ATOM aAtom2, 
                                FLAGS fFlags );
extern void             AtomRemoveBond( ATOM aAtom1, ATOM aAtom2 );
extern bool             bAtomBondedTo( ATOM aAtom1, ATOM aAtom2 );
extern int              iAtomFindBondOrder( ATOM aAtom, ATOM aNeighbor );
extern void             AtomFindSetBondOrder( ATOM aAtom, ATOM aNeighbor, 
                                int iOrder );
extern FLAGS            fAtomFindBondFlags( ATOM aAtom, ATOM aNeighbor );
extern ATOM             aAtomDuplicate( ATOM aAtom );
extern void             AtomResetPointers( ATOM aAtom );
extern void             AtomCheck( ATOM aAtom, int *iPErrors, int *iPWarnings );
extern void             AtomYouAreBeingRemoved( ATOM aAtom );
extern void             AtomIAmBeingRemoved( ATOM aAtom, CONTAINER cRemoved );
extern void             AtomSetAttribute( ATOM aAtom, STRING sAttr, 
                                OBJEKT oAttr );
extern OBJEKT           oAtomGetAttribute( ATOM aAtom, STRING sAttr );
extern int              iAtomHybridization( ATOM aAtom );
extern int              iAtomBondOrderFromName( char *sName );
extern bool             bAtomSpaceConflict( ATOM aAtom1, ATOM aAtom2 );
extern double           dAtomVanderWaals( ATOM aAtom );
extern bool             bAtomSetTmpRadius( ATOM aAtom );

extern void             AtomBondTo( ATOM aAtom1, ATOM aAtom2 );

#endif /* ATOM_H */
