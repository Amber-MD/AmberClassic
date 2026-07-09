/*
 *      File: unit.h
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
 *              UNIT
 *      Superclass: 
 *              CONTAINER
 *
 *      Description:
 *
 *              UNITs can contain molecules residues and/or atoms
 *              also parameters.
 */
#ifndef UNIT_H
#define UNIT_H

/*
 *-----------------------------------------------------------------------
 *
 *       Define object typedefs here.
 *        
 *        Object typedef MUST include the superclass object as
 *        its first structure element.
 *        
 *        UNITs have a mode, in UNITNORMAL mode, the unit is defined
 *        by the contents of the container and the pointers between
 *        the contents.  In UNITNORMAL mode the coordinates of the atoms
 *        are EXTERNAL coordinates.
 *        In UNITTABLES mode the UNIT is represented by the tables
 *        that have been read in, or generated from the contents.
 *        In UNITINTERNAL the UNIT is represented by the contents and 
 *        pointers between the contents, the coordinates of the atoms
 *        are represented by INTERNAL coordinates attached to the
 *        atoms.
 */


#include        "atom.h"
#include        "bag.h"
#include        "restraint.h"
#include        "dictionary.h"
#include        "parmLib.h"


#define UNITNORMAL      1
#define UNITTABLES      2
#define UNITINTERNAL    3


        /* UNIT flags */

#define UNITALLFLAGS            0xFFFFFFFF
#define UNITUSEBOUNDINGBOX      0x00000001
#define UNITBEINGEDITED         0x00000002
#define UNITUSESOLVENTCAP       0x00000004
#define UNITBOXOCT              0x00000008




typedef struct  {
        CONTAINERt      cHeader;
        OBJEKT          aHead;
        OBJEKT          aTail;
        PARMSET         psParameters;
        FLAGS           fFlags;
        int             iMode;
        double          dAlpha, dBeta, dGamma;
        double          dXWidth;
        double          dYWidth;
        double          dZWidth;
        VECTOR          vCapOrigin;
        double          dCapRadius;

        STRING          sDescription;
        DICTIONARY      dHeterogens;

                        /* Dictionary of named groups (LISTs) of ATOMs */
                        
        DICTIONARY      dAtomGroups;

                        /* Stuff to maintain RESTRAINTs */

        BAGLOOP         blRestraintLoop;
        BAG             bRestraints;
                
                        /* The following are PRIVATE, they are used */
                        /* to store interactions and pointers to parameters */
                        /* for writing UNITs to DATABASEs and for writing */
                        /* UNITs with parameters for SPASMS */
        int             iMoleculeCount;
        int             iCapTempInt;
        VARARRAY        vaAtoms;       // SAVEATOMt
        VARARRAY        vaBonds;       // SAVEBONDt
        VARARRAY        vaAngles;      // SAVEANGLEt
        VARARRAY        vaTorsions;    // SAVETORSIONt
        VARARRAY        vaConnectivity;// SAVECONNECTIVITYt
        VARARRAY        vaRestraints;  // SAVERESTRAINTt
        VARARRAY        vaResidues;    // SAVERESIDUEt
        VARARRAY        vaMolecules;   // SAVEMOLECULEt
        VARARRAY        vaHierarchy;   // SAVEHIERARCHYt (unused?)
        VARARRAY        vaConnect;     // int
        VARARRAY        vaGroupNames;  // STRING
        VARARRAY        vaGroupAtoms;  // SAVEGROUPSt

// from  zUnitIOFindAndCountMolecules()
        int             iFirstSolvent;
        VARARRAY        vaAtomsPerMolecule; // int (was also called vaMolecules and not part of UNIT struct)
} UNITt;

typedef UNITt   *UNIT;


#ifdef DEBUG
static inline CONTAINER container_from_unit(UNIT u) { return u ? &(u->cHeader) : NULL; }
static inline OBJEKT objekt_from_unit(UNIT u) { return u ? &(u->cHeader.oHeader) : NULL; }
static inline UNIT unit_from_objekt(OBJEKT o)
  { return o ? (assert(iObjectType(o)==UNITid), (UNIT)o) : NULL; }
static inline UNIT unit_from_container(CONTAINER c)
  { return c ? (assert(iObjectType(&(c->oHeader))==UNITid), (UNIT)c) : NULL; }
static inline UNIT unit_from_genp(void *p)
  { return p ? (assert(iObjectType((OBJEKT)p)==UNITid),(UNIT)p) : NULL; }
static inline UNIT unit_from_unit(UNIT u) { return u; }
#define UNIT_from(x) _Generic((x), \
    OBJEKT: unit_from_objekt, \
    CONTAINER: unit_from_container, \
    UNIT: unit_from_unit, \
    GENP: unit_from_genp \
)(x)
#else
#define UNIT_from(x) ((UNIT)(x))
#endif

/*
======================================================================

        Define object messages here.
        
        There must be at least a Create, Destroy, and Describe message.
        Hook into the messages of the superclasses so that
        when the message is sent to the most primitive superclass
        of this class that it will eventually make it into these routines.
*/


/*      Define Create, Destroy, Describe methods */

extern UNIT     uUnitCreate(void);
extern void     UnitDelete(UNIT *uPUnit);
extern void     UnitDescribe(UNIT uUnit);
extern UNIT     uUnitDuplicate(UNIT uOld);
extern UNIT     uUnitCopy(UNIT uOld);
extern void     UnitResetPointers(UNIT uUnit);

extern void     UnitJoin(UNIT uA, UNIT uB);
extern void     UnitSequence(UNIT uA, UNIT uB);


extern void     UnitSave(UNIT uUnit, DATABASE db, PARMLIB plParameters);
extern UNIT     uUnitLoad(DATABASE db);
extern void     UnitCheck(UNIT uUnit, int *iPErrors, int *iPWarnings);
extern void     UnitCheckForParms( UNIT uUnit, PARMLIB plParms, 
                        PARMSET psParmSet );

extern void     UnitSaveAmberParmFile(UNIT uUnit, char *prmtopName, char *crdName,
                        PARMLIB plParms, bool bPolar, bool bPert, bool bNetcdf);

extern void     UnitYouAreBeingRemoved(UNIT uUnit);
extern void     UnitIAmBeingRemoved(UNIT uUnit, CONTAINER cRemoved);

extern bool     bUnitCanBePerturbed(UNIT uUnit);

                /* Restraint add/remove/loop */

extern void     UnitAddRestraint(UNIT uUnit, RESTRAINT rRest );
extern bool     bUnitRemoveRestraint(UNIT uUnit, RESTRAINT rRest);
extern void     UnitLoopRestraints(UNIT uUnit);
extern RESTRAINT        rUnitNextRestraint(UNIT uUnit);

extern int      iUnitRestraintTypeCount(UNIT uUnit, int iType);

extern void     UnitSetAttribute(UNIT uUnit, STRING sAttr, OBJEKT oAttr);
extern OBJEKT   oUnitGetAttribute(UNIT uUnit, STRING sAttr);

extern bool     bUnitCapContainsAtom(UNIT uUnit, ATOM aAtom);
extern bool     bUnitCapContainsContainer(UNIT uUnit, CONTAINER cCont);

extern bool     bUnitGroupCreate(UNIT uUnit, char *cPName);
extern LIST     lUnitGroup(UNIT uUnit, char *sGroup);
extern bool     bUnitGroupAddAtom(UNIT uUnit, char *sGroup, ATOM aAtom);
extern bool     bUnitGroupFindAtom(UNIT uUnit, char *sGroup, ATOM aAtom, 
                        bool *bPFound);
extern bool     bUnitGroupRemoveAtom(UNIT uUnit, char *sGroup, ATOM aAtom);
extern bool     bUnitGroupDestroy(UNIT uUnit, char *sGroup);
extern void     UnitFindBoundingBox(UNIT uUnit, VECTOR *vPLower, 
                        VECTOR *vPUpper);

extern void     UnitSetUseBox(UNIT uUnit, bool b);
extern void     UnitSetBoxOct(UNIT uUnit, bool b);
extern void     UnitSetUseSolventCap(UNIT uUnit, bool b);
extern void     UnitDestroy( UNIT *uPUnit );
extern bool     zbUnitIgnoreHwHwOwAngle( STRING sA, STRING sB, STRING sC );
extern bool     zbUnitIgnoreAngle( STRING sA, STRING sB, STRING sC );

#define uCopyUnit(u) UNIT_from(oCopy(OBJEKT_from(u)))

#define sUnitDescription( u )     ( UNIT_from(u)->sDescription )
#define UnitSetDescription( u, s ) ( StringCopyMax( UNIT_from(u)->sDescription, s,\
                                        sizeof(STRING) ) )
#define UnitUseParameters( u, p ) \
        ( UNIT_from(u)->psParameters = (p),CDU(u) )
#define psUnitParameters( u )     ( UNIT_from(u)->psParameters )

#define bUnitHeadUsed(u)                ( UNIT_from(u)->aHead != NULL )
#define aUnitHead(u)                    ( ATOM_from(UNIT_from(u)->aHead) )
#define UnitSetHead( u, a )             ( UNIT_from(u)->aHead = (OBJEKT)(a),CDU(u) )

#define bUnitTailUsed(u)        ( UNIT_from(u)->aTail != NULL )
#define aUnitTail(u)            ( ATOM_from(UNIT_from(u)->aTail) )
#define UnitSetTail( u, a )     ( UNIT_from(u)->aTail = (OBJEKT)(a),CDU(u) )

#define UnitSetMode( u, m )             (UNIT_from(u)->iMode = m,CDU(u) )
#define iUnitMode( u )                  (UNIT_from(u)->iMode)
#define dUnitAtomGroups(u)              (UNIT_from(u)->dAtomGroups)
#define UnitGetBox( u, xP, yP, zP ) (\
        *(xP) = UNIT_from(u)->dXWidth,\
        *(yP) = UNIT_from(u)->dYWidth,\
        *(zP) = UNIT_from(u)->dZWidth )
#define UnitSetBox( u, x, y, z ) (\
        UNIT_from(u)->dXWidth = (x),\
        UNIT_from(u)->dYWidth = (y),\
        UNIT_from(u)->dZWidth = (z),CDU(u) )
#define UnitGetCell( u, xP, yP, zP, aP, bP, gP ) (\
        *(xP) = UNIT_from(u)->dXWidth,\
        *(yP) = UNIT_from(u)->dYWidth,\
        *(zP) = UNIT_from(u)->dZWidth,\
        *(aP) = UNIT_from(u)->dAlpha,\
        *(bP) = UNIT_from(u)->dBeta,\
        *(gP) = UNIT_from(u)->dGamma)
#define UnitSetCell( u, x, y, z, a, b, g ) (\
        UNIT_from(u)->dXWidth = (x),\
        UNIT_from(u)->dYWidth = (y),\
        UNIT_from(u)->dZWidth = (z),\
        UNIT_from(u)->dAlpha = (a),\
        UNIT_from(u)->dBeta = (b),\
        UNIT_from(u)->dGamma = (g), CDU(u) )
#define UnitGetSolventCap( u, xP, yP, zP, rP ) {\
        (*xP) = dVX(&(u->vCapOrigin));\
        (*yP) = dVY(&(u->vCapOrigin));\
        (*zP) = dVZ(&(u->vCapOrigin));\
        (*rP) = u->dCapRadius;}
#define UnitSetSolventCap( u, x, y, z, r ) {\
        VectorDef( &(u->vCapOrigin), x, y, z );\
        u->dCapRadius = r;CDU(u); }
#define UnitSetBeta( u, d )     ( UNIT_from(u)->dBeta = d, \
        UNIT_from(u)->dAlpha = UNIT_from(u)->dGamma = 90.0*DEGTORAD, CDU(u))
#define dUnitBeta( u )          ( UNIT_from(u)->dBeta )
#define UnitDefineFlags( u, f ) ( UNIT_from(u)->fFlags = (f), CDU(u) )
#define UnitSetFlags( u, f )    { UNIT_from(u)->fFlags |= (f); CDU(u); }
#define UnitResetFlags( u, f )  { UNIT_from(u)->fFlags &= ~(f); CDU(u); }
#define bUnitFlagsSet( u, f )   ( (UNIT_from(u)->fFlags & (f)) == (f) )
#define bUnitUseBox(u)          (bUnitFlagsSet(u,UNITUSEBOUNDINGBOX))
#define bUnitBoxOct(u)          (bUnitFlagsSet(u,UNITBOXOCT))
#define bUnitUseSolventCap(u)           bUnitFlagsSet(u,UNITUSESOLVENTCAP)


#endif /* UNIT_H */
