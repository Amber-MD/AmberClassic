/*
 *      File: residue.h
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
 *              RESIDUE
 *      Superclass: 
 *              CONTAINER
 *
 *      Description:
 *
 *              Residues can contain atoms
 */
 
#ifndef RESIDUE_H
#define RESIDUE_H

/*
-----------------------------------------------------------------------

        Define object typedefs here.
        
        Object typedef MUST include the superclass object as
        its first structure element.
*/


                /* Some containers, UNITs and RESIDUEs have connect */
                /* pointers which point to atoms to which connections */
                /* can be made easily.  Other atoms can be connected to */
                /* but they have to be refered to directly */
#define MAXCONNECT      6

#define NOEND           -1

#define	CONNECT0	0
#define	CONNECT1	1
#define	CONNECT2	2
#define	CONNECT3	3
#define	CONNECT4	4
#define	CONNECT5	5


#define FIRSTEND        CONNECT0	/* This is the end that is bare in */
                                        /* the first RESIDUE */
#define LASTEND         CONNECT1	/* This is the end that is bare in */
                                        /* the last RESIDUE */

#define NEND            FIRSTEND        /* N terminus */
#define CEND            LASTEND         /* C terminus */
#define SEND            CONNECT2	/* Disulphide bridges */

		/* Residue types       related PDB class */

#define	RESTYPEUNDEFINED	'?'
#define	RESTYPESOLVENT		'w'    // HETAS
#define	RESTYPEPROTEIN		'p'    // ATOMP
#define	RESTYPENUCLEIC		'n'    // ATOMN
#define	RESTYPESACCHARIDE	's'    // ATOMS
#define	RESTYPEION              'i'    // HETAI
#define	RESTYPELIGAND           'l'    // HETAD=drug,HETAIN=inhibitor,HETAC=coenzyme
#define	RESTYPEOTHER            'o'    // non-standard polymer, anything else
#define RESTYPECLASSPOLYMER     "nops"
#define RESTYPECLASSNONPOLYMER  "ilw"

                /* Residue flags */

#define	RESIDUEPERM	0x0000FFFF
#define RESIDUEUNKNOWN  0x00000001
#define RESIDUEBULKSOLVENT 0x00000002

// These are not used in residue, for fGetPdbResMapped()
#define RESIDUENOEND    0x00000008
#define RESIDUEFIRSTEND 0x00000010
#define RESIDUELASTEND  0x00000020

#define	RESIDUETEMP	0xFFFF0000
#define	RESIDUEINCAP	0x00010000

typedef struct {
	char	sName1[ATOMNAMELEN];
	char	sName2[ATOMNAMELEN];
	char	sName3[ATOMNAMELEN];
	char	sName4[ATOMNAMELEN];
} IMPROPERt;

typedef struct RESIDUESTRUCT {
	CONTAINERt      cHeader;
	FLAGS           fFlags;
	char		cResType;
	OBJEKT          aaConnect[MAXCONNECT];
	OBJEKT		aSolventImagingAtom;	/* ATOM */
	int		iPdbResSeq;
	char		sChainId[3], cICode;
	VARARRAY	vaImpropers;
	double		dTemp;
	int		iTemp; // do we need this? only used for default resSeq
                               // parent already has iTempInt
	int		iMolecule;
} RESIDUEt;

typedef RESIDUEt	*RESIDUE;

#ifdef DEBUG
static inline OBJEKT objekt_from_residue(RESIDUE r) { return r ? &(r->cHeader.oHeader) : NULL; }
static inline CONTAINER container_from_residue(RESIDUE r) { return r ? &(r->cHeader) : NULL; }
static inline RESIDUE residue_from_objekt(OBJEKT o)
  { return o ? (assert (iObjectType(o)==RESIDUEid), (RESIDUE)o) : NULL; }
static inline RESIDUE residue_from_container(CONTAINER c)
  { return c ? (assert(iObjectType(&(c->oHeader))==RESIDUEid), (RESIDUE)c) : NULL; }
static inline RESIDUE residue_from_residue(RESIDUE r) { return r; }
#define RESIDUE_from(x) _Generic((x), \
    OBJEKT: residue_from_objekt, \
    CONTAINER: residue_from_container, \
    RESIDUE: residue_from_residue \
)(x)
#else
#define RESIDUE_from(x) ((RESIDUE)(x))
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

extern RESIDUE		rResidueCreate(void);
extern void		ResidueDelete(RESIDUE *rPResidue);
extern void		ResidueDescribe(RESIDUE rResidue);
extern void             ResidueDestroy(RESIDUE *rPResidue);
extern RESIDUE		rResidueDuplicate(RESIDUE rOld);
extern void		ResidueResetPointers(RESIDUE rRes);
extern void		ResidueCheck(RESIDUE rRes, 
				int *iPErrors, int *iPWarnings);
extern void		ResidueMutate(RESIDUE rNew, RESIDUE rOld);
extern char		*sResidueTypeNameFromChar(char c);

extern RESIDUE		rResidueConnected(RESIDUE rRes, int iConnect);
extern int		iResidueConnectFromName(char *sName);
extern void		ResidueSetAttribute( RESIDUE rRes, 
				STRING sAttr, OBJEKT oAttr );
extern OBJEKT		oResidueGetAttribute( RESIDUE rRes,
				STRING sAttr);

extern bool	bResidueCrossLink(RESIDUE rA, int iConnectA,
			RESIDUE rB, int iConnectB, int iOrder);
extern void	ResidueYouAreBeingRemoved(RESIDUE rRes);
extern void	ResidueIAmBeingRemoved(RESIDUE rRes, CONTAINER cRemoved);

#define rCopyResidue(r) RESIDUE_from(oCopy(OBJEKT_from(r)))

#define bResidueConnectUsed(r,c)     (RESIDUE_from(r)->aaConnect[c]!=NULL)
#define ResidueSetConnectAtom(r,c,a) (RESIDUE_from(r)->aaConnect[c]=(OBJEKT)(a),\
					CDU(r))
#define aResidueConnectAtom(r,c)        (ATOM)(RESIDUE_from(r)->aaConnect[c])
#define ResidueSetDescription(r,s) (strcpy( RESIDUE_from(r)->sDescription,s),\
					CDU(r))
#define sResidueDescription(r)          (RESIDUE_from(r)->sDescription)
#define bResidueFlagsSet(r,f)           ((RESIDUE_from(r)->fFlags & f)!= 0)
#define ResidueSetFlags(r,f)            (RESIDUE_from(r)->fFlags |= f,CDU(r) )
#define ResidueDefineFlags(r,f)  (RESIDUE_from(r)->fFlags = f,CDU(r))
#define ResidueResetFlags(r,f)    (RESIDUE_from(r)->fFlags &= ~f,CDU(r))
#define	ResidueSetType(r,c)	(RESIDUE_from(r)->cResType = (c),CDU(r))
#define	cResidueType(r)		(RESIDUE_from(r)->cResType)
#define ResidueSetImagingAtom(r,a) (RESIDUE_from(r)->aSolventImagingAtom=(OBJEKT)(a),CDU(r))
#define aResidueImagingAtom(r)     (ATOM)(RESIDUE_from(r)->aSolventImagingAtom)
#define	iResiduePdbSequence(r)	(RESIDUE_from(r)->iPdbResSeq)
#define	ResidueSetPdbSequence(r,i) (RESIDUE_from(r)->iPdbResSeq=(i))
#define	sResidueChainId(r) (RESIDUE_from(r)->sChainId)
#define	ResidueSetChainId(r,s) do { \
        if ((s)[0]==' ' && (s)[1]==0) strcpy(RESIDUE_from(r)->sChainId,""); \
        else strcpy(RESIDUE_from(r)->sChainId,s); \
     } while (0);

#endif /* RESIDUE_H */
