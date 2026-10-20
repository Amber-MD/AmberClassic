/*
 *	File:	octree.h
 *
 ************************************************************************
 *                                                                      *
 *     Designed by:    Bill Ross					*
 *     Author:         Bill Ross					*
 *                                                                      *
 *     VERSION: 1.0                                                     *
 *     Programmers:                                                     *
 *             Bill Ross						*
 *                                                                      *
 *     Principal Investigator: Peter A. Kollman                         *
 *                                                                      *
 ************************************************************************
 *
 *	Description:
 *		Handle octrees
 *
 */

# ifndef OCTREE_H
# define	OCTREE_H

typedef struct oc {
	int		iNodeNum;
	int		iStatus;
	int		iDepth;
	VECTOR		vCorner;
	int             iAtoms;
	ATOM		*PaAtomList;
	float		*PfCharges;	/* float to save space */
	struct oc	*PonChildren;
} OCTNODEt;

typedef OCTNODEt	*OCTNODE;

typedef	struct	{
	int		iType;	
	int		iMaxDepth;
	int		iTreePoints;
	int		iDielectric;
	double		dGridSize;
	int		*PiDensities;
	float		*PfCharges;	/* float to save space */
	double		*PdHalfEdges;
	double		*PdHalfDiagonals;
	VARARRAY        vaAtoms;
	OCTNODEt	onHead;
} OCTREEt;

typedef OCTREEt		*OCTREE;

/* OCTREEt.iOctType types */

#define OCT_SHELL 		1
#define OCT_INTERIOR_SOLUTE 	2
#define OCT_INTERIOR_SOLVENT	3

/* dielectric options */

#define DIEL_R2		1
#define DIEL_R		2

/* types of printing */

#define COLOR_CUT	1
#define COLOR_RANGE	2
#define COLOR_DEPTH	3
#define COLOR_NONE	4

extern OCTREE	octOctTreeCreate(UNIT uUnit, int iType, double dGridSpace,
			double dAddExtent, double dShellExtent, bool bIncludeSolvent);
extern void	OctTreeDestroy( OCTREE *PoctTree );

extern void	OctTreeInitCharges(OCTREE octTree, int iDielectric, 
			double dCutDist, VECTOR *PvMin, VECTOR *PvMax);	

extern void	OctTreeDescribe();	/* ( OCTREE ) */
extern void	OctTreeDeleteSphere(OCTREE octTree, VECTOR *PvPoint, double dRadius);
extern void	OctTreeUpdateCharge(OCTREE octTree, VECTOR *PvNewPoint, 
			float fCharge, double dCutDist, VECTOR *PvMax, VECTOR *PvMin);
extern RESIDUE	rOctTreeCheckSolvent(OCTREE octTree, VECTOR *PvPoint);

extern void	OctTreePrintGrid(OCTREE octTree, char *sFileName, int iColor);

# endif
