
typedef struct {
	double	dDSearchDistance;
	double	dDESPGridSpace;
	double	dDESPBoxSize;
	float	fDESPDielectric;
	int	iDESPConstant;
	int	pdbwritecharges;
	int nocenter;
	int reorder_residues;
	double	dGridSpace;
	double	dShellExtent;
	int	iDielectricFlag;
	int	iGBparm;
	int 	iOldPrmtopFormat;
	int		iGibbs;
	int 	iCharmm; 
	int	iResidueImpropers;
	int	iDeleteExtraPointAngles;
	int     bPdbHybrid36;
	int     bPdbKeepChainId;
	int     iPdbReadBioMT;
	int     iFlexibleWater;
	double  dDipoleDampFactor;
	double  dSceeScaleFactor;
	double  dScnbScaleFactor;
	int     iCMAP;
	int     iIPOL;    
	int     iIPOLset;    /* indicate IPOL set in frcmod/parm.dat */
} defaultstruct ;

extern defaultstruct GDefaults;

