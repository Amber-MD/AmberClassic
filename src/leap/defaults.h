#ifndef DEFAULTS_H
# define DEFAULTS_H

typedef struct {
	double	dDSearchDistance;
	double	dDESPGridSpace;
	double	dDESPBoxSize;
	float	fDESPDielectric;
	int	iDESPConstant;
	int	pdbwritecharges;
	int nocenter;
	int reorder_residues;
	int reverse_lists;
	double	dGridSpace;
	double	dShellExtent;
	int	iDielectricFlag;
	int	iGBparm;
	int	iGibbs;
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

        double  dPdbLinkCovalentCutoff;
        double  dPdbCrosslinkCovalentCutoff;
        int     iPdbIgnoreNonConnect;
        int     iPdbReadModel;
        int     iPdbConvertResName;
        BOOL    bPdbAutoMatch;
        BOOL    bPdbAutoLink;
        BOOL    bPdbAutoLoadRes;
        BOOL    bPdbUseLinkRecords;
        BOOL    bPdbUseConect;
        BOOL    bPdbLinkIons;
        BOOL    bPdbResetChainID;
        BOOL    bPdbExpandBioMt;
        BOOL    bPdbExpandNCSMt;
        BOOL    bPdbExpandSymm;
        BOOL    bCIFReadAuth;
        char    cPdbAltLocSelect;
        STRING  sPdbPatchFilename;
} defaultstruct ;

extern defaultstruct GDefaults;

#endif
