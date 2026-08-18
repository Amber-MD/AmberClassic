#ifndef DEFAULTS_H
# define DEFAULTS_H
#include "objekt.h"

typedef struct {
	bool	bCompatible;
	double	dDSearchDistance;
	double	dDESPGridSpace;
	double	dDESPBoxSize;
	float	fDESPDielectric;
	int	iDESPConstant;
	bool	pdbwritecharges;
	bool    bNoCenter;
	bool    bReorderResidues;
	bool    bReorderMolecules;
	bool    bReverseLists;
	bool    bOrigCMAPOrder;
        bool    bRandomOrientation;
        bool    bOldPrmtopFormat;
        double  dPrmtopFormat;
	double	dGridSpace;
	double	dShellExtent;
	int	iDielectricFlag;
	int	iGBparm;
	bool	bGibbs;
	bool 	bCharmm; 
	bool	bResidueImpropers;
	bool	bDeleteExtraPointAngles;
	bool    bPdbHybrid36;
        bool    bKeepInputSolvent;
	bool    bPdbKeepChainId;
	int     iPdbReadBioMT;
	bool    bFlexibleWater;
	double  dDipoleDampFactor;
	double  dSceeScaleFactor;
	double  dScnbScaleFactor;
	bool    bCMAP;
	int     iIPOL;    
	int     iIPOLset;    /* indicate IPOL set in frcmod/parm.dat */

        double  dPdbLinkCovalentCutoff;
        double  dPdbCrosslinkCovalentCutoff;
        double  dPdbBulkQ;
        int     iPdbIgnoreNonConnect;
        int     iPdbReadModel;
        int     iPdbConvertResName;
        bool    bPrmtopNetcdf;
        bool    bPdbAutoMatch;
        bool    bPdbExactMatch;
        bool    bPdbAutoLink;
        bool    bPdbAutoLoadRes;
        bool    bPdbUseLinkRecords;
        bool    bPdbUseConect;
        bool    bPdbLinkIons;
        bool    bPdbResetChainID;
        bool    bPdbExpandNCSMt;
        bool    bPdbExpandSymm;
        bool    bCIFReadAuth;
        bool    bMaskPDBMode;
        bool    bTiming;
        char    cPdbAltLocSelect;
        STRING  sPdbPatchFilename;
} defaultstruct ;

extern defaultstruct GDefaults;

extern void InitializeDefaults(void);
extern OBJEKT SetDefault(char *sParam, OBJEKT oValue);
extern void PrintDefaultSettings(char *cPParam);

extern char *PBRadii_options[], *PBRadii_optionDesc[];

#define GBPARM_BONDI 0
#define GBPARM_AMBER6 1
#define GBPARM_MBONDI 2
#define GBPARM_MBONDI2 3
#define GBPARM_PARSE 4
#define GBPARM_MBONDI3 5
#define GBPARM_MAX 5

#endif
