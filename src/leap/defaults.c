#include <stdio.h>
#include "basics.h"
#include "defaults.h"
#include "classes.h"
#include "objekt.h"
#include "container.h"
#include "unit.h"
#include "residue.h"
#include "vector.h"
#include "octree.h" // DIEL_R*
#include "tools.h" // DEFAULT_DISTANCE_SEARCH

defaultstruct GDefaults;

typedef struct {
    int  iType; // B=bool, D=double, I=integer, S=string-enum
    char *sNameLower, *sName;
    void *pValue;
    union {
        double real;
        int integer;
    } defval; // integer for boolean, character, and string-enum types as well
    // direct string default (currently) is always empty. If this changes, may need: char *defval.string
    char **options;    // setting tokens
    char **optionDesc; // full descriptive name
} DefaultSetting;

static char opt_null[]="null";
char *PBRadii_options[] =
     {"bondi", "amber6", "mbondi", opt_null /*"pbamber"*/ ,opt_null,opt_null, "mbondi2", "parse", "mbondi3", NULL};
char *PBRadii_optionDesc[] =
     {"Bondi radii","Amber6 modified Bondi radii","Modified Bondi radii","","","",
       "H(N)-modified Bondi radii", "ArgH and AspGluO modified Bondi2 radii", NULL};
static char *Dielectric_options[] = {opt_null,"constant","distance",NULL};
static char *Dielectric_optionDesc[] = {"Undefined","Constant","Distance",NULL};
static char *PdbConvertResname_options[] = {"none","keep","variant","standard",NULL};
static char *PdbConvertResname_optionDesc[] = {"PdbResMap","Retain all names",
                                      "Variant names","Standard names",NULL };

// NOTE: uninitialzed values get zero, and most defaults are zero
static DefaultSetting
zSDefaultSettings[] = {
    { 'B', "pdbwritecharges", "PdbWriteCharges", &GDefaults.pdbwritecharges },
    { 'B', "nocenter", "NoCenter", &GDefaults.nocenter },
    { 'B', "reorder_residues", "Reorder_Residues", &GDefaults.reorder_residues, .defval.integer = 1 },
    { 'B', "reverse_lists", "Reverse_Lists", &GDefaults.reverse_lists },
    //{ 'B', "oldprmtopformat", "OldPrmtopFormat", &GDefaults.iOldPrmtopFormat },
    { 'B', "gibbs", "Gibbs", &GDefaults.iGibbs },
    { 'B', "hybrid36", "Hybrid36", &GDefaults.bPdbHybrid36, .defval.integer=1 },
    { 'B', "keep_chainid", "Keep_chainId", &GDefaults.bPdbKeepChainId },
    { 'I', "pdbreadbiomt", "PdbReadBioMT", &GDefaults.iPdbReadBioMT },
    { 'B', "charmm", "Charmm", &GDefaults.iCharmm },
    { 'B', "flexiblewater", "FlexibleWater", &GDefaults.iFlexibleWater },
    { 'B', "deleteextrapointangles", "DeleteExtraPointAngles", &GDefaults.iDeleteExtraPointAngles,
               .defval.integer=1 },
    { 'S', "pbradii", "PB Radii", &GDefaults.iGBparm, .defval.integer=2,
               .options = PBRadii_options, .optionDesc = PBRadii_optionDesc },
    { 'D', "searchdistance", "SearchDistance", &GDefaults.dDSearchDistance,
               .defval.real=DEFAULT_DISTANCE_SEARCH },
    { 'D', "gridspace", "GridSpace", &GDefaults.dGridSpace, .defval.real=1.0 },
    { 'D', "shellextent", "ShellExtent", &GDefaults.dShellExtent, .defval.real=4.0 },
    { 'D', "dipole_damp_factor", "Dipole Damping Factor", &GDefaults.dDipoleDampFactor },
    { 'D', "scee", "SCEE 1-4 Scale Factor", &GDefaults.dSceeScaleFactor, .defval.real=1.2 },
    { 'D', "scnb", "SCNB 1-4 Scale Factor", &GDefaults.dScnbScaleFactor, .defval.real=2.0 },
    { 'B', "cmap", "CMAP", &GDefaults.iCMAP },
    { 'I', "ipol", "IPOL", &GDefaults.iIPOL },
    { 'S', "dielectric", "Dielectric", &GDefaults.iDielectricFlag, .defval.integer=DIEL_R2,
               .options = Dielectric_options, .optionDesc = Dielectric_optionDesc },
    { 'B', "residueimpropers", "ResidueImpropers", &GDefaults.iResidueImpropers },
    { 'B', "pdb_auto_match","PDB_Auto_Match", &GDefaults.bPdbAutoMatch },
    { 'B', "pdb_auto_link","PDB_Auto_Link", &GDefaults.bPdbAutoLink },
    { 'D', "pdb_link_cutoff","PDB_Link_Cutoff", &GDefaults.dPdbLinkCovalentCutoff, .defval.real=1.2 },
    { 'D', "pdb_crosslink_cutoff","PDB_CrossLink_Cutoff", &GDefaults.dPdbCrosslinkCovalentCutoff,
               .defval.real=1.1 },
    { 'B', "pdb_auto_load","PDB_Auto_Load", &GDefaults.bPdbAutoLoadRes },
    { 'C', "pdb_altloc","PDB_AltLoc", &GDefaults.cPdbAltLocSelect, .defval.integer = 'A' },
    { 'B', "pdb_use_link","PDB_Use_LINK", &GDefaults.bPdbUseLinkRecords },
    { 'B', "pdb_use_conect","PDB_Use_CONECT", &GDefaults.bPdbUseConect },
    { 'B', "pdb_link_ions","PDB_Link_Ions", &GDefaults.bPdbLinkIons },
    { 'B', "pdb_reset_chainids","PDB_Reset_ChainIds", &GDefaults.bPdbResetChainID },
    { 'S', "pdb_convert_resname","PDB_Convert_ResName", &GDefaults.iPdbConvertResName,
               .options = PdbConvertResname_options, .optionDesc = PdbConvertResname_optionDesc  },
    { 'B', "pdb_ignore_nonconnect","PDB_Ignore_NonConnect", &GDefaults.iPdbIgnoreNonConnect },
    { 'I', "pdb_read_model","PDB_read_Model", &GDefaults.iPdbReadModel, .defval.integer=-1 },
    { 'B', "pdb_expand_biomt","PDB_Expand_BioMT", &GDefaults.bPdbExpandBioMt },
    { 'B', "pdb_expand_ncs","PDB_Expand_NCS", &GDefaults.bPdbExpandNCSMt, .defval.integer=1 },
    { 'B', "pdb_expand_symmetry","PDB_Expand_Symmetry", &GDefaults.bPdbExpandSymm },
    { 'B', "cif_read_auth","CIF_Read_auth", &GDefaults.bCIFReadAuth, .defval.integer=1 },
    { 'S', "pdb_patch_file","PDB_Patch_File", &GDefaults.sPdbPatchFilename },
    { 'B', "prmtop_netcdf","PrmTop_NetCDF", &GDefaults.bPrmtopNetcdf },
    { 0 }
};

void InitializeDefaults(void) {
    for (int i=0; zSDefaultSettings[i].sName; i++) {
        switch (zSDefaultSettings[i].iType) {
            case 'B':
                *((bool*)zSDefaultSettings[i].pValue) = zSDefaultSettings[i].defval.integer;
                break;
            case 'I':
                *((int*)zSDefaultSettings[i].pValue) = zSDefaultSettings[i].defval.integer;
                break;
            case 'C':
                *((char*)zSDefaultSettings[i].pValue) = zSDefaultSettings[i].defval.integer;
                break;
            case 'D':
                *((double*)zSDefaultSettings[i].pValue) = zSDefaultSettings[i].defval.real;
                break;
            case 'S':
                if (zSDefaultSettings[i].options)
                    *((int*)zSDefaultSettings[i].pValue) = zSDefaultSettings[i].defval.integer;
                else
                    (*((STRING*)zSDefaultSettings[i].pValue))[0] = 0;
                break;
        }
    }
}

OBJEKT
SetDefault(char *sParam, OBJEKT oValue) {
STRING sParamLower;
int iOptIndex;
        strcpy(sParamLower,sParam);
        StringLower(sParamLower);
        iOptIndex = -1;
        for (int i=0; zSDefaultSettings[i].sName; i++) {
            if (!strcmp(zSDefaultSettings[i].sNameLower,sParamLower) ||
                // special case for scee, scnb
                (sParamLower[0]=='s' && sParamLower[1]=='c' && !strncmp(zSDefaultSettings[i].sNameLower,sParamLower,4))) {
                iOptIndex=i;
                break;
            }
        }
        if (iOptIndex < 0) {
            VPFATAL("Default setting \"%s\" not found.\n",sParam);
            return NULL;
        }
        int iType = iObjectType(oValue);
        double dValue = ( iType == ODOUBLEid ) ? dODouble(oValue) : 0.0;
        double intpart;
        modf(dValue, &intpart);
        int iValue = (int)intpart;
        int iOptionValue = -1;
        char *sValue = NULL;
        bool bValue = FALSE;
        if (iType == OSTRINGid) {
            sValue = sOString(oValue);
            if (!strcmp("?",sValue)) {
                OBJEKT oResult;
                switch (zSDefaultSettings[iOptIndex].iType) {
                case 'B':
                    dValue = (double) *((bool*)zSDefaultSettings[iOptIndex].pValue);
                    break;
                case 'I':
                    dValue = (double) *((int*)zSDefaultSettings[iOptIndex].pValue);
                    break;
                case 'D':
                    dValue = *((double*)zSDefaultSettings[iOptIndex].pValue);
                    break;
                case 'C':
                    oResult = (OBJEKT) oCreate(OSTRINGid);
                    {
                        char str[2];
                        str[0] = *((char*)zSDefaultSettings[iOptIndex].pValue);
                        str[1] = 0;
                        OStringDefine( (OSTRING)oResult, strdup(str) );
                    }
                    return oResult;
                    break;
                case 'S':
                    iValue = *((int*)zSDefaultSettings[iOptIndex].pValue);
                    oResult = (OBJEKT) oCreate(OSTRINGid);
                    OStringDefine( (OSTRING)oResult, zSDefaultSettings[iOptIndex].options[iValue] );
                    return oResult;
                    break;
                }
                oResult = (OBJEKT) oCreate(ODOUBLEid);
                ODoubleSet( (ODOUBLE)oResult, dValue );
                return oResult;
            }
        }

        // Special case: IPOL can only be set once
        if ( !strcmp( sParamLower,"ipol") && GDefaults.iIPOLset > 0 ) {
            VPWARN("IPOL has already been set to %i in frcmod/parm.dat.\n",
                            GDefaults.iIPOL);
            VP0("Please change the setting in frcmod/parm.dat.\n");
            return(NULL);
        }


        switch (zSDefaultSettings[iOptIndex].iType) {
        case 'B':
            if ( iType == OSTRINGid && (!strcasecmp( sValue, "on" ) || !strcasecmp(sValue,"true")))
                bValue = TRUE;
            else if ( iType == OSTRINGid && (!strcasecmp( sValue, "off" ) || !strcasecmp(sValue,"false")))
                bValue = FALSE;
            else if ( iType == ODOUBLEid && dValue != 0.0 )
                bValue = TRUE;
            else if ( iType == ODOUBLEid && dValue == 0.0 )
                bValue = FALSE;
            else {
                VPFATAL("Set %s: value must be 'on'/'true'/1 or 'off'/'false'/0\n",zSDefaultSettings[iOptIndex].sName);
                return NULL;
            }
            *((bool*)zSDefaultSettings[iOptIndex].pValue) = bValue;
            break;
        case 'D':
            if ( iType != ODOUBLEid ) {
                VPFATAL("Set %s: value must be a number\n",zSDefaultSettings[iOptIndex].sName);
                return NULL;
            }
            if (dValue < 0.0) {
                dValue = zSDefaultSettings[iOptIndex].defval.real;
                VPWARN("set default: %s value must be greater than 0; resetting to %g\n",
                       zSDefaultSettings[iOptIndex].sName, dValue);
            }
            printf("Set value = %g\n",dValue);
            *((double*)zSDefaultSettings[iOptIndex].pValue) = dValue;
            break;
        case 'I':
            if ( iType != ODOUBLEid || dValue != (double)iValue ) {
                VPFATAL("Set %s: value must be an integer\n",zSDefaultSettings[iOptIndex].sName);
                return NULL;
            }
            if ( !strcmp( sParamLower, "ipol") && ( iValue < 0 || iValue > 4 ) ) {
                 iValue = zSDefaultSettings[iOptIndex].defval.integer;
                 VPWARN("Only IPOL = 0 to 4 is supported; resetting IPOL to %d.\n",iValue);
            }
            *((int*)zSDefaultSettings[iOptIndex].pValue) = iValue;
            break;
        case 'C':
            if ( iType != OSTRINGid || (iType == OSTRINGid && strlen(sValue) != 1) ) {
                VPFATAL("Set %s: value must be a single character\n",zSDefaultSettings[iOptIndex].sName);
                return NULL;
            }
            *((char*)zSDefaultSettings[iOptIndex].pValue) = sValue[0];
            break;
        case 'S':
            if ( iType != OSTRINGid) {
                VPFATAL("Set %s: value must be an string\n",zSDefaultSettings[iOptIndex].sName);
                return NULL;
            }
            // Can be raw string or named enum
            if (zSDefaultSettings[iOptIndex].options) {
                for (int i=0;zSDefaultSettings[iOptIndex].options[i];i++) {
                    if (!strcasecmp(zSDefaultSettings[iOptIndex].options[i],sValue)) {
                        iOptionValue = i;
                        break;
                    }
                }
                if (iOptionValue < 0) {
                    STRING sLine;
                    sprintf(sLine,"Set %s: value must be one of '%s'",
                            zSDefaultSettings[iOptIndex].sName, zSDefaultSettings[iOptIndex].options[0]);
                    for (int i=1;zSDefaultSettings[iOptIndex].options[i];i++) {
                        sprintf(sLine+strlen(sLine),", '%s'",zSDefaultSettings[iOptIndex].options[i]);
                    }
                    VPFATAL("%s",sLine);
                    return NULL;
                }
                *((int*)zSDefaultSettings[iOptIndex].pValue) = iOptionValue;
                VP0("%s: Using %s\n", zSDefaultSettings[iOptIndex].sName,
                        zSDefaultSettings[iOptIndex].optionDesc[iOptionValue]);
            } else {
                strcpy(*((STRING*)zSDefaultSettings[iOptIndex].pValue), sValue);
            }
            break;
        }
        return(NULL);
}

void
zPrintDefaultSetting(int index) {
    VP0("%25s :    ",zSDefaultSettings[index].sName);
    if (zSDefaultSettings[index].iType == 'B') {
        VP0("%s\n", *((bool*)zSDefaultSettings[index].pValue) ?  "on/true" : "off/false");
    } else if (zSDefaultSettings[index].iType == 'I') {
        VP0("%i\n",*((int*)zSDefaultSettings[index].pValue));
    } else if (zSDefaultSettings[index].iType == 'S') {
        if (!zSDefaultSettings[index].optionDesc) {
            VP0("\"%s\"\n", (char*)zSDefaultSettings[index].pValue);
            return;
        }
        int ioption = *((int*)zSDefaultSettings[index].pValue);
        // We loop only because list length is defined by NULL
        for (int i=0;zSDefaultSettings[index].optionDesc[i]; i++) {
            if (ioption == i) {
                VP0("%s\n",zSDefaultSettings[index].optionDesc[ioption] );
                return;
            }
        }
        VP0("ERROR: Undfined otion index: %d\n",ioption );
    } else if (zSDefaultSettings[index].iType == 'D') {
        VP0("%lf\n",*((double*)zSDefaultSettings[index].pValue) );
    } else if (zSDefaultSettings[index].iType == 'C') {
        char c = *((char*)zSDefaultSettings[index].pValue);
        if (c) VP0("'%c'\n",*((char*)zSDefaultSettings[index].pValue) );
        else VP0("''\n");
    }
}


void
PrintDefaultSettings(char *cPParam) {
    if (cPParam) {
        for (int i=0;zSDefaultSettings[i].sName;i++) {
            if (!strcasecmp(zSDefaultSettings[i].sName,cPParam)) {
                zPrintDefaultSetting(i);
                return;
            }
        }
        VPWARN("Default parameter %s not found\n", cPParam );
        return;
    }
    for (int i=0;zSDefaultSettings[i].sName;i++)
        zPrintDefaultSetting(i);
}

