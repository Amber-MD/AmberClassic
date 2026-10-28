/*
 *      File:   amber.c
 *
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
 *      Description:
 *              Read AMBER prep and parameter database files into
 *              UNITs and PARMSETs.
 */

#include        "basics.h"
#include        "vector.h"
#include        "classes.h"
#include        "leap.h"
#include        "tools.h"
#include        "amber.h"
#include        "defaults.h"
#include        "sort.h"
#include        "zMatrix.h"
#include        "cmap.h"

#define FGETS(s,f)      { (s)[0]=0; fgets(s,sizeof(s),f); }

/*
 *      Mapping of AMBER types to ELEMENT names used by iAtomTypesToElement
 */

typedef struct  {
        char    sType[ATOMTYPELEN];
        int     iElement;
        int     iHybridization;
} TYPEELEMENTt;

static VARARRAY         vaAtomTypes = NULL;

// NONBONDt and MASSt are merged and passed to iParmSetAddAtom()
// and get storderd as ATOMPARMt
// Loop over NONBONDt, and skipped if MASSt absent.
typedef struct  {
        char            sNonBondType[ATOMNAMELEN];
        double          dRStar;
        double          dDepth;
        double          dRStar14;
        double          dDepth14;
        DESCRIPTION     sDesc;
} NONBONDt;

typedef struct  {
        char            sType[ATOMTYPELEN];
        double          dMass;
        double          dPolar;
        double          dScreenF;
        int             iElement;
        int             iHybridization;
} MASSt;

#define MAXEQUIV        20
typedef char    SHORTt[ATOMNAMELEN];
typedef struct  {
        char            sName[ATOMNAMELEN];
        int             iEquivs;
        SHORTt          saEquivs[MAXEQUIV];
} EQUIVt;

/*
 *----------------------------------------------------
 *
 *      Private routines
 */

/*
 * These are intended to support pseudo free formating by removing dashes from
 * CO-CX-2C-C    1     0.00000    0.0 -4.0  CASH/CASN chi1'
 * The format is at most a2,1x,a2,1x,a2,1x,a2; hence, the magic 9.
 */
#define REPLACEDASHES(s) {StringRTrim(sLine); for(int z=0;z<9;z++) if(s[z]=='-')s[z]=' '; else if(s[z]=='\n') break;}
#define NODASHES(s)     {for(int z=0;z<9;z++) if(s[z]=='-')s[z]=' ';}
#define NODASHESWL(s)     {for(int z=0;z<strlen(s);z++) if(s[z]=='-')s[z]=' ';}


/*
 *      zAmberConvertWildCard
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      If the atom type is an "X", an old AMBER wild card
 *      character then convert it to the new wild card character '?'.
 */
static inline void
zAmberConvertWildCard( char *sType )
{
    if ( strcmp( sType, AMBER_WILD_CARD ) == 0 ) {
        strcpy( sType, WILD_CARD_TYPE );
    }
}



/*
 *      zAmberReadMasses
 *
 *      Read the atomic masses into a VARARRAY
 */
static void
zAmberReadParmSetMasses( VARARRAY *vaPMasses, FILE *fIn )
{
STRING          sLine;
int             iRead;
VARARRAY        vaMasses;
MASSt           mMass;
double          dPolar;
double          dScreenF;
int             iElement;
STRING          sHybridization;

                /* Read the atom masses into a VARARRAY until we */
                /* read the non-bond parameters */

    memset(&mMass, 0, sizeof(mMass));                   /* for Purify */
    MESSAGE("Reading MASSes.\n" );
    vaMasses = vaVarArrayCreate( sizeof(MASSt) );

    while (1) {
        FGETS( sLine, fIn );
        iRead = sscanf( sLine, "%s %lf %lf %lf %d %s",
                        mMass.sType, &mMass.dMass, &dPolar, &dScreenF,
                        &iElement, sHybridization );
        if (iRead <= 0)
            break;
        MESSAGE("Read %d MASS terms from: %s", iRead, sLine );
        if ( iRead == 1 ) {
            VPWARN("%s: mass not read - omitting\n", mMass.sType );
            continue;
        }
        if ( mMass.dMass < 0.0 ) {
            VPWARN("%s: mass %lf - omitting\n", mMass.sType, mMass.dMass );
            continue;
        }
        mMass.dPolar = -1.0;    /* flag that default=0.0 will be used */
        mMass.dScreenF = 0.0;
        mMass.iElement=-10;   /* NOELEMENT==EP==-1 ; -10 seems safe */
        mMass.iHybridization=HUNDEFINED;
        switch( iRead ){
            default:
            case 6:
                StringLower( sHybridization );
                if ( !strcmp( sHybridization, "sp2") )
                    mMass.iHybridization = HSP2;
                else if ( !strcmp( sHybridization, "sp3") )
                    mMass.iHybridization = HSP3;
                else if ( !strcmp( sHybridization, "sp") )
                    mMass.iHybridization = HSP1;
                else if ( !strcmp( sHybridization, "sp0") )
                    mMass.iHybridization = 0;
                else {
                    mMass.iHybridization = HUNDEFINED;
                    VPWARN("atom type %s - unknown hybridization %s\n",
                             mMass.sType, sHybridization );
                }
                /* fall through */
            case 5:
                mMass.iElement = iElement;
                /* fall through */
            case 4:
                mMass.dScreenF = dScreenF;
                /* fall through */
            case 3:
                if ( dPolar < 0.0 ) {
                    if ( dPolar < -0.99 && dPolar > -1.01 )
                        break;
                    VPWARN("%s: polarization %lf - omitting\n",
                                                mMass.sType, dPolar );
                    continue;
                }
                if ( dPolar > 15.0 ){
                    VPWARN("%s: polarization %lf\n", mMass.sType, dPolar );
                }
                mMass.dPolar = dPolar;
                /* fall through */
            case 2:
                break;
        };

        VarArrayAdd( vaMasses, (GENP)&mMass );
    }
    MESSAGE("Done reading MASSes.\n" );
    *vaPMasses = vaMasses;
}





//---------------------------------------------------------------------------------------------
// zAmberReadParmSetBonds: read the parameters for bonds between atom types from an Amber
//                         format parameter file (i.e. parm10.dat).
//
// Arguments:
//   psParms:    the developing parameter set
//   fIn:        the input Amber parameters file
//---------------------------------------------------------------------------------------------
static void zAmberReadParmSetBonds(PARMSET psParms, FILE *fIn)
{
    int    iRead, iDesc1, iDesc;
    STRING sLine;
    STRING saStr[2];
    double dKb, dR0, dKpull, dRpull0, dKpress, dRpress0;

    // For Purify
    memset(saStr, 0, sizeof(saStr));
    MESSAGE("Reading BONDs.\n" );
    while (1) {
        FGETS(sLine, fIn);
        REPLACEDASHES(sLine);
        StringRTrim(sLine);
        iRead = sscanf(sLine, "%s %s %lf %lf %n%lf %lf %lf %lf %n", saStr[0], saStr[1],
                       &dKb, &dR0, &iDesc1, &dKpull, &dRpull0, &dKpress, &dRpress0, &iDesc);
        if (iRead < 4) // XXX was <=0
            break;
        if (iRead != 8) {
            iDesc = iDesc1;
            dKpull   =   0.0;
            dKpress  =   0.0;
            dRpull0  = 100.0;
            dRpress0 =   0.0;
        }
        MESSAGE("Read %d BOND terms from: %s", iRead, sLine );
        iParmSetAddBond(psParms, saStr[0], saStr[1], dKb, dR0, dKpull, dRpull0,
                        dKpress, dRpress0, &sLine[iDesc]);
    }
    MESSAGE("Done reading BONDs.\n" );
}

/*
 *      zAmberReadParmSetAngles
 *
 *      Read the angle parameter terms.
 */
static void
zAmberReadParmSetAngles( PARMSET psParms, FILE *fIn )
{
STRING          sLine;
int             iRead, iDesc;
STRING          saStr[3];
double          dKt, dT0, dTkub, dRkub, zero;

    zero = 0.0;
    memset(saStr, 0, sizeof(saStr));                /* for Purify */
    MESSAGE("Reading ANGLEs.\n" );
    while (1) {
        FGETS( sLine, fIn );
        REPLACEDASHES(sLine);
        if( GDefaults.bCharmm ){
            iRead = sscanf( sLine, "%s %s %s %lf %lf %lf %lf %n",
                            saStr[0], saStr[1], saStr[2], &dKt, &dT0, &dTkub, &dRkub, &iDesc );
            if ( iRead != 7 )
                break;
            MESSAGE("Read 7 ANGLE terms from: %s", sLine );
            iParmSetAddAngle( psParms, saStr[0], saStr[1], saStr[2],
                              dKt, dT0*DEGTORAD, dTkub, dRkub, &sLine[iDesc] );
        } else {
            iRead = sscanf( sLine, "%s %s %s %lf %lf %n",
                            saStr[0], saStr[1], saStr[2], &dKt, &dT0, &iDesc );
            if ( iRead != 5 )
                break;
            MESSAGE("Read 5 ANGLE terms from: %s", sLine );
            iParmSetAddAngle( psParms, saStr[0], saStr[1], saStr[2],
                              dKt, dT0*DEGTORAD, zero, zero, &sLine[iDesc] );
        }
    }
    if ( iRead > 0 ){
        VPWARN("Incomplete Angle line:\n%s", sLine );
    }
    MESSAGE("Done reading ANGLEs.\n" );
}

static void
zAmberReadParmSetIPOL( FILE *fIn )
{
STRING          sLine;
int             iRead,myiPOL;
    while (1) {
        FGETS( sLine, fIn );
        NODASHES(sLine);
        iRead = sscanf( sLine, "%i", &myiPOL );
        if ( iRead != 1 ) break;
        MESSAGE("Read: %s\n", sLine );
        if (GDefaults.iIPOLset == 1 && GDefaults.iIPOL != myiPOL ) {
           VPWARN("Conflicting IPOL. %i %i\n", myiPOL, GDefaults.iIPOL);
        } else {
           GDefaults.iIPOL = myiPOL;
           GDefaults.iIPOLset = 1;
        }
    }
}

/*
 *  zAmberReadParmSetCMAP
 *
 *  Read CMAP parameter records from frcmod/dat file into psParms.
 *
 *  Token ordering rules:
 *    CMAP_TITLE     — always first in a record; commits any prior record
 *    CMAP_PARAMETER — always last; commits current record, resets bReading
 *    CMAP_COUNT     — either total count or per-record 1-based index;
 *                     validated as max, not exact match
 *    All others     — may appear in any order between TITLE and PARAMETER
 *    Any token other than COUNT or TITLE before bReading is set → warning
 */
static void
zAmberReadParmSetCMAP(PARMSET psParms, FILE *fIn)
{
    STRING sLine;
    char tok1[32], tok2[32];
    int ipos;
    int iMaxCount = 0;       /* max seen in CMAP_COUNT */
    bool bReading = false;   /* true between CMAP_TITLE and CMAP_PARAMETER */
    CMAPt cmap;

    while (1) {
        FGETS(sLine, fIn);
        if (sLine[0] == '\0') break;
        StringRTrim(sLine);
        if (sLine[0] != '%') {
            if (strlen(sLine))
                VPWARN("Unrecognized line parsing CMAP data: %s\n", sLine);
            continue;
        }
        sscanf(sLine, "%31s%31s%n", tok1, tok2, &ipos);

        if (strcmp("%FLAG", tok1) != 0) {
            if (strcmp("%COMMENT", tok1) == 0) {
                // No bReading check: COMMENT lines preceding TITLE drop into the next cmap
                /* accumulate comment into current record */
                char *cmnt = &sLine[9];
                size_t used = strnlen(cmap.comments, sizeof(cmap.comments));
                if (used < sizeof(cmap.comments) - 1) {
                    if (used > 0) cmap.comments[used++] = '\n';
                    size_t avail = sizeof(cmap.comments) - used - 1;
                    size_t n = strnlen(cmnt, avail);
                    memcpy(cmap.comments + used, cmnt, n);
                    cmap.comments[used + n] = '\0';
                }
            }
            else VPWARN("Unrecognized line parsing CMAP data: %s\n", sLine);
            continue;
        }

        /* Invalid record guard: only COUNT and TITLE allowed before bReading ── */
        if (!bReading &&
            strcmp("CMAP_COUNT", tok2) != 0 &&
            strcmp("CMAP_TITLE", tok2) != 0) {
            VPWARN("CMAP %%FLAG %s before %%FLAG CMAP_TITLE — skipping\n", tok2);
            continue;
        }

        if (strcmp("CMAP_COUNT", tok2) == 0) {
            int iCount = atoi(&sLine[ipos]);
            if (iCount > iMaxCount) iMaxCount = iCount;
            if (iCount > iParmSetTotalCMAPParms(psParms) + 1)
                VPWARN("CMAP_COUNT=%d but only %d records read so far\n",
                       iCount, iParmSetTotalCMAPParms(psParms));

        } else if (strcmp("CMAP_TITLE", tok2) == 0) {
            /* PARAMETER must have closed the previous record */
            if (bReading)
                VPWARN("CMAP record '%s' had no CMAP_PARAMETER — discarding\n",
                       cmap.title);
            /* discard incomplete record — free any allocated inner pointers */
            if (bReading && cmap.reslist) { free(cmap.reslist); cmap.reslist = NULL; }
            if (bReading && cmap.map)     { free(cmap.map);     cmap.map     = NULL; }



            memset(&cmap, 0, sizeof(CMAPt));
            bReading = true;

            FGETS(sLine, fIn);
            StringTrim(sLine);
            sLine[sizeof(cmap.title) - 1] = '\0';
            strcpy(cmap.title, sLine);
            VP2("Reading CMAP: %s\n", cmap.title);

            /* defaults — standard backbone phi/psi */
            cmap.residx[0] = -1;
            cmap.residx[4] =  1;
            cmap.iRefResIndex = 2; // default reference residue
            strcpy(cmap.atmname[0], "C");
            strcpy(cmap.atmname[1], "N");
            strcpy(cmap.atmname[2], "CA");
            strcpy(cmap.atmname[3], "C");
            strcpy(cmap.atmname[4], "N");

        } else if (strcmp("CMAP_RESLIST", tok2) == 0) {
            if (cmap.reslist) free(cmap.reslist);
            cmap.nres = atoi(&sLine[ipos]);
            cmap.reslist = malloc(sizeof(WRD) * cmap.nres);
            for (int i = 0; i < cmap.nres; i++)
                fscanf(fIn, "%7s", cmap.reslist[i]);

        } else if (strcmp("CMAP_RESOLUTION", tok2) == 0) {
            cmap.resolution = atoi(&sLine[ipos]);

        } else if (strcmp("CMAP_PARAMETER", tok2) == 0) {
            if (cmap.resolution <= 0) {
                VPWARN("CMAP_PARAMETER before CMAP_RESOLUTION in '%s' — skipping\n",
                       cmap.title);
            } else {
                if (cmap.map) {
                    VPWARN("CMAP record '%s': duplicate PARAMETER\n", cmap.title);
                    free(cmap.map);
                }
                cmap.map = malloc(sizeof(double) *
                              cmap.resolution * cmap.resolution);
                int m = cmap.resolution * cmap.resolution;
                for (int i = 0; i < m; i++)
                    fscanf(fIn, "%lf", &cmap.map[i]);
            }
            /* PARAMETER is always last — commit and reset */
            if (cmap.resolution <= 0 || !cmap.map)
                VPWARN("CMAP record '%s': missing map data — skipping\n", cmap.title);
            else
                iParmSetAddCMAP(psParms, &cmap);
            bReading = false;
        } else if (strcmp("CMAP_ATMLIST", tok2) == 0) {
            NODASHESWL(sLine);
            sscanf(&sLine[ipos], "%7s%7s%7s%7s%7s",
                   cmap.atmname[0], cmap.atmname[1], cmap.atmname[2],
                   cmap.atmname[3], cmap.atmname[4]);

        } else if (strcmp("CMAP_RESIDX", tok2) == 0) {
            NODASHESWL(sLine);
            sscanf(&sLine[ipos], "%d%d%d%d%d",
                   &cmap.residx[0], &cmap.residx[1], &cmap.residx[2],
                   &cmap.residx[3], &cmap.residx[4]);
            for (int i=0;i<4;i++)
                if (!cmap.residx[i]) { cmap.iRefResIndex = i; break; }

        } else if (strcmp("CMAP_TERMMAP", tok2) == 0) {
            int iTermMap = atoi(&sLine[ipos]);
            if (iTermMap != 0)
                VPFATALEXIT("CMAP_TERMMAP=%d not supported\n", iTermMap);

        } else {
            VPWARN("Unknown CMAP %%FLAG: %s\n", tok2);
        }
    }

    /* end of file — if still reading, CMAP_PARAMETER was never seen */
    if (bReading) {
        VPWARN("CMAP record '%s' had no CMAP_PARAMETER — discarding\n",
               cmap.title);
        if (cmap.reslist) free(cmap.reslist);
        if (cmap.map)     free(cmap.map);
    }
}

/*
 *      zAmberReadParmSetPropers
 *
 *      Read the proper torsion parameter terms.
 */
static void
zAmberReadParmSetPropers( PARMSET psParms, FILE *fIn )
{
    STRING          sLine;
    int             iRead, iN, iDesc, iDesc1, iDesc2;
    STRING          saStr[4], taStr[4];
    double          dDivisions, dKp, dP0, dN, dScEE, dScNB;
    char            *cScEE, *cScNB;

    memset(saStr, 0, sizeof(saStr));                    /* for Purify */
    memset(taStr, 0, sizeof(taStr));
    MESSAGE("Reading DIHEDRALs.\n" );
    while (1) {
        FGETS( sLine, fIn );
        REPLACEDASHES(sLine);

// beginning of a group.
/* Read in parmxx.dat format
 *
 *        DI-DJ-DK-DL division kp PHI periodicity SCEE SCNB
 */
        iRead = sscanf( sLine, "%s %s %s %s %lf %lf %lf %lf %n%lf %n%lf %n",
                                saStr[0], saStr[1], saStr[2], saStr[3],
                                &dDivisions, &dKp, &dP0, &dN, &iDesc1,
                                &dScEE, &iDesc2, &dScNB, &iDesc );
        if ( iRead <= 0 )
            break;
        if ( sLine[0] == ' ' && sLine[1] == ' ') {
            VPWARN("Bad proper torsion definition (skipping):\n (%s)\n",
                                                                sLine );
            continue;
        }
/* Try GLYCAM format
 *
 *        DI-DJ-DK-DL division kp PHI periodicity SCEE=scee SCNB=scnb
 */
        if (iRead < 9) {
            cScEE = strstr(sLine, "SCEE");
            if(cScEE!=NULL) {
                iRead += sscanf( cScEE, "SCEE=%lf", &dScEE);
            }
            cScNB = strstr(sLine, "SCNB");
            if(cScNB!=NULL) {
                iRead += sscanf( cScNB, "SCNB=%lf", &dScNB);
            }
        }

/*
 *  Neither SCEE nor SCNB has been set. Mark dScEE = 0.0, dScNB = 0.0 "not read".
 */
        if (iRead < 9) {
            iDesc = iDesc1;
            dScEE = -1.0;
            dScNB = -1.0;
        } else if (iRead < 10) {
            iDesc = iDesc2;
            dScNB = -1.0;
        }

/* For some reason, we must divide dKp by dDivisions */

        if ( dDivisions == 0.0 )
            dDivisions = 1.0;
        dKp /= dDivisions;
        iN = (int)floor(dN+0.5);
        zAmberConvertWildCard( saStr[0] );
        zAmberConvertWildCard( saStr[1] ); // FIXME: can these be wildcard in proper?
        zAmberConvertWildCard( saStr[2] ); // FIXME: can these be wildcard ever?
        zAmberConvertWildCard( saStr[3] );

        MESSAGE("Read %d DIHEDRAL terms from: %s", iRead, sLine );
        iParmSetAddProperTerm( psParms,
                                saStr[0], saStr[1], saStr[2], saStr[3],
                                abs(iN), dKp, dP0*DEGTORAD, dScEE,
                                dScNB, &sLine[iDesc] );
        while( iN < 0 ) {
            FGETS( sLine, fIn );
            NODASHES(sLine);

// The same group of torsions share identical SCEE, SCNB

            iRead = sscanf( &sLine[11], "%lf %lf %lf %lf",
                            &dDivisions, &dKp, &dP0, &dN);
            if ( iRead<=0 ) break;

            if ( dDivisions == 0.0 )
                dDivisions = 1.0;
            dKp /= dDivisions;
            iN = (int)floor(dN+0.5);
            MESSAGE("Read %d DIHEDRAL extra terms from: %s", iRead, sLine );
            iParmSetAddProperTerm( psParms, saStr[0], saStr[1], saStr[2], saStr[3],
                                   abs(iN), dKp, dP0*DEGTORAD, dScEE, dScNB, "MULT" );
        }
        if ( iRead <= 0 )
            break;
    }
    MESSAGE("Done reading DIHEDRALs.\n" );
}


/*
 *      zAmberReadParmSetImpropers
 *
 *      Read the improper torsion parameter terms.
 */
static void
zAmberReadParmSetImpropers( PARMSET psParms, FILE *fIn )
{
    STRING          sLine;
    int             iRead, iN, iDesc;
    STRING          saStr[4];
    double          dKp, dP0, dN, dScEE, dScNB;
    bool            bPrintLine;

    memset(saStr, 0, sizeof(saStr));                    /* for Purify */
    MESSAGE("Reading IMPROPERs.\n" );
    while (1) {
        FGETS( sLine, fIn );
        REPLACEDASHES(sLine);
        /*
         *  get atoms & values, skipping possible
         *      extraneous idivf value allowed by
         *      the input spec for PARM
         */
        iRead = sscanf( sLine, "%s %s %s %s %n",
                                saStr[0], saStr[1], saStr[2], saStr[3], &iDesc );
        if ( iRead != 4 )
            break;
        MESSAGE("Read %d IMPROPER strings;\n", iRead );
        iRead = sscanf( sLine+15, "%lf %lf %lf", &dKp, &dP0, &dN );
        if ( iRead != 3 )
            break;
        MESSAGE("read %d IMPROPER floats from: %s", iRead, sLine );
        dN = (int)floor(dN+0.5);
        zAmberConvertWildCard( saStr[0] );
        zAmberConvertWildCard( saStr[1] );
        zAmberConvertWildCard( saStr[2] );
        zAmberConvertWildCard( saStr[3] );
        iN = (int)dN;
        dScEE = 0.0;
        dScNB = 0.0;
        /*
         *  check everything in case a format or other user error
         *      led to wrong values (e.g. IDIVF offset)
         */
        bPrintLine = false;
        if ( dKp < 0.0 ) {
            VPWARN("Expected Improper Torsion PK>=0 (%f)\n", dKp );
            bPrintLine = true;
        }
        if ( dP0 < 179.999  ||  dP0 > 180.001 ) {
            VPWARN("Expected Improper Torsion PHASE=180 (%f)\n", dP0 );
            bPrintLine = true;
        }
        if ( iN < 1  ||  iN > 6 ) {
            VPWARN("Unexpected Improper Torsion PN term (%d)\n", iN );
            bPrintLine = true;
        }
        if ( bPrintLine )
            VP0("Here is the Improper Torsion line in question:\n%s", sLine );

        iParmSetAddImproperTerm( psParms,
                       saStr[0], saStr[1], saStr[2], saStr[3],
                       iN, dKp, dP0*DEGTORAD, dScEE, dScNB, &sLine[iDesc] );
    }
    if ( iRead > 0 ){
        VPWARN("Incomplete Improper Torsion line:\n%s", sLine );
    }
    MESSAGE("Done reading IMPROPERs.\n" );
}




/*
 *      zAmberReadParmSetHBonds
 *
 *      Read hbond parameters.
 */
static void
zAmberReadParmSetHBonds( PARMSET psParms, FILE *fIn )
{
STRING          sLine;
int             iRead, iDesc;
STRING          saStr[2];
double          dA, dB;

    memset(saStr, 0, sizeof(saStr));                    /* for Purify */
    MESSAGE("Reading H-BONDs.\n" );
    while (1) {
        FGETS( sLine, fIn );
        iRead = sscanf( sLine, "%s %s %lf %lf %n", saStr[0], saStr[1], &dA, &dB, &iDesc );
        if ( iRead <= 0 )
            break;
        MESSAGE("Read %d H-BOND terms from: %s", iRead, sLine );
        iParmSetAddHBond( psParms, saStr[0], saStr[1], dA, dB, &sLine[iDesc] );
    }
    MESSAGE("Done reading H-BONDs.\n" );
}


/*
 *      zAmberReadParmSetNonBonds
 *
 *      Read non-bond parameters.
 */
static void
zAmberReadParmSetNonBonds( VARARRAY *vaPNonBonds, FILE *fIn )
{
    STRING          sLine;
    int             iRead, iDesc;
    STRING          saStr;
    double          dRStar, dDepth, dRStar14, dDepth14;
    NONBONDt        nbNonBond;
    VARARRAY        vaNonBonds;

    memset(saStr, 0, sizeof(saStr));                        /* for Purify */
    memset(&nbNonBond, 0, sizeof(nbNonBond));               /* for Purify */
    vaNonBonds = vaVarArrayCreate(sizeof(NONBONDt));
    MESSAGE("Reading NON-BONDs.\n" );
    while (1) {
        FGETS( sLine, fIn );
        if( GDefaults.bCharmm ){
            iRead = sscanf( sLine, "%s %lf %lf %lf %lf %n",
                    saStr, &dRStar, &dDepth, &dRStar14, &dDepth14, &iDesc );
            if ( iRead <= 0 )
                break;
            strcpy( nbNonBond.sNonBondType, saStr );
            nbNonBond.dRStar = dRStar;
            nbNonBond.dDepth = dDepth;
            nbNonBond.dRStar14 = dRStar14;
            nbNonBond.dDepth14 = dDepth14;
            strcpy( nbNonBond.sDesc, &sLine[iDesc] );
            VarArrayAdd( vaNonBonds, (GENP)&nbNonBond );
        } else {
            iRead = sscanf( sLine, "%s %lf %lf %n", saStr, &dRStar, &dDepth, &iDesc );
            if ( iRead <= 0 )
                break;
            strcpy( nbNonBond.sNonBondType, saStr );
            nbNonBond.dRStar = dRStar;
            nbNonBond.dDepth = dDepth;
            nbNonBond.dRStar14 = dRStar;
            nbNonBond.dDepth14 = dDepth;
            strcpy( nbNonBond.sDesc, &sLine[iDesc] );
            VarArrayAdd( vaNonBonds, (GENP)&nbNonBond );
        }
        MESSAGE("Read %d NON-BOND terms from: %s", iRead, sLine );
    }
    MESSAGE("Done reading NON-BONDs.\n" );
    *vaPNonBonds = vaNonBonds;
}


/*
 *      zAmberReadParmSetNBPairEdits
 *
 *      Author: David Cerutti (Case group, 2013)
 *
 *      Read non-bonded parameter pair changes specific to a force field
 *      which break the standard Lennard-Jones combining rules.
 */
static void
zAmberReadParmSetNBPairEdits( PARMSET psParms, FILE *fIn, int segfound )
{
    STRING          sLine;
    int             iRead, iDesc;
    STRING          saStr[2];
    double          dEI, dEJ, dRI, dRJ;

    memset(saStr, 0, sizeof(saStr));                        /* for Purify */
    MESSAGE("Reading LJEDITs.\n" );

    /* Scan the file for the card */
    while ( segfound == 0) {
        if (fgets( sLine, MAXSTRINGLENGTH, fIn ) == NULL) {
            break;
        }
        iRead = sscanf( sLine, "%s", saStr[0]);
        if (iRead > 0) {
            segfound = strcmp(saStr[0], "LJEDIT") == 0;
        }
    }

    /* Scan the file for Lennard-Jones pair edits */
    while (segfound) {
        if (fgets( sLine, MAXSTRINGLENGTH, fIn ) == NULL) {
            break;
        }
        iRead = sscanf( sLine, "%s %s %lf %lf %lf %lf %n", saStr[0], saStr[1],
                        &dRI, &dEI, &dRJ, &dEJ, &iDesc);
        if ( iRead <= 0 ) {
            break;
        } else if ( iRead < 6 ) {
            VPWARN("Incomplete LJEDIT line (skipping):\n%s", sLine );
            continue;
        } else if ( iRead > 6 ) {
            /* Silently ignore trailing characters */
        }
        MESSAGE("Read %d LJEDIT terms from: %s", iRead, sLine );
        iParmSetAddNBEdit( psParms, saStr[0], saStr[1], dEI, dEJ, dRI, dRJ, &sLine[iDesc] );
    }
    MESSAGE("Done reading LJEDITs.\n" );
}


/*
 *      zbAmberDetermineParmSetFrcModType
 *
 *      Determine whether this is an old AMBER parameter set
 *      or if it is one of the new FRCMOD parameter sets.
 */
static bool
zbAmberDetermineParmSetFrcModType( FILE *fIn, bool *bPMass, bool *bPNonBond )
{
STRING          sLine;
bool            bNew;

    bNew = false;
    *bPMass = false;
    *bPNonBond = false;
    do {
        FGETS( sLine, fIn );
        if ( strncmp( sLine, "MASS", 4 ) == 0 ) {
            bNew = true;
            *bPMass = true;
        } else if ( strncmp( sLine, "BOND", 4 ) == 0 ) {
            bNew = true;
        } else if ( strncmp( sLine, "ANGL", 4 ) == 0 ) {
            bNew = true;
        } else if ( strncmp( sLine, "DIHE", 4 ) == 0 ) {
            bNew = true;
        } else if ( strncmp( sLine, "IMPR", 4 ) == 0 ) {
            bNew = true;
        } else if ( strncmp( sLine, "HBON", 4 ) == 0 ) {
            bNew = true;
        } else if ( strncmp( sLine, "NONB", 4 ) == 0 ) {
            bNew = true;
            *bPNonBond = true;
        } else if ( strncmp( sLine, "LJEDIT", 4 ) == 0 ) {
            bNew = true;
        } else if ( strncmp( sLine, "NBEDIT", 4 ) == 0 ) {
            bNew = true;
        } else if ( strncmp( sLine, "CMAP", 4 ) == 0 ) {
            bNew = true;
        } else if ( strncmp( sLine, "LJEDIT", 4 ) == 0 ) {
            bNew = true;
        }
    } while ( !feof(fIn) );
    fseek( fIn, 0, 0 );
    return bNew;
}

/*
 *      iAtomTypeToElement
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Convert a type to an element number.
 *      If the type is unknown then return NOELEMENT.
 */
static int
iAtomTypeToElement( char *sType )
{
int             iElement, i, iCount;
STRING          sLowerType, sLowerTable;
TYPEELEMENTt    *tP;

    iCount = iVarArrayElementCount( vaAtomTypes );
    if ( iCount == 0 ) {
        VP1(" (UNKNOWN ATOM TYPE: %s [no types loaded])\n", sType );
        return  NOELEMENT ;
    }

    iElement = NOELEMENT;
    strcpy( sLowerType, sType );
    StringLower( sLowerType );
    tP = PVAI( vaAtomTypes, TYPEELEMENTt, 0 );
    for (i=0; i<iCount; i++, tP++) {
        strcpy( sLowerTable, tP->sType );
        StringLower( sLowerTable );
        if ( strcmp( sLowerTable, sLowerType ) == 0 ) {
            iElement = tP->iElement;
            break;
        }
    }
    if ( iElement == NOELEMENT ) {
        VP1("(UNKNOWN ATOM TYPE: %s)\n", sType );
    }
    return iElement;
}




/*
 *      zpsAmberReadParmSetOld
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Read a PARMSET from an AMBER parameter database file.
 */
static PARMSET
zpsAmberReadParmSetOld( FILE *fIn )
{
STRING          sLine;
STRING          saStr;
VARARRAY        vaMasses, vaEquivs, vaNonBonds;
EQUIVt          eEquiv;
PARMSET         psParms;
NONBONDt        *nbNonBondP;
int             i, j, k, iRead, iSet;

    memset(saStr, 0, sizeof(saStr));                    /* for Purify */
    psParms = (PARMSET)oCreate(PARMSETid);

                /* Read and write the title */

    FGETS( psParms->sTitle, fIn );
    StringTrim(psParms->sTitle);
    VP0("Reading title:\n%s\n", psParms->sTitle );

                /* Read the masses */

    zAmberReadParmSetMasses( &vaMasses, fIn );


                /* Skip hydrophilicity stuff */
    while (1) {
        FGETS( sLine, fIn );
        iRead = sscanf( sLine,
                "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s",
                        saStr, saStr, saStr, saStr, saStr,
                        saStr, saStr, saStr, saStr, saStr,
                        saStr, saStr, saStr, saStr, saStr,
                        saStr, saStr, saStr, saStr, saStr );
        if ( iRead < 20 ) break;
    }

                /* Read the bond parameters */

    MESSAGE("Reading bond parameters\n" );
    zAmberReadParmSetBonds( psParms, fIn );



                /* Read the angle parameters */

    MESSAGE("Reading angle parameters\n" );
    zAmberReadParmSetAngles( psParms, fIn );



                /* Read the torsion parameters */

    MESSAGE("Reading torsion parameters\n" );
    zAmberReadParmSetPropers( psParms, fIn );



                /* Read the improper parameters */

    MESSAGE("Reading improper parameters\n" );
    zAmberReadParmSetImpropers( psParms, fIn );




                /* Read the H-bond parameters */

    MESSAGE("Reading H-Bond parameters\n" );
    zAmberReadParmSetHBonds( psParms, fIn );






                /* Read the equivalences for non-bond parameters */

    MESSAGE("Reading non-bond equivalences\n" );
    vaEquivs = vaVarArrayCreate(sizeof(EQUIVt));
    memset(&eEquiv, 0, sizeof(eEquiv));                 /* for Purify */
    while (1) {
        FGETS( sLine, fIn );
        if ( strlen(sLine)==0 )
                break;
                        /* Remove the \n at the end */
        sLine[strlen(sLine)-1] = '\0';
        sRemoveLeadingSpaces( sLine );
        if ( strlen(sLine)==0 )
                break;
        MESSAGE("Read: %s\n", sLine );
        sRemoveFirstString( sLine, eEquiv.sName );
        for (i=0; i<MAXEQUIV; i++) {
            sRemoveLeadingSpaces(sLine);
            if ( strlen(sLine)==0 )
                break;
            sRemoveFirstString( sLine, eEquiv.saEquivs[i] );
        }
        if ( strlen(sLine)!=0 ){
             VPWARN("Parmset: equiv %s: Unrecognized equivs (> %d)!!\n",
                        eEquiv.sName, MAXEQUIV );
        }
        eEquiv.iEquivs = i;
        VarArrayAdd( vaEquivs, (GENP)&eEquiv );
    }

                /* Read in the Non-Bond parameters, look them up in the */
                /* Masses table and then add them and their equivalences */
                /* to the PARMSET */

                /* Skip a line */
    FGETS( sLine, fIn );
    MESSAGE("Reading non-bond parameters\n" );
    zAmberReadParmSetNonBonds( &vaNonBonds, fIn );

    /*
     *  correlate nb's & masses
     */
    nbNonBondP = PVAI( vaNonBonds, NONBONDt, 0 );
    for ( i=0; i<iVarArrayElementCount(vaNonBonds); i++, nbNonBondP++ ) {
        double          dMass=0.0, dPolar=0.0, dDepth, dRStar, dDepth14, dRStar14;
        double          dScreenF=0.0;
        int             iElement=-10, iHybridization=HUNDEFINED;
        MASSt           *MassP;
        EQUIVt          *EquivP;

        strcpy( saStr, nbNonBondP->sNonBondType );
        dRStar = nbNonBondP->dRStar;
        dDepth = nbNonBondP->dDepth;
        dRStar14 = nbNonBondP->dRStar14;
        dDepth14 = nbNonBondP->dDepth14;

                        /* Lookup the MASS if it is defined */
        /*
         *  look up mass/polar param TODO - use avl instead of vararray
         *      AND CHECK FOR DUPLICATE DEFN'S!
         */
        iSet = 0;
        MassP = PVAI(vaMasses,MASSt,0);
        for ( j=0; j<iVarArrayElementCount(vaMasses); j++, MassP++ ) {
            if ( strcmp( saStr, MassP->sType )==0 ) {
                dMass = MassP->dMass;
                dPolar = MassP->dPolar;
                dScreenF = MassP->dScreenF;
                iElement = MassP->iElement;
                iHybridization = MassP->iHybridization;
                iSet = 1;
            }
        }
        if ( iSet == 0 ) {
            VPWARN("No mass was defined for non-bond atom type: %s - ignoring\n",
                        saStr );
        } else {
                MESSAGE("Adding %d atom type: %s  %lf %lf %lf %lf\n",
                        i, saStr, dMass, dPolar, dRStar, dDepth );
                iParmSetAddAtom( psParms, saStr, dMass, dPolar,
                                dDepth, dRStar, dDepth14, dRStar14,
                                dScreenF,
                                iElement!=-10?iElement
                                              :iAtomTypeToElement(saStr),
                                iHybridization!=HUNDEFINED?iHybridization
                                              :iAtomTypeHybridization(saStr),
                                nbNonBondP->sDesc );
        }
        EquivP = PVAI(vaEquivs,EQUIVt,0);
        for ( j=0; j<iVarArrayElementCount(vaEquivs); j++, EquivP++ ) {
            if ( strcmp( saStr, EquivP->sName )==0 ) {
                for ( k=0; k<EquivP->iEquivs; k++ ) {
                    if ( iSet == 0 ) {
                        VPWARN("No mass equivalenced type: %s - skipping\n",
                                    EquivP->saEquivs[k] );
                        continue;
                    }
                    iParmSetAddAtom( psParms, EquivP->saEquivs[k],
                       dMass, dPolar, dDepth, dRStar, dDepth14, dRStar14,
                       dScreenF,
                       iElement!=-10?iElement
                                    :iAtomTypeToElement(saStr),
                       iHybridization!=HUNDEFINED?iHybridization
                                    :iAtomTypeHybridization(saStr),
                       nbNonBondP->sDesc );
                }
            }
        }
    }

    /* Read in non-bonded pair interaction edits, check their */
    /* validity against existing atom types, and add them to  */
    /* the parmset.                                           */
    MESSAGE("Seeking non-bond pair edits\n" );
    zAmberReadParmSetNBPairEdits( psParms, fIn, 0 );

    VarArrayDestroy( &vaMasses );
    VarArrayDestroy( &vaEquivs );
    VarArrayDestroy( &vaNonBonds );

    return psParms;
}




/*
 *      zpsAmberReadParmSetFrcMod
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Read a PARMSET from an AMBER parameter database file.
 */
static PARMSET
zpsAmberReadParmSetFrcMod( FILE *fIn )
{
    STRING          sLine, saStr;
    VARARRAY        vaMasses, vaNonBonds;
    PARMSET         psParms;
    int             i, j;

    VPTRACEENTER("zpsAmberReadParmSetFrcMod" );
    memset(saStr, 0, sizeof(saStr));                    /* for Purify */
    psParms = (PARMSET)oCreate(PARMSETid);

                /* Read and write the title */

    FGETS( psParms->sTitle, fIn );
    StringTrim(psParms->sTitle);
    VP0("Reading title:\n%s\n", psParms->sTitle );

    vaMasses = NULL;
    vaNonBonds = NULL;

    // Mass and NonBond are deferred; collected int VARARRAYs
    // and then saved in ParmSet as a common Atom parameter set
    do {
        FGETS( sLine, fIn );
        if ( strncmp( sLine, "MASS", 4 ) == 0 ) {
            zAmberReadParmSetMasses( &vaMasses, fIn );
        } else if ( strncmp( sLine, "BOND", 4 ) == 0 ) {
            zAmberReadParmSetBonds( psParms, fIn );
        } else if ( strncmp( sLine, "ANGL", 4 ) == 0 ) {
            zAmberReadParmSetAngles( psParms, fIn );
        } else if ( strncmp( sLine, "DIHE", 4 ) == 0 ) {
            zAmberReadParmSetPropers( psParms, fIn );
        } else if ( strncmp( sLine, "IMPR", 4 ) == 0 ) {
            zAmberReadParmSetImpropers( psParms, fIn );
        } else if ( strncmp( sLine, "HBON", 4 ) == 0 ) {
            zAmberReadParmSetHBonds( psParms, fIn );
        } else if ( strncmp( sLine, "NONB", 4 ) == 0 ) {
            zAmberReadParmSetNonBonds( &vaNonBonds, fIn );
        } else if ( strncmp( sLine, "CMAP", 4 ) == 0 ) {
            zAmberReadParmSetCMAP( psParms, fIn );
        } else if ( strncmp( sLine, "IPOL", 4 ) == 0 ) {
            zAmberReadParmSetIPOL( fIn );
        } else if ( strncmp( sLine, "LJEDIT", 6 ) == 0 ||
                    strncmp( sLine, "NBEDIT", 6 ) == 0 ) {
            zAmberReadParmSetNBPairEdits( psParms, fIn, 1 );
        } else if ( strncmp( sLine, "END", 3 ) == 0 ) {
            continue;
        } else {
            char *s = sLine;
            while ( isspace( *s ) ) ++s ;
            if ( strlen(s) > 0 ) {
                VPWARN("Unknown keyword (%s) in parameter file. Perhaps "
                        "a format issue?\n", sLine );
            }
        }
    } while ( !feof(fIn) );

                /* If there are NonBond parameters then use them */
                /* this implies that there are MASSES */

    if ( vaNonBonds ) {
        NONBONDt        *nbP = PVAI( vaNonBonds, NONBONDt, 0);

        for ( i=0; i<iVarArrayElementCount(vaNonBonds); i++, nbP++ ) {
            double      dMass, dPolar, dRStar, dDepth, dRStar14, dDepth14;
            double      dScreenF;
            int         iElement, iHybridization;
            int         iSet = 0;
            MASSt       *mP;

            strcpy( saStr, nbP->sNonBondType );
            dRStar = nbP->dRStar;
            dDepth = nbP->dDepth;
            dRStar14 = nbP->dRStar14;
            dDepth14 = nbP->dDepth14;

                            /* Lookup the MASS if it is defined */
            mP = PVAI(vaMasses,MASSt,0);
            for ( j=0; j<iVarArrayElementCount(vaMasses); j++, mP++ ) {
                if ( strcmp( saStr, mP->sType )==0 ) {
                        dMass = mP->dMass;
                        dPolar = mP->dPolar;
                        dScreenF = mP->dScreenF;
                        iElement = mP->iElement;
                        iHybridization = mP->iHybridization;
                        iSet = 1;
                        break;
                }
            }
            if ( iSet == 0 ) {
                VPWARN("No mass defined for non-bond atom type: %s - skipping\n",
                            saStr );
                continue;
            }
            MESSAGE("Adding %d atom type: %s  %lf %lf %lf %lf\n",
                            i, saStr, dMass, dPolar, dRStar, dDepth );
            iParmSetAddAtom( psParms, saStr, dMass, dPolar, dDepth, dRStar,
                    dDepth14, dRStar14, dScreenF,
                                iElement!=-10?iElement
                                             :iAtomTypeToElement(saStr),
                                iHybridization!=HUNDEFINED?iHybridization
                                             :iAtomTypeHybridization(saStr),
                                nbP->sDesc );
        }
    }

    if ( vaMasses )
        VarArrayDestroy( &vaMasses );
    if ( vaNonBonds )
        VarArrayDestroy( &vaNonBonds );

    VPTRACEEXIT("zpsAmberReadParmSetFrcMod" );
    return psParms;
}



/*
 *      uAmberReadUnitFromPrep
 *
 *      Author: Christian Schafmeister (1991)
 *
 *TODO: Currently this routine does not do a distance search to
 *TODO: create bonds between atoms a specified distance apart.
 *TODO: The PREP file contains the distance.
 *
 *      Read a UNIT from a PREP file.
 *      If the first line in the file is 'STOP' then return NULL.
 */
typedef struct {
        VECTOR  vPos;
        ATOM    aAtom;
        int     iUnfilled;
        int     iBack;
} TREENODEt;

                /* Define a special type for sorting prep format lines */
typedef struct  {
        int     iLineNumber;
        STRING  sLine;
} INPUTLINEt;

static UNIT
uAmberReadUnitFromPrep( FILE *fIn )
{
TREENODEt       *cPCoor;
TREENODEt       cCoor;
STRING          sLine, sChain, sDesc, sResName;
STRING          saStr[12];
int             iRead, iIndex;
UNIT            uUnit = NULL;
RESIDUE         rRes;
ATOM            aAtom, aAtom2;
VARARRAY        vaCoor = NULL;
VARARRAY        vaLines = NULL;
STRING          sName, sType, sSave;
double          dX, dY, dZ, dCharge;
int             iBond, iAngle, iTorsion;
double          dBond, dAngle, dTorsion;
VECTOR          vPos, vPos1, vPos2, vPos3;
ATOM            aMain0, aMain1;
bool            bFirstTime;
double          da[5];
int             i, iLen, iErrors, iCharge, iCharges = 0, iChgWarn = 0, iAtoms = 0;
double          dCutoff;
INPUTLINEt      ilLine;
int             iaTreeStack[500];
int             iTreeStackTop;
bool            bUseFirstColumn;
#define         TREEPOP()       ( iTreeStackTop-- )
#define         TREEPUSH( i )   ( iaTreeStack[++iTreeStackTop] = i )
#define         iTREETOP()      ( iaTreeStack[iTreeStackTop] )
#define         iTREESTACKPOS() ( iTreeStackTop )
#define         bTREEATTOP()    ( iTreeStackTop == 0 )

                /* Read the long name of the UNIT and null-terminate it */

    FGETS( sDesc, fIn );
    if ( strncmp( sDesc, "STOP", 4 ) == 0 ) {
        return NULL;
    }
    if ( feof( fIn ) ) {
        VPWARN("Unexpected EOF: Assuming STOP\n" );
        return NULL;
    }
    StringTrim(sDesc);

    uUnit = (UNIT)oCreate(UNITid);
    rRes  = (RESIDUE)oCreate(RESIDUEid);
    ContainerAdd( (CONTAINER)uUnit, (OBJEKT)rRes );
    DEREF(rRes);

                /* Skip a line (this is the mainchain filename, if present) */
    FGETS( sLine, fIn );
    if ( feof( fIn ) ) {
        VPERROR("Unexpected EOF: discarding residue (%s)\n",  sDesc );
        goto ERROR;
    }

                /* Define an array to store the coordinates generated */
                /* from the tree-matrix */

    vaCoor = vaVarArrayCreate( sizeof(TREENODEt) );

                /* Read in the short name, INT/XYZ, 0/1 */

    FGETS( sLine, fIn );
    if ( feof( fIn ) ) {
        VPERROR("Unexpected EOF: discarding residue (%s)\n", sDesc );
        goto ERROR;
    }
    sscanf( sLine, "%s %s %s", sResName, saStr[2], saStr[3] );
    ContainerSetName( uUnit, sResName );
    ContainerSetName( rRes, sResName );
    StringTrim(sDesc);
    if (!strcmp(sDesc,"This is a remark line")) sDesc[0]=0;
    if ( strlen(sDesc) > 0 ) {
        UnitSetDescription( uUnit, sDesc );
        VP2("Residue %s Desc:\n%s\n", sResName, sDesc );
    }

                /* Read the line: CORR/CHANGE, OMIT/NOMIT, DU, ALL/BEG */

    FGETS( sLine, fIn );
    if ( feof( fIn ) ) goto EOFERROR;
    sscanf( sLine, "%s", saStr[1] );
    bUseFirstColumn = true;
    if ( strncmp( saStr[1], "CORR", 4 ) == 0 ) {
        bUseFirstColumn = false;
    }

                /* Read the CUT line */

    FGETS( sLine, fIn );
    if ( feof( fIn ) ) goto EOFERROR;
    sscanf( sLine, "%lf", &dCutoff );


                /* Read all of the lines defining the atoms so that */
                /* we can sort then by the first number */

    vaLines = vaVarArrayCreate(sizeof(INPUTLINEt));
    memset(ilLine.sLine, 0, sizeof(ilLine.sLine));              /* for Purify */
    for (;;) {
        FGETS( sLine, fIn );
        if ( feof( fIn ) ) goto EOFERROR;
        iRead = sscanf( sLine, "%s %s %s %s %s %s %s %s %s %s %s",
                        saStr[1], saStr[2], saStr[3], saStr[4], saStr[5],
                        saStr[6], saStr[7], saStr[8], saStr[9], saStr[10],
                        saStr[11] );
        if ( iRead <= 0 )
                break;
        ilLine.iLineNumber = atoi(saStr[1]);
        strcpy( ilLine.sLine, sLine );
        MESSAGE("iRead = %d   Adding: |%s|\n", iRead, sLine );
        VarArrayAdd( vaLines, (GENP)&ilLine );
    }

        /* Save the last line read, because it will be a command */
        /* that has to be parsed after all the atoms are defined */

    strcpy( sSave, sLine );

#ifdef  DEBUG
MESSAGE("------------\n" );
for ( i=0; i<iVarArrayElementCount(vaLines); i++ ) {
    MESSAGE("%s", PVAI(vaLines,INPUTLINEt,i)->sLine );
}
MESSAGE("============\n" );
#endif
    if ( bUseFirstColumn ) {
        INPUTLINEt      *Pil = PVAI(vaLines,INPUTLINEt,0);

        SortByInteger( (GENP) Pil,
                        iVarArrayElementCount(vaLines),
                        sizeof(INPUTLINEt),
                        (GENP) &Pil->iLineNumber, true );
    }
#ifdef  DEBUG
MESSAGE("------------\n" );
for ( i=0; i<iVarArrayElementCount(vaLines); i++ ) {
    MESSAGE("%s", PVAI(vaLines,INPUTLINEt,i)->sLine );
}
MESSAGE("============\n" );
#endif


                /* Read in ATOM lines until the first string is
                 * the keyword LOOP, IMPROPER, DONE, CHARGE
                 * Create atoms for every entry except for those of
                 * type DU
                 * Read the tree type records and build a tree using
                 * the back indices (iBack) in the vaCoor array
                 * and the index stack iaTreeStack/iTreeStackTop
                 * Keep adding back indices to a node until the
                 * iUnfilled field is zero, then pop the stack */

    iTreeStackTop = -1;
    TREEPUSH(-1);
    MESSAGE("Reading atom records\n" );
    aMain0 = NULL;
    aMain1 = NULL;
    iIndex = 0;
    iAngle = -1;
    iTorsion = -1;
    memset(&cCoor, 0, sizeof(cCoor));                   /* for Purify */
    for ( i=0; i<iVarArrayElementCount(vaLines); i++ ) {
        strcpy( sLine, PVAI( vaLines, INPUTLINEt, i )->sLine );
        iRead = sscanf( sLine, "%s %s %s %s %s %s %s %s %s %s %s",
                        saStr[1], saStr[2], saStr[3], saStr[4], saStr[5],
                        saStr[6], saStr[7], saStr[8], saStr[9], saStr[10],
                        saStr[11] );
        if ( iRead <= 0 ) {
            if ( feof(fIn) ) goto DONE;
            continue;
        }
        MESSAGE("Read atom: %s\n", saStr[2] );
        iAtoms++;

                /* Create an entry in the vaCoor array */

        cCoor.iBack = 0;
        cCoor.iUnfilled = 0;
        VarArrayAdd( vaCoor, (GENP)&cCoor );
        cPCoor = PVAI( vaCoor, TREENODEt, iIndex );

        iErrors = 0;
        strcpy( sName, saStr[2] );
        iLen = strlen( sName );
        if ( iLen > 4 ) {
            VPERROR("residue %s atom %s: name must be <= 4 characters\n", sResName, sName );
            iErrors++;
        }

        strcpy( sType, saStr[3] );
        iLen = strlen( sType );
        if ( iLen > 2 ) {
            VPERROR("residue %s atom %s type %s: type must be 1 or 2 characters\n", sResName, sName, sType );
            iErrors++;
        }

        strcpy( sChain, saStr[4] );
        iLen = strlen( sChain );
        if ( iLen > 1 ) {
            VPERROR("residue %s atom %s: chain type [%s] must be 1 character\n", sResName, sName, sChain );
            iErrors++;
        }

        if ( iErrors )
            goto ERROR;


                /* Define the node type */

        switch ( sChain[0] ) {
            case 'M':
                cPCoor->iUnfilled = 10000;
                break;
            case 'B':
                cPCoor->iUnfilled = 2;
                break;
            case '3':
                cPCoor->iUnfilled = 3;
                break;
            case 'S':
                cPCoor->iUnfilled = 1;
                break;
            case 'E':
                cPCoor->iUnfilled = 0;
                break;
            case '4':
                cPCoor->iUnfilled = 4;
                break;
            case '5':
                cPCoor->iUnfilled = 5;
                break;
            case '6':
                cPCoor->iUnfilled = 6;
                break;
            default:
                VPFATAL("Atom %s: Invalid chain specifier [%c] in PREP file.\n",
                        saStr[2], sChain[0] );
                VPFATALDELAYEDEXIT("%s\n", sLine );
                goto ERROR;
        }

                /* Point the new node back to the top node on the stack */
                /* If the top node is not accepting any more */
                /* connections then pop it. */

        cPCoor->iBack = iTREETOP();
        if ( iTREETOP() >= 0 ) {
            PVAI( vaCoor, TREENODEt, iTREETOP() )->iUnfilled--;
            if ( PVAI( vaCoor, TREENODEt, iTREETOP() )->iUnfilled == 0 ) {
                TREEPOP();
                if ( iTREESTACKPOS() <= 0 ) DFATAL("Stack is empty" );
            }
        }

                /* Now if the current node is unfilled then push it on */
                /* the tree stack */

        if ( cPCoor->iUnfilled != 0 ) TREEPUSH( iIndex );

                /* Get the coordinates for the atom */

        iBond = cPCoor->iBack;
        if ( iRead == 7 ) {
            bStringToDouble( saStr[5], &dX );
            bStringToDouble( saStr[6], &dY );
            bStringToDouble( saStr[7], &dZ );
            dCharge = 0.0;
        } else if ( iRead == 8 ) {
            bStringToDouble( saStr[5], &dX );
            bStringToDouble( saStr[6], &dY );
            bStringToDouble( saStr[7], &dZ );
            bStringToDouble( saStr[8], &dCharge );
            iCharges++;
        } else {
            if ( iBond >= 0 ) {
                iAngle= PVAI( vaCoor, TREENODEt, iBond )->iBack;
                if ( iAngle >= 0 ) {
                    iTorsion = PVAI( vaCoor, TREENODEt, iAngle )->iBack;
                }
            }
            bStringToDouble( saStr[8], &dBond );
            bStringToDouble( saStr[9], &dAngle );
            dAngle *= DEGTORAD;
            bStringToDouble( saStr[10], &dTorsion );
            dTorsion *= DEGTORAD;
            if ( iBond == -1 ) {
                ZMatrixNothing( &vPos );
            } else if ( iAngle == -1 ) {
                vPos1 = PVAI( vaCoor, TREENODEt, iBond )->vPos;
                ZMatrixBond( &vPos, &vPos1, dBond );
            } else if ( iTorsion == -1 ) {
                vPos1 = PVAI( vaCoor, TREENODEt, iBond )->vPos;
                vPos2 = PVAI( vaCoor, TREENODEt, iAngle )->vPos;
                ZMatrixBondAngle( &vPos, &vPos1, &vPos2, dBond, dAngle );
            } else {
                vPos1 = PVAI( vaCoor, TREENODEt, iBond )->vPos;
                vPos2 = PVAI( vaCoor, TREENODEt, iAngle )->vPos;
                vPos3 = PVAI( vaCoor, TREENODEt, iTorsion )->vPos;
                ZMatrixBondAngleTorsion( &vPos, &vPos1, &vPos2, &vPos3,
                                dBond, dAngle, dTorsion );
            }
            dX = dVX(&vPos);
            dY = dVY(&vPos);
            dZ = dVZ(&vPos);
            if (iRead < 11)
                dCharge = 0.0;
            else {
                bStringToDouble( saStr[11], &dCharge );
                iCharges++;
            }
/* printf("--- %f [%d]\n", dCharge, iRead);*/
/* %%% */
        }
        if ( iIndex >= iVarArrayElementCount(vaCoor) ) {
            VarArraySetSize( vaCoor, iIndex+1 );
        }
        cPCoor = PVAI( vaCoor, TREENODEt, iIndex );
        VectorDef( &(cPCoor->vPos), dX, dY, dZ );
        aAtom = NULL;
        if ( strcmp( sType, "DU" ) != 0 ) {
            MESSAGE("Creating atom: %s\n", sName );
            aAtom = (ATOM)oCreate(ATOMid);
            ContainerSetName( aAtom, sName );
            AtomSetElement( aAtom, iAtomTypeToElement(sType) );
            AtomSetType( aAtom, sType );
            VectorDef( &vPos, dX, dY, dZ );
            AtomSetPosition( aAtom, vPos );
            AtomSetCharge( aAtom, dCharge );
            ContainerAdd( (CONTAINER)rRes, (OBJEKT)aAtom );
            DEREF(aAtom);
            if ( iBond != -1 ) {
                aAtom2 = (PVAI( vaCoor, TREENODEt, iBond))->aAtom;
                if ( aAtom2 != NULL ) {
                    MESSAGE("Created bond between: %s - %s\n",
                                sContainerName(aAtom),
                                sContainerName(aAtom2) );
                    AtomBondTo( aAtom, aAtom2 );
                }
            }
            if ( (strcmp( sChain, "M" )==0) && (strcmp( sType, "DU" )!=0) ) {
                if ( aMain0 == NULL ) aMain0 = aAtom;
                aMain1 = aAtom;
            }
        }
        else {
            if ( strncmp( sName, "DUMM", 4 ) != 0 ) {
                VPNOTE("Entry of type DUmmy; atom %s not created\n", sName );
            }
        }
        cPCoor->aAtom = aAtom;
        iIndex++;
    }

    if ( iCharges == 0 ) {
        VPNOTE("(no charges read on atoms lines in %s)\n", sResName );
    } else {
        if ( iCharges != iAtoms ){
            VPWARN("%d of %d atoms missing charges on atoms lines: %s\n",
                    iAtoms-iCharges, iAtoms, sResName );
        }
        iChgWarn = 1;
        iCharges = 0;
    }

        /* Get the last line read which is the first command */

    strcpy( sLine, sSave );

    MESSAGE("Parsing commands after atom records\n" );

                /* Parse the commands that follow the atom records */
    bFirstTime = true;
    while (1) {
        if ( bFirstTime ) bFirstTime = false;
        else            FGETS( sLine, fIn );
        if ( feof( fIn ) ) {
                VPERROR("Unexpected EOF: discarding residue (%s)\n", sResName );
                return NULL;
        }
        iRead = sscanf( sLine, "%s", saStr[1] );
        if ( iRead <= 0 ) continue;
        MESSAGE("Parsed command: %s\n", saStr[1] );
        if ( strcmp( saStr[1], "DONE" ) == 0 ) break;
        else if ( strcmp( saStr[1], "LOOP" ) == 0 ) {
            int problem = 0;
            while (1) {
                FGETS( sLine, fIn );
                if ( feof( fIn ) ) {
                        VPERROR("Unexpected EOF: discarding (%s)\n", sResName );
                        goto ERROR;
                }
                iRead = sscanf( sLine, "%s %s", saStr[1], saStr[2] );
                if ( iRead <= 0 ) break;
                if ( iRead != 2 ) {
                    VPERROR("** only got 1 atom (%s) from LOOP line - skipping\n",
                                        sLine );
                    problem++;
                    continue;
                }
                aAtom = (ATOM)cContainerFindName( (CONTAINER)rRes,
                                                        ATOMid, saStr[1] );
                aAtom2= (ATOM)cContainerFindName( (CONTAINER)rRes,
                                                        ATOMid, saStr[2] );
                if ( aAtom == NULL ) {
                    VPERROR("** LOOP atom %s not found - bond not formed\n",
                            saStr[1] );
                    problem++;
                    continue;
                }
                if ( aAtom2 == NULL ) {
                    VPERROR("** LOOP atom %s not found - bond not formed\n",
                            saStr[2] );
                    problem++;
                    continue;
                }
                if (bAtomBondedTo(aAtom, aAtom2)) {
                    VPWARN("%s:  LOOP: redundant bond  %s--%s  ignored\n", sResName, saStr[1], saStr[2] );
                } else
                        AtomBondTo( aAtom, aAtom2 );
            }
            if (problem) {
                VPERROR("Discarding residue (%s) to EOF\n", sResName );
                goto ERROR;
            }
        } else if ( strcmp( saStr[1], "CHARGE" )==0 ) {
            if ( iChgWarn ){
                VPWARN("Per-line charges being overridden by CHARGE block in %s\n", sResName );
            }
            iCharge = 3;
            while (1) {
                FGETS( sLine, fIn );
                if ( feof( fIn ) ) {
                        VPERROR("Unexpected EOF: discarding residue (%s)\n", sResName );
                        goto ERROR;
                }
                iRead = sscanf( sLine, "%lf %lf %lf %lf %lf",
                                &da[0], &da[1], &da[2],
                                &da[3], &da[4] );
                if ( iRead<=0 ) break;
                for ( i=0; i<iRead; i++ ) {
                    if ( iCharge >= iIndex ) {
                        VPWARN("Skipping extra charge in: %s\n", sResName );
                        continue;
                    }
                    AtomSetCharge( PVAI( vaCoor, TREENODEt, iCharge )->aAtom,
                                        da[i] );
                    iCharge++;
                }
            }
            if ( iCharge < iAtoms ) {
                    VPERROR("Incomplete CHARGE block in: %s: %d/%d\n", sResName, iCharge, iAtoms );
                    goto ERROR;
            }
        } else if ( strcmp( saStr[1], "IMPROPER" )==0 ) {
            int problem = 0;
            IMPROPERt   Imp;

                /* maybe will use this info someday; for now this is
                   the 1st step towards saving it
                        (debate continues over whether it has value)  */

            rRes->vaImpropers = vaVarArrayCreate( sizeof(IMPROPERt) );
            while (1) {
                FGETS( sLine, fIn );
                if ( feof( fIn ) ) {
                        VPERROR("Unexpected EOF: discarding residue (%s)\n", sResName );
                        goto ERROR;
                }
                iRead = sscanf( sLine, "%s %s %s %s",
                                saStr[1], saStr[2], saStr[3], saStr[4] );
                if ( iRead <= 0 ) break;
                if ( iRead != 4 ) {
                    VPNOTE("** got %d atom%s from IMPROPER line (%s) %s\n",
                                iRead, (iRead>1 ? "s":""), sLine,
                                " - expected 4 (skipping line)" );
                    problem++;
                    continue;
                }
                for ( i=1; i<=4; i++ ) {
                    if ( strcmp(saStr[i], "-M") && strcmp(saStr[i], "+M") &&
                             !cContainerFindName( (CONTAINER)rRes, ATOMid, saStr[i] ) ) {
                        VPNOTE("** %s: 'IMPROPER' atom %s not found\n",  sResName, saStr[i] );
                        problem++;
                    }
                }
                if ( !problem ) {
                    memset(&Imp, ' ', sizeof(Imp));
                    memcpy(Imp.sName1, saStr[1], strlen(saStr[1]));
                    memcpy(Imp.sName2, saStr[2], strlen(saStr[2]));
                    memcpy(Imp.sName3, saStr[3], strlen(saStr[3]));
                    memcpy(Imp.sName4, saStr[4], strlen(saStr[4]));
                    VarArrayAdd( rRes->vaImpropers, (GENP)&Imp );
                }
            }
            if (problem && 0) {
                VPERROR("Discarding residue (%s) to EOF\n", sResName );
                goto ERROR;
            }
        }
    }

    UnitSetHead( uUnit, aMain0 );
    ResidueSetConnectAtom( rRes, NEND, aMain0 );
    UnitSetTail( uUnit, aMain1 );
    ResidueSetConnectAtom( rRes, CEND, aMain1 );

                /* Check if there was a cutoff value, if there was */
                /* use it in a Distance Search to create bonds */

    if ( bUseFirstColumn  &&  dCutoff > 0.01 ) {
        VP1("Distance search to create bonds for: %s  distance: %10.5lf\n",
                sContainerName(uUnit), dCutoff );
        iToolDistanceSearch( (CONTAINER)uUnit, dCutoff, true, /* Absolute distance */
                                DISTANCE_SEARCH_CREATE_BONDS );
    }

DONE:
    VarArrayDestroy( &vaCoor );
    VarArrayDestroy( &vaLines );

    return uUnit;

EOFERROR:
    VPERROR("Unexpected EOF: discarding residue (%s)\n", sResName );
ERROR:
    VarArrayDestroy( &vaCoor );
    VarArrayDestroy( &vaLines );

    Destroy( (OBJEKT *)&uUnit );

    return NULL;
}


/*
 *-------------------------------
 *       Public routines
 */


/*
 *      AmberAddAtomTypes
 *
 *      Author: Bill Ross (1996)
 *
 *      Add to array of atom types / hybridizations
 */
void
AmberAddAtomTypes( LIST lEntries )
{
LISTLOOP        llEntries, llEntry;
ASSOC           aEntry1, aEntry2;
LIST            lEntry;
OBJEKT          oEntry2;
TYPEELEMENTt    tTypeInfo;

    if ( iObjectType(lEntries) != LISTid ) {
        DFATAL("AmberAddAtomTypes: Need LIST" );
    }
    if ( vaAtomTypes == NULL )
        vaAtomTypes = vaVarArrayCreate( sizeof( TYPEELEMENTt ) );

    llEntries = llListLoop(lEntries);
    while ( (aEntry1 = (ASSOC)oListNext(&llEntries)) ) {
        /*
         *  prepare to go over the Type, Element, Hybridization
         *      items - list items w/in the overall list
         */
        lEntry = (LIST)oAssocObject(aEntry1);   /* view current entry as list */
        llEntry = llListLoop(lEntry);           /* set up to loop over it */
        aEntry2 = (ASSOC)oListNext(&llEntry);   /* get 1st thing */

        if ( aEntry2 == NULL ) {
                VPWARN("addAtomTypes: null atom type entry\n" );
                continue;
        }
        oEntry2 = (OBJEKT)oAssocObject(aEntry2); /* get the objekt it is */
        if ( iObjectType(oEntry2) != OSTRINGid ) {
                VPWARN("addAtomTypes: atom type - expected string\n" );
                continue;
        }
#if 0
        /*  Joung-Cheatham use ion types with three characters,
            so skip this warning  */
        if ( strlen( sOString( oEntry2 ) ) > 2 ) {
                VPWARN("addAtomTypes: type %s - max length is 2\n",
                                        sOString( oEntry2 ) );
                continue;
        }
#endif
        memset(&tTypeInfo, 0, sizeof(tTypeInfo));       /* for Purify */
        strcpy( tTypeInfo.sType, sOString( oEntry2 ) );

        aEntry2 = (ASSOC)oListNext(&llEntry);   /* get next thing */
        if ( aEntry2 == NULL ) {
                VPWARN("addAtomTypes: incomplete type entry (%s)\n",
                                        tTypeInfo.sType );
                continue;
        }
        oEntry2 = (OBJEKT)oAssocObject(aEntry2); /* get the objekt it is */
        if ( iObjectType(oEntry2) != OSTRINGid ) {
                VPWARN("addAtomTypes: %s - expected string for element\n",
                                                tTypeInfo.sType );
                continue;
        }
        if ( strlen( sOString( oEntry2 ) ) == 0 ) {
                /* NULL ("") entry for dummy atoms */
                tTypeInfo.iElement = NOELEMENT;
        } else {
                tTypeInfo.iElement = iElementNumber( sOString( oEntry2 ) );
                if ( tTypeInfo.iElement == NOELEMENT ) {
                        VPWARN("atom type %s - unknown element %s\n",
                                                tTypeInfo.sType,
                                                sOString( oEntry2 ) );
                        continue;
                }
        }
        aEntry2 = (ASSOC)oListNext(&llEntry);   /* get next thing */
        if ( aEntry2 == NULL ) {
                VPWARN("incomplete atomtype entry (%s)\n", tTypeInfo.sType );
                continue;
        }
        oEntry2 = (OBJEKT)oAssocObject(aEntry2); /* get the objekt it is */
        if ( iObjectType(oEntry2) == OSTRINGid ) {
                StringLower( sOString( oEntry2 ) );
                if ( !strcmp( sOString( oEntry2 ), "sp2") )
                        tTypeInfo.iHybridization = HSP2;
                else if ( !strcmp( sOString( oEntry2 ), "sp3") )
                        tTypeInfo.iHybridization = HSP3;
                else {
                        VPWARN("atom type %s - unknown hybridization %s\n",
                                                tTypeInfo.sType, sOString( oEntry2 ) );
                        continue;
                }
        } else if ( iObjectType(oEntry2) == OINTEGERid ) {
                tTypeInfo.iHybridization = iOInteger( oEntry2 );
                if (tTypeInfo.iHybridization != HSP2  &&
                    tTypeInfo.iHybridization != HSP3  &&
                    tTypeInfo.iHybridization != 0) {
                        VPWARN("atom type %s - unknown hybridization %d\n",
                                                tTypeInfo.iHybridization );
                        continue;
                }
        } else {
          VPWARN("atom type %s - expected string or integer for hybridization\n",
                                                tTypeInfo.sType );
          continue;
        }
        VarArrayAdd( vaAtomTypes, (GENP)&tTypeInfo );
    }

}

/*
 *      dAmberReadPrepFile
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Read an AMBER PREP file into a dictionary.
 */
DICTIONARY
dAmberReadPrepFile( char *sFilename )
{
DICTIONARY      dUnits;
UNIT            uUnit;
FILE            *fIn;
STRING          sTemp;

    fIn = FOPENCOMPLAIN( sFilename, "r" );
    if ( fIn == NULL )
        return NULL;
    VP0("Loading Prep file: %s\n", GsBasicsFullName );

    dUnits = dDictionaryCreate();

                /* Skip the two header lines */
    FGETS( sTemp, fIn );
    FGETS( sTemp, fIn );

                /* Read each UNIT in the PREP file */

    while ( !feof(fIn) && (uUnit = uAmberReadUnitFromPrep(fIn)) != NULL ) {
        DictionaryAdd( dUnits, sContainerName(uUnit), (GENP)uUnit );
    }

    fclose( fIn );
    return dUnits;
}


// Standard POSIX cksum CRC32 lookup table
static const uint32_t crc32_table[256] = {
    0x00000000, 0x04C11DB7, 0x09823B6E, 0x0D4326D9, 0x130476DC, 0x17C56B6B, 0x1A864DB2, 0x1E475005,
    0x2608EDB8, 0x22C9F00F, 0x2F8AD6D6, 0x2B4BCB61, 0x350C9B64, 0x31CD86D3, 0x3C8EA00A, 0x384FBDBD,
    0x4C11DB70, 0x48D0C6C7, 0x4593E01E, 0x4152FDA9, 0x5F15ABAC, 0x5BD4B61B, 0x569790C2, 0x52568D75,
    0x6A1936C8, 0x6ED82B7F, 0x639B0DA6, 0x675A1011, 0x791D4014, 0x7DDC5DA3, 0x709F7B7A, 0x745E66CD,
    0x9823B6E0, 0x9CE2AB57, 0x91A18D8E, 0x95609039, 0x8B27C03C, 0x8FE6DD8B, 0x82A5FB52, 0x8664E6E5,
    0xBE2B5B58, 0xBAEA46EF, 0xB7A96036, 0xB3687D81, 0xAD2F2D84, 0xA9EE3033, 0xA4AD16EA, 0xA06C0B5D,
    0xD4326D90, 0xD0F37027, 0xDDB056FE, 0xD9714B49, 0xC7361B4C, 0xC3F706FB, 0xCEB42022, 0xCA753D95,
    0xF23A8028, 0xF6FB9D9F, 0xFBB8BB46, 0xFF79A6F1, 0xE13EF6F4, 0xE5FFEB43, 0xE8BCCD9A, 0xEC7DD02D,
    0x34867077, 0x30476DC0, 0x3D044B19, 0x39C556AE, 0x278206AB, 0x23431B1C, 0x2E003DC5, 0x2AC12072,
    0x128E9DCF, 0x164F8078, 0x1B0CA6A1, 0x1FCDBB16, 0x018AEB13, 0x054BF6A4, 0x0808D07D, 0x0CC9CDCA,
    0x7897AB07, 0x7C56B6B0, 0x71159069, 0x75D48DDE, 0x6B93DDDB, 0x6F52C06C, 0x6211E6B5, 0x66D0FB02,
    0x5E9F46BF, 0x5A5E5B08, 0x571D7DD1, 0x53DC6066, 0x4D9B3063, 0x495A2DD4, 0x44190B0D, 0x40D816BA,
    0xACA5C697, 0xA864DB20, 0xA527FDF9, 0xA1E6E04E, 0xBFA1B04B, 0xBB60ADFC, 0xB6238B25, 0xB2E29692,
    0x8AAD2B2F, 0x8E6C3698, 0x832F1041, 0x87EE0DF6, 0x99A95DF3, 0x9D684044, 0x902B669D, 0x94EA7B2A,
    0xE0B41DE7, 0xE4750050, 0xE9362689, 0xEDF73B3E, 0xF3B06B3B, 0xF771768C, 0xFA325055, 0xFEF34DE2,
    0xC6BCF05F, 0xC27DEDE8, 0xCF3ECB31, 0xCBFFD686, 0xD5B88683, 0xD1799B34, 0xDC3ABDED, 0xD8FBA05A,
    0x690CE0EE, 0x6DCDFD59, 0x608EDB80, 0x644FC637, 0x7A089632, 0x7EC98B85, 0x738AAD5C, 0x774BB0EB,
    0x4F040D56, 0x4BC510E1, 0x46863638, 0x42472B8F, 0x5C007B8A, 0x58C1663D, 0x558240E4, 0x51435D53,
    0x251D3B9E, 0x21DC2629, 0x2C9F00F0, 0x285E1D47, 0x36194D42, 0x32D850F5, 0x3F9B762C, 0x3B5A6B9B,
    0x0315D626, 0x07D4CB91, 0x0A97ED48, 0x0E56F0FF, 0x1011A0FA, 0x14D0BD4D, 0x1993DB94, 0x1D52C623,
    0xF12F560E, 0xF5EE4BB9, 0xF8AD6D60, 0xFC6C70D7, 0xE22B20D2, 0xE6EA3D65, 0xEBA91BBC, 0xEF68060B,
    0xD727BBB6, 0xD3E6A601, 0xDEA580D8, 0xDA649D6F, 0xC423CD6A, 0xC0E2D0DD, 0xCDA1F604, 0xC960EBB3,
    0xBD3E8D7E, 0xB9FF90C9, 0xB4BCB610, 0xB07DD0A7, 0xAE3AFBA2, 0xAAFBE615, 0xA7B8C0CC, 0xA379DD7B,
    0x9B3660C6, 0x9FF77D71, 0x92B45BA8, 0x9675461F, 0x8832161A, 0x8CF30BAD, 0x81B02D74, 0x857130C3,
    0x5D8A7079, 0x594B6DCE, 0x54084B17, 0x50C956A0, 0x4E8E06A5, 0x4A4F1B12, 0x470C3DCB, 0x43CD207C,
    0x7B829D11, 0x7F4380A6, 0x7200A67F, 0x76C1BB38, 0x6886EBB7, 0x6C47F600, 0x6104D0D2, 0x65C5CD65,
    0x118BBB09, 0x154AA6BE, 0x18098067, 0x1CC89DD0, 0x028FCDDD, 0x064EFE6A, 0x0B0BEDB3, 0x0FC2F004,
    0x378356B1, 0x33424B06, 0x3E016DD0, 0x3AC07067, 0x2487206D, 0x20463DCB, 0x2D051B03, 0x29C406B4,
    0xC5A94B99, 0xC168562E, 0xCC2B70F7, 0xC8EA6D40, 0xD6AD3D45, 0xD26C20F2, 0xDF2F062B, 0xDBEE1B9C,
    0xE3A1A621, 0xE760BB96, 0xEA23ED4F, 0xEE62F0F8, 0xF0A5D0FD, 0xF464CD4A, 0xF927E593, 0xFDE6F824,
    0x89B89D09, 0x8D7980BE, 0x803A7667, 0x84FB6BDA, 0x9ABCEBD5, 0x9E7DF662, 0x933ED0BB, 0x97FFCD0C,
    0xAFB070B1, 0xAB716D06, 0xA6324BDF, 0xA2F35668, 0xBCB4066D, 0xB8751BDF, 0xB5363D03, 0xB1F720B4
};

#define CHUNK_SIZE 8192
void generate_cksum_fingerprint(FILE *file, char *dest_str, size_t dest_len) {
    uint32_t crc = 0;
    uint64_t total_bytes = 0;
    uint8_t buffer[CHUNK_SIZE];
    size_t bytes_read;

    // 1. Process file chunks
    while ((bytes_read = fread(buffer, 1, CHUNK_SIZE, file)) > 0) {
        total_bytes += bytes_read;
        for (size_t i = 0; i < bytes_read; i++) {
            crc = (crc << 8) ^ crc32_table[((crc >> 24) ^ buffer[i]) & 0xFF];
        }
    }

    // 2. POSIX specification requires appending the file length to the CRC math
    uint64_t len_tmp = total_bytes;
    while (len_tmp > 0) {
        crc = (crc << 8) ^ crc32_table[((crc >> 24) ^ (len_tmp & 0xFF)) & 0xFF];
        len_tmp >>= 8;
    }

    // 3. Finalize CRC (bitwise NOT)
    crc = ~crc;

    // 4. Concatenate <sum>-<size> as a decimal string
    snprintf(dest_str, dest_len, "%u-%llu", crc, (unsigned long long)total_bytes);
    fseek( file, 0, 0 );
}




/*
 *      psAmberReadParmSet
 *
 *      Read either an old AMBER parm set or a new frcmod parmset.
 *      Determine which type.
 */
PARMSET
psAmberReadParmSet( FILE *fIn, char *sFilename )
{
bool            bFrcMod;
bool            bMass, bNonBonds;
PARMSET         psParms;

    bFrcMod = zbAmberDetermineParmSetFrcModType( fIn, &bMass, &bNonBonds );
    if ( !bFrcMod ) {
        psParms = zpsAmberReadParmSetOld( fIn );
    } else {
        /*
         *  frcmod
         */
        if ( bMass != bNonBonds ) {
                VPERROR("Modified force field files must have both MASS "
                    "and NONBOND entries or neither.\n"
                    "Could not load parameter set from %s.\n",
                    sFilename );
                fclose(fIn);
                return NULL;
        }
        VP0("Reading force field modification type file (frcmod)\n" );
        psParms = zpsAmberReadParmSetFrcMod( fIn );
    }

    fseek(fIn,0,0);
    generate_cksum_fingerprint(fIn, psParms->cksum, sizeof(psParms->cksum));
    fclose(fIn);
    strcpy( sParmName(psParms), sFilename );
    return psParms;
}





/*
 *      iAtomTypeHybridization
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the element type of an AMBER atom type.
 *      Return HUNKNOWN if it is undefined.
 */
int
iAtomTypeHybridization( char *sType )
{
int             i, iCount;
TYPEELEMENTt    *tP;

    iCount = iVarArrayElementCount( vaAtomTypes );
    if ( iCount == 0 )
        return  HUNKNOWN ;

    tP = PVAI( vaAtomTypes, TYPEELEMENTt, 0 );
    for (i=0; i<iCount; i++, tP++) {
        if ( strcmp( sType, tP->sType ) == 0 )
            return tP->iHybridization;
    }
    return HUNKNOWN;
}

