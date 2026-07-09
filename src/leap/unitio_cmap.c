
#include        "basics.h"
#include        "classes.h"
#include        "avl.h"
#include        "defaults.h"
#include        "fortran.h"
#include        "mathop.h"
#include        "sort.h"
#include        "cmap.h"
#include        "unitio.h"

#ifdef BINTRAJ
#  include      "netcdf.h"

/* ----------------------------------------------------------------
   NC_CHECK: error handler for NetCDF API calls.
   ---------------------------------------------------------------- */
#define NC_CHECK(call) \
    do { int _e = (call); if (_e != NC_NOERR) { \
        fprintf(stderr, "NetCDF error %s: %s (%s:%d)\n", #call, nc_strerror(_e), __FILE__,__LINE__); \
        FREE(mapflag); FREE(mapidx); FREE(phipsi); FREE(stdpt0); \
        return _e; \
    }} while(0)

#endif

#define COMMENT_DIM(name, n) \
    sprintf(sComment, "%%COMMENT DIMENSION=%s; SIZE=%d", (name), (n)); \
    FortranWriteString(sComment)

#define COMMENT_SIZE(n) \
    sprintf(sComment, "%%COMMENT SIZE=%d", (n)); \
    FortranWriteString(sComment)

#define AMBERINDEX(i)   3*(i-1)

typedef struct {
    SAVETORSIONt *tp;
} SAVETORSIONtp;

CMAP *Gcmap=NULL;
CMAPLST *Gcmaplst=NULL;
int GiCmapNum=0;

static void get4atoms(UNIT u, SAVETORSIONt * pt, SAVEATOMt * sa4[])
{
    sa4[0] = PVAI(u->vaAtoms, SAVEATOMt, pt->iAtom1 - 1);
    sa4[1] = PVAI(u->vaAtoms, SAVEATOMt, pt->iAtom2 - 1);
    sa4[2] = PVAI(u->vaAtoms, SAVEATOMt, pt->iAtom3 - 1);
    sa4[3] = PVAI(u->vaAtoms, SAVEATOMt, pt->iAtom4 - 1);
}

static int copyatoms(int atoms[], SAVETORSIONt * sa4, SAVETORSIONt * sb4,
                     int k1, int k2)
{
    int atma[4], atmb[4], i;

    if (k1 > 0) {
        atma[0] = AMBERINDEX(sa4->iAtom1);
        atma[1] = AMBERINDEX(sa4->iAtom2);
        atma[2] = AMBERINDEX(sa4->iAtom3);
        atma[3] = AMBERINDEX(sa4->iAtom4);
    } else {
        atma[3] = AMBERINDEX(sa4->iAtom1);
        atma[2] = AMBERINDEX(sa4->iAtom2);
        atma[1] = AMBERINDEX(sa4->iAtom3);
        atma[0] = AMBERINDEX(sa4->iAtom4);
    }

    if (k2 > 0) {
        atmb[0] = AMBERINDEX(sb4->iAtom1);
        atmb[1] = AMBERINDEX(sb4->iAtom2);
        atmb[2] = AMBERINDEX(sb4->iAtom3);
        atmb[3] = AMBERINDEX(sb4->iAtom4);
    } else {
        atmb[3] = AMBERINDEX(sb4->iAtom1);
        atmb[2] = AMBERINDEX(sb4->iAtom2);
        atmb[1] = AMBERINDEX(sb4->iAtom3);
        atmb[0] = AMBERINDEX(sb4->iAtom4);
    }
// check overlap
    for (i = 0; i < 3; i++)
        if (atmb[i] != atma[i + 1]) {
            VP0("Atom indices mismatch?? in copyatoms %i %i\n", atma[i + 1],
                 atmb[i]);
            return -1;
        } else {
            atoms[i] = atma[i];
            atoms[i + 2] = atmb[i + 1];
        }
    return 4;

}

static int copyatoms_netcdf(int atoms[], SAVETORSIONt * sa4, SAVETORSIONt * sb4,
                     int k1, int k2)
{
    int atma[4], atmb[4], i;

    if (k1 > 0) {
        atma[0] = sa4->iAtom1;
        atma[1] = sa4->iAtom2;
        atma[2] = sa4->iAtom3;
        atma[3] = sa4->iAtom4;
    } else {
        atma[3] = sa4->iAtom1;
        atma[2] = sa4->iAtom2;
        atma[1] = sa4->iAtom3;
        atma[0] = sa4->iAtom4;
    }

    if (k2 > 0) {
        atmb[0] = sb4->iAtom1;
        atmb[1] = sb4->iAtom2;
        atmb[2] = sb4->iAtom3;
        atmb[3] = sb4->iAtom4;
    } else {
        atmb[3] = sb4->iAtom1;
        atmb[2] = sb4->iAtom2;
        atmb[1] = sb4->iAtom3;
        atmb[0] = sb4->iAtom4;
    }
// check overlap
    for (i = 0; i < 3; i++)
        if (atmb[i] != atma[i + 1]) {
            VP0("Atom indices mismatch?? in copyatoms %i %i\n", atma[i + 1],
                 atmb[i]);
            return -1;
        } else {
            atoms[i] = atma[i];
            atoms[i + 2] = atmb[i + 1];
        }
    return 4;

}

static int cmpresname1(UNIT u, SAVEATOMt * sa4, WRD reslist[], int nres)
{
    int i;
    char *sname;

    sname = PVAI(u->vaResidues, SAVERESIDUEt, sa4->iResidueIndex - 1)->sName;
    for (i = 0; i < nres; i++) {
        if (strcmp(sname, reslist[i]) == 0)
            return (i + 1);
    }

    return -1;

}

static int cmpresname4(UNIT u, SAVEATOMt * sa4[], WRD reslist[], int nres)
{
    int i, j;

    for (j = 0; j < 4; j++) {
        if ((i = cmpresname1(u, sa4[j], reslist, nres)) > 0)
            return (i + 1);
    }

    return -1;

}

static int cmp4vs4(SAVEATOMt * sa4[], WRD atm4[])
{
    int i, l1, l2;
    l1 = 0;
    l2 = 0;
    for (i = 0; i < 4; i++)
        if (strcmp(sa4[i]->sName, atm4[i]) == 0)
            l1++;
    if (l1 == 4)
        return 4;
    for (i = 0; i < 4; i++)
        if (strcmp(sa4[3 - i]->sName, atm4[i]) == 0)
            l2++;
    if (l2 == 4)
        return -4;

    return 0;

}

static int cmp_residx(SAVEATOMt * sa4[], SAVEATOMt * sb4[], int *residx)
{
    int i, l1, l2;
    int idx0, idx1;

    idx0 = sa4[0]->iResidueIndex;
    idx1 = residx[0];
    l1 = 0;
    l2 = 0;
    for (i = 0; i < 4; i++) {
        if ((sa4[i]->iResidueIndex - idx0) == (residx[i] - idx1))
            l1++;
    }
    for (i = 0; i < 4; i++) {
        if ((sb4[i]->iResidueIndex - idx0) == (residx[i + 1] - idx1))
            l2++;
    }

    if (l1 == 4 && l2 == 4)
        return 1;

    return 0;

}

void SaveAmberParmCMAP(UNIT uUnit, FILE * fOut)
{
    //
    // CMAP parameters, Mengjuei Hsieh and Yong Duan
    //
    int i, j, k, l;
    int mapid, maptypes;
//    int mapcount;
    int *mapflag, *mapidx;
    int iNumDIH;
    int nprospect, ires;
    SAVEATOMt *sa4[4], *sb4[4];
    SAVETORSIONt *stPTorsion2, *stPTorsion;
    SAVETORSIONtp *stdpt0;
    STRING sTmp;
    PHIPSI *phipsi;
    CMAPLST *cmaplstt;
    int k1, k2;
    CMAP *Gcmap;
    int maxmap;
    STRING sComment;

    if (!GDefaults.iCMAP)
        return;
    if (GiCmapNum <= 0)
        return;

    for (cmaplstt = Gcmaplst; cmaplstt->next != NULL; cmaplstt = cmaplstt->next) {
        Gcmap = cmaplstt->cmap;
    }
    iNumDIH = iVarArrayElementCount(uUnit->vaTorsions);
    if (iNumDIH > 0) {
        stPTorsion = PVAI(uUnit->vaTorsions, SAVETORSIONt, 0);
    } else {
        return;
    }

//    mapcount = 0;
    MALLOC(mapflag, int *, sizeof(int) * (GiCmapNum + 1));
    MALLOC(mapidx, int *, sizeof(int) * (GiCmapNum + 1));

    i = 0;
    for (i = 0; i <= GiCmapNum; i++) {
        mapflag[i] = 0;
        mapidx[i] = 0;
    }

    maxmap = iNumDIH;
    MALLOC(phipsi, PHIPSI *, sizeof(PHIPSI) * (maxmap));
    MALLOC(stdpt0, SAVETORSIONtp *, sizeof(SAVETORSIONtp) * (iNumDIH));
    nprospect = 0;
//    mapcount = 0;
    // Loop over dihedral list
// pre-filter removes the irrelevant torsions first ...
    for (i = 0; i < iNumDIH; i++, stPTorsion++) {

//            if (stPTorsion->bCalc14 == 0) continue; //cycle

        get4atoms(uUnit, stPTorsion, sa4);

        for (cmaplstt = Gcmaplst; cmaplstt->next != NULL;
             cmaplstt = cmaplstt->next) {
            Gcmap = cmaplstt->cmap;
            if (cmpresname4(uUnit, sa4, Gcmap->reslist, Gcmap->nres) > 0) // potential match
            {
                if (abs(cmp4vs4(sa4, (WRD *) (&Gcmap->atmname[0]))) == 4 || abs(cmp4vs4(sa4, (WRD *) (&Gcmap->atmname[1]))) == 4) {       // keep in the list
                    stdpt0[nprospect].tp = stPTorsion;
                    nprospect++;
                    break;
                }
            } else if (Gcmap->termmap > 0) {
                if (cmpresname4(uUnit, sa4, Gcmap->creslist, Gcmap->nres) > 0) {
                    if (abs(cmp4vs4(sa4, (WRD *) (&Gcmap->catmname[0]))) == 4
                        || abs(cmp4vs4(sa4, (WRD *) (&Gcmap->catmname[1]))) ==
                        4) {
                        stdpt0[nprospect].tp = stPTorsion;
                        nprospect++;
                        break;
                    }
                } else if (cmpresname4(uUnit, sa4, Gcmap->nreslist, Gcmap->nres) >
                           0)
                    if (abs(cmp4vs4(sa4, (WRD *) (&Gcmap->natmname[0]))) == 4
                        || abs(cmp4vs4(sa4, (WRD *) (&Gcmap->natmname[1]))) ==
                        4) {
                        stdpt0[nprospect].tp = stPTorsion;
                        nprospect++;
                        break;
                    }
            }
        }
    }
//
    k = 0;
    int k1c, k2c, k1n, k2n;
    for (i = 0; i < nprospect; i++) {
        stPTorsion = stdpt0[i].tp;
        get4atoms(uUnit, stPTorsion, sa4);
        mapid = 0;
        for (cmaplstt = Gcmaplst; cmaplstt->next != NULL;
             cmaplstt = cmaplstt->next) {
            Gcmap = cmaplstt->cmap;
            mapid++;
            k1 = cmp4vs4(sa4, (WRD *) (&Gcmap->atmname[0]));     // check atmnames 0..3
            k1c = cmp4vs4(sa4, (WRD *) (&Gcmap->catmname[0]));   // check atmnames 0..3
            k1n = cmp4vs4(sa4, (WRD *) (&Gcmap->natmname[0]));   // check atmnames 0..3

            if (Gcmap->termmap == 0) {
                k1c = 0;
                k1n = 0;
            }

            if (abs(k1) == 4 || abs(k1c) == 4 || abs(k1n) == 4) {       // match the atom names
                int iresc=0, iresn=0;
                // "0" in residx[] marks "present residue"
                for (l = 0; l < 5 && Gcmap->residx[l] != 0; l++);
                if (l < 4) {
                    ires =
                        cmpresname1(uUnit, sa4[l], Gcmap->reslist, Gcmap->nres);
                    iresc =
                        cmpresname1(uUnit, sa4[l], Gcmap->creslist, Gcmap->nres);
                    iresn =
                        cmpresname1(uUnit, sa4[l], Gcmap->nreslist, Gcmap->nres);

                    if (abs(k1) != 4)
                        ires = 0;
                    if (abs(k1c) != 4)
                        iresc = 0;
                    if (abs(k1n) != 4)
                        iresn = 0;
                }

                if (ires > 0 || iresc > 0 || iresn > 0 ||       // found on the reslist of the cmap
                    l == 4) {   // or the "present residue" pointer is the last one.
                    for (j = 0; j < nprospect; j++) {
                        stPTorsion2 = stdpt0[j].tp;
                        get4atoms(uUnit, stPTorsion2, sb4);
                        if (l == 4) {
                            ires =
                                cmpresname1(uUnit, sb4[3], Gcmap->reslist,
                                            Gcmap->nres);
                            iresc =
                                cmpresname1(uUnit, sb4[3], Gcmap->creslist,
                                            Gcmap->nres);
                            iresn =
                                cmpresname1(uUnit, sb4[3], Gcmap->nreslist,
                                            Gcmap->nres);
                        }

                        if (Gcmap->termmap == 0) {
                            iresc = 0;
                            iresn = 0;
                        }

                        k2 = cmp4vs4(sb4, (WRD *) (&Gcmap->atmname[1])); // check atmnames 1..4
                        k2c = cmp4vs4(sb4, (WRD *) (&Gcmap->catmname[1]));       // check atmnames 1..4
                        k2n = cmp4vs4(sb4, (WRD *) (&Gcmap->natmname[1]));       // check atmnames 1..4

                        if (abs(k2) == 4 && abs(k1) == 4 && ires > 0)
                            // check residue indicies, 0: present, -1: residue before, +1: residue after
                            if (cmp_residx(sa4, sb4, Gcmap->residx))     // we found two torsions and the Gcmap
                            {
                                copyatoms(phipsi[k].atoms, stPTorsion,
                                          stPTorsion2, k1, k2);
                                phipsi[k].mapid = mapid;
                                // mark as used
                                mapflag[mapid - 1] = 1;
                                k++;
                            }

                        if (abs(k2c) == 4 && abs(k1c) == 4 && iresc > 0)
                            if (cmp_residx(sa4, sb4, Gcmap->cresidx)) {
                                copyatoms(phipsi[k].atoms, stPTorsion,
                                          stPTorsion2, k1c, k2c);
                                phipsi[k].mapid = mapid;
                                // mark as used
                                mapflag[mapid - 1] = 1;
                                k++;
                            }

                        if (abs(k2n) == 4 && abs(k1n) == 4 && iresn > 0)
                            if (cmp_residx(sa4, sb4, Gcmap->nresidx)) {
                                copyatoms(phipsi[k].atoms, stPTorsion,
                                          stPTorsion2, k1n, k2n);
                                phipsi[k].mapid = mapid;
                                // mark as used
                                mapflag[mapid - 1] = 1;
                                k++;
                            }
                    }
                }
            }
        }
    }

    int mk = 0;
    for (i = 0; i < k; i++) {
        int flag = 1;

// Remove the redundant torsions.

        int l;
        if (i > 0)
            for (l = 0; l < i; l++) {
                if (phipsi[i].atoms[0] == phipsi[l].atoms[0]
                    && phipsi[i].atoms[1] == phipsi[l].atoms[1]
                    && phipsi[i].atoms[2] == phipsi[l].atoms[2]
                    && phipsi[i].atoms[3] == phipsi[l].atoms[3]
                    && phipsi[i].atoms[4] == phipsi[l].atoms[4]
                    && phipsi[i].mapid == phipsi[l].mapid)
                    flag = -1;
            }
//
        for (j = 0; j < 5; j++)
            if (phipsi[i].atoms[j] < 0)
                flag = -1;
        if (phipsi[i].mapid <= 0)
            flag = -1;
        if (flag > 0)
            mk++;
    }

    //wmap
    maptypes = 0;

    for (i = 0; i < GiCmapNum; i++) {
        if (mapflag[i] == 1) {
            mapidx[i] = maptypes;
            maptypes++;
        }
    }
    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG CMAP_COUNT");
    FortranWriteString("%COMMENT SIZE=2; MEMBERS=CMAP_KINDS,CMAP_TYPES");
    //FortranWriteString("%COMMENT");
//        for (cmntt=cmnt0; cmntt->next != NULL; cmntt=cmntt->next){
//            sTmp[0]='\0';
//            strcat(sTmp,"%COMMENT");
//            strcat(sTmp,cmntt->record);
//            FortranWriteString(sTmp);
    //fprintf(fpout,"%%COMMENT%s",cmntt->record);
//        }
    FortranWriteString("%FORMAT(2I8)");
//        sprintf(sTmp,"%8d%8d",mapcount,GiCmapNum);
//        sprintf(sTmp,"%8d%8d",k,maptypes);
    int CMAP_KINDS=mk, CMAP_TYPES=maptypes;
    sprintf(sTmp, "%8d%8d", mk, maptypes);
    FortranWriteString(sTmp);
    //FortranEndLine();

    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG CMAP_RESOLUTION");
    COMMENT_DIM("CMAP_TYPES", CMAP_TYPES);
    FortranWriteString("%FORMAT(20I4)");
    FortranFormat(20, "%4d");
    mapid = 0;
    for (cmaplstt = Gcmaplst; cmaplstt->next != NULL; cmaplstt = cmaplstt->next) {
        Gcmap = cmaplstt->cmap;
        if (mapflag[mapid])
            FortranWriteInt(Gcmap->resolution);
        mapid++;
    }
    FortranEndLine();

    mapid = 0;
    for (cmaplstt = Gcmaplst; cmaplstt->next != NULL; cmaplstt = cmaplstt->next) {
        Gcmap = cmaplstt->cmap;
        int msize;
        if (mapflag[mapid] == 1) {
            sprintf(sTmp, "%%FLAG CMAP_PARAMETER_%02d", mapidx[mapid] + 1);
            FortranFormat(1, "%-80s");
            FortranWriteString(sTmp);
            msize = Gcmap->resolution * Gcmap->resolution;
            sprintf(sComment,
                  "%%COMMENT DIMENSION=CMAP_RESOLUTION(%d),CMAP_RESOLUTION(%d); SIZE=%d",
                  mapidx[mapid]+1,mapidx[mapid]+1,msize);
            FortranWriteString(sComment);
            sprintf(sTmp, "%%COMMENT TITLE=%s", Gcmap->title);
            FortranWriteString(sTmp);
            FortranWriteString("%FORMAT(8F9.5)");
            FortranFormat(8, "%9.5lf");
            for (j = 0; j < msize; j += 8) {
                int l;
                for (l = j; l < msize && l < j + 8; l++) {
                    FortranWriteDouble(Gcmap->map[l]);
                }
            }
            FortranEndLine();
        }
        mapid++;
    }

    FortranFormat(1, "%-80s");
    FortranWriteString("%FLAG CMAP_INDEX");
    sprintf(sComment, "%%COMMENT DIMENSION=CMAP_KINDS; STRIDE=6; SIZE=%d; MEMBERS:ATOM_I,ATOM_J,ATOM_K,ATOM_L,ATOM_M,PARM_INDEX", CMAP_KINDS*6);
    FortranWriteString(sComment);
    FortranWriteString("%FORMAT(6I8)");
    FortranFormat(6, "%8d");
    for (i = 0; i < k; i++) {
        int flag = 1;

// Remove the redundant torsions.

        int l;
        if (i > 0)
            for (l = 0; l < i; l++) {
                if (phipsi[i].atoms[0] == phipsi[l].atoms[0]
                    && phipsi[i].atoms[1] == phipsi[l].atoms[1]
                    && phipsi[i].atoms[2] == phipsi[l].atoms[2]
                    && phipsi[i].atoms[3] == phipsi[l].atoms[3]
                    && phipsi[i].atoms[4] == phipsi[l].atoms[4]
                    && phipsi[i].mapid == phipsi[l].mapid)
                    flag = -1;
            }
//
        for (j = 0; j < 5; j++)
            if (phipsi[i].atoms[j] < 0)
                flag = -1;
        if (phipsi[i].mapid <= 0)
            flag = -1;
        if (flag > 0) {
            for (j = 0; j < 5; j++)
                FortranWriteInt(phipsi[i].atoms[j] / 3 + 1);
            FortranWriteInt(mapidx[phipsi[i].mapid - 1] + 1);
        }
    }
    FortranEndLine();
    FREE(mapflag);
    FREE(mapidx);
    FREE(phipsi);
    FREE(stdpt0);
}

#ifdef BINTRAJ
/* ================================================================
   6. SaveAmberParmCMAPNetcdf()
   NetCDF equivalent of SaveAmberParmCMAP().
   Differences from Fortran path:
   - CMAP_COUNT: 2-element NC_INT array, not formatted string
   - CMAP_PARAMETER_nn: 2D NC_DOUBLE grid [res][res]
   - CMAP_INDEX: raw atom indices, no /3+1 AMBERINDEX transform;
     reader applies AMBERINDEX at runtime same as torsions
   ================================================================ */
/* ================================================================
   SaveAmberParmCMAPNetcdf()

   Dimensions:
     CMAP_KINDS  — number of unique phi/psi interaction pairs (mk)
     CMAP_TYPES  — number of unique map parameter sets (maptypes)
     CMAP_RESOLUTION_nn  — scalar int, grid size (res x res) for maptype nn

   Index arrays:
     CMAP_INDEX_ATOMS(CMAP_KINDS, cmap_atom_quint)
       5 raw 1-based atom indices defining two adjacent torsions.
       Atoms 1-4 = phi torsion, atoms 2-5 = psi torsion.
       Shared central bond atoms 2-3 defines the phi/psi junction.
       No /3+1 AMBERINDEX transform — reader applies at runtime.
     CMAP_INDEX_MAP(CMAP_KINDS)
       1-based index into /cmap/map_nn/ subgroup.
   Parameter arrays:
     CMAP_PARAMETER_nn(CMAP_RESOLUTION_nn,CMAP_RESOLUTION_nn) — energy surface in kcal/mol
   ================================================================ */
// TODO don't duplicate all of the CMAP indexing code
int SaveAmberParmCMAPNetcdf(UNIT uUnit, int ncid)
{
    int i, j, k, l, mapid, maptypes, mk;
    int *mapflag, *mapidx;
    int iNumDIH, nprospect;
    SAVEATOMt *sa4[4], *sb4[4];
    SAVETORSIONt *stPTorsion, *stPTorsion2;
    SAVETORSIONt **stdpt0;
    PHIPSI *phipsi;
    CMAPLST *cmaplstt;
    CMAP *cmap;
    int k1, k2, k1c, k2c, k1n, k2n;
    int maxmap, ires, iresc, iresn;

    if (!GDefaults.iCMAP) return NC_NOERR;
    if (GiCmapNum <= 0)       return NC_NOERR;

    iNumDIH = iVarArrayElementCount(uUnit->vaTorsions);
    if (iNumDIH <= 0) return NC_NOERR;
    stPTorsion = PVAI(uUnit->vaTorsions, SAVETORSIONt, 0);

    MALLOC(mapflag, int *, sizeof(int) * (GiCmapNum + 1));
    MALLOC(mapidx,  int *, sizeof(int) * (GiCmapNum + 1));
    for (i = 0; i <= GiCmapNum; i++) { mapflag[i] = 0; mapidx[i] = 0; }

    maxmap = iNumDIH;
    MALLOC(phipsi, PHIPSI *, sizeof(PHIPSI) * maxmap);
    MALLOC(stdpt0, SAVETORSIONt **, sizeof(SAVETORSIONt *) * iNumDIH);
    nprospect = 0;

    /* pre-filter — identical to original */
    for (i = 0; i < iNumDIH; i++, stPTorsion++) {
        get4atoms(uUnit, stPTorsion, sa4);
        for (cmaplstt = Gcmaplst; cmaplstt->next != NULL; cmaplstt = cmaplstt->next) {
            cmap = cmaplstt->cmap;
            if (cmpresname4(uUnit, sa4, cmap->reslist, cmap->nres) > 0) {
                if (abs(cmp4vs4(sa4,(WRD*)(&cmap->atmname[0])))==4 ||
                    abs(cmp4vs4(sa4,(WRD*)(&cmap->atmname[1])))==4) {
                    stdpt0[nprospect++] = stPTorsion; break;
                }
            } else if (cmap->termmap > 0) {
                if (cmpresname4(uUnit, sa4, cmap->creslist, cmap->nres) > 0) {
                    if (abs(cmp4vs4(sa4,(WRD*)(&cmap->catmname[0])))==4 ||
                        abs(cmp4vs4(sa4,(WRD*)(&cmap->catmname[1])))==4) {
                        stdpt0[nprospect++] = stPTorsion; break;
                    }
                } else if (cmpresname4(uUnit, sa4, cmap->nreslist, cmap->nres) > 0)
                    if (abs(cmp4vs4(sa4,(WRD*)(&cmap->natmname[0])))==4 ||
                        abs(cmp4vs4(sa4,(WRD*)(&cmap->natmname[1])))==4) {
                        stdpt0[nprospect++] = stPTorsion; break;
                    }
            }
        }
    }

    /* match pairs — identical to original */
    k = 0; mapid = 0;
    for (cmaplstt = Gcmaplst; cmaplstt->next != NULL; cmaplstt = cmaplstt->next) {
        cmap = cmaplstt->cmap;
        mapid++;
        for (i = 0; i < nprospect; i++) {
            stPTorsion = stdpt0[i];
            get4atoms(uUnit, stPTorsion, sa4);
            k1  = cmp4vs4(sa4,(WRD*)(&cmap->atmname[0]));
            k1c = cmp4vs4(sa4,(WRD*)(&cmap->catmname[0]));
            k1n = cmp4vs4(sa4,(WRD*)(&cmap->natmname[0]));
            if (cmap->termmap == 0) { k1c = 0; k1n = 0; }
            if (abs(k1)!=4 && abs(k1c)!=4 && abs(k1n)!=4) continue;
            iresc = iresn = 0;
            for (l = 0; l < 5 && cmap->residx[l] != 0; l++);
            if (l < 4) {
                ires  = cmpresname1(uUnit, sa4[l], cmap->reslist,  cmap->nres);
                iresc = cmpresname1(uUnit, sa4[l], cmap->creslist, cmap->nres);
                iresn = cmpresname1(uUnit, sa4[l], cmap->nreslist, cmap->nres);
                if (abs(k1)  != 4) ires  = 0;
                if (abs(k1c) != 4) iresc = 0;
                if (abs(k1n) != 4) iresn = 0;
            } else ires = 1;
            if (ires <= 0 && iresc <= 0 && iresn <= 0 && l != 4) continue;
            for (j = 0; j < nprospect; j++) {
                stPTorsion2 = stdpt0[j];
                get4atoms(uUnit, stPTorsion2, sb4);
                if (l == 4) {
                    ires  = cmpresname1(uUnit, sb4[3], cmap->reslist,  cmap->nres);
                    iresc = cmpresname1(uUnit, sb4[3], cmap->creslist, cmap->nres);
                    iresn = cmpresname1(uUnit, sb4[3], cmap->nreslist, cmap->nres);
                }
                if (cmap->termmap == 0) { iresc = 0; iresn = 0; }
                k2  = cmp4vs4(sb4,(WRD*)(&cmap->atmname[1]));
                k2c = cmp4vs4(sb4,(WRD*)(&cmap->catmname[1]));
                k2n = cmp4vs4(sb4,(WRD*)(&cmap->natmname[1]));
                if (abs(k2)==4  && abs(k1)==4  && ires>0  && cmp_residx(sa4,sb4,cmap->residx))
                    { copyatoms_netcdf(phipsi[k].atoms,stPTorsion,stPTorsion2,k1,k2);
                      phipsi[k].mapid=mapid; mapflag[mapid-1]=1; k++; }
                if (abs(k2c)==4 && abs(k1c)==4 && iresc>0 && cmp_residx(sa4,sb4,cmap->cresidx))
                    { copyatoms_netcdf(phipsi[k].atoms,stPTorsion,stPTorsion2,k1c,k2c);
                      phipsi[k].mapid=mapid; mapflag[mapid-1]=1; k++; }
                if (abs(k2n)==4 && abs(k1n)==4 && iresn>0 && cmp_residx(sa4,sb4,cmap->nresidx))
                    { copyatoms_netcdf(phipsi[k].atoms,stPTorsion,stPTorsion2,k1n,k2n);
                      phipsi[k].mapid=mapid; mapflag[mapid-1]=1; k++; }
            }
        }
    }

    /* ── count unique (mk) and build mapidx ── */
    mk = 0;
    for (i = 0; i < k; i++) {
        int flag = 1;
        if (i > 0)
            for (l = 0; l < i; l++)
                if (phipsi[i].atoms[0]==phipsi[l].atoms[0] &&
                    phipsi[i].atoms[1]==phipsi[l].atoms[1] &&
                    phipsi[i].atoms[2]==phipsi[l].atoms[2] &&
                    phipsi[i].atoms[3]==phipsi[l].atoms[3] &&
                    phipsi[i].atoms[4]==phipsi[l].atoms[4] &&
                    phipsi[i].mapid  ==phipsi[l].mapid)      flag = -1;
        for (j = 0; j < 5; j++) if (phipsi[i].atoms[j] < 0) flag = -1;
        if (phipsi[i].mapid <= 0) flag = -1;
        if (flag > 0) mk++;
    }
    maptypes = 0;
    for (i = 0; i < GiCmapNum; i++)
        if (mapflag[i] == 1) { mapidx[i] = maptypes++; }

    if (mk == 0 || maptypes == 0) {
        FREE(mapflag); FREE(mapidx); FREE(phipsi); FREE(stdpt0);
        return NC_NOERR;
    }

    /* ── write CMAP_KINDS and CMAP_TYPES scalars ── */
    {

    /* ── CMAP_INDEX_ATOMS and CMAP_INDEX_MAP ── */

        int dimid_types, dimid_kinds, dimid_quint;
        int vid_atoms, vid_map;
        int dims2[2];
        int *atoms_buf = malloc(mk * 5 * sizeof(int));
        int *map_buf   = malloc(mk * sizeof(int));
        int idx = 0;

        NC_CHECK(nc_redef(ncid));

        NC_CHECK(nc_def_dim(ncid, "CMAP_TYPES", maptypes, &dimid_types));
        NC_CHECK(nc_def_dim(ncid, "CMAP_KINDS",     mk, &dimid_kinds));
        NC_CHECK(nc_def_dim(ncid, "atom_quint", 5, &dimid_quint));

        dims2[0] = dimid_kinds; dims2[1] = dimid_quint;
        NC_CHECK(nc_def_var(ncid, "CMAP_INDEX_ATOMS", NC_INT, 2, dims2, &vid_atoms));
        NC_CHECK(nc_put_att_text(ncid, vid_atoms, "long_name", 50,
                                 "5 raw 1-based atom indices: atoms 1-4=phi, 2-5=psi"));
        NC_CHECK(nc_put_att_text(ncid, vid_atoms, "reference_dimension", 6, "NTOTAT"));

        NC_CHECK(nc_def_var(ncid, "CMAP_INDEX_MAP", NC_INT, 1, &dimid_kinds, &vid_map));
        NC_CHECK(nc_put_att_text(ncid, vid_map, "long_name", 79,
                "1-based index into CMAP_RESOLUTION_<nn> and CMAP_PARAMETER_<nn> CMAP_TITLE_<nn>"));
        NC_CHECK(nc_put_att_text(ncid, vid_atoms, "reference_dimenson", 10, "CMAP_TYPES"));
        NC_CHECK(nc_enddef(ncid));

        for (i = 0; i < k; i++) {
            int flag = 1;
            if (i > 0)
                for (l = 0; l < i; l++)
                    if (phipsi[i].atoms[0]==phipsi[l].atoms[0] &&
                        phipsi[i].atoms[1]==phipsi[l].atoms[1] &&
                        phipsi[i].atoms[2]==phipsi[l].atoms[2] &&
                        phipsi[i].atoms[3]==phipsi[l].atoms[3] &&
                        phipsi[i].atoms[4]==phipsi[l].atoms[4] &&
                        phipsi[i].mapid  ==phipsi[l].mapid)      flag = -1;
            for (j = 0; j < 5; j++) if (phipsi[i].atoms[j] < 0) flag = -1;
            if (phipsi[i].mapid <= 0) flag = -1;
            if (flag > 0) {
                for (j = 0; j < 5; j++)
                    atoms_buf[idx*5 + j] = phipsi[i].atoms[j]; /* raw, no /3+1 */
                map_buf[idx] = mapidx[phipsi[i].mapid - 1] + 1;
                idx++;
            }
        }
        NC_CHECK(nc_put_var_int(ncid, vid_atoms, atoms_buf));
        NC_CHECK(nc_put_var_int(ncid, vid_map,   map_buf));
        free(atoms_buf); free(map_buf);
    }

    /* ── subgroups — one per active map type ──
       Each contains coordinate variables phi and psi (angle values
       in degrees) plus the GRID energy surface.                     */
    mapid = 0;
    for (cmaplstt = Gcmaplst; cmaplstt->next != NULL; cmaplstt = cmaplstt->next) {
        cmap = cmaplstt->cmap;
        if (mapflag[mapid] == 1) {
            int dimid_res, dims2[2];
            int vid_grid;
            int res = cmap->resolution;

            NC_CHECK(nc_redef(ncid));

            /* RESOLUTION_nn dimension shared by both phi and psi axes. */
            STRING sResolution;
            sprintf(sResolution,"CMAP_RESOLUTION_%02d",mapid+1);
            NC_CHECK(nc_def_dim(ncid, sResolution, res, &dimid_res));

            /* GRID(RESOLUTION, RESOLUTION) — same dim for both axes */
            dims2[0] = dims2[1] = dimid_res;
            STRING sParam;
            sprintf(sParam,"CMAP_PARAMETER_%02d",mapid+1);
            NC_CHECK(nc_def_var(ncid, sParam, NC_DOUBLE, 2, dims2, &vid_grid));
            NC_CHECK(nc_put_att_text(ncid, vid_grid, "units",      8, "kcal/mol"));
            NC_CHECK(nc_put_att_text(ncid, vid_grid, "long_name",  36,
                                     "CMAP energy grid: rows=phi, cols=psi"));
            NC_CHECK(nc_put_att_text(ncid, vid_grid, "coordinates", 18, sResolution));
            NC_CHECK(nc_put_att_text(ncid, vid_grid, "title", strlen(Gcmap->title), Gcmap->title));
            NC_CHECK(nc_enddef(ncid));

            NC_CHECK(nc_put_var_double(ncid, vid_grid, cmap->map));
        }
        mapid++;
    }

    FREE(mapflag); FREE(mapidx); FREE(phipsi); FREE(stdpt0);
    return NC_NOERR;
}
#endif
