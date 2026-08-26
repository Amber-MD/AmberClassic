/*
************************************************************************
*           All Copyright Reserved!                                    *
*                                                                      *
*  Prog:    parmchk2                                                   *
*  Version: version 2.0                                                *
*  Author:  Junmei Wang                                                *
*                                                                      *
*  Department of Pharmaceutical Chemistry                              *
*  School of Pharmacy                                                  *
*  University of California                                            *
*  San Francisco   CA 94143                                            *
*  October, 2001                                                       *
************************************************************************
* Refactored 2026 by Juno Krahn
* TODO: this is really too large for a single monolithic source file
* and the other soruce files should be compiles not inlined.
* bool + true/false should be used where relevant e.g. aliph, saturate, improper, ifreeze
*/

#include <unistd.h>
#include <stdbool.h>
#include "common.h"
#include "define.h"
#include "atom.h"

char *amberhome;
#include "utility.c"
#include "common.c"
#include "rotate.c"
#include "ac.c"
#include "mol2.c"
#include "prep.c"

#define MAXPARM           250
#define MAXEQUA            25
#define MAXCORR            50
#define MAX_FF_ATOMTYPE   250
#define MAX_FF_VDW        250
#define MAX_FF_BOND      2000
#define MAX_FF_ANGLE     5000
#define MAX_FF_TORSION   1500
#define MAX_FF_IMPROPER   500
#define MAX_EQU_VDW        20
#define MAX_EQU_VDW_NUM    50
#define INITSCORE      999999
#define NAMELEN             8

#ifdef DEBUG
int debug = 0;
#else
#define debug 0
#endif

#define MAXBEST 5
typedef struct { int idx; double score; } RankedMatch;

/* Running insertion into a score-sorted (ascending) array of at most
 * maxbest entries.  One step of insertion sort per call — O(MAXBEST).  */
#define RANKED_INSERT(best, nbest, new_idx, new_score)              \
    do {                                                            \
        if ((nbest) < maxbest ||                                    \
            (new_score) < (best)[(nbest)-1].score) {               \
            int _pos = ((nbest) < maxbest) ? (nbest) : (nbest)-1;  \
            while (_pos > 0 &&                                      \
                   (new_score) < (best)[_pos-1].score) {            \
                (best)[_pos] = (best)[_pos-1];                      \
                _pos--;                                             \
            }                                                       \
            (best)[_pos].idx   = (new_idx);                         \
            (best)[_pos].score = (new_score);                       \
            if ((nbest) < maxbest) (nbest)++;                       \
        }                                                           \
    } while(0)

/*------------------------------------------------------------------------
 * GROW_ARRAY: type-plastic realloc with zeroing of new chunk.
 * memset matches calloc behaviour. Allocations in chunk number of objects.
 *----------------------------------------------------------------------*/
#define GROW_ARRAY(ptr, count, max, chunk, type)                        \
    do {                                                                \
        if ((count) >= (max)) {                                         \
            int _old = (max);                                           \
            (max) += (chunk);                                           \
            (ptr) = (type*)realloc((ptr), sizeof(type)*(max));         \
            if (!(ptr)) {                                               \
                fprintf(stdout, "memory allocation error for *%s in %s:%d\n", \
                        #ptr, __FILE__, __LINE__);                      \
                exit(1);                                                \
            }                                                           \
            memset((ptr)+_old, 0, sizeof(type)*(chunk));               \
        }                                                               \
    } while(0)

#define VLOG(level, ...) do { if (verbose>=(level)) fprintf(stderr, __VA_ARGS__); } while(0)

/*========================================================================
 * Enumerations
 *======================================================================*/
typedef enum {
    FILETYPE_NONE=0, FILETYPE_PREPI, FILETYPE_PREPC, FILETYPE_AC,
    FILETYPE_MOL2, FILETYPE_FRCMOD, FILETYPE_LEAPLOG
} FILETYPE;

typedef enum {
    FFTYPE_NONE=0, FFTYPE_GAFF, FFTYPE_GAFF2, FFTYPE_PARM99,
    FFTYPE_PARM10, FFTYPE_LIPID14, FFTYPE_GLYCAM
} FFTYPE;

typedef enum {
    FRCMOD_COUNT=0, FRCMOD_READ_R, FRCMOD_READ_MAIN
} FrcmodMode;

/*========================================================================
 * Structs
 *======================================================================*/
typedef struct { int atid1, atid2, atid3, atid4; } IMPROPERID;
typedef struct { char atomtype[NAMELEN]; int pid; } EQUA;

typedef struct {
    char   atomtype[NAMELEN];
    int    pid, type;
    double bl, blf, ba, baf, cba, cbaf, ctor, tor, ps, improper;
} CORR;

typedef struct {
    char   atomtype[NAMELEN];
    int    group, equtype, improper, ncorr, nequa, atomicnum;
    double mass;
    CORR   corr[MAXCORR];
    EQUA   equa[MAXEQUA];
} PARM;

typedef struct { char name[NAMELEN]; double mass, pol; int attn; } ATOMTYPE;

typedef struct {
    char name1[NAMELEN], name2[NAMELEN], name3[NAMELEN];
    double angle, force; int attn;
} ANGLE;

typedef struct {
    char name1[NAMELEN], name2[NAMELEN], name3[NAMELEN], name4[NAMELEN];
    int mul; double fterm, phase, force; int attn;
} TORSION;

typedef struct {
    char name1[NAMELEN], name2[NAMELEN], name3[NAMELEN], name4[NAMELEN];
    int mul, numX; double fterm, phase, force; int attn;
} IMPROPER;

typedef struct { char name[NAMELEN]; double radius, pot; int attn; } VDW;
typedef struct { char name[MAX_EQU_VDW_NUM][NAMELEN]; int num; } EQU_VDW;

typedef struct {
    double BL, BLF, BA, BAF, X, X3, BA_CTR, TOR_CTR, IMPROPER, GROUP, EQUTYPE;
} WT;

typedef struct {
    double BL, BLF, BA, BAF, BA_CTR, BAF_CTR, TOR, TOR_CTR, FRACT1, FRACT2;
} DEFAULT_VALUE;

/* Query structs — atid* are atom indices (v1 topology path) or -1 (v2). */
typedef struct {
    char name1[NAMELEN], name2[NAMELEN];
    double length; int attn; int atid1, atid2;
} BondQuery;

typedef struct {
    char name1[NAMELEN], name2[NAMELEN], name3[NAMELEN];
    double bl1, bl2, angle; int attn; int atid1, atid2, atid3;
} AngleQuery;

typedef struct {
    char name1[NAMELEN], name2[NAMELEN], name3[NAMELEN], name4[NAMELEN];
    int attn; int atid1, atid2, atid3, atid4;
} TorsionQuery;

typedef struct {
    char name1[NAMELEN], name2[NAMELEN], name3[NAMELEN], name4[NAMELEN];
    int numX, attn; int atid1, atid2, atid3, atid4;
} ImproperQuery;

/*========================================================================
 * Globals
 *======================================================================*/
MOLINFO minfo; CONTROLINFO cinfo;
ATOM *atom; BOND *bond_array;
WT wt; DEFAULT_VALUE dv;

int atomnum=0, bondnum=0;
IMPROPERID *improper;

char ifilename[MAXCHAR]="", ofilename[MAXCHAR]="", pfilename[MAXCHAR]="";
char cfilename[MAXCHAR]="", atcfilename[MAXCHAR]="";
char additional_frcmod_file[MAXCHAR]="", frcmod_str[MAXCHAR]="";

FFTYPE ffset=FFTYPE_GAFF;
int gaff_set=0, impropernum=0, output_improper_flag=1, attn_opt=1;
int allparm_flag=0, pformat=1;
FILETYPE iformat=FILETYPE_NONE;
int fc_opt=1;
int verbose=0;           /* 0=off, 1=summary, 2=full step tracing */
int equa_corr_inherit=0; /* -icor: inherit CORR from EQUA types   */
int maxbest=1;           /* -multi N: top N CORR matches           */

int maxatomtype, maxvdwparm, maxbondparm, maxangleparm;
int maxtorsionparm, maximproperparm;
int atomtypenum=0, vdwparmnum=0, bondparmnum=0;
int angleparmnum=0, torsionparmnum=0, improperparmnum=0;
int atomtypenum2=0, vdwparmnum2=0, bondparmnum2=0;
int angleparmnum2=0, torsionparmnum2=0, improperparmnum2=0;
int last_atomtypenum=0, last_vdwparmnum=0, last_bondparmnum=0;
int last_angleparmnum=0, last_torsionparmnum=0, last_improperparmnum=0;
int afrcmod_atomtypenum_beg, afrcmod_bondparmnum_beg;
int afrcmod_angleparmnum_beg, afrcmod_torsionparmnum_beg;
int afrcmod_improperparmnum_beg, afrcmod_vdwparmnum_beg;
int afrcmod_atomtypenum_end, afrcmod_bondparmnum_end;
int afrcmod_angleparmnum_end, afrcmod_torsionparmnum_end;
int afrcmod_improperparmnum_end, afrcmod_vdwparmnum_end;

double THRESHOLD_BA, equtype_penalty_score;

int maxparmnum, parmnum=0, parmnum2=0;
int *parmid;
PARM *parm; ATOMTYPE *atomtype;
BOND_FF *bondparm; ANGLE *angleparm;
TORSION *torsionparm; IMPROPER *improperparm; VDW *vdwparm;

EQU_VDW equ_vdw[MAX_EQU_VDW]; int equ_vdw_num=0;

ATOMTYPE *r_atomtype; BOND_FF *r_bondparm; ANGLE *r_angleparm;
TORSION *r_torsionparm; IMPROPER *r_improperparm; VDW *r_vdwparm;
int r_atomtype_num=0, r_bondparm_num=0, r_angleparm_num=0;
int r_torsionparm_num=0, r_improperparm_num=0, r_vdwparm_num=0;

FILE *fpout;
ANGLE bestba;
int bestblid=-1, bestbaid=-1, besttorid=-1, bestimproperid=-1;
int nblf_parm=0;
char *remark_line;

/*========================================================================
 * Small helpers
 *======================================================================*/
static inline void copy_type(char *dest, const char *src) {
    dest[0]=src[0];
    if (src[1]==' '||src[1]=='\0'||src[1]==')'||src[1]=='-') dest[1]='\0';
    else { dest[1]=src[1]; dest[2]='\0'; }
}

static void strip_ansi_and_print(const char *s) {
    while (*s) {
        if (*s=='\033'&&*(s+1)=='[') { s+=2; while(*s&&*s!='m')s++; if(*s)s++; }
        else putchar(*s++);
    }
}

static void print_help(void) {
    static const char help[]=
        "\033[31mUsage: parmchk2 -i   \033[0m input file name\n"
        "\033[31m                -o   \033[0m frcmod file name\n"
        "\033[31m                -f   \033[0m input file format (prepi, prepc, ac, mol2, frcmod, leaplog)\n"
        "\033[31m                -s   \033[0m ff parm set, suppressed by \"-p\"\n"
        "\033[34m                      1 or gaff:   \033[0m gaff (the default)\n"
        "\033[34m                      2 or gaff2:  \033[0m gaff2\n"
        "\033[34m                      3 or parm99: \033[0m parm99\n"
        "\033[34m                      4 or parm10: \033[0m parm10\n"
        "\033[34m                      5 or lipid14:\033[0m lipid14\n"
        "\033[34m                      6 or glycam: \033[0m GLYCAM_06j\n"
        "\033[31m                -frc \033[0m frcmod files to load; supported names:\n"
        "\033[31m                     \033[0m ff99SB, ff14SB, ff03 for proteins\n"
        "\033[31m                     \033[0m bsc1, ol15, ol3 for DNA; yil for RNA\n"
        "\033[31m                     \033[0m eg. ff14SB+bsc1+yil, ff99SB+bsc1\n"
        "\033[31m                -p   \033[0m parmfile, suppresses -s, optional\n"
        "\033[31m                -pf  \033[0m parmfile format\n"
        "\033[34m                      1:\033[0m amber FF data file (the default)\n"
        "\033[34m                      2:\033[0m additional frcmod-format parameter file\n"
        "\033[31m                -afrc\033[0m additional frcmod file, optional\n"
        "\033[31m                -c   \033[0m atom type corresponding score file (default: PARMCHK.DAT)\n"
        "\033[31m                -atc \033[0m additional atom type corresponding score file, optional\n"
        "\033[31m                     \033[0m type 'parmchk2 -l' for format details\n"
        "\033[31m                -a   \033[0m print all FF parameters including those already in parmfile\n"
        "\033[31m                     \033[0m Y/N (default N)\n"
        "\033[31m                -w   \033[0m print improper parameters matched via X wildcard\n"
        "\033[31m                     \033[0m Y/N (default Y)\n"
        "\033[31m                -fc  \033[0m force constant calculation for frcmod/leaplog input\n"
        "\033[34m                      1:\033[0m default behaviour\n"
        "\033[34m                      2:\033[0m empirical calculation before corresponding atom types\n"
        "\033[31m                -att \033[0m scope of parmchk for frcmod input\n"
        "\033[34m                      1:\033[0m all parameters (default)\n"
        "\033[34m                      2:\033[0m only parameters flagged ATTN\n"
        "New:\n"
        "\033[31m                -icor\033[0m inherit CORR tables from EQUA equivalent types\n"
        "\033[31m                     \033[0m scores adjusted by WEIGHT_EQUTYPE from PARMCHK.DAT\n"
        "\033[31m                     \033[0m useful for adding simple EQUA-only tables; off by default\n"
        "\033[31m                -multi\033[0m report top N CORR matches (default 1)\n"
        "\033[31m                      \033[0m winner written as parameter, rest as # alt[] comments\n"
        "\033[31m                -v   \033[0m verbose logging to stderr\n"
        "\033[34m                      0:\033[0m off (default)\n"
        "\033[34m                      1:\033[0m summary/brief\n"
        "\033[34m                      2:\033[0m verbose info\n";

    int color=(!strcmp(COLORTEXT,"YES")||!strcmp(COLORTEXT,"yes"));
    if (color) printf("%s",help); else strip_ansi_and_print(help);
}

static FILETYPE parse_iformat(const char *s) {
    if (!strcmp(s,"prepi"))                        return FILETYPE_PREPI;
    if (!strcmp(s,"prepc"))                        return FILETYPE_PREPC;
    if (!strcmp(s,"ac"))                           return FILETYPE_AC;
    if (!strcmp(s,"mol2"))                         return FILETYPE_MOL2;
    if (!strcmp(s,"frcmod")||!strcmp(s,"frc"))     return FILETYPE_FRCMOD;
    if (!strcmp(s,"leaplog")||!strcmp(s,"leap"))   return FILETYPE_LEAPLOG;
    fprintf(stderr,"Unknown input format: %s\n",s); exit(1);
}

static FFTYPE parse_ffset(const char *s) {
    char low[64]; int len=(int)strlen(s); if(len>63)len=63;
    for(int ii=0;ii<len;ii++) low[ii]=tolower((unsigned char)s[ii]);
    low[len]='\0';
    if(!strcmp(low,"gaff")   ||!strcmp(s,"1")) return FFTYPE_GAFF;
    if(!strcmp(low,"gaff2")  ||!strcmp(s,"2")) return FFTYPE_GAFF2;
    if(!strcmp(low,"parm99") ||!strcmp(s,"3")) return FFTYPE_PARM99;
    if(!strcmp(low,"parm10") ||!strcmp(s,"4")) return FFTYPE_PARM10;
    if(!strcmp(low,"lipid14")||!strcmp(s,"5")) return FFTYPE_LIPID14;
    if(!strcmp(low,"glycam") ||!strcmp(s,"6")) return FFTYPE_GLYCAM;
    fprintf(stderr,"Unknown force field set: %s\n",s); exit(1);
}

static void update_last_counts(void) {
    last_atomtypenum=atomtypenum; last_bondparmnum=bondparmnum;
    last_angleparmnum=angleparmnum; last_torsionparmnum=torsionparmnum;
    last_improperparmnum=improperparmnum; last_vdwparmnum=vdwparmnum;
}

static void init_corr(CORR *c, const char *at, int type, int pid) {
    strcpy(c->atomtype,at);
    c->bl=c->blf=c->ba=c->baf=0.0;
    c->cba=c->cbaf=c->tor=c->ctor=c->improper=c->ps=0.0;
    c->type=type; c->pid=pid;
}

static int lookup_parmid(const char *at) {
    for (int ii=0;ii<parmnum;ii++)
        if (!strcmp(parm[ii].atomtype,at)) return ii;
    fprintf(stderr,"Atom type %s does not exist in %s\n",at,cfilename);
    if (atcfilename[0]) fprintf(stderr,"\t\t or in %s\n",atcfilename);
    exit(1);
}

static inline int get_pid(int atid, const char *name) {
    return (atid>=0) ? parmid[atid] : lookup_parmid(name);
}

/*========================================================================
 * build_ext_corr: extend CORR tables for user-supplied (parmnum2..parmnum)
 * atom types that have only EQUA entries, by appending the CORR tables of
 * their EQUA-matched pairs, with scores adjusted by wt.EQUTYPE.
 * Only called when -icor is active.  Assumes EQUA-only -atc input table.
 *======================================================================*/
static void build_ext_corr(void) {
    for (int pid=parmnum2; pid<parmnum; pid++) {
        int n=parm[pid].ncorr;
        for (int e=1; e<parm[pid].nequa && n<MAXCORR; e++) {
            int epid=parm[pid].equa[e].pid;
            if (!parm[epid].ncorr) continue;
            VLOG(2,"  %s inherits CORR from EQUA atom %s:",
                 parm[pid].atomtype,parm[epid].atomtype);
            for (int j=0; j<parm[epid].ncorr && n<MAXCORR; j++) {
                parm[pid].corr[n]=parm[epid].corr[j];
                parm[pid].corr[n].bl      +=wt.EQUTYPE;
                parm[pid].corr[n].blf     +=wt.EQUTYPE;
                parm[pid].corr[n].ba      +=wt.EQUTYPE;
                parm[pid].corr[n].baf     +=wt.EQUTYPE;
                parm[pid].corr[n].cba     +=wt.EQUTYPE;
                parm[pid].corr[n].cbaf    +=wt.EQUTYPE;
                parm[pid].corr[n].tor     +=wt.EQUTYPE;
                parm[pid].corr[n].ctor    +=wt.EQUTYPE;
                parm[pid].corr[n].improper+=wt.EQUTYPE;
                VLOG(2," %s",parm[parm[epid].corr[j].pid].atomtype);
                n++;
            }
            VLOG(2,"\n");
        }
        parm[pid].ncorr=n;
    }
}

/*========================================================================
 * improper_id1 / improper_id2
 *======================================================================*/
void improper_id1(char *filename) {
    FILE  *fpin;
    int    readindex = 0;
    char   t1[MAXCHAR], t2[MAXCHAR], t3[MAXCHAR], t4[MAXCHAR];
    char line[MAXCHAR];

    if ((fpin = fopen(filename, "r")) == NULL) {
        fprintf(stdout, "Cannot open %s in improper_id1(), exit\n", filename);
        exit(1);
    }
    while (fgets(line, MAXCHAR, fpin)) {
        sscanf(line, "%s", t1);
        if (!strcmp("IMPROPER", t1)) { readindex = 1; continue; }
        if (readindex == 1) {
            if (spaceline(line) == 1      || !strcmp(t1, "DONE")   ||
                !strcmp(t1, "STOP")       || !strcmp(t1, "CHARGE") ||
                !strcmp(t1, "LOOP")       || !strcmp(t1, "IMPROPER"))
                { readindex = 0; continue; }
            sscanf(line, "%s%s%s%s", t1, t2, t3, t4);
            if (!strcmp(t1,"-M")||!strcmp(t2,"-M")||!strcmp(t3,"-M")||!strcmp(t4,"-M")) continue;
            if (!strcmp(t1,"+M")||!strcmp(t2,"+M")||!strcmp(t3,"+M")||!strcmp(t4,"+M")) continue;
            printf("IMPR(%s,%s,%s,%s) FIXME!!!\n", t1, t2, t3, t4);
            for (int i = 0; i < atomnum; i++) {
                if (!strcmp(t1, atom[i].name)) { strcpy(t1, atom[i].ambername); improper[impropernum].atid1 = i; continue; }
                if (!strcmp(t2, atom[i].name)) { strcpy(t2, atom[i].ambername); improper[impropernum].atid2 = i; continue; }
                if (!strcmp(t3, atom[i].name)) { strcpy(t3, atom[i].ambername); improper[impropernum].atid3 = i; continue; }
                if (!strcmp(t4, atom[i].name)) { strcpy(t4, atom[i].ambername); improper[impropernum].atid4 = i; continue; }
            }
            impropernum++;
        }
    }
    fclose(fpin);
}
void improper_id2(void) {
    impropernum = 0;
    for (int i = 0; i < atomnum; i++)
        if (atom[i].improper == 1) {
            if (atom[i].con[0] < 0 || atom[i].con[1] < 0 || atom[i].con[2] < 0)
                continue;
            improper[impropernum].atid3 = i;
            improper[impropernum].atid1 = atom[i].con[0];
            improper[impropernum].atid2 = atom[i].con[1];
            improper[impropernum].atid4 = atom[i].con[2];
            impropernum++;
        }
}

/*========================================================================
 * assign_parmid
 *======================================================================*/
void assign_parmid(void) {
    for (int i=0; i<atomnum; i++) {
        parmid[i] = lookup_parmid(atom[i].ambername);
        if (parm[parmid[i]].improper == 1)
            atom[i].improper = 1;
    }
}

/*========================================================================
 * equtype_penalty
 *======================================================================*/
double equtype_penalty(int id1, int id2, int id3, int id4) {
    int n1=parm[id1].equtype, n2=parm[id2].equtype;
    int n3=parm[id3].equtype, n4=parm[id4].equtype;
    if (n1==0 && n2==0) return 0.0;
    int tn1=abs(n1)+abs(n2), tn2=abs(n3)+abs(n4);
    if (tn2==3 && tn1!=3) return wt.EQUTYPE;
    if (tn1==3 && tn2!=3) return wt.EQUTYPE;
    if ((n1+n2)==0 && n3<0 && n4<0) return 0.5*wt.EQUTYPE;
    return 0.0;
}

/*========================================================================
 * check_bond_length / empbond / empangle
 *======================================================================*/
double check_bond_length(const char *name1, const char *name2) {
    for (int kk=0; kk<bondparmnum; kk++)
        if ((!strcmp(bondparm[kk].name1,name1)&&!strcmp(bondparm[kk].name2,name2))||
            (!strcmp(bondparm[kk].name1,name2)&&!strcmp(bondparm[kk].name2,name1)))
            return bondparm[kk].length;
    return 0.0;
}

static void ensure_blba_loaded(void) {
    char blba_parmfile[MAXCHAR];
    static int iread_blbaparm = 0;
    if (iread_blbaparm) return;
    iread_blbaparm = 1;
    if (gaff_set==2)
        build_dat_path(blba_parmfile,"PARM_BLBA_GAFF2.DAT",sizeof blba_parmfile,0);
    else
        build_dat_path(blba_parmfile,"PARM_BLBA_GAFF.DAT", sizeof blba_parmfile,0);
    read_blba_parm(blba_parmfile, blf_parm, &nblf_parm, baf_parm);
}

int empbond(const char *type1, const char *type2,
            const char *name1, const char *name2,
            int id1, int id2, double input_bl) {
    int    flag=0, pid1=-1, pid2=-1;
    double bl=0.0, force, kref=-1, rref=-1, delta_r1, delta_r2;

    ensure_blba_loaded();

    bl = input_bl;
    if (bl<=0)
        for (int kk=0; kk<bondparmnum; kk++)
            if ((!strcmp(bondparm[kk].name1,type1)&&!strcmp(bondparm[kk].name2,type2))||
                (!strcmp(bondparm[kk].name1,type2)&&!strcmp(bondparm[kk].name2,type1)))
                { bl=bondparm[kk].length; break; }
    if (bl<=0) return 0;

    for (int kk=0; kk<nblf_parm; kk++)
        if ((blf_parm[kk].id1==id1&&blf_parm[kk].id2==id2)||
            (blf_parm[kk].id1==id2&&blf_parm[kk].id2==id1))
            { kref=blf_parm[kk].bfkai; rref=blf_parm[kk].refbondlength; break; }

    if (kref<0 || rref<0) {
        for (int kk=0;kk<nblf_parm;kk++) if(blf_parm[kk].id1==id1&&blf_parm[kk].id2==id1){pid1=kk;break;}
        for (int kk=0;kk<nblf_parm;kk++) if(blf_parm[kk].id1==id2&&blf_parm[kk].id2==id2){pid2=kk;break;}
        if (pid1<0||pid2<0) return 0;
        if (rref<0) { rref=0.5*(blf_parm[pid1].refbondlength+blf_parm[pid2].refbondlength); flag=1; }
        delta_r1=fabs(rref-blf_parm[pid1].refbondlength);
        delta_r2=fabs(rref-blf_parm[pid2].refbondlength);
        kref=(blf_parm[pid1].bfkai*delta_r1+blf_parm[pid2].bfkai*delta_r2)/(delta_r1+delta_r2);
    }
    if (kref<=0) return 0;

    force = exp(kref - blf_exp_const*log(bl));
    fprintf(fpout,"%-2s-%-2s%8.2lf%8.3lf", name1, name2, force, bl);
    if (!flag)
        fprintf(fpout,"       calculated for %2s-%2s\n", name1, name2);
    else
        fprintf(fpout,"Warning: reference bond length of %s-%s (%8.4lf) calculated "
                "using %s-%s (%7.4lf) and %s-%s (%7.4lf)\n",
                blf_parm[pid1].elem1, blf_parm[pid2].elem1, rref,
                blf_parm[pid1].elem1, blf_parm[pid1].elem2, blf_parm[pid1].refbondlength,
                blf_parm[pid2].elem1, blf_parm[pid2].elem2, blf_parm[pid2].refbondlength);

    strcpy(bondparm[bondparmnum].name1,name1);
    strcpy(bondparm[bondparmnum].name2,name2);
    bondparm[bondparmnum].force  = force;
    bondparm[bondparmnum].length = bl;
    bondparmnum++;
    GROW_ARRAY(bondparm,bondparmnum,maxbondparm,MAX_FF_BOND,BOND_FF);
    return 1;
}

int empangle(const char *type1, const char *type2, const char *type3,
             const char *name1, const char *name2, const char *name3,
             int id1, int id2, int id3,
             double input_bl1, double input_bl2, double input_angle, int flag) {
    int    num1=-1, num2=-1;
    double bl1=0.0, bl2=0.0, cparm, dparm, zparm1, zparm2, angle, force;

    ensure_blba_loaded();

    if (input_angle<=0) {
        for (int kk=0;kk<angleparmnum;kk++)
            if (!strcmp(angleparm[kk].name1,type1)&&!strcmp(angleparm[kk].name2,type2)&&
                !strcmp(angleparm[kk].name3,type1)) { num1=kk; break; }
        if (num1==-1) return 0;
        for (int kk=0;kk<angleparmnum;kk++)
            if (!strcmp(angleparm[kk].name1,type3)&&!strcmp(angleparm[kk].name2,type2)&&
                !strcmp(angleparm[kk].name3,type3)) { num2=kk; break; }
        if (num2==-1) return 0;
        angle=0.5*(angleparm[num1].angle+angleparm[num2].angle);
    } else {
        angle=input_angle;
    }
    if (angle<=0) return 0;

    bl1=(input_bl1>0)?input_bl1:check_bond_length(type1,type2);
    if (bl1==0.0) return 0;
    bl2=(input_bl2>0)?input_bl2:check_bond_length(type2,type3);
    if (bl2==0.0) return 0;

    zparm1=baf_parm[id1].anglez;
    zparm2=baf_parm[id3].anglez;
    cparm =baf_parm[id2].anglec;
    dparm =(bl1-bl2)*(bl1-bl2)/((bl1+bl2)*(bl1+bl2));
    force =143.9*zparm1*cparm*zparm2*exp(-2*dparm)/(bl1+bl2);
    force/=sqrt(angle*3.1415926/180.0);

    if (flag==0) {
        fprintf(fpout,"%-2s-%-2s-%-2s%9.3lf%12.3lf", name1,name2,name3,force,angle);
        fprintf(fpout,"   Calculated with empirical approach for %s-%s-%s\n",type1,type2,type3);
        strcpy(angleparm[angleparmnum].name1,name1);
        strcpy(angleparm[angleparmnum].name2,name2);
        strcpy(angleparm[angleparmnum].name3,name3);
        angleparm[angleparmnum].angle=angle;
        angleparm[angleparmnum].force=force;
        angleparmnum++;
        GROW_ARRAY(angleparm,angleparmnum,maxangleparm,MAX_FF_ANGLE,ANGLE);
    } else {
        bestba.angle=angle; bestba.force=force;
        strcpy(bestba.name1,type1); strcpy(bestba.name2,type2); strcpy(bestba.name3,type3);
        bestbaid=999999;
    }
    return 1;
}

/*========================================================================
 * Query collectors
 *======================================================================*/
static int collect_bond_queries(BondQuery **out) {
    int idx;
    char tmp[NAMELEN];
    if (iformat<=FILETYPE_MOL2) {
        int n=0;
        for (int ii=0;ii<atomnum;ii++)
            for (int jj=ii+1;jj<atomnum;jj++)
                for (int kk=0;kk<6;kk++) if(atom[ii].con[kk]==jj){n++;break;}
        BondQuery *q=malloc(sizeof(BondQuery)*n);
        idx=0;
        for (int ii=0;ii<atomnum;ii++)
            for (int jj=ii+1;jj<atomnum;jj++)
                for (int kk=0;kk<6;kk++) if(atom[ii].con[kk]==jj){
                    strcpy(q[idx].name1,atom[ii].ambername);
                    strcpy(q[idx].name2,atom[jj].ambername);
                    if(strcmp(q[idx].name1,q[idx].name2)>0){
                        strcpy(tmp,q[idx].name2);strcpy(q[idx].name2,q[idx].name1);strcpy(q[idx].name1,tmp);
                    }
                    q[idx].length=0; q[idx].attn=1;
                    q[idx].atid1=ii; q[idx].atid2=jj;
                    idx++; break;
                }
        *out=q; return idx;
    } else {
        BondQuery *q=malloc(sizeof(BondQuery)*r_bondparm_num);
        for (int ii=0;ii<r_bondparm_num;ii++){
            strcpy(q[ii].name1,r_bondparm[ii].name1);
            strcpy(q[ii].name2,r_bondparm[ii].name2);
            if(strcmp(q[ii].name1,q[ii].name2)>0){
                strcpy(tmp,q[ii].name2);strcpy(q[ii].name2,q[ii].name1);strcpy(q[ii].name1,tmp);
            }
            q[ii].length=r_bondparm[ii].length;
            q[ii].attn  =r_bondparm[ii].attn;
            q[ii].atid1=q[ii].atid2=-1;
        }
        *out=q; return r_bondparm_num;
    }
}

static int collect_angle_queries(AngleQuery **out) {
    int idx;
    char tmp[NAMELEN];
    if (iformat<=FILETYPE_MOL2) {
        int n=0;
        for (int ii=0;ii<atomnum;ii++) for (int jj=0;jj<atomnum;jj++) for (int kk=0;kk<atomnum;kk++){
            if(ii==kk) continue;
            int ij=0,jk=0;
            for(int c=0;c<6;c++) if(atom[ii].con[c]==jj){ij=1;break;}
            for(int c=0;c<6;c++) if(atom[jj].con[c]==kk){jk=1;break;}
            if(ij&&jk) n++;
        }
        AngleQuery *q=malloc(sizeof(AngleQuery)*n);
        idx=0;
        for (int ii=0;ii<atomnum;ii++) for (int jj=0;jj<atomnum;jj++) for (int kk=0;kk<atomnum;kk++){
            if(ii==kk) continue;
            int ij=0,jk=0;
            for(int c=0;c<6;c++) if(atom[ii].con[c]==jj){ij=1;break;}
            for(int c=0;c<6;c++) if(atom[jj].con[c]==kk){jk=1;break;}
            if(!ij||!jk) continue;
            strcpy(q[idx].name1,atom[ii].ambername);
            strcpy(q[idx].name2,atom[jj].ambername);
            strcpy(q[idx].name3,atom[kk].ambername);
            if(strcmp(q[idx].name1,q[idx].name3)>0){
                strcpy(tmp,q[idx].name3);strcpy(q[idx].name3,q[idx].name1);strcpy(q[idx].name1,tmp);
            }
            /* look up bond lengths now while we have atom names — used by empangle */
            q[idx].bl1  = check_bond_length(atom[ii].ambername, atom[jj].ambername);
            q[idx].bl2  = check_bond_length(atom[jj].ambername, atom[kk].ambername);
            q[idx].angle= 0; q[idx].attn=1;
            q[idx].atid1=ii; q[idx].atid2=jj; q[idx].atid3=kk;
            idx++;
        }
        *out=q; return idx;
    } else {
        AngleQuery *q=malloc(sizeof(AngleQuery)*r_angleparm_num);
        for (int ii=0;ii<r_angleparm_num;ii++){
            strcpy(q[ii].name1,r_angleparm[ii].name1);
            strcpy(q[ii].name2,r_angleparm[ii].name2);
            strcpy(q[ii].name3,r_angleparm[ii].name3);
            q[ii].bl1  =check_bond_length(r_angleparm[ii].name1,r_angleparm[ii].name2);
            q[ii].bl2  =check_bond_length(r_angleparm[ii].name2,r_angleparm[ii].name3);
            q[ii].angle=r_angleparm[ii].angle;
            q[ii].attn =r_angleparm[ii].attn;
            q[ii].atid1=q[ii].atid2=q[ii].atid3=-1;
        }
        *out=q; return r_angleparm_num;
    }
}

static int collect_torsion_queries(TorsionQuery **out) {
    int idx;
    char tmp[NAMELEN];
    if (iformat<=FILETYPE_MOL2) {
        int n=0;
        for (int ii=0;ii<atomnum;ii++) for (int jj=0;jj<atomnum;jj++) for (int kk=jj+1;kk<atomnum;kk++) for (int ll=0;ll<atomnum;ll++){
            if(ii==kk||ll==jj) continue;
            int ij=0,jk=0,lk=0;
            for(int c=0;c<6;c++) if(atom[ii].con[c]==jj){ij=1;break;}
            for(int c=0;c<6;c++) if(atom[jj].con[c]==kk){jk=1;break;}
            for(int c=0;c<6;c++) if(atom[ll].con[c]==kk){lk=1;break;}
            if(ij&&jk&&lk) n++;
        }
        TorsionQuery *q=malloc(sizeof(TorsionQuery)*n);
        idx=0;
        for (int ii=0;ii<atomnum;ii++) for (int jj=0;jj<atomnum;jj++) for (int kk=jj+1;kk<atomnum;kk++) for (int ll=0;ll<atomnum;ll++){
            if(ii==kk||ll==jj) continue;
            int ij=0,jk=0,lk=0;
            for(int c=0;c<6;c++) if(atom[ii].con[c]==jj){ij=1;break;}
            for(int c=0;c<6;c++) if(atom[jj].con[c]==kk){jk=1;break;}
            for(int c=0;c<6;c++) if(atom[ll].con[c]==kk){lk=1;break;}
            if(!ij||!jk||!lk) continue;
            strcpy(q[idx].name1,atom[ii].ambername);
            strcpy(q[idx].name2,atom[jj].ambername);
            strcpy(q[idx].name3,atom[kk].ambername);
            strcpy(q[idx].name4,atom[ll].ambername);
            if(strcmp(q[idx].name2,q[idx].name3)>0){
                strcpy(tmp,q[idx].name3);strcpy(q[idx].name3,q[idx].name2);strcpy(q[idx].name2,tmp);
                strcpy(tmp,q[idx].name4);strcpy(q[idx].name4,q[idx].name1);strcpy(q[idx].name1,tmp);
            } else if(!strcmp(q[idx].name2,q[idx].name3)&&strcmp(q[idx].name1,q[idx].name4)>0){
                strcpy(tmp,q[idx].name4);strcpy(q[idx].name4,q[idx].name1);strcpy(q[idx].name1,tmp);
            }
            q[idx].attn=1;
            q[idx].atid1=ii; q[idx].atid2=jj; q[idx].atid3=kk; q[idx].atid4=ll;
            idx++;
        }
        *out=q; return idx;
    } else {
        TorsionQuery *q=malloc(sizeof(TorsionQuery)*r_torsionparm_num);
        for (int ii=0;ii<r_torsionparm_num;ii++){
            strcpy(q[ii].name1,r_torsionparm[ii].name1);
            strcpy(q[ii].name2,r_torsionparm[ii].name2);
            strcpy(q[ii].name3,r_torsionparm[ii].name3);
            strcpy(q[ii].name4,r_torsionparm[ii].name4);
            q[ii].attn=r_torsionparm[ii].attn;
            q[ii].atid1=q[ii].atid2=q[ii].atid3=q[ii].atid4=-1;
        }
        *out=q; return r_torsionparm_num;
    }
}

static int collect_improper_queries(ImproperQuery **out) {
    char tmp[NAMELEN];
    if (iformat<=FILETYPE_MOL2) {
        ImproperQuery *q=malloc(sizeof(ImproperQuery)*impropernum);
        for (int ii=0;ii<impropernum;ii++){
            int a1=improper[ii].atid1,a2=improper[ii].atid2,
                a3=improper[ii].atid3,a4=improper[ii].atid4;
            strcpy(q[ii].name1,atom[a1].ambername);
            strcpy(q[ii].name2,atom[a2].ambername);
            strcpy(q[ii].name3,atom[a3].ambername);
            strcpy(q[ii].name4,atom[a4].ambername);
            if(strcmp(q[ii].name1,q[ii].name2)>0){strcpy(tmp,q[ii].name2);strcpy(q[ii].name2,q[ii].name1);strcpy(q[ii].name1,tmp);}
            if(strcmp(q[ii].name1,q[ii].name4)>0){strcpy(tmp,q[ii].name4);strcpy(q[ii].name4,q[ii].name1);strcpy(q[ii].name1,tmp);}
            if(strcmp(q[ii].name2,q[ii].name4)>0){strcpy(tmp,q[ii].name4);strcpy(q[ii].name4,q[ii].name2);strcpy(q[ii].name2,tmp);}
            q[ii].numX=0; q[ii].attn=1;
            q[ii].atid1=a1; q[ii].atid2=a2; q[ii].atid3=a3; q[ii].atid4=a4;
        }
        *out=q; return impropernum;
    } else {
        ImproperQuery *q=malloc(sizeof(ImproperQuery)*r_improperparm_num);
        for (int ii=0;ii<r_improperparm_num;ii++){
            strcpy(q[ii].name1,r_improperparm[ii].name1);
            strcpy(q[ii].name2,r_improperparm[ii].name2);
            strcpy(q[ii].name3,r_improperparm[ii].name3);
            strcpy(q[ii].name4,r_improperparm[ii].name4);
            q[ii].numX=r_improperparm[ii].numX;
            q[ii].attn=r_improperparm[ii].attn;
            q[ii].atid1=q[ii].atid2=q[ii].atid3=q[ii].atid4=-1;
        }
        *out=q; return r_improperparm_num;
    }
}
/*========================================================================
 * rleaplog
 *======================================================================*/
void rleaplog(char *filename, int read_mode) {
/*
Building topology.
Building atom parameters.
For atom (R<NPO | I:1 | 77>.A<OH 10>) could not find vdW (or other) parameters for type (Os)
Building bond parameters.
Could not find bond parameter for atom types: N3 - P
Building angle parameters.
Could not find angle parameter for atom types: OT - P - OH
Building proper torsion parameters.
 ** No torsion terms for atom types: CT-N3-P-O2
Building improper torsion parameters.
 ** Warning: No sp2 improper torsion term for  CA-OV-C-CA
*/
    typedef enum {
        SECT_NONE=0, SECT_TOPO, SECT_ATOM,
        SECT_BOND, SECT_ANGLE, SECT_TORS, SECT_IMPR
    } SECT;

    SECT  rindex = SECT_NONE;
    int   suc;
    char *pos, *tok1, *tok2, *tok3, *tok4;
    char  line[MAXCHAR];
    FILE *fp;

    if ((fp = fopen(filename, "r")) == NULL) {
        fprintf(stdout, "Cannot open %s in rleaplog(), exit\n", filename);
        exit(1);
    }
    if (debug && read_mode)
        printf("\n-- Begin read leaplog file: %s --\n\n", filename);

    while (fgets(line, MAXCHAR, fp)) {

        if      (strstr(line, "Building topology"))         { rindex = SECT_TOPO;  continue; }
        else if (strstr(line, "Building atom parameters"))  { rindex = SECT_ATOM;  continue; }
        else if (strstr(line, "Building bond parameters"))  { rindex = SECT_BOND;  continue; }
        else if (strstr(line, "Building angle parameters")) { rindex = SECT_ANGLE; continue; }
        else if (strstr(line, "Building proper torsion"))   { rindex = SECT_TORS;  continue; }
        else if (strstr(line, "Building improper torsion")) { rindex = SECT_IMPR;  continue; }

        if (rindex == SECT_ATOM) {
            if (strncmp(line, "For atom",8)) continue;
            if (!(pos=strstr(line, "for type"))) continue;
            pos += strlen("for type");
            tok1 = strtok(pos, " :()");

            if (!read_mode) {
                r_atomtype_num++;
                r_vdwparm_num++;
            } else {
                suc = 0;
                for (int i = 0; i < r_atomtype_num; i++)
                    if (!strcmp(r_atomtype[i].name, tok1)) { suc = 1; break; }
                if (!suc) {
                    strcpy(r_atomtype[r_atomtype_num].name, tok1);
                    r_atomtype_num++;
                    strcpy(r_vdwparm[r_vdwparm_num].name, tok1);
                    r_vdwparm_num++;
                }
            }
        }
        else if (rindex == SECT_BOND) {
            if (!strstr(line, "Could not find bond")) continue;
            pos = strchr(line+20,':') + 1;
            if (!(tok1 = strtok(pos,  " -"))) continue;
            if (!(tok2 = strtok(NULL, " -\n\r"))) continue;

            if (!read_mode) {
                r_bondparm_num++;
            } else {
                suc = 0;
                for (int i = 0; i < r_bondparm_num; i++) {
                    if (!strcmp(r_bondparm[i].name1,tok1)&&!strcmp(r_bondparm[i].name2,tok2)){suc=1;break;}
                    if (!strcmp(r_bondparm[i].name1,tok2)&&!strcmp(r_bondparm[i].name2,tok1)){suc=1;break;}
                }
                if (!suc) {
                    strcpy(r_bondparm[r_bondparm_num].name1, tok1);
                    strcpy(r_bondparm[r_bondparm_num].name2, tok2);
                    r_bondparm_num++;
                }
            }
        }
        else if (rindex == SECT_ANGLE) {
            if (!strstr(line, "Could not find angle")) continue;
            pos = strchr(line+20,':') + 1;
            if (!(tok1 = strtok(pos,  " -"))) continue;
            if (!(tok2 = strtok(NULL, " -"))) continue;
            if (!(tok3 = strtok(NULL, " -\n\r"))) continue;

            if (!read_mode) {
                r_angleparm_num++;
            } else {
                suc = 0;
                for (int i = 0; i < r_angleparm_num; i++) {
                    if (!strcmp(r_angleparm[i].name2, tok2)) {
                        if (!strcmp(r_angleparm[i].name1,tok1)&&!strcmp(r_angleparm[i].name3,tok3)){suc=1;break;}
                        if (!strcmp(r_angleparm[i].name1,tok3)&&!strcmp(r_angleparm[i].name3,tok1)){suc=1;break;}
                    }
                }
                if (!suc) {
                    strcpy(r_angleparm[r_angleparm_num].name1, tok1);
                    strcpy(r_angleparm[r_angleparm_num].name2, tok2);
                    strcpy(r_angleparm[r_angleparm_num].name3, tok3);
                    r_angleparm_num++;
                }
            }
        }
        else if (rindex == SECT_TORS) {
            if (!strstr(line, "No torsion terms")) continue;
            pos = strchr(line+20,':') + 1;
            if (!(tok1 = strtok(pos,  " -"))) continue;
            if (!(tok2 = strtok(NULL, " -"))) continue;
            if (!(tok3 = strtok(NULL, " -"))) continue;
            if (!(tok4 = strtok(NULL, " -\n\r"))) continue;

            if (!read_mode) {
                r_torsionparm_num++;
            } else {
                suc = 0;
                for (int i = 0; i < r_torsionparm_num; i++) {
                    if (!strcmp(r_torsionparm[i].name1,tok1)&&!strcmp(r_torsionparm[i].name2,tok2)&&
                        !strcmp(r_torsionparm[i].name3,tok3)&&!strcmp(r_torsionparm[i].name4,tok4)){suc=1;break;}
                    if (!strcmp(r_torsionparm[i].name4,tok1)&&!strcmp(r_torsionparm[i].name3,tok2)&&
                        !strcmp(r_torsionparm[i].name2,tok3)&&!strcmp(r_torsionparm[i].name1,tok4)){suc=1;break;}
                }
                if (!suc) {
                    strcpy(r_torsionparm[r_torsionparm_num].name1, tok1);
                    strcpy(r_torsionparm[r_torsionparm_num].name2, tok2);
                    strcpy(r_torsionparm[r_torsionparm_num].name3, tok3);
                    strcpy(r_torsionparm[r_torsionparm_num].name4, tok4);
                    r_torsionparm_num++;
                }
            }
        }
        else if (rindex == SECT_IMPR) {
            /* ** Warning: No sp2 improper torsion term for  CA-OV-C-CA */
            /* skipped for now */
            continue;
        }
    }
    fclose(fp);

    if (debug && read_mode) {
        printf("\nMASS\n");
        for (int i=0;i<r_atomtype_num;i++)
            printf("\n%-2s %9.4lf %9.4lf", r_atomtype[i].name, r_atomtype[i].mass, r_atomtype[i].pol);
        printf("\n\nBOND\n");
        for (int i=0;i<r_bondparm_num;i++)
            printf("\n%-2s %-2s %9.4lf %9.4lf", r_bondparm[i].name1, r_bondparm[i].name2,
                   r_bondparm[i].force, r_bondparm[i].length);
        printf("\n\nANGLE\n");
        for (int i=0;i<r_angleparm_num;i++)
            printf("\n%-2s %-2s %-2s %9.4lf %9.4lf", r_angleparm[i].name1, r_angleparm[i].name2,
                   r_angleparm[i].name3, r_angleparm[i].force, r_angleparm[i].angle);
        printf("\n\nTORSION\n");
        for (int i=0;i<r_torsionparm_num;i++)
            printf("\n%-2s %-2s %-2s %-2s %9.4lf %9.4lf %9.4lf",
                   r_torsionparm[i].name1, r_torsionparm[i].name2,
                   r_torsionparm[i].name3, r_torsionparm[i].name4,
                   r_torsionparm[i].phase, r_torsionparm[i].force, r_torsionparm[i].fterm);
        printf("\n\nVDW\n");
        for (int i=0;i<r_vdwparm_num;i++)
            printf("\n%-2s %9.4lf %9.4lf", r_vdwparm[i].name,
                   r_vdwparm[i].radius, r_vdwparm[i].pot);
        printf("\n-- Finished reading leaplog file: %s --\n\n", filename);
    }
}

/*========================================================================
 * rleaplog
 *======================================================================*/
void __rleaplog(char *filename, int read_mode) {
/*
Building topology.
Building atom parameters.
For atom (R<NPO | I:1 | 77>.A<OH 10>) could not find vdW (or other) parameters for type (Os)
Building bond parameters.
Could not find bond parameter for atom types: N3 - P
Building angle parameters.
Could not find angle parameter for atom types: OT - P - OH
Building proper torsion parameters.
 ** No torsion terms for atom types: CT-N3-P-O2
Building improper torsion parameters.
 ** Warning: No sp2 improper torsion term for  CA-OV-C-CA
*/
    typedef enum {
        SECT_NONE=0, SECT_TOPO, SECT_ATOM,
        SECT_BOND, SECT_ANGLE, SECT_TORS, SECT_IMPR
    } SECT;

    SECT  rindex = SECT_NONE;
    int   suc;
    char *pos;
    char  tok1[MAXCHAR], tok2[MAXCHAR], tok3[MAXCHAR], tok4[MAXCHAR];
    char line[MAXCHAR];
    FILE *fp;

    if ((fp = fopen(filename, "r")) == NULL) {
        fprintf(stdout, "Cannot open %s in rleaplog(), exit\n", filename);
        exit(1);
    }
    if (debug && read_mode)
        printf("\n-- Begin read leaplog file: %s --\n\n", filename);

    while (fgets(line, MAXCHAR, fp)) {

        if      (strstr(line, "Building topology"))         { rindex = SECT_TOPO;  continue; }
        else if (strstr(line, "Building atom parameters"))  { rindex = SECT_ATOM;  continue; }
        else if (strstr(line, "Building bond parameters"))  { rindex = SECT_BOND;  continue; }
        else if (strstr(line, "Building angle parameters")) { rindex = SECT_ANGLE; continue; }
        else if (strstr(line, "Building proper torsion"))   { rindex = SECT_TORS;  continue; }
        else if (strstr(line, "Building improper torsion")) { rindex = SECT_IMPR;  continue; }

        switch (rindex) {
        case SECT_ATOM:
            /* For atom (R<NPO | I:1 | 77>.A<OH 10>) could not find vdW (or other) parameters for type (Oq) */
            if ((pos = strstr(line, "for type (")) == NULL) continue;
            pos += strlen("for type (");
            copy_type(tok1, pos);

            if (!read_mode) {
                r_atomtype_num++;
                r_vdwparm_num++;
            } else {
                suc = 0;
                for (int i = 0; i < r_atomtype_num; i++)
                    if (!strcmp(r_atomtype[i].name, tok1)) { suc = 1; break; }
                if (!suc) {
                    strcpy(r_atomtype[r_atomtype_num].name, tok1);
                    r_atomtype_num++;
                    strcpy(r_vdwparm[r_vdwparm_num].name, tok1);
                    r_vdwparm_num++;
                }
            }
            break;
        case SECT_BOND:
            /* Could not find bond parameter for atom types: N3 - P */
            if ((pos = strstr(line, "Could not find bond parameter for atom types:")) == NULL)
                continue;
            pos += strlen("Could not find bond parameter for atom types:");
            if (sscanf(pos, "%s - %s", tok1, tok2) != 2) continue;

            if (!read_mode) {
                r_bondparm_num++;
            } else {
                suc = 0;
                for (int i = 0; i < r_bondparm_num; i++) {
                    if (!strcmp(r_bondparm[i].name1,tok1)&&!strcmp(r_bondparm[i].name2,tok2)){suc=1;break;}
                    if (!strcmp(r_bondparm[i].name1,tok2)&&!strcmp(r_bondparm[i].name2,tok1)){suc=1;break;}
                }
                if (!suc) {
                    strcpy(r_bondparm[r_bondparm_num].name1, tok1);
                    strcpy(r_bondparm[r_bondparm_num].name2, tok2);
                    r_bondparm_num++;
                }
            }
            break;
        case SECT_ANGLE:
            /* Could not find angle parameter for atom types: OT - P - OH */
            if ((pos = strstr(line, "Could not find angle parameter for atom types:")) == NULL)
                continue;
            pos += strlen("Could not find angle parameter for atom types:");
            if (sscanf(pos, "%s - %s - %s", tok1, tok2, tok3) != 3) continue;

            if (!read_mode) {
                r_angleparm_num++;
            } else {
                suc = 0;
                for (int i = 0; i < r_angleparm_num; i++) {
                    if (!strcmp(r_angleparm[i].name2, tok2)) {
                        if (!strcmp(r_angleparm[i].name1,tok1)&&!strcmp(r_angleparm[i].name3,tok3)){suc=1;break;}
                        if (!strcmp(r_angleparm[i].name1,tok3)&&!strcmp(r_angleparm[i].name3,tok1)){suc=1;break;}
                    }
                }
                if (!suc) {
                    strcpy(r_angleparm[r_angleparm_num].name1, tok1);
                    strcpy(r_angleparm[r_angleparm_num].name2, tok2);
                    strcpy(r_angleparm[r_angleparm_num].name3, tok3);
                    r_angleparm_num++;
                }
            }
            break;
        case SECT_TORS:
            /* ** No torsion terms for atom types: CT-N3-P-O2 */
            if ((pos = strstr(line, "No torsion terms for atom types:")) == NULL)
                continue;
            pos += strlen("No torsion terms for atom types:");
            while (*pos == ' ') pos++;
            if (sscanf(pos, "%[^-]-%[^-]-%[^-]-%s", tok1, tok2, tok3, tok4) != 4) continue;
            char *end = tok4 + strlen(tok4) - 1;
            while (end > tok4 && isspace((unsigned char)*end)) *end-- = '\0';

            if (!read_mode) {
                r_torsionparm_num++;
            } else {
                suc = 0;
                for (int i = 0; i < r_torsionparm_num; i++) {
                    if (!strcmp(r_torsionparm[i].name1,tok1)&&!strcmp(r_torsionparm[i].name2,tok2)&&
                        !strcmp(r_torsionparm[i].name3,tok3)&&!strcmp(r_torsionparm[i].name4,tok4)){suc=1;break;}
                    if (!strcmp(r_torsionparm[i].name4,tok1)&&!strcmp(r_torsionparm[i].name3,tok2)&&
                        !strcmp(r_torsionparm[i].name2,tok3)&&!strcmp(r_torsionparm[i].name1,tok4)){suc=1;break;}
                }
                if (!suc) {
                    strcpy(r_torsionparm[r_torsionparm_num].name1, tok1);
                    strcpy(r_torsionparm[r_torsionparm_num].name2, tok2);
                    strcpy(r_torsionparm[r_torsionparm_num].name3, tok3);
                    strcpy(r_torsionparm[r_torsionparm_num].name4, tok4);
                    r_torsionparm_num++;
                }
            }
            break;
        case SECT_IMPR:
            /* skipped for now */
            break;
        default:
            break;
        }
    }
    fclose(fp);

    if (debug && read_mode) {
        printf("\nMASS\n");
        for (int i=0;i<r_atomtype_num;i++)
            printf("\n%-2s %9.4lf %9.4lf", r_atomtype[i].name, r_atomtype[i].mass, r_atomtype[i].pol);
        printf("\n\nBOND\n");
        for (int i=0;i<r_bondparm_num;i++)
            printf("\n%-2s %-2s %9.4lf %9.4lf", r_bondparm[i].name1, r_bondparm[i].name2,
                   r_bondparm[i].force, r_bondparm[i].length);
        printf("\n\nANGLE\n");
        for (int i=0;i<r_angleparm_num;i++)
            printf("\n%-2s %-2s %-2s %9.4lf %9.4lf", r_angleparm[i].name1, r_angleparm[i].name2,
                   r_angleparm[i].name3, r_angleparm[i].force, r_angleparm[i].angle);
        printf("\n\nTORSION\n");
        for (int i=0;i<r_torsionparm_num;i++)
            printf("\n%-2s %-2s %-2s %-2s %9.4lf %9.4lf %9.4lf",
                   r_torsionparm[i].name1, r_torsionparm[i].name2,
                   r_torsionparm[i].name3, r_torsionparm[i].name4,
                   r_torsionparm[i].phase, r_torsionparm[i].force, r_torsionparm[i].fterm);
        printf("\n\nVDW\n");
        for (int i=0;i<r_vdwparm_num;i++)
            printf("\n%-2s %9.4lf %9.4lf", r_vdwparm[i].name,
                   r_vdwparm[i].radius, r_vdwparm[i].pot);
        printf("\n-- Finished reading leaplog file: %s --\n\n", filename);
    }
}

/*========================================================================
 * frcmod_parse: unified parser for rfrcmod (r_* arrays) and readfrcmod
 * (main arrays).  Three modes:
 *   FRCMOD_COUNT    — first pass, increment r_*_num counters only
 *   FRCMOD_READ_R   — second pass, fill r_* arrays
 *   FRCMOD_READ_MAIN— fill main atomtype/bondparm/... arrays
 *======================================================================*/
static void frcmod_parse(const char *filename, FrcmodMode mode) {
    FILE *fpin;
    int   mindex=0, bindex=0, aindex=0, tindex=0, iindex=0, vindex=0;
    int   pos_tor=0, tmpnum=0;
    char  lline[MAXCHAR];
    char  type[MAXCHAR], type1[MAXCHAR], type2[MAXCHAR];
    char  type3[MAXCHAR], type4[MAXCHAR], tmpchar[MAXCHAR];

    if ((fpin = fopen(filename, "r")) == NULL) {
        fprintf(stdout, "Cannot open %s in frcmod_parse(), exit\n", filename);
        exit(1);
    }

    while (fgets(lline, MAXCHAR, fpin)) {

        if      (!strncmp(lline,"MASS",4))     { mindex=1; continue; }
        else if (!strncmp(lline,"BOND",4))     { bindex=1; continue; }
        else if (!strncmp(lline,"ANGL",4))     { aindex=1; continue; }
        else if (!strncmp(lline,"DIHE",4))     {
            if (mode==FRCMOD_READ_MAIN) pos_tor = ftell(fpin);
            tindex=1; continue;
        }
        else if (!strncmp(lline,"IMPROPER",8)) { iindex=1; continue; }
        else if (!strncmp(lline,"NONBON",6))   { vindex=1; continue; }

        if (strlen(lline) <= 2) {
            if (mindex==1) mindex=0;
            if (bindex==1) bindex=0;
            if (aindex==1) aindex=0;
            if (tindex==1) {
                if (mode==FRCMOD_READ_MAIN) { fseek(fpin,pos_tor,0); tindex=2; }
                else tindex=0;
                continue;
            }
            if (tindex==2) tindex=0;
            if (iindex==1) iindex=0;
            if (vindex==1) vindex=0;
            continue;
        }

        /* MASS */
        if (mindex==1) {
            if (mode==FRCMOD_COUNT) {
                r_atomtype_num++;
            } else if (mode==FRCMOD_READ_R) {
                r_atomtype[r_atomtype_num].attn = 0;
                sscanf(lline, "%s%lf%lf", r_atomtype[r_atomtype_num].name,
                       &r_atomtype[r_atomtype_num].mass, &r_atomtype[r_atomtype_num].pol);
                if (strstr(lline,"ATTN")) r_atomtype[r_atomtype_num].attn=1;
                r_atomtype_num++;
            } else {
                sscanf(lline, "%s%lf%lf", atomtype[atomtypenum].name,
                       &atomtype[atomtypenum].mass, &atomtype[atomtypenum].pol);
                atomtypenum++;
                GROW_ARRAY(atomtype,atomtypenum,maxatomtype,MAX_FF_ATOMTYPE,ATOMTYPE);
            }
        }

        /* BOND */
        if (bindex==1) {
            if (mode==FRCMOD_COUNT) {
                r_bondparm_num++;
            } else {
                BOND_FF *bp = (mode==FRCMOD_READ_R) ? &r_bondparm[r_bondparm_num]
                                                     : &bondparm[bondparmnum];
                bp->attn = 0;
                copy_type(bp->name1, lline);
                copy_type(bp->name2, lline+3);
                if (strcmp(bp->name1,bp->name2)>0) {
                    strcpy(tmpchar,bp->name1); strcpy(bp->name1,bp->name2); strcpy(bp->name2,tmpchar);
                }
                sscanf(&lline[5], "%lf%lf", &bp->force, &bp->length);
                if (strstr(lline,"ATTN")) bp->attn=1;
                if (mode==FRCMOD_READ_R) r_bondparm_num++;
                else {
                    bondparmnum++;
                    GROW_ARRAY(bondparm,bondparmnum,maxbondparm,MAX_FF_BOND,BOND_FF);
                }
            }
        }

        /* ANGLE */
        if (aindex==1) {
            if (mode==FRCMOD_COUNT) {
                r_angleparm_num++;
            } else {
                ANGLE *ap = (mode==FRCMOD_READ_R) ? &r_angleparm[r_angleparm_num]
                                                   : &angleparm[angleparmnum];
                ap->attn = 0;
                copy_type(ap->name1, lline);
                copy_type(ap->name2, lline+3);
                copy_type(ap->name3, lline+6);
                if (strcmp(ap->name1,ap->name3)>0) {
                    strcpy(tmpchar,ap->name1); strcpy(ap->name1,ap->name3); strcpy(ap->name3,tmpchar);
                }
                sscanf(&lline[8], "%lf%lf", &ap->force, &ap->angle);
                if (strstr(lline,"ATTN")) ap->attn=1;
                if (mode==FRCMOD_READ_R) r_angleparm_num++;
                else {
                    angleparmnum++;
                    GROW_ARRAY(angleparm,angleparmnum,maxangleparm,MAX_FF_ANGLE,ANGLE);
                }
            }
        }

        /* DIHE */
        if (tindex==1 || tindex==2) {
            int hasX = (lline[0]=='X'||lline[3]=='X'||lline[6]=='X'||lline[9]=='X');
            if (tindex==1 &&  hasX) continue;
            if (tindex==2 && !hasX) continue;

            if (mode==FRCMOD_COUNT) {
                r_torsionparm_num++;
            } else {
                TORSION *tp = (mode==FRCMOD_READ_R) ? &r_torsionparm[r_torsionparm_num]
                                                     : &torsionparm[torsionparmnum];
                tp->attn = 0;
                copy_type(tp->name1, lline);
                copy_type(tp->name2, lline+3);
                copy_type(tp->name3, lline+6);
                copy_type(tp->name4, lline+9);
                if (strcmp(tp->name2,tp->name3)>0) {
                    strcpy(tmpchar,tp->name2); strcpy(tp->name2,tp->name3); strcpy(tp->name3,tmpchar);
                    strcpy(tmpchar,tp->name1); strcpy(tp->name1,tp->name4); strcpy(tp->name4,tmpchar);
                } else if (!strcmp(tp->name2,tp->name3) && strcmp(tp->name1,tp->name4)>0) {
                    strcpy(tmpchar,tp->name1); strcpy(tp->name1,tp->name4); strcpy(tp->name4,tmpchar);
                }
                sscanf(&lline[11], "%d%lf%lf%lf",
                       &tp->mul, &tp->force, &tp->phase, &tp->fterm);
                if (strstr(lline,"ATTN")) tp->attn=1;
                if (mode==FRCMOD_READ_R) r_torsionparm_num++;
                else {
                    torsionparmnum++;
                    GROW_ARRAY(torsionparm,torsionparmnum,maxtorsionparm,MAX_FF_TORSION,TORSION);
                }
            }
        }

        /* IMPROPER */
        if (iindex==1) {
            if (mode==FRCMOD_COUNT) {
                r_improperparm_num++;
            } else {
                IMPROPER *ip = (mode==FRCMOD_READ_R) ? &r_improperparm[r_improperparm_num]
                                                      : &improperparm[improperparmnum];
                ip->attn = 0;
                tmpnum = 0;
                copy_type(type1, lline);
                copy_type(type2, lline+3);
                copy_type(type3, lline+6);
                copy_type(type4, lline+9);
                if (lline[0]=='X') tmpnum++;
                if (lline[3]=='X') tmpnum++;
                if (lline[6]=='X') tmpnum++;
                if (lline[9]=='X') tmpnum++;
                if (tmpnum==0) {
                    if(strcmp(type1,type2)>0){strcpy(type,type2);strcpy(type2,type1);strcpy(type1,type);}
                    if(strcmp(type1,type4)>0){strcpy(type,type4);strcpy(type4,type1);strcpy(type1,type);}
                    if(strcmp(type2,type4)>0){strcpy(type,type4);strcpy(type4,type2);strcpy(type2,type);}
                }
                if (tmpnum==1 || (tmpnum==2 && lline[6]=='X')) {
                    if(lline[0]=='X'&&strcmp(type2,type4)>0){strcpy(type,type4);strcpy(type4,type2);strcpy(type2,type);}
                    if(lline[3]=='X'&&strcmp(type1,type4)>0){strcpy(type,type4);strcpy(type4,type1);strcpy(type1,type);}
                    if(lline[9]=='X'&&strcmp(type1,type2)>0){strcpy(type,type1);strcpy(type1,type2);strcpy(type2,type);}
                }
                strcpy(ip->name1,type1); strcpy(ip->name2,type2);
                strcpy(ip->name3,type3); strcpy(ip->name4,type4);
                sscanf(&lline[11], "%lf%lf%lf", &ip->force, &ip->phase, &ip->fterm);
                ip->mul  = 1;
                ip->numX = tmpnum;
                if (strstr(lline,"ATTN")) ip->attn=1;
                if (mode==FRCMOD_READ_R) r_improperparm_num++;
                else {
                    improperparmnum++;
                    GROW_ARRAY(improperparm,improperparmnum,maximproperparm,MAX_FF_IMPROPER,IMPROPER);
                }
            }
        }

        /* NONBON */
        if (vindex==1) {
            if (mode==FRCMOD_COUNT) {
                r_vdwparm_num++;
            } else if (mode==FRCMOD_READ_R) {
                r_vdwparm[r_vdwparm_num].attn = 0;
                sscanf(lline, "%s%lf%lf", r_vdwparm[r_vdwparm_num].name,
                       &r_vdwparm[r_vdwparm_num].radius, &r_vdwparm[r_vdwparm_num].pot);
                if (strstr(lline,"ATTN")) r_vdwparm[r_vdwparm_num].attn=1;
                r_vdwparm_num++;
            } else {
                sscanf(lline, "%s%lf%lf", vdwparm[vdwparmnum].name,
                       &vdwparm[vdwparmnum].radius, &vdwparm[vdwparmnum].pot);
                vdwparmnum++;
                GROW_ARRAY(vdwparm,vdwparmnum,maxvdwparm,MAX_FF_VDW,VDW);
            }
        }
    }
    fclose(fpin);

    if (debug && mode==FRCMOD_READ_R) {
        printf("\n-- frcmod r_* arrays after reading %s --\n", filename);
        printf("\nMASS\n");
        for (int i=0;i<r_atomtype_num;i++)
            printf("\n%-2s %9.4lf %9.4lf %5d", r_atomtype[i].name,
                   r_atomtype[i].mass, r_atomtype[i].pol, r_atomtype[i].attn);
        printf("\n\nBOND\n");
        for (int i=0;i<r_bondparm_num;i++)
            printf("\n%-2s %-2s %9.4lf %9.4lf %5d", r_bondparm[i].name1, r_bondparm[i].name2,
                   r_bondparm[i].force, r_bondparm[i].length, r_bondparm[i].attn);
        printf("\n\nANGLE\n");
        for (int i=0;i<r_angleparm_num;i++)
            printf("\n%-2s %-2s %-2s %9.4lf %9.4lf %5d", r_angleparm[i].name1,
                   r_angleparm[i].name2, r_angleparm[i].name3,
                   r_angleparm[i].force, r_angleparm[i].angle, r_angleparm[i].attn);
        printf("\n\nTORSION\n");
        for (int i=0;i<r_torsionparm_num;i++)
            printf("\n%-2s %-2s %-2s %-2s %9.4lf %9.4lf %9.4lf %5d",
                   r_torsionparm[i].name1, r_torsionparm[i].name2,
                   r_torsionparm[i].name3, r_torsionparm[i].name4,
                   r_torsionparm[i].phase, r_torsionparm[i].force,
                   r_torsionparm[i].fterm, r_torsionparm[i].attn);
        printf("\n\nVDW\n");
        for (int i=0;i<r_vdwparm_num;i++)
            printf("\n%-2s %9.4lf %9.4lf %5d", r_vdwparm[i].name,
                   r_vdwparm[i].radius, r_vdwparm[i].pot, r_vdwparm[i].attn);
        printf("\n-- done --\n\n");
    }
    if (debug && mode==FRCMOD_READ_MAIN) {
        printf("\n-- main arrays after reading frcmod %s --\n", filename);
        printf("\nMASS\n");
        for (int i=last_atomtypenum;i<atomtypenum;i++)
            printf("\n%-2s %9.4lf %9.4lf", atomtype[i].name, atomtype[i].mass, atomtype[i].pol);
        printf("\n\nBOND\n");
        for (int i=last_bondparmnum;i<bondparmnum;i++)
            printf("\n%-2s %-2s %9.4lf %9.4lf", bondparm[i].name1, bondparm[i].name2,
                   bondparm[i].force, bondparm[i].length);
        printf("\n\nANGLE\n");
        for (int i=last_angleparmnum;i<angleparmnum;i++)
            printf("\n%-2s %-2s %-2s %9.4lf %9.4lf", angleparm[i].name1,
                   angleparm[i].name2, angleparm[i].name3,
                   angleparm[i].force, angleparm[i].angle);
        printf("\n\nTORSION\n");
        for (int i=last_torsionparmnum;i<torsionparmnum;i++)
            printf("\n%-2s %-2s %-2s %-2s %9.4lf %9.4lf %9.4lf",
                   torsionparm[i].name1, torsionparm[i].name2,
                   torsionparm[i].name3, torsionparm[i].name4,
                   torsionparm[i].phase, torsionparm[i].force, torsionparm[i].fterm);
        printf("\n\nIMPROPER\n");
        for (int i=last_improperparmnum;i<improperparmnum;i++)
            printf("\n%-2s %-2s %-2s %-2s %9.4lf %9.4lf %9.4lf",
                   improperparm[i].name1, improperparm[i].name2,
                   improperparm[i].name3, improperparm[i].name4,
                   improperparm[i].phase, improperparm[i].force, improperparm[i].fterm);
        printf("\n\nVDW\n");
        for (int i=last_vdwparmnum;i<vdwparmnum;i++)
            printf("\n%-2s %9.4lf %9.4lf", vdwparm[i].name,
                   vdwparm[i].radius, vdwparm[i].pot);
        printf("\n-- done --\n\n");
    }
}

void rfrcmod(char *filename, int read_mode) {
    frcmod_parse(filename, read_mode ? FRCMOD_READ_R : FRCMOD_COUNT);
}
void readfrcmod(char *filename) {
    frcmod_parse(filename, FRCMOD_READ_MAIN);
}

/*========================================================================
 * readparm
 *======================================================================*/
void readparm(char *filename) {
    int   mindex=-1, bindex=0, aindex=0, tindex=0, iindex=0, vindex=0;
    int   num=0, tmpnum, vdwparmnum_old;
    long  pos_tor=0;
    FILE *fpin;
    char  lline[MAXCHAR];
    char  type[MAXCHAR], type1[MAXCHAR], type2[MAXCHAR];
    char  type3[MAXCHAR], type4[MAXCHAR];

    if ((fpin = fopen(filename, "r")) == NULL) {
        fprintf(stdout, "Cannot open %s in readparm(), exit\n", filename);
        exit(1);
    }
    while (fgets(lline, MAXCHAR, fpin)) {
        num++;
        if (mindex==-1 && num==1)             { mindex=1; continue; }
        if (mindex==1  && spaceline(lline)==1) { mindex=0; num=0; bindex=-1; continue; }
        if (bindex==-1 && num==1)             { bindex=1; continue; }
        if (bindex==1  && spaceline(lline)==1) { bindex=0; aindex=1; continue; }
        if (aindex==1  && spaceline(lline)==1) { aindex=0; tindex=1; pos_tor=ftell(fpin); continue; }
        if (tindex==1  && spaceline(lline)==1) { tindex=2; fseek(fpin,pos_tor,0); continue; }
        if (tindex==2  && spaceline(lline)==1) { tindex=0; iindex=1; continue; }
        if (iindex==1  && spaceline(lline)==1) { iindex=0; vindex=-1; num=0; continue; }
        if (vindex==-1 && num==2)             { vindex=1; continue; }
        if (vindex==2  && spaceline(lline)==1) { vindex=0; continue; }
        if (vindex==1  && !strncmp(lline,"MOD4",4)) { vindex=2; continue; }

        if (mindex==1) {
            sscanf(lline, "%s%lf%lf", atomtype[atomtypenum].name,
                   &atomtype[atomtypenum].mass, &atomtype[atomtypenum].pol);
            atomtypenum++;
            GROW_ARRAY(atomtype,atomtypenum,maxatomtype,MAX_FF_ATOMTYPE,ATOMTYPE);
        }
        if (bindex==1) {
            copy_type(bondparm[bondparmnum].name1, lline);
            copy_type(bondparm[bondparmnum].name2, lline+3);
            sscanf(&lline[5], "%lf%lf", &bondparm[bondparmnum].force, &bondparm[bondparmnum].length);
            bondparmnum++;
            GROW_ARRAY(bondparm,bondparmnum,maxbondparm,MAX_FF_BOND,BOND_FF);
        }
        if (aindex==1) {
            copy_type(angleparm[angleparmnum].name1, lline);
            copy_type(angleparm[angleparmnum].name2, lline+3);
            copy_type(angleparm[angleparmnum].name3, lline+6);
            sscanf(&lline[8], "%lf%lf", &angleparm[angleparmnum].force, &angleparm[angleparmnum].angle);
            angleparmnum++;
            GROW_ARRAY(angleparm,angleparmnum,maxangleparm,MAX_FF_ANGLE,ANGLE);
        }
        if (tindex==1 || tindex==2) {
            int hasX = (lline[0]=='X'||lline[3]=='X'||lline[6]=='X'||lline[9]=='X');
            if (tindex==1 &&  hasX) continue;
            if (tindex==2 && !hasX) continue;
            if (lline[0]==' ') continue;
            copy_type(torsionparm[torsionparmnum].name1, lline);
            copy_type(torsionparm[torsionparmnum].name2, lline+3);
            copy_type(torsionparm[torsionparmnum].name3, lline+6);
            copy_type(torsionparm[torsionparmnum].name4, lline+9);
            sscanf(&lline[11], "%d%lf%lf%lf",
                   &torsionparm[torsionparmnum].mul,  &torsionparm[torsionparmnum].force,
                   &torsionparm[torsionparmnum].phase, &torsionparm[torsionparmnum].fterm);
            torsionparmnum++;
            GROW_ARRAY(torsionparm,torsionparmnum,maxtorsionparm,MAX_FF_TORSION,TORSION);
        }
        if (iindex==1) {
            tmpnum=0;
            copy_type(type1,lline); copy_type(type2,lline+3);
            copy_type(type3,lline+6); copy_type(type4,lline+9);
            if (lline[0]=='X') tmpnum++;
            if (lline[3]=='X') tmpnum++;
            if (lline[6]=='X') tmpnum++;
            if (lline[9]=='X') tmpnum++;
            if (tmpnum==0) {
                if(strcmp(type1,type2)>0){strcpy(type,type2);strcpy(type2,type1);strcpy(type1,type);}
                if(strcmp(type1,type4)>0){strcpy(type,type4);strcpy(type4,type1);strcpy(type1,type);}
                if(strcmp(type2,type4)>0){strcpy(type,type4);strcpy(type4,type2);strcpy(type2,type);}
            }
            if (tmpnum==1 || (tmpnum==2 && lline[6]=='X')) {
                if(lline[0]=='X'&&strcmp(type2,type4)>0){strcpy(type,type4);strcpy(type4,type2);strcpy(type2,type);}
                if(lline[3]=='X'&&strcmp(type1,type4)>0){strcpy(type,type4);strcpy(type4,type1);strcpy(type1,type);}
                if(lline[9]=='X'&&strcmp(type1,type2)>0){strcpy(type,type1);strcpy(type1,type2);strcpy(type2,type);}
            }
            strcpy(improperparm[improperparmnum].name1,type1);
            strcpy(improperparm[improperparmnum].name2,type2);
            strcpy(improperparm[improperparmnum].name3,type3);
            strcpy(improperparm[improperparmnum].name4,type4);
            sscanf(&lline[11], "%lf%lf%lf",
                   &improperparm[improperparmnum].force,
                   &improperparm[improperparmnum].phase,
                   &improperparm[improperparmnum].fterm);
            improperparm[improperparmnum].mul  = 1;
            improperparm[improperparmnum].numX = tmpnum;
            improperparmnum++;
            GROW_ARRAY(improperparm,improperparmnum,maximproperparm,MAX_FF_IMPROPER,IMPROPER);
        }
        if (vindex==1) {
            if (spaceline(lline)==1) continue;
            equ_vdw[equ_vdw_num].num = 0;
            char *tok = strtok(lline, " \t\n\r");
            while (tok) {
                if (equ_vdw[equ_vdw_num].num >= MAX_EQU_VDW_NUM) {
                    printf("\nError: equivalent vdw atoms exceed MAX_EQU_VDW_NUM\n"); exit(1);
                }
                strcpy(equ_vdw[equ_vdw_num].name[equ_vdw[equ_vdw_num].num++], tok);
                tok = strtok(NULL, " \t\n\r");
            }
            if (equ_vdw[equ_vdw_num].num >= 2) equ_vdw_num++;
            if (equ_vdw_num >= MAX_EQU_VDW) {
                printf("\nError: equivalent vdw parameters exceed MAX_EQU_VDW\n"); exit(1);
            }
        }
        if (vindex==2) {
            sscanf(lline, "%s%lf%lf", vdwparm[vdwparmnum].name,
                   &vdwparm[vdwparmnum].radius, &vdwparm[vdwparmnum].pot);
            vdwparmnum++;
            GROW_ARRAY(vdwparm,vdwparmnum,maxvdwparm,MAX_FF_VDW,VDW);
        }
    }
    fclose(fpin);

    if (equ_vdw_num > 0) {
        vdwparmnum_old = vdwparmnum;
        for (int i=0; i<equ_vdw_num; i++)
            for (int j=1; j<equ_vdw[i].num; j++) {
                if (strlen(equ_vdw[i].name[j]) < 1) continue;
                for (int k=0; k<vdwparmnum_old; k++)
                    if (!strcmp(vdwparm[k].name, equ_vdw[i].name[0])) {
                        strcpy(vdwparm[vdwparmnum].name, equ_vdw[i].name[j]);
                        vdwparm[vdwparmnum].radius = vdwparm[k].radius;
                        vdwparm[vdwparmnum].pot    = vdwparm[k].pot;
                        vdwparmnum++;
                        GROW_ARRAY(vdwparm,vdwparmnum,maxvdwparm,MAX_FF_VDW,VDW);
                        break;
                    }
            }
    }

    if (debug) {
        printf("\n-- main arrays after reading parm file %s --\n", filename);
        printf("\nMASS\n");
        for (int i=last_atomtypenum;i<atomtypenum;i++)
            printf("\n%-2s %9.4lf %9.4lf", atomtype[i].name, atomtype[i].mass, atomtype[i].pol);
        printf("\n\nBOND\n");
        for (int i=last_bondparmnum;i<bondparmnum;i++)
            printf("\n%-2s %-2s %9.4lf %9.4lf", bondparm[i].name1, bondparm[i].name2,
                   bondparm[i].force, bondparm[i].length);
        printf("\n\nANGLE\n");
        for (int i=last_angleparmnum;i<angleparmnum;i++)
            printf("\n%-2s %-2s %-2s %9.4lf %9.4lf", angleparm[i].name1,
                   angleparm[i].name2, angleparm[i].name3,
                   angleparm[i].force, angleparm[i].angle);
        printf("\n\nTORSION\n");
        for (int i=last_torsionparmnum;i<torsionparmnum;i++)
            printf("\n%-2s %-2s %-2s %-2s %9.4lf %9.4lf %9.4lf",
                   torsionparm[i].name1, torsionparm[i].name2,
                   torsionparm[i].name3, torsionparm[i].name4,
                   torsionparm[i].phase, torsionparm[i].force, torsionparm[i].fterm);
        printf("\n\nIMPROPER\n");
        for (int i=last_improperparmnum;i<improperparmnum;i++)
            printf("\n%-2s %-2s %-2s %-2s %9.4lf %9.4lf %9.4lf",
                   improperparm[i].name1, improperparm[i].name2,
                   improperparm[i].name3, improperparm[i].name4,
                   improperparm[i].phase, improperparm[i].force, improperparm[i].fterm);
        printf("\n\nVDW\n");
        for (int i=last_vdwparmnum;i<vdwparmnum;i++)
            printf("\n%-2s %9.4lf %9.4lf", vdwparm[i].name,
                   vdwparm[i].radius, vdwparm[i].pot);
        printf("\n-- done --\n\n");
    }
}

/*========================================================================
 * read_parmchk_parm
 *======================================================================*/
static void read_parmchk_parm(const char *filename, int parse_weights) {
    FILE  *fpin;
    char   lline[MAXCHAR];
    char   at[NAMELEN], equaname[NAMELEN], corrname[NAMELEN];
    int    group, improper_flag, ncorr, nequa, equtype, atomicnum;
    double bl,blf,ba,baf,cba,cbaf,tor,ctor,ps,mass;
    int    loop_start = parmnum;

    if ((fpin = fopen(filename, "r")) == NULL) {
        fprintf(stdout, "Cannot open %s in read_parmchk_parm_impl(), exit\n", filename);
        exit(1);
    }
    while (fgets(lline, MAXCHAR, fpin)) {
        if (!strncmp("PARM", lline, 4)) {
            equtype = 0;
            sscanf(&lline[4], "%s%d%d%lf%d%d",
                   at, &improper_flag, &group, &mass, &equtype, &atomicnum);
            strcpy(parm[parmnum].atomtype,  at);
            parm[parmnum].improper  = improper_flag;
            parm[parmnum].group     = group;
            parm[parmnum].mass      = mass;
            parm[parmnum].equtype   = equtype;
            parm[parmnum].atomicnum = atomicnum;
            strcpy(parm[parmnum].equa[0].atomtype, at);
            init_corr(&parm[parmnum].corr[0], at, 0, parmnum);
            nequa = ncorr = 1;
            parm[parmnum].nequa = parm[parmnum].ncorr = 1;
            parmnum++;
            GROW_ARRAY(parm,parmnum,maxparmnum,MAXPARM,PARM);
        }
        if (!strncmp("EQUA", lline, 4)) {
            sscanf(&lline[4], "%s", equaname);
            strcpy(parm[parmnum-1].equa[nequa].atomtype, equaname);
            parm[parmnum-1].nequa++;
            nequa++;
            if (parm[parmnum-1].nequa > MAXEQUA) {
                fprintf(stdout,"Too many equal atom types for %s\n",parm[parmnum-1].atomtype);
                exit(1);
            }
            init_corr(&parm[parmnum-1].corr[ncorr], equaname, 1, -1);
            parm[parmnum-1].ncorr++;
            ncorr++;
            if (parm[parmnum-1].ncorr > MAXCORR) {
                fprintf(stdout,"Too many corresponding atom types for %s\n",parm[parmnum-1].atomtype);
                exit(1);
            }
        }
        if (!strncmp("CORR", lline, 4)) {
            bl=blf=ba=baf=cba=cbaf=tor=ctor=ps=0.0;
            sscanf(&lline[4], "%s%lf%lf%lf%lf%lf%lf%lf%lf%lf",
                   corrname,&bl,&blf,&cba,&cbaf,&ba,&baf,&ctor,&tor,&ps);
            CORR *c = &parm[parmnum-1].corr[ncorr];
            strcpy(c->atomtype, corrname);
            c->bl=bl; c->blf=blf; c->ba=ba; c->baf=baf;
            c->cba=cba; c->cbaf=cbaf; c->tor=tor; c->ctor=ctor;
            c->improper=ps; c->ps=ps; c->pid=-1; c->type=2;
            parm[parmnum-1].ncorr++;
            ncorr++;
            if (parm[parmnum-1].ncorr > MAXCORR) {
                fprintf(stdout,"Too many corresponding atom types for %s\n",parm[parmnum-1].atomtype);
                exit(1);
            }
        }
        if (parse_weights) {
            if (!strncmp("WEIGHT_BL",      lline, 9))  sscanf(&lline[9],  "%lf",&wt.BL);
            if (!strncmp("WEIGHT_BLF",     lline,10))  sscanf(&lline[10], "%lf",&wt.BLF);
            if (!strncmp("WEIGHT_BA",      lline, 9))  sscanf(&lline[9],  "%lf",&wt.BA);
            if (!strncmp("WEIGHT_BAF",     lline,10))  sscanf(&lline[10], "%lf",&wt.BAF);
            if (!strncmp("WEIGHT_X",       lline, 8))  sscanf(&lline[8],  "%lf",&wt.X);
            if (!strncmp("WEIGHT_X3",      lline, 9))  sscanf(&lline[9],  "%lf",&wt.X3);
            if (!strncmp("WEIGHT_BA_CTR",  lline,13))  sscanf(&lline[13], "%lf",&wt.BA_CTR);
            if (!strncmp("WEIGHT_TOR_CTR", lline,14))  sscanf(&lline[14], "%lf",&wt.TOR_CTR);
            if (!strncmp("WEIGHT_IMPROPER",lline,15))  sscanf(&lline[15], "%lf",&wt.IMPROPER);
            if (!strncmp("WEIGHT_GROUP",   lline,12))  sscanf(&lline[12], "%lf",&wt.GROUP);
            if (!strncmp("WEIGHT_EQUTYPE", lline,14))  sscanf(&lline[14], "%lf",&wt.EQUTYPE);
            if (!strncmp("THRESHOLD_BA",   lline,12))  sscanf(&lline[12], "%lf",&THRESHOLD_BA);
            if (!strncmp("DEFAULT_BL",     lline,10))  sscanf(&lline[10], "%lf",&dv.BL);
            if (!strncmp("DEFAULT_BLF",    lline,11))  sscanf(&lline[11], "%lf",&dv.BLF);
            if (!strncmp("DEFAULT_BA",     lline,10))  sscanf(&lline[10], "%lf",&dv.BA);
            if (!strncmp("DEFAULT_BAF",    lline,11))  sscanf(&lline[11], "%lf",&dv.BAF);
            if (!strncmp("DEFAULT_BA_CTR", lline,14))  sscanf(&lline[14], "%lf",&dv.BA_CTR);
            if (!strncmp("DEFAULT_BAF_CTR",lline,15))  sscanf(&lline[15], "%lf",&dv.BAF_CTR);
            if (!strncmp("DEFAULT_TOR",    lline,11))  sscanf(&lline[11], "%lf",&dv.TOR);
            if (!strncmp("DEFAULT_TOR_CTR",lline,15))  sscanf(&lline[15], "%lf",&dv.TOR_CTR);
            if (!strncmp("DEFAULT_FRACT1", lline,14))  sscanf(&lline[14], "%lf",&dv.FRACT1);
            if (!strncmp("DEFAULT_FRACT2", lline,14))  sscanf(&lline[14], "%lf",&dv.FRACT2);
        }
    }
    fclose(fpin);

    for (int i=loop_start; i<parmnum; i++) {
        for (int j=0; j<parm[i].ncorr; j++) {
            CORR *c = &parm[i].corr[j];
            if (c->bl   < 0) c->bl   = dv.BL;
            if (c->blf  < 0) c->blf  = dv.BLF;
            if (c->ba   < 0) c->ba   = dv.BA;
            if (c->baf  < 0) c->baf  = dv.BAF;
            if (c->cba  < 0) c->cba  = dv.BA_CTR;
            if (c->cbaf < 0) c->cbaf = dv.BAF_CTR;
            if (c->tor  < 0) c->tor  = dv.TOR;
            if (c->ctor < 0) c->ctor = dv.TOR_CTR;
            c->ctor = c->ctor * dv.FRACT1 + c->ps * dv.FRACT2;
        }
        for (int j=0; j<parm[i].nequa; j++)
            for (int k=0; k<parmnum; k++)
                if (!strcmp(parm[i].equa[j].atomtype, parm[k].atomtype))
                    { parm[i].equa[j].pid = k; break; }
        for (int j=0; j<parm[i].ncorr; j++)
            for (int k=0; k<parmnum; k++)
                if (!strcmp(parm[i].corr[j].atomtype, parm[k].atomtype))
                    { parm[i].corr[j].pid = k; break; }
    }
}
/*========================================================================
 * cleanup_frcmod
 *======================================================================*/
void cleanup_frcmod(char *filename) {
    FILE *fp1, *fp2;
    char  command[MAXCHAR]="cp -f ";
    int   status;
    typedef struct { char name[NAMELEN]; }                               ATOMTYPE_CL;
    typedef struct { char name1[NAMELEN], name2[NAMELEN]; }              BOND_CL;
    typedef struct { char name1[NAMELEN], name2[NAMELEN], name3[NAMELEN]; } ANGLE_CL;
    typedef struct { char name1[NAMELEN],name2[NAMELEN],name3[NAMELEN],name4[NAMELEN]; double fterm; } TORSION_CL;

    ATOMTYPE_CL at[MAX_FF_ATOMTYPE];   int atnum=0;
    BOND_CL     bond[MAX_FF_BOND];     int bondnum_cl=0;
    ANGLE_CL    angle[MAX_FF_ANGLE];   int anglenum=0;
    TORSION_CL  tor[MAX_FF_TORSION];   int tornum=0;
    TORSION_CL  improp[MAX_FF_IMPROPER]; int impropnum=0;
    ATOMTYPE_CL vdw[MAX_FF_VDW];       int vdwnum=0;

    int    writeflag=1, typeflag=0, tmpint;
    char   name1[NAMELEN],name2[NAMELEN],name3[NAMELEN],name4[NAMELEN];
    char   lline[MAXCHAR], line2[MAXCHAR], tmpchar[MAXCHAR];
    double tmpf1, tmpf2, tmpf3;

    strcat(command,filename); strcat(command," ANTECHAMBER.FRCMOD");
    status=system(command);
    if (status!=0){
        fprintf(stdout,"Error: cannot run \"%s\" in cleanup_frcmod, exit\n",command); exit(1);
    }
    if ((fp1=fopen("ANTECHAMBER.FRCMOD","r"))==NULL){
        fprintf(stdout,"Error: Cannot open ANTECHAMBER.FRCMOD for cleanup\n"); exit(1);
    }
    if ((fp2=fopen(filename,"w"))==NULL){
        fprintf(stdout,"Error: Cannot open %s for writing\n",filename); exit(1);
    }

    while (fgets(lline,MAXCHAR,fp1)){
        sscanf(lline,"%s",tmpchar);
        if (!strcmp("MASS",    tmpchar)) typeflag=1;
        if (!strcmp("BOND",    tmpchar)) typeflag=2;
        if (!strcmp("ANGLE",   tmpchar)) typeflag=3;
        if (!strcmp("DIHE",    tmpchar)) typeflag=4;
        if (!strcmp("IMPROPER",tmpchar)) typeflag=5;
        if (!strcmp("NONBON",  tmpchar)) typeflag=6;
        if (strlen(lline)<=1){ typeflag=0; writeflag=1; }

        if (typeflag==1){
            writeflag=1;
            for (int ii=0;ii<atnum;ii++) if(!strcmp(tmpchar,at[ii].name)){writeflag=0;break;}
            if (writeflag) strcpy(at[atnum++].name,tmpchar);
        }
        if (typeflag==2){
            writeflag=1;
            strcpy(line2,lline);
            for (int ii=0;ii<6;ii++) if(line2[ii]=='-') line2[ii]=' ';
            sscanf(line2,"%s%s",name1,name2);
            for (int ii=0;ii<bondnum_cl;ii++){
                if(!strcmp(name1,bond[ii].name1)&&!strcmp(name2,bond[ii].name2)){writeflag=0;break;}
                if(!strcmp(name1,bond[ii].name2)&&!strcmp(name2,bond[ii].name1)){writeflag=0;break;}
            }
            if (writeflag){ strcpy(bond[bondnum_cl].name1,name1); strcpy(bond[bondnum_cl++].name2,name2); }
        }
        if (typeflag==3){
            writeflag=1;
            strcpy(line2,lline);
            for (int ii=0;ii<9;ii++) if(line2[ii]=='-') line2[ii]=' ';
            sscanf(line2,"%s%s%s",name1,name2,name3);
            for (int ii=0;ii<anglenum;ii++){
                if(!strcmp(name1,angle[ii].name1)&&!strcmp(name2,angle[ii].name2)&&!strcmp(name3,angle[ii].name3)){writeflag=0;break;}
                if(!strcmp(name1,angle[ii].name3)&&!strcmp(name2,angle[ii].name2)&&!strcmp(name3,angle[ii].name1)){writeflag=0;break;}
            }
            if (writeflag){
                strcpy(angle[anglenum].name1,name1);
                strcpy(angle[anglenum].name2,name2);
                strcpy(angle[anglenum++].name3,name3);
            }
        }
        if (typeflag==4){
            writeflag=1;
            strcpy(line2,lline);
            for (int ii=0;ii<12;ii++) if(line2[ii]=='-') line2[ii]=' ';
            sscanf(line2,"%s%s%s%s%d%lf%lf%lf",name1,name2,name3,name4,&tmpint,&tmpf1,&tmpf2,&tmpf3);
            for (int ii=0;ii<tornum;ii++){
                if(!strcmp(name1,tor[ii].name1)&&!strcmp(name2,tor[ii].name2)&&
                   !strcmp(name3,tor[ii].name3)&&!strcmp(name4,tor[ii].name4)&&tor[ii].fterm==tmpf3){writeflag=0;break;}
                if(!strcmp(name1,tor[ii].name4)&&!strcmp(name2,tor[ii].name3)&&
                   !strcmp(name3,tor[ii].name2)&&!strcmp(name4,tor[ii].name1)&&tor[ii].fterm==tmpf3){writeflag=0;break;}
            }
            if (writeflag){
                strcpy(tor[tornum].name1,name1); strcpy(tor[tornum].name2,name2);
                strcpy(tor[tornum].name3,name3); strcpy(tor[tornum].name4,name4);
                tor[tornum++].fterm=tmpf3;
            }
        }
        if (typeflag==5){
            writeflag=1;
            strcpy(line2,lline);
            for (int ii=0;ii<12;ii++) if(line2[ii]=='-') line2[ii]=' ';
            sscanf(line2,"%s%s%s%s",name1,name2,name3,name4);
            for (int ii=0;ii<impropnum;ii++){
                if(!strcmp(name1,improp[ii].name1)&&!strcmp(name2,improp[ii].name2)&&
                   !strcmp(name3,improp[ii].name3)&&!strcmp(name4,improp[ii].name4)){writeflag=0;break;}
                if(!strcmp(name1,improp[ii].name4)&&!strcmp(name2,improp[ii].name3)&&
                   !strcmp(name3,improp[ii].name2)&&!strcmp(name4,improp[ii].name1)){writeflag=0;break;}
            }
            if (writeflag){
                strcpy(improp[impropnum].name1,name1); strcpy(improp[impropnum].name2,name2);
                strcpy(improp[impropnum].name3,name3); strcpy(improp[impropnum++].name4,name4);
            }
        }
        if (typeflag==6){
            writeflag=1;
            for (int ii=0;ii<vdwnum;ii++) if(!strcmp(tmpchar,vdw[ii].name)){writeflag=0;break;}
            if (writeflag) strcpy(vdw[vdwnum++].name,tmpchar);
        }
        if (writeflag) fprintf(fp2,"%s",lline);
    }
    fclose(fp1); fclose(fp2);
}



/*========================================================================
 * chk_atomtype_core / chk_atomtype
 *======================================================================*/
static void chk_atomtype_core(const char *aname, int pid) {
    int suc=0, pid2; char type1[NAMELEN];
    for (int jj=0;jj<atomtypenum;jj++)
        if (!strcmp(aname,atomtype[jj].name)){
            suc=1;
            if (allparm_flag)
                fprintf(fpout,"%-2s %-8.3lf   %8.3lf\n",
                        aname,atomtype[jj].mass,atomtype[jj].pol);
            break;
        }
    if (!suc && parm[pid].nequa>1)
        for (int jj=1;jj<parm[pid].nequa&&!suc;jj++){
            pid2=parm[pid].equa[jj].pid; strcpy(type1,parm[pid2].atomtype);
            for (int kk=0;kk<atomtypenum2&&!suc;kk++)
                if (!strcmp(type1,atomtype[kk].name)){
                    suc=1;
                    fprintf(fpout,"%-2s %-8.3lf   %8.3lf               %s %-3s\n",
                            aname,parm[pid].mass,atomtype[kk].pol,"same as",atomtype[kk].name);
                    atomtype[atomtypenum]=atomtype[kk];
                    strcpy(atomtype[atomtypenum].name,aname);
                    atomtypenum++;
                    GROW_ARRAY(atomtype,atomtypenum,maxatomtype,MAX_FF_ATOMTYPE,ATOMTYPE);
                }
        }
    if (!suc && parm[pid].ncorr>1)
        for (int jj=1;jj<parm[pid].ncorr&&!suc;jj++){
            if (parm[pid].corr[jj].type<=1) continue;
            pid2=parm[pid].corr[jj].pid; strcpy(type1,parm[pid2].atomtype);
            for (int kk=0;kk<atomtypenum2&&!suc;kk++)
                if (!strcmp(type1,atomtype[kk].name)){
                    suc=1;
                    fprintf(fpout,"%-2s %-8.3lf   %8.3lf               %s %-3s\n",
                            aname,parm[pid].mass,atomtype[kk].pol,"same as",atomtype[kk].name);
                    atomtype[atomtypenum]=atomtype[kk];
                    strcpy(atomtype[atomtypenum].name,aname);
                    atomtypenum++;
                    GROW_ARRAY(atomtype,atomtypenum,maxatomtype,MAX_FF_ATOMTYPE,ATOMTYPE);
                }
        }
    if (!suc){
        fprintf(fpout,"%-2s %-8.3lf   %8.3lf               %s\n",
                aname,parm[pid].mass,0.0,"ATTN, no polarizability parameter");
        strcpy(atomtype[atomtypenum].name,aname);
        atomtypenum++;
        GROW_ARRAY(atomtype,atomtypenum,maxatomtype,MAX_FF_ATOMTYPE,ATOMTYPE);
    }
}

void chk_atomtype(void) {
    fprintf(fpout,"%s\n","MASS");
    if (iformat<=FILETYPE_MOL2) {
        for (int ii=0;ii<atomnum;ii++)
            chk_atomtype_core(atom[ii].ambername,parmid[ii]);
    } else {
        for (int ii=0;ii<r_atomtype_num;ii++){
            if (attn_opt==2 && !r_atomtype[ii].attn){
                fprintf(fpout,"%-2s %-8.3lf   %8.3lf\n",
                        r_atomtype[ii].name,r_atomtype[ii].mass,r_atomtype[ii].pol);
                continue;
            }
            chk_atomtype_core(r_atomtype[ii].name,lookup_parmid(r_atomtype[ii].name));
        }
    }
}

/*========================================================================
 * vdw helpers
 *======================================================================*/
static int vdw_lookup(const char *at_name, const char *corr_name, int nparm) {
    for (int ii=0;ii<nparm;ii++)
        if (!strcmp(vdwparm[ii].name,corr_name)){
            fprintf(fpout,"  %-2s%16.4lf%8.4lf             %s %-3s\n",
                    at_name,vdwparm[ii].radius,vdwparm[ii].pot,"same as",vdwparm[ii].name);
            vdwparm[vdwparmnum]=vdwparm[ii];
            strcpy(vdwparm[vdwparmnum].name,at_name);
            vdwparmnum++;
            GROW_ARRAY(vdwparm,vdwparmnum,maxvdwparm,MAX_FF_VDW,VDW);
            return 1;
        }
    return 0;
}

static void chk_vdw_core(const char *aname, int pid) {
    int suc=0, pid2;
    for (int jj=0;jj<vdwparmnum&&!suc;jj++)
        if (!strcmp(vdwparm[jj].name,aname)){
            suc=1;
            if (allparm_flag)
                fprintf(fpout,"  %-2s%16.4lf%8.4lf\n",
                        vdwparm[jj].name,vdwparm[jj].radius,vdwparm[jj].pot);
        }
    if (!suc && parm[pid].nequa>1)
        for (int jj=1;jj<parm[pid].nequa&&!suc;jj++){
            pid2=parm[pid].equa[jj].pid;
            suc=vdw_lookup(aname,parm[pid2].atomtype,vdwparmnum2);
        }
    if (!suc && parm[pid].ncorr>1)
        for (int jj=1;jj<parm[pid].ncorr&&!suc;jj++){
            if (parm[pid].corr[jj].type<=1) continue;
            pid2=parm[pid].corr[jj].pid;
            suc=vdw_lookup(aname,parm[pid2].atomtype,vdwparmnum2);
        }
    if (!suc){
        fprintf(fpout,"  %-2s%16.4lf%8.4lf             %s\n",
                aname,0.0,0.0,"ATTN, need revision");
        strcpy(vdwparm[vdwparmnum].name,aname);
        vdwparmnum++;
        GROW_ARRAY(vdwparm,vdwparmnum,maxvdwparm,MAX_FF_VDW,VDW);
    }
}

void chk_vdw(void) {
    fprintf(fpout,"\n%s\n","NONBON");
    if (iformat<=FILETYPE_MOL2) {
        for (int ii=0;ii<atomnum;ii++) chk_vdw_core(atom[ii].ambername,parmid[ii]);
    } else {
        for (int ii=0;ii<r_vdwparm_num;ii++){
            if (attn_opt==2 && !r_vdwparm[ii].attn){
                fprintf(fpout,"  %-2s%16.4lf%8.4lf\n",
                        r_vdwparm[ii].name,r_vdwparm[ii].radius,r_vdwparm[ii].pot);
                continue;
            }
            chk_vdw_core(r_vdwparm[ii].name,lookup_parmid(r_vdwparm[ii].name));
        }
    }
    fprintf(fpout,"\n\n\n");
}

/*========================================================================
 * Pure finder functions — no side effects, return index or -1.
 *======================================================================*/
static int bond_find(const char *c1, const char *c2, int nparm) {
    for (int kk=0;kk<nparm;kk++)
        if ((!strcmp(bondparm[kk].name1,c1)&&!strcmp(bondparm[kk].name2,c2))||
            (!strcmp(bondparm[kk].name1,c2)&&!strcmp(bondparm[kk].name2,c1)))
            return kk;
    return -1;
}

static int angle_find(const char *c1, const char *c2, const char *c3, int nparm) {
    for (int kk=0;kk<nparm;kk++)
        if ((!strcmp(angleparm[kk].name1,c1)&&!strcmp(angleparm[kk].name2,c2)&&
              !strcmp(angleparm[kk].name3,c3))||
            (!strcmp(angleparm[kk].name3,c1)&&!strcmp(angleparm[kk].name2,c2)&&
              !strcmp(angleparm[kk].name1,c3)))
            return kk;
    return -1;
}

static int torsion_find(const char *c1, const char *c2,
                        const char *c3, const char *c4, int nparm) {
    for (int kk=0;kk<nparm;kk++)
        if ((!strcmp(torsionparm[kk].name1,c1)&&!strcmp(torsionparm[kk].name2,c2)&&
              !strcmp(torsionparm[kk].name3,c3)&&!strcmp(torsionparm[kk].name4,c4))||
            (!strcmp(torsionparm[kk].name4,c1)&&!strcmp(torsionparm[kk].name3,c2)&&
              !strcmp(torsionparm[kk].name2,c3)&&!strcmp(torsionparm[kk].name1,c4)))
            return kk;
    return -1;
}

/* print_torsion: output one torsion entry, handles fterm<0 chains */
static void print_torsion(int tid, const char *n1, const char *n2,
                           const char *n3, const char *n4, double sc) {
    while (torsionparm[tid].fterm<0){
        fprintf(fpout,"%-2s-%-2s-%-2s-%-2s%4d%9.3lf%14.3lf%16.3lf"
                "      same as %-2s-%-2s-%-2s-%-2s\n",
                n1,n2,n3,n4,torsionparm[tid].mul,torsionparm[tid].force,
                torsionparm[tid].phase,torsionparm[tid].fterm,
                torsionparm[tid].name1,torsionparm[tid].name2,
                torsionparm[tid].name3,torsionparm[tid].name4);
        torsionparm[torsionparmnum]=torsionparm[tid];
        strcpy(torsionparm[torsionparmnum].name1,n1);
        strcpy(torsionparm[torsionparmnum].name2,n2);
        strcpy(torsionparm[torsionparmnum].name3,n3);
        strcpy(torsionparm[torsionparmnum].name4,n4);
        torsionparmnum++;
        GROW_ARRAY(torsionparm,torsionparmnum,maxtorsionparm,MAX_FF_TORSION,TORSION);
        tid++;
    }
    fprintf(fpout,"%-2s-%-2s-%-2s-%-2s%4d%9.3lf%14.3lf%16.3lf"
            "      same as %-2s-%-2s-%-2s-%-2s, penalty score=%5.1lf\n",
            n1,n2,n3,n4,torsionparm[tid].mul,torsionparm[tid].force,
            torsionparm[tid].phase,torsionparm[tid].fterm,
            torsionparm[tid].name1,torsionparm[tid].name2,
            torsionparm[tid].name3,torsionparm[tid].name4,sc);
    torsionparm[torsionparmnum]=torsionparm[tid];
    strcpy(torsionparm[torsionparmnum].name1,n1);
    strcpy(torsionparm[torsionparmnum].name2,n2);
    strcpy(torsionparm[torsionparmnum].name3,n3);
    strcpy(torsionparm[torsionparmnum].name4,n4);
    torsionparmnum++;
    GROW_ARRAY(torsionparm,torsionparmnum,maxtorsionparm,MAX_FF_TORSION,TORSION);
}

/*------------------------------------------------------------------------
 * VLOG_SETS: print the candidate atom type sets for a CORR/EQUA step.
 * pid_arr: array of pids, n: count, label: "1", "2", "X", etc.
 * For CORR entries, resolve pid=-1 by name lookup before printing.
 *----------------------------------------------------------------------*/
#define VLOG_CORR_SET(atom_pid, ncorr, label) \
    do { \
        fprintf(stderr,"    %s =", label); \
        for(int _x=0;_x<(ncorr);_x++){ \
            int _p=parm[atom_pid].corr[_x].pid; \
            if(_p<0) _p=lookup_parmid(parm[atom_pid].corr[_x].atomtype); \
            fprintf(stderr," %s",parm[_p].atomtype); \
        } \
        fprintf(stderr,"\n"); \
    } while(0)

#define VLOG_EQUA_SET(atom_pid, label) \
    do { \
        fprintf(stderr,"    %s =", label); \
        for(int _x=0;_x<parm[atom_pid].nequa;_x++) \
            fprintf(stderr," %s",parm[parm[atom_pid].equa[_x].pid].atomtype); \
        fprintf(stderr,"\n"); \
    } while(0)

/*========================================================================
 * chk_bond_core
 * Step 1: exact match
 * Step 2: EQUA single-winner
 * Step 3: CORR ranked → empirical
 * epid1/2: EQUA substitute pids; cpid1/2: CORR result pids
 *======================================================================*/
void chk_bond_core(BondQuery *queries, int nqueries) {
    int suc, epid1, epid2, cpid1, cpid2;
    char type1[NAMELEN], type2[NAMELEN];
    double score, score1, score2, bestscore;

    fprintf(fpout,"\n%s\n","BOND");
    for (int ii=0;ii<nqueries;ii++){
        BondQuery *q=&queries[ii];
        const char *name1=q->name1, *name2=q->name2;

        VLOG(1,"BOND %s-%s:\n",name1,name2);

        if (iformat>=FILETYPE_FRCMOD && attn_opt==2 && !q->attn){
            fprintf(fpout,"%-2s-%-2s%8.2lf%8.3lf\n",
                    name1,name2,r_bondparm[ii].force,q->length);
            continue;
        }

        /* Step 1: exact match */
        VLOG(2,"  Step 1 exact: %s-%s\n",name1,name2);
        {
            int ei=bond_find(name1,name2,bondparmnum);
            if (ei>=0){
                if (allparm_flag)
                    fprintf(fpout,"%-2s-%-2s%8.2lf%8.3lf\n",
                            name1,name2,bondparm[ei].force,bondparm[ei].length);
                VLOG(1,"  -> %s-%s exact match\n",name1,name2);
                continue;
            }
        }

        int pid1=get_pid(q->atid1,name1);
        int pid2=get_pid(q->atid2,name2);

        /* Step 2: EQUA single-winner */
        if (verbose>=2){
            fprintf(stderr,"  Step 2 EQUA: %d x %d\n",
                    parm[pid1].nequa,parm[pid2].nequa);
            VLOG_EQUA_SET(pid1,"1");
            VLOG_EQUA_SET(pid2,"2");
        }
        suc=0; bestblid=-1;
        for (int m=0;m<parm[pid1].nequa&&!suc;m++){
            epid1=parm[pid1].equa[m].pid; strcpy(type1,parm[epid1].atomtype);
            for (int n=0;n<parm[pid2].nequa&&!suc;n++){
                if (m==0&&n==0) continue;
                epid2=parm[pid2].equa[n].pid; strcpy(type2,parm[epid2].atomtype);
                int idx=bond_find(type1,type2,bondparmnum2);
                if (idx>=0){ bestblid=idx; bestscore=0; suc=1;
                    VLOG(1,"  -> %s-%s EQUA match penalty=0\n",type1,type2); }
            }
        }

        /* Step 3: CORR ranked search */
        if (!suc){
            int nc1=parm[pid1].ncorr, nc2=parm[pid2].ncorr;
            if (verbose>=2){
                fprintf(stderr,"  Step 3 CORR: %d x %d%s\n",nc1,nc2,
                        equa_corr_inherit?" (incl. EQUA-inherited)":"");
                VLOG_CORR_SET(pid1,nc1,"1");
                VLOG_CORR_SET(pid2,nc2,"2");
            }
            RankedMatch best[MAXBEST]; int nbest=0;
            for (int m=0;m<nc1;m++){
                cpid1=parm[pid1].corr[m].pid;
                if (cpid1<0) cpid1=lookup_parmid(parm[pid1].corr[m].atomtype);
                strcpy(type1,parm[cpid1].atomtype);
                score1=parm[pid1].corr[m].bl*wt.BL+parm[pid1].corr[m].blf*wt.BLF;
                for (int n=0;n<nc2;n++){
                    if (parm[pid1].corr[m].type<=1&&parm[pid2].corr[n].type<=1) continue;
                    cpid2=parm[pid2].corr[n].pid;
                    if (cpid2<0) cpid2=lookup_parmid(parm[pid2].corr[n].atomtype);
                    strcpy(type2,parm[cpid2].atomtype);
                    score2=parm[pid2].corr[n].bl*wt.BL+parm[pid2].corr[n].blf*wt.BLF;
                    score=score1+score2;
                    if (parm[cpid1].group!=parm[cpid2].group) score+=wt.GROUP;
                    score+=equtype_penalty(pid1,pid2,cpid1,cpid2);
                    int idx=bond_find(type1,type2,bondparmnum2);
                    if (idx>=0) RANKED_INSERT(best,nbest,idx,score);
                }
            }
            if (nbest>0){
                suc=1;
                int wi=best[0].idx;
                VLOG(1,"  -> %s-%s CORR match penalty=%.1f\n",
                     bondparm[wi].name1,bondparm[wi].name2,best[0].score);
                for (int r=1;r<nbest;r++)
                    VLOG(1,"     %s-%s penalty=%.1f (alt %d)\n",
                         bondparm[best[r].idx].name1,bondparm[best[r].idx].name2,
                         best[r].score,r);
                fprintf(fpout,"%-2s-%-2s%8.2lf%8.3lf       same as %2s-%2s,"
                        " penalty score=%5.1lf\n",
                        name1,name2,bondparm[wi].force,bondparm[wi].length,
                        bondparm[wi].name1,bondparm[wi].name2,best[0].score);
                bondparm[bondparmnum]=bondparm[wi];
                strcpy(bondparm[bondparmnum].name1,name1);
                strcpy(bondparm[bondparmnum].name2,name2);
                bondparmnum++;
                GROW_ARRAY(bondparm,bondparmnum,maxbondparm,MAX_FF_BOND,BOND_FF);
                for (int r=1;r<nbest;r++){
                    int ai=best[r].idx;
                    fprintf(fpout,"# alt[%d]: %-2s-%-2s%8.2lf%8.3lf, score=%5.1lf\n",
                            r,bondparm[ai].name1,bondparm[ai].name2,
                            bondparm[ai].force,bondparm[ai].length,best[r].score);
                }
                continue;
            }
        }

        /* EQUA winner output */
        if (suc&&bestblid>=0){
            fprintf(fpout,"%-2s-%-2s%8.2lf%8.3lf       same as %2s-%2s,"
                    " penalty score=%5.1lf\n",
                    name1,name2,bondparm[bestblid].force,bondparm[bestblid].length,
                    bondparm[bestblid].name1,bondparm[bestblid].name2,bestscore);
            bondparm[bondparmnum]=bondparm[bestblid];
            strcpy(bondparm[bondparmnum].name1,name1);
            strcpy(bondparm[bondparmnum].name2,name2);
            bondparmnum++;
            GROW_ARRAY(bondparm,bondparmnum,maxbondparm,MAX_FF_BOND,BOND_FF);
            continue;
        }

        /* Step 4: empirical */
        {
            int an1=parm[pid1].atomicnum, an2=parm[pid2].atomicnum;
            VLOG(2,"  Step 4 empirical: Z=%d-Z=%d\n",an1,an2);
            suc=empbond(name1,name2,name1,name2,an1,an2,q->length);
            if (suc){
                VLOG(1,"  -> %s-%s empirical\n",name1,name2);
            } else {
                if (q->length<=0)
                    VLOG(2,"    no bond length (topology-only input)\n");
                else
                    VLOG(2,"    failed (bl=%.4f, missing Z=%d-Z=%d?)\n",
                         q->length,an1,an2);
            }
        }
        if (!suc){
            VLOG(1,"  -> %s-%s ATTN\n",name1,name2);
            fprintf(fpout,"%-2s-%-2s%8.2lf%8.3lf       %s\n",
                    name1,name2,0.0,0.0,"ATTN, need revision");
            strcpy(bondparm[bondparmnum].name1,name1);
            strcpy(bondparm[bondparmnum].name2,name2);
            bondparmnum++;
            GROW_ARRAY(bondparm,bondparmnum,maxbondparm,MAX_FF_BOND,BOND_FF);
        }
    }
}

/*========================================================================
 * chk_angle_core
 * epid1/2/3: EQUA substitutes; cpid1/2/3: CORR result pids
 *======================================================================*/
void chk_angle_core(AngleQuery *queries, int nqueries) {
    int suc, suc2, epid1, epid2, epid3, cpid1, cpid2, cpid3;
    char type1[NAMELEN], type2[NAMELEN], type3[NAMELEN];
    double score, score1, score2, score3, bestscore;
    int an1, an2, an3;

    fprintf(fpout,"\n%s\n","ANGLE");
    for (int ii=0;ii<nqueries;ii++){
        AngleQuery *q=&queries[ii];
        const char *name1=q->name1,*name2=q->name2,*name3=q->name3;

        VLOG(1,"ANGLE %s-%s-%s:\n",name1,name2,name3);

        if (iformat>=FILETYPE_FRCMOD && attn_opt==2 && !q->attn){
            fprintf(fpout,"%-2s-%-2s-%-2s%9.3lf%12.3lf\n",
                    name1,name2,name3,r_angleparm[ii].force,r_angleparm[ii].angle);
            continue;
        }

        /* Step 1: exact match */
        VLOG(2,"  Step 1 exact: %s-%s-%s\n",name1,name2,name3);
        {
            int ei=angle_find(name1,name2,name3,angleparmnum);
            if (ei>=0){
                if (allparm_flag)
                    fprintf(fpout,"%-2s-%-2s-%-2s%9.3lf%12.3lf\n",
                            name1,name2,name3,angleparm[ei].force,angleparm[ei].angle);
                VLOG(1,"  -> %s-%s-%s exact match\n",name1,name2,name3);
                continue;
            }
        }

        int pid1=get_pid(q->atid1,name1);
        int pid2=get_pid(q->atid2,name2);
        int pid3=get_pid(q->atid3,name3);
        int nc1=parm[pid1].ncorr, nc2=parm[pid2].ncorr, nc3=parm[pid3].ncorr;

        /* Step 2: EQUA single-winner — central atom first */
        if (verbose>=2){
            fprintf(stderr,"  Step 2 EQUA: %d x %d x %d\n",
                    parm[pid1].nequa,parm[pid2].nequa,parm[pid3].nequa);
            VLOG_EQUA_SET(pid1,"1");
            VLOG_EQUA_SET(pid2,"2");
            VLOG_EQUA_SET(pid3,"3");
        }
        suc=0; bestbaid=-1;
        for (int m=0;m<parm[pid2].nequa&&!suc;m++){
            epid2=parm[pid2].equa[m].pid; strcpy(type2,parm[epid2].atomtype);
            for (int n=0;n<parm[pid1].nequa&&!suc;n++){
                epid1=parm[pid1].equa[n].pid; strcpy(type1,parm[epid1].atomtype);
                for (int o=0;o<parm[pid3].nequa&&!suc;o++){
                    if (m==0&&n==0&&o==0) continue;
                    epid3=parm[pid3].equa[o].pid; strcpy(type3,parm[epid3].atomtype);
                    int idx=angle_find(type1,type2,type3,angleparmnum2);
                    if (idx>=0){ bestba=angleparm[idx]; bestbaid=idx;
                        bestscore=0; suc=1;
                        VLOG(1,"  -> %s-%s-%s EQUA match penalty=0\n",
                             type1,type2,type3); }
                }
            }
        }
        if (suc&&bestbaid>=0){
            fprintf(fpout,"%-2s-%-2s-%-2s%9.3lf%12.3lf   same as %-2s-%-2s-%-2s,"
                    " penalty score=%5.1lf\n",
                    name1,name2,name3,bestba.force,bestba.angle,
                    bestba.name1,bestba.name2,bestba.name3,bestscore);
            angleparm[angleparmnum]=bestba;
            strcpy(angleparm[angleparmnum].name1,name1);
            strcpy(angleparm[angleparmnum].name2,name2);
            strcpy(angleparm[angleparmnum].name3,name3);
            angleparmnum++;
            GROW_ARRAY(angleparm,angleparmnum,maxangleparm,MAX_FF_ANGLE,ANGLE);
            bestbaid=-1; continue;
        }

        /* Step 3: CORR ranked search */
        if (!suc){
            if (verbose>=2){
                fprintf(stderr,"  Step 3 CORR: %d x %d x %d%s\n",nc1,nc2,nc3,
                        equa_corr_inherit?" (incl. EQUA-inherited)":"");
                VLOG_CORR_SET(pid1,nc1,"1");
                VLOG_CORR_SET(pid2,nc2,"2");
                VLOG_CORR_SET(pid3,nc3,"3");
            }
            RankedMatch best[MAXBEST]; int nbest=0;
            for (int m=0;m<nc1;m++){
                cpid1=parm[pid1].corr[m].pid;
                if (cpid1<0) cpid1=lookup_parmid(parm[pid1].corr[m].atomtype);
                strcpy(type1,parm[cpid1].atomtype);
                score1=parm[pid1].corr[m].ba*wt.BA+parm[pid1].corr[m].baf*wt.BAF;
                for (int n=0;n<nc2;n++){
                    cpid2=parm[pid2].corr[n].pid;
                    if (cpid2<0) cpid2=lookup_parmid(parm[pid2].corr[n].atomtype);
                    strcpy(type2,parm[cpid2].atomtype);
                    score2=(parm[pid2].corr[n].cba*wt.BA+
                            parm[pid2].corr[n].cbaf*wt.BAF)*wt.BA_CTR;
                    for (int o=0;o<nc3;o++){
                        if (parm[pid1].corr[m].type<=1&&parm[pid2].corr[n].type<=1&&
                            parm[pid3].corr[o].type<=1) continue;
                        cpid3=parm[pid3].corr[o].pid;
                        if (cpid3<0) cpid3=lookup_parmid(parm[pid3].corr[o].atomtype);
                        strcpy(type3,parm[cpid3].atomtype);
                        score3=parm[pid3].corr[o].ba*wt.BA+parm[pid3].corr[o].baf*wt.BAF;
                        score=score1+score2+score3+wt.GROUP;
                        if (parm[cpid1].group==parm[cpid2].group&&
                            parm[cpid1].group==parm[cpid3].group) score-=wt.GROUP;
                        if (score<=THRESHOLD_BA){
                            int idx=angle_find(type1,type2,type3,angleparmnum2);
                            if (idx>=0) RANKED_INSERT(best,nbest,idx,score);
                        }
                    }
                }
            }
            if (nbest>0){
                suc=1;
                int wi=best[0].idx;
                VLOG(1,"  -> %s-%s-%s CORR match penalty=%.1f\n",
                     angleparm[wi].name1,angleparm[wi].name2,angleparm[wi].name3,
                     best[0].score);
                for (int r=1;r<nbest;r++)
                    VLOG(1,"     %s-%s-%s penalty=%.1f (alt %d)\n",
                         angleparm[best[r].idx].name1,angleparm[best[r].idx].name2,
                         angleparm[best[r].idx].name3,best[r].score,r);
                fprintf(fpout,"%-2s-%-2s-%-2s%9.3lf%12.3lf   same as %-2s-%-2s-%-2s,"
                        " penalty score=%5.1lf\n",
                        name1,name2,name3,angleparm[wi].force,angleparm[wi].angle,
                        angleparm[wi].name1,angleparm[wi].name2,angleparm[wi].name3,
                        best[0].score);
                angleparm[angleparmnum]=angleparm[wi];
                strcpy(angleparm[angleparmnum].name1,name1);
                strcpy(angleparm[angleparmnum].name2,name2);
                strcpy(angleparm[angleparmnum].name3,name3);
                angleparmnum++;
                GROW_ARRAY(angleparm,angleparmnum,maxangleparm,MAX_FF_ANGLE,ANGLE);
                for (int r=1;r<nbest;r++){
                    int ai=best[r].idx;
                    fprintf(fpout,"# alt[%d]: %-2s-%-2s-%-2s%9.3lf%12.3lf, score=%5.1lf\n",
                            r,angleparm[ai].name1,angleparm[ai].name2,angleparm[ai].name3,
                            angleparm[ai].force,angleparm[ai].angle,best[r].score);
                }
                bestbaid=-1; continue;
            }
        }

        /* Step 4: empirical self */
        {
            an1=parm[pid1].atomicnum; an2=parm[pid2].atomicnum; an3=parm[pid3].atomicnum;
            VLOG(2,"  Step 4 empirical self: Z=%d-Z=%d-Z=%d (bl1=%.4f bl2=%.4f)\n",
                 an1,an2,an3,q->bl1,q->bl2);
            if (q->bl1<=0||q->bl2<=0){
                VLOG(2,"    no bond lengths (topology-only input)\n");
                goto angle_attn;
            }
            suc=empangle(name1,name2,name3,name1,name2,name3,
                         an1,an2,an3,q->bl1,q->bl2,q->angle,0);
            if (suc){ VLOG(1,"  -> %s-%s-%s empirical\n",name1,name2,name3); continue; }
            VLOG(2,"    empirical self failed\n");
        }

        /* Step 5: EQUA empirical single-winner */
        if (!suc){
            if (verbose>=2){
                fprintf(stderr,"  Step 5 EQUA empirical: %d x %d x %d\n",
                        parm[pid1].nequa,parm[pid2].nequa,parm[pid3].nequa);
                VLOG_EQUA_SET(pid1,"1");
                VLOG_EQUA_SET(pid2,"2");
                VLOG_EQUA_SET(pid3,"3");
            }
            bestbaid=-1;
            for (int m=0;m<parm[pid2].nequa&&!suc;m++){
                epid2=parm[pid2].equa[m].pid; strcpy(type2,parm[epid2].atomtype);
                for (int n=0;n<parm[pid1].nequa&&!suc;n++){
                    epid1=parm[pid1].equa[n].pid; strcpy(type1,parm[epid1].atomtype);
                    for (int o=0;o<parm[pid3].nequa&&!suc;o++){
                        if (m==0&&n==0&&o==0) continue;
                        epid3=parm[pid3].equa[o].pid; strcpy(type3,parm[epid3].atomtype);
                        an1=parm[epid1].atomicnum;
                        an2=parm[epid2].atomicnum;
                        an3=parm[epid3].atomicnum;
                        suc2=empangle(type1,type2,type3,name1,name2,name3,
                                      an1,an2,an3,q->bl1,q->bl2,q->angle,1);
                        if (suc2){ bestbaid=999999; bestscore=0; suc=1;
                            VLOG(1,"  -> %s-%s-%s EQUA empirical penalty=0\n",
                                 type1,type2,type3); }
                    }
                }
            }
        }
        if (suc&&bestbaid>=0){
            fprintf(fpout,"%-2s-%-2s-%-2s%9.3lf%12.3lf   Calculated using %-2s-%-2s-%-2s,"
                    " penalty score=%5.1lf\n",
                    name1,name2,name3,bestba.force,bestba.angle,
                    bestba.name1,bestba.name2,bestba.name3,bestscore);
            angleparm[angleparmnum]=bestba;
            strcpy(angleparm[angleparmnum].name1,name1);
            strcpy(angleparm[angleparmnum].name2,name2);
            strcpy(angleparm[angleparmnum].name3,name3);
            angleparmnum++;
            GROW_ARRAY(angleparm,angleparmnum,maxangleparm,MAX_FF_ANGLE,ANGLE);
            bestbaid=-1; continue;
        }

        /* Step 6: CORR empirical ranked */
        if (!suc){
            if (verbose>=2){
                fprintf(stderr,"  Step 6 CORR empirical: %d x %d x %d%s\n",nc1,nc2,nc3,
                        equa_corr_inherit?" (incl. EQUA-inherited)":"");
                VLOG_CORR_SET(pid1,nc1,"1");
                VLOG_CORR_SET(pid2,nc2,"2");
                VLOG_CORR_SET(pid3,nc3,"3");
            }
            RankedMatch best[MAXBEST]; int nbest=0;
            bestscore=INITSCORE;
            for (int m=0;m<nc1;m++){
                cpid1=parm[pid1].corr[m].pid;
                if (cpid1<0) cpid1=lookup_parmid(parm[pid1].corr[m].atomtype);
                strcpy(type1,parm[cpid1].atomtype);
                score1=parm[pid1].corr[m].ba*wt.BA+parm[pid1].corr[m].baf*wt.BAF;
                for (int n=0;n<nc2;n++){
                    cpid2=parm[pid2].corr[n].pid;
                    if (cpid2<0) cpid2=lookup_parmid(parm[pid2].corr[n].atomtype);
                    strcpy(type2,parm[cpid2].atomtype);
                    score2=(parm[pid2].corr[n].cba*wt.BA+
                            parm[pid2].corr[n].cbaf*wt.BAF)*wt.BA_CTR;
                    for (int o=0;o<nc3;o++){
                        if (parm[pid1].corr[m].type<=1&&parm[pid2].corr[n].type<=1&&
                            parm[pid3].corr[o].type<=1) continue;
                        cpid3=parm[pid3].corr[o].pid;
                        if (cpid3<0) cpid3=lookup_parmid(parm[pid3].corr[o].atomtype);
                        strcpy(type3,parm[cpid3].atomtype);
                        score3=parm[pid3].corr[o].ba*wt.BA+parm[pid3].corr[o].baf*wt.BAF;
                        score=score1+score2+score3+wt.GROUP;
                        if (parm[cpid1].group==parm[cpid2].group&&
                            parm[cpid1].group==parm[cpid3].group) score-=wt.GROUP;
                        if (score<bestscore){
                            an1=parm[cpid1].atomicnum;
                            an2=parm[cpid2].atomicnum;
                            an3=parm[cpid3].atomicnum;
                            int ok=empangle(type1,type2,type3,
                                            name1,name2,name3,
                                            an1,an2,an3,q->bl1,q->bl2,q->angle,1);
                            if (ok){ RANKED_INSERT(best,nbest,-1,score);
                                if (nbest==1) bestscore=score; }
                        }
                    }
                }
            }
            if (nbest>0&&bestbaid>=0){
                suc=1;
                VLOG(1,"  -> %s-%s-%s CORR empirical penalty=%.1f\n",
                     bestba.name1,bestba.name2,bestba.name3,best[0].score);
                fprintf(fpout,"%-2s-%-2s-%-2s%9.3lf%12.3lf   Calculated using %-2s-%-2s-%-2s,"
                        " penalty score=%5.1lf\n",
                        name1,name2,name3,bestba.force,bestba.angle,
                        bestba.name1,bestba.name2,bestba.name3,best[0].score);
                angleparm[angleparmnum]=bestba;
                strcpy(angleparm[angleparmnum].name1,name1);
                strcpy(angleparm[angleparmnum].name2,name2);
                strcpy(angleparm[angleparmnum].name3,name3);
                angleparmnum++;
                GROW_ARRAY(angleparm,angleparmnum,maxangleparm,MAX_FF_ANGLE,ANGLE);
                if (nbest>1)
                    fprintf(fpout,"# %d additional empirical candidates"
                            " (averaging planned)\n",nbest-1);
                bestbaid=-1; continue;
            }
        }

        angle_attn:
        VLOG(1,"  -> %s-%s-%s ATTN\n",name1,name2,name3);
        fprintf(fpout,"%-2s-%-2s-%-2s%9.3lf  %10.3lf   %s\n",
                name1,name2,name3,0.0,0.0,"ATTN, need revision");
        strcpy(angleparm[angleparmnum].name1,name1);
        strcpy(angleparm[angleparmnum].name2,name2);
        strcpy(angleparm[angleparmnum].name3,name3);
        angleparmnum++;
        GROW_ARRAY(angleparm,angleparmnum,maxangleparm,MAX_FF_ANGLE,ANGLE);
    }
}

/*========================================================================
 * chk_torsion_core
 * epid1-4: EQUA substitutes; cpid1-4: CORR result pids
 *======================================================================*/
void chk_torsion_core(TorsionQuery *queries, int nqueries) {
    int suc, epid1, epid2, epid3, epid4, cpid1, cpid2, cpid3, cpid4;
    char type1[NAMELEN], type2[NAMELEN], type3[NAMELEN], type4[NAMELEN];
    double score, score1, score2, score3, score4, bestscore;

    fprintf(fpout,"\n%s\n","DIHE");
    for (int ii=0;ii<nqueries;ii++){
        TorsionQuery *q=&queries[ii];
        const char *name1=q->name1,*name2=q->name2,
                   *name3=q->name3,*name4=q->name4;

        VLOG(1,"DIHE %s-%s-%s-%s:\n",name1,name2,name3,name4);

        if (iformat>=FILETYPE_FRCMOD && attn_opt==2 && !q->attn){
            fprintf(fpout,"%-2s-%-2s-%-2s-%-2s%4d%9.3lf%14.3lf%16.3lf\n",
                    name1,name2,name3,name4,
                    r_torsionparm[ii].mul,r_torsionparm[ii].force,
                    r_torsionparm[ii].phase,r_torsionparm[ii].fterm);
            continue;
        }

        int pid1=get_pid(q->atid1,name1);
        int pid2=get_pid(q->atid2,name2);
        int pid3=get_pid(q->atid3,name3);
        int pid4=get_pid(q->atid4,name4);
        int nc1=parm[pid1].ncorr, nc2=parm[pid2].ncorr;
        int nc3=parm[pid3].ncorr, nc4=parm[pid4].ncorr;

        /* Step 1: exact match */
        VLOG(2,"  Step 1 exact: %s-%s-%s-%s\n",name1,name2,name3,name4);
        {
            int idx=torsion_find(name1,name2,name3,name4,torsionparmnum);
            if (idx>=0){
                if (allparm_flag){
                    int nn=idx;
                    while(torsionparm[nn].fterm<0){
                        fprintf(fpout,"%-2s-%-2s-%-2s-%-2s%4d%9.3lf%14.3lf%16.3lf\n",
                                name1,name2,name3,name4,
                                torsionparm[nn].mul,torsionparm[nn].force,
                                torsionparm[nn].phase,torsionparm[nn].fterm);
                        nn++;
                    }
                    fprintf(fpout,"%-2s-%-2s-%-2s-%-2s%4d%9.3lf%14.3lf%16.3lf\n",
                            name1,name2,name3,name4,
                            torsionparm[nn].mul,torsionparm[nn].force,
                            torsionparm[nn].phase,torsionparm[nn].fterm);
                }
                VLOG(1,"  -> %s-%s-%s-%s exact match\n",name1,name2,name3,name4);
                continue;
            }
        }

        /* Step 2: EQUA specific single-winner */
        if (verbose>=2){
            fprintf(stderr,"  Step 2 EQUA specific: %d x %d x %d x %d\n",
                    parm[pid1].nequa,parm[pid2].nequa,parm[pid3].nequa,parm[pid4].nequa);
            VLOG_EQUA_SET(pid1,"1");
            VLOG_EQUA_SET(pid2,"2");
            VLOG_EQUA_SET(pid3,"3");
            VLOG_EQUA_SET(pid4,"4");
        }
        suc=0; besttorid=-1;
        for (int m=0;m<parm[pid2].nequa&&!suc;m++){
            epid2=parm[pid2].equa[m].pid; strcpy(type2,parm[epid2].atomtype);
            for (int n=0;n<parm[pid3].nequa&&!suc;n++){
                epid3=parm[pid3].equa[n].pid; strcpy(type3,parm[epid3].atomtype);
                for (int p=0;p<parm[pid1].nequa&&!suc;p++){
                    epid1=parm[pid1].equa[p].pid; strcpy(type1,parm[epid1].atomtype);
                    for (int qi=0;qi<parm[pid4].nequa&&!suc;qi++){
                        if (m==0&&n==0&&p==0&&qi==0) continue;
                        epid4=parm[pid4].equa[qi].pid; strcpy(type4,parm[epid4].atomtype);
                        int idx=torsion_find(type1,type2,type3,type4,torsionparmnum2);
                        if (idx>=0){ besttorid=idx; bestscore=0; suc=1;
                            VLOG(1,"  -> %s-%s-%s-%s EQUA specific penalty=0\n",
                                 type1,type2,type3,type4); }
                    }
                }
            }
        }
        if (suc&&besttorid>=0){
            print_torsion(besttorid,name1,name2,name3,name4,bestscore);
            besttorid=-1; continue;
        }

        /* Step 3: X exact — literal library match, silent unless allparm_flag */
        VLOG(2,"  Step 3 X exact: X-%s-%s-X\n",name2,name3);
        {
            int idx=torsion_find("X",name2,name3,"X",torsionparmnum2);
            if (idx>=0){
                VLOG(1,"  -> X-%s-%s-X exact match\n",name2,name3);
                if (allparm_flag)
                    print_torsion(idx,name1,name2,name3,name4,0.0);
                else {
                    torsionparm[torsionparmnum]=torsionparm[idx];
                    strcpy(torsionparm[torsionparmnum].name1,name1);
                    strcpy(torsionparm[torsionparmnum].name2,name2);
                    strcpy(torsionparm[torsionparmnum].name3,name3);
                    strcpy(torsionparm[torsionparmnum].name4,name4);
                    torsionparmnum++;
                    GROW_ARRAY(torsionparm,torsionparmnum,maxtorsionparm,MAX_FF_TORSION,TORSION);
                }
                continue;
            }
        }

        /* Step 4: X EQUA — resolved by equivalence, print with penalty=0 */
        if (verbose>=2){
            fprintf(stderr,"  Step 4 X EQUA: %d x %d\n",
                    parm[pid2].nequa,parm[pid3].nequa);
            fprintf(stderr,"    1 = X\n");
            VLOG_EQUA_SET(pid2,"2");
            VLOG_EQUA_SET(pid3,"3");
            fprintf(stderr,"    4 = X\n");
        }
        besttorid=-1;
        for (int n=0;n<parm[pid2].nequa&&!suc;n++){
            epid2=parm[pid2].equa[n].pid; strcpy(type2,parm[epid2].atomtype);
            for (int p=0;p<parm[pid3].nequa&&!suc;p++){
                if (n==0&&p==0) continue;
                epid3=parm[pid3].equa[p].pid; strcpy(type3,parm[epid3].atomtype);
                int idx=torsion_find("X",type2,type3,"X",torsionparmnum2);
                if (idx>=0){ besttorid=idx; bestscore=0; suc=1;
                    VLOG(1,"  -> X-%s-%s-X EQUA match penalty=0\n",type2,type3); }
            }
        }
        if (suc&&besttorid>=0){
            print_torsion(besttorid,name1,name2,name3,name4,bestscore);
            besttorid=-1; continue;
        }

        /* Step 5: CORR specific ranked */
        if (!suc){
            if (verbose>=2){
                fprintf(stderr,"  Step 5 CORR specific: %d x %d x %d x %d%s\n",
                        nc1,nc2,nc3,nc4,
                        equa_corr_inherit?" (incl. EQUA-inherited)":"");
                VLOG_CORR_SET(pid1,nc1,"1");
                VLOG_CORR_SET(pid2,nc2,"2");
                VLOG_CORR_SET(pid3,nc3,"3");
                VLOG_CORR_SET(pid4,nc4,"4");
            }
            RankedMatch best[MAXBEST]; int nbest=0;
            for (int m=0;m<nc1;m++){
                cpid1=parm[pid1].corr[m].pid;
                if (cpid1<0) cpid1=lookup_parmid(parm[pid1].corr[m].atomtype);
                strcpy(type1,parm[cpid1].atomtype);
                score1=parm[pid1].corr[m].tor;
                for (int n=0;n<nc2;n++){
                    cpid2=parm[pid2].corr[n].pid;
                    if (cpid2<0) cpid2=lookup_parmid(parm[pid2].corr[n].atomtype);
                    strcpy(type2,parm[cpid2].atomtype);
                    score2=parm[pid2].corr[n].ctor*wt.TOR_CTR;
                    for (int p=0;p<nc3;p++){
                        cpid3=parm[pid3].corr[p].pid;
                        if (cpid3<0) cpid3=lookup_parmid(parm[pid3].corr[p].atomtype);
                        strcpy(type3,parm[cpid3].atomtype);
                        score3=parm[pid3].corr[p].ctor*wt.TOR_CTR;
                        for (int qi=0;qi<nc4;qi++){
                            if (parm[pid1].corr[m].type<=1&&parm[pid2].corr[n].type<=1&&
                                parm[pid3].corr[p].type<=1&&parm[pid4].corr[qi].type<=1) continue;
                            cpid4=parm[pid4].corr[qi].pid;
                            if (cpid4<0) cpid4=lookup_parmid(parm[pid4].corr[qi].atomtype);
                            strcpy(type4,parm[cpid4].atomtype);
                            score4=parm[pid4].corr[qi].tor;
                            score=score1+score2+score3+score4+wt.GROUP;
                            score+=equtype_penalty(pid2,pid3,cpid2,cpid3);
                            if (parm[cpid1].group==parm[cpid2].group&&
                                parm[cpid1].group==parm[cpid3].group&&
                                parm[cpid1].group==parm[cpid4].group) score-=wt.GROUP;
                            int idx=torsion_find(type1,type2,type3,type4,torsionparmnum2);
                            if (idx>=0) RANKED_INSERT(best,nbest,idx,score);
                        }
                    }
                }
            }
            if (nbest>0){
                suc=1;
                VLOG(1,"  -> %s-%s-%s-%s CORR specific penalty=%.1f\n",
                     torsionparm[best[0].idx].name1,torsionparm[best[0].idx].name2,
                     torsionparm[best[0].idx].name3,torsionparm[best[0].idx].name4,
                     best[0].score);
                for (int r=1;r<nbest;r++)
                    VLOG(1,"     %s-%s-%s-%s penalty=%.1f (alt %d)\n",
                         torsionparm[best[r].idx].name1,torsionparm[best[r].idx].name2,
                         torsionparm[best[r].idx].name3,torsionparm[best[r].idx].name4,
                         best[r].score,r);
                print_torsion(best[0].idx,name1,name2,name3,name4,best[0].score);
                for (int r=1;r<nbest;r++){
                    int ai=best[r].idx;
                    fprintf(fpout,"# alt[%d]: %2s-%2s-%2s-%2s%4d%9.3lf%14.3lf%16.3lf, score=%5.1lf\n",
                            r,torsionparm[ai].name1,torsionparm[ai].name2,
                            torsionparm[ai].name3,torsionparm[ai].name4,
                            torsionparm[ai].mul,torsionparm[ai].force,
                            torsionparm[ai].phase,torsionparm[ai].fterm,
                            best[r].score);
                }
                besttorid=-1; continue;
            }
        }

        /* Step 6: CORR X ranked */
        if (!suc){
            if (verbose>=2){
                fprintf(stderr,"  Step 6 CORR X: %d x %d%s\n",nc2,nc3,
                        equa_corr_inherit?" (incl. EQUA-inherited)":"");
                fprintf(stderr,"    1 = X\n");
                VLOG_CORR_SET(pid2,nc2,"2");
                VLOG_CORR_SET(pid3,nc3,"3");
                fprintf(stderr,"    4 = X\n");
            }
            RankedMatch best[MAXBEST]; int nbest=0;
            for (int n=0;n<nc2;n++){
                cpid2=parm[pid2].corr[n].pid;
                if (cpid2<0) cpid2=lookup_parmid(parm[pid2].corr[n].atomtype);
                strcpy(type2,parm[cpid2].atomtype);
                score2=parm[pid2].corr[n].ctor*wt.TOR_CTR;
                for (int p=0;p<nc3;p++){
                    if (parm[pid2].corr[n].type<=1&&parm[pid3].corr[p].type<=1) continue;
                    cpid3=parm[pid3].corr[p].pid;
                    if (cpid3<0) cpid3=lookup_parmid(parm[pid3].corr[p].atomtype);
                    strcpy(type3,parm[cpid3].atomtype);
                    score3=parm[pid3].corr[p].ctor*wt.TOR_CTR;
                    score=score2+score3;
                    score+=equtype_penalty(pid2,pid3,cpid2,cpid3);
                    if (parm[cpid2].group!=parm[cpid3].group) score+=wt.GROUP;
                    int idx=torsion_find("X",type2,type3,"X",torsionparmnum2);
                    if (idx>=0) RANKED_INSERT(best,nbest,idx,score);
                }
            }
            if (nbest>0){
                suc=1;
                VLOG(1,"  -> X-%s-%s-X CORR match penalty=%.1f\n",
                     torsionparm[best[0].idx].name2,torsionparm[best[0].idx].name3,
                     best[0].score);
                for (int r=1;r<nbest;r++)
                    VLOG(1,"     X-%s-%s-X penalty=%.1f (alt %d)\n",
                         torsionparm[best[r].idx].name2,torsionparm[best[r].idx].name3,
                         best[r].score,r);
                print_torsion(best[0].idx,name1,name2,name3,name4,best[0].score);
                for (int r=1;r<nbest;r++){
                    int ai=best[r].idx;
                    fprintf(fpout,"# alt[%d]: X-%2s-%2s-X%4d%9.3lf%14.3lf%16.3lf, score=%5.1lf\n",
                            r,torsionparm[ai].name2,torsionparm[ai].name3,
                            torsionparm[ai].mul,torsionparm[ai].force,
                            torsionparm[ai].phase,torsionparm[ai].fterm,
                            best[r].score);
                }
                besttorid=-1; continue;
            }
        }

        VLOG(1,"  -> %s-%s-%s-%s ATTN\n",name1,name2,name3,name4);
        fprintf(fpout,"%-2s-%-2s-%-2s-%-2s%4d%9.3lf%14.3lf%16.3lf      %s\n",
                name1,name2,name3,name4,1,0.0,0.0,2.0,"ATTN, need revision");
        strcpy(torsionparm[torsionparmnum].name1,name1);
        strcpy(torsionparm[torsionparmnum].name2,name2);
        strcpy(torsionparm[torsionparmnum].name3,name3);
        strcpy(torsionparm[torsionparmnum].name4,name4);
        torsionparmnum++;
        GROW_ARRAY(torsionparm,torsionparmnum,maxtorsionparm,MAX_FF_TORSION,TORSION);
    }
}

/*========================================================================
 * chk_improper_core
 *======================================================================*/
void chk_improper_core(ImproperQuery *queries, int nqueries) {
    int suc;
    int pid1, pid2, pid3, pid4, pid5, pid6, pid7, pid8;
    char type1[NAMELEN], type2[NAMELEN], type3[NAMELEN], type4[NAMELEN];
    double score, score1, score2, score3, score4, bestscore;

    fprintf(fpout,"\n%s\n","IMPROPER");
    for (int ii=0;ii<nqueries;ii++){
        ImproperQuery *q=&queries[ii];
        const char *name1=q->name1,*name2=q->name2,*name3=q->name3,*name4=q->name4;

        if (iformat>=FILETYPE_FRCMOD && attn_opt==2 && !q->attn){
            fprintf(fpout,"%-2s-%-2s-%-2s-%-2s %11.1lf%15.1lf%12.1lf\n",
                    name1,name2,name3,name4,
                    r_improperparm[ii].force,r_improperparm[ii].phase,
                    r_improperparm[ii].fterm);
            continue;
        }

        pid1=get_pid(q->atid1,name1); pid2=get_pid(q->atid2,name2);
        pid3=get_pid(q->atid3,name3); pid4=get_pid(q->atid4,name4);
        int ipm_lim=(iformat<=FILETYPE_MOL2)?improperparmnum:improperparmnum2;

/*	step 1 check directly	*/
        suc=0;
        for (int jj=0;jj<improperparmnum;jj++)
            if (!strcmp(improperparm[jj].name1,name1)&&
                !strcmp(improperparm[jj].name2,name2)&&
                !strcmp(improperparm[jj].name3,name3)&&
                !strcmp(improperparm[jj].name4,name4)){
                suc=1;
                if (allparm_flag)
                    fprintf(fpout,"%-2s-%-2s-%-2s-%-2s %11.1lf%15.1lf%12.1lf\n",
                            name1,name2,name3,name4,
                            improperparm[jj].force,improperparm[jj].phase,
                            improperparm[jj].fterm);
                break;
            }
        if (suc) continue;

/* 	Step 2 check special improper torsional parameters using equal atom types */
        bestimproperid=-1;
        for (int m=0;m<parm[pid1].nequa&&!suc;m++){
            pid5=parm[pid1].equa[m].pid; strcpy(type1,parm[pid5].atomtype);
            for (int n=0;n<parm[pid2].nequa&&!suc;n++){
                pid6=parm[pid2].equa[n].pid; strcpy(type2,parm[pid6].atomtype);
                for (int p=0;p<parm[pid3].nequa&&!suc;p++){
                    pid7=parm[pid3].equa[p].pid; strcpy(type3,parm[pid7].atomtype);
                    for (int qi=0;qi<parm[pid4].nequa&&!suc;qi++){
                        if (m==0&&n==0&&p==0&&qi==0) continue;
                        pid8=parm[pid4].equa[qi].pid; strcpy(type4,parm[pid8].atomtype);
                        for (int kk=0;kk<ipm_lim;kk++){
                            if (improperparm[kk].numX>0) continue;
                            if (!strcmp(improperparm[kk].name1,type1)&&
                                !strcmp(improperparm[kk].name2,type2)&&
                                !strcmp(improperparm[kk].name3,type3)&&
                                !strcmp(improperparm[kk].name4,type4))
                                { bestimproperid=kk; bestscore=0; suc=1; break; }
                        }
                    }
                }
            }
        }
        if (suc&&bestimproperid>=0) goto improper_output;

/*	step 3 considering general atom types */
        if (!suc){
            bestimproperid=-1; bestscore=INITSCORE;
            for (int jj=0;jj<ipm_lim;jj++){
                if (!improperparm[jj].numX) continue;
                if (strcmp(improperparm[jj].name1,name1)!=0&&strcmp(improperparm[jj].name1,"X")!=0) continue;
                if (strcmp(improperparm[jj].name2,name2)!=0&&strcmp(improperparm[jj].name2,"X")!=0) continue;
                if (strcmp(improperparm[jj].name3,name3)!=0&&strcmp(improperparm[jj].name3,"X")!=0) continue;
                if (strcmp(improperparm[jj].name4,name4)!=0&&strcmp(improperparm[jj].name4,"X")!=0) continue;
                score1=(!strcmp(improperparm[jj].name1,"X"))?wt.X:0;
                score2=(!strcmp(improperparm[jj].name2,"X"))?wt.X:0;
                score3=(!strcmp(improperparm[jj].name3,"X"))?wt.X3:0;
                score4=(!strcmp(improperparm[jj].name4,"X"))?wt.X:0;
                score=score1+score2+score3+score4;
                if (score<bestscore){ bestimproperid=jj; bestscore=score; suc=1; }
            }
        }
        if (suc&&bestimproperid>=0) goto improper_output;

/*	Step 4 considering both equal atom types and general terms	*/
        if (!suc){
            bestscore=INITSCORE; bestimproperid=-1;
            for (int m=0;m<parm[pid1].nequa;m++){
                pid5=parm[pid1].equa[m].pid; strcpy(type1,parm[pid5].atomtype);
                for (int n=0;n<parm[pid2].nequa;n++){
                    pid6=parm[pid2].equa[n].pid; strcpy(type2,parm[pid6].atomtype);
                    for (int p=0;p<parm[pid3].nequa;p++){
                        pid7=parm[pid3].equa[p].pid; strcpy(type3,parm[pid7].atomtype);
                        for (int qi=0;qi<parm[pid4].nequa;qi++){
                            if (m==0&&n==0&&p==0&&qi==0) continue;
                            pid8=parm[pid4].equa[qi].pid; strcpy(type4,parm[pid8].atomtype);
                            for (int kk=0;kk<ipm_lim;kk++){
                                if (!improperparm[kk].numX) continue;
                                if (strcmp(improperparm[kk].name1,type1)!=0&&strcmp(improperparm[kk].name1,"X")!=0) continue;
                                if (strcmp(improperparm[kk].name2,type2)!=0&&strcmp(improperparm[kk].name2,"X")!=0) continue;
                                if (strcmp(improperparm[kk].name3,type3)!=0&&strcmp(improperparm[kk].name3,"X")!=0) continue;
                                if (strcmp(improperparm[kk].name4,type4)!=0&&strcmp(improperparm[kk].name4,"X")!=0) continue;
                                score1=(!strcmp(improperparm[kk].name1,"X"))?wt.X:0;
                                score2=(!strcmp(improperparm[kk].name2,"X"))?wt.X:0;
                                score3=(!strcmp(improperparm[kk].name3,"X"))?wt.X3:0;
                                score4=(!strcmp(improperparm[kk].name4,"X"))?wt.X:0;
                                score=score1+score2+score3+score4;
                                if (score<bestscore){ bestimproperid=kk; bestscore=score; suc=1; }
                                break;
                            }
                        }
                    }
                }
            }
        }
        if (suc&&bestimproperid>=0) goto improper_output;

/*	Step 5 considering corresponding atom types for specific improper parameters	*/
        if (!suc){
            bestscore=INITSCORE; bestimproperid=-1;
            for (int m=0;m<parm[pid1].ncorr;m++){
                pid5=parm[pid1].corr[m].pid; strcpy(type1,parm[pid5].atomtype);
                score1=parm[pid1].corr[m].improper;
                for (int n=0;n<parm[pid2].ncorr;n++){
                    pid6=parm[pid2].corr[n].pid; strcpy(type2,parm[pid6].atomtype);
                    score2=parm[pid2].corr[n].improper;
                    for (int p=0;p<parm[pid3].ncorr;p++){
                        pid7=parm[pid3].corr[p].pid; strcpy(type3,parm[pid7].atomtype);
                        score3=parm[pid3].corr[p].improper*wt.IMPROPER;
                        for (int qi=0;qi<parm[pid4].ncorr;qi++){
                            if (parm[pid1].corr[m].type<=1&&parm[pid2].corr[n].type<=1&&
                                parm[pid3].corr[p].type<=1&&parm[pid4].corr[qi].type<=1) continue;
                            pid8=parm[pid4].corr[qi].pid; strcpy(type4,parm[pid8].atomtype);
                            score4=parm[pid4].corr[qi].improper;
                            score=score1+score2+score3+score4+wt.GROUP;
                            if (parm[pid5].group==parm[pid6].group&&
                                parm[pid5].group==parm[pid7].group&&
                                parm[pid5].group==parm[pid8].group) score-=wt.GROUP;
                            if (score<bestscore)
                                for (int kk=0;kk<ipm_lim;kk++){
                                    if (improperparm[kk].numX>0) continue;
                                    if (!strcmp(improperparm[kk].name1,type1)&&
                                        !strcmp(improperparm[kk].name2,type2)&&
                                        !strcmp(improperparm[kk].name3,type3)&&
                                        !strcmp(improperparm[kk].name4,type4))
                                        { bestimproperid=kk; bestscore=score; suc=1; break; }
                                }
                        }
                    }
                }
            }
        }
        if (suc&&bestimproperid>=0) goto improper_output;

/*	Step 6  considering corresponding atom types for general improper parameters	*/
        if (!suc){
            bestscore=INITSCORE; bestimproperid=-1;
            for (int m=0;m<parm[pid1].ncorr;m++){
                pid5=parm[pid1].corr[m].pid; strcpy(type1,parm[pid5].atomtype);
                score1=parm[pid1].corr[m].improper;
                for (int n=0;n<parm[pid2].ncorr;n++){
                    pid6=parm[pid2].corr[n].pid; strcpy(type2,parm[pid6].atomtype);
                    score2=parm[pid2].corr[n].improper;
                    for (int p=0;p<parm[pid3].ncorr;p++){
                        pid7=parm[pid3].corr[p].pid; strcpy(type3,parm[pid7].atomtype);
                        score3=parm[pid3].corr[p].improper+wt.IMPROPER;
                        for (int qi=0;qi<parm[pid4].ncorr;qi++){
                            if (parm[pid1].corr[m].type<=1&&parm[pid2].corr[n].type<=1&&
                                parm[pid3].corr[p].type<=1&&parm[pid4].corr[qi].type<=1) continue;
                            pid8=parm[pid4].corr[qi].pid; strcpy(type4,parm[pid8].atomtype);
                            score4=parm[pid4].corr[qi].improper;
                            score=score1+score2+score3+score4+wt.GROUP;
                            if (parm[pid5].group==parm[pid6].group&&
                                parm[pid5].group==parm[pid7].group&&
                                parm[pid5].group==parm[pid8].group) score-=wt.GROUP;
                            for (int kk=0;kk<ipm_lim;kk++){
                                if (!improperparm[kk].numX) continue;
                                if (strcmp(improperparm[kk].name1,type1)!=0&&strcmp(improperparm[kk].name1,"X")!=0) continue;
                                if (strcmp(improperparm[kk].name2,type2)!=0&&strcmp(improperparm[kk].name2,"X")!=0) continue;
                                if (strcmp(improperparm[kk].name3,type3)!=0&&strcmp(improperparm[kk].name3,"X")!=0) continue;
                                if (strcmp(improperparm[kk].name4,type4)!=0&&strcmp(improperparm[kk].name4,"X")!=0) continue;
                                if (!strcmp(improperparm[kk].name1,"X")) score+=wt.X;
                                if (!strcmp(improperparm[kk].name2,"X")) score+=wt.X;
                                if (!strcmp(improperparm[kk].name3,"X")) score+=wt.X3;
                                if (!strcmp(improperparm[kk].name4,"X")) score+=wt.X;
                                if (score<bestscore){ bestimproperid=kk; bestscore=score; suc=1; }
                                break;
                            }
                        }
                    }
                }
            }
        }
        if (suc&&bestimproperid>=0) goto improper_output;

        fprintf(fpout,"%-2s-%-2s-%-2s-%-2s %11.1lf%15.1lf%12.1lf"
                "          Using the default value\n",
                name1,name2,name3,name4,1.1,180.0,2.0);
        strcpy(improperparm[improperparmnum].name1,name1);
        strcpy(improperparm[improperparmnum].name2,name2);
        strcpy(improperparm[improperparmnum].name3,name3);
        strcpy(improperparm[improperparmnum].name4,name4);
        improperparm[improperparmnum].phase=180.0;
        improperparm[improperparmnum].fterm=2.0;
        improperparm[improperparmnum].force=1.1;
        improperparmnum++;
        GROW_ARRAY(improperparm,improperparmnum,maximproperparm,MAX_FF_IMPROPER,IMPROPER);
        continue;

    improper_output:
        if (output_improper_flag||allparm_flag) {
            int impr_general= (
                    !strcmp(improperparm[bestimproperid].name1,"X") ||
                    !strcmp(improperparm[bestimproperid].name2,"X") ||
                    !strcmp(improperparm[bestimproperid].name3,"X") ||
                    !strcmp(improperparm[bestimproperid].name4,"X") );
            fprintf(fpout,"%-2s-%-2s-%-2s-%-2s %11.1lf%15.1lf%12.1lf          "
                          "Same as %-2s-%-2s-%-2s-%-2s, penalty score=%5.1lf%s\n",
                    name1,name2,name3,name4,
                    improperparm[bestimproperid].force,
                    improperparm[bestimproperid].phase,
                    improperparm[bestimproperid].fterm,
                    improperparm[bestimproperid].name1,
                    improperparm[bestimproperid].name2,
                    improperparm[bestimproperid].name3,
                    improperparm[bestimproperid].name4,
                    bestscore, impr_general ? " (general term)":"");
        }
        strcpy(improperparm[improperparmnum].name1,name1);
        strcpy(improperparm[improperparmnum].name2,name2);
        strcpy(improperparm[improperparmnum].name3,name3);
        strcpy(improperparm[improperparmnum].name4,name4);
        improperparm[improperparmnum].phase=improperparm[bestimproperid].phase;
        improperparm[improperparmnum].fterm=improperparm[bestimproperid].fterm;
        improperparm[improperparmnum].force=improperparm[bestimproperid].force;
        improperparmnum++;
        GROW_ARRAY(improperparm,improperparmnum,maximproperparm,MAX_FF_IMPROPER,IMPROPER);
        bestimproperid=-1;
    }
}

char *assemble_command_line(int argc, char *argv[]) {
    const size_t max_len = 1024; // must be <= leap MAXSTRINGLENGTH
    char *prefix = "ParmChk2("  ") ";

    if (argc <= 0 || argv == NULL) return NULL;

    size_t prefix_len = strlen(prefix);
    size_t total_size = prefix_len;

    // 1. Calculate required size with overflow check
    for (int i = 1; i < argc; i++) {
        size_t arg_len = strlen(argv[i]);
        // Add length of argument + 1 for space or null terminator + 2 for quotes
        if (arg_len > max_len - total_size - 1) return NULL;
        total_size += arg_len + 3;
    }

    // Allocate memory for the combined string
    char *result = malloc(total_size);
    if (!result) return NULL; // Allocation failed

    strcpy(result, prefix);
    size_t offset = prefix_len;

    // 2. Concatenate safely
    for (int i = 1; i < argc; i++) {
        bool need_quotes = strchr( argv[i], ' ');
        int written = snprintf(result + offset, total_size - offset, "%s%s%s%s",
                               need_quotes ? "\"" : "",
                               argv[i],
                               need_quotes ? "\"" : "",
                               (i < argc - 1) ? " " : "");

        if (written < 0 || (size_t)written >= total_size - offset) {
            free(result);
            return NULL;
        }
        offset += (size_t)written;
    }

    return result;
}


/*========================================================================
 * main
 *======================================================================*/
typedef enum {
    OPT_I, OPT_O, OPT_F, OPT_S, OPT_P, OPT_PF, OPT_C,
    OPT_A, OPT_W, OPT_FC, OPT_ATC, OPT_ATT, OPT_FRC, OPT_AFRC,
    OPT_V, OPT_ICOR, OPT_MULTI, OPT_UNKNOWN=-1
} OptKey;

typedef struct { const char *flag; OptKey key; } OptEntry;
static const OptEntry opt_table[]={
    {"-i",OPT_I},{"-o",OPT_O},{"-f",OPT_F},{"-s",OPT_S},{"-p",OPT_P},
    {"-pf",OPT_PF},{"-c",OPT_C},{"-a",OPT_A},{"-w",OPT_W},{"-fc",OPT_FC},
    {"-atc",OPT_ATC},{"-att",OPT_ATT},{"-frc",OPT_FRC},{"-afrc",OPT_AFRC},
    {"-v",OPT_V},{"-icor",OPT_ICOR},{"-multi",OPT_MULTI},
    {NULL,OPT_UNKNOWN}
};
static OptKey find_opt(const char *arg) {
    for (const OptEntry *e=opt_table;e->flag;e++)
        if (!strcmp(arg,e->flag)) return e->key;
    return OPT_UNKNOWN;
}

static const struct { const char *key; const char *file; } frcmod_table[]={
    {"ff14SB","frcmod.ff14SB"},{"ff99SB","frcmod.ff99SB"},{"ff03","frcmod.ff03"},
    {"bsc1","frcmod.parmbsc1"},{"ol15","frcmod.DNA.OL15"},
    {"ol3","frcmod.RNA.OL3"},{"yil","frcmod.parmCHI_YIL"},
    {NULL,NULL}
};

int main(int argc, char *argv[]) {
    int overflow_flag=0;
    char frcmod_filename[MAXCHAR];

    default_cinfo(&cinfo); default_minfo(&minfo);
    amberhome=(char*)getenv("AMBERHOME");
    if (!amberhome){ fprintf(stdout,"AMBERHOME is not set!\n"); exit(1); }
    minfo.connect_file[0]='\0';
    build_dat_path(minfo.connect_file,"CONNECT.TPL",sizeof minfo.connect_file,0);

    if (argc==2&&(!strcmp(argv[1],"-l")||!strcmp(argv[1],"-L"))){
        printf("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
        printf("Format of Additional PARMCHK ATCS file\n");
        printf("# Each atom type starts with \"PARM\" in the first four columns.\n");
        printf("# equivalent_flag: 1 for cc/ce/cg/nc/ne/pc/pe\n");
        printf("#                  2 for cd/cf/ch/nd/nf/pd/pf\n");
        printf("#                  0 for others\n");
        printf("# Example: np defined as equivalent to n3\n");
        printf("---------------------------------------------------------------------------\n");
        printf("#   atomtype  improper_flag  group_id  mass  equivalent_flag  atomic_num\n");
        printf("PARM    np      1               1       14.01   0       7\n");
        printf("EQUA    n3\n");
        printf("---------------------------------------------------------------------------\n");
        exit(0);
    }
    if (argc<7||argc%2==0){ print_help(); exit(1); }

    for (int ii=1;ii<argc;ii+=2) {
        if (ii+1>=argc){ fprintf(stderr,"Missing value for %s\n",argv[ii]); exit(1); }
        const char *val=argv[ii+1];
        switch (find_opt(argv[ii])){
            case OPT_I:    strcpy(ifilename,val); break;
            case OPT_O:    strcpy(ofilename,val); break;
            case OPT_P:    strcpy(pfilename,val); break;
            case OPT_PF:   pformat=atoi(val); break;
            case OPT_C:    strcpy(cfilename,val); break;
            case OPT_FC:   fc_opt=atoi(val); break;
            case OPT_ATT:  attn_opt=atoi(val); break;
            case OPT_ATC:  strcpy(atcfilename,val); break;
            case OPT_FRC:  strcpy(frcmod_str,val); break;
            case OPT_AFRC: strcpy(additional_frcmod_file,val); break;
            case OPT_A:    allparm_flag=(val[0]=='Y'||val[0]=='y'); break;
            case OPT_W:    output_improper_flag=(val[0]=='Y'||val[0]=='y'); break;
            case OPT_F:    iformat=parse_iformat(val); break;
            case OPT_S:    ffset=parse_ffset(val); break;
            case OPT_V:    verbose=atoi(val); break;
            case OPT_ICOR: equa_corr_inherit=(val[0]=='Y'||val[0]=='y'); break;
            case OPT_MULTI: maxbest=atoi(val);
                            if (maxbest<1) maxbest=1;
                            if (maxbest>MAXBEST) maxbest=MAXBEST;
                            break;
            case OPT_UNKNOWN:
                if (!strcmp(argv[ii],"-at")){
                    fprintf(stderr,"Error: -at does not exist;"
                            " did you mean -atc or -att?\n"); exit(1);
                }
                fprintf(stderr,"Unknown option: %s\n",argv[ii]); exit(1);
        }
    }

    if (!pfilename[0]){
        gaff_set=1;
        const char *dat;
        switch (ffset){
            case FFTYPE_GAFF:    dat="gaff.dat";      break;
            case FFTYPE_GAFF2:   dat="gaff2.dat"; gaff_set=2; break;
            case FFTYPE_PARM99:  dat="parm99.dat";    break;
            case FFTYPE_PARM10:  dat="parm10.dat";    break;
            case FFTYPE_LIPID14: dat="lipid14.dat";   break;
            case FFTYPE_GLYCAM:  dat="GLYCAM06j.dat"; break;
            default:             dat="gaff.dat";      break;
        }
        build_path(pfilename,"/dat/leap/parm/",dat,sizeof pfilename,0);
        pformat=1;
    }
    if (pformat!=1&&pformat!=2) pformat=1;
    if (attn_opt!=1&&attn_opt!=2) attn_opt=1;
    if (fc_opt!=1&&fc_opt!=2) fc_opt=1;

    if (!cfilename[0]){
        build_dat_path(cfilename,"PARMCHK.DAT",sizeof cfilename,0);
        if (access(cfilename,R_OK)){
            if (!access("PARMCHK.DAT",R_OK)) strcpy(cfilename,"PARMCHK.DAT");
            else { fprintf(stderr,"PARMCHK.DAT not found; exiting.\n"); exit(1); }
        }
    }

    maxparmnum=MAXPARM; maxatomtype=MAX_FF_ATOMTYPE; maxvdwparm=MAX_FF_VDW;
    maxbondparm=MAX_FF_BOND; maxangleparm=MAX_FF_ANGLE;
    maxtorsionparm=MAX_FF_TORSION; maximproperparm=MAX_FF_IMPROPER;

    atom        =(ATOM*)      malloc(sizeof(ATOM)     *cinfo.maxatom);
    bond_array  =(BOND*)      malloc(sizeof(BOND)     *cinfo.maxbond);
    parm        =(PARM*)      calloc(maxparmnum,        sizeof(PARM));
    parmid      =(int*)       malloc(sizeof(int)       *cinfo.maxatom);
    improper    =(IMPROPERID*)malloc(sizeof(IMPROPERID)*cinfo.maxatom);
    atomtype    =(ATOMTYPE*)  calloc(maxatomtype,       sizeof(ATOMTYPE));
    bondparm    =(BOND_FF*)   calloc(maxbondparm,       sizeof(BOND_FF));
    angleparm   =(ANGLE*)     calloc(maxangleparm,      sizeof(ANGLE));
    torsionparm =(TORSION*)   calloc(maxtorsionparm,    sizeof(TORSION));
    improperparm=(IMPROPER*)  calloc(maximproperparm,   sizeof(IMPROPER));
    vdwparm     =(VDW*)       calloc(maxvdwparm,        sizeof(VDW));

    if (!atom||!bond_array||!parm||!parmid||!improper||
        !atomtype||!bondparm||!angleparm||!torsionparm||!improperparm||!vdwparm){
        fprintf(stdout,"memory allocation error\n"); exit(1);
    }
    for (int ii=0;ii<cinfo.maxbond;ii++) bond_array[ii].jflag=-1;

    switch (iformat){
        case FILETYPE_PREPI:
            overflow_flag=rprepi(ifilename,&atomnum,atom,&bondnum,bond_array,&cinfo,&minfo); break;
        case FILETYPE_PREPC:
            overflow_flag=rprepc(ifilename,&atomnum,atom,&cinfo,&minfo); break;
        case FILETYPE_AC:
            overflow_flag=rac(ifilename,&atomnum,atom,&bondnum,bond_array,&cinfo,&minfo); break;
        case FILETYPE_MOL2:
            overflow_flag=rmol2(ifilename,&atomnum,atom,&bondnum,bond_array,&cinfo,&minfo,1); break;
        case FILETYPE_FRCMOD:
        case FILETYPE_LEAPLOG:
            if (iformat==FILETYPE_FRCMOD) rfrcmod(ifilename,0);
            else                          rleaplog(ifilename,0);
            r_atomtype    =(ATOMTYPE*) calloc(r_atomtype_num,    sizeof(ATOMTYPE));
            r_bondparm    =(BOND_FF*)  calloc(r_bondparm_num,    sizeof(BOND_FF));
            r_angleparm   =(ANGLE*)    calloc(r_angleparm_num,   sizeof(ANGLE));
            r_torsionparm =(TORSION*)  calloc(r_torsionparm_num, sizeof(TORSION));
            r_improperparm=(IMPROPER*) calloc(r_improperparm_num,sizeof(IMPROPER));
            r_vdwparm     =(VDW*)      calloc(r_vdwparm_num,     sizeof(VDW));
            if (!r_atomtype||!r_bondparm||!r_angleparm||
                !r_torsionparm||!r_improperparm||!r_vdwparm){
                fprintf(stdout,"memory allocation error for r_* arrays\n"); exit(1);
            }
            r_atomtype_num=r_bondparm_num=r_angleparm_num=0;
            r_torsionparm_num=r_improperparm_num=r_vdwparm_num=0;
            if (iformat==FILETYPE_FRCMOD) rfrcmod(ifilename,1);
            else                          rleaplog(ifilename,1);
            break;
        default:
            fprintf(stderr,"Error: invalid input file format\n"); exit(1);
    }

    if (overflow_flag){
        cinfo.maxatom=atomnum+10; cinfo.maxbond=bondnum+4*atomnum+10;
        free(atom); free(bond_array);
        atom      =(ATOM*)malloc(sizeof(ATOM)*cinfo.maxatom);
        bond_array=(BOND*)malloc(sizeof(BOND)*cinfo.maxbond);
        if (!atom||!bond_array){ fprintf(stdout,"memory allocation error\n"); exit(1); }
        for (int ii=0;ii<cinfo.maxbond;ii++) bond_array[ii].jflag=-1;
        if (iformat==FILETYPE_PREPI)
            overflow_flag=rprepi(ifilename,&atomnum,atom,&bondnum,bond_array,&cinfo,&minfo);
        if (iformat==FILETYPE_PREPC)
            overflow_flag=rprepc(ifilename,&atomnum,atom,&cinfo,&minfo);
        if (iformat==FILETYPE_AC)
            overflow_flag=rac(ifilename,&atomnum,atom,&bondnum,bond_array,&cinfo,&minfo);
        if (iformat==FILETYPE_MOL2)
            overflow_flag=rmol2(ifilename,&atomnum,atom,&bondnum,bond_array,&cinfo,&minfo,1);
    }

    if (iformat==FILETYPE_PREPI){ atomicnum(atomnum,atom); adjustatomname(atomnum,atom,1); }
    if (iformat==FILETYPE_PREPC){
        atomicnum(atomnum,atom);
        overflow_flag=connect(minfo.connect_file,atomnum,atom,&bondnum,bond_array,cinfo.maxbond);
        adjustatomname(atomnum,atom,1);
    }
    if (iformat==FILETYPE_AC)   atomicnum(atomnum,atom);
    if (iformat==FILETYPE_MOL2){ atomicnum(atomnum,atom); default_inf(atomnum,atom,0); }

    read_parmchk_parm(cfilename, 1);
    parmnum2 = parmnum;
    if (atcfilename[0]) read_parmchk_parm(atcfilename, 0);
    assign_parmid();

    if (equa_corr_inherit){
        VLOG(1,"EQUA-CORR inheritance enabled (WEIGHT_EQUTYPE=%.2f)\n",wt.EQUTYPE);
        build_ext_corr();
    }

    if (frcmod_str[0])
        for (int fi=0;frcmod_table[fi].key;fi++)
            if (strstr(frcmod_str,frcmod_table[fi].key)){
                build_path(frcmod_filename,"/dat/leap/parm/",
                           frcmod_table[fi].file,sizeof frcmod_filename,0);
                readfrcmod(frcmod_filename);
                update_last_counts();
            }

    if (additional_frcmod_file[0]){
        afrcmod_atomtypenum_beg=atomtypenum; afrcmod_bondparmnum_beg=bondparmnum;
        afrcmod_angleparmnum_beg=angleparmnum; afrcmod_torsionparmnum_beg=torsionparmnum;
        afrcmod_improperparmnum_beg=improperparmnum; afrcmod_vdwparmnum_beg=vdwparmnum;
        readfrcmod(additional_frcmod_file);
        afrcmod_atomtypenum_end=atomtypenum; afrcmod_bondparmnum_end=bondparmnum;
        afrcmod_angleparmnum_end=angleparmnum; afrcmod_torsionparmnum_end=torsionparmnum;
        afrcmod_improperparmnum_end=improperparmnum; afrcmod_vdwparmnum_end=vdwparmnum;
    }
    update_last_counts();

    if (pformat==1) readparm(pfilename);
    if (pformat==2) readfrcmod(pfilename);

    if (iformat==FILETYPE_PREPI||iformat==FILETYPE_PREPC) improper_id1(ifilename);
    if (iformat==FILETYPE_AC   ||iformat==FILETYPE_MOL2)  improper_id2();

    atomtypenum2=atomtypenum; vdwparmnum2=vdwparmnum;
    bondparmnum2=bondparmnum; angleparmnum2=angleparmnum;
    torsionparmnum2=torsionparmnum; improperparmnum2=impropernum;

    if ((fpout=fopen(ofilename,"w"))==NULL){
        fprintf(stdout,"Cannot open %s for writing, exit\n",ofilename); exit(1);
    }

    {
        BondQuery    *bq; int nb=collect_bond_queries(&bq);
        AngleQuery   *aq; int na=collect_angle_queries(&aq);
        TorsionQuery *tq; int nt=collect_torsion_queries(&tq);
        ImproperQuery*iq; int ni=collect_improper_queries(&iq);

        if (iformat>FILETYPE_MOL2) allparm_flag=1;
        char *remark = assemble_command_line(argc, argv);
        if (!remark) remark = "Remark line goes here";
        fprintf(fpout,"%s\n",remark);
        chk_atomtype();
        chk_bond_core(bq,nb); chk_angle_core(aq,na);
        chk_torsion_core(tq,nt); chk_improper_core(iq,ni);
        chk_vdw();
        free(bq); free(aq); free(tq); free(iq);
    }

    fclose(fpout);
    if (allparm_flag&&iformat<=FILETYPE_MOL2) cleanup_frcmod(ofilename);
    return 0;
}
