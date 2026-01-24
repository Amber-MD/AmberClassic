/*
************************************************************************
*           All Copyright Reserved!                                    *
*                                                                      *
*  Prog:    parmchk                                                    *
*  Version: version 1.0                                                *
*  Author:  Junmei Wang                                                *
*                                                                      *
*  Department of Pharmaceutical Chemistry                              *
*  School of Pharmacy                                                  *
*  University of California                                            *
*  San Francisco   CA 94143                                            *
*  Octomber, 2001                                                      *
************************************************************************
*/
char *amberhome;
# include "common.h"
# include "define.h"
# include "atom.h"
# include "utility.c"
# include "common.c"
# include "rotate.c"
# include "ac.c"
# include "mol2.c"
# include "prep.c"

#define MAXPARM 250   /*maximum atom types */
#define MAXEQUA 25    /*maximum equal atom types for each atom type */
#define MAXCORR 50    /*maximum corresponding atom types for each atom type */
#define MAX_FF_ATOMTYPE 250
#define MAX_FF_VDW 250
#define MAX_FF_BOND 2000
#define MAX_FF_ANGLE 5000
#define MAX_FF_TORSION 1500
#define MAX_FF_IMPROPER 500
#define MAX_EQU_VDW 20
#define MAX_EQU_VDW_NUM 50
#define INITSCORE       999999
#define debug 0

typedef struct {
	int atid1;
	int atid2;
	int atid3;
	int atid4;
} IMPROPERID;

typedef struct {
	char atomtype[10];
	int  pid;
} EQUA;

typedef struct {
	char atomtype[10];
	int  pid;
	int  type;
	double bl;  	/*penalty of bond length */
	double blf; 	/*penalty of bond stretching force constant */
	double ba;  	/*penalty of bond angle */
	double baf; 	/*penalty of bond bending force constant */
	double cba; 	/*penalty of bond angle when replace angle central atom */
	double cbaf;	/*penalty of bond bending force constant of replacing angle central atom */
	double ctor;	/*penalty of torsional angle when replace central atoms */
	double tor; 	/*penalty of torsional angle when replace torsional atoms */
	double ps;      /*overall penalty score*/
	double improper;/* penalty of improper angle when replace improper atoms */
} CORR;

typedef struct {
	char atomtype[10];
	int  group;
	int  equtype;  /*equivalent atom type, cc/ce/cg/nc/ne/pc/pe belong to TYPE 1, cd/cf/ch/nd/nf/pd/pf belong to TYPE 2 and others are 0*/
	int  improper;
	int  ncorr;
	int  nequa;
	int  atomicnum;
	double mass;
	CORR corr[MAXCORR];
	EQUA equa[MAXEQUA];
} PARM;

typedef struct {
	char name[10];
	double mass;
	double pol;
	int attn; 
} ATOMTYPE;

typedef struct {
	char name1[10];
	char name2[10];
	char name3[10];
	double angle;
	double force;
	int attn; 
} ANGLE;

typedef struct {
	char name1[10];
	char name2[10];
	char name3[10];
	char name4[10];
	int mul;
	double fterm;
	double phase;
	double force;
	int attn; 
} TORSION;

typedef struct {
	char name1[10];
	char name2[10];
	char name3[10];
	char name4[10];
	int mul;
	int numX;
	double fterm;
	double phase;
	double force;
	int attn; 
} IMPROPER;

typedef struct {
	char name[10];
	double radius;
	double pot;
	int attn; 
} VDW;

typedef struct {
	char name[MAX_EQU_VDW_NUM][10];
	int num;
} EQU_VDW;

typedef struct {
	double BL;
	double BLF;
	double BA;
	double BAF;
	double X;
	double X3;
	double BA_CTR;
	double TOR_CTR;
	double IMPROPER;
	double GROUP;
	double EQUTYPE;
} WT;

typedef struct {
	double BL;
	double BLF;
	double BA;
	double BAF;
	double BA_CTR;
	double BAF_CTR;
	double TOR;
	double TOR_CTR;
	double FRACT1;
	double FRACT2;
} DEFAULT_VALUE;

MOLINFO minfo;
CONTROLINFO cinfo;
ATOM *atom;
BOND *bond_array;
WT   wt;
DEFAULT_VALUE dv;

int atomnum = 0;
int bondnum = 0;
IMPROPERID *improper;
char ifilename[MAXCHAR];
char ofilename[MAXCHAR];
char pfilename[MAXCHAR];
char cfilename[MAXCHAR];
char atcfilename[MAXCHAR];
char line[MAXCHAR];
char frcmod_str[MAXCHAR]="";
int  ifrcmod_str = 0;
int  ffset = 0;
int  gaff_set = 0;
int  ipfilename = 0;
int  cindex = 0;
int  iatc = 0;
int impropernum = 0;
int output_improper_flag = 1;
int attn_opt = 1;
int i, j, k, l;
int maxatomtype = 0;
int maxvdwparm = 0;
int maxbondparm = 0;
int maxangleparm = 0;
int maxtorsionparm = 0;
int maximproperparm = 0;
int atomtypenum = 0;
int parmnum = 0;
int parmnum2 = 0;
int maxparmnum;
int vdwparmnum = 0;
int bondparmnum = 0;
int angleparmnum = 0;
int torsionparmnum = 0;
int improperparmnum = 0;

double THRESHOLD_BA; 
double equtype_penalty_score;
/* 
  atomtypenum2 etc are numbers of parameters from force field files
  not including those newly found in the file
*/
int atomtypenum2 = 0; 
int vdwparmnum2 = 0;
int bondparmnum2 = 0;
int angleparmnum2 = 0;
int torsionparmnum2 = 0;
int improperparmnum2 = 0;

/* for tracking ff parameters from different sources*/
int last_atomtypenum = 0; 
int last_vdwparmnum = 0;
int last_bondparmnum = 0;
int last_angleparmnum = 0;
int last_torsionparmnum = 0;
int last_improperparmnum = 0;

int afrcmod_atomtypenum_beg ;
int afrcmod_bondparmnum_beg ;
int afrcmod_angleparmnum_beg ;
int afrcmod_torsionparmnum_beg ;
int afrcmod_improperparmnum_beg ;
int afrcmod_vdwparmnum_beg ;

int afrcmod_atomtypenum_end ;
int afrcmod_bondparmnum_end ;
int afrcmod_angleparmnum_end ;
int afrcmod_torsionparmnum_end ;
int afrcmod_improperparmnum_end ;
int afrcmod_vdwparmnum_end ;

/*H-1, C-2, N-3, O-4, F-5, Cl-6, Br-7, I-8, S-9, P-10*/

int *parmid; 
PARM *parm;
ATOMTYPE *atomtype;
BOND_FF *bondparm;
ANGLE *angleparm;
TORSION *torsionparm;
IMPROPER *improperparm;
VDW *vdwparm;
/* for equivalent vdw types */
EQU_VDW equ_vdw[MAX_EQU_VDW];
int equ_vdw_num = 0;

/* for frcmod and leaplog formats*/
ATOMTYPE *r_atomtype;
BOND_FF *r_bondparm;
ANGLE *r_angleparm;
TORSION *r_torsionparm;
IMPROPER *r_improperparm;
VDW *r_vdwparm;
int r_atomtype_num = 0;
int r_bondparm_num = 0;
int r_angleparm_num = 0;
int r_torsionparm_num = 0;
int r_improperparm_num = 0;
int r_vdwparm_num = 0;

FILE *fp, *fpout;

int allparm_flag = 0;
int pformat = 1;
int iformat = -1;
int fc_opt = 1;

ANGLE   bestba;
int	bestblid = -1;
int 	bestbaid = -1;
int 	besttorid= -1;
int 	bestimproperid=-1;
double bestscore;
/* additional frcmod file*/
int iadditional_frcmod = 0;
char additional_frcmod_file[MAXCHAR];

/* reading parameters for calculating bl and ba parameters*/
int iread_blbaparm = 0;
int nblf_parm = 0;
char blba_parmfile[MAXCHAR];

/*
parameter sets
ff03		parm99.dat + frcmod.ff03
ff99SB		parm99.dat + frcmod.ff99SB
ff14SB		parm10.dat + frcmod.ff14SB
lipid14 	lipid14.dat
DNA.bsc1 	parm10.dat frcmod.parmbsc1
DNA.OL15 	parm10.dat frcmod.DNA.OL15
RNA.OL3  	parm10.dat 
RNA.YIL  	parm10.dat frcmod.parmCHI_YIL
the following parameter sets are not supported due to duplicated atom type definitions
glycam_06EP 	GLYCAM_06EPb.dat 
glycam_06   	GLYCAM_06j.dat  
*/

/* set up improper atoms for prepi/prepc file */
void improper_id1(char *filename) {
	FILE *fpin;
	int i;
	int readindex;
	char tmpchar1[MAXCHAR];
	char tmpchar2[MAXCHAR];
	char tmpchar3[MAXCHAR];
	char tmpchar4[MAXCHAR];

	if ((fpin = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open file %s in improper_id1(), exit\n", filename);
		exit(1);
	}
	readindex = 0;
	for (;;) {
		if (fgets(line, MAXCHAR, fpin) == NULL)
			break;
		sscanf(line, "%s", tmpchar1);
		if (strcmp("IMPROPER", tmpchar1) == 0) {
			readindex = 1;
			continue;
		}
		if (spaceline(line) == 1 && readindex == 1)
			readindex = 0;
		if (strcmp(tmpchar1, "DONE") == 0 && readindex == 1)
			readindex = 0;
		if (strcmp(tmpchar1, "STOP") == 0 && readindex == 1)
			readindex = 0;
		if (strcmp(tmpchar1, "CHARGE") == 0 && readindex == 1)
			readindex = 0;
		if (strcmp(tmpchar1, "LOOP") == 0 && readindex == 1)
			readindex = 0;
		if (strcmp(tmpchar1, "IMPROPER") == 0 && readindex == 1)
			readindex = 0;

		if (readindex == 1) {
			sscanf(line, "%s%s%s%s", tmpchar1, tmpchar2, tmpchar3,
				   tmpchar4);
			if (strcmp(tmpchar1, "-M") == 0)
				continue;
			if (strcmp(tmpchar2, "-M") == 0)
				continue;
			if (strcmp(tmpchar3, "-M") == 0)
				continue;
			if (strcmp(tmpchar4, "-M") == 0)
				continue;
			if (strcmp(tmpchar1, "+M") == 0)
				continue;
			if (strcmp(tmpchar2, "+M") == 0)
				continue;
			if (strcmp(tmpchar3, "+M") == 0)
				continue;
			if (strcmp(tmpchar4, "+M") == 0)
				continue;
			for (i = 0; i < atomnum; i++) {
				if (strcmp(tmpchar1, atom[i].name) == 0) {
					strcpy(tmpchar1, atom[i].ambername);
					improper[impropernum].atid1 = i;
					continue;
				}
				if (strcmp(tmpchar2, atom[i].name) == 0) {
					strcpy(tmpchar2, atom[i].ambername);
					improper[impropernum].atid2 = i;
					continue;
				}
				if (strcmp(tmpchar3, atom[i].name) == 0) {
					strcpy(tmpchar3, atom[i].ambername);
					improper[impropernum].atid3 = i;
					continue;
				}
				if (strcmp(tmpchar4, atom[i].name) == 0) {
					strcpy(tmpchar4, atom[i].ambername);
					improper[impropernum].atid4 = i;
					continue;
				}
			}
			impropernum++;
		}
	}
	fclose(fpin);
}

/* set up improper atoms for ac or mol2 file */
void improper_id2() {
	int i;
	impropernum = 0;
	for (i = 0; i < atomnum; i++) 
		if (atom[i].improper == 1) {
			improper[impropernum].atid3 = i;
			improper[impropernum].atid1 = atom[i].con[0];
			improper[impropernum].atid2 = atom[i].con[1];
			improper[impropernum].atid4 = atom[i].con[2];
                        if(atom[i].con[0] <0 || atom[i].con[1] <0 || atom[i].con[2] <0)
                                continue;
			impropernum++;
		}
}

void rleaplog(char *filename, int mode) {
int rindex = 0;
int suc = 0;
int i;
char tmpc1[MAXCHAR];
char tmpc2[MAXCHAR];
char tmpc3[MAXCHAR];
char tmpc4[MAXCHAR];
char tmpc5[MAXCHAR];
char tmpc6[MAXCHAR];
char tmpc7[MAXCHAR];
char tmpc8[MAXCHAR];
char tmpc9[MAXCHAR];
char tmpc10[MAXCHAR];
char tmpc11[MAXCHAR];
char tmpc12[MAXCHAR];
char tmpc13[MAXCHAR];
char tmpc14[MAXCHAR];
	if ((fp = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open file %s in rleaplog(), exit\n", filename);
		exit(1);
	}
/*
Building topology.
Building atom parameters.
For atom: .R<YCM 281>.A<CB 4> Could not find vdW (or other) parameters for type: Cx
Building bond parameters.
Could not find bond parameter for: S - H1
Could not find bond parameter for: S - H1
Building angle parameters.
Could not find angle parameter: H1 - S - H1
Could not find angle parameter: S - S - H1
Could not find angle parameter: S - S - H1
Could not find angle parameter: S - CT - C
Could not find angle parameter: S - CT - C
Could not find angle parameter: CT - S - H1
Could not find angle parameter: CT - S - H1
Could not find angle parameter: N - CT - S
Building proper torsion parameters.
 ** No torsion terms for  CT-S-S-H1
 ** No torsion terms for  CT-S-S-H1
*/
for (;;) {
	if (fgets(line, MAXCHAR, fp) == NULL)
		break;
	sscanf(line, "%s%s%s%s%s%s%s%s%s%s%s%s%s%s", tmpc1, tmpc2, tmpc3, tmpc4, tmpc5, tmpc6, tmpc7, tmpc8, tmpc9, tmpc10,
                     tmpc11, tmpc12, tmpc13, tmpc14);		
	if(strcmp(tmpc1, "Building")  == 0 && strcmp(tmpc2, "topology.") == 0) rindex = 1;	
	if(rindex != 0) {
		if(strcmp(tmpc1, "Building")  == 0 && strcmp(tmpc2, "atom") == 0 && strcmp(tmpc3, "parameters.")  == 0) rindex = 2;
		if(strcmp(tmpc1, "Building")  == 0 && strcmp(tmpc2, "bond") == 0 && strcmp(tmpc3, "parameters.")  == 0) rindex = 3;
		if(strcmp(tmpc1, "Building")  == 0 && strcmp(tmpc2, "angle") == 0 && strcmp(tmpc3, "parameters.")  == 0) rindex = 4;
		if(strcmp(tmpc1, "Building")  == 0 && strcmp(tmpc2, "proper") == 0 && strcmp(tmpc3, "torsion")  == 0) rindex = 5;
		if(strcmp(tmpc1, "Building")  == 0 && strcmp(tmpc2, "improper") == 0 && strcmp(tmpc3, "torsion")  == 0) rindex = 6;
	}	
	if(rindex == 2)  {
		if(mode == 0) {
			r_atomtype_num ++;
			r_vdwparm_num ++;
		}
		if(mode  == 1) {
			if(strcmp(tmpc1, "For") == 0 && strcmp(tmpc2, "atom:") == 0) {
				suc = 0;
				for(i=0; i< r_atomtype_num; i++) 
					if(strcmp(r_atomtype[i].name, tmpc14) == 0) {
						suc =1;
						break;
					}
				if(suc == 0) {
					strcpy(r_atomtype[r_atomtype_num].name, tmpc14);
					r_atomtype_num ++;
					strcpy(r_atomtype[r_vdwparm_num].name, tmpc14);
					r_vdwparm_num ++;
				}
			}
		}
	}
	if(rindex == 3) {
		if(mode == 0) r_bondparm_num ++;
		if(mode == 1) {
			if(strcmp(tmpc1, "Could") == 0 && strcmp(tmpc2, "not") == 0) {
				suc = 0; 
				for(i=0; i< r_bondparm_num; i++) {
					if(strcmp(r_bondparm[i].name1, tmpc7) == 0 && strcmp(r_bondparm[i].name2, tmpc9) == 0) {
						suc =1;
						break;
					}
					if(strcmp(r_bondparm[i].name2, tmpc7) == 0 && strcmp(r_bondparm[i].name1, tmpc9) == 0) {
						suc =1;
						break;
					}
				}
				if(suc == 0) {
					strcpy(r_bondparm[r_bondparm_num].name1, tmpc7);
					strcpy(r_bondparm[r_bondparm_num].name2, tmpc9);
					r_bondparm_num ++;
				}
			}
		}
	}
	if(rindex == 4) {
		if(mode == 0) r_angleparm_num ++;
		if(mode == 1) {
			if(strcmp(tmpc1, "Could") == 0 && strcmp(tmpc2, "not") == 0) {
				suc = 0; 
				for(i=0; i< r_angleparm_num; i++) {
					if(strcmp(r_angleparm[i].name1, tmpc6) == 0 && strcmp(r_angleparm[i].name2, tmpc8) == 0 && strcmp(r_angleparm[i].name3, tmpc10) == 0 ) {
						suc =1;
						break;
					}
					if(strcmp(r_angleparm[i].name3, tmpc6) == 0 && strcmp(r_angleparm[i].name2, tmpc8) == 0 && strcmp(r_angleparm[i].name1, tmpc10) == 0 ) {
						suc =1;
						break;
					}
				}
				if(suc == 0) {
					strcpy(r_angleparm[r_angleparm_num].name1, tmpc6);
					strcpy(r_angleparm[r_angleparm_num].name2, tmpc8);
					strcpy(r_angleparm[r_angleparm_num].name3, tmpc10);
					r_angleparm_num ++;
				}
			}
		}
	}

        if(rindex == 5) {
		if(mode == 0) r_torsionparm_num ++;
		if(mode == 1) {
                	if(strcmp(tmpc2, "No") == 0 && strcmp(tmpc3, "torsion") == 0) {
				for(i=0; i<strlen(line); i++) {
					if(line[i] == '-') line[i] = ' ';
				}
				sscanf(line, "%s%s%s%s%s%s%s%s%s%s%s%s%s%s", tmpc1, tmpc2, tmpc3, tmpc4, tmpc5, tmpc6, tmpc7, tmpc8, tmpc9, tmpc10,
                     			tmpc11, tmpc12, tmpc13, tmpc14);		
			
                       		 suc = 0;
                       		 for(i=0; i< r_torsionparm_num; i++) {
                       		         if(strcmp(r_torsionparm[i].name1, tmpc6) == 0 && strcmp(r_torsionparm[i].name2, tmpc7) == 0 && 
                               		    	strcmp(r_torsionparm[i].name3, tmpc8) == 0 && strcmp(r_torsionparm[i].name4, tmpc9) ==0) {
                                       	 	suc =1;
                                        	break;
                                	}
                                	if(strcmp(r_torsionparm[i].name4, tmpc6) == 0 && strcmp(r_torsionparm[i].name3, tmpc7) == 0 && 
                                  		strcmp(r_torsionparm[i].name2, tmpc8) == 0 && strcmp(r_torsionparm[i].name1, tmpc9) ==0) {
                                        	suc =1;
                                        	break;
                                	}
                        	}
                        	if(suc == 0) {
                                	strcpy(r_torsionparm[r_torsionparm_num].name1, tmpc6);
                                	strcpy(r_torsionparm[r_torsionparm_num].name2, tmpc7);
                                	strcpy(r_torsionparm[r_torsionparm_num].name3, tmpc8);
                                	strcpy(r_torsionparm[r_torsionparm_num].name4, tmpc9);
                                	r_torsionparm_num ++;
                        	}
               	 	}
		}
        }
}
fclose(fp);
if(debug == 1 && mode == 1) {
	printf("\n-- Begin read leaplog file : %s --\n\n", filename);
        printf("\nMASS\n");
        for(i=0;i<r_atomtype_num;i++)
                printf("\n%-2s %9.4lf %9.4lf", r_atomtype[i].name, r_atomtype[i].mass,atomtype[i].pol);

        printf("\n\nBOND\n");
        for(i=0;i<r_bondparm_num;i++)
                printf("\n%-2s %-2s %9.4lf %9.4lf", r_bondparm[i].name1, r_bondparm[i].name2,
                r_bondparm[i].force, r_bondparm[i].length);

        printf("\n\nANGLE\n");
        for(i=0;i<r_angleparm_num;i++)
                printf("\n%-2s %-2s %-2s %9.4lf %9.4lf", r_angleparm[i].name1, r_angleparm[i].name2,
                r_angleparm[i].name3, r_angleparm[i].force, r_angleparm[i].angle);

        printf("\n\nTORSION\n");
        for(i=0;i<r_torsionparm_num;i++)
                printf("\n%-2s %-2s %-2s %-2s %9.4lf %9.4lf %9.4lf", r_torsionparm[i].name1,
                r_torsionparm[i].name2, r_torsionparm[i].name3, r_torsionparm[i].name4,
                r_torsionparm[i].phase, r_torsionparm[i].force, r_torsionparm[i].fterm);

        printf("\n\nIMPROPER\n");
        for(i=0;i<r_improperparm_num;i++)
                printf("\n%-2s %-2s %-2s %-2s %9.4lf %9.4lf %9.4lf", r_improperparm[i].name1,
                r_improperparm[i].name2, r_improperparm[i].name3, r_improperparm[i].name4,
                r_improperparm[i].phase, r_improperparm[i].force, r_improperparm[i].fterm);
        printf("\n\nVDW\n");
        for(i=0;i<r_vdwparm_num;i++)
                printf("\n%-2s %9.4lf %9.4lf", r_vdwparm[i].name, r_vdwparm[i].radius, r_vdwparm[i].pot);

	printf("\n-- Finished reading leaplog file : %s --\n\n", filename);
}
}


void rfrcmod (char *filename, int mode) {
FILE *fp; 
int mindex = 0;
int bindex = 0;
int aindex = 0;
int tindex = 0;
int iindex = 0;
int vindex = 0;
int pos_tor = 0;
int tmpnum  = 0;
char line[MAXCHAR];
char tmpc[MAXCHAR];
char tmpc1[MAXCHAR];
char tmpc2[MAXCHAR];
char tmpc3[MAXCHAR];
char tmpc4[MAXCHAR];
char tmpchar[MAXCHAR];
	if ((fp = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open file %s in rfrcmod(), exit\n", filename);
		exit(1);
	}
	for (;;) {
		if (fgets(line, MAXCHAR, fp) == NULL)
			break;
		if (strncmp(line, "MASS", 4) == 0) {
			mindex = 1;
			continue;
		}
		if (strncmp(line, "BOND", 4) == 0) { 
			bindex = 1;
			continue;
		}
		if (strncmp(line, "ANGL", 4) == 0){ 
			aindex = 1;
			continue;
		}
		if (strncmp(line, "DIHE", 4) == 0) { 
			pos_tor= ftell(fp);
			tindex = 1;
			continue;
		}
		if (strncmp(line, "IMPROPER", 8) == 0) { 
			iindex = 1;
			continue;
		}
		if (strncmp(line, "NONBON", 6) == 0) { 
			vindex = 1;
			continue;
		}
		if (strlen(line) <= 2) {
			if(mindex == 1) mindex = 0;
			if(bindex == 1) bindex = 0;
			if(aindex == 1) aindex = 0;
			if(tindex == 1) tindex = 0;
			if(iindex == 1) iindex = 0;
			if(vindex == 1) vindex = 0;
			continue;
		}
                if (mindex == 1) {
			if(mode == 1) {
			        r_atomtype[r_atomtype_num].attn = 0;
                        	sscanf(line, "%s%lf%lf", r_atomtype[r_atomtype_num].name,
                                   	&r_atomtype[r_atomtype_num].mass,
                                   	&r_atomtype[r_atomtype_num].pol);
				if(strstr(line, "ATTN") != NULL)  r_atomtype[r_atomtype_num].attn = 1;
			}
                       	r_atomtype_num++;
                }
                if (bindex == 1) {
			if(mode == 1) {
			        r_bondparm[r_bondparm_num].attn = 0;
                        	r_bondparm[r_bondparm_num].name1[0] = line[0];
				if(line[1]== ' ') 
                        		r_bondparm[r_bondparm_num].name1[1] = '\0';
				else {
					r_bondparm[r_bondparm_num].name1[1] = line[1];
					r_bondparm[r_bondparm_num].name1[2] = '\0';
				}
                        	r_bondparm[r_bondparm_num].name2[0] = line[3];
				if(line[4]== ' ') 
                        		r_bondparm[r_bondparm_num].name2[1] = '\0';
				else {
					r_bondparm[r_bondparm_num].name2[1] = line[4];
					r_bondparm[r_bondparm_num].name2[2] = '\0';
				}
				if(strcmp(r_bondparm[r_bondparm_num].name1, r_bondparm[r_bondparm_num].name2) > 0) {
					strcpy(tmpchar, r_bondparm[r_bondparm_num].name1);
					strcpy(r_bondparm[r_bondparm_num].name1, r_bondparm[r_bondparm_num].name2);
					strcpy(r_bondparm[r_bondparm_num].name2, tmpchar);
				}	
				
                        	sscanf(&line[5], "%lf%lf", &r_bondparm[r_bondparm_num].force,
                                   	&r_bondparm[r_bondparm_num].length);
				if(strstr(line, "ATTN")  != NULL)
			        	r_bondparm[r_bondparm_num].attn = 1;
			}
                        r_bondparm_num++;
                }
                if (aindex == 1) {
			if(mode == 1) {
		        	r_angleparm[r_angleparm_num].attn = 0;
                        	r_angleparm[r_angleparm_num].name1[0] = line[0];
				if(line[1]== ' ') 
                        		r_angleparm[r_angleparm_num].name1[1] = '\0';
				else {
					r_angleparm[r_angleparm_num].name1[1] = line[1];
					r_angleparm[r_angleparm_num].name1[2] = '\0';
				}
                        	r_angleparm[r_angleparm_num].name2[0] = line[3];
				if(line[4]== ' ') 
                        		r_angleparm[r_angleparm_num].name2[1] = '\0';
				else {
					r_angleparm[r_angleparm_num].name2[1] = line[4];
					r_angleparm[r_angleparm_num].name2[2] = '\0';
				}
                        	r_angleparm[r_angleparm_num].name3[0] = line[6];
				if(line[7]== ' ') 
                        		r_angleparm[r_angleparm_num].name3[1] = '\0';
				else {
					r_angleparm[r_angleparm_num].name3[1] = line[7];
					r_angleparm[r_angleparm_num].name3[2] = '\0';
				}
				if(strcmp(r_angleparm[r_angleparm_num].name1, r_angleparm[r_angleparm_num].name3) > 0) {
					strcpy(tmpchar, r_angleparm[r_angleparm_num].name1);
					strcpy(r_angleparm[r_angleparm_num].name1, r_angleparm[r_angleparm_num].name3);
					strcpy(r_angleparm[r_angleparm_num].name3, tmpchar);
				}	
                        	sscanf(&line[8], "%lf%lf", &r_angleparm[r_angleparm_num].force,
                                   	&r_angleparm[r_angleparm_num].angle);
				if(strstr(line, "ATTN")  != NULL)
		        		r_angleparm[r_angleparm_num].attn = 1;
			}
                        r_angleparm_num++;
                }
                if (tindex == 1) { 
/*
                        if(line[0] == 'X' || line[3] == 'X' || line[6] == 'X' || line[9] == 'X')
                                continue;
*/
			if(mode == 1) {
		       		r_torsionparm[r_torsionparm_num].attn = 0;
                        	r_torsionparm[r_torsionparm_num].name1[0] = line[0];
				if(line[1]== ' ') 
                        		r_torsionparm[r_torsionparm_num].name1[1] = '\0';
				else {
					r_torsionparm[r_torsionparm_num].name1[1] = line[1];
					r_torsionparm[r_torsionparm_num].name1[2] = '\0';
				}
                        	r_torsionparm[r_torsionparm_num].name2[0] = line[3];
				if(line[4]== ' ') 
                        		r_torsionparm[r_torsionparm_num].name2[1] = '\0';
				else {
					r_torsionparm[r_torsionparm_num].name2[1] = line[4];
					r_torsionparm[r_torsionparm_num].name2[2] = '\0';
				}
                        	r_torsionparm[r_torsionparm_num].name3[0] = line[6];
				if(line[7]== ' ') 
                        		r_torsionparm[r_torsionparm_num].name3[1] = '\0';
				else {
					r_torsionparm[r_torsionparm_num].name3[1] = line[7];
					r_torsionparm[r_torsionparm_num].name3[2] = '\0';
				}
                        	r_torsionparm[r_torsionparm_num].name4[0] = line[9];
				if(line[10]== ' ') 
                        		r_torsionparm[r_torsionparm_num].name4[1] = '\0';
				else {
					r_torsionparm[r_torsionparm_num].name4[1] = line[10];
					r_torsionparm[r_torsionparm_num].name4[2] = '\0';
				}
				if(strcmp(r_torsionparm[r_torsionparm_num].name2, r_torsionparm[r_torsionparm_num].name3) > 0) {
					strcpy(tmpchar, r_torsionparm[r_torsionparm_num].name2);
					strcpy(r_torsionparm[r_torsionparm_num].name2, r_torsionparm[r_torsionparm_num].name3);
					strcpy(r_torsionparm[r_torsionparm_num].name3, tmpchar);

					strcpy(tmpchar, r_torsionparm[r_torsionparm_num].name1);
					strcpy(r_torsionparm[r_torsionparm_num].name1, r_torsionparm[r_torsionparm_num].name4);
					strcpy(r_torsionparm[r_torsionparm_num].name4, tmpchar);
				}	
				if(strcmp(r_torsionparm[r_torsionparm_num].name2, r_torsionparm[r_torsionparm_num].name3) == 0) {
					if(strcmp(r_torsionparm[r_torsionparm_num].name1, r_torsionparm[r_torsionparm_num].name4) == 0) {
						strcpy(tmpchar, r_torsionparm[r_torsionparm_num].name1);
						strcpy(r_torsionparm[r_torsionparm_num].name1, r_torsionparm[r_torsionparm_num].name4);
						strcpy(r_torsionparm[r_torsionparm_num].name4, tmpchar);
					}
				}
                        	sscanf(&line[11], "%d%lf%lf%lf",
                                   	&r_torsionparm[r_torsionparm_num].mul,
                                   	&r_torsionparm[r_torsionparm_num].force,
                                   	&r_torsionparm[r_torsionparm_num].phase,
                                   	&r_torsionparm[r_torsionparm_num].fterm);
				if(strstr(line, "ATTN")  != NULL)
		        		r_torsionparm[r_torsionparm_num].attn = 1;
				}
                        r_torsionparm_num++;
                }
                if (iindex == 1) {
			if(mode == 1) {
                        	tmpnum = 0;
                        	tmpc1[0] = line[0];
				if(line[1] == ' ')
                        		tmpc1[1] = '\0';
				else {
                        		tmpc1[1] = line[1];
                        		tmpc1[2] = '\0';
				}
                        	tmpc2[0] = line[3];
				if(line[4] == ' ')
                        		tmpc2[1] = '\0';
				else {
                        		tmpc2[1] = line[4];
                        		tmpc2[2] = '\0';
				}
                        	tmpc3[0] = line[6];
				if(line[7] == ' ')
                        		tmpc3[1] = '\0';
				else {
                        		tmpc3[1] = line[7];
                        		tmpc3[2] = '\0';
				}
                        	tmpc4[0] = line[9];
				if(line[10] == ' ')
                        		tmpc4[1] = '\0';
				else {
                        		tmpc4[1] = line[10];
                        		tmpc4[2] = '\0';
				}
                        	if(line[0] == 'X') tmpnum++;
                        	if(line[3] == 'X') tmpnum++;
                        	if(line[6] == 'X') tmpnum++;
                        	if(line[9] == 'X') tmpnum++;
/*	we should not change the positions of 'X' */
				if(tmpnum == 0) {
                        		if(strcmp(tmpc1, tmpc2) > 0)  {
                                		strcpy(tmpc, tmpc2);
                                		strcpy(tmpc2, tmpc1);
                                		strcpy(tmpc1, tmpc);
                        		}
                        		if(strcmp(tmpc1, tmpc4) > 0)  {
                                		strcpy(tmpc, tmpc4);
                                		strcpy(tmpc4, tmpc1);
                                		strcpy(tmpc1, tmpc);
                        		}
                        		if(strcmp(tmpc2, tmpc4) > 0)  {
                                		strcpy(tmpc, tmpc4);
                                		strcpy(tmpc4, tmpc2);
                                		strcpy(tmpc2, tmpc);
                        		}
				}
				if(tmpnum == 1 || (tmpnum == 2 && line[6] == 'X')) {
					if(line[0] == 'X') 
                        			if(strcmp(tmpc2, tmpc4) > 0)  {
                                			strcpy(tmpc, tmpc4);
                                			strcpy(tmpc4, tmpc2);
                                			strcpy(tmpc2, tmpc);
                        			}
					if(line[3] == 'X') 
                        			if(strcmp(tmpc1, tmpc4) > 0)  {
                                			strcpy(tmpc, tmpc4);
                                			strcpy(tmpc4, tmpc1);
                                			strcpy(tmpc1, tmpc);
                        			}
					if(line[9] == 'X') 
                        			if(strcmp(tmpc1, tmpc2) > 0)  {
                                			strcpy(tmpc, tmpc1);
                                			strcpy(tmpc1, tmpc2);
                                			strcpy(tmpc2, tmpc);
                        			}
				}
		       		r_improperparm[r_improperparm_num].attn = 0;
                        	strcpy(r_improperparm[r_improperparm_num].name1, tmpc1);
                        	strcpy(r_improperparm[r_improperparm_num].name2, tmpc2);
                        	strcpy(r_improperparm[r_improperparm_num].name3, tmpc3);
                        	strcpy(r_improperparm[r_improperparm_num].name4, tmpc4);
                        	sscanf(&line[11], "%lf%lf%lf",
                                   	&r_improperparm[r_improperparm_num].force,
                                   	&r_improperparm[r_improperparm_num].phase,
                                   	&r_improperparm[r_improperparm_num].fterm);
                        	r_improperparm[r_improperparm_num].mul = 1; 
                        	r_improperparm[r_improperparm_num].numX = tmpnum; 
				if(strstr(line, "ATTN")  != NULL)
		        		r_improperparm[r_improperparm_num].attn = 1;
				}
                        r_improperparm_num++;
                }
                if (vindex == 1) {
			if(mode == 1) {
		       		r_vdwparm[r_vdwparm_num].attn = 0;
                        	sscanf(line, "%s%lf%lf", r_vdwparm[r_vdwparm_num].name,
                                   	&r_vdwparm[r_vdwparm_num].radius, &r_vdwparm[r_vdwparm_num].pot);
				if(strstr(line, "ATTN")  != NULL)
		        		r_vdwparm[r_vdwparm_num].attn = 1;
			}
                        r_vdwparm_num++;
                }

}
fclose(fp);
if(debug == 1 && mode == 1) {
	printf("\n-- Begin read frcmod file as input: %s --\n\n", filename);
        printf("\nMASS\n");
        for(i=0;i<r_atomtype_num;i++)
                printf("\n%-2s %9.4lf %9.4lf %5d", r_atomtype[i].name, r_atomtype[i].mass,atomtype[i].pol, 
                                                   r_atomtype[i].attn);

        printf("\n\nBOND\n");
        for(i=0;i<r_bondparm_num;i++)
                printf("\n%-2s %-2s %9.4lf %9.4lf %5d", r_bondparm[i].name1, r_bondparm[i].name2,
                r_bondparm[i].force, r_bondparm[i].length, r_bondparm[i].attn);

        printf("\n\nANGLE\n");
        for(i=0;i<r_angleparm_num;i++)
                printf("\n%-2s %-2s %-2s %9.4lf %9.4lf %5d", r_angleparm[i].name1, r_angleparm[i].name2,
                r_angleparm[i].name3, r_angleparm[i].force, r_angleparm[i].angle, r_angleparm[i].attn);

        printf("\n\nTORSION\n");
        for(i=0;i<r_torsionparm_num;i++)
                printf("\n%-2s %-2s %-2s %-2s %9.4lf %9.4lf %9.4lf %5d", r_torsionparm[i].name1,
                r_torsionparm[i].name2, r_torsionparm[i].name3, r_torsionparm[i].name4,
                r_torsionparm[i].phase, r_torsionparm[i].force, r_torsionparm[i].fterm,
                r_torsionparm[i].attn);

        printf("\n\nIMPROPER\n");
        for(i=0;i<r_improperparm_num;i++)
                printf("\n%-2s %-2s %-2s %-2s %9.4lf %9.4lf %9.4lf %5d", r_improperparm[i].name1,
                r_improperparm[i].name2, r_improperparm[i].name3, r_improperparm[i].name4,
                r_improperparm[i].phase, r_improperparm[i].force, r_improperparm[i].fterm,
                r_improperparm[i].attn);
        printf("\n\nVDW\n");
        for(i=0;i<r_vdwparm_num;i++)
                printf("\n%-2s %9.4lf %9.4lf %5d", r_vdwparm[i].name, r_vdwparm[i].radius, r_vdwparm[i].pot,
                                                   r_vdwparm[i].attn);
	printf("\n-- Finished reading frcmod file as input: %s --\n\n", filename);
}

}



void readfrcmod(char *filename) {
FILE *fp;
int mindex = 0;
int bindex = 0;
int aindex = 0;
int tindex = 0;
int iindex = 0;
int vindex = 0;
int pos_tor = 0;
int tmpnum  = 0;
char line[MAXCHAR];
char tmpc[MAXCHAR];
char tmpc1[MAXCHAR];
char tmpc2[MAXCHAR];
char tmpc3[MAXCHAR];
char tmpc4[MAXCHAR];
	if ((fp = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open file %s in readfrcmod(), exit\n", filename);
		exit(1);
	}
	for (;;) {
		if (fgets(line, MAXCHAR, fp) == NULL)
			break;
		if (strncmp(line, "MASS", 4) == 0) {
			mindex = 1;
			continue;
		}
		if (strncmp(line, "BOND", 4) == 0) { 
			bindex = 1;
			continue;
		}
		if (strncmp(line, "ANGL", 4) == 0){ 
			aindex = 1;
			continue;
		}
		if (strncmp(line, "DIHE", 4) == 0) { 
			pos_tor= ftell(fp);
			tindex = 1;
			continue;
		}
		if (strncmp(line, "IMPROPER", 8) == 0) { 
			iindex = 1;
			continue;
		}
		if (strncmp(line, "NONBON", 6) == 0) { 
			vindex = 1;
			continue;
		}
		if (strlen(line) <= 2) {
			if(mindex == 1) mindex = 0;
			if(bindex == 1) bindex = 0;
			if(aindex == 1) aindex = 0;
			if(tindex == 1) {
				fseek(fp,pos_tor,0);
				tindex = 2;
				continue;
			}
			if(tindex == 2) tindex = 0;
			if(iindex == 1) iindex = 0;
			if(vindex == 1) vindex = 0;
			continue;
		}
                if (mindex == 1) {
                        sscanf(line, "%s%lf%lf", atomtype[atomtypenum].name,
                                   &atomtype[atomtypenum].mass,
                                   &atomtype[atomtypenum].pol);
                        atomtypenum++;
                        if (atomtypenum >= maxatomtype) {
                                maxatomtype += MAX_FF_ATOMTYPE;
                                atomtype =
                                        (ATOMTYPE *) realloc(atomtype,
                                                                                 sizeof(ATOMTYPE) * maxatomtype);
                                if (atomtype == NULL) {
                                        fprintf(stdout,
                                                        "memory allocation error for *atomtype\n");
                                        exit(1);
                                }
                        }
                }
                if (bindex == 1) {
                        bondparm[bondparmnum].name1[0] = line[0];
			if(line[1]== ' ') 
                        	bondparm[bondparmnum].name1[1] = '\0';
			else {
				bondparm[bondparmnum].name1[1] = line[1];
				bondparm[bondparmnum].name1[2] = '\0';
			}
                        bondparm[bondparmnum].name2[0] = line[3];
			if(line[4]== ' ') 
                        	bondparm[bondparmnum].name2[1] = '\0';
			else {
				bondparm[bondparmnum].name2[1] = line[4];
				bondparm[bondparmnum].name2[2] = '\0';
			}
				
                        sscanf(&line[5], "%lf%lf", &bondparm[bondparmnum].force,
                                   &bondparm[bondparmnum].length);
                        bondparmnum++;
                        if (bondparmnum >= maxbondparm) {
                                maxbondparm += MAX_FF_BOND;
                                bondparm =
                                        (BOND_FF *) realloc(bondparm,
                                                                                sizeof(BOND_FF) * maxbondparm);
                                if (bondparm == NULL) {
                                        fprintf(stdout,
                                                        "memory allocation error for *bondparm\n");
                                        exit(1);
                                }
                        }
                }
                if (aindex == 1) {
                        angleparm[angleparmnum].name1[0] = line[0];
			if(line[1]== ' ') 
                        	angleparm[angleparmnum].name1[1] = '\0';
			else {
				angleparm[angleparmnum].name1[1] = line[1];
				angleparm[angleparmnum].name1[2] = '\0';
			}
                        angleparm[angleparmnum].name2[0] = line[3];
			if(line[4]== ' ') 
                        	angleparm[angleparmnum].name2[1] = '\0';
			else {
				angleparm[angleparmnum].name2[1] = line[4];
				angleparm[angleparmnum].name2[2] = '\0';
			}
                        angleparm[angleparmnum].name3[0] = line[6];
			if(line[7]== ' ') 
                        	angleparm[angleparmnum].name3[1] = '\0';
			else {
				angleparm[angleparmnum].name3[1] = line[7];
				angleparm[angleparmnum].name3[2] = '\0';
			}
                        sscanf(&line[8], "%lf%lf", &angleparm[angleparmnum].force,
                                   &angleparm[angleparmnum].angle);
                        angleparmnum++;
                        if (angleparmnum >= maxangleparm) {
                                maxangleparm += MAX_FF_ANGLE;
                                angleparm =
                                        (ANGLE *) realloc(angleparm,
                                                                          sizeof(ANGLE) * maxangleparm);
                                if (angleparm == NULL) {
                                        fprintf(stdout,
                                                        "memory allocation error for *angleparm\n");
                                        exit(1);
                                }
                        }
                }
                if (tindex == 1) { /*we first only read special torsional angle parameters*/
                        if(line[0] == 'X' || line[3] == 'X' || line[6] == 'X' || line[9] == 'X')
                                continue;
                        torsionparm[torsionparmnum].name1[0] = line[0];
			if(line[1]== ' ') 
                        	torsionparm[torsionparmnum].name1[1] = '\0';
			else {
				torsionparm[torsionparmnum].name1[1] = line[1];
				torsionparm[torsionparmnum].name1[2] = '\0';
			}
                        torsionparm[torsionparmnum].name2[0] = line[3];
			if(line[4]== ' ') 
                        	torsionparm[torsionparmnum].name2[1] = '\0';
			else {
				torsionparm[torsionparmnum].name2[1] = line[4];
				torsionparm[torsionparmnum].name2[2] = '\0';
			}
                        torsionparm[torsionparmnum].name3[0] = line[6];
			if(line[7]== ' ') 
                        	torsionparm[torsionparmnum].name3[1] = '\0';
			else {
				torsionparm[torsionparmnum].name3[1] = line[7];
				torsionparm[torsionparmnum].name3[2] = '\0';
			}
                        torsionparm[torsionparmnum].name4[0] = line[9];
			if(line[10]== ' ') 
                        	torsionparm[torsionparmnum].name4[1] = '\0';
			else {
				torsionparm[torsionparmnum].name4[1] = line[10];
				torsionparm[torsionparmnum].name4[2] = '\0';
			}
                        sscanf(&line[11], "%d%lf%lf%lf",
                                   &torsionparm[torsionparmnum].mul,
                                   &torsionparm[torsionparmnum].force,
                                   &torsionparm[torsionparmnum].phase,
                                   &torsionparm[torsionparmnum].fterm);
                        torsionparmnum++;
                        if (torsionparmnum >= maxtorsionparm) {
                                maxtorsionparm += MAX_FF_TORSION;
                                torsionparm =
                                        (TORSION *) realloc(torsionparm,
                                                                                sizeof(TORSION) * maxtorsionparm);
                                if (torsionparm == NULL) {
                                        fprintf(stdout,
                                                        "memory allocation error for *torsionparm\n");
                                        exit(1);
                                }
                        }
                }
                if (tindex == 2) {
                        if(line[0] != 'X' && line[3] != 'X' && line[6] != 'X' && line[9] != 'X')
                                continue;
                        torsionparm[torsionparmnum].name1[0] = line[0];
			if(line[1]== ' ') 
                        	torsionparm[torsionparmnum].name1[1] = '\0';
			else {
				torsionparm[torsionparmnum].name1[1] = line[1];
				torsionparm[torsionparmnum].name1[2] = '\0';
			}
                        torsionparm[torsionparmnum].name2[0] = line[3];
			if(line[4]== ' ') 
                        	torsionparm[torsionparmnum].name2[1] = '\0';
			else {
				torsionparm[torsionparmnum].name2[1] = line[4];
				torsionparm[torsionparmnum].name2[2] = '\0';
			}
                        torsionparm[torsionparmnum].name3[0] = line[6];
			if(line[7]== ' ') 
                        	torsionparm[torsionparmnum].name3[1] = '\0';
			else {
				torsionparm[torsionparmnum].name3[1] = line[7];
				torsionparm[torsionparmnum].name3[2] = '\0';
			}
                        torsionparm[torsionparmnum].name4[0] = line[9];
			if(line[10]== ' ') 
                        	torsionparm[torsionparmnum].name4[1] = '\0';
			else {
				torsionparm[torsionparmnum].name4[1] = line[10];
				torsionparm[torsionparmnum].name4[2] = '\0';
			}
                        sscanf(&line[11], "%d%lf%lf%lf",
                                   &torsionparm[torsionparmnum].mul,
                                   &torsionparm[torsionparmnum].force,
                                   &torsionparm[torsionparmnum].phase,
                                   &torsionparm[torsionparmnum].fterm);
                        torsionparmnum++;
                        if (torsionparmnum >= maxtorsionparm) {
                                maxtorsionparm += MAX_FF_TORSION;
                                torsionparm =
                                        (TORSION *) realloc(torsionparm,
                                                                                sizeof(TORSION) * maxtorsionparm);
                                if (torsionparm == NULL) {
                                        fprintf(stdout,
                                                        "memory allocation error for *torsionparm\n");
                                        exit(1);
                                }
                        }
                }
                if (iindex == 1) {
                        tmpnum = 0;
                        tmpc1[0] = line[0];
			if(line[1] == ' ')
                        	tmpc1[1] = '\0';
			else {
                        	tmpc1[1] = line[1];
                        	tmpc1[2] = '\0';
			}
                        tmpc2[0] = line[3];
			if(line[4] == ' ')
                        	tmpc2[1] = '\0';
			else {
                        	tmpc2[1] = line[4];
                        	tmpc2[2] = '\0';
			}
                        tmpc3[0] = line[6];
			if(line[7] == ' ')
                        	tmpc3[1] = '\0';
			else {
                        	tmpc3[1] = line[7];
                        	tmpc3[2] = '\0';
			}
                        tmpc4[0] = line[9];
			if(line[10] == ' ')
                        	tmpc4[1] = '\0';
			else {
                        	tmpc4[1] = line[10];
                        	tmpc4[2] = '\0';
			}
                        if(line[0] == 'X') tmpnum++;
                        if(line[3] == 'X') tmpnum++;
                        if(line[6] == 'X') tmpnum++;
                        if(line[9] == 'X') tmpnum++;
/*	we should not change the positions of 'X' */
			if(tmpnum == 0) {
                        	if(strcmp(tmpc1, tmpc2) > 0)  {
                                	strcpy(tmpc, tmpc2);
                                	strcpy(tmpc2, tmpc1);
                                	strcpy(tmpc1, tmpc);
                        	}
                        	if(strcmp(tmpc1, tmpc4) > 0)  {
                                	strcpy(tmpc, tmpc4);
                                	strcpy(tmpc4, tmpc1);
                                	strcpy(tmpc1, tmpc);
                        	}
                        	if(strcmp(tmpc2, tmpc4) > 0)  {
                                	strcpy(tmpc, tmpc4);
                                	strcpy(tmpc4, tmpc2);
                                	strcpy(tmpc2, tmpc);
                        	}
			}
			if(tmpnum == 1 || (tmpnum == 2 && line[6] == 'X')) {
				if(line[0] == 'X') 
                        		if(strcmp(tmpc2, tmpc4) > 0)  {
                                		strcpy(tmpc, tmpc4);
                                		strcpy(tmpc4, tmpc2);
                                		strcpy(tmpc2, tmpc);
                        		}
				if(line[3] == 'X') 
                        		if(strcmp(tmpc1, tmpc4) > 0)  {
                                		strcpy(tmpc, tmpc4);
                                		strcpy(tmpc4, tmpc1);
                                		strcpy(tmpc1, tmpc);
                        		}
				if(line[9] == 'X') 
                        		if(strcmp(tmpc1, tmpc2) > 0)  {
                                		strcpy(tmpc, tmpc1);
                                		strcpy(tmpc1, tmpc2);
                                		strcpy(tmpc2, tmpc);
                        		}
			}
                        strcpy(improperparm[improperparmnum].name1, tmpc1);
                        strcpy(improperparm[improperparmnum].name2, tmpc2);
                        strcpy(improperparm[improperparmnum].name3, tmpc3);
                        strcpy(improperparm[improperparmnum].name4, tmpc4);
                        sscanf(&line[11], "%lf%lf%lf",
                                   &improperparm[improperparmnum].force,
                                   &improperparm[improperparmnum].phase,
                                   &improperparm[improperparmnum].fterm);
                        improperparm[improperparmnum].mul = 1; 
                        improperparm[improperparmnum].numX = tmpnum; 
                        improperparmnum++;
                        if (improperparmnum >= maximproperparm) {
                                maximproperparm += MAX_FF_IMPROPER;
                                improperparm =
                                        (IMPROPER *) realloc(improperparm,
                                                                                 sizeof(IMPROPER) *
                                                                                 maximproperparm);
                                if (improperparm == NULL) {
                                        fprintf(stdout,
                                                        "memory allocation error for *improperparm\n");
                                        exit(1);
                                }
                        }
                }
                if (vindex == 1) {
                        sscanf(line, "%s%lf%lf", vdwparm[vdwparmnum].name,
                                   &vdwparm[vdwparmnum].radius, &vdwparm[vdwparmnum].pot);
                        vdwparmnum++;
                        if (vdwparmnum >= maxvdwparm) {
                                maxvdwparm += MAX_FF_VDW;
                                vdwparm =
                                        (VDW *) realloc(vdwparm, sizeof(VDW) * maxvdwparm);
                                if (vdwparm == NULL) {
                                        fprintf(stdout,
                                                        "memory allocation error for *vdwparm\n");
                                        exit(1);
                                }
                        }
                }

}
fclose(fp);
if(debug == 1) {
	printf("\n-- Begin read frcmod file : %s --\n\n", filename);
        printf("\nMASS\n");
        for(i=last_atomtypenum;i<atomtypenum;i++)
                printf("\n%-2s %9.4lf %9.4lf", atomtype[i].name, atomtype[i].mass,atomtype[i].pol);

        printf("\n\nBOND\n");
        for(i=last_bondparmnum;i<bondparmnum;i++)
                printf("\n%-2s %-2s %9.4lf %9.4lf", bondparm[i].name1, bondparm[i].name2,
                bondparm[i].force, bondparm[i].length);

        printf("\n\nANGLE\n");
        for(i=last_angleparmnum;i<angleparmnum;i++)
                printf("\n%-2s %-2s %-2s %9.4lf %9.4lf", angleparm[i].name1, angleparm[i].name2,
                angleparm[i].name3, angleparm[i].force, angleparm[i].angle);

        printf("\n\nTORSION\n");
        for(i=last_torsionparmnum;i<torsionparmnum;i++)
                printf("\n%-2s %-2s %-2s %-2s %9.4lf %9.4lf %9.4lf", torsionparm[i].name1,
                torsionparm[i].name2, torsionparm[i].name3,torsionparm[i].name4,
                torsionparm[i].phase, torsionparm[i].force, torsionparm[i].fterm);

        printf("\n\nIMPROPER\n");
        for(i=last_improperparmnum;i<improperparmnum;i++)
                printf("\n%-2s %-2s %-2s %-2s %9.4lf %9.4lf %9.4lf", improperparm[i].name1,
                improperparm[i].name2, improperparm[i].name3,improperparm[i].name4,
                improperparm[i].phase, improperparm[i].force, improperparm[i].fterm);
        printf("\n\nVDW\n");
        for(i=last_vdwparmnum;i<vdwparmnum;i++)
                printf("\n%-2s %9.4lf %9.4lf", vdwparm[i].name, vdwparm[i].radius, vdwparm[i].pot);

	printf("\n-- Finish reading frcmod file : %s -- \n\n", filename);
}

}

void readparm(char *filename)
{
	int mindex = -1;
	int bindex = 0;
	int aindex = 0;
	int tindex = 0;
	int iindex = 0;
	int vindex = 0;
	int num = 0;
	int tmpnum;
	int i, j, k;
	int flag;
	int vdwparmnum_old;
	int pos_tor = 0;
	FILE *fp;
	char line[MAXCHAR];
	char tmpc[MAXCHAR];
	char tmpc1[MAXCHAR];
	char tmpc2[MAXCHAR];
	char tmpc3[MAXCHAR];
	char tmpc4[MAXCHAR];
	if ((fp = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open file %s in readparm(), exit\n", filename);
		exit(1);
	}
	for (;;) {
		if (fgets(line, MAXCHAR, fp) == NULL)
			break;
		num++;
		if (mindex == -1 && num == 1) {
			mindex = 1;
			continue;
		}
		if (mindex == 1 && spaceline(line) == 1) {
			mindex = 0;
			num = 0;
			bindex = -1;
			continue;
		}
		if (bindex == -1 && num == 1) {
			bindex = 1;
			continue;
		}
		if (bindex == 1 && spaceline(line) == 1) {
			bindex = 0;
			aindex = 1;
			continue;
		}
		if (aindex == 1 && spaceline(line) == 1) {
			aindex = 0;
			tindex = 1;
			pos_tor= ftell(fp);
			continue;
		}
		if (tindex == 1 && spaceline(line) == 1) {
			tindex = 2;
			fseek(fp,pos_tor,0);
			continue;
		}
		if (tindex == 2 && spaceline(line) == 1) {
			tindex = 0;
			iindex = 1;
			continue;
		}
		if (iindex == 1 && spaceline(line) == 1) {
			iindex = 0;
			vindex = -1;
			num = 0;
			continue;
		}
		if (vindex == -1 && num == 2) {
			vindex = 1;
			continue;
		}
		if (vindex == 2 && spaceline(line) == 1) {
			vindex = 0;
			continue;
		}
		if (vindex == 1)
			if (strncmp(line, "MOD4", 4) == 0) {
				vindex = 2;
				continue;
			}
		if (mindex == 1) {
			sscanf(line, "%s%lf%lf", atomtype[atomtypenum].name,
				   &atomtype[atomtypenum].mass,
				   &atomtype[atomtypenum].pol);
			atomtypenum++;
			if (atomtypenum >= maxatomtype) {
				maxatomtype += MAX_FF_ATOMTYPE;
				atomtype =
					(ATOMTYPE *) realloc(atomtype,
										 sizeof(ATOMTYPE) * maxatomtype);
				if (atomtype == NULL) {
					fprintf(stdout,
							"memory allocation error for *atomtype\n");
					exit(1);
				}
			}
		}
		if (bindex == 1) {
			bondparm[bondparmnum].name1[0] = line[0];
                        if(line[1]== ' ')
                                bondparm[bondparmnum].name1[1] = '\0';
                        else {
                                bondparm[bondparmnum].name1[1] = line[1];
                                bondparm[bondparmnum].name1[2] = '\0';
                        }
			bondparm[bondparmnum].name2[0] = line[3];
                        if(line[4]== ' ')
                                bondparm[bondparmnum].name2[1] = '\0';
                        else {
                                bondparm[bondparmnum].name2[1] = line[4];
                                bondparm[bondparmnum].name2[2] = '\0';
                        }
			sscanf(&line[5], "%lf%lf", &bondparm[bondparmnum].force,
				   &bondparm[bondparmnum].length);
			bondparmnum++;
			if (bondparmnum >= maxbondparm) {
				maxbondparm += MAX_FF_BOND;
				bondparm =
					(BOND_FF *) realloc(bondparm,
										sizeof(BOND_FF) * maxbondparm);
				if (bondparm == NULL) {
					fprintf(stdout,
							"memory allocation error for *bondparm\n");
					exit(1);
				}
			}
		}
		if (aindex == 1) {
			angleparm[angleparmnum].name1[0] = line[0];
                        if(line[1]== ' ')
                                angleparm[angleparmnum].name1[1] = '\0';
                        else {
                                angleparm[angleparmnum].name1[1] = line[1];
                                angleparm[angleparmnum].name1[2] = '\0';
                        }
			angleparm[angleparmnum].name2[0] = line[3];
                        if(line[4]== ' ')
                                angleparm[angleparmnum].name2[1] = '\0';
                        else {
                                angleparm[angleparmnum].name2[1] = line[4];
                                angleparm[angleparmnum].name2[2] = '\0';
                        }
			angleparm[angleparmnum].name3[0] = line[6];
                        if(line[7]== ' ')
                                angleparm[angleparmnum].name3[1] = '\0';
                        else {
                                angleparm[angleparmnum].name3[1] = line[7];
                                angleparm[angleparmnum].name3[2] = '\0';
                        }
			sscanf(&line[8], "%lf%lf", &angleparm[angleparmnum].force,
				   &angleparm[angleparmnum].angle);
			angleparmnum++;
			if (angleparmnum >= maxangleparm) {
				maxangleparm += MAX_FF_ANGLE;
				angleparm =
					(ANGLE *) realloc(angleparm,
									  sizeof(ANGLE) * maxangleparm);
				if (angleparm == NULL) {
					fprintf(stdout,
							"memory allocation error for *angleparm\n");
					exit(1);
				}
			}
		}
		if (tindex == 1) { /*we first only read special torsional angle parameters*/
			if(line[0] == 'X' || line[3] == 'X' || line[6] == 'X' || line[9] == 'X')
				continue;
			torsionparm[torsionparmnum].name1[0] = line[0];
                        if(line[1]== ' ')
                                torsionparm[torsionparmnum].name1[1] = '\0';
                        else {
                                torsionparm[torsionparmnum].name1[1] = line[1];
                                torsionparm[torsionparmnum].name1[2] = '\0';
                        }
			torsionparm[torsionparmnum].name2[0] = line[3];
                        if(line[4]== ' ')
                                torsionparm[torsionparmnum].name2[1] = '\0';
                        else {
                                torsionparm[torsionparmnum].name2[1] = line[4];
                                torsionparm[torsionparmnum].name2[2] = '\0';
                        }
			torsionparm[torsionparmnum].name3[0] = line[6];
                        if(line[7]== ' ')
                                torsionparm[torsionparmnum].name3[1] = '\0';
                        else {
                                torsionparm[torsionparmnum].name3[1] = line[7];
                                torsionparm[torsionparmnum].name3[2] = '\0';
                        }
			torsionparm[torsionparmnum].name4[0] = line[9];
                        if(line[10]== ' ')
                                torsionparm[torsionparmnum].name4[1] = '\0';
                        else {
                                torsionparm[torsionparmnum].name4[1] = line[10];
                                torsionparm[torsionparmnum].name4[2] = '\0';
                        }
			sscanf(&line[11], "%d%lf%lf%lf",
				   &torsionparm[torsionparmnum].mul,
				   &torsionparm[torsionparmnum].force,
				   &torsionparm[torsionparmnum].phase,
				   &torsionparm[torsionparmnum].fterm);
			torsionparmnum++;
			if (torsionparmnum >= maxtorsionparm) {
				maxtorsionparm += MAX_FF_TORSION;
				torsionparm =
					(TORSION *) realloc(torsionparm,
										sizeof(TORSION) * maxtorsionparm);
				if (torsionparm == NULL) {
					fprintf(stdout,
							"memory allocation error for *torsionparm\n");
					exit(1);
				}
			}
		}
                if (tindex == 2) {
			if(line[0] != 'X' && line[3] != 'X' && line[6] != 'X' && line[9] != 'X')
				continue;
                        torsionparm[torsionparmnum].name1[0] = line[0];
                        if(line[1]== ' ')
                                torsionparm[torsionparmnum].name1[1] = '\0';
                        else {
                                torsionparm[torsionparmnum].name1[1] = line[1];
                                torsionparm[torsionparmnum].name1[2] = '\0';
                        }
                        torsionparm[torsionparmnum].name2[0] = line[3];
                        if(line[4]== ' ')
                                torsionparm[torsionparmnum].name2[1] = '\0';
                        else {
                                torsionparm[torsionparmnum].name2[1] = line[4];
                                torsionparm[torsionparmnum].name2[2] = '\0';
                        }
                        torsionparm[torsionparmnum].name3[0] = line[6];
                        if(line[7]== ' ') 
                                torsionparm[torsionparmnum].name3[1] = '\0';
                        else {
                                torsionparm[torsionparmnum].name3[1] = line[7];
                                torsionparm[torsionparmnum].name3[2] = '\0';
                        }
                        torsionparm[torsionparmnum].name4[0] = line[9];
                        if(line[10]== ' ') 
                                torsionparm[torsionparmnum].name4[1] = '\0';
                        else {
                                torsionparm[torsionparmnum].name4[1] = line[10];
                                torsionparm[torsionparmnum].name4[2] = '\0';
                        }
                        sscanf(&line[11], "%d%lf%lf%lf",
                                   &torsionparm[torsionparmnum].mul,
                                   &torsionparm[torsionparmnum].force,
                                   &torsionparm[torsionparmnum].phase,
                                   &torsionparm[torsionparmnum].fterm);
                        torsionparmnum++;
                        if (torsionparmnum >= maxtorsionparm) {
                                maxtorsionparm += MAX_FF_TORSION;
                                torsionparm =
                                        (TORSION *) realloc(torsionparm,
                                                                                sizeof(TORSION) * maxtorsionparm);
                                if (torsionparm == NULL) {
                                        fprintf(stdout,
                                                        "memory allocation error for *torsionparm\n");
                                        exit(1);
                                }
                        }
                }
		if (iindex == 1) {
			tmpnum = 0;
			tmpc1[0] = line[0];
                        if(line[1] == ' ') 
                                tmpc1[1] = '\0';
                        else {
                                tmpc1[1] = line[1];
                                tmpc1[2] = '\0';
                        }
			tmpc2[0] = line[3];
                        if(line[4] == ' ') 
                                tmpc2[1] = '\0';
                        else {
                                tmpc2[1] = line[4];
                                tmpc2[2] = '\0';
                        }
			tmpc3[0] = line[6];
                        if(line[7] == ' ') 
                                tmpc3[1] = '\0';
                        else {
                                tmpc3[1] = line[7];
                                tmpc3[2] = '\0';
                        }
			tmpc4[0] = line[9];
                        if(line[10] == ' ')
                                tmpc4[1] = '\0';
                        else {
                                tmpc4[1] = line[10];
                                tmpc4[2] = '\0';
                        }
			if(line[0] == 'X') tmpnum++;
			if(line[3] == 'X') tmpnum++;
			if(line[6] == 'X') tmpnum++;
			if(line[9] == 'X') tmpnum++;

/*      we should not change the position of 'X'*/
                        if(tmpnum == 0) {
                                if(strcmp(tmpc1, tmpc2) > 0)  {
                                        strcpy(tmpc, tmpc2);
                                        strcpy(tmpc2, tmpc1);
                                        strcpy(tmpc1, tmpc);
                                }
                                if(strcmp(tmpc1, tmpc4) > 0)  {
                                        strcpy(tmpc, tmpc4);
                                        strcpy(tmpc4, tmpc1);
                                        strcpy(tmpc1, tmpc);
                                }
                                if(strcmp(tmpc2, tmpc4) > 0)  {
                                        strcpy(tmpc, tmpc4);
                                        strcpy(tmpc4, tmpc2);
                                        strcpy(tmpc2, tmpc);
                                }
                        }
                        if(tmpnum == 1 || (tmpnum == 2 && line[6] == 'X')) {
                                if(line[0] == 'X')
                                        if(strcmp(tmpc2, tmpc4) > 0)  {
                                                strcpy(tmpc, tmpc4);
                                                strcpy(tmpc4, tmpc2);
                                                strcpy(tmpc2, tmpc);
                                        }
                                if(line[3] == 'X') 
                                        if(strcmp(tmpc1, tmpc4) > 0)  {
                                                strcpy(tmpc, tmpc4);
                                                strcpy(tmpc4, tmpc1);
                                                strcpy(tmpc1, tmpc);
                                        }
                                if(line[9] == 'X')
                                        if(strcmp(tmpc1, tmpc2) > 0)  {
                                                strcpy(tmpc, tmpc1);
                                                strcpy(tmpc1, tmpc2);
                                                strcpy(tmpc2, tmpc);
                                        }
                        }
                	strcpy(improperparm[improperparmnum].name1, tmpc1);
                	strcpy(improperparm[improperparmnum].name2, tmpc2);
                	strcpy(improperparm[improperparmnum].name3, tmpc3);
                	strcpy(improperparm[improperparmnum].name4, tmpc4);
			sscanf(&line[11], "%lf%lf%lf",
				   &improperparm[improperparmnum].force,
				   &improperparm[improperparmnum].phase,
				   &improperparm[improperparmnum].fterm);
			improperparm[improperparmnum].mul = 1; 
			improperparm[improperparmnum].numX = tmpnum; 
			improperparmnum++;
			if (improperparmnum >= maximproperparm) {
				maximproperparm += MAX_FF_IMPROPER;
				improperparm =
					(IMPROPER *) realloc(improperparm,
										 sizeof(IMPROPER) *
										 maximproperparm);
				if (improperparm == NULL) {
					fprintf(stdout,
							"memory allocation error for *improperparm\n");
					exit(1);
				}
			}
		}
		if (vindex == 1) {
			if (spaceline(line) == 1)
				continue;
			equ_vdw[equ_vdw_num].num = 0;
			tmpnum = 0;
			flag = 1;
			while (flag) {
				flag = 0;
				sscanf(&line[tmpnum], "%s",
					   equ_vdw[equ_vdw_num].name[equ_vdw[equ_vdw_num].
												 num++]);
				if (equ_vdw[equ_vdw_num].num >= MAX_EQU_VDW_NUM) {
					printf
						("\nError: number of equivalent vdw atoms exceeds MAX_EQU_VDW_NUM, exit\n");
					exit(1);
				}
				for (i = tmpnum; i < strlen(line) - 3; i++) {
					if (line[i - 1] != ' ' && i >= 1 && line[i] != ' ')
						continue;
					if (line[i] != ' ' && line[i + 1] != ' ') {
						tmpnum = i + 2;
						flag = 1;
						break;
					}
					if (line[i] != ' ' && line[i + 1] == ' ') {
						tmpnum = i + 1;
						flag = 1;
						break;
					}
				}
			}
			if (equ_vdw[equ_vdw_num].num >= 2)
				equ_vdw_num++;
			if (equ_vdw_num >= MAX_EQU_VDW_NUM) {
				printf
					("\nError: number of equivalent vdw parameters exceeds MAX_EQU_VDW, exit\n");
				exit(1);
			}
		}
		if (vindex == 2) {
			sscanf(line, "%s%lf%lf", vdwparm[vdwparmnum].name,
				   &vdwparm[vdwparmnum].radius, &vdwparm[vdwparmnum].pot);
			vdwparmnum++;
			if (vdwparmnum >= maxvdwparm) {
				maxvdwparm += MAX_FF_VDW;
				vdwparm =
					(VDW *) realloc(vdwparm, sizeof(VDW) * maxvdwparm);
				if (vdwparm == NULL) {
					fprintf(stdout,
							"memory allocation error for *vdwparm\n");
					exit(1);
				}
			}
		}
	}
	fclose(fp);

	if (equ_vdw_num > 0) {
		vdwparmnum_old = vdwparmnum;
		for (i = 0; i < equ_vdw_num; i++)
			for (j = 1; j < equ_vdw[i].num; j++) {
				if (strlen(equ_vdw[i].name[j]) < 1)
					continue;
				for (k = 0; k < vdwparmnum_old; k++)
					if (strcmp(vdwparm[k].name, equ_vdw[i].name[0]) == 0) {
						strcpy(vdwparm[vdwparmnum].name,
							   equ_vdw[i].name[j]);
						vdwparm[vdwparmnum].radius = vdwparm[k].radius;
						vdwparm[vdwparmnum].pot = vdwparm[k].pot;
						vdwparmnum++;
						break;
					}
			}
	}
	if(debug == 1) {
		printf("\n-- Begin read force field parameter file : %s -- \n\n", filename);
                printf("\nMASS\n");
                for(i=last_atomtypenum;i<atomtypenum;i++)
                        printf("\n%-2s %9.4lf %9.4lf", atomtype[i].name, atomtype[i].mass,atomtype[i].pol);

                printf("\n\nBOND\n");
                for(i=last_bondparmnum;i<bondparmnum;i++)
                        printf("\n%-2s %-2s %9.4lf %9.4lf", bondparm[i].name1, bondparm[i].name2,
                        bondparm[i].force, bondparm[i].length);

                printf("\n\nANGLE\n");
                for(i=last_angleparmnum;i<angleparmnum;i++)
                        printf("\n%-2s %-2s %-2s %9.4lf %9.4lf", angleparm[i].name1, angleparm[i].name2,
                        angleparm[i].name3, angleparm[i].force, angleparm[i].angle);

                printf("\n\nTORSION\n");
                for(i=last_torsionparmnum;i<torsionparmnum;i++)
                        printf("\n%-2s %-2s %-2s %-2s %9.4lf %9.4lf %9.4lf", torsionparm[i].name1,
                        torsionparm[i].name2, torsionparm[i].name3,torsionparm[i].name4,
                        torsionparm[i].phase, torsionparm[i].force, torsionparm[i].fterm);

                printf("\n\nIMPROPER\n");
                for(i=last_improperparmnum;i<improperparmnum;i++)
                        printf("\n%-2s %-2s %-2s %-2s %9.4lf %9.4lf %9.4lf", improperparm[i].name1,
                        improperparm[i].name2, improperparm[i].name3,improperparm[i].name4,
                        improperparm[i].phase, improperparm[i].force, improperparm[i].fterm);
                printf("\n\nVDW\n");
                for(i=last_vdwparmnum;i<vdwparmnum;i++)
                        printf("\n%-2s %9.4lf %9.4lf", vdwparm[i].name, vdwparm[i].radius, vdwparm[i].pot);
		printf("\n-- Finish reading force field parameter file %s --\n\n", filename);
	}
}

void read_parmchk_parm(char *filename)
{
	FILE *fp;
	char line[MAXCHAR];
	char atomtype[10];
	char equaname[10];
	char corrname[10];
	int  i,j, k;
	int  group;
	int  improper;
	int  ncorr;
	int  nequa;
	int  equtype;
	int  atomicnum;
	double bl,blf,ba,baf, cba, cbaf, tor, ctor, itor, ps;
	double mass;
	if ((fp = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open file %s in read_parmchk_parm(), exit\n", filename);
		exit(1);
	}

	parmnum = 0;

	for (;;) {
		if (fgets(line, MAXCHAR, fp) == NULL)
			break;
		if (strncmp("PARM", &line[0], 4) == 0) {
			equtype = 0;
			sscanf(&line[4], "%s%d%d%lf%d%d", atomtype, &improper, &group, &mass, &equtype, &atomicnum) ;
			strcpy(parm[parmnum].atomtype, atomtype);
			parm[parmnum].improper = improper;
			parm[parmnum].group    = group;
			parm[parmnum].mass     = mass;
			parm[parmnum].equtype  = equtype;
			parm[parmnum].atomicnum  = atomicnum;
/*	Initialization */
/*	for the sake of simplicity in the coding, an atom type equals to and corresponds to itself*/
/*      Equal atom types must also be corresponding atom types */
                        strcpy(parm[parmnum].equa[0].atomtype, atomtype);
                        strcpy(parm[parmnum].corr[0].atomtype, atomtype);
                        parm[parmnum].corr[0].bl   = 0;
                        parm[parmnum].corr[0].blf  = 0;
                        parm[parmnum].corr[0].ba   = 0;
                        parm[parmnum].corr[0].baf  = 0;
                        parm[parmnum].corr[0].cba  = 0;
                        parm[parmnum].corr[0].cbaf = 0;
                        parm[parmnum].corr[0].tor  = 0;
                        parm[parmnum].corr[0].ctor = 0;
                        parm[parmnum].corr[0].improper = 0;
                        parm[parmnum].corr[0].ps    = 0;
                        parm[parmnum].corr[0].type = 0;
                        parm[parmnum].corr[0].pid  = parmnum;
			nequa = 1;
			ncorr = 1;
			parm[parmnum].nequa = 1;
			parm[parmnum].ncorr = 1;
			parmnum++;
			if (parmnum >= maxparmnum) {
				maxparmnum += MAXPARM;
				parm = (PARM *) realloc(parm, sizeof(PARM) * maxparmnum);
				if (parm == NULL) {
					fprintf(stdout, "memory allocation error for *parm\n");
					exit(1);
				}
			}
		}
		if (strncmp("EQUA", &line[0], 4) == 0) {
			sscanf(&line[4], "%s",  equaname); 
			strcpy(parm[parmnum-1].equa[nequa].atomtype, equaname);
			parm[parmnum-1].nequa ++;
			nequa ++;
			if(parm[parmnum-1].nequa > MAXEQUA) {
				fprintf(stdout, "Too many equal atom types for Atom Type %d (%s)\n", parmnum, parm[parmnum-1].atomtype);
				exit(1);
			}
			strcpy(parm[parmnum-1].corr[ncorr].atomtype, equaname);
			parm[parmnum-1].corr[ncorr].bl  =0;
			parm[parmnum-1].corr[ncorr].blf =0;
			parm[parmnum-1].corr[ncorr].ba  =0;
			parm[parmnum-1].corr[ncorr].baf =0;
			parm[parmnum-1].corr[ncorr].cba =0;
			parm[parmnum-1].corr[ncorr].cbaf=0;
			parm[parmnum-1].corr[ncorr].tor =0;
			parm[parmnum-1].corr[ncorr].ctor=0;
			parm[parmnum-1].corr[ncorr].ps  =0;
			parm[parmnum-1].corr[ncorr].improper=0;
			parm[parmnum-1].corr[ncorr].pid = -1;
			parm[parmnum-1].corr[ncorr].type = 1;
			parm[parmnum-1].ncorr ++;
			ncorr ++;
			if(parm[parmnum-1].ncorr > MAXCORR) {
				fprintf(stdout, "Too many corresponding atom types for Atom Type %d (%s)\n", parmnum, parm[parmnum-1].atomtype);
				exit(1);
			}
		}
		if (strncmp("CORR", &line[0], 4) == 0) {
			bl  = 0;
			blf = 0;
			ba  = 0;
			baf = 0;
			cba = 0;
			cbaf= 0;
			tor = 0;
			ctor= 0;
			itor= 0;
			ps=0;
			sscanf(&line[4], "%s%lf%lf%lf%lf%lf%lf%lf%lf%lf",  
               corrname, &bl, &blf, &cba, &cbaf, &ba, &baf, &ctor, &tor, &ps); 
			strcpy(parm[parmnum-1].corr[ncorr].atomtype, corrname);
			parm[parmnum-1].corr[ncorr].bl  =bl;
			parm[parmnum-1].corr[ncorr].blf =blf;
			parm[parmnum-1].corr[ncorr].ba  =ba;
			parm[parmnum-1].corr[ncorr].baf =baf;
			parm[parmnum-1].corr[ncorr].cba =cba;
			parm[parmnum-1].corr[ncorr].cbaf=cbaf;
			parm[parmnum-1].corr[ncorr].tor =tor;
			parm[parmnum-1].corr[ncorr].ctor=ctor;
			parm[parmnum-1].corr[ncorr].improper=ps; /*assign general score to improper parameter*/
			parm[parmnum-1].corr[ncorr].ps=ps;
			parm[parmnum-1].corr[ncorr].pid = -1;
			parm[parmnum-1].corr[ncorr].type = 2;
			parm[parmnum-1].ncorr ++;
			ncorr ++;
			if(parm[parmnum-1].ncorr > MAXCORR) {
				fprintf(stdout, "Too many corresponding atom types for Atom Type %d (%s)\n", parmnum, parm[parmnum-1].atomtype);
				exit(1);
			}
		}
		if (strncmp("WEIGHT_BL", &line[0], 9) == 0)   	sscanf(&line[9], "%lf", &wt.BL);
		if (strncmp("WEIGHT_BLF", &line[0], 10) == 0) 	sscanf(&line[10], "%lf", &wt.BLF);
		if (strncmp("WEIGHT_BA", &line[0], 9) == 0)   	sscanf(&line[9], "%lf", &wt.BA);
		if (strncmp("WEIGHT_BAF", &line[0], 10) == 0) 	sscanf(&line[10], "%lf", &wt.BAF);
		if (strncmp("WEIGHT_X", &line[0], 8) == 0)   	sscanf(&line[8], "%lf", &wt.X);
		if (strncmp("WEIGHT_X3", &line[0], 9) == 0)  	sscanf(&line[9], "%lf", &wt.X3);
		if (strncmp("WEIGHT_BA_CTR",  &line[0], 13) == 0)  sscanf(&line[13], "%lf", &wt.BA_CTR);
		if (strncmp("WEIGHT_TOR_CTR", &line[0], 14) == 0)  sscanf(&line[14], "%lf", &wt.TOR_CTR);
		if (strncmp("WEIGHT_IMPROPER", &line[0], 15) == 0) sscanf(&line[15], "%lf", &wt.IMPROPER);
		if (strncmp("WEIGHT_GROUP", &line[0], 12) == 0)	   sscanf(&line[12], "%lf", &wt.GROUP);
		if (strncmp("WEIGHT_EQUTYPE", &line[0], 14) == 0)  sscanf(&line[14], "%lf", &wt.EQUTYPE);
		if (strncmp("THRESHOLD_BA", &line[0], 12) == 0)    sscanf(&line[12], "%lf", &THRESHOLD_BA);
		
		if (strncmp("DEFAULT_BL", &line[0], 10) == 0)   sscanf(&line[10], "%lf", &dv.BL);
		if (strncmp("DEFAULT_BLF", &line[0], 11) == 0) 	sscanf(&line[11], "%lf", &dv.BLF);
		if (strncmp("DEFAULT_BA", &line[0], 10) == 0)   sscanf(&line[10], "%lf", &dv.BA);
		if (strncmp("DEFAULT_BAF", &line[0], 11) == 0) 	sscanf(&line[11], "%lf", &dv.BAF);
		if (strncmp("DEFAULT_BA_CTR", &line[0], 14) == 0)  	sscanf(&line[14], "%lf", &dv.BA_CTR);
		if (strncmp("DEFAULT_BAF_CTR", &line[0], 15) == 0)  	sscanf(&line[15], "%lf", &dv.BAF_CTR);
		if (strncmp("DEFAULT_TOR", &line[0], 11) == 0) 		sscanf(&line[11], "%lf", &dv.TOR);
		if (strncmp("DEFAULT_TOR_CTR", &line[0], 15) == 0) 	sscanf(&line[15], "%lf", &dv.TOR_CTR);
		if (strncmp("DEFAULT_FRACT1", &line[0], 14) == 0) 	sscanf(&line[14], "%lf", &dv.FRACT1);
		if (strncmp("DEFAULT_FRACT2", &line[0], 14) == 0) 	sscanf(&line[14], "%lf", &dv.FRACT2);
	}
	fclose(fp);
	for(i=0;i<parmnum;i++) 
		for(j=0;j<parm[i].ncorr;j++) {
			if(parm[i].corr[j].bl   < 0) parm[i].corr[j].bl  = dv.BL; 
			if(parm[i].corr[j].blf  < 0) parm[i].corr[j].blf = dv.BLF; 
			if(parm[i].corr[j].ba   < 0) parm[i].corr[j].ba  = dv.BA ;
			if(parm[i].corr[j].baf  < 0) parm[i].corr[j].baf = dv.BAF;
			if(parm[i].corr[j].cba  < 0) parm[i].corr[j].cba = dv.BA_CTR;
			if(parm[i].corr[j].cbaf < 0) parm[i].corr[j].cbaf= dv.BAF_CTR;
			if(parm[i].corr[j].tor  < 0) parm[i].corr[j].tor = dv.TOR;
			if(parm[i].corr[j].ctor < 0) parm[i].corr[j].ctor= dv.TOR_CTR;
			parm[i].corr[j].ctor = parm[i].corr[j].ctor * dv.FRACT1  + parm[i].corr[j].ps * dv.FRACT2; 
		}
/*	now finding the pid of equal/corresponding atom types */
	for(i=0;i<parmnum;i++) {
		for(j=0;j<parm[i].nequa;j++) 
			for(k=0;k<parmnum;k++) 
				if(strcmp(parm[i].equa[j].atomtype, parm[k].atomtype) == 0) {
					parm[i].equa[j].pid = k;
					break;
				}
		for(j=0;j<parm[i].ncorr;j++) 
			for(k=0;k<parmnum;k++) 
				if(strcmp(parm[i].corr[j].atomtype, parm[k].atomtype) == 0) {
					parm[i].corr[j].pid = k;
					break;
				}
	}
	if(debug == 1) {
		printf("\n-- Begin read parmchk parm file %s --\n\n", filename);
		for(i=0;i<parmnum;i++) {
			printf("PARM %5d %5s %5d %5d %9.4lf %5d\n", i+1, parm[i].atomtype, parm[i].improper, parm[i].group, parm[i].mass, parm[i].equtype);
			for(j=0;j<parm[i].nequa;j++)
				printf("EQUA %5d %5s %5d \n", 
				j+1, parm[i].equa[j].atomtype, parm[i].equa[j].pid + 1); 
			for(j=0;j<parm[i].ncorr;j++)
				printf("CORR %5d %5s %5d %5.1lf %5.1lf %5.1lf %5.1lf %5.1lf %5.1lf %5.1lf %5.1lf %5.1lf %5.1lf %5d\n", 
				j+1, parm[i].corr[j].atomtype, parm[i].corr[j].pid + 1, 
                                parm[i].corr[j].bl, parm[i].corr[j].blf, parm[i].corr[j].cba, parm[i].corr[j].cbaf,
				parm[i].corr[j].ba, parm[i].corr[j].baf, parm[i].corr[j].ctor, parm[i].corr[j].tor, 
				parm[i].corr[j].improper, parm[i].corr[j].ps, parm[i].corr[j].type); 
			printf("\n");
		}
		printf("WEIGHTS\n");
		printf("WT for BL: %9.2lf\n", wt.BL);
		printf("WT for BLF: %9.2lf\n", wt.BLF);
		printf("WT for BA: %9.2lf\n", wt.BA);
		printf("WT for BAF: %9.2lf\n", wt.BAF);
		printf("WT for X: %9.2lf\n", wt.X);
		printf("WT for X3: %9.2lf\n", wt.X3);
		printf("WT for BA_CTR: %9.2lf\n", wt.BA_CTR);
		printf("WT for TOR_CTR: %9.2lf\n", wt.TOR_CTR);
		printf("WT for IMPROPER: %9.2lf\n", wt.IMPROPER);
		printf("WT for GROUP: %9.2lf\n", wt.GROUP);
		printf("WT for EQUTYPE: %9.2lf\n", wt.EQUTYPE);
		printf("\nDEFAULT VALUES\n");
		printf("DEFAULT for BL: %9.2lf\n", dv.BL);
		printf("DEFAULT for BLF: %9.2lf\n", dv.BLF);
		printf("DEFAULT for BA: %9.2lf\n", dv.BA);
		printf("DEFAULT for BAF: %9.2lf\n", dv.BAF);
		printf("DEFAULT for BA_CTR: %9.2lf\n", dv.BA_CTR);
		printf("DEFAULT for TOR: %9.2lf\n", dv.TOR);
		printf("DEFAULT for TOR_CTR: %9.2lf\n", dv.TOR_CTR);
		printf("DEFAULT for FRACT1: %9.2lf\n", dv.FRACT1);
		printf("DEFAULT for FRACT2: %9.2lf\n", dv.FRACT2);
		printf("THRESHOLD for BA: %9.2lf\n", THRESHOLD_BA);

		printf("\n-- Finish reading parmchk parm file %s --\n\n", filename);
	}
}


void read_parmchk_parm_v2(char *filename)
{
	FILE *fp;
	char line[MAXCHAR];
	char atomtype[10];
	char equaname[10];
	char corrname[10];
	int  i,j, k;
	int  group;
	int  improper;
	int  ncorr;
	int  nequa;
	int  equtype;
	int  atomicnum;
	double bl,blf,ba,baf, cba, cbaf, tor, ctor, itor, ps;
	double mass;
	if ((fp = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open file %s in read_parmchk_parm_v2(), exit\n", filename);
		exit(1);
	}

	for (;;) {
		if (fgets(line, MAXCHAR, fp) == NULL)
			break;
		if (strncmp("PARM", &line[0], 4) == 0) {
			equtype = 0;
			sscanf(&line[4], "%s%d%d%lf%d%d", atomtype, &improper, &group, &mass, &equtype, &atomicnum) ;
			strcpy(parm[parmnum].atomtype, atomtype);
			parm[parmnum].improper = improper;
			parm[parmnum].group    = group;
			parm[parmnum].mass     = mass;
			parm[parmnum].equtype  = equtype;
			parm[parmnum].atomicnum  = atomicnum;
/*	Initialization */
/*	for the sake of simplicity in the coding, an atom type equals to and corresponds to itself*/
/*      Equal atom types must also be corresponding atom types */
                        strcpy(parm[parmnum].equa[0].atomtype, atomtype);
                        strcpy(parm[parmnum].corr[0].atomtype, atomtype);
                        parm[parmnum].corr[0].bl   = 0;
                        parm[parmnum].corr[0].blf  = 0;
                        parm[parmnum].corr[0].ba   = 0;
                        parm[parmnum].corr[0].baf  = 0;
                        parm[parmnum].corr[0].cba  = 0;
                        parm[parmnum].corr[0].cbaf = 0;
                        parm[parmnum].corr[0].tor  = 0;
                        parm[parmnum].corr[0].ctor = 0;
                        parm[parmnum].corr[0].improper = 0;
                        parm[parmnum].corr[0].ps    = 0;
                        parm[parmnum].corr[0].type = 0;
                        parm[parmnum].corr[0].pid  = parmnum;
			nequa = 1;
			ncorr = 1;
			parm[parmnum].nequa = 1;
			parm[parmnum].ncorr = 1;
			parmnum++;
			if (parmnum >= maxparmnum) {
				maxparmnum += MAXPARM;
				parm = (PARM *) realloc(parm, sizeof(PARM) * maxparmnum);
				if (parm == NULL) {
					fprintf(stdout, "memory allocation error for *parm\n");
					exit(1);
				}
			}
		}
		if (strncmp("EQUA", &line[0], 4) == 0) {
			sscanf(&line[4], "%s",  equaname); 
			strcpy(parm[parmnum-1].equa[nequa].atomtype, equaname);
			parm[parmnum-1].nequa ++;
			nequa ++;
			if(parm[parmnum-1].nequa > MAXEQUA) {
				fprintf(stdout, "Too many equal atom types for Atom Type %d (%s)\n", parmnum, parm[parmnum-1].atomtype);
				exit(1);
			}
			strcpy(parm[parmnum-1].corr[ncorr].atomtype, equaname);
			parm[parmnum-1].corr[ncorr].bl  =0;
			parm[parmnum-1].corr[ncorr].blf =0;
			parm[parmnum-1].corr[ncorr].ba  =0;
			parm[parmnum-1].corr[ncorr].baf =0;
			parm[parmnum-1].corr[ncorr].cba =0;
			parm[parmnum-1].corr[ncorr].cbaf=0;
			parm[parmnum-1].corr[ncorr].tor =0;
			parm[parmnum-1].corr[ncorr].ctor=0;
			parm[parmnum-1].corr[ncorr].ps  =0;
			parm[parmnum-1].corr[ncorr].improper=0;
			parm[parmnum-1].corr[ncorr].pid = -1;
			parm[parmnum-1].corr[ncorr].type = 1;
			parm[parmnum-1].ncorr ++;
			ncorr ++;
			if(parm[parmnum-1].ncorr > MAXCORR) {
				fprintf(stdout, "Too many corresponding atom types for Atom Type %d (%s)\n", parmnum, parm[parmnum-1].atomtype);
				exit(1);
			}
		}
		if (strncmp("CORR", &line[0], 4) == 0) {
			bl  = 0;
			blf = 0;
			ba  = 0;
			baf = 0;
			cba = 0;
			cbaf= 0;
			tor = 0;
			ctor= 0;
			itor= 0;
			ps=0;
			sscanf(&line[4], "%s%lf%lf%lf%lf%lf%lf%lf%lf%lf",  
              corrname, &bl, &blf, &cba, &cbaf, &ba, &baf, &ctor, &tor, &ps); 
			strcpy(parm[parmnum-1].corr[ncorr].atomtype, corrname);
			parm[parmnum-1].corr[ncorr].bl  =bl;
			parm[parmnum-1].corr[ncorr].blf =blf;
			parm[parmnum-1].corr[ncorr].ba  =ba;
			parm[parmnum-1].corr[ncorr].baf =baf;
			parm[parmnum-1].corr[ncorr].cba =cba;
			parm[parmnum-1].corr[ncorr].cbaf=cbaf;
			parm[parmnum-1].corr[ncorr].tor =tor;
			parm[parmnum-1].corr[ncorr].ctor=ctor;
			parm[parmnum-1].corr[ncorr].improper=ps; /*assign general score to improper parameter*/
			parm[parmnum-1].corr[ncorr].ps=ps;
			parm[parmnum-1].corr[ncorr].pid = -1;
			parm[parmnum-1].corr[ncorr].type = 2;
			parm[parmnum-1].ncorr ++;
			ncorr ++;
			if(parm[parmnum-1].ncorr > MAXCORR) {
				fprintf(stdout, "Too many corresponding atom types for Atom Type %d (%s)\n", parmnum, parm[parmnum-1].atomtype);
				exit(1);
			}
		}
	}
	fclose(fp);
	for(i=parmnum2;i<parmnum;i++) 
		for(j=0;j<parm[i].ncorr;j++) {
			if(parm[i].corr[j].bl   < 0) parm[i].corr[j].bl  = dv.BL; 
			if(parm[i].corr[j].blf  < 0) parm[i].corr[j].blf = dv.BLF; 
			if(parm[i].corr[j].ba   < 0) parm[i].corr[j].ba  = dv.BA ;
			if(parm[i].corr[j].baf  < 0) parm[i].corr[j].baf = dv.BAF;
			if(parm[i].corr[j].cba  < 0) parm[i].corr[j].cba = dv.BA_CTR;
			if(parm[i].corr[j].cbaf < 0) parm[i].corr[j].cbaf= dv.BAF_CTR;
			if(parm[i].corr[j].tor  < 0) parm[i].corr[j].tor = dv.TOR;
			if(parm[i].corr[j].ctor < 0) parm[i].corr[j].ctor= dv.TOR_CTR;
			parm[i].corr[j].ctor = parm[i].corr[j].ctor * dv.FRACT1  + parm[i].corr[j].ps * dv.FRACT2; 
		}
/*	now finding the pid of equal/corresponding atom types */
	for(i=parmnum2;i<parmnum;i++) {
		for(j=0;j<parm[i].nequa;j++) 
			for(k=0;k<parmnum;k++) 
				if(strcmp(parm[i].equa[j].atomtype, parm[k].atomtype) == 0) {
					parm[i].equa[j].pid = k;
					break;
				}
		for(j=0;j<parm[i].ncorr;j++) 
			for(k=0;k<parmnum;k++) 
				if(strcmp(parm[i].corr[j].atomtype, parm[k].atomtype) == 0) {
					parm[i].corr[j].pid = k;
					break;
				}
	}
	if(debug == 1) {
		printf("-- Begin read additional parmchk parameter file: %s\n\n", filename);
		for(i=parmnum2;i<parmnum;i++) {
			printf("PARM %5d %5s %5d %5d %9.4lf %5d\n", i+1, parm[i].atomtype, parm[i].improper, parm[i].group, parm[i].mass, parm[i].equtype);
			for(j=0;j<parm[i].nequa;j++)
				printf("EQUA %5d %5s %5d \n", 
				j+1, parm[i].equa[j].atomtype, parm[i].equa[j].pid + 1); 
			for(j=0;j<parm[i].ncorr;j++)
				printf("CORR %5d %5s %5d %5.1lf %5.1lf %5.1lf %5.1lf %5.1lf %5.1lf %5.1lf %5.1lf %5.1lf %5.1lf %5d\n", 
				j+1, parm[i].corr[j].atomtype, parm[i].corr[j].pid + 1, 
                                parm[i].corr[j].bl, parm[i].corr[j].blf, parm[i].corr[j].cba, parm[i].corr[j].cbaf,
				parm[i].corr[j].ba, parm[i].corr[j].baf, parm[i].corr[j].ctor, parm[i].corr[j].tor, 
				parm[i].corr[j].improper, parm[i].corr[j].ps, parm[i].corr[j].type); 
			printf("\n");
		}
		printf("-- Finish reading additional parmchk parameter file: %s\n\n", filename);
	}
}
void assign_parmid(void)
{
	int i, j;
	for (i = 0; i < atomnum; i++) {
		parmid[i] = -1;
		for (j = 0; j < parmnum; j++)
			if (strcmp(parm[j].atomtype, atom[i].ambername) == 0) {
				parmid[i] = j;
				if(parm[j].improper == 1) 
					atom[i].improper = 1;
				break;
			}
                if(parmid[i] < 0) {
                        fprintf(stderr, "Atom type of %s does not exist in PARMCHK.DAT\n", atom[i].ambername);
                        exit(1);
                }
	}
}

int assign_parmid_v2(char *at) {
int pid;
	pid = -1;
	for (j = 0; j < parmnum; j++) {
		if (strcmp(parm[j].atomtype, at) == 0) {
			pid = j;
			break;
		}
	}
        if(pid < 0) {
                fprintf(stderr, "Atom type of %s does not exist in PARMCHK.DAT\n", at);
		exit(1);
	}
	return pid;
}

double check_bond_length(char *name1, char *name2) {
double bl = 0;
int i;
	for (i = 0; i < bondparmnum; i++) {
		if ((	strcmp(bondparm[i].name1, name1) == 0 &&
                        strcmp(bondparm[i].name2, name2) == 0)) {
			bl = bondparm[i].length;
			break;
		}
		if ((	strcmp(bondparm[i].name1, name2) == 0 &&
                     	strcmp(bondparm[i].name2, name1) == 0)) {
			bl = bondparm[i].length;
			break;
		}
	}
	return bl;
}
int empbond(char *tmpc1, char *tmpc2,  char *name1, char *name2, int id1, int id2, double input_bl) {
	double bl = 0.0;
	double force;
	double kref = -1;
	double rref = -1;
	double delta_r1, delta_r2;
	int flag = 0;
	int pid1 = -1, pid2 = -1;

	if(iread_blbaparm == 0) {
		iread_blbaparm = 1;
/*allow other force field parameter sets, such as ff99SB and ff14SB to use PARM_BLBA_GAFF.DAT*/
		if(gaff_set == 2) 
    			build_dat_path(blba_parmfile, "PARM_BLBA_GAFF2.DAT",
    				sizeof blba_parmfile, 0);
		else
    			build_dat_path(blba_parmfile, "PARM_BLBA_GAFF.DAT",
    				sizeof blba_parmfile, 0);
		read_blba_parm (blba_parmfile, blf_parm, &nblf_parm, baf_parm); 
		if(debug == 1) {
			for(i=0;i<nblf_parm;i++)
				printf("BL %5d %5s %5d %5s %5d %9.4lf %9.4lf\n", i+1, blf_parm[i].elem1, 
					blf_parm[i].id1, blf_parm[i].elem2, blf_parm[i].id2,
					blf_parm[i].refbondlength, blf_parm[i].bfkai);
			for(i=0;i<120;i++)
				printf("BA %5d %5s %5d %9.4lf %9.4lf\n", i+1, baf_parm[i].elem, 
					baf_parm[i].id, baf_parm[i].anglec, baf_parm[i].anglez); 
			printf("PC %9.4lf\n", blf_exp_const);
		}
	}
	bl = input_bl;
	if(bl <= 0) {
		for (i = 0; i < bondparmnum; i++) {
			if ((strcmp(bondparm[i].name1, tmpc1) == 0 &&
                     	     strcmp(bondparm[i].name2, tmpc2) == 0)) {
				bl = bondparm[i].length;
				break;
			}
			if ((strcmp(bondparm[i].name1, tmpc2) == 0 &&
                     	     strcmp(bondparm[i].name2, tmpc1) == 0)) {
				bl = bondparm[i].length;
				break;
			}
		}
	}
	if( bl <= 0) return 0; /*later we may use an empirical rule to predict bond length*/

	for(i=0;i<nblf_parm;i++) {
		if(blf_parm[i].id1 == id1 && blf_parm[i].id2 == id2) {
			kref = blf_parm[i].bfkai;		
			rref = blf_parm[i].refbondlength;
			break;
		} 
		if(blf_parm[i].id1 == id2 && blf_parm[i].id2 == id1) {
			kref = blf_parm[i].bfkai;		
			rref = blf_parm[i].refbondlength;
			break;
		} 
	}
	if( kref < 0 || rref < 0) {
		for(i=0;i<nblf_parm;i++) 
			if(blf_parm[i].id1 == id1 && blf_parm[i].id2 == id1) {
				pid1 = i;
				break;
			}
		for(i=0;i<nblf_parm;i++) 
			if(blf_parm[i].id1 == id2 && blf_parm[i].id2 == id2) {
				pid2 = i;
				break;
			}
		if(pid1 < 0 || pid2 < 0) return 0; 
		if(rref < 0) {
			rref = 0.5 * (blf_parm[pid1].refbondlength + blf_parm[pid2].refbondlength); /*need to verify this assumption*/
			flag = 1;
		}
                delta_r1 = fabs(rref - blf_parm[pid1].refbondlength); 
                delta_r2 = fabs(rref - blf_parm[pid2].refbondlength); 
                kref =blf_parm[pid1].bfkai * delta_r1 ;
                kref+=blf_parm[pid2].bfkai * delta_r2 ;
                kref/=(delta_r1 + delta_r2);
	}
	if( kref <= 0) return 0;

	force = exp(kref - blf_exp_const * log(bl));		
        fprintf(fpout, "%-2s-%-2s%8.2lf%8.3lf", name1, name2, force, bl);
        fprintf(fpout, "       calculated for %2s-%2s", name1, name2); 
	if(flag == 0)
		fprintf(fpout, "\n");
	else
		fprintf(fpout, "Warning: reference bond length of %s-%s (%8.4lf) was calculated using those of %s-%s (%7.4lf) and %s-%s (%7.4lf)\n", 
			blf_parm[pid1].elem1, blf_parm[pid2].elem1, rref,
			blf_parm[pid1].elem1, blf_parm[pid1].elem2, blf_parm[pid1].refbondlength,
			blf_parm[pid2].elem1, blf_parm[pid2].elem2, blf_parm[pid2].refbondlength);

        strcpy(bondparm[bondparmnum].name1, name1);
        strcpy(bondparm[bondparmnum].name2, name2);
	bondparm[bondparmnum].force = force;
	bondparm[bondparmnum].length = bl;
        bondparmnum++;
        if (bondparmnum >= maxbondparm) {
        	maxbondparm += MAX_FF_BOND;
                bondparm =
               		(BOND_FF *) realloc(bondparm,
                        sizeof(BOND_FF) * maxbondparm);
                if (bondparm == NULL) {
                	fprintf(stdout,
                        	"memory allocation error for *bondparm\n");
                        exit(1);
                }
         }
	return 1;
}

int empangle(char *tmpc1, char *tmpc2, char *tmpc3, char *name1,
	     char *name2, char *name3, int id1, int id2, int id3, 
	     double input_bl1, double input_bl2, double input_angle, int flag)
{
	int num1 = -1, num2 = -1;
	double bondlength1 = 0.0;
	double bondlength2 = 0.0;
	double cparm, dparm, zparm1, zparm2;
	double angle;
	double force;

	if(iread_blbaparm == 0) {
		iread_blbaparm = 1;
/*allow other force field parameter sets, such as ff99SB and ff14SB to use PARM_BLBA_GAFF.DAT*/
		if(gaff_set == 2) 
    			build_dat_path(blba_parmfile, "PARM_BLBA_GAFF2.DAT",
    				sizeof blba_parmfile, 0);
		else
    			build_dat_path(blba_parmfile, "PARM_BLBA_GAFF.DAT",
    				sizeof blba_parmfile, 0);
		read_blba_parm (blba_parmfile, blf_parm, &nblf_parm, baf_parm); 
		if(debug == 1) {
			for(i=0;i<nblf_parm;i++)
				printf("BL %5d %5s %5d %5s %5d %9.4lf %9.4lf\n", i+1, blf_parm[i].elem1, 
					blf_parm[i].id1, blf_parm[i].elem2, blf_parm[i].id2,
					blf_parm[i].refbondlength, blf_parm[i].bfkai);
			for(i=0;i<120;i++)
				printf("BA %5d %5s %5d %9.4lf %9.4lf\n", i+1, baf_parm[i].elem, 
					baf_parm[i].id, baf_parm[i].anglec, baf_parm[i].anglez); 
			printf("PC %9.4lf\n", blf_exp_const);
		}
	}

	if(input_angle <= 0) {
		for (i = 0; i < angleparmnum; i++)
		    if (strcmp(angleparm[i].name1, tmpc1) == 0 && 
		    	strcmp(angleparm[i].name2, tmpc2) == 0 && 
		    	strcmp(angleparm[i].name3, tmpc1) == 0) {
				num1 = i;
				break;
			}
		if (num1 == -1)
			return 0;
		for (i = 0; i < angleparmnum; i++)
			if (strcmp(angleparm[i].name1, tmpc3) == 0 && 
		    	strcmp(angleparm[i].name2, tmpc2) == 0 && 
		    	strcmp(angleparm[i].name3, tmpc3) == 0) {
				num2 = i;
				break;
			}
		if (num2 == -1)
			return 0;

		angle = 0.5 * (angleparm[num1].angle + angleparm[num2].angle);
	}
	else
		angle = input_angle;
	if( angle <= 0) return 0;

	if(input_bl1 <= 0) 
		bondlength1=check_bond_length(tmpc1, tmpc2);
	else
		bondlength1 = input_bl1;
	if (bondlength1 == 0.0)
		return 0;

	if(input_bl2 <= 0) 
		bondlength2=check_bond_length(tmpc2, tmpc3);
	else
		bondlength2 = input_bl2;
	if (bondlength2 == 0.0)
		return 0;

	/* calculate the bond angle force  */
	zparm1 = baf_parm[id1].anglez;
	zparm2 = baf_parm[id3].anglez;
	cparm  = baf_parm[id2].anglec;
	dparm = (bondlength1 - bondlength2) * (bondlength1 - bondlength2);
	dparm =
		dparm / ((bondlength1 + bondlength2) * (bondlength1 + bondlength2));
	force =
		143.9 * zparm1 * cparm * zparm2 * exp(-2 * dparm) / (bondlength1 + bondlength2);
	force /= sqrt(angle * 3.1415926 / 180.0);
	if(flag == 0) {
		fprintf(fpout, "%-2s-%-2s-%-2s%9.3lf%12.3lf", name1, name2, name3, force, angle);
		fprintf(fpout, "   Calculated with empirical approach for %s-%s-%s\n", tmpc1, tmpc2, tmpc3);
		strcpy(angleparm[angleparmnum].name1, name1);
		strcpy(angleparm[angleparmnum].name2, name2);
		strcpy(angleparm[angleparmnum].name3, name3);
		angleparm[angleparmnum].angle = angle;
		angleparm[angleparmnum].force = force;
		angleparmnum++;
		if (angleparmnum >= maxangleparm) {
			maxangleparm += MAX_FF_ANGLE;
			angleparm =
				(ANGLE *) realloc(angleparm, sizeof(ANGLE) * maxangleparm);
			if (angleparm == NULL) {
				fprintf(stdout, "memory allocation error for *angleparm\n");
				exit(1);
			}
		}
	}
	if(flag == 1) {
		bestba.angle = angle;
		bestba.force = force;
		strcpy(bestba.name1, tmpc1);
		strcpy(bestba.name2, tmpc2);
		strcpy(bestba.name3, tmpc3);
		bestbaid = 999999;
	}
	return 1;
}

double equtype_penalty(int id1, int id2, int id3, int id4){
int tn1; 
int tn2;
int n1, n2, n3, n4;
	n1 = parm[id1].equtype;
	n2 = parm[id2].equtype;
	n3 = parm[id3].equtype;
	n4 = parm[id4].equtype;
	if(n1 == 0 && n2 == 0) return 0.0;

	tn1 = fabs(n1) + fabs(n2);
	tn2 = fabs(n3) + fabs(n4);
	if(tn2 == 3 && tn1 != 3) return wt.EQUTYPE;
	if(tn1 == 3 && tn2 != 3) return wt.EQUTYPE;
/*	e.g. for X-cd-cf-X, both X-cd-cd-X and X-cf-cf are possible replacement.  
 	To avoid ambiguity, the replacement of X-cf-cf-X now has a penalty score of 0.5 *wt.EQUTYPE 	
	Equtype: cc -1
                 cd -2
                 ce  1
                 cf  2
		 ...
*/ 
	if((n1 + n2) == 0) 
		if(n3 < 0 && n4 < 0) return 0.5*wt.EQUTYPE; 
	return 0.0;
}

void chk_atomtype(void)
{
/*	in this version, only polarizability parameter is replaced, as mass must be present in the PARMCHK.PARM file */
	int i, j;
	int suc;
	int pid;
	int pid2;
	char tmpc[5];
	fprintf(fpout, "%s\n", "Remark line goes here");
	fprintf(fpout, "%s\n", "MASS");
	for (i = 0; i < atomnum; i++) {
		suc = 0;
		pid=parmid[i];
		for (j = 0; j < atomtypenum; j++)
			if (strcmp(atom[i].ambername, atomtype[j].name) == 0) {
				suc = 1;
				if(allparm_flag == 1) {
					fprintf(fpout, "%-2s %-8.3lf   %8.3lf\n",
						atom[i].ambername, atomtype[j].mass,
						atomtype[j].pol);
				}
				break;
			}
/*	check equal atom types	*/	
/*	equal atom types	*/
		if (suc == 0 && parm[pid].nequa > 1) 
                        for(j=1;j<parm[pid].nequa;j++) {
                                pid2 = parm[pid].equa[j].pid;
                                strcpy(tmpc, parm[pid2].atomtype);
                                for (k = 0; k < atomtypenum2; k++)
                                        if (strcmp(tmpc, atomtype[k].name) == 0) {
                                                suc = 1;
                                                fprintf(fpout, "%-2s %-8.3lf   %8.3lf",
                                                                atom[i].ambername, parm[pid].mass,
                                                                atomtype[k].pol);
                                                fprintf(fpout, "               %s %-3s\n", "same as", atomtype[k].name);
                                                atomtype[atomtypenum] = atomtype[k];
                                                strcpy(atomtype[atomtypenum].name, atom[i].ambername);
                                                atomtypenum++;
                                                if (atomtypenum >= maxatomtype) {
                                                        maxatomtype += MAX_FF_ATOMTYPE;
                                                        atomtype =
                                                                (ATOMTYPE *) realloc(atomtype,
                                                                                                 sizeof(ATOMTYPE) *
                                                                                                 maxatomtype);
                                                        if (atomtype == NULL) {
                                                                fprintf(stdout,
                                                                                "memory allocation error for *atomtype\n");
                                                                exit(1);
                                                        }
                                                }
                                                break;
                                        }
                                if(suc == 1) break;
                        }
/*	corresponding atom types	*/
		if (suc == 0 && parm[pid].ncorr > 1) 
			for(j=1;j<parm[pid].ncorr;j++) {
				if(parm[pid].corr[j].type <= 1) continue;
				pid2 = parm[pid].corr[j].pid;
				strcpy(tmpc, parm[pid2].atomtype);
				for (k = 0; k < atomtypenum2; k++)
					if (strcmp(tmpc, atomtype[k].name) == 0) {
						suc = 1;
						fprintf(fpout, "%-2s %-8.3lf   %8.3lf",
								atom[i].ambername, parm[pid].mass,
								atomtype[k].pol);
						fprintf(fpout, "               %s %-3s\n", "same as", atomtype[k].name);
						atomtype[atomtypenum] = atomtype[k];
						strcpy(atomtype[atomtypenum].name, atom[i].ambername);
						atomtypenum++;
						if (atomtypenum >= maxatomtype) {
							maxatomtype += MAX_FF_ATOMTYPE;
							atomtype =
								(ATOMTYPE *) realloc(atomtype,
												 sizeof(ATOMTYPE) *
												 maxatomtype);
							if (atomtype == NULL) {
								fprintf(stdout,
										"memory allocation error for *atomtype\n");
								exit(1);
							}
						}
						break;
					}
				if(suc == 1) break;
			}
		if (suc == 0) {
			fprintf(fpout, "%-2s %-8.3lf   %8.3lf", atom[i].ambername, parm[pid].mass, 0.0);
			fprintf(fpout, "               %5s\n", "ATTN, no polarizability parameter");
			strcpy(atomtype[atomtypenum].name, atom[i].ambername);
			atomtypenum++;
			if (atomtypenum >= maxatomtype) {
				maxatomtype += MAX_FF_ATOMTYPE;
				atomtype =
					(ATOMTYPE *) realloc(atomtype,
										 sizeof(ATOMTYPE) * maxatomtype);
				if (atomtype == NULL) {
					fprintf(stdout,
							"memory allocation error for *atomtype\n");
					exit(1);
				}
			}
		}
	}
}

void chk_atomtype_v2 (void){
/*	in this version, only polarizability parameter is replaced, as mass must be present in the PARMCHK.PARM file */
	int i, j;
	int suc;
	int pid;
	int pid2;
	char tmpc[5];
	fprintf(fpout, "%s\n", "Remark line goes here");
	fprintf(fpout, "%s\n", "MASS");
	for (i = 0; i < r_atomtype_num; i++) {
		if(iformat == 4 && attn_opt == 2 && r_atomtype[i].attn == 0) {
			fprintf(fpout, "%-2s %-8.3lf   %8.3lf\n", r_atomtype[i].name,  
                                r_atomtype[i].mass, r_atomtype[i].pol); 
			continue;
		}
		suc = 0;
                pid = assign_parmid_v2(r_atomtype[i].name);
		for (j = 0; j < atomtypenum; j++)
			if (strcmp(r_atomtype[i].name, atomtype[j].name) == 0) {
				suc = 1;
					fprintf(fpout, "%-2s %-8.3lf   %8.3lf\n",
                                               r_atomtype[i].name,  parm[pid].mass, atomtype[j].pol); 
				break;
			}
/*	check equal atom types	*/	
/*	equal atom types	*/
		if (suc == 0 && parm[pid].nequa > 1) 
                        for(j=1;j<parm[pid].nequa;j++) {
                                pid2 = parm[pid].equa[j].pid;
                                strcpy(tmpc, parm[pid2].atomtype);
                                for (k = 0; k < atomtypenum2; k++)
                                        if (strcmp(tmpc, atomtype[k].name) == 0) {
                                                suc = 1;
                                                fprintf(fpout, "%-2s %-8.3lf   %8.3lf",
                                                                r_atomtype[i].name, parm[pid].mass,
                                                                atomtype[k].pol);
                                                fprintf(fpout, "               %s %-3s\n", "same as", atomtype[k].name);
                                                atomtype[atomtypenum] = atomtype[k];
                                                strcpy(atomtype[atomtypenum].name, r_atomtype[i].name);
                                                atomtypenum++;
                                                if (atomtypenum >= maxatomtype) {
                                                        maxatomtype += MAX_FF_ATOMTYPE;
                                                        atomtype =
                                                                (ATOMTYPE *) realloc(atomtype,
                                                                                                 sizeof(ATOMTYPE) *
                                                                                                 maxatomtype);
                                                        if (atomtype == NULL) {
                                                                fprintf(stdout,
                                                                                "memory allocation error for *atomtype\n");
                                                                exit(1);
                                                        }
                                                }
                                                break;
                                        }
                                if(suc == 1) break;
                        }
/*	corresponding atom types	*/
		if (suc == 0 && parm[pid].ncorr > 1) 
			for(j=1;j<parm[pid].ncorr;j++) {
				if(parm[pid].corr[j].type <= 1) continue;
				pid2 = parm[pid].corr[j].pid;
				strcpy(tmpc, parm[pid2].atomtype);
				for (k = 0; k < atomtypenum2; k++)
					if (strcmp(tmpc, atomtype[k].name) == 0) {
						suc = 1;
						fprintf(fpout, "%-2s %-8.3lf   %8.3lf",
								r_atomtype[i].name, parm[pid].mass,
								atomtype[k].pol);
						fprintf(fpout, "               %s %-3s\n", "same as", atomtype[k].name);
						atomtype[atomtypenum] = atomtype[k];
						strcpy(atomtype[atomtypenum].name, r_atomtype[i].name);
						atomtypenum++;
						if (atomtypenum >= maxatomtype) {
							maxatomtype += MAX_FF_ATOMTYPE;
							atomtype =
								(ATOMTYPE *) realloc(atomtype,
												 sizeof(ATOMTYPE) *
												 maxatomtype);
							if (atomtype == NULL) {
								fprintf(stdout,
										"memory allocation error for *atomtype\n");
								exit(1);
							}
						}
						break;
					}
				if(suc == 1) break;
			}
		if (suc == 0) {
			fprintf(fpout, "%-2s %-8.3lf   %8.3lf", r_atomtype[i].name, parm[pid].mass, 0.0);
			fprintf(fpout, "               %5s\n", "ATTN, no polarizability parameter");
			strcpy(atomtype[atomtypenum].name, r_atomtype[i].name);
			atomtypenum++;
			if (atomtypenum >= maxatomtype) {
				maxatomtype += MAX_FF_ATOMTYPE;
				atomtype =
					(ATOMTYPE *) realloc(atomtype,
										 sizeof(ATOMTYPE) * maxatomtype);
				if (atomtype == NULL) {
					fprintf(stdout,
							"memory allocation error for *atomtype\n");
					exit(1);
				}
			}
		}
	}

}

int vdw(char *at_name, char *corr_name, int nparm)
{
	int suc = 0;
	int i;
	for (i = 0; i < nparm; i++)
		if (strcmp(vdwparm[i].name, corr_name) == 0) {
			suc = 1;
			fprintf(fpout, "  %-2s%16.4lf%8.4lf",
					at_name, vdwparm[i].radius, vdwparm[i].pot);
			fprintf(fpout, "             %s %-3s\n", "same as",
					vdwparm[i].name);
			vdwparm[vdwparmnum] = vdwparm[i];
			strcpy(vdwparm[vdwparmnum].name, at_name);
			vdwparmnum++;
			if (vdwparmnum >= maxvdwparm) {
				maxvdwparm += MAX_FF_VDW;
				vdwparm =
					(VDW *) realloc(vdwparm, sizeof(VDW) * maxvdwparm);
				if (vdwparm == NULL) {
					fprintf(stdout,
							"memory allocation error for *vdwparm\n");
					exit(1);
				}
			}
			break;
		}
	return suc;
}

void chk_vdw(void)
{
	int i, j;
	int pid, pid2;
	int suc;
	fprintf(fpout, "\n%s\n", "NONBON");
	for (i = 0; i < atomnum; i++) {
		suc = 0;
		pid = parmid[i];
		for (j = 0; j < vdwparmnum; j++)
			if (strcmp(vdwparm[j].name, atom[i].ambername) == 0) {
				suc = 1;
				if(allparm_flag == 1) {
					fprintf(fpout, "  %-2s%16.4lf%8.4lf\n",
						vdwparm[j].name, vdwparm[j].radius, vdwparm[j].pot);
				}
				break;
			}
		if (suc == 0 && parm[pid].nequa > 1) 
			for(j=1; j< parm[pid].nequa; j++) {
				pid2 = parm[pid].equa[j].pid;	
				suc = vdw(atom[i].ambername, parm[pid2].atomtype, vdwparmnum2);
				if(suc == 1) break;
			}
		if (suc == 0 && parm[pid].ncorr > 1) 
			for(j=1; j< parm[pid].ncorr; j++) {
				if(parm[pid].corr[j].type <= 1) continue;
				pid2 = parm[pid].corr[j].pid;	
				suc = vdw(atom[i].ambername, parm[pid2].atomtype, vdwparmnum2);
				if(suc == 1) break;
			}
		if (suc == 0) {
			fprintf(fpout, "  %-2s%16.4lf%8.4lf", atom[i].ambername, 0.0,
					0.0);
			fprintf(fpout, "             %s\n", "ATTN, need revision");
			strcpy(vdwparm[vdwparmnum].name, atom[i].ambername);
			vdwparmnum++;
			if (vdwparmnum >= maxvdwparm) {
				maxvdwparm += MAX_FF_VDW;
				vdwparm =
					(VDW *) realloc(vdwparm, sizeof(VDW) * maxvdwparm);
				if (vdwparm == NULL) {
					fprintf(stdout,
							"memory allocation error for *vdwparm\n");
					exit(1);
				}
			}
		}
	}
	fprintf(fpout, "\n\n\n");
}

void chk_vdw_v2 (void){
	int i, j;
	int pid, pid2;
	int suc;
	fprintf(fpout, "\n%s\n", "NONBON");
	for (i = 0; i < r_vdwparm_num; i++) {
		if(iformat == 4 && attn_opt == 2 && r_vdwparm[i].attn == 0) {
			fprintf(fpout, "  %-2s%16.4lf%8.4lf\n",
				r_vdwparm[i].name, r_vdwparm[i].radius, r_vdwparm[i].pot);
			continue;
		}
		suc = 0;
                pid = assign_parmid_v2(r_vdwparm[i].name);
		for (j = 0; j < vdwparmnum; j++)
			if (strcmp(vdwparm[j].name, r_vdwparm[i].name) == 0) {
				suc = 1;
				if(allparm_flag == 1) {
					fprintf(fpout, "  %-2s%16.4lf%8.4lf\n",
						vdwparm[j].name, vdwparm[j].radius, vdwparm[j].pot);
				}
				break;
			}
		if (suc == 0 && parm[pid].nequa > 1) 
			for(j=1; j< parm[pid].nequa; j++) {
				pid2 = parm[pid].equa[j].pid;	
				suc = vdw(r_vdwparm[i].name, parm[pid2].atomtype, vdwparmnum2);
				if(suc == 1) break;
			}
		if (suc == 0 && parm[pid].ncorr > 1) 
			for(j=1; j< parm[pid].ncorr; j++) {
				if(parm[pid].corr[j].type <= 1) continue;
				pid2 = parm[pid].corr[j].pid;	
				suc = vdw(r_vdwparm[i].name, parm[pid2].atomtype, vdwparmnum2);
				if(suc == 1) break;
			}
		if (suc == 0) {
			fprintf(fpout, "  %-2s%16.4lf%8.4lf", r_vdwparm[i].name, 0.0,
					0.0);
			fprintf(fpout, "             %s\n", "ATTN, need revision");
			strcpy(vdwparm[vdwparmnum].name, r_vdwparm[i].name);
			vdwparmnum++;
			if (vdwparmnum >= maxvdwparm) {
				maxvdwparm += MAX_FF_VDW;
				vdwparm =
					(VDW *) realloc(vdwparm, sizeof(VDW) * maxvdwparm);
				if (vdwparm == NULL) {
					fprintf(stdout,
							"memory allocation error for *vdwparm\n");
					exit(1);
				}
			}
		}
	}
	fprintf(fpout, "\n\n\n");
}

int bond(char *at_name1, char *at_name2, char *corr_name1,
		 char *corr_name2, int nparm, int index)
{
	int suc = 0;
	int k;
	for (k = 0; k < nparm; k++)
		if ((strcmp(bondparm[k].name1, corr_name1) == 0 &&
		     strcmp(bondparm[k].name2, corr_name2) == 0) ||
                    (strcmp(bondparm[k].name1, corr_name2) == 0 &&
		     strcmp(bondparm[k].name2, corr_name1) == 0)) {
			suc = 1;
			if (index == 0 && allparm_flag == 1) {
				fprintf(fpout, "%-2s-%-2s%8.2lf%8.3lf\n", 
						at_name1, at_name2, 
						bondparm[k].force, bondparm[k].length);
			}
			if (index == 1) {
				bestblid = k;
				break;
			}
		}
	return suc;
}

void chk_bond(void)
{
	int i, j, m,n;
	int suc, suc2;
	int pid1, pid2;
	int pid3, pid4;
	int atomicnum1, atomicnum2;
	char tmpc[5];
	char tmpc1[5], tmpc2[5];
	char name1[5], name2[5];
	double score, score1, score2;
	fprintf(fpout, "\n%s\n", "BOND");

	for (i = 0; i < atomnum; i++)
		for (j = i + 1; j < atomnum; j++)
			if (atom[i].con[0] == j || atom[i].con[1] == j
				|| atom[i].con[2] == j || atom[i].con[3] == j
				|| atom[i].con[4] == j || atom[i].con[5] == j) {
				suc = 0;
				strcpy(tmpc1, atom[i].ambername);
				strcpy(tmpc2, atom[j].ambername);
                                if(strcmp(tmpc1, tmpc2) > 0)  {
                                        strcpy(tmpc, tmpc2);
                                        strcpy(tmpc2, tmpc1);
                                        strcpy(tmpc1, tmpc);
                                }
				strcpy(name1, tmpc1);
				strcpy(name2, tmpc2);

				suc = bond(name1, name2, tmpc1, tmpc2, bondparmnum, 0);
				if(suc == 1) continue;
/*	for equal atom types	*/
                                if(suc == 0) {
					bestblid = -1;
                                        pid1 = parmid[i];
                                        pid2 = parmid[j];
                                        for(m=0;m<parm[pid1].nequa; m++) {
						if(suc == 1) break;
                                                pid3 = parm[pid1].equa[m].pid;
                                                strcpy(tmpc1, parm[pid3].atomtype);
                                                for(n=0;n<parm[pid2].nequa; n++) {
							if(suc == 1) break;
                                                        if(m==0 && n== 0) continue;  
                                                        pid4 = parm[pid2].equa[n].pid;
                                                        strcpy(tmpc2, parm[pid4].atomtype);
                                                        suc2 = bond(name1, name2, tmpc1, tmpc2, bondparmnum2, 1);
							if(suc2 == 1) {
								bestscore = 0;
								suc = 1;
							}
                                        	}
					}
				}
/*	for corresponding atom types	*/
				if(suc == 0) {
					pid1 = parmid[i];
					pid2 = parmid[j];
					bestblid = -1;
					bestscore = INITSCORE;	
					for(m=0;m<parm[pid1].ncorr; m++) {
						pid3 = parm[pid1].corr[m].pid;
						strcpy(tmpc1, parm[pid3].atomtype);	
						score1= parm[pid1].corr[m].bl * wt.BL +
						        parm[pid1].corr[m].blf* wt.BLF;
						for(n=0;n<parm[pid2].ncorr; n++) {
							if(parm[pid1].corr[m].type <= 1 && 
                                                           parm[pid2].corr[n].type <= 1) 
								continue;	
							pid4 = parm[pid2].corr[n].pid;
							strcpy(tmpc2, parm[pid4].atomtype);	
							score2= parm[pid2].corr[n].bl * wt.BL +
						        	parm[pid2].corr[n].blf* wt.BLF;
							score = score1 + score2;
							if(parm[pid3].group != parm[pid4].group) score += wt.GROUP;
							equtype_penalty_score = equtype_penalty(pid1, pid2, pid3, pid4);	
							score += equtype_penalty_score;
							if(score < bestscore) {
								suc2= bond(name1, name2, tmpc1, tmpc2, bondparmnum2, 1);
								if(suc2== 1) {
									bestscore = score;
									suc = 1;
								}
							}
						}
					}
				}
				if(suc==1 && bestblid >= 0) {
                               		fprintf(fpout, "%-2s-%-2s%8.2lf%8.3lf", name1,
                                               		name2, bondparm[bestblid].force, bondparm[bestblid].length);
                               		fprintf(fpout, "       same as %2s-%2s, penalty score=%5.1lf\n",
							bondparm[bestblid].name1, bondparm[bestblid].name2, bestscore);
                               		bondparm[bondparmnum] = bondparm[bestblid]; 
                               		strcpy(bondparm[bondparmnum].name1, name1);
                               		strcpy(bondparm[bondparmnum].name2, name2);
                               		bondparmnum++;
                               		if (bondparmnum >= maxbondparm) {
                                       		maxbondparm += MAX_FF_BOND;
                                       		bondparm =
                                               		(BOND_FF *) realloc(bondparm,
                                                                                       sizeof(BOND_FF) * maxbondparm);
                                       		if (bondparm == NULL) {
                                               		fprintf(stdout,
                                                      	         	"memory allocation error for *bondparm\n");
                                               		exit(1);
                                       		}
                               		}
					bestbaid = -1;
					continue;
                        	}
			
				if(suc == 0) {
/*direct calculation*/
					strcpy(tmpc1, atom[i].ambername);
					strcpy(tmpc2, atom[j].ambername);
                                        pid1 = parmid[i];
                                        pid2 = parmid[j];
					atomicnum1 = parm[pid1].atomicnum;
					atomicnum2 = parm[pid2].atomicnum;
					suc = empbond(tmpc1, tmpc2, name1, name2, atomicnum1, atomicnum2, 0);
				}
			
				if(suc == 0) {
					fprintf(fpout, "%-2s-%-2s%8.2lf%8.3lf", name1, name2, 0.0, 0.0);
					fprintf(fpout, "       %s\n", "ATTN, need revision");
					strcpy(bondparm[bondparmnum].name1, name1);
					strcpy(bondparm[bondparmnum].name2, name2);
					bondparmnum++;
					if (bondparmnum >= maxbondparm) {
						maxbondparm += MAX_FF_BOND;
						bondparm =
							(BOND_FF *) realloc(bondparm,
												sizeof(BOND_FF) *
												maxbondparm);
						if (bondparm == NULL) {
							fprintf(stdout,
									"memory allocation error for *bondparm\n");
							exit(1);
						}
					}
				}
		}
}

void chk_bond_v2 (void){
	int i, j, m,n;
	int suc, suc2;
	int pid1, pid2;
	int pid3, pid4;
	int atomicnum1, atomicnum2;
	char tmpc[5];
	char tmpc1[5], tmpc2[5];
	char name1[5], name2[5];
	double score, score1, score2;
	fprintf(fpout, "\n%s\n", "BOND");

	for (i = 0; i < r_bondparm_num; i++)  {
		if(iformat == 4 && attn_opt == 2 && r_bondparm[i].attn == 0) {
			fprintf(fpout, "%-2s-%-2s%8.2lf%8.3lf\n", 
				r_bondparm[i].name1, r_bondparm[i].name2, 
				r_bondparm[i].force, r_bondparm[i].length); 
			continue;
		}
		suc = 0;
		strcpy(tmpc1, r_bondparm[i].name1);
		strcpy(tmpc2, r_bondparm[i].name2);
                if(strcmp(tmpc1, tmpc2) > 0)  {
                	strcpy(tmpc, tmpc2);
                        strcpy(tmpc2, tmpc1);
                        strcpy(tmpc1, tmpc);
                }
		strcpy(name1, tmpc1);
		strcpy(name2, tmpc2);
               	pid1 = assign_parmid_v2(r_bondparm[i].name1);
               	pid2 = assign_parmid_v2(r_bondparm[i].name2);

		suc = bond(name1, name2, tmpc1, tmpc2, bondparmnum, 0);
		if(suc == 1) continue;

/*	for equal atom types	*/
                if(suc == 0) {
			bestblid = -1;
                        for(m=0;m<parm[pid1].nequa; m++) {
				if(suc == 1) break;
                                pid3 = parm[pid1].equa[m].pid;
                                strcpy(tmpc1, parm[pid3].atomtype);
                                for(n=0;n<parm[pid2].nequa; n++) {
					if(suc == 1) break;
                                        if(m==0 && n== 0) continue;  
                                        pid4 = parm[pid2].equa[n].pid;
                                        strcpy(tmpc2, parm[pid4].atomtype);
                                        suc2 = bond(name1, name2, tmpc1, tmpc2, bondparmnum2, 1);
					if(suc2 == 1) {
						bestscore = 0;
						suc = 1;
					}
                           	}
			}
		}

		if(suc == 0 && fc_opt == 2 && r_bondparm[i].length > 0) {
			atomicnum1 = parm[pid1].atomicnum;
			atomicnum2 = parm[pid2].atomicnum;
			empbond(r_bondparm[i].name1, r_bondparm[i].name2, name1, name2, atomicnum1, atomicnum2, r_bondparm[i].length);
			continue;
		}

/*	for corresponding atom types	*/
		if(suc == 0) {
			bestblid = -1;
			bestscore = INITSCORE;	
			for(m=0;m<parm[pid1].ncorr; m++) {
				pid3 = parm[pid1].corr[m].pid;
				strcpy(tmpc1, parm[pid3].atomtype);	
				score1= parm[pid1].corr[m].bl * wt.BL +
				        parm[pid1].corr[m].blf* wt.BLF;
				for(n=0;n<parm[pid2].ncorr; n++) {
					if(parm[pid1].corr[m].type <= 1 && 
                                           parm[pid2].corr[n].type <= 1) 
						continue;	
					pid4 = parm[pid2].corr[n].pid;
					strcpy(tmpc2, parm[pid4].atomtype);	
					score2= parm[pid2].corr[n].bl * wt.BL +
				        	parm[pid2].corr[n].blf* wt.BLF;
					score = score1 + score2;
					if(parm[pid3].group != parm[pid4].group) score += wt.GROUP;
					equtype_penalty_score = equtype_penalty(pid1, pid2, pid3, pid4);	
					score += equtype_penalty_score;
					if(score < bestscore) {
						suc2= bond(name1, name2, tmpc1, tmpc2, bondparmnum2, 1);
						if(suc2== 1) {
							bestscore = score;
							suc = 1;
						}
					}
				}
			}
		}
		if(suc==1 && bestblid >= 0) {
            		fprintf(fpout, "%-2s-%-2s%8.2lf%8.3lf", name1,
                              		name2, bondparm[bestblid].force, bondparm[bestblid].length);
                    	fprintf(fpout, "       same as %2s-%2s, penalty score=%5.1lf\n",
					bondparm[bestblid].name1, bondparm[bestblid].name2, bestscore);
                             		bondparm[bondparmnum] = bondparm[bestblid]; 
                               		strcpy(bondparm[bondparmnum].name1, name1);
                               		strcpy(bondparm[bondparmnum].name2, name2);
                        bondparmnum++;
                        if (bondparmnum >= maxbondparm) {
                         	maxbondparm += MAX_FF_BOND;
                             	bondparm =
                               		(BOND_FF *) realloc(bondparm,
                                        	sizeof(BOND_FF) * maxbondparm);
                                      	if (bondparm == NULL) {
                                       		fprintf(stdout,
                                               		"memory allocation error for *bondparm\n");
                                       		exit(1);
                                       	}
                        } 
			bestbaid = -1;
			continue;
                }
		if(suc == 0) {
			atomicnum1 = parm[pid1].atomicnum;
			atomicnum2 = parm[pid2].atomicnum;
			suc=empbond(r_bondparm[i].name1, r_bondparm[i].name2, name1, name2, atomicnum1, atomicnum2, r_bondparm[i].length);
		}
		if(suc == 0) {
			fprintf(fpout, "%-2s-%-2s%8.2lf%8.3lf", name1, name2, 0.0, 0.0);
			fprintf(fpout, "       %s\n", "ATTN, need revision");
			strcpy(bondparm[bondparmnum].name1, name1);
			strcpy(bondparm[bondparmnum].name2, name2);
			bondparmnum++;
			if (bondparmnum >= maxbondparm) {
				maxbondparm += MAX_FF_BOND;
				bondparm =
					(BOND_FF *) realloc(bondparm,
							sizeof(BOND_FF) * maxbondparm);
				if (bondparm == NULL) {
					fprintf(stdout,
							"memory allocation error for *bondparm\n");
					exit(1);
				}
			}
		}
	}
}

int angle(char *at_name1, char *at_name2, char *at_name3, char *corr_name1,
		  char *corr_name2, char *corr_name3, int nparm, int index)
{
	int suc = 0;
	int l=0;
	for (l = 0; l < nparm; l++)
		if ((strcmp(angleparm[l].name1, corr_name1) == 0 && 
                     strcmp(angleparm[l].name2, corr_name2) == 0 &&
                     strcmp(angleparm[l].name3, corr_name3) == 0) ||
		    (strcmp(angleparm[l].name3, corr_name1) == 0 && 
                     strcmp(angleparm[l].name2, corr_name2) == 0 &&
                     strcmp(angleparm[l].name1, corr_name3) == 0)) {
			suc = 1;
			if (allparm_flag == 1 && index == 0) {
				fprintf(fpout,
						"%-2s-%-2s-%-2s%9.3lf%12.3lf\n",
						at_name1, at_name2, at_name3,
						angleparm[l].force, angleparm[l].angle);
			}
			if (index == 1) {
				bestba = angleparm[l];
				bestbaid = l;
				break;
			}
		}
	return suc;
}



void chk_angle(void)
{

	int i, j, k, m, n, o;
	int suc, suc2;
	int pid1, pid2, pid3;
	int pid4, pid5, pid6;
	int atomicnum1, atomicnum2, atomicnum3;
	char tmpc[5];
	char tmpc1[5], tmpc2[5], tmpc3[5];
	char tmpc4[5], tmpc5[5], tmpc6[5];
	char name1[5], name2[5], name3[5];
	double score, score1, score2, score3;
	fprintf(fpout, "\n%s\n", "ANGLE");

	/* NB: non-standard indentation in next four lines; for readability  */
	for (i = 0; i < atomnum; i++) {
		for (j = 0; j < atomnum; j++) {
			for (k = 0; k < atomnum; k++) {
				if (i != k) {
					if (atom[i].con[0] == j || atom[i].con[1] == j
						|| atom[i].con[2] == j || atom[i].con[3] == j
						|| atom[i].con[4] == j || atom[i].con[5] == j) {
						if (atom[j].con[0] == k || atom[j].con[1] == k
							|| atom[j].con[2] == k || atom[j].con[3] == k
							|| atom[j].con[4] == k
							|| atom[j].con[5] == k) {
							suc = 0;
							strcpy(tmpc1, atom[i].ambername);
							strcpy(tmpc2, atom[j].ambername);
							strcpy(tmpc3, atom[k].ambername);
                                			if(strcmp(tmpc1, tmpc3) > 0)  {
                                        			strcpy(tmpc, tmpc3);
                                        			strcpy(tmpc3, tmpc1);
                                        			strcpy(tmpc1, tmpc);
                                			}
							strcpy(name1, tmpc1);
							strcpy(name2, tmpc2);
							strcpy(name3, tmpc3);
							suc = angle(name1, name2, name3, tmpc1, tmpc2, tmpc3, angleparmnum, 0);
							if(suc == 1) continue;

/* for equal atom types */
							if(suc == 0) {
								bestbaid = -1;
                                        			pid1 = parmid[j]; /* it is j not i, as we do not want to replace atom type for central atom */
                                        			pid2 = parmid[i]; /* it is i not j */
                                        			pid3 = parmid[k];
                                                                for(m=0;m<parm[pid1].nequa; m++) {
									if(suc == 1) break;
                                                                        pid4 = parm[pid1].equa[m].pid;
                                                                        strcpy(tmpc4, parm[pid4].atomtype);
                                                                        for(n=0;n<parm[pid2].nequa; n++) {
										if(suc == 1) break;
                                                                                pid5 = parm[pid2].equa[n].pid;
                                                                                strcpy(tmpc5, parm[pid5].atomtype);
                                                                                for(o=0;o<parm[pid3].nequa; o++) {
											if(suc == 1) break;
                                                                                        if(m== 0 && n==0 && o==0) continue;
                                                                                        pid6 = parm[pid3].corr[o].pid;
                                                                                        strcpy(tmpc6, parm[pid6].atomtype);
                                                                                        suc2 = angle(name1, name2, name3, tmpc5, tmpc4, tmpc6, angleparmnum2, 1);
											if(suc2 == 1) {
												bestscore = 0;
												suc = 1;
											}
                                                                                }
                                                                        }
                                                                }
							}
/*	for corresponding atom types */
                                			if(suc == 0) {
                                                                pid1 = parmid[i];
                                                                pid2 = parmid[j];
                                                                pid3 = parmid[k];
								bestbaid = -1;
                                        			bestscore = INITSCORE;
                                        			for(m=0;m<parm[pid1].ncorr; m++) {
                                                			pid4 = parm[pid1].corr[m].pid;
                                                			strcpy(tmpc4, parm[pid4].atomtype);
                                                			score1= parm[pid1].corr[m].ba * wt.BA +
                                                        			parm[pid1].corr[m].baf* wt.BAF;
                                                			for(n=0;n<parm[pid2].ncorr; n++) {
                                                        			pid5 = parm[pid2].corr[n].pid;
                                                        			strcpy(tmpc5, parm[pid5].atomtype);
                                                        			score2= parm[pid2].corr[n].cba * wt.BA +
                                                                			parm[pid2].corr[n].cbaf* wt.BAF;
										score2 *= wt.BA_CTR;
                                                				for(o=0;o<parm[pid3].ncorr; o++) {
											if(parm[pid1].corr[m].type <= 1 && 
                                                           				   parm[pid2].corr[n].type <= 1 &&
                                                           				   parm[pid3].corr[o].type <= 1) 
												continue;	
                                                        				pid6 = parm[pid3].corr[o].pid;
                                                        				strcpy(tmpc6, parm[pid6].atomtype);
                                                        				score3= parm[pid3].corr[o].ba * wt.BA +
                                                                				parm[pid3].corr[o].baf* wt.BAF;
											score = score1 + score2 + score3 + wt.GROUP;
											if(parm[pid4].group == parm[pid5].group && 
											   parm[pid4].group == parm[pid6].group)
                                                                                        	score -= wt.GROUP;
                                                        				if(score < bestscore && score <=THRESHOLD_BA) {
                                                                				suc2= angle(name1, name2, name3, tmpc4, tmpc5, tmpc6, angleparmnum2, 1);
                                                                				if(suc2== 1) {
													bestscore = score;
													suc = 1;
												}
                                                        				}
										}
                                                			}
								}
							}
							if(suc == 1 && bestbaid >= 0) {
                                				fprintf(fpout, "%-2s-%-2s-%-2s%9.3lf%12.3lf",
                                                				name1, name2, name3,
										bestba.force, bestba.angle);
                                				fprintf(fpout, "   same as %-2s-%-2s-%-2s, penalty score=%5.1lf\n",
                                                				bestba.name1, bestba.name2, bestba.name3, bestscore); 
                                				angleparm[angleparmnum] = bestba;
                                				strcpy(angleparm[angleparmnum].name1, name1);
                                				strcpy(angleparm[angleparmnum].name2, name2);
                                				strcpy(angleparm[angleparmnum].name3, name3);
                                				angleparmnum++;
                                				if (angleparmnum >= maxangleparm) {
                                        				maxangleparm += MAX_FF_ANGLE;
                                        				angleparm = (ANGLE *) realloc(angleparm,
                                                                                 	sizeof(ANGLE) * maxangleparm);
                                        				if (angleparm == NULL) {
                                                				fprintf(stdout,
                                                               				"memory allocation error for *angleparm\n");
                                                				exit(1);
                                        				}
								}
								bestbaid = -1;
								continue;
                                        		}
/* from here estimate the bond angle parameters with empirical method for corresponding names */
							if (suc == 0) {
								atomicnum1 = parm[pid1].atomicnum;
								atomicnum2 = parm[pid2].atomicnum;
								atomicnum3 = parm[pid3].atomicnum;
								suc = empangle(name1, name2, name3, name1, name2, name3,
											atomicnum1, atomicnum2, atomicnum3, 0, 0, 0, 0);
								if(suc == 1) continue;
							}

/* for equal atom types*/
                                                        if (suc == 0) {
								bestbaid = -1;
                                                                pid1 = parmid[j]; /* it is j, not i*/
                                                                pid2 = parmid[i]; /* it is i, not j*/
                                                                pid3 = parmid[k];
                                                                for(m=0;m<parm[pid1].nequa; m++) {
									if(suc == 1) break;
                                                                        pid4 = parm[pid1].equa[m].pid;
                                                                        strcpy(tmpc4, parm[pid4].atomtype);
                                                                        for(n=0;n<parm[pid2].nequa; n++) {
										if(suc == 1) break;
                                                                                pid5 = parm[pid2].equa[n].pid;
                                                                                strcpy(tmpc5, parm[pid5].atomtype);
                                                                                for(o=0;o<parm[pid3].nequa; o++) {
											if(suc == 1) break;
                                                                                        if(m==0 && n== 0 && o==0) continue;
                                                                                        pid6 = parm[pid3].equa[o].pid;
                                                                                        strcpy(tmpc6, parm[pid6].atomtype);
											atomicnum1 = parm[pid4].atomicnum;
											atomicnum2 = parm[pid5].atomicnum;
											atomicnum3 = parm[pid6].atomicnum;
                                                                                        suc2= empangle(tmpc5, tmpc4, tmpc6, name1, name2, name3,
													atomicnum2, atomicnum1, atomicnum3, 0, 0, 0, 1);
											if(suc2== 1) {
												bestscore = 0;
												suc = 1;
											}
                                                                                }
                                                                        }
                                                                }
                                                        }
/*	for corresponding atom types*/
							if (suc == 0) {
                                                                pid1 = parmid[i];
                                                                pid2 = parmid[j];
                                                                pid3 = parmid[k];
                                                                bestscore = INITSCORE;
								bestbaid  = -1;
                                                                for(m=0;m<parm[pid1].ncorr; m++) {
                                                                        pid4 = parm[pid1].corr[m].pid;
                                                                        strcpy(tmpc4, parm[pid4].atomtype);
                                                                        score1= parm[pid1].corr[m].ba * wt.BA +
                                                                                parm[pid1].corr[m].baf* wt.BAF;
                                                                        for(n=0;n<parm[pid2].ncorr; n++) {
                                                                                pid5 = parm[pid2].corr[n].pid;
                                                                                strcpy(tmpc5, parm[pid5].atomtype);
                                                                                score2= parm[pid2].corr[n].cba * wt.BA +
                                                                                        parm[pid2].corr[n].cbaf* wt.BAF;
                                                                                score2 *= wt.BA_CTR; 
                                                                                for(o=0;o<parm[pid3].ncorr; o++) {
											if(parm[pid1].corr[m].type <= 1 && 
                                                           				   parm[pid2].corr[n].type <= 1 &&
                                                           				   parm[pid3].corr[o].type <= 1) 
												continue;	
                                                                                        pid6 = parm[pid3].corr[o].pid;
                                                                                        strcpy(tmpc6, parm[pid6].atomtype);
                                                                                        score3= parm[pid3].corr[o].ba * wt.BA +
                                                                                                parm[pid3].corr[o].baf* wt.BAF;
                                                                                        score = score1 + score2 + score3 + wt.GROUP;
											if(parm[pid4].group == parm[pid5].group &&
                                                                                           parm[pid4].group == parm[pid6].group)
                                                                                        	score -= wt.GROUP;
                                                                                        if(score < bestscore) {
												atomicnum1 = parm[pid4].atomicnum;
												atomicnum2 = parm[pid5].atomicnum;
												atomicnum3 = parm[pid6].atomicnum;
												suc2= empangle(tmpc4, tmpc5, tmpc6, name1, name2, name3,
													atomicnum1, atomicnum2, atomicnum3, 0,0,0, 1);
                                                                                                if(suc2== 1) {
													bestscore = score;
													suc = 1;
												}
                                                                                        }       
                                                                                }       
                                                                        }       
								}
							}
							if(suc==1 && bestbaid >= 0) {
                                				fprintf(fpout, "%-2s-%-2s-%-2s%9.3lf%12.3lf",
                                                				name1, name2, name3,
										bestba.force, bestba.angle);
                                				fprintf(fpout, "   Calculated using %2s-%2s-%2s, penalty score=%5.1lf\n",
                                                				bestba.name1, bestba.name2, bestba.name3, bestscore); 
                                				angleparm[angleparmnum] = bestba;
                                				strcpy(angleparm[angleparmnum].name1, name1);
                                				strcpy(angleparm[angleparmnum].name2, name2);
                                				strcpy(angleparm[angleparmnum].name3, name3);
                                				angleparmnum++;
                                				if (angleparmnum >= maxangleparm) {
                                        				maxangleparm += MAX_FF_ANGLE;
                                        				angleparm = (ANGLE *) realloc(angleparm,
                                                                                 	sizeof(ANGLE) * maxangleparm);
                                        				if (angleparm == NULL) {
                                                				fprintf(stdout,
                                                               				"memory allocation error for *angleparm\n");
                                                				exit(1);
                                        				}
								}
								bestbaid = -1;
								continue;
							}
							else {
								fprintf(fpout,
										"%-2s-%-2s-%-2s%9.3lf  %10.3lf",
										name1, name2, name3, 0.0, 0.0);
								fprintf(fpout, "   %s\n",
										"ATTN, need revision");
								strcpy(angleparm[angleparmnum].name1,
									   name1);
								strcpy(angleparm[angleparmnum].name2,
									   name2);
								strcpy(angleparm[angleparmnum].name3,
									   name3);
								angleparmnum++;
								if (angleparmnum >= maxangleparm) {
									maxangleparm += MAX_FF_ANGLE;
									angleparm = (ANGLE *) realloc(angleparm, sizeof(ANGLE) * maxangleparm);
									if (angleparm == NULL) {
										fprintf(stdout,
												"memory allocation error for *angleparm\n");
										exit(1);
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void chk_angle_v2 (void){
	int i, m, n, o;
	int suc, suc2;
	int pid1, pid2, pid3;
	int pid4, pid5, pid6;
	int atomicnum1, atomicnum2, atomicnum3;
	char tmpc[5];
	char tmpc1[5], tmpc2[5], tmpc3[5];
	char tmpc4[5], tmpc5[5], tmpc6[5];
	char name1[5], name2[5], name3[5];
	double score, score1, score2, score3;
	double bl1, bl2;
	fprintf(fpout, "\n%s\n", "ANGLE");

	/* NB: non-standard indentation in next four lines; for readability  */
	for (i = 0; i < r_angleparm_num; i++) {
		if(iformat == 4 && attn_opt == 2 && r_angleparm[i].attn == 0) {
			fprintf(fpout,
					"%-2s-%-2s-%-2s%9.3lf%12.3lf\n",
					r_angleparm[i].name1, r_angleparm[i].name2, r_angleparm[i].name3,
					r_angleparm[i].force, r_angleparm[i].angle); 
			continue;
		}
		suc = 0;
		strcpy(tmpc1, r_angleparm[i].name1);	
		strcpy(tmpc2, r_angleparm[i].name2);	
		strcpy(tmpc3, r_angleparm[i].name3);	
		strcpy(name1, tmpc1);
		strcpy(name2, tmpc2);
		strcpy(name3, tmpc3);
		suc = angle(name1, name2, name3, tmpc1, tmpc2, tmpc3, angleparmnum, 0);
		if(suc == 1) continue;

/* for equal atom types */
		bl1 = check_bond_length(r_angleparm[i].name1, r_angleparm[i].name2);
		bl2 = check_bond_length(r_angleparm[i].name2, r_angleparm[i].name3);

		if(suc == 0) {
               		pid1 = assign_parmid_v2(r_angleparm[i].name2);
               		pid2 = assign_parmid_v2(r_angleparm[i].name1);
               		pid3 = assign_parmid_v2(r_angleparm[i].name3);
			bestbaid = -1;
                        for(m=0;m<parm[pid1].nequa; m++) {
				if(suc == 1) break;
                                pid4 = parm[pid1].equa[m].pid;
                                strcpy(tmpc4, parm[pid4].atomtype);
                                for(n=0;n<parm[pid2].nequa; n++) {
					if(suc == 1) break;
                                        pid5 = parm[pid2].equa[n].pid;
                                        strcpy(tmpc5, parm[pid5].atomtype);
                                        for(o=0;o<parm[pid3].nequa; o++) {
						if(suc == 1) break;
                                                if(m== 0 && n==0 && o==0) continue;
                                                pid6 = parm[pid3].corr[o].pid;
                                                strcpy(tmpc6, parm[pid6].atomtype);
						if(bl1 <= 0) bl1=check_bond_length(tmpc4, tmpc5);
						if(bl2 <= 0) bl1=check_bond_length(tmpc4, tmpc6);
                                                suc2 = angle(name1, name2, name3, tmpc5, tmpc4, tmpc6, angleparmnum2, 1);
						if(suc2 == 1) {
							bestscore = 0;
							suc = 1;
						}
                                       	 }
                                 }
                        }
	       	}
/*do empirical calculations if fc_opt == 2*/
		if(suc == 0 && fc_opt == 2 && r_angleparm[i].angle > 0) {
               		pid1 = assign_parmid_v2(r_angleparm[i].name1);
               		pid2 = assign_parmid_v2(r_angleparm[i].name2);
               		pid3 = assign_parmid_v2(r_angleparm[i].name3);
			atomicnum1 = parm[pid1].atomicnum;
			atomicnum2 = parm[pid2].atomicnum;
			atomicnum3 = parm[pid3].atomicnum;
			suc= empangle(r_angleparm[i].name1, r_angleparm[i].name2, r_angleparm[i].name3, name1, name2, name3,
				atomicnum1, atomicnum2, atomicnum3, bl1, bl2, r_angleparm[i].angle, 0);
			if(suc == 1) continue;

		}
/* for corresponding atom types */
            	if(suc == 0) {
               		pid1 = assign_parmid_v2(r_angleparm[i].name1);
               		pid2 = assign_parmid_v2(r_angleparm[i].name2);
               		pid3 = assign_parmid_v2(r_angleparm[i].name3);
			bestbaid = -1;
                       	bestscore = INITSCORE;
                       	for(m=0;m<parm[pid1].ncorr; m++) {
                       		pid4 = parm[pid1].corr[m].pid;
                               	strcpy(tmpc4, parm[pid4].atomtype);
                               	score1= parm[pid1].corr[m].ba * wt.BA + 
					parm[pid1].corr[m].baf* wt.BAF;
                                for(n=0;n<parm[pid2].ncorr; n++) {
                                	pid5 = parm[pid2].corr[n].pid;
                                       	strcpy(tmpc5, parm[pid5].atomtype);
                                       	score2= parm[pid2].corr[n].cba * wt.BA +
                                       		parm[pid2].corr[n].cbaf* wt.BAF;
					score2 *= wt.BA_CTR;
                                       	for(o=0;o<parm[pid3].ncorr; o++) {
						if(parm[pid1].corr[m].type <= 1 && 
                                       	   		parm[pid2].corr[n].type <= 1 &&
                                                        parm[pid3].corr[o].type <= 1) 
						continue;	
                                               	pid6 = parm[pid3].corr[o].pid;
                                              	strcpy(tmpc6, parm[pid6].atomtype);
                                              	score3= parm[pid3].corr[o].ba * wt.BA +
                                               		parm[pid3].corr[o].baf* wt.BAF;
						score = score1 + score2 + score3 + wt.GROUP;
						if(parm[pid4].group == parm[pid5].group && 
					   	   parm[pid4].group == parm[pid6].group)
                                                  	score -= wt.GROUP;
                                               	if(score < bestscore && score <=THRESHOLD_BA) {
                                               		suc2= angle(name1, name2, name3, tmpc4, tmpc5, tmpc6, angleparmnum2, 1);
                                                       	if(suc2== 1) {
								bestscore = score;
								suc = 1;
							}
                                               	}
					}
                                 }
			}
		}
		if(suc == 1 && bestbaid >= 0) {
               		fprintf(fpout, "%-2s-%-2s-%-2s%9.3lf%12.3lf", name1, name2, name3,
				bestba.force, bestba.angle);
                        fprintf(fpout, "   same as %-2s-%-2s-%-2s, penalty score=%5.1lf\n",
                        	bestba.name1, bestba.name2, bestba.name3, bestscore); 
                        angleparm[angleparmnum] = bestba;
                        strcpy(angleparm[angleparmnum].name1, name1);
                        strcpy(angleparm[angleparmnum].name2, name2);
                        strcpy(angleparm[angleparmnum].name3, name3);
                        angleparmnum++;
                        if (angleparmnum >= maxangleparm) {
                        	maxangleparm += MAX_FF_ANGLE;
                               	angleparm = (ANGLE *) realloc(angleparm,
                                     	sizeof(ANGLE) * maxangleparm);
                            	if (angleparm == NULL) {
                               		fprintf(stdout,
                               			"memory allocation error for *angleparm\n");
                               		exit(1);
                               	}
			}
			bestbaid = -1;
			continue;
           	}
/* from here estimate the bond angle parameters with empirical method for corresponding names */
		if (suc == 0) {
			atomicnum1 = parm[pid1].atomicnum;
			atomicnum2 = parm[pid2].atomicnum;
			atomicnum3 = parm[pid3].atomicnum;
			suc = empangle(name1, name2, name3, name1, name2, name3,
				       atomicnum1, atomicnum2, atomicnum3, bl1, bl2, r_angleparm[i].angle, 0);
			if(suc == 1) continue;
		}

/* for equal atom types*/
                if (suc == 0) {
			bestbaid = -1;
               		pid1 = assign_parmid_v2(r_angleparm[i].name2);
               		pid2 = assign_parmid_v2(r_angleparm[i].name1);
               		pid3 = assign_parmid_v2(r_angleparm[i].name3);
                        for(m=0;m<parm[pid1].nequa; m++) {
				if(suc == 1) break;
                                pid4 = parm[pid1].equa[m].pid;
                                strcpy(tmpc4, parm[pid4].atomtype);
                                for(n=0;n<parm[pid2].nequa; n++) {
					if(suc == 1) break;
                                        pid5 = parm[pid2].equa[n].pid;
                                        strcpy(tmpc5, parm[pid5].atomtype);
                                        for(o=0;o<parm[pid3].nequa; o++) {
						if(suc == 1) break;
                                                if(m==0 && n== 0 && o==0) continue;
                                                pid6 = parm[pid3].equa[o].pid;
                                                strcpy(tmpc6, parm[pid6].atomtype);
						atomicnum1 = parm[pid4].atomicnum;
						atomicnum2 = parm[pid5].atomicnum;
						atomicnum3 = parm[pid6].atomicnum;
                                                suc2= empangle(tmpc5, tmpc4, tmpc6, name1, name2, name3,
						               atomicnum2, atomicnum1, atomicnum3, bl1, bl2, r_angleparm[i].angle, 1);
						if(suc2== 1) {
							bestscore = 0;
							suc = 1;
						}
                                       	}
                             	}
                      	}
              	}
/*	for corresponding atom types*/
		if (suc == 0) {
               		pid1 = assign_parmid_v2(r_angleparm[i].name1);
               		pid2 = assign_parmid_v2(r_angleparm[i].name2);
               		pid3 = assign_parmid_v2(r_angleparm[i].name3);
                        bestscore = INITSCORE;
			bestbaid  = -1;
                        for(m=0;m<parm[pid1].ncorr; m++) {
                        	pid4 = parm[pid1].corr[m].pid;
                                strcpy(tmpc4, parm[pid4].atomtype);
                               	score1= parm[pid1].corr[m].ba * wt.BA +
                                        parm[pid1].corr[m].baf* wt.BAF;
                                for(n=0;n<parm[pid2].ncorr; n++) {
                                	pid5 = parm[pid2].corr[n].pid;
                                        strcpy(tmpc5, parm[pid5].atomtype);
                                        score2= parm[pid2].corr[n].cba * wt.BA +
                                                parm[pid2].corr[n].cbaf* wt.BAF;
                                        score2 *= wt.BA_CTR; 
                                        for(o=0;o<parm[pid3].ncorr; o++) {
						if(parm[pid1].corr[m].type <= 1 && 
                                                 	parm[pid2].corr[n].type <= 1 &&
                                                        parm[pid3].corr[o].type <= 1) 
								continue;	
                                               	pid6 = parm[pid3].corr[o].pid;
                                               	strcpy(tmpc6, parm[pid6].atomtype);
                                                score3= parm[pid3].corr[o].ba * wt.BA +
                                                        parm[pid3].corr[o].baf* wt.BAF;
                                                score = score1 + score2 + score3 + wt.GROUP;
						if(parm[pid4].group == parm[pid5].group &&
                                                   parm[pid4].group == parm[pid6].group)
                                                     	score -= wt.GROUP;
                                                if(score < bestscore) {
							atomicnum1 = parm[pid4].atomicnum;
							atomicnum2 = parm[pid5].atomicnum;
							atomicnum3 = parm[pid6].atomicnum;
							suc2= empangle(tmpc4, tmpc5, tmpc6, name1, name2, name3,
									atomicnum1, atomicnum2, atomicnum3, bl1, bl2, r_angleparm[i].angle, 1);
                                                      	if(suc2== 1) {
								bestscore = score;
								suc = 1;
							}
                                              	}       
                                          }       
                                }      
			}
		}
		if(suc==1 && bestbaid >= 0) {
       			fprintf(fpout, "%-2s-%-2s-%-2s%9.3lf%12.3lf",
                       		name1, name2, name3,
				bestba.force, bestba.angle);
                       	fprintf(fpout, "   Calculated using %2s-%2s-%2s, penalty score=%5.1lf\n",
                       		bestba.name1, bestba.name2, bestba.name3, bestscore); 
                       		angleparm[angleparmnum] = bestba;
                               strcpy(angleparm[angleparmnum].name1, name1);
                               strcpy(angleparm[angleparmnum].name2, name2);
                               strcpy(angleparm[angleparmnum].name3, name3);
                               angleparmnum++;
                         if (angleparmnum >= maxangleparm) {
                         	maxangleparm += MAX_FF_ANGLE;
                               	angleparm = (ANGLE *) realloc(angleparm,
                               		sizeof(ANGLE) * maxangleparm);
                               	if (angleparm == NULL) {
                               		fprintf(stdout,
                               			"memory allocation error for *angleparm\n");
                                       	exit(1);
                               	}
			}
			bestbaid = -1;
			continue;
		}
		else {
			fprintf(fpout,
				"%-2s-%-2s-%-2s%9.3lf  %10.3lf",
				name1, name2, name3, 0.0, 0.0);
			fprintf(fpout, "   %s\n",
				"ATTN, need revision");
			strcpy(angleparm[angleparmnum].name1, name1);
			strcpy(angleparm[angleparmnum].name2, name2);
			strcpy(angleparm[angleparmnum].name3, name3);
			angleparmnum++;
			if (angleparmnum >= maxangleparm) {
				maxangleparm += MAX_FF_ANGLE;
				angleparm = (ANGLE *) realloc(angleparm, sizeof(ANGLE) * maxangleparm);
				if (angleparm == NULL) {
					fprintf(stdout, "memory allocation error for *angleparm\n");
					exit(1);
				}
			}
		}
	}

}

int torsion(char *at_name1, char *at_name2, char *at_name3, char *at_name4,
			char *corr_name1, char *corr_name2, char *corr_name3,
			char *corr_name4, int nparm, int index)
{
	int suc = 0;
	int m, n;
	for (m = 0; m < nparm; m++)
		if ((strcmp(torsionparm[m].name1, corr_name1) == 0 &&
                     strcmp(torsionparm[m].name2, corr_name2) == 0 &&
                     strcmp(torsionparm[m].name3, corr_name3) == 0 &&
                     strcmp(torsionparm[m].name4, corr_name4) == 0) ||
		    (strcmp(torsionparm[m].name4, corr_name1) == 0 &&
                     strcmp(torsionparm[m].name3, corr_name2) == 0 &&
                     strcmp(torsionparm[m].name2, corr_name3) == 0 &&
                     strcmp(torsionparm[m].name1, corr_name4) == 0)) {
			suc = 1;
			n = m;
			if (allparm_flag == 1 && index == 0) {
				while (torsionparm[n].fterm < 0) {
                        		fprintf(fpout,
                                		"%-2s-%-2s-%-2s-%-2s%4d%9.3lf%14.3lf%16.3lf\n",
                                        	at_name1, at_name2, at_name3, at_name4, 
						torsionparm[n].mul, torsionparm[n].force,
                                        	torsionparm[n].phase, torsionparm[n].fterm);
					n++;
				}
				fprintf(fpout,
						"%-2s-%-2s-%-2s-%-2s%4d%9.3lf%14.3lf%16.3lf\n",
						at_name1, at_name2, at_name3, at_name4, 
						torsionparm[n].mul, torsionparm[n].force,
						torsionparm[n].phase, torsionparm[n].fterm);
			}
			if (index == 1) {
				besttorid = m;
				break;
			}
		}
	return suc;
}

void print_torsion(int besttorid, char *name1, char *name2, char *name3, char *name4) {

	while (torsionparm[besttorid].fterm < 0) {
		fprintf(fpout, "%-2s-%-2s-%-2s-%-2s%4d%9.3lf%14.3lf%16.3lf",
                       	name1, name2, name3, name4, 
			torsionparm[besttorid].mul,   torsionparm[besttorid].force,
                       	torsionparm[besttorid].phase, torsionparm[besttorid].fterm);
        	fprintf(fpout, "      same as %-2s-%-2s-%-2s-%-2s\n",
                       	torsionparm[besttorid].name1,
                       	torsionparm[besttorid].name2,
                       	torsionparm[besttorid].name3,
                       	torsionparm[besttorid].name4);
        	torsionparm[torsionparmnum] = torsionparm[besttorid];
        	strcpy(torsionparm[torsionparmnum].name1, name1);
        	strcpy(torsionparm[torsionparmnum].name2, name2);
        	strcpy(torsionparm[torsionparmnum].name3, name3);
        	strcpy(torsionparm[torsionparmnum].name4, name4);
        	torsionparmnum++;
        	if (torsionparmnum >= maxtorsionparm) {
        		maxtorsionparm += MAX_FF_TORSION;
                	torsionparm = (TORSION *) realloc(torsionparm, sizeof(TORSION) * maxtorsionparm);
                	if (torsionparm == NULL) {
                		fprintf(stdout, "memory allocation error for *torsionparm\n");
                        	exit(1);
                	}
        	}
        	besttorid ++;
	}
	fprintf(fpout, "%-2s-%-2s-%-2s-%-2s%4d%9.3lf%14.3lf%16.3lf",
               	name1, name2, name3, name4, 
		torsionparm[besttorid].mul,   torsionparm[besttorid].force,
                torsionparm[besttorid].phase, torsionparm[besttorid].fterm);
	fprintf(fpout, "      same as %-2s-%-2s-%-2s-%-2s, penalty score=%5.1lf\n",
               	torsionparm[besttorid].name1,
               	torsionparm[besttorid].name2,
               	torsionparm[besttorid].name3,
               	torsionparm[besttorid].name4, bestscore);
	torsionparm[torsionparmnum] = torsionparm[besttorid];
	strcpy(torsionparm[torsionparmnum].name1, name1);
	strcpy(torsionparm[torsionparmnum].name2, name2);
	strcpy(torsionparm[torsionparmnum].name3, name3);
	strcpy(torsionparm[torsionparmnum].name4, name4);
	torsionparmnum++;
	if (torsionparmnum >= maxtorsionparm) {
        	maxtorsionparm += MAX_FF_TORSION;
                torsionparm = (TORSION *) realloc(torsionparm, sizeof(TORSION) * maxtorsionparm);
                if (torsionparm == NULL) {
                	fprintf(stdout, "memory allocation error for *torsionparm\n");
                         exit(1);
                }
         }
}

void chk_torsion(void)
{

	int i, j, k, l;
	int m, n, p, q;
	int pid1, pid2, pid3, pid4;
	int pid5, pid6, pid7, pid8;
	int suc, suc2;
	char tmpc[5];
	char tmpc1[5], tmpc2[5], tmpc3[5], tmpc4[5];
	char tmpc5[5], tmpc6[5], tmpc7[5], tmpc8[5];
	char name1[5], name2[5], name3[5], name4[5];
	double score, score1, score2, score3, score4;
	fprintf(fpout, "\n%s\n", "DIHE");

	/* NB: non-standard indentation in next four lines; for readability  */
	for (i = 0; i < atomnum; i++) {
		for (j = 0; j < atomnum; j++) {
			for (k = j + 1; k < atomnum; k++) {
				for (l = 0; l < atomnum; l++) {
					if (i != k && l != j) {
						if (atom[i].con[0] == j || atom[i].con[1] == j
							|| atom[i].con[2] == j || atom[i].con[3] == j
							|| atom[i].con[4] == j
							|| atom[i].con[5] == j) {
							if (atom[j].con[0] == k || atom[j].con[1] == k
								|| atom[j].con[2] == k
								|| atom[j].con[3] == k
								|| atom[j].con[4] == k
								|| atom[j].con[5] == k) {
								if (atom[l].con[0] == k
									|| atom[l].con[1] == k
									|| atom[l].con[2] == k
									|| atom[l].con[3] == k
									|| atom[l].con[4] == k
									|| atom[l].con[5] == k) {
									suc = 0;
									strcpy(tmpc1, atom[i].ambername);
									strcpy(tmpc2, atom[j].ambername);
									strcpy(tmpc3, atom[k].ambername);
									strcpy(tmpc4, atom[l].ambername);

                                					if(strcmp(tmpc2, tmpc3) > 0)  {
                                        					strcpy(tmpc, tmpc3);
                                        					strcpy(tmpc3, tmpc2);
                                        					strcpy(tmpc2, tmpc);
                                        					strcpy(tmpc, tmpc4);
                                        					strcpy(tmpc4, tmpc1);
                                        					strcpy(tmpc1, tmpc);
                                					} else if (strcmp(tmpc2, tmpc3) == 0) {
										if(strcmp(tmpc1, tmpc4) > 0) {
                                        						strcpy(tmpc, tmpc4);
                                        						strcpy(tmpc4, tmpc1);
                                        						strcpy(tmpc1, tmpc);
										}
									}
/*
printf("%d_%d_%d_%d\n", i+1,j+1,k+1,l+1);
fflush(stdout);
*/
									strcpy(name1, tmpc1);
									strcpy(name2, tmpc2);
									strcpy(name3, tmpc3);
									strcpy(name4, tmpc4);
/* Step 1 check if the special torsional parameter exists or not */
									suc = torsion(name1, name2, name3, name4,
											tmpc1, tmpc2, tmpc3, tmpc4, torsionparmnum, 0);
									if(suc == 1) continue;
/* Step 2 check special torsional parameters using equal atom types */
                                                                        if (suc == 0) {
                                                                        	besttorid = -1;
                                                                                pid1 = parmid[j]; /*central atoms in the first two layers of loops */
                                                                                pid2 = parmid[k];
                                                                                pid3 = parmid[i];
                                                                                pid4 = parmid[l];
                                                                                for(m=0; m<parm[pid1].nequa; m++) {
											if(suc == 1) break;
                                                                                        pid5 = parm[pid1].equa[m].pid;
                                                                                        strcpy(tmpc5, parm[pid5].atomtype);
                                                                                        for(n=0;n<parm[pid2].nequa; n++) {
												if(suc == 1) break;
                                                                                                pid6 = parm[pid2].equa[n].pid;
                                                                                                strcpy(tmpc6, parm[pid6].atomtype);
                                                                                                for(p=0;p<parm[pid3].nequa; p++) {
													if(suc == 1) break;
                                                                                                        pid7 = parm[pid3].equa[p].pid;
                                                                                                        strcpy(tmpc7, parm[pid7].atomtype);
                                                                                                        for(q=0;q<parm[pid4].nequa; q++) {
														if(suc == 1) break;
                                                                                                                if(m==0 && n==0 && p== 0 && q==0) continue;
                                                                                                                pid8 = parm[pid4].equa[q].pid;
                                                                                                                strcpy(tmpc8, parm[pid8].atomtype);
                                                                                                                suc2= torsion(name1, name2, name3, name4, tmpc7, tmpc5, tmpc6, tmpc8, torsionparmnum2, 1);
														if(suc2== 1) {
															bestscore = 0;
															suc = 1;
														}
                                                                                                        }
                                                                                                }
                                                                                        }
                                                                                }
										if(suc == 1 && besttorid >=0) {
											print_torsion(besttorid, name1, name2, name3, name4);
											besttorid = -1;
											continue;
										}
                                                                        }
									
/* Step 3. check general torsional angle terms*/
									if(suc == 0) {
										suc = torsion(name1, name2, name3, name4,
											"X", tmpc2, tmpc3, "X", torsionparmnum2, 0);
										if(suc == 1) continue;
									}

/* Step 4 check general torsional angle terms using equal atom types*/
                                                                        if (suc == 0) {
                                                                        	besttorid = -1;
                                                                                pid2 = parmid[j];
                                                                                pid3 = parmid[k];
                                                                                for(n=0;n<parm[pid2].nequa; n++) {
                                                                                        pid6 = parm[pid2].equa[n].pid;
                                                                                        strcpy(tmpc6, parm[pid6].atomtype);
                                                                                        for(p=0;p<parm[pid3].nequa; p++) {
												if(n==0 && p==0) continue;
                                                                                                pid7 = parm[pid3].equa[p].pid;
                                                                                                strcpy(tmpc7, parm[pid7].atomtype);
                                                                                                suc2= torsion(name1, name2, name3, name4, "X", tmpc6, tmpc7, "X", torsionparmnum2, 1);
												if(suc2== 1) {
													bestscore = 0;
													suc = 1;
												}
                                                                                        }
                                                                                }
										if(suc == 1 && besttorid >=0) {
											print_torsion(besttorid, name1, name2, name3, name4);
											besttorid = -1;
											continue;
										}
                                                                        }
/* Step 5. check special torsional parameters using corresponding atom types */
									if (suc == 0) {
										pid1 = parmid[i];	
										pid2 = parmid[j];	
										pid3 = parmid[k];	
										pid4 = parmid[l];	
                                                                		bestscore = INITSCORE;
										besttorid = -1;
										for(m=0; m<parm[pid1].ncorr; m++) {
											pid5 = parm[pid1].corr[m].pid;
                                                                        		strcpy(tmpc5, parm[pid5].atomtype);
                                                                        		score1= parm[pid1].corr[m].tor ; 
                                                                        		for(n=0;n<parm[pid2].ncorr; n++) {
												pid6 = parm[pid2].corr[n].pid;
                                                                        			strcpy(tmpc6, parm[pid6].atomtype);
                                                                        			score2= parm[pid2].corr[n].ctor * wt.TOR_CTR; 
                                                                        			for(p=0;p<parm[pid3].ncorr; p++) {
													pid7 = parm[pid3].corr[p].pid;
                                                                        				strcpy(tmpc7, parm[pid7].atomtype);
                                                                        				score3= parm[pid3].corr[p].ctor * wt.TOR_CTR; 
                                                                        				for(q=0;q<parm[pid4].ncorr; q++) {
														if(parm[pid1].corr[m].type <= 1 && 
                                                           				   			   parm[pid2].corr[n].type <= 1 &&
                                                           				   			   parm[pid3].corr[p].type <= 1 && 
                                                           				   			   parm[pid4].corr[q].type <= 1) 
															continue;	
														pid8 = parm[pid4].corr[q].pid;
                                                                        					strcpy(tmpc8, parm[pid8].atomtype);
                                                                        					score4= parm[pid4].corr[q].tor; 
														score = score1 + score2 + score3 + score4;
														score += wt.GROUP;
														equtype_penalty_score = equtype_penalty(pid2, pid3, pid6, pid7);	
														score += equtype_penalty_score;
														if(parm[pid5].group == parm[pid6].group && 
                                                                                                                   parm[pid5].group == parm[pid7].group &&
                                                                                                                   parm[pid5].group == parm[pid8].group)
															score -= wt.GROUP;
														if(score < bestscore) {
															suc2= torsion(name1, name2, name3, name4, tmpc5, tmpc6, tmpc7, tmpc8, torsionparmnum2, 1);
															if(suc2== 1) {
																bestscore = score;
																suc = 1;
															}
														}
													}
												}
											}
										}
										if(suc == 1 && besttorid >=0) {
											print_torsion(besttorid, name1, name2, name3, name4);
											besttorid = -1;
											continue;
										}
									}
/* Step 6. check general torsional parameters using corresponding atom types */
                                                                        if (suc == 0) {
                                                                                pid2 = parmid[j];
                                                                                pid3 = parmid[k];
                                                                                bestscore = INITSCORE;
										besttorid = -1;
                                                                                for(n=0;n<parm[pid2].ncorr; n++) {
                                                                                        pid6 = parm[pid2].corr[n].pid;
                                                                                        strcpy(tmpc6, parm[pid6].atomtype);
                                                                                        score2= parm[pid2].corr[n].ctor * wt.TOR_CTR;
                                                                                        for(p=0;p<parm[pid3].ncorr; p++) {
                                                                                                pid7 = parm[pid3].corr[p].pid;
												if(parm[pid2].corr[n].type <= 1 && 
                                                           				   	   parm[pid3].corr[p].type <= 1) 
													continue;	
                                                                                                strcpy(tmpc7, parm[pid7].atomtype);
                                                                                                score3= parm[pid3].corr[p].ctor * wt.TOR_CTR;
                                                                                                score = score2 + score3;
												equtype_penalty_score = equtype_penalty(pid2, pid3, pid6, pid7);	
												score += equtype_penalty_score;
												if(parm[pid6].group != parm[pid7].group) 
													score += wt.GROUP;
                                                                                                if(score < bestscore) {
                                                                                                        suc2= torsion(name1, name2, name3, name4, "X", tmpc6, tmpc7, "X", torsionparmnum2, 1);
													if(suc2== 1) {
														bestscore = score;
														suc = 1;
													}
                                                                                                }
                                                                                        }
                                                                                }
										if(suc == 1 && besttorid >=0) {
											print_torsion(besttorid, name1, name2, name3, name4);
											besttorid = -1;
											continue;
										}
                                                                        }

									if (suc == 0) {
										fprintf(fpout,
												"%-2s-%-2s-%-2s-%-2s%4d%9.3lf%14.3lf%16.3lf",
												name1, name2, name3, name4,
												1, 0.0, 0.0, 2.0);
										fprintf(fpout, "      %s\n",
												"ATTN, need revision");
										strcpy(torsionparm[torsionparmnum].
											   name1, name1);
										strcpy(torsionparm[torsionparmnum].
											   name2, name2);
										strcpy(torsionparm[torsionparmnum].
											   name3, name3);
										strcpy(torsionparm[torsionparmnum].
											   name4, name4);
										torsionparmnum++;
										if (torsionparmnum >=
											maxtorsionparm) {
											maxtorsionparm +=
												MAX_FF_TORSION;
											torsionparm = (TORSION *)
												realloc(torsionparm,
														sizeof(TORSION) *
														maxtorsionparm);
											if (torsionparm == NULL) {
												fprintf(stdout,
														"memory allocation error for *torsionparm\n");
												exit(1);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
void chk_torsion_v2 (void){
	int i;
	int m, n, p, q;
	int pid1, pid2, pid3, pid4;
	int pid5, pid6, pid7, pid8;
	int suc, suc2;
	char tmpc[5];
	char tmpc1[5], tmpc2[5], tmpc3[5], tmpc4[5];
	char tmpc5[5], tmpc6[5], tmpc7[5], tmpc8[5];
	char name1[5], name2[5], name3[5], name4[5];
	double score, score1, score2, score3, score4;
	fprintf(fpout, "\n%s\n", "DIHE");

	/* NB: non-standard indentation in next four lines; for readability  */
	for (i = 0; i < r_torsionparm_num; i++) {
		if(iformat == 4 && attn_opt == 2 && r_torsionparm[i].attn == 0) {
			fprintf(fpout, "%-2s-%-2s-%-2s-%-2s%4d%9.3lf%14.3lf%16.3lf\n",
				r_torsionparm[i].name1, r_torsionparm[i].name2, r_torsionparm[i].name3, r_torsionparm[i].name4, 
				r_torsionparm[i].mul,   r_torsionparm[i].force,
				r_torsionparm[i].phase, r_torsionparm[i].fterm);
			continue;
		}
		suc = 0;
                strcpy(tmpc1, r_torsionparm[i].name1);
                strcpy(tmpc2, r_torsionparm[i].name2);
                strcpy(tmpc3, r_torsionparm[i].name3);
                strcpy(tmpc4, r_torsionparm[i].name4);

		strcpy(name1, tmpc1);
		strcpy(name2, tmpc2);
		strcpy(name3, tmpc3);
		strcpy(name4, tmpc4);
/* Step 1 check if the special torsional parameter exists or not */
		suc = torsion(name1, name2, name3, name4,
				tmpc1, tmpc2, tmpc3, tmpc4, torsionparmnum, 0);
		if(suc == 1) continue;
/* Step 2 check special torsional parameters using equal atom types */
                if (suc == 0) {
                	besttorid = -1;
                       	pid1 = assign_parmid_v2(name2);
                       	pid2 = assign_parmid_v2(name3);
                        pid3 = assign_parmid_v2(name1); 
                        pid4 = assign_parmid_v2(name4);
                        for(m=0; m<parm[pid1].nequa; m++) {
				if(suc == 1) break;
                                pid5 = parm[pid1].equa[m].pid;
                                strcpy(tmpc5, parm[pid5].atomtype);
                                for(n=0;n<parm[pid2].nequa; n++) {
					if(suc == 1) break;
                                        pid6 = parm[pid2].equa[n].pid;
                                        strcpy(tmpc6, parm[pid6].atomtype);
                                        for(p=0;p<parm[pid3].nequa; p++) {
						if(suc == 1) break;
                                                pid7 = parm[pid3].equa[p].pid;
                                                strcpy(tmpc7, parm[pid7].atomtype);
                                                for(q=0;q<parm[pid4].nequa; q++) {
							if(suc == 1) break;
                                                        if(m==0 && n==0 && p== 0 && q==0) continue;
                                                        pid8 = parm[pid4].equa[q].pid;
                                                        strcpy(tmpc8, parm[pid8].atomtype);
                                                        suc2= torsion(name1, name2, name3, name4, tmpc7, tmpc5, tmpc6, tmpc8, torsionparmnum2, 1);
							if(suc2== 1) {
								bestscore = 0;
								suc = 1;
							}
                                                }
                                         }
                                }
                         }
			if(suc == 1 && besttorid >=0) {
				print_torsion(besttorid, name1, name2, name3, name4);
				besttorid = -1;
				continue;
			}
               	}
									
/* Step 3. check general torsional angle terms*/
		if(suc == 0) {
			suc = torsion(name1, name2, name3, name4,
				"X", tmpc2, tmpc3, "X", torsionparmnum2, 0);
			if(suc == 1) continue;
		}

/* Step 4 check general torsional angle terms using equal atom types*/
                if (suc == 0) {
                     	besttorid = -1;
                        pid2 = assign_parmid_v2(name2);
                        pid3 = assign_parmid_v2(name3); 
                        for(n=0;n<parm[pid2].nequa; n++) {
                        	pid6 = parm[pid2].equa[n].pid;
                                strcpy(tmpc6, parm[pid6].atomtype);
                                for(p=0;p<parm[pid3].nequa; p++) {
					if(n==0 && p==0) continue;
                                        pid7 = parm[pid3].equa[p].pid;
                                        strcpy(tmpc7, parm[pid7].atomtype);
                                        suc2= torsion(name1, name2, name3, name4, "X", tmpc6, tmpc7, "X", torsionparmnum2, 1);
					if(suc2== 1) {
						bestscore = 0;
						suc = 1;
					}
                                }
                        }
			if(suc == 1 && besttorid >=0) {
				print_torsion(besttorid, name1, name2, name3, name4);
				besttorid = -1;
				continue;
			}
              	}
/* Step 5. check special torsional parameters using corresponding atom types */
		if (suc == 0) {
                        pid1 = assign_parmid_v2(name1);
                        pid2 = assign_parmid_v2(name2);
                        pid3 = assign_parmid_v2(name3);
                        pid4 = assign_parmid_v2(name4);
               		bestscore = INITSCORE;
			besttorid = -1;
			for(m=0; m<parm[pid1].ncorr; m++) {
				pid5 = parm[pid1].corr[m].pid;
                       		strcpy(tmpc5, parm[pid5].atomtype);
                       		score1= parm[pid1].corr[m].tor ; 
                       		for(n=0;n<parm[pid2].ncorr; n++) {
					pid6 = parm[pid2].corr[n].pid;
                       			strcpy(tmpc6, parm[pid6].atomtype);
                                        score2= parm[pid2].corr[n].ctor * wt.TOR_CTR; 
                                        for(p=0;p<parm[pid3].ncorr; p++) {
						pid7 = parm[pid3].corr[p].pid;
                                               	strcpy(tmpc7, parm[pid7].atomtype);
                                               	score3= parm[pid3].corr[p].ctor * wt.TOR_CTR; 
                                               	for(q=0;q<parm[pid4].ncorr; q++) {
							if(parm[pid1].corr[m].type <= 1 && 
                                                       		parm[pid2].corr[n].type <= 1 &&
                                                              	parm[pid3].corr[p].type <= 1 && 
                                                              	parm[pid4].corr[q].type <= 1) 
								continue;	
							pid8 = parm[pid4].corr[q].pid;
                                                       	strcpy(tmpc8, parm[pid8].atomtype);
                                                       	score4= parm[pid4].corr[q].tor; 
							score = score1 + score2 + score3 + score4;
							score += wt.GROUP;
							equtype_penalty_score = equtype_penalty(pid2, pid3, pid6, pid7);	
							score += equtype_penalty_score;
							if(parm[pid5].group == parm[pid6].group && 
                                                        	parm[pid5].group == parm[pid7].group &&
                                                                parm[pid5].group == parm[pid8].group)
							score -= wt.GROUP;
							if(score < bestscore) {
								suc2= torsion(name1, name2, name3, name4, tmpc5, tmpc6, tmpc7, tmpc8, torsionparmnum2, 1);
								if(suc2== 1) {
									bestscore = score;
									suc = 1;
								}
							}
						}
					}
				}
			}
			if(suc == 1 && besttorid >=0) {
				print_torsion(besttorid, name1, name2, name3, name4);
				besttorid = -1;
				continue;
			}
		}
/* Step 6. check general torsional parameters using corresponding atom types */
                if (suc == 0) {
                        pid2 = assign_parmid_v2(name2);
                        pid3 = assign_parmid_v2(name3);
                        bestscore = INITSCORE;
			besttorid = -1;
                        for(n=0;n<parm[pid2].ncorr; n++) {
                        	pid6 = parm[pid2].corr[n].pid;
                                strcpy(tmpc6, parm[pid6].atomtype);
                                score2= parm[pid2].corr[n].ctor * wt.TOR_CTR;
                                for(p=0;p<parm[pid3].ncorr; p++) {
                                	pid7 = parm[pid3].corr[p].pid;
					if(parm[pid2].corr[n].type <= 1 && 
                                           parm[pid3].corr[p].type <= 1) 
						continue;	
                                         strcpy(tmpc7, parm[pid7].atomtype);
                                         score3= parm[pid3].corr[p].ctor * wt.TOR_CTR;
                                         score = score2 + score3;
					 equtype_penalty_score = equtype_penalty(pid2, pid3, pid6, pid7);	
					 score += equtype_penalty_score;
					 if(parm[pid6].group != parm[pid7].group) 
						score += wt.GROUP;
                                         if(score < bestscore) {
                                              	suc2= torsion(name1, name2, name3, name4, "X", tmpc6, tmpc7, "X", torsionparmnum2, 1);
						if(suc2== 1) {
							bestscore = score;
							suc = 1;
						}
                                         }
                                 }
                        }
			if(suc == 1 && besttorid >=0) {
				print_torsion(besttorid, name1, name2, name3, name4);
				besttorid = -1;
				continue;
			}
                 }

		if (suc == 0) {
			fprintf(fpout,
				"%-2s-%-2s-%-2s-%-2s%4d%9.3lf%14.3lf%16.3lf",
				name1, name2, name3, name4,
				1, 0.0, 0.0, 2.0);
			fprintf(fpout, "      %s\n",
				"ATTN, need revision");
			strcpy(torsionparm[torsionparmnum].name1, name1);
			strcpy(torsionparm[torsionparmnum].name2, name2);
			strcpy(torsionparm[torsionparmnum].name3, name3);
			strcpy(torsionparm[torsionparmnum].name4, name4);
			torsionparmnum++;
			if (torsionparmnum >=
				maxtorsionparm) {
				maxtorsionparm += MAX_FF_TORSION;
				torsionparm = (TORSION *)
					realloc(torsionparm, sizeof(TORSION) * maxtorsionparm);
					if (torsionparm == NULL) {
						fprintf(stdout,
							"memory allocation error for *torsionparm\n");
							exit(1);
						}
					}
		}
	}
}


void chk_improper(void)
{
	int i,j,k;
	int m,n,p,q;
	int suc;
	int tmpnum;
	int index1, index2, index4;
        char tmpc[5];
        char tmpc1[5], tmpc2[5], tmpc3[5], tmpc4[5];
        char name1[5], name2[5], name3[5], name4[5];
        int pid1, pid2, pid3, pid4;
        int pid5, pid6, pid7, pid8;
	double score, score1, score2, score3, score4;
	fprintf(fpout, "\n%s\n", "IMPROPER");

	for (i = 0; i < impropernum; i++) {
		suc = 0;
		strcpy(tmpc1, atom[improper[i].atid1].ambername);
		strcpy(tmpc2, atom[improper[i].atid2].ambername);
		strcpy(tmpc3, atom[improper[i].atid3].ambername);
		strcpy(tmpc4, atom[improper[i].atid4].ambername);

		if(strcmp(tmpc1, tmpc2) > 0)  {
			strcpy(tmpc, tmpc2);
			strcpy(tmpc2, tmpc1);
			strcpy(tmpc1, tmpc);
		}
		if(strcmp(tmpc1, tmpc4) > 0)  {
			strcpy(tmpc, tmpc4);
			strcpy(tmpc4, tmpc1);
			strcpy(tmpc1, tmpc);
		}
		if(strcmp(tmpc2, tmpc4) > 0)  {
			strcpy(tmpc, tmpc4);
			strcpy(tmpc4, tmpc2);
			strcpy(tmpc2, tmpc);
		}
		strcpy(name1, tmpc1);
		strcpy(name2, tmpc2);
		strcpy(name3, tmpc3);
		strcpy(name4, tmpc4);
/*	step 1 check directly	*/
		for (j = 0; j < improperparmnum; j++)
			if(strcmp(improperparm[j].name1, name1) == 0 && 
			   strcmp(improperparm[j].name2, name2) == 0 && 
			   strcmp(improperparm[j].name3, name3) == 0 && 
			   strcmp(improperparm[j].name4, name4) == 0) {
				suc = 1;
				if(allparm_flag == 1) {
                                	fprintf(fpout, "%-2s-%-2s-%-2s-%-2s %11.1lf%15.1lf%12.1lf\n", name1,
                                        	name2, name3, name4, improperparm[j].force, improperparm[j].phase,
                                                improperparm[j].fterm);
				}
				break;
			}
		if(suc == 1) continue;
/* 	Step 2 check special improper torsional parameters using equal atom types */
                if(suc == 0) {                                         
                        pid1 = parmid[improper[i].atid1]; 
                        pid2 = parmid[improper[i].atid2];
                        pid3 = parmid[improper[i].atid3];
                        pid4 = parmid[improper[i].atid4];
                	bestimproperid = -1;
                        for(m=0; m<parm[pid1].nequa; m++) {
				if(suc == 1) break;
                                pid5 = parm[pid1].equa[m].pid;
                                strcpy(tmpc1, parm[pid5].atomtype);
                                for(n=0;n<parm[pid2].nequa; n++) {
					if(suc == 1) break;
                                        pid6 = parm[pid2].equa[n].pid;
                                        strcpy(tmpc2, parm[pid6].atomtype);
                                        for(p=0;p<parm[pid3].nequa; p++) {
						if(suc == 1) break;
                                                pid7 = parm[pid3].equa[p].pid;
                                                strcpy(tmpc3, parm[pid7].atomtype);
                                                for(q=0;q<parm[pid4].nequa; q++) {
							if(suc == 1) break;
                                                        if(m==0 && n==0 && p== 0 && q==0) continue;
                                                        pid8 = parm[pid4].equa[q].pid;
                                                        strcpy(tmpc4, parm[pid8].atomtype);
                                                        for (k = 0; k < improperparmnum2; k++) {
                                                                if (improperparm[k].numX > 0) continue;
                                                                if(strcmp(improperparm[k].name1, tmpc1) == 0 &&
                                                                   strcmp(improperparm[k].name2, tmpc2) == 0 &&
                                                                   strcmp(improperparm[k].name3, tmpc3) == 0 &&
                                                                   strcmp(improperparm[k].name4, tmpc4) == 0) {
                                                                        bestimproperid = k;
                                                                        bestscore = 0;
                                                                        suc = 1;
                                                                        break;
                                                        	}
                                                	}
                                        	}
                                	}
                        	}
                	}
                	if(suc == 1 && bestimproperid >= 0) {
                        	if(output_improper_flag == 1 || allparm_flag == 1) {
                                	fprintf(fpout, "%-2s-%-2s-%-2s-%-2s %11.1lf%15.1lf%12.1lf", name1,
                                                	name2, name3, name4, improperparm[bestimproperid].force,
                                                	improperparm[bestimproperid].phase,
                                                	improperparm[bestimproperid].fterm);
                                	fprintf(fpout, "          Same as %-2s-%-2s-%-2s-%-2s, penalty score=%5.1lf)\n",
                                                	improperparm[bestimproperid].name1,
                                                	improperparm[bestimproperid].name2,
                                                	improperparm[bestimproperid].name3,
                                                	improperparm[bestimproperid].name4, bestscore);
                        	}               
                        	strcpy(improperparm[improperparmnum].name1, name1);
                        	strcpy(improperparm[improperparmnum].name2, name2);    
                        	strcpy(improperparm[improperparmnum].name3, name3);    
                        	strcpy(improperparm[improperparmnum].name4, name4);
                        	improperparm[improperparmnum].phase = improperparm[bestimproperid].phase;
                        	improperparm[improperparmnum].fterm = improperparm[bestimproperid].fterm;
                        	improperparm[improperparmnum].force = improperparm[bestimproperid].force;
                        	improperparmnum++;
                        	if (improperparmnum >= maximproperparm) {
                                	maximproperparm += MAX_FF_IMPROPER;
                                	improperparm = (IMPROPER *) realloc(improperparm, sizeof(IMPROPER) * maximproperparm);
                                	if (improperparm == NULL) {
                                        	fprintf(stdout, "memory allocation error for *improperparm\n");
                                        	exit(1);
                                	} 
                        	}
				bestimproperid = -1;
				continue;
			}
                }    
/*	step 3 considering general atom types */
		if (suc == 0) {
			bestimproperid = -1;
                        bestscore = INITSCORE;
			for (j = 0; j < improperparmnum2; j++) {
				if(improperparm[j].numX == 0) continue;
				if(strcmp(improperparm[j].name1, name1) != 0 && strcmp(improperparm[j].name1, "X") != 0) 
					continue;	
				if(strcmp(improperparm[j].name2, name2) != 0 && strcmp(improperparm[j].name2, "X") != 0) 
					continue;	
				if(strcmp(improperparm[j].name3, name3) != 0 && strcmp(improperparm[j].name3, "X") != 0) 
					continue;	
				if(strcmp(improperparm[j].name4, name4) != 0 && strcmp(improperparm[j].name4, "X") != 0) 
					continue;	
				score1 = 0;
				score2 = 0;
				score3 = 0;
				score4 = 0;
				if(strcmp(improperparm[j].name1, "X") == 0) score1 = wt.X;
				if(strcmp(improperparm[j].name2, "X") == 0) score2 = wt.X;
				if(strcmp(improperparm[j].name3, "X") == 0) score3 = wt.X3;
				if(strcmp(improperparm[j].name4, "X") == 0) score4 = wt.X;
				score = score1 + score2 + score3 + score4;
				if(score < bestscore) {
					bestimproperid = j;
					bestscore = score;
					suc = 1;
				}
			}
                	if(suc == 1 && bestimproperid >= 0) {
                        	if(output_improper_flag == 1 || allparm_flag == 1) {
                               		fprintf(fpout, "%-2s-%-2s-%-2s-%-2s %11.1lf%15.1lf%12.1lf", name1,
                                        	      	name2, name3, name4, improperparm[bestimproperid].force,
                                              		improperparm[bestimproperid].phase,
                                              		improperparm[bestimproperid].fterm);
                               		fprintf(fpout, "          Using general improper torsional angle %2s-%2s-%2s-%2s, penalty score=%5.1lf)\n",
                                       	        	improperparm[bestimproperid].name1,
                                       		       	improperparm[bestimproperid].name2,
                                               		improperparm[bestimproperid].name3,
                                               		improperparm[bestimproperid].name4, bestscore);
                        	}
                       	 	strcpy(improperparm[improperparmnum].name1, name1);
                        	strcpy(improperparm[improperparmnum].name2, name2);
                        	strcpy(improperparm[improperparmnum].name3, name3);
                        	strcpy(improperparm[improperparmnum].name4, name4);
                        	improperparm[improperparmnum].phase = improperparm[bestimproperid].phase;
                        	improperparm[improperparmnum].fterm = improperparm[bestimproperid].fterm;
                        	improperparm[improperparmnum].force = improperparm[bestimproperid].force;
                        	improperparmnum++;
                        	if (improperparmnum >= maximproperparm) {
                                	maximproperparm += MAX_FF_IMPROPER;
                                	improperparm = (IMPROPER *) realloc(improperparm, sizeof(IMPROPER) * maximproperparm);
                                	if (improperparm == NULL) {
                                        	fprintf(stdout, "memory allocation error for *improperparm\n");
                                        	exit(1);
                                	}
                        	}
				bestimproperid = -1;
				continue;
               		}
		}
/*	Step 4 considering both equal atom types and general terms	*/
                if(suc == 0) {          
                        pid1 = parmid[improper[i].atid1];
                        pid2 = parmid[improper[i].atid2];
                        pid3 = parmid[improper[i].atid3];
                        pid4 = parmid[improper[i].atid4];
                        bestscore = INITSCORE;  
                        bestimproperid = -1;            
                        for(m=0; m<parm[pid1].nequa; m++) {
                                pid5 = parm[pid1].equa[m].pid;
                                strcpy(tmpc1, parm[pid5].atomtype);
                                for(n=0;n<parm[pid2].nequa; n++) {
                                        pid6 = parm[pid2].equa[n].pid;
                                        strcpy(tmpc2, parm[pid6].atomtype);
                                        for(p=0;p<parm[pid3].nequa; p++) {
                                                pid7 = parm[pid3].equa[p].pid; 
                                                strcpy(tmpc3, parm[pid7].atomtype);
                                                for(q=0;q<parm[pid4].nequa; q++) {
							if(m==0 && n==0 && p==0 && q==0) continue;  
                                                        pid8 = parm[pid4].equa[q].pid;
                                                        strcpy(tmpc4, parm[pid8].atomtype);
                                                        for (k = 0; k < improperparmnum2; k++){
                                                                if (improperparm[k].numX <= 0) continue;
                                                                if(strcmp(improperparm[k].name1, tmpc1) != 0 && strcmp(improperparm[k].name1, "X") != 0)
                                                                        continue;
                                                                if(strcmp(improperparm[k].name2, tmpc2) != 0 && strcmp(improperparm[k].name2, "X") != 0)
                                                                        continue;
                                                                if(strcmp(improperparm[k].name3, tmpc3) != 0 && strcmp(improperparm[k].name3, "X") != 0)
                                                                        continue;       
                                                                if(strcmp(improperparm[k].name4, tmpc4) != 0 && strcmp(improperparm[k].name4, "X") != 0)
                                                                        continue;  
								score1 = 0;
								score2 = 0;
								score3 = 0;
								score4 = 0;
                                                                if(strcmp(improperparm[k].name1, "X") == 0) score1= wt.X;
                                                                if(strcmp(improperparm[k].name2, "X") == 0) score2= wt.X;
                                                                if(strcmp(improperparm[k].name3, "X") == 0) score3= wt.X3;
                                                                if(strcmp(improperparm[k].name4, "X") == 0) score4= wt.X;
								score = score1 + score2 + score3 + score4;
                                                                if(score < bestscore) {
                                                                        bestimproperid = k;
                                                                        bestscore = score;
                                                                        suc = 1;
                                                                }
                                                                break;
                                                        }
                                                }
                                        }
                                }
                        }
                	if(suc == 1 && bestimproperid >= 0) {
                        	if(output_improper_flag == 1 || allparm_flag == 1) { 
                                	fprintf(fpout, "%-2s-%-2s-%-2s-%-2s %11.1lf%15.1lf%12.1lf", name1,
                                                	name2, name3, name4, improperparm[bestimproperid].force,
                                                	improperparm[bestimproperid].phase,
                                                	improperparm[bestimproperid].fterm);
                                	fprintf(fpout, "          Same as %-2s-%-2s-%-2s-%-2s, penalty score=%5.1lf (use general term))\n",
                                                	improperparm[bestimproperid].name1,
                                                	improperparm[bestimproperid].name2,
                                                	improperparm[bestimproperid].name3,
                                                	improperparm[bestimproperid].name4, bestscore);
                        	}       
                        	strcpy(improperparm[improperparmnum].name1, name1); 
                        	strcpy(improperparm[improperparmnum].name2, name2); 
                        	strcpy(improperparm[improperparmnum].name3, name3); 
                        	strcpy(improperparm[improperparmnum].name4, name4); 
                        	improperparm[improperparmnum].phase = improperparm[bestimproperid].phase;
                        	improperparm[improperparmnum].fterm = improperparm[bestimproperid].fterm;
                        	improperparm[improperparmnum].force = improperparm[bestimproperid].force;
                        	improperparmnum++;
                        	if (improperparmnum >= maximproperparm) {
                                	maximproperparm += MAX_FF_IMPROPER;
                                	improperparm = (IMPROPER *) realloc(improperparm, sizeof(IMPROPER) * maximproperparm);
                                	if (improperparm == NULL) { 
                                        	fprintf(stdout, "memory allocation error for *improperparm\n");
                                        	exit(1);
                                	}       
                        	}       
				bestimproperid = -1;
				continue;
                	}     
		}

/*	Step 5 considering corresponding atom types for specific improper parameters	*/
		if(suc == 0) {
                	pid1 = parmid[improper[i].atid1];
                	pid2 = parmid[improper[i].atid2];
                	pid3 = parmid[improper[i].atid3];
                	pid4 = parmid[improper[i].atid4];
                        bestscore = INITSCORE;
			bestimproperid = -1;
                        for(m=0; m<parm[pid1].ncorr; m++) {
                                pid5 = parm[pid1].corr[m].pid;
                                strcpy(tmpc1, parm[pid5].atomtype);
                                score1= parm[pid1].corr[m].improper;
                                for(n=0;n<parm[pid2].ncorr; n++) {
                                        pid6 = parm[pid2].corr[n].pid;
                                        strcpy(tmpc2, parm[pid6].atomtype);
                                	score2= parm[pid2].corr[n].improper;
                                        for(p=0;p<parm[pid3].ncorr; p++) {
                                                pid7 = parm[pid3].corr[p].pid;
                                                strcpy(tmpc3, parm[pid7].atomtype);
                                		score3= parm[pid3].corr[p].improper * wt.IMPROPER;
                                                for(q=0;q<parm[pid4].ncorr; q++) {
                                                	if(parm[pid1].corr[m].type <= 1 &&
                                                           parm[pid2].corr[n].type <= 1 &&
                                                           parm[pid3].corr[p].type <= 1 &&
                                                           parm[pid4].corr[q].type <= 1)
                                                        	continue;
                                                        pid8 = parm[pid4].corr[q].pid;
                                                        strcpy(tmpc4, parm[pid8].atomtype);
                                			score4= parm[pid4].corr[q].improper;
                                                        score = score1 + score2 + score3 + score4;
							score += wt.GROUP;
							if(parm[pid5].group == parm[pid6].group && 
                                   			   parm[pid5].group == parm[pid7].group &&
                                   			   parm[pid5].group == parm[pid8].group)
			   	       			   score -= wt.GROUP;
                                                        if(score < bestscore)  
								for (k = 0; k < improperparmnum2; k++) {
									if (improperparm[k].numX > 0) continue; 
									if(strcmp(improperparm[k].name1, tmpc1) == 0 && 
			   						   strcmp(improperparm[k].name2, tmpc2) == 0 && 
			   						   strcmp(improperparm[k].name3, tmpc3) == 0 && 
			   						   strcmp(improperparm[k].name4, tmpc4) == 0) {
										bestimproperid = k;
                                                                		bestscore = score;
										suc = 1;
										break;
									}
								}
                                                }
                                        }
                                }
                        }
                	if(suc == 1 && bestimproperid >= 0) {
                        	if(output_improper_flag == 1 || allparm_flag == 1) {
                               		fprintf(fpout, "%-2s-%-2s-%-2s-%-2s %11.1lf%15.1lf%12.1lf", name1,
                                       		       	name2, name3, name4, improperparm[bestimproperid].force,
                                              		improperparm[bestimproperid].phase,
                                              		improperparm[bestimproperid].fterm);
                               		fprintf(fpout, "          Same as %-2s-%-2s-%-2s-%-2s, penalty score=%5.1lf)\n",
                                               		improperparm[bestimproperid].name1,
                                               		improperparm[bestimproperid].name2,
                                               		improperparm[bestimproperid].name3,
                                               		improperparm[bestimproperid].name4, bestscore);
                        	}
                        	strcpy(improperparm[improperparmnum].name1, name1);
                        	strcpy(improperparm[improperparmnum].name2, name2);
                        	strcpy(improperparm[improperparmnum].name3, name3);
                        	strcpy(improperparm[improperparmnum].name4, name4);
                        	improperparm[improperparmnum].phase = improperparm[bestimproperid].phase;
                        	improperparm[improperparmnum].fterm = improperparm[bestimproperid].fterm;
                        	improperparm[improperparmnum].force = improperparm[bestimproperid].force;
                        	improperparmnum++;
                        	if (improperparmnum >= maximproperparm) {
                               		maximproperparm += MAX_FF_IMPROPER;
                                	improperparm = (IMPROPER *) realloc(improperparm, sizeof(IMPROPER) * maximproperparm);
                                	if (improperparm == NULL) {
                                        	fprintf(stdout, "memory allocation error for *improperparm\n");
                                        	exit(1);
                                	}
                        	}
				bestimproperid = -1;
				continue;
                	}
		}
/*	Step6  considering corresponding atom types for general improper parameters	*/
		if(suc == 0) {
                	pid1 = parmid[improper[i].atid1];
                	pid2 = parmid[improper[i].atid2];
                	pid3 = parmid[improper[i].atid3];
                	pid4 = parmid[improper[i].atid4];
                        bestscore = INITSCORE;
			bestimproperid = -1;
                        for(m=0; m<parm[pid1].ncorr; m++) {
                                pid5 = parm[pid1].corr[m].pid;
                                strcpy(tmpc1, parm[pid5].atomtype);
                                score1= parm[pid1].corr[m].improper ; 
                                for(n=0;n<parm[pid2].ncorr; n++) {
                                        pid6 = parm[pid2].corr[n].pid;
                                        strcpy(tmpc2, parm[pid6].atomtype);
                                        score2= parm[pid2].corr[n].improper;
                                        for(p=0;p<parm[pid3].ncorr; p++) {
                                                pid7 = parm[pid3].corr[p].pid;
                                                strcpy(tmpc3, parm[pid7].atomtype);
                                                score3= parm[pid3].corr[p].improper + wt.IMPROPER;
                                                for(q=0;q<parm[pid4].ncorr; q++) {
                                                	if(parm[pid1].corr[m].type <= 1 &&
                                                           parm[pid2].corr[n].type <= 1 &&
                                                           parm[pid3].corr[p].type <= 1 &&
                                                           parm[pid4].corr[q].type <= 1)
                                                        	continue;
                                                        pid8 = parm[pid4].corr[q].pid;
                                                        strcpy(tmpc4, parm[pid8].atomtype);
                                                        score4= parm[pid4].corr[q].improper; 
                                                        score = score1 + score2 + score3 + score4;
							score += wt.GROUP;
							if(parm[pid5].group == parm[pid6].group && 
                                   			   parm[pid5].group == parm[pid7].group &&
                                   			   parm[pid5].group == parm[pid8].group)
			   	       			   score -= wt.GROUP;
							for (k = 0; k < improperparmnum2; k++){
								if (improperparm[k].numX <= 0) continue; 
								if(strcmp(improperparm[k].name1, tmpc1) != 0 && strcmp(improperparm[k].name1, "X") != 0) 
									continue;	
								if(strcmp(improperparm[k].name2, tmpc2) != 0 && strcmp(improperparm[k].name2, "X") != 0) 
									continue;	
								if(strcmp(improperparm[k].name3, tmpc3) != 0 && strcmp(improperparm[k].name3, "X") != 0) 
									continue;	
								if(strcmp(improperparm[k].name4, tmpc4) != 0 && strcmp(improperparm[k].name4, "X") != 0) 
									continue;	
								if(strcmp(improperparm[k].name1, "X") == 0) score += wt.X;
								if(strcmp(improperparm[k].name2, "X") == 0) score += wt.X;
								if(strcmp(improperparm[k].name3, "X") == 0) score += wt.X3;
								if(strcmp(improperparm[k].name4, "X") == 0) score += wt.X;
								if(score < bestscore) {
									bestimproperid = k;
                                                               		bestscore = score;
									suc = 1;
								}
								break;
							}
                                                }
                                        }
                                }
                        }
                	if(suc == 1 && bestimproperid >= 0) {
                        	if(output_improper_flag == 1 || allparm_flag == 1) {
                               		fprintf(fpout, "%-2s-%-2s-%-2s-%-2s %11.1lf%15.1lf%12.1lf", name1,
                                       	        name2, name3, name4, improperparm[bestimproperid].force,
                                       	        improperparm[bestimproperid].phase,
                                      	        improperparm[bestimproperid].fterm);
                               		fprintf(fpout, "          Same as %-2s-%-2s-%-2s-%-2s, penalty score=%5.1lf (use general term))\n",
                                       	        improperparm[bestimproperid].name1,
                                       	        improperparm[bestimproperid].name2,
                                       	        improperparm[bestimproperid].name3,
                                       	        improperparm[bestimproperid].name4, bestscore);
                        	}
                        	strcpy(improperparm[improperparmnum].name1, name1);
                        	strcpy(improperparm[improperparmnum].name2, name2);
                        	strcpy(improperparm[improperparmnum].name3, name3);
                        	strcpy(improperparm[improperparmnum].name4, name4);
                        	improperparm[improperparmnum].phase = improperparm[bestimproperid].phase;
                        	improperparm[improperparmnum].fterm = improperparm[bestimproperid].fterm;
                        	improperparm[improperparmnum].force = improperparm[bestimproperid].force;
                        	improperparmnum++;
                        	if (improperparmnum >= maximproperparm) {
                                	maximproperparm += MAX_FF_IMPROPER;
                                	improperparm = (IMPROPER *) realloc(improperparm, sizeof(IMPROPER) * maximproperparm);
                                	if (improperparm == NULL) {
                                        	fprintf(stdout, "memory allocation error for *improperparm\n");
                                        	exit(1);
                                	}
                        	}
				bestimproperid = -1;
				continue;
                	}
		}
		if (suc == 0) {
			fprintf(fpout, "%-2s-%-2s-%-2s-%-2s %11.1lf%15.1lf%12.1lf", name1,
					name2, name3, name4, 1.1, 180.0, 2.0);
			fprintf(fpout, "          Using the default value\n");
			strcpy(improperparm[improperparmnum].name1, name1);
			strcpy(improperparm[improperparmnum].name2, name2);
			strcpy(improperparm[improperparmnum].name3, name3);
			strcpy(improperparm[improperparmnum].name4, name4);
			improperparm[improperparmnum].phase = 180.0;
			improperparm[improperparmnum].fterm = 2.0;
			improperparm[improperparmnum].force = 1.1;
			improperparmnum++;
			if (improperparmnum >= maximproperparm) {
				maximproperparm += MAX_FF_IMPROPER;
				improperparm =
					(IMPROPER *) realloc(improperparm,
										 sizeof(IMPROPER) *
										 maximproperparm);
				if (improperparm == NULL) {
					fprintf(stdout,
							"memory allocation error for *improperparm\n");
					exit(1);
				}
			}

		}
	}
}
void chk_improper_v2 (void){
	int i,j,k;
	int m,n,p,q;
	int suc;
	int tmpnum;
	int index1, index2, index4;
        char tmpc[5];
        char tmpc1[5], tmpc2[5], tmpc3[5], tmpc4[5];
        char name1[5], name2[5], name3[5], name4[5];
        int pid1, pid2, pid3, pid4;
        int pid5, pid6, pid7, pid8;
	double score, score1, score2, score3, score4;
	fprintf(fpout, "\n%s\n", "IMPROPER");

	for (i = 0; i < r_improperparm_num; i++) {
		if(iformat == 4 && attn_opt == 2 && r_improperparm[i].attn == 0) {
                    	fprintf(fpout, "%-2s-%-2s-%-2s-%-2s %11.1lf%15.1lf%12.1lf\n", 
				r_improperparm[i].name1, r_improperparm[i].name2, 
                                r_improperparm[i].name3, r_improperparm[i].name4, 
				r_improperparm[i].force, r_improperparm[i].phase,
                                r_improperparm[i].fterm);
			continue;
		}
		suc = 0;
		strcpy(tmpc1, r_improperparm[i].name1); 
		strcpy(tmpc2, r_improperparm[i].name2); 
		strcpy(tmpc3, r_improperparm[i].name3); 
		strcpy(tmpc4, r_improperparm[i].name4); 

		strcpy(name1, tmpc1);
		strcpy(name2, tmpc2);
		strcpy(name3, tmpc3);
		strcpy(name4, tmpc4);
/*	step 1 check directly	*/
		for (j = 0; j < improperparmnum; j++)
			if(strcmp(improperparm[j].name1, name1) == 0 && 
			   strcmp(improperparm[j].name2, name2) == 0 && 
			   strcmp(improperparm[j].name3, name3) == 0 && 
			   strcmp(improperparm[j].name4, name4) == 0) {
				suc = 1;
				if(allparm_flag == 1) {
                                	fprintf(fpout, "%-2s-%-2s-%-2s-%-2s %11.1lf%15.1lf%12.1lf\n", name1,
                                        	name2, name3, name4, improperparm[j].force, improperparm[j].phase,
                                                improperparm[j].fterm);
				}
				break;
			}
		if(suc == 1) continue;
/* 	Step 2 check special improper torsional parameters using equal atom types */
                if(suc == 0) {                                         
			pid1 = assign_parmid_v2(name1);
			pid2 = assign_parmid_v2(name2);
			pid3 = assign_parmid_v2(name3);
			pid4 = assign_parmid_v2(name4);
                	bestimproperid = -1;
                        for(m=0; m<parm[pid1].nequa; m++) {
				if(suc == 1) break;
                                pid5 = parm[pid1].equa[m].pid;
                                strcpy(tmpc1, parm[pid5].atomtype);
                                for(n=0;n<parm[pid2].nequa; n++) {
					if(suc == 1) break;
                                        pid6 = parm[pid2].equa[n].pid;
                                        strcpy(tmpc2, parm[pid6].atomtype);
                                        for(p=0;p<parm[pid3].nequa; p++) {
						if(suc == 1) break;
                                                pid7 = parm[pid3].equa[p].pid;
                                                strcpy(tmpc3, parm[pid7].atomtype);
                                                for(q=0;q<parm[pid4].nequa; q++) {
							if(suc == 1) break;
                                                        if(m==0 && n==0 && p== 0 && q==0) continue;
                                                        pid8 = parm[pid4].equa[q].pid;
                                                        strcpy(tmpc4, parm[pid8].atomtype);
                                                        for (k = 0; k < improperparmnum2; k++) {
                                                                if (improperparm[k].numX > 0) continue;
                                                                if(strcmp(improperparm[k].name1, tmpc1) == 0 &&
                                                                   strcmp(improperparm[k].name2, tmpc2) == 0 &&
                                                                   strcmp(improperparm[k].name3, tmpc3) == 0 &&
                                                                   strcmp(improperparm[k].name4, tmpc4) == 0) {
                                                                        bestimproperid = k;
                                                                        bestscore = 0;
                                                                        suc = 1;
                                                                        break;
                                                        	}
                                                	}
                                        	}
                                	}
                        	}
                	}
                	if(suc == 1 && bestimproperid >= 0) {
                        	if(output_improper_flag == 1 || allparm_flag == 1) {
                                	fprintf(fpout, "%-2s-%-2s-%-2s-%-2s %11.1lf%15.1lf%12.1lf", name1,
                                                	name2, name3, name4, improperparm[bestimproperid].force,
                                                	improperparm[bestimproperid].phase,
                                                	improperparm[bestimproperid].fterm);
                                	fprintf(fpout, "          Same as %-2s-%-2s-%-2s-%-2s, penalty score=%5.1lf)\n",
                                                	improperparm[bestimproperid].name1,
                                                	improperparm[bestimproperid].name2,
                                                	improperparm[bestimproperid].name3,
                                                	improperparm[bestimproperid].name4, bestscore);
                        	}               
                        	strcpy(improperparm[improperparmnum].name1, name1);
                        	strcpy(improperparm[improperparmnum].name2, name2);    
                        	strcpy(improperparm[improperparmnum].name3, name3);    
                        	strcpy(improperparm[improperparmnum].name4, name4);
                        	improperparm[improperparmnum].phase = improperparm[bestimproperid].phase;
                        	improperparm[improperparmnum].fterm = improperparm[bestimproperid].fterm;
                        	improperparm[improperparmnum].force = improperparm[bestimproperid].force;
                        	improperparmnum++;
                        	if (improperparmnum >= maximproperparm) {
                                	maximproperparm += MAX_FF_IMPROPER;
                                	improperparm = (IMPROPER *) realloc(improperparm, sizeof(IMPROPER) * maximproperparm);
                                	if (improperparm == NULL) {
                                        	fprintf(stdout, "memory allocation error for *improperparm\n");
                                        	exit(1);
                                	} 
                        	}
				bestimproperid = -1;
				continue;
			}
                }    
/*	step 3 considering general atom types */
		if (suc == 0) {
			bestimproperid = -1;
                        bestscore = INITSCORE;
			for (j = 0; j < improperparmnum2; j++) {
				if(improperparm[j].numX == 0) continue;
				if(strcmp(improperparm[j].name1, name1) != 0 && strcmp(improperparm[j].name1, "X") != 0) 
					continue;	
				if(strcmp(improperparm[j].name2, name2) != 0 && strcmp(improperparm[j].name2, "X") != 0) 
					continue;	
				if(strcmp(improperparm[j].name3, name3) != 0 && strcmp(improperparm[j].name3, "X") != 0) 
					continue;	
				if(strcmp(improperparm[j].name4, name4) != 0 && strcmp(improperparm[j].name4, "X") != 0) 
					continue;	
				score1 = 0;
				score2 = 0;
				score3 = 0;
				score4 = 0;
				if(strcmp(improperparm[j].name1, "X") == 0) score1 = wt.X;
				if(strcmp(improperparm[j].name2, "X") == 0) score2 = wt.X;
				if(strcmp(improperparm[j].name3, "X") == 0) score3 = wt.X3;
				if(strcmp(improperparm[j].name4, "X") == 0) score4 = wt.X;
				score = score1 + score2 + score3 + score4;
				if(score < bestscore) {
					bestimproperid = j;
					bestscore = score;
					suc = 1;
				}
			}
                	if(suc == 1 && bestimproperid >= 0) {
                        	if(output_improper_flag == 1 || allparm_flag == 1) {
                               		fprintf(fpout, "%-2s-%-2s-%-2s-%-2s %11.1lf%15.1lf%12.1lf", name1,
                                        	      	name2, name3, name4, improperparm[bestimproperid].force,
                                              		improperparm[bestimproperid].phase,
                                              		improperparm[bestimproperid].fterm);
                               		fprintf(fpout, "          Using general improper torsional angle %2s-%2s-%2s-%2s, penalty score=%5.1lf)\n",
                                       	        	improperparm[bestimproperid].name1,
                                       		       	improperparm[bestimproperid].name2,
                                               		improperparm[bestimproperid].name3,
                                               		improperparm[bestimproperid].name4, bestscore);
                        	}
                       	 	strcpy(improperparm[improperparmnum].name1, name1);
                        	strcpy(improperparm[improperparmnum].name2, name2);
                        	strcpy(improperparm[improperparmnum].name3, name3);
                        	strcpy(improperparm[improperparmnum].name4, name4);
                        	improperparm[improperparmnum].phase = improperparm[bestimproperid].phase;
                        	improperparm[improperparmnum].fterm = improperparm[bestimproperid].fterm;
                        	improperparm[improperparmnum].force = improperparm[bestimproperid].force;
                        	improperparmnum++;
                        	if (improperparmnum >= maximproperparm) {
                                	maximproperparm += MAX_FF_IMPROPER;
                                	improperparm = (IMPROPER *) realloc(improperparm, sizeof(IMPROPER) * maximproperparm);
                                	if (improperparm == NULL) {
                                        	fprintf(stdout, "memory allocation error for *improperparm\n");
                                        	exit(1);
                                	}
                        	}
				bestimproperid = -1;
				continue;
               		}
		}
/*	Step 4 considering both equal atom types and general terms	*/
                if(suc == 0) {          
			pid1 = assign_parmid_v2(name1);
			pid2 = assign_parmid_v2(name2);
			pid3 = assign_parmid_v2(name3);
			pid4 = assign_parmid_v2(name4);
                        bestscore = INITSCORE;  
                        bestimproperid = -1;            
                        for(m=0; m<parm[pid1].nequa; m++) {
                                pid5 = parm[pid1].equa[m].pid;
                                strcpy(tmpc1, parm[pid5].atomtype);
                                for(n=0;n<parm[pid2].nequa; n++) {
                                        pid6 = parm[pid2].equa[n].pid;
                                        strcpy(tmpc2, parm[pid6].atomtype);
                                        for(p=0;p<parm[pid3].nequa; p++) {
                                                pid7 = parm[pid3].equa[p].pid; 
                                                strcpy(tmpc3, parm[pid7].atomtype);
                                                for(q=0;q<parm[pid4].nequa; q++) {
							if(m==0 && n==0 && p==0 && q==0) continue;  
                                                        pid8 = parm[pid4].equa[q].pid;
                                                        strcpy(tmpc4, parm[pid8].atomtype);
                                                        for (k = 0; k < improperparmnum2; k++){
                                                                if (improperparm[k].numX <= 0) continue;
                                                                if(strcmp(improperparm[k].name1, tmpc1) != 0 && strcmp(improperparm[k].name1, "X") != 0)
                                                                        continue;
                                                                if(strcmp(improperparm[k].name2, tmpc2) != 0 && strcmp(improperparm[k].name2, "X") != 0)
                                                                        continue;
                                                                if(strcmp(improperparm[k].name3, tmpc3) != 0 && strcmp(improperparm[k].name3, "X") != 0)
                                                                        continue;       
                                                                if(strcmp(improperparm[k].name4, tmpc4) != 0 && strcmp(improperparm[k].name4, "X") != 0)
                                                                        continue;  
								score1 = 0;
								score2 = 0;
								score3 = 0;
								score4 = 0;
                                                                if(strcmp(improperparm[k].name1, "X") == 0) score1= wt.X;
                                                                if(strcmp(improperparm[k].name2, "X") == 0) score2= wt.X;
                                                                if(strcmp(improperparm[k].name3, "X") == 0) score3= wt.X3;
                                                                if(strcmp(improperparm[k].name4, "X") == 0) score4= wt.X;
								score = score1 + score2 + score3 + score4;
                                                                if(score < bestscore) {
                                                                        bestimproperid = k;
                                                                        bestscore = score;
                                                                        suc = 1;
                                                                }
                                                                break;
                                                        }
                                                }
                                        }
                                }
                        }
                	if(suc == 1 && bestimproperid >= 0) {
                        	if(output_improper_flag == 1 || allparm_flag == 1) { 
                                	fprintf(fpout, "%-2s-%-2s-%-2s-%-2s %11.1lf%15.1lf%12.1lf", name1,
                                                	name2, name3, name4, improperparm[bestimproperid].force,
                                                	improperparm[bestimproperid].phase,
                                                	improperparm[bestimproperid].fterm);
                                	fprintf(fpout, "          Same as %-2s-%-2s-%-2s-%-2s, penalty score=%5.1lf (use general term))\n",
                                                	improperparm[bestimproperid].name1,
                                                	improperparm[bestimproperid].name2,
                                                	improperparm[bestimproperid].name3,
                                                	improperparm[bestimproperid].name4, bestscore);
                        	}       
                        	strcpy(improperparm[improperparmnum].name1, name1); 
                        	strcpy(improperparm[improperparmnum].name2, name2); 
                        	strcpy(improperparm[improperparmnum].name3, name3); 
                        	strcpy(improperparm[improperparmnum].name4, name4); 
                        	improperparm[improperparmnum].phase = improperparm[bestimproperid].phase;
                        	improperparm[improperparmnum].fterm = improperparm[bestimproperid].fterm;
                        	improperparm[improperparmnum].force = improperparm[bestimproperid].force;
                        	improperparmnum++;
                        	if (improperparmnum >= maximproperparm) {
                                	maximproperparm += MAX_FF_IMPROPER;
                                	improperparm = (IMPROPER *) realloc(improperparm, sizeof(IMPROPER) * maximproperparm);
                                	if (improperparm == NULL) { 
                                        	fprintf(stdout, "memory allocation error for *improperparm\n");
                                        	exit(1);
                                	}       
                        	}       
				bestimproperid = -1;
				continue;
                	}     
		}

/*	Step 5 considering corresponding atom types for specific improper parameters	*/
		if(suc == 0) {
			pid1 = assign_parmid_v2(name1);
			pid2 = assign_parmid_v2(name2);
			pid3 = assign_parmid_v2(name3);
			pid4 = assign_parmid_v2(name4);

                        bestscore = INITSCORE;
			bestimproperid = -1;
                        for(m=0; m<parm[pid1].ncorr; m++) {
                                pid5 = parm[pid1].corr[m].pid;
                                strcpy(tmpc1, parm[pid5].atomtype);
                                score1= parm[pid1].corr[m].improper;
                                for(n=0;n<parm[pid2].ncorr; n++) {
                                        pid6 = parm[pid2].corr[n].pid;
                                        strcpy(tmpc2, parm[pid6].atomtype);
                                	score2= parm[pid2].corr[n].improper;
                                        for(p=0;p<parm[pid3].ncorr; p++) {
                                                pid7 = parm[pid3].corr[p].pid;
                                                strcpy(tmpc3, parm[pid7].atomtype);
                                		score3= parm[pid3].corr[p].improper * wt.IMPROPER;
                                                for(q=0;q<parm[pid4].ncorr; q++) {
                                                	if(parm[pid1].corr[m].type <= 1 &&
                                                           parm[pid2].corr[n].type <= 1 &&
                                                           parm[pid3].corr[p].type <= 1 &&
                                                           parm[pid4].corr[q].type <= 1)
                                                        	continue;
                                                        pid8 = parm[pid4].corr[q].pid;
                                                        strcpy(tmpc4, parm[pid8].atomtype);
                                			score4= parm[pid4].corr[q].improper;
                                                        score = score1 + score2 + score3 + score4;
							score += wt.GROUP;
							if(parm[pid5].group == parm[pid6].group && 
                                   			   parm[pid5].group == parm[pid7].group &&
                                   			   parm[pid5].group == parm[pid8].group)
			   	       			   score -= wt.GROUP;
                                                        if(score < bestscore)  
								for (k = 0; k < improperparmnum2; k++) {
									if (improperparm[k].numX > 0) continue; 
									if(strcmp(improperparm[k].name1, tmpc1) == 0 && 
			   						   strcmp(improperparm[k].name2, tmpc2) == 0 && 
			   						   strcmp(improperparm[k].name3, tmpc3) == 0 && 
			   						   strcmp(improperparm[k].name4, tmpc4) == 0) {
										bestimproperid = k;
                                                                		bestscore = score;
										suc = 1;
										break;
									}
								}
                                                }
                                        }
                                }
                        }
                	if(suc == 1 && bestimproperid >= 0) {
                        	if(output_improper_flag == 1 || allparm_flag == 1) {
                               		fprintf(fpout, "%-2s-%-2s-%-2s-%-2s %11.1lf%15.1lf%12.1lf", name1,
                                       		       	name2, name3, name4, improperparm[bestimproperid].force,
                                              		improperparm[bestimproperid].phase,
                                              		improperparm[bestimproperid].fterm);
                               		fprintf(fpout, "          Same as %-2s-%-2s-%-2s-%-2s, penalty score=%5.1lf)\n",
                                               		improperparm[bestimproperid].name1,
                                               		improperparm[bestimproperid].name2,
                                               		improperparm[bestimproperid].name3,
                                               		improperparm[bestimproperid].name4, bestscore);
                        	}
                        	strcpy(improperparm[improperparmnum].name1, name1);
                        	strcpy(improperparm[improperparmnum].name2, name2);
                        	strcpy(improperparm[improperparmnum].name3, name3);
                        	strcpy(improperparm[improperparmnum].name4, name4);
                        	improperparm[improperparmnum].phase = improperparm[bestimproperid].phase;
                        	improperparm[improperparmnum].fterm = improperparm[bestimproperid].fterm;
                        	improperparm[improperparmnum].force = improperparm[bestimproperid].force;
                        	improperparmnum++;
                        	if (improperparmnum >= maximproperparm) {
                               		maximproperparm += MAX_FF_IMPROPER;
                                	improperparm = (IMPROPER *) realloc(improperparm, sizeof(IMPROPER) * maximproperparm);
                                	if (improperparm == NULL) {
                                        	fprintf(stdout, "memory allocation error for *improperparm\n");
                                        	exit(1);
                                	}
                        	}
				bestimproperid = -1;
				continue;
                	}
		}
/*	Step6  considering corresponding atom types for general improper parameters	*/
		if(suc == 0) {
			pid1 = assign_parmid_v2(name1);
			pid2 = assign_parmid_v2(name2);
			pid3 = assign_parmid_v2(name3);
			pid4 = assign_parmid_v2(name4);
                        bestscore = INITSCORE;
			bestimproperid = -1;
                        for(m=0; m<parm[pid1].ncorr; m++) {
                                pid5 = parm[pid1].corr[m].pid;
                                strcpy(tmpc1, parm[pid5].atomtype);
                                score1= parm[pid1].corr[m].improper ; 
                                for(n=0;n<parm[pid2].ncorr; n++) {
                                        pid6 = parm[pid2].corr[n].pid;
                                        strcpy(tmpc2, parm[pid6].atomtype);
                                        score2= parm[pid2].corr[n].improper;
                                        for(p=0;p<parm[pid3].ncorr; p++) {
                                                pid7 = parm[pid3].corr[p].pid;
                                                strcpy(tmpc3, parm[pid7].atomtype);
                                                score3= parm[pid3].corr[p].improper + wt.IMPROPER;
                                                for(q=0;q<parm[pid4].ncorr; q++) {
                                                	if(parm[pid1].corr[m].type <= 1 &&
                                                           parm[pid2].corr[n].type <= 1 &&
                                                           parm[pid3].corr[p].type <= 1 &&
                                                           parm[pid4].corr[q].type <= 1)
                                                        	continue;
                                                        pid8 = parm[pid4].corr[q].pid;
                                                        strcpy(tmpc4, parm[pid8].atomtype);
                                                        score4= parm[pid4].corr[q].improper; 
                                                        score = score1 + score2 + score3 + score4;
							score += wt.GROUP;
							if(parm[pid5].group == parm[pid6].group && 
                                   			   parm[pid5].group == parm[pid7].group &&
                                   			   parm[pid5].group == parm[pid8].group)
			   	       			   score -= wt.GROUP;
							for (k = 0; k < improperparmnum2; k++){
								if (improperparm[k].numX <= 0) continue; 
								if(strcmp(improperparm[k].name1, tmpc1) != 0 && strcmp(improperparm[k].name1, "X") != 0) 
									continue;	
								if(strcmp(improperparm[k].name2, tmpc2) != 0 && strcmp(improperparm[k].name2, "X") != 0) 
									continue;	
								if(strcmp(improperparm[k].name3, tmpc3) != 0 && strcmp(improperparm[k].name3, "X") != 0) 
									continue;	
								if(strcmp(improperparm[k].name4, tmpc4) != 0 && strcmp(improperparm[k].name4, "X") != 0) 
									continue;	
								if(strcmp(improperparm[k].name1, "X") == 0) score += wt.X;
								if(strcmp(improperparm[k].name2, "X") == 0) score += wt.X;
								if(strcmp(improperparm[k].name3, "X") == 0) score += wt.X3;
								if(strcmp(improperparm[k].name4, "X") == 0) score += wt.X;
								if(score < bestscore) {
									bestimproperid = k;
                                                               		bestscore = score;
									suc = 1;
								}
								break;
							}
                                                }
                                        }
                                }
                        }
                	if(suc == 1 && bestimproperid >= 0) {
                        	if(output_improper_flag == 1 || allparm_flag == 1) {
                               		fprintf(fpout, "%-2s-%-2s-%-2s-%-2s %11.1lf%15.1lf%12.1lf", name1,
                                       	        name2, name3, name4, improperparm[bestimproperid].force,
                                       	        improperparm[bestimproperid].phase,
                                      	        improperparm[bestimproperid].fterm);
                               		fprintf(fpout, "          Same as %-2s-%-2s-%-2s-%-2s, penalty score=%5.1lf (use general term))\n",
                                       	        improperparm[bestimproperid].name1,
                                       	        improperparm[bestimproperid].name2,
                                       	        improperparm[bestimproperid].name3,
                                       	        improperparm[bestimproperid].name4, bestscore);
                        	}
                        	strcpy(improperparm[improperparmnum].name1, name1);
                        	strcpy(improperparm[improperparmnum].name2, name2);
                        	strcpy(improperparm[improperparmnum].name3, name3);
                        	strcpy(improperparm[improperparmnum].name4, name4);
                        	improperparm[improperparmnum].phase = improperparm[bestimproperid].phase;
                        	improperparm[improperparmnum].fterm = improperparm[bestimproperid].fterm;
                        	improperparm[improperparmnum].force = improperparm[bestimproperid].force;
                        	improperparmnum++;
                        	if (improperparmnum >= maximproperparm) {
                                	maximproperparm += MAX_FF_IMPROPER;
                                	improperparm = (IMPROPER *) realloc(improperparm, sizeof(IMPROPER) * maximproperparm);
                                	if (improperparm == NULL) {
                                        	fprintf(stdout, "memory allocation error for *improperparm\n");
                                        	exit(1);
                                	}
                        	}
				bestimproperid = -1;
				continue;
                	}
		}
		if (suc == 0) {
			fprintf(fpout, "%-2s-%-2s-%-2s-%-2s %11.1lf%15.1lf%12.1lf", name1,
					name2, name3, name4, 1.1, 180.0, 2.0);
			fprintf(fpout, "          Using the default value\n");
			strcpy(improperparm[improperparmnum].name1, name1);
			strcpy(improperparm[improperparmnum].name2, name2);
			strcpy(improperparm[improperparmnum].name3, name3);
			strcpy(improperparm[improperparmnum].name4, name4);
			improperparm[improperparmnum].phase = 180.0;
			improperparm[improperparmnum].fterm = 2.0;
			improperparm[improperparmnum].force = 1.1;
			improperparmnum++;
			if (improperparmnum >= maximproperparm) {
				maximproperparm += MAX_FF_IMPROPER;
				improperparm =
					(IMPROPER *) realloc(improperparm,
										 sizeof(IMPROPER) *
										 maximproperparm);
				if (improperparm == NULL) {
					fprintf(stdout,
							"memory allocation error for *improperparm\n");
					exit(1);
				}
			}

		}
	}

}

void cleanup_frcmod(char * filename) {
	FILE *fp1, *fp2;
	char command[MAXCHAR] = "cp -f ";
	int status;
	typedef struct {
        	char name[10];
	} ATOMTYPE_CL;
	typedef struct {
        	char name1[10];
        	char name2[10];
	} BOND_CL;
	typedef struct {
        	char name1[10];
        	char name2[10];
        	char name3[10];
	} ANGLE_CL;
	typedef struct {
        	char name1[10];
        	char name2[10];
        	char name3[10];
        	char name4[10];
		double fterm;
	} TORSION_CL;
	ATOMTYPE_CL at[MAX_FF_ATOMTYPE];
	int atnum = 0;
	BOND_CL bond[MAX_FF_BOND];
	int bondnum = 0;
	ANGLE_CL angle[MAX_FF_ANGLE];
	int anglenum = 0;
	TORSION_CL tor[MAX_FF_TORSION];
	int tornum = 0;
	TORSION_CL improp[MAX_FF_IMPROPER];
	int impropnum = 0;
	ATOMTYPE_CL vdw[MAX_FF_VDW];	
	int vdwnum = 0;
	int writeflag = 1;
	int typeflag = 0;
	int tmpint;
	char name1[10];
	char name2[10];
	char name3[10];
	char name4[10];
	char line[MAXCHAR];
	char line2[MAXCHAR];
	char tmpchar[MAXCHAR];
	double tmpf1, tmpf2, tmpf3;

	strcat(command, filename);
	strcat(command, " ANTECHAMBER.FRCMOD");
	status = system(command);
        if(status != 0) {
                fprintf(stdout, "Error: cannot run \"%s\" in cleanup_frcmod function properly, exit\n", command);
                exit(1);
        }
	
        if ((fp1 = fopen("ANTECHAMBER.FRCMOD", "r")) == NULL) {
		fprintf(stdout, "Error: Cannot open the redundant frcmod file - ANTECHAMBER.FRCMOD to do clean up\n"); 
		exit(1);
        }
        if ((fp2 = fopen(filename, "w")) == NULL) {
		fprintf(stdout, "Error: Cannot open the frcmod file %s to write out\n", filename); 
		exit(1);
        }
	for(;;) {
                if (fgets(line, MAXCHAR, fp1) == NULL) break;
                sscanf(line, "%s", tmpchar);
                if (strcmp("MASS", tmpchar) == 0) 
			typeflag = 1;
                if (strcmp("BOND", tmpchar) == 0) 
			typeflag = 2;
                if (strcmp("ANGLE", tmpchar) == 0) 
			typeflag = 3;
                if (strcmp("DIHE", tmpchar) == 0) 
			typeflag = 4;
                if (strcmp("IMPROPER", tmpchar) == 0) 
			typeflag = 5;
                if (strcmp("NONBON", tmpchar) == 0) 
			typeflag = 6;
		if(strlen(line) <= 1) {
			typeflag = 0;
			writeflag = 1;
		}
		if(typeflag == 1) {
			writeflag = 1;
			for(i=0;i<atnum;i++)  {
				if(strcmp(tmpchar, at[i].name) == 0) {
					writeflag = 0;
					break;
				}
			}
			if(writeflag == 1) {
				strcpy(at[atnum].name, tmpchar);
				atnum++;
			}
		}

		if(typeflag == 2) {
			writeflag = 1;
			strcpy(line2, line);
			for(i=0;i<6;i++)  /* we do not want to delete other '-' after the 6th columns */
				if(line2[i] == '-') line2[i] = ' ';
			sscanf(line2, "%s%s", name1, name2);
			for(i=0;i<bondnum;i++)  {
				if(strcmp(name1, bond[i].name1) == 0 && strcmp(name2, bond[i].name2) == 0) {
					writeflag = 0;
					break;
				}
				if(strcmp(name1, bond[i].name2) == 0 && strcmp(name2, bond[i].name1) == 0) {
					writeflag = 0;
					break;
				}
			}
			if(writeflag == 1) {
				strcpy(bond[bondnum].name1, name1);
				strcpy(bond[bondnum].name2, name2);
				bondnum++;
			}
		}

                if(typeflag == 3) {
                        writeflag = 1;
			strcpy(line2, line);
			for(i=0;i<9;i++)  /* we do not want to delete other '-' after the 9th columns */
				if(line2[i] == '-') line2[i] = ' ';
                        sscanf(line2, "%s%s%s", name1, name2,name3);
                        for(i=0;i<anglenum;i++)  {
                                if(strcmp(name1, angle[i].name1) == 0 && strcmp(name2, angle[i].name2) == 0 && strcmp(name3, angle[i].name3) == 0 ) {
                                        writeflag = 0;
                                        break;
                                }
                                if(strcmp(name1, angle[i].name3) == 0 && strcmp(name2, angle[i].name2) == 0 && strcmp(name3, angle[i].name1) == 0 ) {
                                        writeflag = 0;
                                        break;
                                }
                        }
                        if(writeflag == 1) {
                                strcpy(angle[anglenum].name1, name1);
                                strcpy(angle[anglenum].name2, name2);
                                strcpy(angle[anglenum].name3, name3);
                                anglenum++;
                        }
                }

                if(typeflag == 4) {
                        writeflag = 1;
			strcpy(line2, line);
			for(i=0;i<12;i++)  /* we do not want to delete '-' after 12th column, esp., '-' before fterm */
				if(line2[i] == '-') line2[i] = ' ';
                        sscanf(line2, "%s%s%s%s%d%lf%lf%lf", name1, name2,name3,name4, &tmpint, &tmpf1, &tmpf2, &tmpf3);
                        for(i=0;i<tornum;i++)  {
                                if(strcmp(name1, tor[i].name1) == 0 && strcmp(name2, tor[i].name2) == 0 && 
				   strcmp(name3, tor[i].name3) == 0 && strcmp(name4, tor[i].name4) == 0 && 
				   tor[i].fterm == tmpf3) {
                                        writeflag = 0;
                                        break;
                                }
                                if(strcmp(name1, tor[i].name4) == 0 && strcmp(name2, tor[i].name3) == 0 && 
				   strcmp(name3, tor[i].name2) == 0 && strcmp(name4, tor[i].name1) == 0 &&
                                   tor[i].fterm == tmpf3) {
                                        writeflag = 0;
                                        break;
                                }
                        }
                        if(writeflag == 1) {
                                strcpy(tor[tornum].name1, name1);
                                strcpy(tor[tornum].name2, name2);
                                strcpy(tor[tornum].name3, name3);
                                strcpy(tor[tornum].name4, name4);
				tor[tornum].fterm = tmpf3;
                                tornum++;
                        }
                }

                if(typeflag == 5) {
                        writeflag = 1;
			strcpy(line2, line);
			for(i=0;i<12;i++)  /* we do not want to delete '-' after 12th column*/
				if(line2[i] == '-') line2[i] = ' ';
                        sscanf(line2, "%s%s%s%s", name1, name2,name3,name4);
                        for(i=0;i<impropnum;i++)  {
                                if(strcmp(name1, improp[i].name1) == 0 && strcmp(name2, improp[i].name2) == 0 &&
                                   strcmp(name3, improp[i].name3) == 0 && strcmp(name4, improp[i].name4) == 0) {
                                        writeflag = 0;
                                        break;
                                }
                                if(strcmp(name1, improp[i].name4) == 0 && strcmp(name2, improp[i].name3) == 0 &&
                                   strcmp(name3, improp[i].name2) == 0 && strcmp(name4, improp[i].name1) == 0) {
                                        writeflag = 0;
                                        break;
                                }
                        }
                        if(writeflag == 1) {
                                strcpy(improp[impropnum].name1, name1);
                                strcpy(improp[impropnum].name2, name2);
                                strcpy(improp[impropnum].name3, name3);
                                strcpy(improp[impropnum].name4, name4);
                                impropnum++;
                        }
                }

		if(typeflag == 6) {
			writeflag = 1;
			for(i=0;i<vdwnum;i++) 
				if(strcmp(tmpchar, vdw[i].name) == 0) {
					writeflag = 0;
					break;
				}
			if(writeflag == 1) {
				strcpy(vdw[vdwnum].name, tmpchar);
				vdwnum++;
			}
		}

		if(writeflag == 1) 
			fprintf(fp2, "%s", line);
	}
	fclose(fp1);
	fclose(fp2);
}

int main(int argc, char *argv[])
{
	int i,j;
	int i_parm10 = 0;
	int i_parm99 = 0;
	FILE *fptmp;
	int overflow_flag = 0;			/*if overflow_flag ==1, reallocate memory */
	char frcmod_filename[MAXCHAR];
	char lowerstr[MAXCHAR];

	default_cinfo(&cinfo);
	default_minfo(&minfo);
    amberhome = (char *) getenv("AMBERCLASSICHOME");
    if( amberhome == NULL ){
       fprintf( stdout, "AMBERCLASSICHOME is not set!\n" );
       exit(1);
    }
    minfo.connect_file[0] = '\0';
    build_dat_path(minfo.connect_file, "CONNECT.TPL",
    	sizeof minfo.connect_file, 0);
	if (argc == 2 && (strcmp(argv[1], "-l") == 0 || strcmp(argv[1], "-L") == 0)) {
		printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
		printf("Format of Additional PARMCHK ATCS (atom type corresponding score) file\n");
		printf("# Each atom type definition starts with \"PARM\" index in the first\n");
		printf("#   four columns. Each atom types may have multiple \"EQUA\" atom types\n");
		printf("#   equivalent_flag:  1 for cc/ce/cg/nc/ne/pc/pe\n");
		printf("#                     2 for cd/cf/ch/nd/nf/pd/pf\n");
		printf("#                     0 for others\n");
		printf("+++++++++In the following example, np is defined which is equal to n3+++++++\n");
		printf("----------------------------------------------------------------------------\n");
		printf("#   atomtype    improper_flag group_id  mass    equivalent_flag  atomic_num\n");
		printf("PARM    np      1               1       14.01   0       7 \n");
		printf("EQUA    n3\n");
		printf("----------------------------------------------------------------------------\n");
		exit(0);
	}
	if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("[31mUsage: parmchk2 -i   [0m input file name\n"
		               "[31m                -o   [0m frcmod file name\n"
			       "[31m                -f   [0m input file format (prepi, prepc, ac, mol2, frcmod, leaplog) \n"
			       "[31m                -s   [0m ff parm set, it is suppressed by \"-p\" option\n"
			       "[34m                      1 or gaff:   [0m gaff (the default)\n"
			       "[34m                      2 or gaff2:  [0m gaff2 \n"
			       "[34m                      3 or parm99: [0m parm99\n"
			       "[34m                      4 or parm10: [0m parm10\n"
			       "[34m                      5 or lipid14:[0m lipid14\n"
			       "[31m                -frc [0m frcmod files to be loaded, the supported frcmods include\n" 
			       "[31m                     [0m ff99SB, ff14SB, ff03 for proteins , bsc1, ol15, ol3 for DNA and yil for RNA\n"
			       "[31m                     [0m eg. ff14SB+bsc1+yil, ff99SB+bsc1\n"
			       "[31m                -p   [0m parmfile, suppress '-s' flag, optional\n"
			       "[31m                -pf  [0m parmfile format \n"
			       "[34m                      1:[0m for amber FF data file (the default)\n"
			       "[34m                      2:[0m for additional force field parameter file\n"
			       "[31m                -afrc[0m additional frcmod file, no matter using -p or not, optional\n" 
			       "[31m                -c   [0m atom type corresponding score file, default is PARMCHK.DAT\n" 
			       "[31m                -atc [0m additional atom type corresponding score file, optional\n" 
			       "[31m                     [0m type 'parmchk2 -l' to learn details\n" 
			       "[31m                -a   [0m print out all force field parameters including those in the parmfile\n"
			       "[31m                     [0m can be 'Y' (yes) or 'N' (no) default is 'N' \n"
			       "[31m                -w   [0m print out parameters that matching improper dihedral parameters\n"
			       "[31m                     [0m that contain 'X' in the force field parameter file, can be 'Y' (yes)\n"
			       "[31m                     [0m or 'N' (no), default is 'Y'\n"
			       "[31m                -fc  [0m option of force constant calculation for '-f frcmod' or '-f leaplog'\n"
			       "[34m                      1:[0m default behavior (the default option)\n"
			       "[34m                      2:[0m do empirical calculation before using corresponding atom types\n"
			       "[31m                -att [0m for the frcmod input format, option of performing parmchk\n"
			       "[34m                      1:[0m for all parameters (the default)\n"
			       "[34m                      2:[0m only for those with ATTN \n");
			exit(1);
		}
		if (argc != 7 && argc != 9  && argc != 11 && argc != 13 && argc != 15 && argc != 17 && argc != 19
                              && argc != 21 && argc != 23 && argc != 25 && argc != 27 && argc != 29) {
			printf("[31mUsage: parmchk2 -i   [0m input file name\n"
		               "[31m                -o   [0m frcmod file name\n"
			       "[31m                -f   [0m input file format (prepi, prepc, ac, mol2, frcmod, leaplog) \n"
			       "[31m                -s   [0m ff parm set, it is suppressed by \"-p\" option\n"
			       "[34m                      1 or gaff:   [0m gaff (the default)\n"
			       "[34m                      2 or gaff2:  [0m gaff2 \n"
			       "[34m                      3 or parm99: [0m parm99\n"
			       "[34m                      4 or parm10: [0m parm10\n"
			       "[34m                      5 or lipid14:[0m lipid14\n"
			       "[31m                -frc [0m frcmod files to be loaded, the supported frcmods include\n" 
			       "[31m                     [0m ff99SB, ff14SB, ff03 for proteins , bsc1, ol15, ol3 for DNA and yil for RNA\n"
			       "[31m                     [0m eg. ff14SB+bsc1+yil, ff99SB+bsc1\n"
			       "[31m                -p   [0m parmfile, suppress '-s' flag, optional\n"
			       "[31m                -pf  [0m parmfile format \n"
			       "[34m                      1:[0m for amber FF data file (the default)\n"
			       "[34m                      2:[0m for additional force field parameter file\n"
			       "[31m                -afrc[0m additional frcmod file, no matter using -p or not, optional\n" 
			       "[31m                -c   [0m atom type corresponding score file, default is PARMCHK.DAT\n" 
			       "[31m                -atc [0m additional atom type corresponding score file, optional\n" 
			       "[31m                     [0m type 'parmchk2 -l' to learn details\n" 
			       "[31m                -a   [0m print out all force field parameters including those in the parmfile\n"
			       "[31m                     [0m can be 'Y' (yes) or 'N' (no) default is 'N' \n"
			       "[31m                -w   [0m print out parameters that matching improper dihedral parameters\n"
			       "[31m                     [0m that contain 'X' in the force field parameter file, can be 'Y' (yes)\n"
			       "[31m                     [0m or 'N' (no), default is 'Y'\n"
			       "[31m                -fc  [0m option of force constant calculation for '-f frcmod' or '-f leaplog'\n"
			       "[34m                      1:[0m default behavior (the default option)\n"
			       "[34m                      2:[0m do empirical calculation before using corresponding atom types\n"
			       "[31m                -att [0m for the frcmod input format, option of performing parmchk\n"
			       "[34m                      1:[0m for all parameters (the default)\n"
			       "[34m                      2:[0m only for those with ATTN \n");
			exit(1);
		}
	} else {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("Usage: parmchk2 -i    input file name\n"
		               "                -o    frcmod file name\n"
			       "                -f    input file format (prepi, prepc, ac, mol2, frcmod, leaplog) \n"
			       "                -s    ff parm set, it is suppressed by \"-p\" option\n"
			       "                      1 or gaff:    gaff (the default)\n"
			       "                      2 or gaff2:   gaff2 \n"
			       "                      3 or parm99:  parm99\n"
			       "                      4 or parm10:  parm10\n"
			       "                      5 or lipid14: lipid14\n"
			       "                -frc  frcmod files to be loaded, the supported frcmods include\n" 
			       "                      ff99SB, ff14SB, ff03 for proteins , bsc1, ol15, ol3 for DNA and yil for RNA\n"
			       "                      eg. ff14SB+bsc1+yil, ff99SB+bsc1\n"
			       "                -p    parmfile, suppress '-s' flag, optional\n"
			       "                -pf   parmfile format \n"
			       "                      1: for amber FF data file (the default)\n"
			       "                      2: for additional force field parameter file\n"
			       "                -afrc additional frcmod file, no matter using -p or not, optional\n" 
			       "                -c    atom type corresponding score file, default is PARMCHK.DAT\n" 
			       "                -atc  additional atom type corresponding score file, optional\n" 
			       "                      type 'parmchk2 -l' to learn details\n" 
			       "                -a    print out all force field parameters including those in the parmfile\n"
			       "                      can be 'Y' (yes) or 'N' (no) default is 'N' \n"
			       "                -w    print out parameters that matching improper dihedral parameters\n"
			       "                      that contain 'X' in the force field parameter file, can be 'Y' (yes)\n"
			       "                      or 'N' (no), default is 'Y'\n"
			       "                -fc   option of force constant calculation for '-f frcmod' or '-f leaplog'\n"
			       "                      1: default behavior (the default option)\n"
			       "                      2: do empirical calculation before using corresponding atom types\n"
			       "                -att  for the frcmod input format, option of performing parmchk\n"
			       "                      1: for all parameters (the default)\n"
			       "                      2: only for those with ATTN \n");
			exit(1);
		}
		if (argc != 7 && argc != 9  && argc != 11 && argc != 13 && argc != 15 && argc != 17 && argc != 19
                              && argc != 21 && argc != 23 && argc != 25 && argc != 27 && argc != 29) {
			printf("Usage: parmchk2 -i    input file name\n"
		               "                -o    frcmod file name\n"
			       "                -f    input file format (prepi, prepc, ac, mol2, frcmod, leaplog) \n"
			       "                -s    ff parm set, it is suppressed by \"-p\" option\n"
			       "                      1 or gaff:    gaff (the default)\n"
			       "                      2 or gaff2:   gaff2 \n"
			       "                      3 or parm99:  parm99\n"
			       "                      4 or parm10:  parm10\n"
			       "                      5 or lipid14: lipid14\n"
			       "                -frc  frcmod files to be loaded, the supported frcmods include\n" 
			       "                      ff99SB, ff14SB, ff03 for proteins , bsc1, ol15, ol3 for DNA and yil for RNA\n"
			       "                      eg. ff14SB+bsc1+yil, ff99SB+bsc1\n"
			       "                -p    parmfile, suppress '-s' flag, optional\n"
			       "                -pf   parmfile format \n"
			       "                      1: for amber FF data file (the default)\n"
			       "                      2: for additional force field parameter file\n"
			       "                -afrc additional frcmod file, no matter using -p or not, optional\n" 
			       "                -c    atom type corresponding score file, default is PARMCHK.DAT\n" 
			       "                -atc  additional atom type corresponding score file, optional\n" 
			       "                      type 'parmchk2 -l' to learn details\n" 
			       "                -a    print out all force field parameters including those in the parmfile\n"
			       "                      can be 'Y' (yes) or 'N' (no) default is 'N' \n"
			       "                -w    print out parameters that matching improper dihedral parameters\n"
			       "                      that contain 'X' in the force field parameter file, can be 'Y' (yes)\n"
			       "                      or 'N' (no), default is 'Y'\n"
			       "                -fc   option of force constant calculation for '-f frcmod' or '-f leaplog'\n"
			       "                      1: default behavior (the default option)\n"
			       "                      2: do empirical calculation before using corresponding atom types\n"
			       "                -att  for the frcmod input format, option of performing parmchk\n"
			       "                      1: for all parameters (the default)\n"
			       "                      2: only for those with ATTN \n");
			exit(1);
		}
	}
	iformat = -1;  /* input file format; -1 is an invalid format */
	cindex = 0;
	for (i = 1; i < argc; i += 2) {
		if (strcmp(argv[i], "-i") == 0)
			strcpy(ifilename, argv[i + 1]);
		if (strcmp(argv[i], "-p") == 0) {
			strcpy(pfilename, argv[i + 1]);
			ipfilename = 1;
		}
		if (strcmp(argv[i], "-pf") == 0)
			pformat = atoi(argv[i+1]);
		if (strcmp(argv[i], "-o") == 0)
			strcpy(ofilename, argv[i + 1]);
		if (strcmp(argv[i], "-c") == 0) {
			strcpy(cfilename, argv[i + 1]);
			cindex = 1;
		}
		if (strcmp(argv[i], "-f") == 0) {
			if (strcmp(argv[i + 1], "prepi") == 0)
				iformat = 0;
			if (strcmp(argv[i + 1], "prepc") == 0)
				iformat = 1;
			if (strcmp(argv[i + 1], "ac") == 0)
				iformat = 2;
			if (strcmp(argv[i + 1], "mol2") == 0)
				iformat = 3;
			if (strcmp(argv[i + 1], "frcmod") == 0)
				iformat = 4;
			if (strcmp(argv[i + 1], "frc") == 0)
				iformat = 4;
			if (strcmp(argv[i + 1], "leaplog") == 0)
				iformat = 5;
			if (strcmp(argv[i + 1], "leap") == 0)
				iformat = 5;
		}
		if (strcmp(argv[i], "-a") == 0) {
			if(argv[i + 1][0] == 'Y' ||argv[i + 1][0] == 'y') 
				allparm_flag = 1;
			if(argv[i + 1][0] == 'N' ||argv[i + 1][0] == 'n') 
				allparm_flag = 0;
		}
		if (strcmp(argv[i], "-w") == 0) {
			if(argv[i + 1][0] == 'Y' ||argv[i + 1][0] == 'y') 
				output_improper_flag = 1;
			if(argv[i + 1][0] == 'N' ||argv[i + 1][0] == 'n') 
				output_improper_flag = 0;
		}
		if (strcmp(argv[i], "-s") == 0) {
			for(j=0; j<strlen(argv[i + 1]); j++)	
				lowerstr[j] = tolower(argv[i + 1][j]);
			if (strcmp(lowerstr, "gaff") == 0 || strcmp(argv[i + 1], "1") == 0) 
				ffset = 1;
			if (strcmp(lowerstr, "gaff2") == 0 || strcmp(argv[i + 1], "2") == 0)
				ffset = 2;
			if (strcmp(lowerstr, "parm99") == 0 || strcmp(argv[i + 1], "3") == 0)
				ffset = 3;
			if (strcmp(lowerstr, "parm10") == 0 || strcmp(argv[i + 1], "4") == 0)
				ffset = 4;
			if (strcmp(lowerstr, "lipid14") == 0 || strcmp(argv[i + 1], "5") == 0)
				ffset = 5;
		}
		if (strcmp(argv[i], "-frc") == 0) {
			strcpy(frcmod_str, argv[i+1]);
			ifrcmod_str = 1;
		}
		if (strcmp(argv[i], "-afrc") == 0) {
			strcpy(additional_frcmod_file, argv[i+1]);
			iadditional_frcmod = 1;
		}
		if (strcmp(argv[i], "-fc") == 0) 
			fc_opt = atoi(argv[i + 1]);
		if (strcmp(argv[i], "-atc") == 0) {
			strcpy(atcfilename, argv[i + 1]);
			iatc = 1;
		} else if (strcmp(argv[i], "-att") == 0) {
			attn_opt = atoi(argv[i + 1]);
		} else if (strcmp(argv[i], "-at") == 0) {
			printf("Error: option -at does not exist !\n"
			"    You are probably looking for -s or mistyping"
			" options -atc or -att \n");
			exit(1);
		}
	}
	if(ipfilename == 0) {
		if(ffset != 1 && ffset != 2 && ffset !=3 && ffset !=4 && ffset !=5) 
			ffset = 1;
		if (ffset == 1) gaff_set = 1;
		if (ffset == 2) gaff_set = 2;
		if(gaff_set == 1) 
			build_path(pfilename, "/dat/leap/parm/", "gaff.dat",
				sizeof pfilename, 0);
		if(gaff_set == 2) 
			build_path(pfilename, "/dat/leap/parm/", "gaff2.dat",
				sizeof pfilename, 0);
		if (ffset == 3) {
			build_path(pfilename, "/dat/leap/parm/", "parm99.dat",
				sizeof pfilename, 0);
			gaff_set = 1;
		}
		if (ffset == 4) {
			build_path(pfilename, "/dat/leap/parm/", "parm10.dat",
				sizeof pfilename, 0);
			gaff_set = 1;
		}
		if (ffset == 5) {
			build_path(pfilename, "/dat/leap/parm/", "lipid14.dat",
				sizeof pfilename, 0);
			gaff_set = 1;
		}
		pformat = 1;
	}
	if(pformat != 1 && pformat != 2) pformat = 1;
	if(attn_opt != 1 && attn_opt !=2) attn_opt = 1;
	if(fc_opt != 1 && fc_opt !=2) fc_opt = 1;

	if (cindex == 0) {
		build_dat_path(cfilename, "PARMCHK.DAT", sizeof cfilename, 0);
		cindex = 1;
	}
	if (cindex == 0) {
		if ((fptmp = fopen("PARMCHK.DAT", "r")) != NULL) {
			strcpy(cfilename, "PARMCHK.DAT");
			cindex = 1;
			fclose(fptmp);
		}
		else {
			fprintf(stderr, "PARMCHK.DAT, the parmchk corresponding"
				" score file, does not exist; exiting.\n");
			exit(1);
		}
	}

	/*  memory allocation */
	/* initialize */
	maxparmnum = MAXPARM;
	maxatomtype = MAX_FF_ATOMTYPE;
	maxvdwparm = MAX_FF_VDW;
	maxbondparm = MAX_FF_BOND;
	maxangleparm = MAX_FF_ANGLE;
	maxtorsionparm = MAX_FF_TORSION;
	maximproperparm = MAX_FF_IMPROPER;

	atomtypenum = 0;
	vdwparmnum = 0;
	bondparmnum = 0;
	angleparmnum = 0;
	torsionparmnum = 0;
	improperparmnum = 0;


	/*read in prep or ac file */
	atom = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
	if (atom == NULL) {
		fprintf(stdout, "memory allocation error for *atom\n");
		exit(1);
	}
	bond_array = (BOND *) malloc(sizeof(BOND) * cinfo.maxbond);
	if (bond_array == NULL) {
		fprintf(stdout, "memory allocation error for *bond_array\n");
		exit(1);
	}
	for (i = 0; i < cinfo.maxbond; ++i) {
		bond_array[i].jflag = -1; /* bond type has not been assigned */
	}

	overflow_flag = 0;

	switch (iformat) {
	case 0:
		overflow_flag = rprepi(ifilename, &atomnum, atom, &bondnum,
			bond_array, &cinfo, &minfo);
		break;
	case 1:
		overflow_flag = rprepc(ifilename, &atomnum, atom, &cinfo, &minfo);
		break;
	case 2:
		overflow_flag = rac(ifilename, &atomnum, atom, &bondnum,
			bond_array, &cinfo, &minfo);
		break;
	case 3:
		overflow_flag = rmol2(ifilename, &atomnum, atom, &bondnum,
			bond_array, &cinfo, &minfo, 1);
		break;
	case 4:
		rfrcmod(ifilename, 0); 
	        r_atomtype = (ATOMTYPE *) calloc(r_atomtype_num, sizeof(ATOMTYPE));
        	if (r_atomtype == NULL) {
                	fprintf(stdout, "memory allocation error for *r_atomtype\n");
                	exit(1);
        	}

        	r_bondparm = (BOND_FF *) calloc(r_bondparm_num, sizeof(BOND_FF));
        	if (r_bondparm == NULL) {
                	fprintf(stdout, "memory allocation error for *r_bondparm\n");
                	exit(1);
        	}

        	r_angleparm = (ANGLE *) calloc(r_angleparm_num, sizeof(ANGLE));
        	if (r_angleparm == NULL) {
                	fprintf(stdout, "memory allocation error for *r_angleparm\n");
                	exit(1);
        	}

        	r_torsionparm = (TORSION *) calloc(r_torsionparm_num, sizeof(TORSION));
        	if (r_torsionparm == NULL) {
                	fprintf(stdout, "memory allocation error for *r_torsionparm\n");
                	exit(1);
        	}

        	r_improperparm = (IMPROPER *) calloc(r_improperparm_num, sizeof(IMPROPER));
        	if (r_improperparm == NULL) {
                	fprintf(stdout, "memory allocation error for *r_improperparm\n");
                	exit(1);
        	}

        	r_vdwparm = (VDW *) calloc(r_vdwparm_num, sizeof(VDW));
        	if (r_vdwparm == NULL) {
                	fprintf(stdout, "memory allocation error for *r_vdwparm\n");
                	exit(1);
        	}
		r_atomtype_num = 0;	
		r_bondparm_num = 0;	
		r_angleparm_num = 0;	
		r_torsionparm_num = 0;	
		r_improperparm_num = 0;	
		r_vdwparm_num = 0;	
		rfrcmod(ifilename, 1); 
		break;
	case 5:
		rleaplog(ifilename, 0); 
                r_atomtype = (ATOMTYPE *) calloc(r_atomtype_num, sizeof(ATOMTYPE));
                if (r_atomtype == NULL) {
                        fprintf(stdout, "memory allocation error for *r_atomtype\n");
                        exit(1);
                }

                r_bondparm = (BOND_FF *) calloc(r_bondparm_num, sizeof(BOND_FF));
                if (r_bondparm == NULL) {
                        fprintf(stdout, "memory allocation error for *r_bondparm\n");
                        exit(1);
                }

                r_angleparm = (ANGLE *) calloc(r_angleparm_num, sizeof(ANGLE));
                if (r_angleparm == NULL) {
                        fprintf(stdout, "memory allocation error for *r_angleparm\n");
                        exit(1);
                }

                r_torsionparm = (TORSION *) calloc(r_torsionparm_num, sizeof(TORSION));
                if (r_torsionparm == NULL) {
                        fprintf(stdout, "memory allocation error for *r_torsionparm\n");
                        exit(1);
                }

                r_improperparm = (IMPROPER *) calloc(r_improperparm_num, sizeof(IMPROPER));
                if (r_improperparm == NULL) {
                        fprintf(stdout, "memory allocation error for *r_improperparm\n");
                        exit(1);
                }

                r_vdwparm = (VDW *) calloc(r_vdwparm_num, sizeof(VDW));
                if (r_vdwparm == NULL) {
                        fprintf(stdout, "memory allocation error for *r_vdwparm\n");
                        exit(1);
                }
		r_atomtype_num = 0;	
		r_bondparm_num = 0;	
		r_angleparm_num = 0;	
		r_torsionparm_num = 0;	
		r_improperparm_num = 0;	
		r_vdwparm_num = 0;	
		rleaplog(ifilename, 1); 
		break;
	default:
		printf("Error: invalid input file format"
			"(valid: prepi, prepc, ac, mol2, frcmod, leaplog)!\n");
		exit(1);
		break;
}

	if (overflow_flag) {
		cinfo.maxatom = atomnum + 10;
                /* warning bondnum is not correct if overflow_flag != 0 */
                /* so add to maxbond a fudge based on the number of atoms */
                cinfo.maxbond = bondnum + 4*atomnum + 10;
		free(atom);
		free(bond_array);
		atom = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
		if (atom == NULL) {
			fprintf(stdout, "memory allocation error for *atom\n");
			exit(1);
		}
		bond_array = (BOND *) malloc(sizeof(BOND) * cinfo.maxbond);
		if (bond_array == NULL) {
			fprintf(stdout, "memory allocation error for *bond_array\n");
			exit(1);
		}
		int i;
		for (i = 0; i < cinfo.maxbond; ++i) {
			bond_array[i].jflag = -1; /* bond type has not been assigned */
		}
		if (iformat == 0)
			overflow_flag =
				rprepi(ifilename, &atomnum, atom, &bondnum, bond_array, &cinfo, &minfo);
		if (iformat == 1)
			overflow_flag =
				rprepc(ifilename, &atomnum, atom, &cinfo, &minfo);
		if (iformat == 2)
			overflow_flag =
				rac(ifilename, &atomnum, atom, &bondnum, bond_array, &cinfo,
					&minfo);
		if (iformat == 3)
			overflow_flag =
				rmol2(ifilename, &atomnum, atom, &bondnum, bond_array,
					  &cinfo, &minfo, 1);

	}

	if (iformat == 0) {
		atomicnum(atomnum, atom);
/*
		overflow_flag =
			connect(minfo.connect_file, atomnum, atom, &bondnum,
					bond_array, cinfo.maxbond);
*/
		adjustatomname(atomnum, atom, 1);
	}
	if (iformat == 1) {
		atomicnum(atomnum, atom);
		overflow_flag =
			connect(minfo.connect_file, atomnum, atom, &bondnum,
					bond_array, cinfo.maxbond);
		adjustatomname(atomnum, atom, 1);
	}
	if (iformat == 2) {
		atomicnum(atomnum, atom);
/*
		overflow_flag =
			connect(minfo.connect_file, atomnum, atom, &bondnum,
					bond_array, cinfo.maxbond);
*/
	}
	if (iformat == 3) {
		atomicnum(atomnum, atom);
		default_inf(atomnum, atom, 0);
	}

        parm = (PARM *) calloc(maxparmnum, sizeof(PARM));
        if (parm == NULL) {
                fprintf(stdout, "memory allocation error for *parm\n");
                exit(1);
        }
        parmid = (int *) malloc(sizeof(int) * atomnum);
        if (parmid == NULL) {
                fprintf(stdout, "memory allocation error for *parmid\n");
                exit(1);
        }
	improper = (IMPROPERID *) malloc(sizeof(IMPROPERID) * atomnum);
	if (improper == NULL) {
		fprintf(stdout, "memory allocation error for *improper\n");
		exit(1);
	}

	atomtype = (ATOMTYPE *) calloc(maxatomtype, sizeof(ATOMTYPE));
	if (atomtype == NULL) {
		fprintf(stdout, "memory allocation error for *atomtype\n");
		exit(1);
	}

	bondparm = (BOND_FF *) calloc(maxbondparm, sizeof(BOND_FF));
	if (bondparm == NULL) {
		fprintf(stdout, "memory allocation error for *bondparm\n");
		exit(1);
	}

	angleparm = (ANGLE *) calloc(maxangleparm, sizeof(ANGLE));
	if (angleparm == NULL) {
		fprintf(stdout, "memory allocation error for *angleparm\n");
		exit(1);
	}

	torsionparm = (TORSION *) calloc(maxtorsionparm, sizeof(TORSION));
	if (torsionparm == NULL) {
		fprintf(stdout, "memory allocation error for *torsionparm\n");
		exit(1);
	}

	improperparm = (IMPROPER *) calloc(maximproperparm, sizeof(IMPROPER));
	if (improperparm == NULL) {
		fprintf(stdout, "memory allocation error for *improperparm\n");
		exit(1);
	}

	vdwparm = (VDW *) calloc(maxvdwparm, sizeof(VDW));
	if (vdwparm == NULL) {
		fprintf(stdout, "memory allocation error for *vdwparm\n");
		exit(1);
	}

	if ((fpout = fopen(ofilename, "w")) == NULL) {
		fprintf(stdout, "Cannot open a file to write: %s, exit\n", ofilename);
		exit(1);
	}

	/* read in parameters */
	if (cindex == 1) {
		read_parmchk_parm(cfilename);	/*atom type corresponding file */
		parmnum2 = parmnum;
	}
	/* read in additional atom type corresponding file*/
	if (iatc == 1) 
		read_parmchk_parm_v2 (atcfilename);
	assign_parmid();
	if(ifrcmod_str == 1) {
		if(strstr(frcmod_str, "ff14SB") != NULL) {
			build_path(frcmod_filename, "/dat/leap/parm/", "frcmod.ff14SB",
				sizeof frcmod_filename, 0);
                        readfrcmod(frcmod_filename);
			last_atomtypenum = atomtypenum;
			last_bondparmnum = bondparmnum;
			last_angleparmnum = angleparmnum;
			last_torsionparmnum = torsionparmnum;
			last_improperparmnum = improperparmnum;
			last_vdwparmnum = vdwparmnum;
		}
		if(strstr(frcmod_str, "ff99SB") != NULL) {
			build_path(frcmod_filename, "/dat/leap/parm/", "frcmod.ff99SB",
				sizeof frcmod_filename, 0);
                        readfrcmod(frcmod_filename);
			last_atomtypenum = atomtypenum;
			last_bondparmnum = bondparmnum;
			last_angleparmnum = angleparmnum;
			last_torsionparmnum = torsionparmnum;
			last_improperparmnum = improperparmnum;
			last_vdwparmnum = vdwparmnum;
		}
		if(strstr(frcmod_str, "ff03") != NULL) {
			build_path(frcmod_filename, "/dat/leap/parm/", "frcmod.ff03",
				sizeof frcmod_filename, 0);
                        readfrcmod(frcmod_filename);
			last_atomtypenum = atomtypenum;
			last_bondparmnum = bondparmnum;
			last_angleparmnum = angleparmnum;
			last_torsionparmnum = torsionparmnum;
			last_improperparmnum = improperparmnum;
			last_vdwparmnum = vdwparmnum;
		}
		if(strstr(frcmod_str, "bsc1") != NULL) {
			build_path(frcmod_filename, "/dat/leap/parm/", "frcmod.parmbsc1",
				sizeof frcmod_filename, 0);
                        readfrcmod(frcmod_filename);
			last_atomtypenum = atomtypenum;
			last_bondparmnum = bondparmnum;
			last_angleparmnum = angleparmnum;
			last_torsionparmnum = torsionparmnum;
			last_improperparmnum = improperparmnum;
			last_vdwparmnum = vdwparmnum;
		}
		if(strstr(frcmod_str, "ol15") != NULL) {
			build_path(frcmod_filename, "/dat/leap/parm/", "frcmod.DNA.OL15",
				sizeof frcmod_filename, 0);
                        readfrcmod(frcmod_filename);
			last_atomtypenum = atomtypenum;
			last_bondparmnum = bondparmnum;
			last_angleparmnum = angleparmnum;
			last_torsionparmnum = torsionparmnum;
			last_improperparmnum = improperparmnum;
			last_vdwparmnum = vdwparmnum;
		}
		if(strstr(frcmod_str, "yil") != NULL) {
			build_path(frcmod_filename, "/dat/leap/parm/", "frcmod.parmCHI_YIL",
				sizeof frcmod_filename, 0);
                        readfrcmod(frcmod_filename);
			last_atomtypenum = atomtypenum;
			last_bondparmnum = bondparmnum;
			last_angleparmnum = angleparmnum;
			last_torsionparmnum = torsionparmnum;
			last_improperparmnum = improperparmnum;
			last_vdwparmnum = vdwparmnum;
		}
	}
	if(iadditional_frcmod == 1) {
		afrcmod_atomtypenum_beg = atomtypenum;
		afrcmod_bondparmnum_beg = bondparmnum;
		afrcmod_angleparmnum_beg = angleparmnum;
		afrcmod_torsionparmnum_beg = torsionparmnum;
		afrcmod_improperparmnum_beg = improperparmnum;
		afrcmod_vdwparmnum_beg = vdwparmnum;
		readfrcmod (additional_frcmod_file);
		afrcmod_atomtypenum_end = atomtypenum;
		afrcmod_bondparmnum_end = bondparmnum;
		afrcmod_angleparmnum_end = angleparmnum;
		afrcmod_torsionparmnum_end = torsionparmnum;
		afrcmod_improperparmnum_end = improperparmnum;
		afrcmod_vdwparmnum_end = vdwparmnum;
	}
	
	last_atomtypenum = atomtypenum;
	last_bondparmnum = bondparmnum;
	last_angleparmnum = angleparmnum;
	last_torsionparmnum = torsionparmnum;
	last_improperparmnum = improperparmnum;
	last_vdwparmnum = vdwparmnum;

	if(pformat == 1) readparm(pfilename);		/*principle parameter file */
	if(pformat == 2) readfrcmod(pfilename);	        /*principle parameter file in frcmod format*/

	if (iformat == 0 || iformat == 1)
		improper_id1(ifilename);
	if (iformat == 2 || iformat == 3)
		improper_id2();

	atomtypenum2 = atomtypenum; 
	vdwparmnum2 = vdwparmnum;
	bondparmnum2 = bondparmnum;
	angleparmnum2 = angleparmnum;
	torsionparmnum2 = torsionparmnum;
	improperparmnum2 = impropernum;

	if(iformat <= 3) {
/*
printf("chk atomtype ...\n");
fflush(stdout);
*/
		chk_atomtype();
/*
printf("chk bond ...\n");
fflush(stdout);
*/
		chk_bond();
/*
printf("chk angle ...\n");
fflush(stdout);
*/
		chk_angle();
/*
printf("chk torsion ...\n");
fflush(stdout);
*/
		chk_torsion();
/*
printf("chk improper ...\n");
fflush(stdout);
*/
		chk_improper();
/*
printf("chk vdw ...\n");
fflush(stdout);
*/
		chk_vdw();
	}
	if(iformat == 4 || iformat == 5) {
		allparm_flag = 1;
		chk_atomtype_v2();
		chk_bond_v2();
		chk_angle_v2();
		chk_torsion_v2();
		chk_improper_v2();
		chk_vdw_v2();
	}
	fclose(fpout);
	if(allparm_flag == 1 && iformat <= 3) 
		cleanup_frcmod(ofilename);
/*
		free(atom);
		free(bond_array);
		free(corrname);
		free(similarname);
		free(similarname2);
		free(corrindex);
		free(similarindex);
		free(similarindex2);
		free(improper);
		free(improperindex);
		free(corr);
		free(atomtype);
		free(bondparm);
		free(angleparm);
		free(torsionparm);
		free(improperparm);
		free(vdwparm);
*/
	return (0);
}
