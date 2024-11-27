/*
************************************************************************
*           All Copyright Reserved!                                    *
*                                                                      *
*  Prog:    bcc                                                        *
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
# include <math.h>
# include "define.h"
# include "atom.h"
# include "eprintf.h"
# include "common.c"
# include "rotate.c"
# include "ac.c"
# include "pdb.c"
# define MAX_BCCPARM 1024
# define MAX_BCCCORR 50
# define debug 0

ATOM *atom;
BOND *bond;
int atomnum = 0;
int bondnum = 0;
MOLINFO minfo;
CONTROLINFO cinfo;

typedef struct {
char ati[10];
char atj[10];
int    bt;
double bcc;
} BCCPARM;

typedef struct {
char at[10];
char ref[10];
} BCCCORR;

FILE *fpparm;
FILE *fpin;
FILE *fpout;

char ifilename[MAXCHAR];
char ofilename[MAXCHAR];
char pfilename[MAXCHAR] = "BCCPARM.DAT";
char line[MAXCHAR];
char am1bcc_type[MAXCHAR]="bcc";
int i, j, k, l;
int intstatus = 1;
BCCPARM bcc_parm[MAX_BCCPARM];
BCCCORR bcc_corr[MAX_BCCCORR];
int num_bcc_parm=0;
int num_bcc_corr=0;
int  i_am1bcc_type = 0;
double totalcharge = 0.0;


void rparm(char *filename)
{
    int i, j, k;
    int tmpint; 
    char tmpchar[MAXCHAR];
    char line[MAXCHAR];
    FILE *fpin;

    fpin = efopen(filename, "r");

    for (;;) {
        if (fgets(line, MAXCHAR, fpin) == NULL)
            break;
	if(strncmp(line, "CORR", 4) == 0) { 
		sscanf(line, "%s%s%s", tmpchar, bcc_corr[num_bcc_corr].at, bcc_corr[num_bcc_corr].ref);
		num_bcc_corr++;
	}
	else {
        	sscanf(line, "%d%s%s%d%lf", &tmpint,  bcc_parm[num_bcc_parm].ati, bcc_parm[num_bcc_parm].atj,
                                                     &bcc_parm[num_bcc_parm].bt, &bcc_parm[num_bcc_parm].bcc); 
		num_bcc_parm++;
	}
    }
    if (debug > 1) {
        for (i = 0; i < num_bcc_corr; i++) 
            printf("\nCORR %5d%5s%5s", i+1, bcc_corr[i].at, bcc_corr[i].ref);
        for (i = 0; i < num_bcc_parm; i++) 
            printf("\nPARM %5d%5s%5s%5d%10.4lf", i+1, bcc_parm[i].ati, bcc_parm[i].atj, bcc_parm[i].bt, bcc_parm[i].bcc);
    }
}

double chk_bcc_parm(char *ati, char *atj, int bt) {
double bcc=-999.9;
int i;
        for(i=0; i<num_bcc_parm; i++) {
                if(strcmp(ati, bcc_parm[i].ati) == 0 &&
                   strcmp(atj, bcc_parm[i].atj) == 0 &&
                   bt == bcc_parm[i].bt) {
                        bcc = bcc_parm[i].bcc;
                        break;
                }
                if(strcmp(ati, bcc_parm[i].atj) == 0 &&
                   strcmp(atj, bcc_parm[i].ati) == 0 &&
                   bt == bcc_parm[i].bt) {
                        bcc =-bcc_parm[i].bcc;
                        break;
                }
        }
return bcc;
}
double bccparm(int ati, int atj, int bt) {
double bcc=0;
char ati_str[10];
char atj_str[10]; 
char ati_corr_str[10];
char atj_corr_str[10]; 
int  icorr_ati = 0;
int  icorr_atj = 0;
int suc = 0;
	strcpy(ati_str, atom[ati].ambername);
	strcpy(atj_str, atom[atj].ambername);
	bcc = chk_bcc_parm(ati_str, atj_str, bt); 
	if(bcc > -999.0) suc = 1;	

	if(suc == 0 && num_bcc_corr > 0) {
		for(i=0;i<num_bcc_corr;i++) 
			if(strcmp(bcc_corr[i].at, ati_str) == 0) {		
				strcpy(ati_corr_str, bcc_corr[i].ref); 
				icorr_ati = 1;
				break;
			}
		for(i=0;i<num_bcc_corr;i++) 
			if(strcmp(bcc_corr[i].at, atj_str) == 0) {		
				strcpy(atj_corr_str, bcc_corr[i].ref); 
				icorr_atj = 1;
				break;
			}
		if(suc == 0 && icorr_ati == 1) {
			bcc = chk_bcc_parm(ati_corr_str, atj_str, bt); 
			if(bcc > -999.0) suc = 1;	
		}
		if(suc == 0 && icorr_atj == 1) {
			bcc = chk_bcc_parm(ati_str, atj_corr_str, bt); 
			if(bcc > -999.0) suc = 1;	
		}
		if(suc == 0 && icorr_ati == 1 && icorr_atj == 1) {
			bcc = chk_bcc_parm(ati_corr_str, atj_corr_str, bt); 
			if(bcc > -999.0) suc = 1;	
		}
	}
	if(suc == 0) bcc = 0;
return bcc;
}
void charge(void)
{
    int i, j, k, l; 
    double bcc;
    if (debug > 0)
        printf
            ("\n bond  at1  at2       pre-bcc       correction(code)         post-bcc\n");
    for (l = 0; l < bondnum; l++) {
        i = bond[l].bondi;
        j = bond[l].bondj;
        k = bond[l].type;
	bcc=bccparm(i,j,k);
	
        if (debug > 0) {
            printf("%4d %4s %4s  %8.4lf %8.4lf   %8.4lf (%02s-%02s-%02d)   %8.4lf %8.4lf\n",
                   l, atom[i].name, atom[j].name, atom[i].charge, atom[j].charge,
                   bcc, atom[i].ambername, atom[j].ambername, k, atom[i].charge + bcc, atom[j].charge - bcc);

        }
        atom[i].charge += bcc;
        atom[j].charge -= bcc;
    }

    if (debug > 1)
        for (i = 0; i < atomnum; i++)
            printf("\n%5d%9.4lf", i + 1, atom[i].charge);
}


int main(int argc, char *argv[])
{
    int i;
    int index;
    int format = 0;
    int judge_flag = 0;
    char command_at[2 * MAXCHAR];
    char command_bt[2 * MAXCHAR];
    int overflow_flag = 0;      /*if overflow_flag ==1, reallocate memory */

    amberhome = egetenv("AMBERCLASSICHOME");
    default_cinfo(&cinfo);
    default_minfo(&minfo);

    if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
        if (argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-H") == 0)) {
            printf("[31mUsage: am1bcc -i[0m input file name in ac format \n"
                   "[31m              -o[0m output file name \n"
                   "[31m              -f[0m output file format(pdb or ac, optional, default is ac)\n"
                   "[31m              -t[0m am1bcc type, default is 'bcc', can also be 'abcg2'\n"
                   "[31m              -p[0m bcc parm file name (optional))\n"
                   "[31m              -s[0m status information, can be 0 (brief), 1 (the default) and 2 (verbose)\n"
                   "[31m              -j[0m atom and bond type judge option, default is 0)\n"
                   "[32m                 0[0m: No judgement\n"
                   "[32m                 1[0m: Atom type\n"
                   "[32m                 2[0m: Full bond type\n"
                   "[32m                 3[0m: Partial bond type\n"
                   "[32m                 4[0m: Atom and full bond type\n"
                   "[32m                 5[0m: Atom and partial bond type\n");
            exit(0);
        }
        if (argc != 7 && argc != 9 && argc != 11 && argc != 13 && argc != 15) {
            printf("[31mUsage: am1bcc -i[0m input file name in ac format \n"
                   "[31m              -o[0m output file name \n"
                   "[31m              -f[0m output file format(pdb or ac, optional, default is ac)\n"
                   "[31m              -t[0m am1bcc type, default is 'bcc', can also be 'abcg2'\n"
                   "[31m              -p[0m bcc parm file name (optional))\n"
                   "[31m              -s[0m status information, can be 0 (brief), 1 (the default) and 2 (verbose)\n"
                   "[31m              -j[0m atom and bond type judge option, default is 0)\n"
                   "[32m                 0[0m: No judgement\n"
                   "[32m                 1[0m: Atom type\n"
                   "[32m                 2[0m: Full bond type\n"
                   "[32m                 3[0m: Partial bond type\n"
                   "[32m                 4[0m: Atom and full bond type\n"
                   "[32m                 5[0m: Atom and partial bond type\n");
            exit(1);
        }
    } else {
        if (argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-H") == 0)) {
            printf("Usage: am1bcc -i input file name in ac format \n"
                   "              -o output file name \n"
                   "              -f output file format(pdb or ac, optional, default is ac)\n"
                   "              -t am1bcc type, default is 'bcc', can also be 'abcg2'\n"
                   "              -p bcc parm file name (optional))\n"
                   "              -s status information, can be 0 (brief), 1 (the default) and 2 (verbose)\n"
                   "              -j atom and bond type judge option, default is 0)\n"
                   "                 0: No judgement\n" "                 1: Atom type\n"
                   "                 2: Full bond type\n"
                   "                 3: Partial bond type\n"
                   "                 4: Atom and full bond type\n"
                   "                 5: Atom and partial bond type\n");
            exit(0);
        }
        if (argc != 7 && argc != 9 && argc != 11 && argc != 13 && argc != 15) {
            printf("Usage: am1bcc -i input file name in ac format \n"
                   "              -o output file name \n"
                   "              -f output file format(pdb or ac, optional, default is ac)\n"
                   "              -t am1bcc type, default is 'bcc', can also be 'abcg2'\n"
                   "              -p bcc parm file name (optional))\n"
                   "              -s status information, can be 0 (brief), 1 (the default) and 2 (verbose)\n"
                   "              -j atom and bond type judge option, default is 0)\n"
                   "                 0: No judgement\n" "                 1: Atom type\n"
                   "                 2: Full bond type\n"
                   "                 3: Partial bond type\n"
                   "                 4: Atom and full bond type\n"
                   "                 5: Atom and partial bond type\n");
            exit(1);
        }

    }

    index = 0;
    for (i = 1; i < argc; i += 2) {
        if (strcmp(argv[i], "-i") == 0)
            strcpy(ifilename, argv[i + 1]);
        if (strcmp(argv[i], "-o") == 0)
            strcpy(ofilename, argv[i + 1]);
        if (strcmp(argv[i], "-t") == 0)
            strcpy(am1bcc_type, argv[i + 1]);
        if (strcmp(argv[i], "-j") == 0)
            judge_flag = atoi(argv[i + 1]);
        if (strcmp(argv[i], "-s") == 0)
            intstatus = atoi(argv[i + 1]);
        if (strcmp(argv[i], "-f") == 0) {
            if (strcmp(argv[i + 1], "PDB") == 0 || strcmp(argv[i + 1], "pdb") == 0) {
                format = 1;
            } else {
                format = 0;
            }
        }
        if (strcmp(argv[i], "-p") == 0) {
            strcpy(pfilename, argv[i + 1]);
            index = 1;
        }
    }

    i_am1bcc_type = 0;
    if(strcmp(am1bcc_type, "abcg2") == 0 || strcmp(am1bcc_type, "ABCG2") == 0) 
    i_am1bcc_type = 1;

    if (index == 0) {
        pfilename[0] = '\0';
	if(i_am1bcc_type == 1)
        	build_dat_path(pfilename, "BCCPARM_ABCG2.DAT", sizeof pfilename, 0);
	else
        	build_dat_path(pfilename, "BCCPARM.DAT", sizeof pfilename, 0);
    }
    command_at[0] = '\0';
    command_bt[0] = '\0';
    build_exe_path(command_at, "atomtype", sizeof command_at, 1);
    if(i_am1bcc_type == 1)
    	strcat(command_at, " -f ac -p abcg2 -o ANTECHAMBER_AM1BCC.AC -i ");
    else
    	strcat(command_at, " -f ac -p bcc -o ANTECHAMBER_AM1BCC.AC -i ");
    build_exe_path(command_bt, "bondtype", sizeof command_bt, 1);
    strcat(command_bt, " -f ac -o ANTECHAMBER_AM1BCC.AC -i ");

    if (judge_flag == 1) {
        strcat(command_at, ifilename);
        if (intstatus == 2)
            fprintf(stdout, "\nRunning: %s\n", command_at);
        esystem(command_at);
    }
    if (judge_flag == 2) {
        strcat(command_bt, ifilename);
        strcat(command_bt, " -j full");
        if (intstatus == 2)
            fprintf(stdout, "\nRunning: %s\n", command_bt);
        esystem(command_bt);
    }
    if (judge_flag == 3) {
        strcat(command_bt, ifilename);
        if (intstatus == 2)
            fprintf(stdout, "\nRunning: %s\n", command_bt);
        esystem(command_bt);
    }
    if (judge_flag == 4) {
        strcat(command_bt, ifilename);
        strcat(command_bt, " -j full");
        if (intstatus == 2)
            fprintf(stdout, "\nRunning: %s\n", command_bt);
        esystem(command_bt);
        strcat(command_at, "ANTECHAMBER_AM1BCC.AC");
        if (intstatus == 2)
            fprintf(stdout, "\nRunning: %s\n", command_at);
        esystem(command_at);
    }
    if (judge_flag == 5) {
        strcat(command_bt, ifilename);
        if (intstatus == 2)
            fprintf(stdout, "\nRunning: %s\n", command_bt);
        esystem(command_bt);
        strcat(command_at, "ANTECHAMBER_AM1BCC.AC");
        if (intstatus == 2)
            fprintf(stdout, "\nRunning: %s\n", command_at);
        esystem(command_at);
    }

    if (judge_flag != 0)
        strcpy(ifilename, "ANTECHAMBER_AM1BCC.AC");
/*allocate memory*/
    atom = (ATOM *) emalloc(sizeof(ATOM) * cinfo.maxatom);
    bond = (BOND *) emalloc(sizeof(BOND) * cinfo.maxbond);
    for (i = 0; i < cinfo.maxbond; ++i) {
        bond[i].jflag = -1;     /* bond type has not been assigned */
    }
/*read ac file*/
    overflow_flag = rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
    if (overflow_flag) {
        cinfo.maxatom = atomnum + 10;
        cinfo.maxbond = bondnum + 10;
        free(atom);
        free(bond);
        atom = (ATOM *) emalloc(sizeof(ATOM) * cinfo.maxatom);
        bond = (BOND *) emalloc(sizeof(BOND) * cinfo.maxbond);
        int i;
        for (i = 0; i < cinfo.maxbond; ++i) {
            bond[i].jflag = -1; /* bond type has not been assigned */
        }
        overflow_flag = rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
    }
    atomicnum(atomnum, atom);
    adjustatomname(atomnum, atom, 1);

    rparm(pfilename);
    charge();
    if (format == 0)
        wac(ofilename, atomnum, atom, bondnum, bond, cinfo, minfo);
    if (format == 1)
        wpdb(ofilename, atomnum, atom);
/*
	 free(atom);
	 free(bond);
*/
    return (0);
}
