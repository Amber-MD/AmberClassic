/*
************************************************************************
*           All Copyright Reserved!                                    *
*                                                                      *
*  Prog:    respgen                                                    *
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
# include "ac.c"
# define MAXCONF 20
# define MAXCHARGE 20
# define MAXSEPBOND 100 /*one one defines more than 100 separating bonds*/

typedef struct {
	char name[20];
	double charge;
	int id;
} CHARGE; 

typedef struct {
	char resat[20];
	char capat[20];
	char atomtype[20];
	char atomtype_cap[20];
	int  resid;
	int  capid;
	int  type;
} SEPBOND; 

ATOM *atom;
BOND *bond;
int atomnum = 0;
int bondnum = 0;
int *flag;
CONTROLINFO cinfo;
MOLINFO minfo;
int i, j, k;
FILE *fpin;
FILE *fpout;
char line[MAXCHAR];
char ifilename[MAXCHAR];

int confnum = 1;
SEPBOND sepbond[MAXSEPBOND];
int sepbondnum = 0;
CHARGE charge[MAXCHARGE];
int chargenum = 0;
int residueatnum ;
double netcharge = 0;
int i_netcharge = 0;
char inputfile[MAXCHAR];
char espfile[MAXCHAR];
char prepfile[MAXCHAR];
char atom_type_set[MAXCHAR] = "default";
char residue_symbol[MAXCHAR] = "MOL";
char residue_file_name[MAXCHAR] = "MOL.res";
char head_atom_name[MAXCHAR]="UN";
char tail_atom_name[MAXCHAR]="UN";
char head_atom_type[MAXCHAR]="DU";
char tail_atom_type[MAXCHAR]="DU";
char pre_head_atom_type[MAXCHAR]="DU";
char post_tail_atom_type[MAXCHAR]="DU";

int *selectelement;
int *selectelement2;
int selectnum = 0;
int selectnum2 = 0;
int pre_head_id = -1;
int post_tail_id = -1;
int head_id = -1;
int tail_id = -1;
int length;
int start_atomid;

void findpath(ATOM atm[], int selectnum, int startnum) {
int i, j;
int start;
int resetindex = -1;
start = -1;
selectelement[selectnum++] = startnum;
for (i = 0; i < atm[startnum].connum; i++) {
	if(atm[startnum].con[i] == -1) return;
	if(atm[startnum].con[i] == tail_id) {
		selectelement[selectnum++] = tail_id;
		for(j = 0; j < selectnum; j++)
                	selectelement2[j] = selectelement[j];
		selectnum2 = selectnum ;
		return;
	}
        start = atm[startnum].con[i];
        for (j = 0; j < selectnum; j++)
                if (start == selectelement[j]) {
                        resetindex = 1;
                        break;
                }
	if(resetindex == 1) {
		resetindex = -1;
		continue;
	}
	if (start == -1) return;
        findpath(atm, selectnum, start);
}
}


void group_atom1(int id) {
        int i, j, num, num_old;
        int atid;
        num = 1;
        num_old = 0;

       	for (i = 0; i<atomnum; i++)
             	flag[i] = -1;
        flag[sepbond[id].capid] = 1;
        while (num_old < num) {
                num_old = num;
                for (i = 0; i<atomnum; i++)
                        if(flag [i] != -1)
                        for(j =0 ; j<6; j++){
                                atid = atom[i].con[j];
                                if(atid == -1) continue;
                               	if(atid == sepbond[id].resid) continue;
				if(flag[atid] != -1) continue;
                                flag[atid] = 1;
                                num ++;
                        }
       }
}

void group_atom2() {
        int i, j, k, num, num_old;
        int atid;
	int break_flag ;
        num = 1;
        num_old = 0;

       	for (i = 0; i<atomnum; i++)
             	flag[i] = -1;
	for(i=0;i<sepbondnum;i++) 
                flag[sepbond[i].resid] = 1;

        while (num_old < num) {
                num_old = num;
                for (i = 0; i<atomnum; i++)
                        if(flag [i] != -1)
                        for(j =0 ; j<6; j++){
                                atid = atom[i].con[j];
                                if(atid == -1) continue;
				break_flag = 0;
				for(k=0;k<sepbondnum;k++)
                                	if(atid == sepbond[k].capid) {
						break_flag = 1; 
						break; 
					}
				if(break_flag == 1) continue;
                                if(flag[atid] != -1) continue;
                                flag[atid] = 1;
                                num ++;
                        }
                }
}

void rcharge(char *filename, int atomnum, ATOM atom[])
{
	FILE *fpcharge;
	char line[MAXCHAR];
	int i,j;
	int flag;
	int break_flag = 0;
	int ncol = 0;
	double chg[10];

	if ((fpcharge = fopen(filename, "r")) == NULL) {
		fprintf(stdout, "Cannot open charge file to read: %s , exit\n", filename);
		exit(1);
	}
	i = 0;
	flag = 0;
	for (;;) {
		if (fgets(line, MAXCHAR, fpcharge) == NULL) {
/*       printf("\nFinished reading file %s", filename); */
			break;
		}
		sscanf(line, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &chg[0], &chg[1], &chg[2], 
			   &chg[3], &chg[4], &chg[5], &chg[6], &chg[7], &chg[8], &chg[9]);
//	check how many colums per line
		if(flag == 0) {
			for(j = 0; j<strlen(line); j++)
				if(line[j] == '.') ncol++;
			flag = 1;
		}
		for(j=0; j< ncol; j++) {
			if(i==atomnum - j) {
				break_flag = 1;
				break;
			}
			atom[i+j].charge = chg[j];
		}
		if(break_flag == 1) 
			break;
		else 	
			i = i + ncol;
	}
	fclose(fpcharge);
}

void adjust_charge() {
double posc= 0, negc = 0, totc = 0, nc = 0;
double diff;
int i;
	for(i=0; i<atomnum; i++)	
		if(flag[i] == 1) {
			nc += atom[i].charge;
			totc += fabs(atom[i].charge);
			if(atom[i].charge > 0) posc+=atom[i].charge;
			if(atom[i].charge < 0) negc+=atom[i].charge;
		}
	if(totc == 0) return;
	if(nc == netcharge) return;
	diff = (nc - netcharge)/totc;
	for(i=0; i<atomnum; i++) 
		if(flag[i] == 1) {
			if(diff > 0) atom[i].charge += atom[i].charge * diff; 	
			if(diff < 0) atom[i].charge -= atom[i].charge * diff; 	
		}
}

void proceed1() {
	char command[MAXCHAR];
	char tmpchar[MAXCHAR];

	printf("[31mStep 1:  Generate RESP input files ... [0m\n");
	printf("\nRun respgen to generate first stage resp input file ... \n");
	build_exe_path(command, "respgen", sizeof command, 1);
	strcat(command, " -i ");
	strcat(command, inputfile);
	sprintf(tmpchar, "%d", confnum);
	strcat(command, " -o RESIDUE_GEN_RESP.INPUT1 -f resp1 -a RESIDUE_GEN_RESPGEN.DAT -n ");
	strcat(command, tmpchar);
	printf("\nCommand: %s\n", command);
	system(command);
	fflush(stdout);

	printf("\nRun respgen to generate second stage resp input file ... \n");
	build_exe_path(command, "respgen", sizeof command, 1);
	strcat(command, " -i ");
	strcat(command, inputfile);
	sprintf(tmpchar, "%d", confnum);
	strcat(command, " -o RESIDUE_GEN_RESP.INPUT2 -f resp2 -a RESIDUE_GEN_RESPGEN.DAT -n ");
	strcat(command, tmpchar);
	printf("\nCommand: %s\n", command);
	system(command);
	fflush(stdout);

	printf("\n\n[31mStep 2:  Run resp to get RESP charges ... [0m\n");
	build_exe_path(command, "resp", sizeof command, 1);
	strcat(command, " -O -i RESIDUE_GEN_RESP.INPUT1 -o RESIDUE_GEN_RESP.OUTPUT1 -q QIN -e ");
	strcat(command, espfile);
	printf("\nCommand: %s\n", command);
	system(command);
	fflush(stdout);

	build_exe_path(command, "resp", sizeof command, 1);
	strcat(command, " -O -i RESIDUE_GEN_RESP.INPUT2 -o RESIDUE_GEN_RESP.OUTPUT2 -q qout -e ");
	strcat(command, espfile);
	printf("\nCommand: %s", command);
	system(command);
	fflush(stdout);

	printf("[31mStep 3:  Read RESP charges ... [0m\n");
	rcharge("qout", atomnum, atom) ;
	wac("RESIDUE_GEN.AC", atomnum, atom, bondnum, bond, cinfo, minfo);

	printf("\n[31mStep 4:  Generate prep input file ... [0m\n");
	build_exe_path(command, "prepgen", sizeof command, 1);
	strcat(command, " -i RESIDUE_GEN.AC -o ");
	strcat(command, prepfile);
	strcat(command, " -m RESIDUE_GEN_MAINCHAIN.DAT -rf "); 
	strcat(command, residue_file_name);
	strcat(command, " -rn "); 
	strcat(command, residue_symbol);
	printf("\nCommand: %s\n", command);
	fflush(stdout);
	system(command);
	printf("\n\n[31mAll Done ![0m\n");
}

void proceed2() {
	char command[MAXCHAR];

	if(i_netcharge == 1) adjust_charge();
	printf("[31mStep 1:  Generate AC file ... [0m\n");
	wac("RESIDUE_GEN.AC", atomnum, atom, bondnum, bond, cinfo, minfo);

	printf("\n[31mStep 2:  Generate prep input file ... [0m\n");
	build_exe_path(command, "prepgen", sizeof command, 1);
	strcat(command, " -i RESIDUE_GEN.AC -o ");
	strcat(command, prepfile);
	strcat(command, " -m RESIDUE_GEN_MAINCHAIN.DAT -rf "); 
	strcat(command, residue_file_name);
	strcat(command, " -rn "); 
	strcat(command, residue_symbol);
	printf("\nCommand: %s\n", command);
	fflush(stdout);
	system(command);
	printf("\n\n[31mAll Done ![0m\n");
}

int main(int argc, char *argv[])
{
	int i,j;
	int index;
	int num_sep_bond_type1 = 0;
	int num_sep_bond_type2 = 0;
	int overflow_flag = 0;			/*if overflow_flag ==1, reallocate memory */
	char command[MAXCHAR];

    amberhome = (char *) getenv("AMBERCLASSICHOME");
    if( amberhome == NULL ){
       fprintf( stdout, "AMBERCLASSICHOME is not set!\n" );
       exit(1);
    }
	if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("[31mUsage: residuegen [0m input_file \n");
			exit(1);
		}
		if (argc != 2) {
			printf("[31mUsage: residuegen [0m input_file \n");
			exit(1);
		}
	} else {
		if (argc == 2
			&& (strcmp(argv[1], "-h") == 0
				|| strcmp(argv[1], "-H") == 0)) {
			printf("Usage residuegen  input_file\n");
			exit(1);
		}
		if (argc != 2) {
			printf("Usage residuegen  input_file\n");
			exit(1);
		}
	}
        if ((fpin = fopen(argv[1], "r")) == NULL) {
                fprintf(stdout, "Cannot open input file: %s, exit\n", argv[1]);
                exit(1);
        }
	
	for(i=0; i<MAXSEPBOND; i++) {
		sepbond[i].type = 0;
		strcpy(sepbond[i].atomtype, "DU"); 
		strcpy(sepbond[i].atomtype_cap, "DU"); 
	}
        for (;;) {
                if (fgets(line, MAXCHAR, fpin) == NULL) break;
                if (strncmp("ESP_FILE", line, 8) == 0) {
			sscanf(&line[9], "%s", espfile); 
                        continue;
                }
                if (strncmp("INPUT_FILE", line, 10) == 0) {
			sscanf(&line[11], "%s", inputfile); 
                        continue;
                }
                if (strncmp("CONF_NUM", line, 8) == 0) {
			sscanf(&line[9], "%d", &confnum); 
                        continue;
                }
                if (strncmp("ATOM_TYPE_SET", line, 13) == 0) {
			sscanf(&line[13], "%s", &atom_type_set); 
                        continue;
                }
                if (strncmp("SEP_BOND", line, 8) == 0) {
			sscanf(&line[9], "%s%s%d%s%s", sepbond[sepbondnum].resat, sepbond[sepbondnum].capat, 
				&sepbond[sepbondnum].type, sepbond[sepbondnum].atomtype, sepbond[sepbondnum].atomtype_cap); 
			sepbondnum++;
			if(sepbondnum >= MAXSEPBOND) {
				fprintf(stderr, "Too many sepbond defined\n");
				exit(1);
			}
                        continue;
                }
                if (strncmp("HEAD_ATOM", line, 9) == 0) {
			sscanf(&line[9], "%s%s%s", head_atom_name, head_atom_type, pre_head_atom_type); 
                        continue;
                }
                if (strncmp("TAIL_ATOM", line, 9) == 0) {
			sscanf(&line[9], "%s%s%s", tail_atom_name, tail_atom_type, post_tail_atom_type); 
                        continue;
                }
                if (strncmp("NET_CHARGE", line, 10) == 0) {
			sscanf(&line[11], "%lf", &netcharge); 
			i_netcharge = 1;
                        continue;
                }
                if (strncmp("PREP_FILE", line, 9) == 0) {
			sscanf(&line[10], "%s", prepfile); 
                        continue;
                }
                if (strncmp("RESIDUE_FILE_NAME", line, 17) == 0) {
			sscanf(&line[18], "%s", residue_file_name); 
                        continue;
                }
                if (strncmp("RESIDUE_SYMBOL", line, 14) == 0) {
			sscanf(&line[15], "%s", residue_symbol); 
                        continue;
                }
                if (strncmp("ATOM_CHARGE", line, 11) == 0) {
			if ( chargenum >= MAXCHARGE ) {
                        fprintf(stdout, "Error:  MAXCHARGE is too small in residuegen.c\n"
                                        "  Increase MAXCHARGE and reinstall antechamber\n"
                               );
                        exit(1);
			}
			sscanf(&line[12], "%s%lf", charge[chargenum].name, &charge[chargenum].charge); 
			chargenum++;
                        continue;
                }
	}
	fclose(fpin);

	printf("\nINPUT_FILE:\t\t%s", inputfile);
	printf("\nESP_FILE:\t\t%s", espfile);
	printf("\nCONF_NUM:\t\t%d", confnum);
	printf("\nNET_CHARGE:\t\t%-9.5lf", netcharge);
	if(strcmp(head_atom_name, "HEAD_ATOM") != 0)
		printf("\nHEAD_ATOM:\t%s\t%s\t%s", head_atom_name, head_atom_type, pre_head_atom_type); 
	if(strcmp(tail_atom_name, "TAIL_ATOM") != 0)
		printf("\nTAIL_ATOM:\t%s\t%s\t%s", tail_atom_name, tail_atom_type, post_tail_atom_type); 
	printf("\nThere are %d separating bonds", sepbondnum);
	for(i=0;i<sepbondnum;i++)
		printf("\nSEPBOND\t\t\t%d\t%s\t%s\t%d\t%s\t%s", i+1, sepbond[i].resat, sepbond[i].capat, sepbond[i].type, sepbond[i].atomtype, sepbond[i].atomtype_cap); 
	printf("\nThere are %d predefined charges", chargenum);
	for(i=0;i<chargenum;i++) 
			printf("\nATOM_CHARGE\t\t%d\t%s\t%9.5lf", i+1, charge[i].name, charge[i].charge);
	printf("\nPREP_FILE: %s\n", prepfile);

/*	read molecule information */
	default_cinfo(&cinfo);
	default_minfo(&minfo);

	atom = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
	if (atom == NULL) {
		fprintf(stdout, "memory allocation error for *atom\n");
		exit(1);
	}

	bond = (BOND *) malloc(sizeof(BOND) * cinfo.maxbond);
	if (bond == NULL) {
		fprintf(stdout, "memory allocation error for *bond\n");
		exit(1);
	}
	for (i = 0; i < cinfo.maxbond; ++i) {
		bond[i].jflag = -1; /* bond type has not been assigned */
	}

	overflow_flag =
		rac(inputfile, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);

	if (overflow_flag) {
		cinfo.maxatom = atomnum + 10;
		cinfo.maxbond = bondnum + 10;
		free(atom);
		free(bond);
		atom = (ATOM *) malloc(sizeof(ATOM) * cinfo.maxatom);
		if (atom == NULL) {
			fprintf(stdout, "memory allocation error for *atom\n");
			exit(1);
		}
		bond = (BOND *) malloc(sizeof(BOND) * cinfo.maxbond);
		if (bond == NULL) {
			fprintf(stdout, "memory allocation error for *bond\n");
			exit(1);
		}
		int i;
		for (i = 0; i < cinfo.maxbond; ++i) {
			bond[i].jflag = -1; /* bond type has not been assigned */
		}
		overflow_flag =
			rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
	}
	atomicnum(atomnum, atom);
	if(strcmp(atom_type_set, "default") != 0) {
		build_exe_path(command, "antechamber", sizeof command, 1);
		strcat(command, " -fi ac -fo ac -i ");
		strcat(command, inputfile);
		strcat(command, " -o RESIDUE_GEN.AC0 -at ");
		strcat(command, atom_type_set);
		printf("\nCommand: %s\n\n", command);
		system(command);
		fflush(stdout);
		rac("RESIDUE_GEN.AC0", &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
	}	

/*	allocate memory*/
	selectelement = (int *) malloc(sizeof(int) * atomnum);
	if (selectelement == NULL) {
       		fprintf(stdout, "memory allocation error for *selectelement\n");
       		exit(1);
	}
	selectelement2 = (int *) malloc(sizeof(int) * atomnum);
	if (selectelement2 == NULL) {
       		fprintf(stdout, "memory allocation error for *selectelement2\n");
       		exit(1);
	}
	flag = (int *) malloc(sizeof(int) * atomnum);
	if (flag == NULL) {
        	fprintf(stdout, "memory allocation error for *flag\n");
        	exit(1);
	}

/*	find atom ids for predefined charges*/
	for(i =0; i< chargenum; i++) { 
		charge[i].id = -1;
		index = 0;
		for(j =0; j< atomnum; j++) 
			if(strcmp(atom[j].name, charge[i].name) == 0) {
				charge[i].id = j;
				index = 1;
				break;
			}
		if(index == 0) {
			fprintf(stdout, "Error: the charge name, %s does not show up in the molecule\n", charge[i].name);
			exit(1);
		}
	}

/*	find atom ids for separted bonds*/
	for(i =0; i<sepbondnum; i++) {
		sepbond[i].resid = -1;
		sepbond[i].capid = -1;
		for(j =0; j< atomnum; j++) {
			if(strcmp(atom[j].name, sepbond[i].resat) == 0) {
				sepbond[i].resid = j; 
				if(strcmp(sepbond[i].atomtype, "DU") != 0) 
					strcpy(atom[j].ambername, sepbond[i].atomtype);
				continue;
			}
			if(strcmp(atom[j].name, sepbond[i].capat) == 0) { 
				sepbond[i].capid = j; 
/*				
 				if(strcmp(sepbond[i].atomtype_cap, "DU") != 0) 
					strcpy(atom[j].ambername, sepbond[i].atomtype_cap);
*/
				continue;
			}
			if(sepbond[i].resid != -1 && sepbond[i].capid != -1) break;
		}
		if(sepbond[i].resid == -1) {
			fprintf(stdout, "Error: the ATOM_NAME1 %s of %d does not show up in the molecule\n", sepbond[i].resat, i+1);
			exit(1);
		}
		if(sepbond[i].capid == -1) {
			fprintf(stdout, "Error: the ATOM_NAME2 %s of %d does not show up in the molecule\n", sepbond[i].capat, i+1);
			exit(1);
		}
	
		if(sepbond[i].type == 0) 
			printf("SEP_BOND  %5d NORMAL %5s %5d %5s CAP: %5s %5d %5s %5s\n", i+1, sepbond[i].resat, sepbond[i].resid, sepbond[i].atomtype, 
                                                                           sepbond[i].capat, sepbond[i].capid, atom[sepbond[i].capid].ambername, sepbond[i].atomtype_cap);
		if(sepbond[i].type == 1) {
			printf("SEP_BOND  %5d HEAD   %5s %5d %5s CAP: %5s %5d %5s %5s\n", i+1, sepbond[i].resat, sepbond[i].resid, sepbond[i].atomtype, 
                                                                           sepbond[i].capat, sepbond[i].capid, atom[sepbond[i].capid].ambername, sepbond[i].atomtype_cap);
			head_id = sepbond[i].resid;
			pre_head_id = sepbond[i].capid;
			strcpy(pre_head_atom_type, sepbond[i].atomtype_cap);
			if(strcmp(pre_head_atom_type, "DU") == 0)
				strcpy(pre_head_atom_type, atom[pre_head_id].ambername);
			num_sep_bond_type1 ++;
		}
		if(sepbond[i].type ==-1) {
			printf("SEP_BOND  %5d TAIL   %5s %5d %5s CAP: %5s %5d %5s %5s\n", i+1, sepbond[i].resat, sepbond[i].resid, sepbond[i].atomtype, 
                                                                           sepbond[i].capat, sepbond[i].capid, atom[sepbond[i].capid].ambername, sepbond[i].atomtype_cap);
			tail_id = sepbond[i].resid;
			post_tail_id = sepbond[i].capid;
			strcpy(post_tail_atom_type, sepbond[i].atomtype_cap);
			if(strcmp(post_tail_atom_type, "DU") == 0)
				strcpy(post_tail_atom_type, atom[post_tail_id].ambername);
			num_sep_bond_type2 ++;
		}
	}

	if(num_sep_bond_type1 > 1) {
		fprintf(stdout, "Error: only one or none SEP_BOND has a type of '1' allowed\n");
		exit(1);
	}
	if(num_sep_bond_type2 > 1) {
		fprintf(stdout, "Error: only one or none SEP_BOND has a type of '-1' allowed\n");
		exit(1);
	}

/* 	overwrite the head atom info obtained from SEPBOND */
        if(strcmp(head_atom_name, "UN") != 0) {
                for(i=0;i<atomnum;i++)
                        if(strcmp(atom[i].name, head_atom_name) == 0) {
                                head_id = i;
                                break;
                        }
		if(head_id != -1) strcpy(atom[head_id].ambername, head_atom_type);	
        }
        if(strcmp(tail_atom_name, "UN") != 0) {
                for(i=0;i<atomnum;i++)
                        if(strcmp(atom[i].name, tail_atom_name) == 0) {
                                tail_id = i;
                                break;
                        }
		if(tail_id != -1) strcpy(atom[tail_id].ambername, tail_atom_type);	
        }

/*	find cap atoms and residue atoms*/
	if(sepbondnum > 0) {
		for(i=0;i<sepbondnum;i++) {
        		printf("\n\nCap Atoms for #%d separating bond (%s - %s) with a type of %d", i+1, sepbond[i].capat, sepbond[i].resat, sepbond[i].type);
        		group_atom1(i);
        		for (j = 0; j<atomnum; j++)
                	if(flag[j] == 1)
                        	printf("\nATOM\t %d\t%s", j+1, atom[j].name);
		}
       		for (i = 0; i<atomnum; i++)
               		flag[i] = -1;
        	printf("\n\nResidue Atoms");
        	group_atom2();
        	residueatnum = 0;
        	for (i = 0; i<atomnum; i++)
                	if(flag[i] == 1) {
                        	residueatnum ++;
                        	printf("\nATOM\t %d\t%s", i+1, atom[i].name);
                	}
	}
	else {
       		for (i = 0; i<atomnum; i++)
               		flag[i] = 1;
	}

/*      find main chain atoms both head_id and tail_id are not -1, otherwise no main chain atoms are determined*/
        if(head_id != -1 && tail_id != -1) {
		printf("\n\nHEAD ATOM : %5s %5s %5s", atom[head_id].name, atom[head_id].ambername, pre_head_atom_type);
		printf("\n\nTAIL ATOM : %5s %5s %5s", atom[tail_id].name, atom[tail_id].ambername, post_tail_atom_type);
                selectnum = 0;
                for (i = 0; i<atomnum; i++) {
                        selectelement[i] = -1;
                        selectelement2[i] = -1;
		}
                findpath(atom, 0, head_id);
                printf("\n\nMain Chain Atoms");
                for (j = 0; j < selectnum2 ; j++)
                        printf("\nATOM\t%d\t%s", j+1, atom[selectelement2[j]].name);
        }

/*	generate charge input file for respgen */
        printf("\nGenerate charge input file for respgen ... \n");
        if ((fpout = fopen("RESIDUE_GEN_RESPGEN.DAT", "w")) == NULL) {
                fprintf(stdout, "Cannot open input file RESIDUE_GEN_RESPGEN.DAT, exit\n");
                exit(1);
        }
        fprintf(fpout, "//predefined charges in a format of (CHARGE partical_charge atom_ID atom_name)");
        for(i=0;i<chargenum;i++)
                fprintf(fpout, "\nCHARGE        %9.6lf\t%d\t%s", charge[i].charge, charge[i].id + 1, charge[i].name);

        fprintf(fpout, "\n//charge groups in a format of (GROUP num_atom net_charge)");
        fprintf(fpout, "\nGROUP\t%d\t%9.5lf\n", residueatnum, netcharge);       
        fprintf(fpout, "//atoms in the group in a format of (ATOM atom_ID atom_name)");
        for (i = 0; i<atomnum; i++)
                if(flag[i] == 1) 
                        fprintf(fpout, "\nATOM\t%d\t%s", i+1, atom[i].name);
        fclose(fpout);

/*	generate main chain file for prepgen */
        printf("\nBegin Real Work ... \n");
        printf("\nGenerate main chain atom file for prepgen ... \n");
        if ((fpout = fopen("RESIDUE_GEN_MAINCHAIN.DAT", "w")) == NULL) {
                fprintf(stdout, "Cannot open input file RESIDUE_GEN_MAINCHAIN.DAT, exit\n");
                exit(1);
        }
        if(head_id != -1 && pre_head_id != -1) 
		fprintf(fpout, "HEAD_NAME\t%s\n", atom[head_id].name);
        if(tail_id != -1 && post_tail_id != -1) 
        	fprintf(fpout, "TAIL_NAME\t%s\n", atom[tail_id].name);
        for (i = 0; i < selectnum2 ; i++) {
		if(head_id != -1 && pre_head_id != -1) {
/*
                	length = strlen(atom[head_id].name);
                	if(strncmp(atom[selectelement2[i]].name, atom[head_id].name, length) == 0) 
*/
			if(strcmp(atom[selectelement2[i]].name, atom[head_id].name) == 0) 
				continue;
		}
		if(tail_id != -1 && post_tail_id != -1) {
/*
                	length = strlen(atom[tail_id].name);
                	if(strncmp(atom[selectelement2[i]].name, atom[tail_id].name, length) == 0) 
*/
                	if(strcmp(atom[selectelement2[i]].name, atom[tail_id].name) == 0) 
				continue;
		}
                fprintf(fpout, "MAIN_CHAIN\t%s\n", atom[selectelement2[i]].name);
        }
        for (i = 0; i<atomnum; i++)
                if(flag[i] == 1)
                        continue;
                else
                        fprintf(fpout, "OMIT_NAME\t%s\n", atom[i].name);
	if(strcmp(pre_head_atom_type, "DU") != 0) {
		fprintf(fpout, "PRE_HEAD_TYPE\t%s\n", pre_head_atom_type);
		printf("PRE_HEAD_TYPE\t%s\n", pre_head_atom_type);
	}
	if(strcmp(post_tail_atom_type, "DU") != 0) {
		fprintf(fpout, "POST_TAIL_TYPE\t%s\n", post_tail_atom_type);
		printf("POST_TAIL_TYPE\t%s\n", post_tail_atom_type);
	}
        fprintf(fpout, "CHARGE\t%9.5lf\n", netcharge);
        fclose(fpout);
	fflush(stdout);

	if(confnum >= 1)
		proceed1();
	else
		proceed2();
	printf("\n");
	return (0);

}
