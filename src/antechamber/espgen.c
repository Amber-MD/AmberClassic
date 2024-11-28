/*
************************************************************************
*           All Copyright Reserved!                                    *
*                                                                      *
*  Prog:    espgen                                                     *
*  Version: version 1.1                                                *
*  Author:  Junmei Wang                                                *
*  GAMESS support added by Robin Betz, 2019                            *
*                                                                      *
*  Department of Pharmaceutical Chemistry                              *
*  School of Pharmacy                                                  *
*  University of California                                            *
*  San Francisco   CA 94143                                            *
*  Octomber, 2001                                                      *
************************************************************************
*/
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include "define.h"
# include "atom.h"
# include "eprintf.h"
# define MAXESP 99999
# define MAX_ATOM_CENTER 1000

typedef struct {
        char atomname[10];
        char atomtype[10];
        int  atomicnum;
        int  iequ;
        double cord_x;
        double cord_y;
        double cord_z;
        double charge_mul;
        double charge_esp;
        double radius_pcm;
        double radius_esp;
} ATOMINFO;

int i, j, k;
int fail = 0;
int method = 0;                 /* output dipole and quadrupole? */
int maxesp = 0;
int format = 1;                 /* 1: log files of G98, G03, g16, 2: G09 ESP file, specia output triggered by iop(6/50=1) */
int formalchrg = 0;
int multiplicity=1;
char line[MAXCHAR];
char ifilename[MAXCHAR];
char ofilename[MAXCHAR];

int    iresname = 0;
char   resname[MAXCHAR]="MOL";
char   atomname[10];
char   atomname2[MAXCHAR];
char   atomtype[10];

ATOMINFO at[MAX_ATOM_CENTER];

ESP *esp;
DM dipole;
QM quadrupole;
FILE *fpin, *fpout;

int Found_Stationary = 0;
int esp_index = 0;
int espvalue_index = 0;
int g09_index = 0;
int opt_index = 0;
int npoint = 0;
int natom = 0;
int nset = 0;

char tmpchar1[MAXCHAR], tmpchar2[MAXCHAR], tmpchar3[MAXCHAR];
char tmpchar4[MAXCHAR], tmpchar5[MAXCHAR], tmpchar6[MAXCHAR];
char tmpchar7[MAXCHAR], tmpchar8[MAXCHAR], tmpchar9[MAXCHAR];
int tmpint0, tmpint1, tmpint2, tmpint3, tmpint4;
double tmpfloat1, tmpfloat2, tmpfloat3, tmpfloat4;

int pflag = 0;
int remark= 0;
double netcharge = 0;

void read_respin(char *filename) {
FILE *fp;
char tmpc[MAXCHAR];
int count;
	if ((fp = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "\n Cannot open the %s file, exit", filename);
		exit (1);
	}
	for (;;) {
		if (fgets(line, MAXCHAR, fp) == NULL)
			break;
		strcpy(tmpc, "");
		tmpc[0]='\0';
		sscanf(line, "%s", tmpc);
		if(strcmp(tmpc, "&end") == 0) {
			fgets(line, MAXCHAR, fp);
			fgets(line, MAXCHAR, fp);
			fgets(line, MAXCHAR, fp);
			count = 0;
			while(count < natom) {
				fgets(line, MAXCHAR, fp);
				sscanf(line, "%d%d", &tmpint1, &tmpint2);
				if(tmpint2 > 0)  
					at[count].iequ = tmpint2;
				count++;	
			}			
		}
	}
	fclose(fp);
}
void assign_atomtype() {
char commandstr[MAXCHAR];
int rflag = -1;
int count ;
int len;
FILE *fp;
	strcpy(commandstr, "antechamber -fi ");
	if(format == 1) 
		strcat(commandstr, " gout -fo ac -i ");
	else
		strcat(commandstr, " gesp -fo ac -i ");
	strcat(commandstr, ifilename); 	
	strcat(commandstr, " -o ESPGEN.TMP -pf y -dr no");
//	fprintf(stdout, "Running command: %s ...\n", commandstr);
	rflag=system(commandstr);
	if(rflag == -1) {
		fprintf(stderr, "Error happens when running %s\n", commandstr);
		exit(1);
	}
	if ((fp = fopen("ESPGEN.TMP", "r")) == NULL) {
		fprintf(stderr, "\n Cannot open the ESPGEN.TMP file, exit");
		exit (1);
	}
	count = 0;
	for (;;) {
		if (fgets(line, MAXCHAR, fp) == NULL)
			break;
		if (strncmp(line, "ATOM", 4) == 0) {
			if(iresname == 0 && count == 0) {
				resname[0] = line[17];
				resname[1] = line[18];
				resname[2] = line[19];
				if(resname[1] == ' ')
					resname[1]='\0';
				else if(resname[2] == ' ')
					resname[2]='\0';
				else
                        		resname[3] = '\0';
			}
                	atomname[0] = line[12];
                	atomname[1] = line[13];
                	atomname[2] = line[14];
                	atomname[3] = line[15];
                	atomname[4] = '\0';
                	sscanf(atomname, "%s", atomname2);
                	if(atomname2[0] >= '0' && atomname2[0] <= '9') {
				len=strlen(atomname2);
                        	for(i=1;i<len;i++)
                                	atomname[i-1] = atomname2[i];
                                atomname[len-1] = atomname2[0];
                                atomname[len] = '\0';
			}
                	else
                        	strcpy(atomname, atomname2);
			sscanf(&line[65], "%s", atomtype);
			strcpy(at[count].atomname, atomname);
			strcpy(at[count].atomtype, atomtype);
			count++;
			if(count >= MAX_ATOM_CENTER) {
				fprintf(stdout,"Error, the number of atomic center exceeds MAX_ATOM_CENTER defined in espgen.c, extend MAX_ATOM_CENTER and recompile the program, tmpint1, MAX_ATOM_CENTER\n");
				exit(1);
			}
		}	
	}
	natom=count;
	fclose(fp);
/* now we generate esp input files*/
        strcpy(commandstr, "respgen -i ESPGEN.TMP -o RESPIN1 -l 20 -f resp1");
        fprintf(stdout, "Running command: %s ...\n", commandstr);
        system(commandstr);
        strcpy(commandstr, "respgen -i ESPGEN.TMP -o RESPIN2 -l 20 -f resp2");
        fprintf(stdout, "Running command: %s ...\n", commandstr);
        system(commandstr);
        read_respin("RESPIN1");
        read_respin("RESPIN2");
}

int judge_atomicnum(char *elem) {
int atomic_num = -1;
switch (elem[0]) {
	case 'A':
		if (elem[1] == 'c' || elem[1] == 'C') 
			atomic_num = 89;
		else if (elem[1] == 'g' || elem[1] == 'G')
			atomic_num = 47;
		else if (elem[1] == 'l' || elem[1] == 'L')
			atomic_num = 13;
		else if (elem[1] == 'm' || elem[1] == 'M')
			atomic_num = 95;
		else if (elem[1] == 'r' || elem[1] == 'R')
			atomic_num = 18;
		else if (elem[1] == 's' || elem[1] == 'S')
			atomic_num = 33;
		else if (elem[1] == 't' || elem[1] == 'T')
			atomic_num = 85;
		else if (elem[1] == 'u' || elem[1] == 'U')
			atomic_num = 79;
		break;
	case 'B':
		if (elem[1] == 'r' || elem[1] == 'R')
			atomic_num = 35;
		else if (elem[1] == 'a' || elem[1] == 'A')
			atomic_num = 56;
		else if (elem[1] == 'e' || elem[1] == 'E')
			atomic_num = 4;
		else if (elem[1] == 'h' || elem[1] == 'H')
			atomic_num = 107;
		else if (elem[1] == 'i' || elem[1] == 'I')
			atomic_num = 83;
		else if (elem[1] == 'k' || elem[1] == 'K')
			atomic_num = 97;
		else
			atomic_num = 5;
		break;
	case 'C':
		if (elem[1] == 'l' || elem[1] == 'L')
			atomic_num = 17;
		else if (elem[1] == 'a')
			atomic_num = 20;
		else if (elem[1] == 'd')
			atomic_num = 48;
		else if (elem[1] == 'e')
			atomic_num = 58;
		else if (elem[1] == 'f')
			atomic_num = 98;
		else if (elem[1] == 'm')
			atomic_num = 96;
		else if (elem[1] == 'o')
			atomic_num = 27;
		else if (elem[1] == 'r')
			atomic_num = 24;
		else if (elem[1] == 's')
			atomic_num = 55;
		else if (elem[1] == 'u')
			atomic_num = 29;
		else
			atomic_num = 6;
		break;
	case 'D':
		if (elem[1] == 'b')
			atomic_num = 105;
		else if (elem[1] == 's' || elem[1] == 'S')
			atomic_num = 110;
		else if (elem[1] == 'y' || elem[1] == 'Y')
			atomic_num = 66;
		else
			atomic_num = 1;
		break;
	case 'E':
		if (elem[1] == 'P')
			atomic_num = 0;
		else if (elem[1] == 'r' || elem[1] == 'R')
			atomic_num = 68;
		else if (elem[1] == 's' || elem[1] == 'S')
			atomic_num = 99;
		else if (elem[1] == 'u' || elem[1] == 'U')
			atomic_num = 63;
		break;
	case 'F':
		if (elem[1] == 'e' || elem[1] == 'E')
			atomic_num = 26;
		else if (elem[1] == 'm' || elem[1] == 'M')
			atomic_num = 100;
		else if (elem[1] == 'r' || elem[1] == 'R')
			atomic_num = 87;
		else
			atomic_num = 9;
		break;
	case 'G':
		if (elem[1] == 'a' || elem[1] == 'A')
			atomic_num = 31;
		else if (elem[1] == 'd' || elem[1] == 'D')
			atomic_num = 64;
		else if (elem[1] == 'e' || elem[1] == 'E')
			atomic_num = 32;
		break;
	case 'H':
		if (elem[1] == 'e')
			atomic_num = 2;
		else if (elem[1] == 'f')
			atomic_num = 72;
		else if (elem[1] == 'g')
			atomic_num = 80;
		else if (elem[1] == 'o')
			atomic_num = 67;
		else if (elem[1] == 's')
			atomic_num = 108;
		else
			atomic_num = 1;
		break;
	case 'I':
		if (elem[1] == 'n' || elem[1] == 'N')
			atomic_num = 49;
		if (elem[1] == 'r' || elem[1] == 'R')
			atomic_num = 77;
		else 
			atomic_num = 53;
		break;
	case 'K':
		if (elem[1] == 'r' || elem[1] == 'R')
			atomic_num = 36;
		else
			atomic_num = 19;
		break;
	case 'l':
		if (elem[1] == 'p')
			atomic_num = 0;
		break;
	case 'L':
		if (elem[1] == 'i' || elem[1] == 'I')
			atomic_num = 3;
		else if (elem[1] == 'a' || elem[1] == 'A')
			atomic_num = 57;
		else if (elem[1] == 'r' || elem[1] == 'R')
			atomic_num = 103;
		else if (elem[1] == 'u' || elem[1] == 'U')
			atomic_num = 71;
		else if (elem[1] == 'P')
			atomic_num = 0;
		break;
	case 'M':
		if (elem[1] == 'n' || elem[1] == 'N')
			atomic_num = 25;
		else if (elem[1] == 'g' || elem[1] == 'G')
			atomic_num = 12;
		else if (elem[1] == 'd' || elem[1] == 'D')
			atomic_num = 101;
		else if (elem[1] == 'o' || elem[1] == 'O')
			atomic_num = 42;
		else if (elem[1] == 't' || elem[1] == 'T')
			atomic_num = 109;
		break;
	case 'N':
		if (elem[1] == 'i')
			atomic_num = 28;
		else if (elem[1] == 'a')
			atomic_num = 11;
		else if (elem[1] == 'b')
			atomic_num = 41;
		else if (elem[1] == 'd')
			atomic_num = 60;
		else if (elem[1] == 'e')
			atomic_num = 10;
		else if (elem[1] == 'o')
			atomic_num = 102;
		else if (elem[1] == 'p')
			atomic_num = 93;
		else
			atomic_num = 7;
		break;
	case 'O':
		if (elem[1] == 's')
			atomic_num = 76;
		else
			atomic_num = 8;
		break;
	case 'P':
		if (elem[1] == 'd')
			atomic_num = 46;
		else if (elem[1] == 't')
			atomic_num = 78;
		else if (elem[1] == 'b')
			atomic_num = 82;
		else if (elem[1] == 'a')
			atomic_num = 91;
		else if (elem[1] == 'm')
			atomic_num = 61;
		else if (elem[1] == 'o')
			atomic_num = 84;
		else if (elem[1] == 'r')
			atomic_num = 59;
		else if (elem[1] == 'u')
			atomic_num = 94;
		else
			atomic_num = 15;
		break;
	case 'R':
		if (elem[1] == 'u' || elem[1] == 'U')
			atomic_num = 44;
		else if (elem[1] == 'h' || elem[1] == 'H')
			atomic_num = 45;
		else if (elem[1] == 'a' || elem[1] == 'A')
			atomic_num = 88;
		else if (elem[1] == 'b' || elem[1] == 'B')
			atomic_num = 37;
		else if (elem[1] == 'e' || elem[1] == 'E')
			atomic_num = 75;
		else if (elem[1] == 'f' || elem[1] == 'F')
			atomic_num = 104;
		else if (elem[1] == 'n' || elem[1] == 'N')
			atomic_num = 86;
		break;
	case 'S':
		if (elem[1] == 'i' || elem[1] == 'I')
			atomic_num = 14;
		else if (elem[1] == 'c')
			atomic_num = 21;
		else if (elem[1] == 'e')
			atomic_num = 34;
		else if (elem[1] == 'r')
			atomic_num = 38;
		else if (elem[1] == 'b')
			atomic_num = 51;
		else if (elem[1] == 'g')
			atomic_num = 106;
		else if (elem[1] == 'm')
			atomic_num = 62;
		else if (elem[1] == 'n')
			atomic_num = 50;
		else 
			atomic_num = 16;
		break;
	case 'T':
		if (elem[1] == 'i' || elem[1] == 'I')
			atomic_num = 22;
		else if (elem[1] == 'l' || elem[1] == 'L')
			atomic_num = 81;
		else if (elem[1] == 'a' || elem[1] == 'A')
			atomic_num = 73;
		else if (elem[1] == 'b' || elem[1] == 'B')
			atomic_num = 65;
		else if (elem[1] == 'c' || elem[1] == 'C')
			atomic_num = 43;
		else if (elem[1] == 'e' || elem[1] == 'E')
			atomic_num = 52;
		else if (elem[1] == 'h' || elem[1] == 'H') 
			atomic_num = 90;
		else if (elem[1] == 'm' || elem[1] == 'M') 
			atomic_num = 69;
		else 
			atomic_num = 1;
		break;
	case 'U':
		atomic_num = 92;
		break;
	case 'V':
		atomic_num = 23;
		break;
	case 'W':
		atomic_num = 74;
		break;
	case 'X':
		if (elem[1] == 'e' || elem[1] == 'E') 
			atomic_num = 54;
		break;
	case 'Y':
		if (elem[1] == 'b' || elem[1] == 'B') 
			atomic_num = 70;
		else 
			atomic_num = 39;
		break;
	case 'Z':
		if (elem[1] == 'n' || elem[1] == 'N') 
			atomic_num = 30;
		else if (elem[1] == 'r' || elem[1] == 'R') 
			atomic_num = 40;
		break;
	default:
		printf("\n Unrecognized atomic name %5s, exit", elem);
}
return atomic_num;
}
void rgamessdat(void)
{
    int atoms = 0, points = 0;
    int i;

    rewind(fpin);

    // Dipole and quadrupole currently unimplemented in this parser
    if (method == 1)
        eprintf("Dipole and quadrupole support unimplemented in GAMESS parser.");

    // Find line with number of atomic centers
    for (;;) {
        if (!fgets(line, MAXCHAR, fpin))
            eprintf("Premature end of file before electrostatic potential");

        // ELECTROSTATIC POTENTIAL COMPUTED FOR <N> ATOMS, TOTAL CHARGE= <N>
        if (strncmp(line, " ELECTROSTATIC POTENTIAL COMPUTED FOR", 37) == 0) {
            sscanf(line, "%*s%*s%*s%*s%d", &atoms);
            if (atoms > MAX_ATOM_CENTER) {
                eprintf("The number of atomic centers (%d) exceeded "
                        "MAX_ATOM_CENTER.\n Increase MAX_ATOM_CENTER in "
                        "espgen.c and rebuild.", atoms);
            }
            break;
        }
    }

    // Read in each atomic center
    for (i = 0; i < atoms; i++) {
        if (!fgets(line, MAXCHAR, fpin))
            eprintf("Premature end of file at atom %d", i);

        // Line has index, atomicnum, x, y, z coordinates
        sscanf(line, "%*d%*f%lf%lf%lf", &at[i].cord_x, &at[i].cord_y, &at[i].cord_z); 
    }

    // Throw out "ELECTROSTATIC POTENTIAL, IPT,X,Y,Z,ELPOTT" line
    if (!fgets(line, MAXCHAR, fpin))
        eprintf("Premature end of file before grid points definition");

    // Read in number of grid points
    if (!fgets(line, MAXCHAR, fpin))
        eprintf("Premature end of file before grid points number");

    if (strncmp(line, " TOTAL NUMBER OF GRID POINTS NPT=", 33))
        eprintf("Grid points line not found as expected");

    sscanf(line, "%*s%*s%*s%*s%*s%*s%d", &points);

    if (points >= maxesp) {
        printf("\nInfo: The number of ESP exceeded MAXESP of %d; "
               "automatically increasing to %d.\n", maxesp, maxesp + MAXESP);
        maxesp += MAXESP;
        esp = (ESP *) erealloc(esp, maxesp * sizeof(ESP));
    }

    // Read in each grid point
    for (i = 0; i < points; i++) {
        if (!fgets(line, MAXCHAR, fpin))
            eprintf("Premature end of file reading ESP point %d", i);

        sscanf(line, "%*d%lf%lf%lf%lf", &esp[i].x, &esp[i].y, &esp[i].z,
               &esp[i].esp);
    }

    // Now write out to the output file
    // First line is number of atoms, number of points, and 0
    fprintf(fpout, "%5d%5d%5d\n", atoms, points, 0);

    // Then all of the atomic coordinates
    // These are already in units of Bohr
    for (i = 0; i < atoms; i++) {
        fprintf(fpout, "%32.7E%16.7E%16.7E\n", at[i].cord_x, at[i].cord_y, at[i].cord_z); 
    }

    // Finally, all of the ESP points, which are already in units of Bohr
    for (i = 0; i < points; i++) {
        fprintf(fpout, "%16.7E%16.7E%16.7E%16.7E\n", esp[i].esp,
                esp[i].x, esp[i].y, esp[i].z);
    }

    printf("Success!");
}

void rglog(void)
{
    int i, j;
    int ipol= 0;

    tmpint1 = 0;
    tmpint2 = 0;
    tmpint3 = 0;
    tmpint4 = 0;
    rewind(fpin);
    for (;;) {
        if (fgets(line, MAXCHAR, fpin) == NULL)
            break;
        strcpy(tmpchar1, "");
        strcpy(tmpchar2, "");
        strcpy(tmpchar3, "");
        strcpy(tmpchar4, "");
        strcpy(tmpchar5, "");
        strcpy(tmpchar6, "");
        tmpchar1[0]='\0';
        tmpchar2[0]='\0';
        tmpchar3[0]='\0';
        tmpchar4[0]='\0';
        tmpchar5[0]='\0';
        tmpchar6[0]='\0';
        sscanf(&line[0], "%s%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3,
                           tmpchar4, tmpchar5, tmpchar6);
        if (strcmp("--", tmpchar1) == 0 
            && strcmp("Stationary", tmpchar2) == 0
            && strcmp("point", tmpchar3) == 0 
            && strcmp("found.", tmpchar4) == 0)
            Found_Stationary = 1;
        if (opt_index == 1 && nset == 2 && Found_Stationary == 0)
            continue;
        if (pflag == 1) {
                if (strcmp("Charge",tmpchar1) == 0 && strcmp("=", tmpchar2) == 0 &&
                    strcmp("Multiplicity",tmpchar4) == 0 && strcmp("=", tmpchar5) == 0 ) {
                                sscanf(&tmpchar3[0], "%d", &formalchrg);
                                sscanf(&tmpchar6[0], "%d", &multiplicity);
                }
                if (strcmp("Input", tmpchar1) == 0 && strcmp("orientation:", tmpchar2) == 0) {
                        tmpint4 = 0;
                        if (fgets(line, MAXCHAR, fpin) == NULL) break;
                        if (fgets(line, MAXCHAR, fpin) == NULL) break;
                        if (fgets(line, MAXCHAR, fpin) == NULL) break;
                        if (fgets(line, MAXCHAR, fpin) == NULL) break;
                        while(1) {
                                fgets(line, MAXCHAR, fpin);
                                if(strncmp(line, " ---------", 10) == 0) break;
                                sscanf(&line[0], "%s%d", tmpchar1, &at[tmpint4].atomicnum);
                                tmpint4++;
                                if(tmpint4 > MAX_ATOM_CENTER) {
                			eprintf("The number of atomic centers (%d) exceeded MAX_ATOM_CENTER.\n"
                        			"Increase MAX_ATOM_CENTER in espgen.c and rebuild.", tmpint1);
                                        exit(1);
                                }
                        }
                }
                if (strcmp(tmpchar1, "Mulliken") == 0 && strcmp(tmpchar2, "atomic") == 0 && strcmp(tmpchar3, "charges:") == 0) {
                        fgets(line, MAXCHAR, fpin);
                        while(1) {
                                fgets(line, MAXCHAR, fpin);
                                strcpy(tmpchar1, "");
                                strcpy(tmpchar2, "");
                                strcpy(tmpchar3, "");
                                tmpchar1[0]='\0';
                                tmpchar2[0]='\0';
                                tmpchar3[0]='\0';
                                sscanf(&line[0], "%s%s%s", tmpchar1, tmpchar2, tmpchar3);
                                if (strcmp(tmpchar1, "Sum") == 0  && strcmp(tmpchar2, "of") ==0 &&
                                    strcmp(tmpchar3, "atomic")==0 && strcmp(tmpchar4, "Mulliken") == 0) {
                                        sscanf(&line[33], "%lf", &netcharge);
                                        break;
                                }
                                tmpint0=atoi(tmpchar1);
                                at[tmpint0-1].charge_mul = atof(tmpchar3);
                                if(tmpint0 == natom) break;
                        }
                }
                if (strcmp(tmpchar1, "Mulliken") == 0 && strcmp(tmpchar2, "charges:") == 0) {
                        fgets(line, MAXCHAR, fpin);
                        while(1) {
                                fgets(line, MAXCHAR, fpin);
                                strcpy(tmpchar1, "");
                                strcpy(tmpchar2, "");
                                strcpy(tmpchar3, "");
                                tmpchar1[0]='\0';
                                tmpchar2[0]='\0';
                                tmpchar3[0]='\0';
                                sscanf(&line[0], "%s%s%s", tmpchar1, tmpchar2, tmpchar3);
                                if (strcmp(tmpchar1, "Sum") == 0 && strcmp(tmpchar2, "of") ==0 &&
                                    strcmp(tmpchar3, "Mulliken") == 0) {
                                        sscanf(&line[26], "%lf", &netcharge);
                                        break;
                                }
                                tmpint0=atoi(tmpchar1);
                                at[tmpint0-1].charge_mul = atof(tmpchar3);
                                if(tmpint0 == natom) break;
                        }
                }

                if (strcmp(tmpchar1, "ESP") == 0 && strcmp(tmpchar2, "charges:") == 0) {
                        fgets(line, MAXCHAR, fpin);
                        while(1) {
                                fgets(line, MAXCHAR, fpin);
                                strcpy(tmpchar1, "");
                                strcpy(tmpchar2, "");
                                strcpy(tmpchar3, "");
                                tmpchar1[0]='\0';
                                tmpchar2[0]='\0';
                                tmpchar3[0]='\0';
                                sscanf(&line[0], "%s%s%s", tmpchar1, tmpchar2, tmpchar3);
                                if (strcmp(tmpchar1, "Sum") == 0 && strcmp(tmpchar2, "of") ==0 &&
                                    strcmp(tmpchar3, "ESP") == 0) {
                                        sscanf(&line[21], "%lf", &netcharge);
                                        break;
                                }
                                tmpint0=atoi(tmpchar1);
                                at[tmpint0-1].charge_esp = atof(tmpchar3);
                                if(tmpint0 == natom) break;
                        }
                }

                if (strcmp(tmpchar1, "Atom") == 0 && strcmp(tmpchar2, "Element") == 0 && strcmp(tmpchar3, "Radius") == 0) {
                        while(1) {
                                fgets(line, MAXCHAR, fpin);

                                if (strcmp(tmpchar1, "Generate") == 0 && strcmp(tmpchar2, "VDW") ==0 && strcmp(tmpchar3, "Surfaces:") == 0) {
                                        pflag = -1;
                                        break;
                                }
                                sscanf(&line[0], "%s%s%s", tmpchar1, tmpchar2, tmpchar3);
                                tmpint0=atoi(tmpchar1);
                                at[tmpint0-1].radius_esp = atof(tmpchar3);
                                if(tmpint0 == natom) break;
                        }
                }
        }
        if (strcmp("Atomic", tmpchar1) == 0 && strcmp("Center", tmpchar2) == 0) {
	    sscanf(&line[32], "%lf%lf%lf", &at[tmpint1].cord_x, &at[tmpint1].cord_y, 
		&at[tmpint1].cord_z);
            tmpint1++;
            if (tmpint1 > MAX_ATOM_CENTER) {
                eprintf("The number of atomic centers (%d) exceeded MAX_ATOM_CENTER.\n"
                        "Increase MAX_ATOM_CENTER in espgen.c and rebuild.", tmpint1);
            }
        }
        if (strcmp("Polarizable", tmpchar1) == 0 && strcmp("Continuum", tmpchar2) == 0 &&
            strcmp("Model", tmpchar3) == 0 && strcmp("(PCM)", tmpchar4) == 0 ) {
                ipol = 1;
                continue;
        }
        if(ipol == 1) {
                if (strcmp("ISph", tmpchar1) == 0 && strcmp("on", tmpchar2) == 0 &&
                    strcmp("Nord", tmpchar3) == 0 && strcmp("Re0", tmpchar4) == 0 ) {
                        while (1) {
                                fgets(line, MAXCHAR, fpin);
                                if (strncmp(line, " ----------", 11) == 0) break;
                                sscanf(&line[0], "%s%s%s%s", tmpchar1, tmpchar2, tmpchar3, tmpchar4);
                                tmpint0=atoi(tmpchar1);
                                at[tmpint0-1].radius_pcm = atof(tmpchar4);
                                if(tmpint0 == natom) break;
                        }
                }
        }

        if (strcmp("ESP", tmpchar1) == 0 && strcmp("Fit", tmpchar2) == 0
            && strcmp("Center", tmpchar3) == 0) {
            sscanf(&line[32], "%lf%lf%lf", &esp[tmpint2].x, &esp[tmpint2].y,
                   &esp[tmpint2].z);
            esp[tmpint2].x /= Bohr;
            esp[tmpint2].y /= Bohr;
            esp[tmpint2].z /= Bohr;
            tmpint2++;
            if (tmpint2 >= maxesp) {
                maxesp += MAXESP;
                printf("\nInfo: The number of ESP exceeded MAXESP; automatically increasing to %d.\n",
                       maxesp);
                esp = (ESP *) erealloc(esp, maxesp * sizeof(ESP));
            }
        }
        if (strncmp("Fit", &line[6], 3) == 0) {
            sscanf(&line[12], "%lf", &esp[tmpint3++].esp);
        }
    }
    if (tmpint3 < 10)           // no ESP at all
        fail = 1;
    if (fail == 0) {
        if(pflag == 1) {
		if(remark == 1) {
                	fprintf(fpout, "# Atomic Unit\n");
                	fprintf(fpout, "# Line 1: #Center, #ESP Sites, Net_Charge, Residue_Name, #Center, #ESP Sites, Multiplicity\n");
                	fprintf(fpout, "# Section 2: Atom X, Y, Z, Atomic_Number, Atom_Type, Atom_Name, Equ_Index, PCM_Radius, ESP_Radius, Mul_Charge, ESP_Charge\n");
		}
                fprintf(fpout, "%5d%5d%5d %-5s %6d %6d %6d\n", tmpint1, tmpint2, formalchrg, resname,tmpint1, tmpint2,multiplicity);
	}
	else
        	fprintf(fpout, "%5d%5d%5d\n", tmpint1, tmpint2, 0);
        for (i = 0; i < tmpint1; i++) {
                if(pflag == 1)
                        fprintf(fpout, "%32.7E%16.7E%16.7E%4d %-3s %-5s %4d %8.4lf%8.4lf%12.6lf%12.6lf\n",
                                at[i].cord_x / Bohr, at[i].cord_y / Bohr, at[i].cord_z / Bohr,
                                at[i].atomicnum,  at[i].atomtype, at[i].atomname, at[i].iequ,
                                at[i].radius_pcm, at[i].radius_esp,
                                at[i].charge_mul, at[i].charge_esp);
                else
                        fprintf(fpout, "%32.7E%16.7E%16.7E\n", at[i].cord_x / Bohr, at[i].cord_y / Bohr, at[i].cord_z / Bohr);
	 }

        for (j = 0; j < tmpint2; j++)
            fprintf(fpout, "%16.7E%16.7E%16.7E%16.7E\n", esp[j].esp, esp[j].x, esp[j].y,
                    esp[j].z);

    }
    if (fail == 0 && method == 1) {
        rewind(fpin);
        Found_Stationary = 0;
        for (;;) {
            if (fgets(line, MAXCHAR, fpin) == NULL)
                break;
            sscanf(&line[0], "%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3, tmpchar4,
                   tmpchar5);
            if (strcmp("--", tmpchar1) == 0 && strcmp("Stationary", tmpchar2) == 0
                && strcmp("point", tmpchar3) == 0 && strcmp("found.", tmpchar4) == 0) {
                Found_Stationary = 1;
                continue;
            }
            if (opt_index == 0 || Found_Stationary == 1) {
                if (strcmp("Dipole", tmpchar1) == 0 && strcmp("moment", tmpchar2) == 0) {
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    sscanf(line, "%s%lf%s%lf%s%lf", tmpchar1, &dipole.x, tmpchar2,
                           &dipole.y, tmpchar3, &dipole.z);

                    continue;
                }
                if (strcmp("Traceless", tmpchar1) == 0
                    && strcmp("Quadrupole", tmpchar2) == 0
                    && strcmp("moment", tmpchar3) == 0) {
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    sscanf(line, "%s%lf%s%lf%s%lf", tmpchar1, &quadrupole.xx, tmpchar2,
                           &quadrupole.yy, tmpchar3, &quadrupole.zz);
                    if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
                    sscanf(line, "%s%lf%s%lf%s%lf", tmpchar1, &quadrupole.xy, tmpchar2,
                           &quadrupole.xz, tmpchar3, &quadrupole.yz);
                    break;
                }
            }
        }
        fprintf(fpout, "%16.7E%16.7E%16.7E\n", dipole.x, dipole.y, dipole.z);
        fprintf(fpout, "%16.7E%16.7E%16.7E\n", quadrupole.xx, quadrupole.yy,
                quadrupole.zz);
        fprintf(fpout, "%16.7E%16.7E%16.7E\n", quadrupole.xy, quadrupole.xz,
                quadrupole.yz);
    }

}

double translate(char *str)
{
    int i;
    int pos = 0;
    int flag = 0;
    char tmpc1[MAXCHAR];
    char tmpc2[MAXCHAR];
    double f, e;

    for (i = 0; i < strlen(str); i++) {
        if (str[i] == 'D' || str[i] == 'd' || str[i] == 'E' || str[i] == 'e') {
            pos = i;
            flag = 1;
            continue;
        }
        if (flag == 0)
            tmpc1[i] = str[i];
        if (flag == 1)
            tmpc2[i - pos - 1] = str[i];

    }

    f = atof(tmpc1);
    e = atof(tmpc2);
    return f * pow(10.0, e);
}

void rgesp()
{
    int rflag = 0;
    int count = 0;
    for (;;) {
        if (fgets(line, MAXCHAR, fpin) == NULL)
            break;
        sscanf(line, "%s%s%s%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3, tmpchar4, tmpchar5,
               tmpchar6, tmpchar7, tmpchar8);
        if (rflag == 3 && strcmp(tmpchar1, "ESP") == 0 && strcmp(tmpchar2, "VALUES") == 0) {
            rflag = 4;
            sscanf(line, "%s%s%s%s%s%s%s%s%d", tmpchar1, tmpchar2, tmpchar3, tmpchar4,
                   tmpchar5, tmpchar6, tmpchar7, tmpchar8, &npoint);
            if (npoint >= maxesp) {
                maxesp = npoint;
                printf("\nInfo: The number of ESP exceeded MAXESP; automatically increasing to %d.",
                       maxesp);
                esp = (ESP *) erealloc(esp, maxesp * sizeof(ESP));
            }
            continue;
        }
        if (rflag == 2 && strcmp(tmpchar1, "TRACELESS") == 0
            && strcmp(tmpchar2, "QUADRUPOLE") == 0) {
            rflag = 3;
            continue;
        }
        if (rflag == 1 && strcmp(tmpchar1, "DIPOLE") == 0) {
            rflag = 2;
            continue;
        }
        if (rflag == 0 && strncmp(line, " ATOMIC COORDINATES", 19) == 0) {
            rflag = 1;
            continue;
        }
        if (rflag == 1) {
            at[natom].cord_x = translate(tmpchar2);
            at[natom].cord_y = translate(tmpchar3);
            at[natom].cord_z = translate(tmpchar4);
            at[natom].charge_esp  = translate(tmpchar5);
            natom++;
            if (natom >= MAX_ATOM_CENTER) {
                eprintf("The number of atomic centers (%d) exceeded MAX_ATOM_CENTER.\n"
                        "Increase MAX_ATOM_CENTER in espgen.c and rebuild.", natom);
            }
        }
        if (rflag == 2) {
            dipole.x = translate(tmpchar2);
            dipole.y = translate(tmpchar4);
            dipole.z = translate(tmpchar6);
            dipole.total = translate(tmpchar8);
        }
        if (rflag == 3) {
            quadrupole.xx = translate(tmpchar2);
            quadrupole.yy = translate(tmpchar4);
            quadrupole.zz = translate(tmpchar6);

            fgets(line, MAXCHAR, fpin);
            sscanf(line, "%s%s%s%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3, tmpchar4,
                   tmpchar5, tmpchar6, tmpchar7, tmpchar8);
            quadrupole.xy = translate(tmpchar2);
            quadrupole.xz = translate(tmpchar4);
            quadrupole.yz = translate(tmpchar6);
        }
        if (rflag == 4) {
            if (count < npoint) {
                esp[count].esp = translate(tmpchar1);
                esp[count].x = translate(tmpchar2);
                esp[count].y = translate(tmpchar3);
                esp[count].z = translate(tmpchar4);
                count++;
            } else
                rflag = -1;
        }

    }
    if(pflag == 1) {
	 if(remark == 1) {
		fprintf(fpout, "# Basic Info: #Center, #ESP Sites, Net_Charge, Multiplicity\n");
        	fprintf(fpout, "%5d%5d%5d%5d\n", natom, npoint, formalchrg, multiplicity);
        	fprintf(fpout, "# ATOM Field: X, Y, Z, Atomic_Number, Atom_Type, Equ_Index, ESP_Charge\n");
    	}
	else 
        	fprintf(fpout, "%5d%5d%5d%5d\n", natom, npoint, formalchrg, multiplicity);
    }
    
    if(pflag == 0) 
        fprintf(fpout, "%5d%5d%5d\n", natom, npoint, 0);
    for (i = 0; i < natom; i++)
        if(pflag == 1)
                fprintf(fpout, "%32.7E%16.7E%16.7E%4d %4s %4d %12.6lf\n", at[i].cord_x, at[i].cord_y, at[i].cord_z,
                        at[i].atomicnum, at[i].atomtype, at[i].iequ, at[i].charge_esp);
        else
                fprintf(fpout, "%32.7E%16.7E%16.7E\n", at[i].cord_x, at[i].cord_y, at[i].cord_z);

    if(pflag == 1 && remark == 1) fprintf(fpout, "# ESP Field: ESP, X, Y, Z\n");
    for (j = 0; j < npoint; j++)
      		fprintf(fpout, "%16.7E%16.7E%16.7E%16.7E\n", esp[j].esp, esp[j].x, esp[j].y, esp[j].z);
    if(method == 1) {
                if(pflag == 1 && remark == 1) fprintf(fpout, "# Dipole momen: X, Y, Z, Total\n");
                fprintf(fpout, "%16.7E%16.7E%16.7E%16.7E\n", dipole.x, dipole.y, dipole.z, dipole.total);
                if(pflag == 1 && remark == 1) fprintf(fpout, "# Quadrupole momen (traceless): XX, YY, ZZ\n");
                fprintf(fpout, "%16.7E%16.7E%16.7E\n", quadrupole.xx, quadrupole.yy, quadrupole.zz);
                if(pflag == 1 && remark == 1) fprintf(fpout, "# Quadrupole momen (traceless): XY, XZ, YZ\n");
                fprintf(fpout, "%16.7E%16.7E%16.7E\n", quadrupole.xy, quadrupole.xz, quadrupole.yz);
    }
}

int main(int argc, char *argv[])
{
    if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
        if (argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-H") == 0)) {
            printf("[31mUsage: espgen -i  [0m input file name \n"
                   "[31m              -o  [0m output file name \n"
                   "[31m              -f  [0m input format:\n"
                   "[34m                   1 [0m Gaussian log file(default)\n" 
                   "[34m                   2 [0m Gaussian ESP file\n" 
                   "[34m                   3 [0m Gamess ESP file\n" 
                   "[31m              -p  [0m generate esp for pGM:\n"
                   "[34m                   0 [0m no, the default)\n" 
                   "[34m                   1 [0m yes \n" 
                   "[31m              -dq [0m print out dipole and quadrupole moments:\n"
                   "[34m                   0 [0m no, the default)\n" 
                   "[34m                   1 [0m yes \n"
                   "[31m              -re [0m print out remark line\n"
                   "[34m                   0 [0m no, the default)\n" 
                   "[34m                   1 [0m yes \n"); 
            exit(0);
        }
        if (argc != 5 && argc != 7 && argc != 9 && argc != 11 && argc != 13) {
            printf("[31mUsage: espgen -i  [0m input file name \n"
                   "[31m              -o  [0m output file name \n"
                   "[31m              -f  [0m input format:\n"
                   "[34m                   1 [0m Gaussian log file(default)\n" 
                   "[34m                   2 [0m Gaussian ESP file\n" 
                   "[34m                   3 [0m Gamess ESP file\n" 
                   "[31m              -p  [0m generate esp for pGM:\n"
                   "[34m                   0 [0m no, the default)\n" 
                   "[34m                   1 [0m yes \n" 
                   "[31m              -dq [0m print out dipole and quadrupole moments:\n"
                   "[34m                   0 [0m no, the default)\n" 
                   "[34m                   1 [0m yes \n"
                   "[31m              -re [0m print out remark line\n"
                   "[34m                   0 [0m no, the default)\n" 
                   "[34m                   1 [0m yes \n"); 
            exit(1);
        }
    } else {
        if (argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-H") == 0)) {
            printf("Usage: espgen -i   input file name \n"
                   "              -o   output file name \n"
                   "              -f   input format:\n"
                   "                   1  Gaussian log file(default)\n" 
                   "                   2  Gaussian ESP file\n" 
                   "                   3  Gamess ESP file\n" 
                   "              -p   generate esp for pGM:\n"
                   "                   0  no, the default)\n" 
                   "                   1  yes \n" 
                   "              -dq  print out dipole and quadrupole moments:\n"
                   "                   0  no, the default)\n" 
                   "                   1  yes \n"
                   "              -re  print out remark line\n"
                   "                   0  no, the default)\n" 
                   "                   1  yes \n"); 
            exit(1);
        }
        if (argc != 5 && argc != 7 && argc != 9 && argc != 11 && argc != 13) {
            printf("Usage: espgen -i   input file name \n"
                   "              -o   output file name \n"
                   "              -f   input format:\n"
                   "                   1  Gaussian log file(default)\n" 
                   "                   2  Gaussian ESP file\n" 
                   "                   3  Gamess ESP file\n" 
                   "              -p   generate esp for pGM:\n"
                   "                   0  no, the default)\n" 
                   "                   1  yes \n" 
                   "              -dq  print out dipole and quadrupole moments:\n"
                   "                   0  no, the default)\n" 
                   "                   1  yes \n"
                   "              -re  print out remark line\n"
                   "                   0  no, the default)\n" 
                   "                   1  yes \n"); 
        }
    }

/* allocate memory for *esp */
    maxesp = MAXESP;
    esp = (ESP *) ecalloc(maxesp, sizeof(ESP));

    method = 0;
    for (i = 1; i < argc; i += 2) {
        if (strcmp(argv[i], "-i") == 0)
            strcpy(ifilename, argv[i + 1]);
        if (strcmp(argv[i], "-o") == 0)
            strcpy(ofilename, argv[i + 1]);
        if (strcmp(argv[i], "-f") == 0)
	    format=atoi(argv[i+1]);
        if (strcmp(argv[i], "-p") == 0)
	    pflag=atoi(argv[i+1]);
        if (strcmp(argv[i], "-re") == 0)
	    remark=atoi(argv[i+1]);
        if (strcmp(argv[i], "-dq") == 0) {
            if (strcmp("YES", argv[i + 1]) == 0 || strcmp("yes", argv[i + 1]) == 0
                || strcmp("1", argv[i + 1]) == 0)
                method = 1;
            if (strcmp("Y", argv[i + 1]) == 0 || strcmp("y", argv[i + 1]) == 0)
                method = 1;
            if (strcmp("NO", argv[i + 1]) == 0 || strcmp("no", argv[i + 1]) == 0
                || strcmp("0", argv[i + 1]) == 0)
                method = 0;
            if (strcmp("N", argv[i + 1]) == 0 || strcmp("n", argv[i + 1]) == 0)
                method = 0;
        }
    }
    if(pflag != 0 && pflag != 1) pflag = 0;
    if(remark != 0 && remark != 1) remark = 0;
    if(pflag == 1 && format != 1) {
	fprintf(stderr, "only the Gaussian log file can be applied to generate pGM-supported ESP\n");
	exit(1);
    }
    fpin = efopen(ifilename, "r");
    fpout = efopen(ofilename, "w");
    nset = 0;
    for (;;) {
        if (fgets(line, MAXCHAR, fpin) == NULL)
            break;
        strcpy(tmpchar1, "");
        strcpy(tmpchar2, "");
        strcpy(tmpchar3, "");
        strcpy(tmpchar4, "");
        strcpy(tmpchar5, "");
        sscanf(&line[0], "%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3, tmpchar4, tmpchar5);

        // First check if it is a GAMESS output file
        if (strncmp(line, " $DATA", 6) == 0) {
            format = 3;
            break;
        }

/*	adding code to detect Gaussion version	*/
/*      adding more code to deal with different Gaussiaon 09 subversions*/
        if (strncmp(line, " ESP FILE -", 11) == 0) {
            format = 2;
            break;
        }
        if (strncmp(line, " Gaussian 09", 12) == 0)
            g09_index = 1;
        if (strcmp("--", tmpchar1) == 0 && strcmp("Stationary", tmpchar2) == 0
            && strcmp("point", tmpchar3) == 0 && strcmp("found.", tmpchar4) == 0)
            opt_index = 1;
        if (strcmp("ESP", tmpchar1) == 0 && strcmp("Fit", tmpchar2) == 0
            && strcmp("Center", tmpchar3) == 0)
            esp_index = 1;
        if (strcmp("Electrostatic", tmpchar1) == 0 && strcmp("Properties", tmpchar2) == 0
            && strcmp("(Atomic", tmpchar3) == 0 && strcmp("Units)", tmpchar4) == 0)
            nset++;
        if (esp_index == 1 && espvalue_index == 0)
            if (strcmp("Fit", tmpchar2) == 0 && tmpchar1[0] >= '0' && tmpchar1[0] <= '9'
                && ((tmpchar3[0] >= '0' && tmpchar3[0] <= '9') || tmpchar3[0] == '-'))
                espvalue_index = 1;
    }

    if (format == 1) {
        if (esp_index == 0) {
            if (g09_index == 0)
                fprintf(stderr,
                        "Error: No ESP fitting centers and fitting values exist, adding 'iop(6/33=2) iop(6/42=6)' to the key word list\n");
            else
                fprintf(stderr,
                        "Error: No ESP fitting centers and fitting values exist, adding 'iop(6/33=2) iop(6/42=6) iop(6/50=1)' to the keyword list\n");
            exit(1);
        }
        if (esp_index == 1 && espvalue_index == 0) {
            fprintf(stderr,
                    "Error: the ESP fitting centers exist, but the fitting values are missing\n");
            if (g09_index == 1)
                fprintf(stderr,
                        "It is recommened to generate esp file for resp fitting from the gesp file generated by adding keyword 'iop(6/50=1) in G09 input'\n");
            exit(1);
        }
    }
    for(i=0;i<MAX_ATOM_CENTER;i++) at[i].iequ = 0;
    if (pflag == 1) assign_atomtype();
    for(i=0;i<MAX_ATOM_CENTER;i++) {
        at[i].charge_mul = 0;
        at[i].charge_esp = 0;
        at[i].radius_esp = 0;
        at[i].radius_pcm = 0;
        at[i].atomicnum = 0;
    }

    if (format == 1)
        rglog();
    if (format == 2)
        rgesp();
    if (format == 3)
        rgamessdat();

    fclose(fpin);
    fclose(fpout);
    if (fail == 1) {
/*		strcat(tmpchar, ofilename);
		esystem(tmpchar);
*/
        eprintf("Cannot run espgen properly.");

    }
/*
	free(esp);
*/
    printf("\n");
    return (0);
}
