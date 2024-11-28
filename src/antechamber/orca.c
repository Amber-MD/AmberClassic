/* ORCA files, modded from GCRT, SSchott */
/* TODO: Adapt missing functionalities from Gaussian */

#include "eprintf.h"


int rorca(char *filename, int *atomnum, ATOM * atom, CONTROLINFO cinfo, MOLINFO minfo)
{
    FILE *fpin;
    int numatom;
    int read_charge = 1;
    int overflow_flag = 0;
    char tmpchar[MAXCHAR];
    char line[MAXCHAR];


    fpin = efopen(filename, "r");
    initial(cinfo.maxatom, atom, minfo.resname);
    numatom = 0;
    for (;;) {
        if (fgets(line, MAXCHAR, fpin) == NULL) {
            if (cinfo.intstatus == 2)
                printf("Info: Finished reading file (%s).\n", filename);
            break;
        }
        sscanf(line, "%s", tmpchar);
        if (tmpchar[0] == '*') {
            sscanf(line, "%*s%s%d%d", tmpchar, &(minfo.icharge), &(minfo.multiplicity));
            if (strcmp("xyz",tmpchar) != 0) {
                sscanf(line, "%*s%d%d", &(minfo.icharge), &(minfo.multiplicity));
                }
            read_charge = 0;
            continue;
        }
        if (strcmp("*",tmpchar) == 0 && read_charge == 0) {
            break;
        }
        if (overflow_flag == 0) {
            if (read_charge == 0) {
                sscanf(line, "%s%lf%lf%lf", atom[numatom].name, &atom[numatom].x,
                       &atom[numatom].y, &atom[numatom].z);
                numatom++;
                }
        }
        if (numatom >= cinfo.maxatom && overflow_flag == 0) {
            printf("Info: The number of atoms (%d) exceeded MAXATOM.\n", numatom);
            overflow_flag = 1;
        }
    }
    if (cinfo.intstatus == 2)
        printf("Info: Finished reading file (%s); atoms read (%d).\n",
               filename, numatom);
    fclose(fpin);

    *atomnum = numatom;
    atomicnum(*atomnum, atom);
/* printf("\n atom number is  %5d", *atomnum); */
    return overflow_flag;
}

void worca(char *filename, int atomnum, ATOM atom[], MOLINFO minfo)
{
// TODO: Left some backbone of GCRT to potentially adapt to ORCA SSchott
    FILE *fpin;
    FILE *fpout;
    char *amberhome;
    char espparm_file[MAXCHAR];
    char line[MAXCHAR];
    char akeyword[MAXCHAR] = "";
    int i, j;
    int iradius_flag;
    int ibs = 0;
    int nespparm = 0;
    int nbasisset = 0;
    double default_radius = 1.7;
    ESPPARM espparm[120];
    BASISSET basisset[100];

    /* int index; */
    fpout = efopen(filename, "w");
    if (strlen(minfo.gn) >= 5)
        fprintf(fpout, "%s\n", minfo.gn);
    if (strlen(minfo.gm) >= 4) {
        fprintf(fpout, "%s\n", minfo.gm);
    }
    else {
        fprintf(fpout, "%s\n", "%MaxCore 4000");
    }

//      when the default gaussian keyword is used, or esp_flag ==1, read ESP.PARM
    amberhome = egetenv("AMBERCLASSICHOME");
    strcpy(espparm_file, amberhome);
    strcat(espparm_file, "/dat/antechamber/ESPPARM.DAT");
    fpin = efopen(espparm_file, "r");
    for (;;) {
        if (fgets(line, MAXCHAR, fpin) == NULL)
            break;
        if (strncmp(line, "DEFAULT RADIUS", 14) == 0)
            sscanf(&line[14], "%lf", &default_radius);
        if (strncmp(line, "BASIS SET", 9) == 0) {
            sscanf(&line[9], "%d%s", &basisset[nbasisset].id, basisset[nbasisset].bs);
            nbasisset++;
        }
        if (strncmp(line, "PARM", 4) == 0) {
            sscanf(&line[4], "%d%s%lf%lf%d%d", &espparm[nespparm].atomicnum,
                   espparm[nespparm].elem, &espparm[nespparm].vdw,
                   &espparm[nespparm].mk, &espparm[nespparm].flag,
                   &espparm[nespparm].bs);
            nespparm++;
        }
    }
    fclose(fpin);

    iradius_flag = 0;
    ibs = 0;
    for (i = 0; i < atomnum; i++) {
        for (j = 0; j < nespparm; j++)
            if (atom[i].atomicnum == espparm[j].atomicnum
                || strcmp(atom[i].element, espparm[j].elem) == 0) {
                if (ibs < espparm[j].bs)
                    ibs = espparm[j].bs;
                if (espparm[j].flag != 0) {
                    iradius_flag = 1;
                    espparm[j].flag = 2;
                }
                break;
            }
    }

    strcpy(minfo.gkeyword, "! HF 6-31G(d)");
    strcat(minfo.gkeyword, " TightSCF TightOpt KeepDens");
//      additional keywords
	if(minfo.gkeyword[0] == '!')
    		fprintf(fpout, "%s\n", minfo.gkeyword);
	else
    		fprintf(fpout, "!%s\n", minfo.gkeyword);
    	if (strlen(akeyword) >= 1) 
        	fprintf(fpout, "!%s\n", akeyword);

    if (minfo.igdsk == 1)
    	fprintf(fpout, "%s\n", minfo.gdsk);
    fprintf(fpout, "%%id \"%s\"\n", "molecule");
    fprintf(fpout, "\n");
    fprintf(fpout, " * xyz %d %d\n", minfo.icharge, minfo.multiplicity);
    initialize_elements_in_atom_to_symbols_upto_atomnum(atomnum, atom);

    for (i = 0; i < atomnum; i++)
//              fprintf(fpout, "%5s%12.4lf    %12.4lf    %12.4lf     \n",
//                              atom[i].element, atom[i].x, atom[i].y, atom[i].z);
        fprintf(fpout, " %-5s %16.10lf %16.10lf %16.10lf \n", atom[i].element,
                atom[i].x, atom[i].y, atom[i].z);
    fprintf(fpout, " * \n");
    fprintf(fpout, "\n");
    fclose(fpout);
}
