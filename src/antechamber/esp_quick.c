#include "eprintf.h"

int rqout(char *filename, int *atomnum, ATOM * atom, CONTROLINFO cinfo, MOLINFO * minfo)
{
    int i;
    int numatom = 0;
    int Standard = 0;
    int read_flag = 0;
    int suc_flag = 0;
    int scan_flag = 1;
    int charge_flag = 1;
    int overflow_flag = 0;
    long input_pos = -9999;
    char tmpchar1[MAXCHAR];
    char tmpchar2[MAXCHAR];
    char tmpchar3[MAXCHAR];
    char tmpchar4[MAXCHAR];
    char tmpchar5[MAXCHAR];
    char line[MAXCHAR];
    FILE *fpin;

    fpin = efopen(filename, "r");
    initial(cinfo.maxatom, atom, (*minfo).resname);
    i = 0;
    for (;;) {
        sscanf(&line[0], "%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3, tmpchar4, tmpchar5);
        if (scan_flag == 1 && charge_flag == 1 && strcmp(tmpchar2, "MOLECULAR") == 0 && strcmp(tmpchar3, "CHARGE") == 0) {
            sscanf(&line[27], "%lf", &(*minfo).dcharge);
            sscanf(&line[65], "%d", &(*minfo).multiplicity);
            (*minfo).icharge = (int) (*minfo).dcharge;
            charge_flag = 0;
            continue;
        }
        if (fgets(line, MAXCHAR, fpin) == NULL) {
            /* now go to input_pos and read coordinates */
            scan_flag = 0;
            if (read_flag == 0 && suc_flag == 0) {
                read_flag = 1;
                if (input_pos >= 0)
                    fseek(fpin, input_pos-100, 0);
                continue;
            }
            if (read_flag == 1 && suc_flag == 0) {
                eprintf("No atom found; the quick output file may not be complete.");
            }
            break;
        }
        if (scan_flag == 1 && strcmp("--",tmpchar1) == 0 && strcmp("INPUT", tmpchar2) == 0
            && strcmp("GEOMETRY", tmpchar3) == 0) {
            input_pos = ftell(fpin);
            continue;
        }
        if (scan_flag == 0 && read_flag == 1) {
            if (suc_flag == 1 && line[17] != '.' && line[30] != '.' && line[43] != '.'){
                break;
            }else if (strlen(line) > 43 && line[17] == '.' && line[30] == '.' && line[43] == '.') {
                suc_flag = 1;
                if (overflow_flag == 0) {
                    sscanf(&line[1], "%s%lf%lf%lf", &atom[numatom].name,
                           &atom[numatom].x, &atom[numatom].y, &atom[numatom].z);
                    atom[numatom].charge = 0.0;
                }
                numatom++;
                if (numatom >= cinfo.maxatom && overflow_flag == 0) {
                    printf("Info: The number of atoms (%d) exceeded MAXATOM.\n", numatom);
                    overflow_flag = 1;
                }
            }
        } 
    }
    *atomnum = numatom;
    fclose(fpin);
    atomicnum(*atomnum, atom);
    return overflow_flag;
}

// This function creates input file to be run by the Quantum chemistry software QUICK
// Input file is a pdb file. Output is a couple of QUICK input files of the specific residue.
// The residue has to be a non-standard residue denoted by HETATM.
// Two files are created with one optimization job and other single point job reading
// data from the optimization job. Thus, the optimization job must be run first.
int pdb_to_quick_input(char *ifilename, char *ofilename, ATOM * atom, MOLINFO * minfo)
{
    FILE *fpout;
    FILE *fpin;
    char tmpchar[25];
    char tmpchar1[25];
    char tmpchar2[25];
    char tmpchar3[25];
    char tmpchar4[25];
    char tmpchar5[25];
    char tmpchar6[25];
    char tmpchar7[25];
    char tmpchar8[25];
    char tmpchar9[25];
    char tmpchar10[25];
    char tmpchar11[25];
    char line[MAXCHAR];
    char file_sp[50];
    char file_data[50];

    printf("ifilename in pdb_to_quick_input is: %s\n",ifilename);
    printf("ofilename in pdb_to_quick_input is: %s\n",ofilename);

    strncpy(file_sp, ofilename, strlen(ofilename)-3);
    strncpy(file_data, ofilename, strlen(ofilename)-3);

    file_sp[strlen(ofilename) - 3] = '\0';
    file_data[strlen(ofilename) - 3] = '\0';

    strcat(file_sp,"_sp.in");
    strcat(file_data,".dat");
    printf("file_sp: %s\n",file_sp);
    printf("file_data: %s\n",file_data);
    printf("multiplicity: %d\n",minfo->multiplicity);
    printf("charge: %d\n", minfo->usercharge);

    fpout = efopen(ofilename, "w");

    if (minfo->multiplicity > -9999){
      fprintf(fpout,"MULT=%d ",minfo->multiplicity);
    }

    if (minfo->usercharge > -9999){
      fprintf(fpout,"CHARGE=%d ",minfo->usercharge);
    }

    if (minfo->multiplicity > -9999 && minfo->multiplicity != 1){
      fprintf(fpout,"opt UDFT pbe0 BASIS=6-311+g(2d,p) write\n\n");
    }else{
      fprintf(fpout,"opt pbe0 BASIS=6-311+g(2d,p) write\n\n");
    };

    fpin = efopen(ifilename, "r");

    for (;;) {
        if (fgets(line, MAXCHAR, fpin) == NULL)
            break;

        sscanf(line,"%s%s%s%s%s%s%s%s%s%s%s%s",tmpchar,tmpchar1,tmpchar2,tmpchar3,tmpchar4,tmpchar5,tmpchar6,tmpchar7,tmpchar8,tmpchar9,tmpchar10,tmpchar11);

        if (strcmp(tmpchar, "HETATM") == 0 && strcmp(tmpchar3, (*minfo).resname) == 0){
           fprintf(fpout,"%-5s %-14s %-14s %-14s\n",tmpchar11,tmpchar6,tmpchar7,tmpchar8);
        }
    }

    fclose(fpin);
    fclose(fpout);

    fpout = efopen(file_sp, "w");

    fprintf(fpout,"$DATA = %-50s\n",file_data);

    if (minfo->multiplicity > -9999){
      fprintf(fpout,"MULT=%d ",minfo->multiplicity);
    }

    if (minfo->usercharge > -9999){
      fprintf(fpout,"CHARGE=%d ",minfo->usercharge);
    }

    if (minfo->multiplicity > -9999 && minfo->multiplicity != 1){
      fprintf(fpout,"UHF BASIS=6-31G* READ_COORD ESP_CHARGE\n\n",file_data);
    }else{
      fprintf(fpout,"HF BASIS=6-31G* READ_COORD ESP_CHARGE\n\n",file_data);
    };

    fclose(fpout);

    return 0;
 
}; // pdb_to_quick_input

// This function creates input file to be run by the Quantum chemistry software QUICK
// Input file is an xyz file. Output is a couple of QUICK input files.
// Two files are created with one optimization job and other single point job reading
// data from the optimization job. Thus, the optimization job must be run first.
int xyz_to_quick_input(char *ifilename, char *ofilename, ATOM * atom, MOLINFO * minfo)
{
    FILE   *fpout;
    FILE   *fpin;
    char   tmpchar[3];
    char   tmpchar1[50];
    char   tmpchar2[50];
    char   tmpchar3[50];
    char   line[MAXCHAR];
    char   file_sp[50];
    char   file_data[50];

    printf("ifilename in xyz_to_quick_input is: %s\n",ifilename);
    printf("ofilename in xyz_to_quick_input is: %s\n",ofilename);

    strncpy(file_sp, ofilename, strlen(ofilename)-3);
    strncpy(file_data, ofilename, strlen(ofilename)-3);

    file_sp[strlen(ofilename) - 3] = '\0';
    file_data[strlen(ofilename) - 3] = '\0';

    strcat(file_sp,"_sp.in");
    strcat(file_data,".dat");
    printf("file_sp: %s\n",file_sp);
    printf("file_data: %s\n",file_data);
    printf("multiplicity: %d\n",minfo->multiplicity);
    printf("charge: %d\n", minfo->usercharge);

    fpout = efopen(ofilename, "w");

    if (minfo->multiplicity > -9999){
      fprintf(fpout,"MULT=%d ",minfo->multiplicity);
    }

    if (minfo->usercharge > -9999){
      fprintf(fpout,"CHARGE=%d ",minfo->usercharge);
    }

    if (minfo->multiplicity > -9999 && minfo->multiplicity != 1){
      fprintf(fpout,"opt UDFT pbe0 BASIS=6-311+g(2d,p) write\n\n");
    }else{
      fprintf(fpout,"opt pbe0 BASIS=6-311+g(2d,p) write\n\n");
    };

    fpin = efopen(ifilename, "r");

    for (int i=0;;i++) {
        if (fgets(line, MAXCHAR, fpin) == NULL)
            break;

        sscanf(line,"%s%s%s%s",tmpchar,tmpchar1,tmpchar2,tmpchar3);

        if (i>1){
          fprintf(fpout,"%-5s    %s    %s    %s\n",tmpchar,tmpchar1,tmpchar2,tmpchar3);
        };

    }

    fclose(fpin);
    fclose(fpout);

    fpout = efopen(file_sp, "w");

    fprintf(fpout,"$DATA = %-50s\n",file_data);

    if (minfo->multiplicity > -9999){
      fprintf(fpout,"MULT=%d ",minfo->multiplicity);
    }

    if (minfo->usercharge > -9999){
      fprintf(fpout,"CHARGE=%d ",minfo->usercharge);
    }

    if (minfo->multiplicity > -9999 && minfo->multiplicity != 1){
      fprintf(fpout,"UHF BASIS=6-31G* READ_COORD ESP_CHARGE\n\n",file_data);
    }else{
      fprintf(fpout,"HF BASIS=6-31G* READ_COORD ESP_CHARGE\n\n",file_data);
    };

    fclose(fpout);

    return 0;
 
}; // xyz_to_quick_input


