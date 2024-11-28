/* divcon CRT */

int rsqmcrt(char *filename, int *atomnum, ATOM atom[], CONTROLINFO cinfo, MOLINFO minfo)
{
    FILE *fpin;
    int read_flag = 0;
    int numatom;
    int overflow_flag = 0;
    int i, j;
    int icharge;
    char line[MAXCHAR];
    char tmpchar1[MAXCHAR];
    char tmpchar2[MAXCHAR];
    double tmpfloat1, tmpfloat2, tmpfloat3;

    fpin = efopen(filename, "r");
    initial(cinfo.maxatom, atom, minfo.resname);
    numatom = 0;
    for (;;) {
        if (fgets(line, MAXCHAR, fpin) == NULL) {
            break;
        }
        sscanf(line, "%s%s", tmpchar1, tmpchar2);
        if (read_flag == 0 && strstr(line, "qmcharge") != 0)
            for (i = 7; i < strlen(line); i++) {
                if (line[i - 7] == 'q' && line[i - 6] == 'm' && line[i - 5] == 'c'
                    && line[i - 4] == 'h' && line[i - 3] == 'a' && line[i - 2] == 'r'
                    && line[i - 1] == 'g' && line[i] == 'e') {
                    for (j = i; j < strlen(line); j++)
                        if (line[j] == '=') {
                            line[j] = ' ';
                            sscanf(&line[j], "%d", &icharge);
                            minfo.icharge = icharge;
                            minfo.dcharge = icharge;
                            break;
                        }
                    break;
                }
            }

        if (strcmp(tmpchar1, "/") == 0) {
            read_flag = 1;
            continue;
        }
        if (read_flag == 1 && strlen(line) < 5)
            break;
        if (read_flag == 1) {
            sscanf(line, "%s%s%lf%lf%lf", tmpchar1, tmpchar2, &tmpfloat1, &tmpfloat2,
                   &tmpfloat3);
            if (overflow_flag == 0) {
                strcpy(atom[numatom].name, tmpchar2);
                atom[numatom].x = tmpfloat1;
                atom[numatom].y = tmpfloat2;
                atom[numatom].z = tmpfloat3;
            }
            numatom++;
            if (numatom >= cinfo.maxatom && overflow_flag == 0) {
                printf("Info: The number of atoms (%d) exceeded MAXATOM.\n", numatom);
                overflow_flag = 1;
            }
        }
    }
    if (cinfo.intstatus == 2)
        printf("Info: Finished reading file (%s); atoms read (%d).\n",
               filename, numatom);
    fclose(fpin);
    *atomnum = numatom;
    return overflow_flag;
}

void wsqmcrt(char *filename, int atomnum, ATOM atom[], MOLINFO minfo)
{
    FILE *fpout;
    int i, nelectrons;

    fpout = efopen(filename, "w");

    /*  write initial keywords and control parameters:  */
    fprintf(fpout, "Run semi-empirical minimization\n");
    fprintf(fpout, " &qmmm\n");
    fprintf(fpout, "  %s  qmcharge=%d,\n /\n", minfo.ekeyword, minfo.icharge);

    /* initialize_elements_in_atom_to_symbols_upto_atomnum(atomnum, atom); */

    nelectrons = 0;
    for (i = 0; i < atomnum; i++) {
        fprintf(fpout, "%4d %5s  %12.4lf  %12.4lf  %12.4lf \n", atom[i].atomicnum,
                atom[i].name, atom[i].x, atom[i].y, atom[i].z);
        nelectrons += atom[i].atomicnum;
    }
    fprintf(fpout, "\n");
    fclose(fpout);
    nelectrons -= minfo.icharge;
    /*  check that the number of electrons is even:   */
    if (nelectrons % 2 != 0) {
        printf("Info: Total number of electrons: %d; net charge: %d\n", 
          nelectrons, minfo.icharge);
        printf("Info: The number of electrons is odd (%d).\n"
               "      Please check the total charge (-nc flag) and spin multiplicity (-m flag).\n",
               nelectrons);
    }
}
