/* Orca OUT */

#include "eprintf.h"


int rorcout(char *filename, int *atomnum, ATOM * atom, CONTROLINFO cinfo, MOLINFO * minfo)
{
    int i;
    int numatom = 0;
    int read_flag = 0;
    int overflow_flag = 0;
    char line[MAXCHAR];
    FILE *fpin;

    fpin = efopen(filename, "r");
    initial(cinfo.maxatom, atom, (*minfo).resname);
    for (;;) {
        if (!fgets(line, MAXCHAR, fpin))
            eprintf("Premature end of file");
        if (read_flag == 0) {
            if ((strncmp(line,"                 *** FINAL ENERGY", 33) == 0) ||
                (strncmp(line,"                       * Single Point ", 38) == 0)) {
//            printf("FOUND_START\n");
            read_flag = 1;
            continue;
            }
            }
        if (read_flag == 1) {
            if (strncmp(line, "CARTESIAN COORDINATES (ANGSTROEM)", 33) == 0 ) {
//                printf("READING COORDS\n");
                read_flag = 2;
            }
            continue;
            }
        if ((strcmp(line, "\n") == 0 || strcmp(line,"\r\n") == 0) && read_flag == 2) {
            read_flag = 3;
        }
        if (read_flag == 2) {
            if (line[0] == '-') {
                continue;
            }
            sscanf(line, "%s%lf%lf%lf", atom[numatom].name, &atom[numatom].x,
                &atom[numatom].y, &atom[numatom].z);
//            printf("%s %d %f %f %f\n",atom[numatom].name,numatom,  atom[numatom].x,  atom[numatom].y,  atom[numatom].z);
            numatom++;
            if (numatom >= cinfo.maxatom) {
                printf("Info: The number of atoms (%d) exceeded MAXATOM (%d).\n",
                       numatom, cinfo.maxatom);
                overflow_flag = 1;
            }
            continue;
        }
        if (read_flag == 3 && strncmp(line," Total Charge           Charge", 30) == 0) {
            sscanf(line, "%*s%*s%*s%*s%d", &minfo->icharge);
            minfo->dcharge = (float) minfo->icharge;
            break;
        }
    }
    *atomnum = numatom;
    fclose(fpin);
    atomicnum(*atomnum, atom);
    initialize_elements_in_atom_to_symbols_upto_atomnum(*atomnum, atom);
    for (i = 0; i < *atomnum; i++)
        strcpy(atom[i].name, atom[i].element);
    return overflow_flag;
}


void worcout()
{
    printf("Warning: Cannot write an Orca output file.\n"
           "         You must run Orca.\n");
}
