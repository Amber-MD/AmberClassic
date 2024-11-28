/* AMBER RST */

#include "eprintf.h"


int rrst(char *filename, int *atomnum, ATOM atom[], CONTROLINFO cinfo)
{
    int i, j;
    int overflow_flag = 0;
    FILE *fpin;
    char line[MAXCHAR];

/* since rst file only has coordinates information, it can only be the additional file */
    fpin = efopen(filename, "r");
/* readin the amber rst file */
    i = 0;
    j = 0;
    for (;;) {
        if (fgets(line, MAXCHAR, fpin) == NULL) {
            if (cinfo.intstatus == 2)
                printf("Info: Finished reading file (%s).\n", filename);
            break;
        }
        j++;
        if (j == 2)
            sscanf(line, "%d", atomnum);
        if (strncmp(".", &line[4], 1) == 0 && strncmp(".", &line[16], 1) == 0
            && strncmp(".", &line[28], 1) == 0) {
            if (overflow_flag == 0)
                sscanf(&line[1], "%lf%lf%lf%lf%lf%lf", &atom[i].x, &atom[i].y, &atom[i].z,
                       &atom[i + 1].x, &atom[i + 1].y, &atom[i + 1].z);
            i = i + 2;
            if (i >= *atomnum)
                break;
            if (i >= cinfo.maxatom && overflow_flag == 0) {
                printf("Info: The number of atoms (%d) exceeded MAXATOM.\n", i);
                overflow_flag = 1;
            }

        }
    }
    fclose(fpin);
    return overflow_flag;

}

void wrst(char *filename, int atomnum, ATOM atom[])
{
    int i, j;
    FILE *fpout;
    fpout = efopen(filename, "w");
    fprintf(fpout, "MOLECULE\n");
    fprintf(fpout, "%5d\n", atomnum);
    i = 0;
    j = 0;
    while (i < atomnum) {
        fprintf(fpout, "%12.7lf%12.7lf%12.7lf", atom[i].x, atom[i].y, atom[i].z);
        i++;
        j++;
        if (j == 2) {
            fprintf(fpout, "\n");
            j = 0;
        }
    }
    fclose(fpout);
}
