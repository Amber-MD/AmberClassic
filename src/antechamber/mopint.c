/* MOPAC INT */

#include "eprintf.h"


int rmopint(char *filename, int *atomnum, ATOM * atom, CONTROLINFO cinfo, MOLINFO minfo)
{
    FILE *fpin;
    int index;
    int numatom;
    int overflow_flag = 0;
    int tmpint1, tmpint2, tmpint3;
    int tmpint4, tmpint5, tmpint6;
    char line[MAXCHAR];
    char tmpchar[20];
    double tmpfloat1, tmpfloat2, tmpfloat3;

    fpin = efopen(filename, "r");
    initial(cinfo.maxatom, atom, minfo.resname);
    index = 0;
    numatom = 0;
    for (;;) {
        if (fgets(line, MAXCHAR, fpin) == NULL) {
            break;
        }
        index++;
        if (index < 4)
            continue;
        sscanf(line, "%s%lf%d%lf%d%lf%d%d%d%d", tmpchar, &tmpfloat1, &tmpint1, &tmpfloat2,
               &tmpint2, &tmpfloat3, &tmpint3, &tmpint4, &tmpint5, &tmpint6);
        if (spaceline(line) == 1)
            break;
        if (overflow_flag == 0) {
            strcpy(atom[numatom].name, tmpchar);
            atom[numatom].bond = tmpfloat1;
            atom[numatom].bondatom = tmpint4 - 1;
            atom[numatom].angle = tmpfloat2;
            atom[numatom].angleatom = tmpint5 - 1;
            atom[numatom].twist = tmpfloat3;
            atom[numatom].twistatom = tmpint6 - 1;

            if (atom[numatom].bondatom >= numatom) {
                eprintf("bond atom ID is larger than ID of current atom (ID: %d, Name: %s).",
                     numatom + 1, atom[numatom].name);
            }
            if (atom[numatom].angleatom >= numatom) {
                eprintf("angle atom ID is larger than ID of current atom (ID: %d, Name: %s).",
                     numatom + 1, atom[numatom].name);
            }
            if (atom[numatom].twistatom >= numatom) {
                eprintf("torsional atom ID is larger than ID of current atom (ID: %d, Name: %s).",
                     numatom + 1, atom[numatom].name);
            }
        }
        numatom++;
        if (numatom >= cinfo.maxatom && overflow_flag == 0) {
            printf("Info: The number of atoms (%d) exceeded MAXATOM.\n",
                   numatom);
            overflow_flag = 1;
        }
    }
    if (cinfo.intstatus == 2)
        printf("Info: Finished reading file (%s); lines read (%d), atoms read (%d).\n",
               filename, index, numatom);
    fclose(fpin);
    *atomnum = numatom;
    return overflow_flag;
}

void wmopint(char *filename, int atomnum, ATOM atom[], MOLINFO minfo)
{
    FILE *fpout;
    int i;

    fpout = efopen(filename, "w");
    intercoord(atomnum, atom, minfo.tor);
    fprintf(fpout, "%s", minfo.ekeyword);
    fprintf(fpout, " CHARGE=%d\n", minfo.icharge);
    fprintf(fpout, "%s\n", "remark line goes here\n");
    initialize_elements_in_atom_to_symbols_upto_atomnum(atomnum, atom);
    for (i = 0; i < atomnum; i++)
        fprintf(fpout, "%5s%12.4lf  1  %12.4lf  1  %12.4lf  1  %5d%5d%5d \n",
                atom[i].element, atom[i].bond, atom[i].angle, atom[i].twist,
                atom[i].bondatom + 1, atom[i].angleatom + 1, atom[i].twistatom + 1);
    fprintf(fpout, "\n");
    fclose(fpout);
}
