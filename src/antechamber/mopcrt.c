/*   MOPac CaRTesian input file    */
/*      Line 1 consists of mopac keywords   */
/*      Line 2 consists of comments         */
/*      Line 3 consists of comments         */
/*      Lines 4 thru the penultimate consist of cartesian coordinates */
/*      The last line is blank.             */

/*****************************************************************************/
/*   Read a MOPac CaRTesian input file    */
/*****************************************************************************/
int rmopcrt(char *filename, int *atomnum, ATOM atom[], CONTROLINFO cinfo, MOLINFO minfo)
{
    FILE *fpin;
    int index;
    int numatom;
    int tmpint1, tmpint2, tmpint3;
    int overflow_flag = 0;
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
        sscanf(line, "%s%lf%d%lf%d%lf%d", tmpchar, &tmpfloat1, &tmpint1, &tmpfloat2,
               &tmpint2, &tmpfloat3, &tmpint3);
        if (spaceline(line) == 1)
            break;
        if (overflow_flag == 0) {
            strcpy(atom[numatom].name, tmpchar);
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
    if (cinfo.intstatus == 2)
        printf("Info: Finished reading file (%s); lines read (%d), atoms read (%d).\n",
               filename, index, numatom);
    fclose(fpin);
    *atomnum = numatom;
    return overflow_flag;
}

/*****************************************************************************/
/*   Write a MOPac CaRTesian input file    */
/*****************************************************************************/
void wmopcrt(char *filename, int atomnum, ATOM atom[], MOLINFO minfo)
{
    FILE *fpout;
    int i, nelectrons;
    char tmpkeyword[MAXCHAR];

    fpout = efopen(filename, "w");

/*      There is a bug for minfo.ekeyword; it consists of
"  qm_theory='AM1', grms_tol=0.0005,\n scfconv=1.d-10, ndiis_attempts=700, "
        which is for sqm not mopac; furthermore the newline creates an
        additional line thus generating an incorrect mopac input file;
        see below for a temporary patch.
*/
    for (i = 0; i < strlen(minfo.ekeyword); i++)
        tmpkeyword[i] = toupper(minfo.ekeyword[i]);

    if (strstr(tmpkeyword, "CHAR") == NULL) {
        fprintf(fpout, "%s", minfo.ekeyword);
        fprintf(fpout, " CHARGE=%d", minfo.icharge);
    }
    if (strstr(tmpkeyword, "DOUBLET") == NULL && minfo.multiplicity == 2)
        fprintf(fpout, " DOUBLET");
    if (strstr(tmpkeyword, "TRIPLET") == NULL && minfo.multiplicity == 3)
        fprintf(fpout, " TRIPLET");
/*
        The original correct statement:
    fprintf(fpout, "\ncreated by wmopcrt() for mopac\n\n");
        Temporary patch for minfo.ekeyword bug above:
*/
    fprintf(fpout, "\n\ncreated by wmopcrt() for mopac\n");

    initialize_elements_in_atom_to_symbols_upto_atomnum(atomnum, atom);
    nelectrons = 0;
    for (i = 0; i < atomnum; i++) {
        fprintf(fpout, "%5s%12.4lf  1  %12.4lf  1  %12.4lf  1   \n", atom[i].element,
                atom[i].x, atom[i].y, atom[i].z);
/*  check that the number of electrons is even:   */
        nelectrons += atom[i].atomicnum;
    }
    fprintf(fpout, "\n");
    fclose(fpout);
    nelectrons -= minfo.icharge;
    printf("Info: Total number of electrons: %d; net charge: %d\n", nelectrons,
            minfo.icharge);
    if (nelectrons % 2 != 0) {
        printf("Info: The number of electrons is odd (%d).\n"
               "      Please check the total charge (-nc flag) and spin multiplicity (-m flag).\n",
               nelectrons);
    }
}
