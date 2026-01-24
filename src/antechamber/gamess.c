/* GAMESS .dat output */

#include "eprintf.h"


int rgamess(char *filename, int *atomnum, ATOM *atom, CONTROLINFO cinfo,
            MOLINFO *minfo)
{
    int numatom = 0, overflow_flag = 0;
    char line[MAXCHAR];
    float atomicnum;
    FILE *fpin;
    int i;

    fpin = efopen(filename, "r");
    initial(cinfo.maxatom, atom, minfo->resname);

    // Find line with number of atomic centers
    for (;;) {
        if (!fgets(line, MAXCHAR, fpin))
            eprintf("Premature end of file before electrostatic potential");

        // ELECTROSTATIC POTENTIAL COMPUTED FOR <N> ATOMS, TOTAL CHARGE= <N>
        if (strncmp(line, " ELECTROSTATIC POTENTIAL COMPUTED FOR", 37) == 0) {
            sscanf(line, "%*s%*s%*s%*s%d%*s%*s%*s%lf", &numatom, &minfo->dcharge);
            minfo->icharge = (int) minfo->dcharge;
            if (numatom >= cinfo.maxatom) {
                printf("Info: The number of atoms (%d) exceeded MAXATOM (%d).\n",
                       numatom, cinfo.maxatom);
                overflow_flag = 1;
            }
            break;
        }
    }

    // Read in each atom
    for (i = 0; i < numatom; i++) {
        if (!fgets(line, MAXCHAR, fpin))
            eprintf("Premature end of file at atom %d", i);

        // Line has index, atomicnum, x, y, z coordinates in units of Bohr
        sscanf(line, "%*d%f%lf%lf%lf", &atomicnum, &atom[i].x,
               &atom[i].y, &atom[i].z);

        atom[i].atomicnum = (int) atomicnum;
        atom[i].charge = 0.0;

        // Units are Bohr in GAMESS so convert to Angstrom
        atom[i].x *= Bohr;
        atom[i].y *= Bohr;
        atom[i].z *= Bohr;
    }

    *atomnum = numatom;
    fclose(fpin);

    // Assign elements from atomic number
    // Set names to be the same as element, for now
    initialize_elements_in_atom_to_symbols_upto_atomnum(numatom, atom);
    for (i = 0; i < numatom; i++) {
        strcpy(atom[i].name, atom[i].element);
    }

    return overflow_flag;
}


void wgamess()
{
    eprintf("Warning: Cannot write a GAMESS output file");
}
