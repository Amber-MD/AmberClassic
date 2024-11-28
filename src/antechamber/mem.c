#include "eprintf.h"


void memory(int flag, int maxatom, int maxbond, int maxring)
{
    if (flag == 0) {
        atom = (ATOM *) emalloc(sizeof(ATOM) * maxatom);
        arom = (AROM *) emalloc(sizeof(AROM) * maxatom);
        bond = (BOND *) emalloc(sizeof(BOND) * maxbond);
        int i;
        for (i = 0; i < maxbond; ++i) {
            bond[i].jflag = -1; /* bond type has not been assigned */
        }
        ring = (RING *) emalloc(sizeof(RING) * maxring);
    }
    if (flag == 3) {
        free(atom);
        free(arom);
        free(bond);
        free(ring);
        atom = (ATOM *) emalloc(sizeof(ATOM) * maxatom);
        arom = (AROM *) emalloc(sizeof(AROM) * maxatom);
        bond = (BOND *) emalloc(sizeof(BOND) * maxbond);
        int i;
        for (i = 0; i < maxbond; ++i) {
            bond[i].jflag = -1; /* bond type has not been assigned */
        }
        ring = (RING *) emalloc(sizeof(RING) * maxring);
    }
    if (flag == 1) {
        free(atom);
        free(arom);
        free(bond);
        atom = (ATOM *) emalloc(sizeof(ATOM) * maxatom);
        arom = (AROM *) emalloc(sizeof(AROM) * maxatom);
        bond = (BOND *) emalloc(sizeof(BOND) * maxbond);
        int i;
        for (i = 0; i < maxbond; ++i) {
            bond[i].jflag = -1; /* bond type has not been assigned */
        }
    }
    if (flag == 2) {
        free(ring);
        ring = (RING *) emalloc(sizeof(RING) * maxring);
    }
}
