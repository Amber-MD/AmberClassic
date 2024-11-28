// I/O functions for SMILES files.
// Assumptions:  requires external variables ifilename and intstatus.

#include <stdio.h>
#include <string.h>
#include "define.h"
#include "eprintf.h"

void rsmiles(char *filename)
{
    extern char *ifilename;
    extern int intstatus;

    char tmpchar[MAXCHAR];
    strcpy(tmpchar, "trigo concord4.02 < ");
    strcat(tmpchar, ifilename);
    if (intstatus == 2)
        fprintf(stdout, "\nRunning: %s\n", tmpchar);
    esystem(tmpchar);
    rmol2("coords.mol2");
}

void wsmiles(char *filename)
{
    printf("\n Sorry, cannot write smiles format at present");
}
