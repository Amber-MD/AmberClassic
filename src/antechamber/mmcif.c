/* Component mmcif files */
/* This file doubles as a header file and as an included (!) source file. */
/* Yes, we are on the road to multiple header files and separate compilation. */

#ifndef MMCIF_H
#define MMCIF_H
#define MAX_CIF_BLOCK_CODE_LENGTH 32
extern char blockId[MAX_CIF_BLOCK_CODE_LENGTH + 1];
#endif                          // MMCIF_H

char blockId[MAX_CIF_BLOCK_CODE_LENGTH + 1] = "Undefined_block_code__use-bk";

/*  Read a "geostd cif" file from elbow  */

#include <assert.h>
#include "../cifparse/cifparse.h"
#include "eprintf.h"

int rmmcif(char *filename, char *blockId, int *atomnum, ATOM atom[], int *bondnum,
           BOND bond[], CONTROLINFO * cinfo, MOLINFO * minfo, int flag)
{
/*if flag =1, read in atom type, if flag ==0, do not read in atom type */

    int i, numatom, numbond, iblock, iCat;
    int col1, col2, col3, col6, col7, col8, col9;
    FILE *fpin;
    char tmps[4];

    fpin = efopen(filename, "r");
    initial((*cinfo).maxatom, atom, (*minfo).resname);
    numatom = 0;
    numbond = 0;

    ndb_cif_init();
    ndb_cif_read_file(fpin);

    /*  Check for target data block.  */
    iblock = -1;
    for (i = 0; i < cifFiles.numDatablock; i++) {
        if (!strcmp(blockId, cifFiles.datablocks[i].datablockName)) {
            iblock = i;
            break;
        }
    }
    if (iblock == -1) {
        eprintf("Target data block (%s) is not in file.", blockId);
    }

    /* get chem_comp_atom atomic info: */
    iCat = get_category_index(iblock, "chem_comp_atom");
    if (iCat == -1) {
        eprintf("Target category (%s) is not in data block (%d).", "chem_comp_atom",
                iblock);
    }
    /* fprintf(stderr, "Read %d rows in category %d\n",
       cifFiles.datablocks[iblock].categories[iCat].numRow, iCat); */

    col1 = get_column_index(iblock, iCat, "type_symbol");   /* element */
    assert(col1 >= 0);
    col2 = get_column_index(iblock, iCat, "atom_id");   /* atom name */
    assert(col2 >= 0);
    col3 = get_column_index(iblock, iCat, "comp_id");   /* res. name */
    assert(col3 >= 0);
    col6 = get_column_index(iblock, iCat, "x");  /* x-coord */
    assert(col6 >= 0);
    col7 = get_column_index(iblock, iCat, "y");  /* y-coord */
    assert(col7 >= 0);
    col8 = get_column_index(iblock, iCat, "z");  /* z-coord */
    assert(col8 >= 0);
    col9 = get_column_index(iblock, iCat, "charge");  /* formal charge */
    if( col9 < 0 )
       printf("no charges in geostd file: assuming net charge of zero\n" );

    /*  --- need some better error processing here!  */
    (*minfo).usercharge = 0;
    for (i = 0; i < cifFiles.datablocks[iblock].categories[iCat].numRow; i++) {

        strcpy(atom[numatom].name,
           cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col2]);
        strcpy(atom[numatom].aa,
           cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col3]);
        atom[numatom].resno = 1;
        atom[numatom].x =
           atof(cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col6]);
        atom[numatom].y =
            atof(cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col7]);
        atom[numatom].z =
            atof(cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col8]);
        if( col9>=0 ) (*minfo).usercharge +=
            atof(cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col9]);

        numatom++;

        strcpy( tmps,
           cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col1]);
        if( strncmp(tmps,"B",1)==0 && strlen(tmps)==1){
           fprintf( stderr, "Error: found boron in residue %s\n",
           cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col3]);
           exit(1);
        }
    }

    fclose(fpin);
    *atomnum = numatom;
    return 0;
}

