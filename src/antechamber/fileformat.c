// File format input processing and validation.

// This file doubles as a header file and as an included (!) source file.
// Yes, we are on the road to multiple header files and separate compilation.

#ifndef FILEFORMAT_H
#define FILEFORMAT_H

extern int checkformat;         // check file format in the manner of acdoctor

#endif                          // FILEFORMAT_H


#ifndef FILEFORMAT_C__YES_THAT_IS_C
#define FILEFORMAT_C__YES_THAT_IS_C

# include <assert.h>
# include <stdio.h>
# include <stdlib.h>
# include "define.h"
# include "atom.h"
// # include "mmcif.h"  // for blockId
# include "checkmolecule.c"  // change to .h when split
# include "eprintf.h"

int checkformat = 1;            // check file format in the manner of acdoctor


// Check the input file syntax based on its file format.

void check_input_file_format(char *filename, char *format)
{
    int imol, iatom, ibond;     /* for mol2 format */
    char line[MAXCHAR];
    char tmpchar[MAXCHAR];
    char tmpchar1[MAXCHAR];
    char tmpchar2[MAXCHAR];
    char tmpchar3[MAXCHAR];
    char FIVEBLANKS[] = "     ";
    char FOURBLANKS[] = "    ";
    char THREEBLANKS[] = "   ";
    int ieps;                   /* for gout format */

    if ( !checkformat) {
        return;
    }
    FILE *fpin;
    if ((fpin = fopen(filename, "r")) == NULL) {
        eprintf("Cannot open the input file (%s).", filename);
    }

    if (strcmp(format, "mol2") == 0) {
        imol = 0;
        iatom = 0;
        ibond = 0;
        printf("-- Check Format for mol2 File --\n");
        for (;;) {
            if (fgets(line, MAXCHAR, fpin) == NULL)
                break;
            sscanf(line, "%s", tmpchar);
            if (strcmp(tmpchar, "@<TRIPOS>MOLECULE") == 0)
                imol = 1;
            if (strcmp(tmpchar, "@<TRIPOS>ATOM") == 0)
                iatom = 1;
            if (strcmp(tmpchar, "@<TRIPOS>BOND") == 0)
                ibond = 1;
            if (imol + iatom + ibond == 3)
                break;
        }
        if (imol == 0) {
            eprintf("No @@<TRIPOS>MOLECULE field.");
        }
        if (iatom == 0) {
            eprintf("No @@<TRIPOS>ATOM field.");
        }
        if (ibond == 0) {
            eprintf("No @@<TRIPOS>BOND field.");
        }

        printf("   Status: pass\n");
    } else if (strcmp(format, "mopint") == 0) {

    } else if (strcmp(format, "mopcrt") == 0) {

    } else if (strcmp(format, "mopout") == 0) {

    } else if (strcmp(format, "gcrt") == 0) {

    } else if (strcmp(format, "orcinp") == 0) {

    } else if (strcmp(format, "orcout") == 0) {

    } else if (strcmp(format, "gzmat") == 0) {

    } else if (strcmp(format, "gout") == 0) {
        printf("-- Check Format for Gaussian Output File --\n");
        ieps = 0;
        for (;;) {
            if (fgets(line, MAXCHAR, fpin) == NULL)
                break;
            sscanf(line, "%s%s%s", tmpchar1, tmpchar2, tmpchar3);
            if (strcmp(tmpchar1, "ESP") == 0 && strcmp(tmpchar2, "Fit") == 0
                && strcmp(tmpchar3, "Center") == 0) {
                ieps = 1;
                break;
            }
        }
        if (ieps == 0) {
            // color on "[33m
            printf("Warning: No ESP information in the Gaussian output file.\n"
                   "         This file cannot be used to generate RESP charges with Gaussian via\n"
                   "         e.g. \"#HF/6-31G* SCF=tight Pop=MK iop(6/33=2) iop(6/42=6) opt\"\n");
        }
        printf("   Status: pass\n");
    } else if (strcmp(format, "jcrt") == 0) {

    } else if (strcmp(format, "jzmat") == 0) {

    } else if (strcmp(format, "jout") == 0) {

    } else if (strcmp(format, "pdb") == 0 || strcmp(format, "mpdb") == 0
               || strcmp(format, "ac") == 0) {
        printf("-- Check Format for %s File --\n", format);
        for (;;) {
            if (fgets(line, MAXCHAR, fpin) == NULL)
                break;
            if (strncmp(line, "ATOM", 4) == 0) {
                if (strncmp(&line[6], FIVEBLANKS, sizeof(FIVEBLANKS) - 1) == 0) {
                    // color on "[31m
                    eprintf("Atom IDs must be in Columns 7-11.");
                }
                if (strncmp(&line[12], FOURBLANKS, sizeof(FOURBLANKS) - 1) == 0) {
                    eprintf("Atom Names must be in Columns 13-16.");
                }
                if (strncmp(&line[17], THREEBLANKS, sizeof(THREEBLANKS) - 1) == 0) {
                    eprintf("Residue Names must be in Columns 18-20.");
                }
                if (strncmp(&line[22], FOURBLANKS, sizeof(FOURBLANKS) - 1) == 0) {
                    eprintf("Residue IDs must be in Columns 23-26.");
                }
                if (line[34] != '.' || line[42] != '.' || line[50] != '.') {
                    eprintf("Coordinates must be in Columns 31-38, 39-46 and 47-54 in %%8.3format.");
                }
            }
        }
        printf("   Status: pass\n");
    } else if (strcmp(format, "csd") == 0) {

    } else if (strcmp(format, "mdl") == 0 || strcmp(format, "sdf") == 0) {
        printf("-- Check Format for %s File --\n", format);
        fgets(line, MAXCHAR, fpin);
        fgets(line, MAXCHAR, fpin);
        if (line[20] == '2' && (line[21] == 'D' || line[21] == 'd')) {
            eprintf("antechamber cannot handle 2D format; try to fill in open valences.");
        }
        fgets(line, MAXCHAR, fpin);
        fgets(line, MAXCHAR, fpin);
        // assume number of atoms is required.
        int atoms = -999;
        // assume number of bonds is not required and defaults to 0.
        int bonds = 0;
        sscanf(line, "%3d%3d", &atoms, &bonds);
        if (atoms < 1 || atoms > 999) {
            eprintf("Invalid number of atoms (%d); MDL SDF supports at most 999 atoms.", atoms);
        }
        if (bonds < 0 || bonds > 999) {
            eprintf("Invalid number of bonds (%d); MDL SDF supports at most 999 bonds.", bonds);
        }
        if (line[6] >= '0' && line[6] <= '9') {
            eprintf("Invalid format: character 7 of line 4 is a digit.");
        }
        printf("   Status: pass\n");
    } else if (strcmp(format, "alc") == 0) {

    } else if (strcmp(format, "hin") == 0) {

    } else if (strcmp(format, "prepi") == 0) {

    } else if (strcmp(format, "prepc") == 0) {

    } else if (strcmp(format, "divcrt") == 0) {

    } else if (strcmp(format, "divout") == 0) {

    // For GAMESS, just check the first line is DATA
    } else if (strcmp(format, "gamess") == 0) {
        fgets(line, MAXCHAR, fpin);
        if (strncmp(line, " $DATA", 6)) {
            eprintf("Invalid format: GAMESS files should begin with ' $DATA'");
        }
    }
}


// Allocate memory for some set of arom, atom, bond, and ring.

void memory(int flag, int maxatom, int maxbond, int maxring)
{
    extern AROM *arom;
    extern ATOM *atom;
    extern BOND *bond;
    extern RING *ring;

    if (flag == 0) {
        atom = (ATOM *) emalloc(sizeof(ATOM) * maxatom);
        arom = (AROM *) emalloc(sizeof(AROM) * maxatom);
        bond = (BOND *) emalloc(sizeof(BOND) * maxbond);
        int i;
        for (i = 0; i < maxbond; ++i) {
            bond[i].jflag = -1; /* bond type has not been assigned */
        }
    }
/*flag = 1  <->atom
       = 2  <->bond
       = 3  <->arom
       = 4  <->atom + bond
       = 5  <->atom + arom 
       = 6  <->bond + arom
       = 7  <->atom + arom +bond
*/
    if (flag == 1 || flag == 4 || flag == 5 || flag == 7) {
        free(atom);
        atom = (ATOM *) emalloc(sizeof(ATOM) * maxatom);
    }
    if (flag == 2 || flag == 4 || flag == 6 || flag == 7) {
        free(bond);
        bond = (BOND *) emalloc(sizeof(BOND) * maxbond);
        int i;
        for (i = 0; i < maxbond; ++i) {
            bond[i].jflag = -1; /* bond type has not been assigned */
        }
    }
    if (flag == 3 || flag == 5 || flag == 6 || flag == 7) {
        free(arom);
        arom = (AROM *) emalloc(sizeof(AROM) * maxatom);
    }
    if (flag == 8) {
        free(ring);
        ring = (RING *) emalloc(sizeof(RING) * maxring);
    }
}


// Set program wide flags and read input file based on file format.

void read_and_validate_input_file()
{
    extern int atomtype_flag;   /*judge atom type? */
    extern int bondtype_flag;   /*judge bond type? */
    extern int default_flag;    /*assign default information? */
    extern int atomname_flag;   /*assign atom name? */
    extern int atomicnum_flag;  /*judge atomic number according to atom name ? */
    extern int adjustatomname_flag;     /*adjust atom name? */
    extern int usr_aan_flag;    /*the user's input for adjusting atom names */
    extern int cartcoord_flag;  /*generate coordinate from internal coordinate ? */
    extern int connect_flag;    /*judge atom connectivity and generate bond ? */
    extern int atomnum;         // 1 more than the highest valid atom array index.
    extern int bondnum;
    extern ATOM *atom;          // the atom array.
    extern BOND *bond;
    extern CONTROLINFO cinfo;
    extern MOLINFO minfo;
    extern char ifilename[];

    int i, j, k;
    int index;
    int overflow_flag = 0;      /*if overflow_flag ==1, reallocate memory */

/*****************************************************************************/
/* Begin big block over input file formats, its big, uh, very big */
/* This very very big block reads in the input file */
/*****************************************************************************/

    if (strcmp("ac", cinfo.intype) == 0 || strcmp("1", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag = rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag =
                rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
        }
        adjustatomname_flag = 1;
        if (usr_aan_flag == 0)
            adjustatomname_flag = 0;
        atomicnum_flag = 1;
        checkbyatomtype = 1;
        checkbybondtype = 1;
    } else if (strcmp("mol2", cinfo.intype) == 0 || strcmp("2", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag =
            rmol2(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo, 0);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag =
                rmol2(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo, 0);
        }
        default_flag = 1;
        atomicnum_flag = 1;
        adjustatomname_flag = 1;
        if (usr_aan_flag == 0)
            adjustatomname_flag = 0;
        checkbyatomtype = 1;
        checkbybondtype = 1;
    } else if (strcmp("ccif", cinfo.intype) == 0 || strcmp("27", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag =
            rmmcif(ifilename, blockId, &atomnum, atom, &bondnum, bond, &cinfo, &minfo, 0);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag =
                rmmcif(ifilename, blockId, &atomnum, atom, &bondnum, bond, &cinfo, &minfo,
                       0);
        }
        default_flag = 1;
        atomicnum_flag = 1;
        adjustatomname_flag = 1;
        if (usr_aan_flag == 0)
            adjustatomname_flag = 0;
        atomtype_flag = 1;
        bondtype_flag = 2;
        connect_flag = 1;
    } else if (strcmp("mopint", cinfo.intype) == 0 || strcmp("9", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag = rmopint(ifilename, &atomnum, atom, cinfo, minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rmopint(ifilename, &atomnum, atom, cinfo, minfo);
        }
        cartcoord_flag = 1;
        atomicnum_flag = 1;
        atomname_flag = 1;
        default_flag = 2;
        connect_flag = 1;
        bondtype_flag = 2;

    } else if (strcmp("mopcrt", cinfo.intype) == 0 || strcmp("10", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag = rmopcrt(ifilename, &atomnum, atom, cinfo, minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rmopcrt(ifilename, &atomnum, atom, cinfo, minfo);
        }
        atomicnum_flag = 1;
        atomname_flag = 1;
        default_flag = 2;
        connect_flag = 1;
        bondtype_flag = 2;
    } else if (strcmp("mopout", cinfo.intype) == 0 || strcmp("12", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag = rmopout(ifilename, &atomnum, atom, &cinfo, &minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rmopout(ifilename, &atomnum, atom, &cinfo, &minfo);
        }
        atomicnum_flag = 1;
        atomname_flag = 1;
        default_flag = 2;
        connect_flag = 1;
        bondtype_flag = 2;
    } else if (strcmp("orcinp", cinfo.intype) == 0 || strcmp("29", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag = rorca(ifilename, &atomnum, atom, cinfo, minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rorca(ifilename, &atomnum, atom, cinfo, minfo);
        }
        atomname_flag = 1;
        default_flag = 2;
        connect_flag = 1;
        bondtype_flag = 2;
    } else if (strcmp("gcrt", cinfo.intype) == 0 || strcmp("8", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag = rgcrt(ifilename, &atomnum, atom, cinfo, minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rgcrt(ifilename, &atomnum, atom, cinfo, minfo);
        }
        atomname_flag = 1;
        default_flag = 2;
        connect_flag = 1;
        bondtype_flag = 2;
    } else if (strcmp("gzmat", cinfo.intype) == 0 || strcmp("7", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag = rgzmat(ifilename, &atomnum, atom, cinfo, minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rgzmat(ifilename, &atomnum, atom, cinfo, minfo);
        }
        cartcoord_flag = 1;
        atomicnum_flag = 1;
        atomname_flag = 1;
        default_flag = 2;
        connect_flag = 1;
        bondtype_flag = 2;

     } else if (strcmp("orcout", cinfo.intype) == 0 || strcmp("30", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag = rorcout(ifilename, &atomnum, atom, cinfo, &minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rorcout(ifilename, &atomnum, atom, cinfo, &minfo);
        }
        atomname_flag = 1;
        default_flag = 2;
        connect_flag = 1;
        bondtype_flag = 2;
   
    } else if (strcmp("gout", cinfo.intype) == 0 || strcmp("11", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag = rgout(ifilename, &atomnum, atom, cinfo, &minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rgout(ifilename, &atomnum, atom, cinfo, &minfo);
        }
        atomname_flag = 1;
        default_flag = 2;
        connect_flag = 1;
        bondtype_flag = 2;

    } else if (strcmp("gamess", cinfo.intype) == 0 || strcmp("27", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag = rgamess(ifilename, &atomnum, atom, cinfo, &minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rgamess(ifilename, &atomnum, atom, cinfo, &minfo);
        }
        atomname_flag = 1;
        default_flag = 2;
        connect_flag = 1;
        bondtype_flag = 2;

    } else if (strcmp("jcrt", cinfo.intype) == 0 || strcmp("18", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag = rjcrt(ifilename, &atomnum, atom, cinfo, minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rjcrt(ifilename, &atomnum, atom, cinfo, minfo);
        }
        default_flag = 2;
        connect_flag = 1;
        bondtype_flag = 2;
    } else if (strcmp("jzmat", cinfo.intype) == 0 || strcmp("19", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag = rjzmat(ifilename, &atomnum, atom, cinfo, minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rjzmat(ifilename, &atomnum, atom, cinfo, minfo);
        }
        cartcoord_flag = 1;
        atomicnum_flag = 1;
        atomname_flag = 0;
        default_flag = 2;
        connect_flag = 1;
        bondtype_flag = 2;
    } else if (strcmp("jout", cinfo.intype) == 0 || strcmp("20", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag = rjout(ifilename, &atomnum, atom, cinfo, &minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rjout(ifilename, &atomnum, atom, cinfo, &minfo);
        }
        default_flag = 2;
        connect_flag = 1;
        bondtype_flag = 2;
    } else if (strcmp("pdb", cinfo.intype) == 0 || strcmp("3", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag = rpdb(ifilename, &atomnum, atom, cinfo, minfo, 0);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rpdb(ifilename, &atomnum, atom, cinfo, minfo, 0);
        }
        adjustatomname_flag = 1;
        if (usr_aan_flag == 0)
            adjustatomname_flag = 0;
        atomicnum_flag = 1;
        default_flag = 2;
        bondtype_flag = 2;
        index = 0;
        for (i = 0; i < atomnum; i++)
            if (atom[i].connum > 0) {
                index = 1;
                break;
            }
        if (index == 1) {
            bondnum = 0;
            for (i = 0; i < atomnum - 1; i++)
                for (j = i + 1; j < atomnum; j++)
                    for (k = 0; k < 6; k++) {
                        if (atom[i].con[k] == -1)
                            break;
                        if (atom[i].con[k] == j) {
                            bond[bondnum].bondi = i;
                            bond[bondnum].bondj = j;
                            bondnum++;
                        }
                    }
        } else
            connect_flag = 1;
    } else if (strcmp("mpdb", cinfo.intype) == 0 || strcmp("4", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag = rpdb(ifilename, &atomnum, atom, cinfo, minfo, 1);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rpdb(ifilename, &atomnum, atom, cinfo, minfo, 1);
        }
        adjustatomname_flag = 1;
        if (usr_aan_flag == 0)
            adjustatomname_flag = 0;
        atomicnum_flag = 1;
        default_flag = 1;
        bondtype_flag = 2;
        index = 0;
        for (i = 0; i < atomnum; i++)
            if (atom[i].connum > 0) {
                index = 1;
                break;
            }
        if (index == 1) {
            bondnum = 0;
            for (i = 0; i < atomnum - 1; i++)
                for (j = i + 1; j < atomnum; j++)
                    for (k = 0; k < 6; k++) {
                        if (atom[i].con[k] == -1)
                            break;
                        if (atom[i].con[k] == j) {
                            bond[bondnum].bondi = i;
                            bond[bondnum].bondj = j;
                            bondnum++;
                        }
                    }
        } else
            connect_flag = 1;

    } else if (strcmp("csd", cinfo.intype) == 0 || strcmp("14", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag = rcsd(ifilename, &atomnum, atom, cinfo, minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rcsd(ifilename, &atomnum, atom, cinfo, minfo);
        }
        atomicnum_flag = 1;
        default_flag = 2;
        connect_flag = 1;
        bondtype_flag = 2;
    } else if (strcmp("mdl", cinfo.intype) == 0 || strcmp("sd", cinfo.intype) == 0
               || strcmp("sdf", cinfo.intype) == 0 || strcmp("15", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag = rmdl(ifilename, &atomnum, atom, &bondnum, bond, cinfo, minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rmdl(ifilename, &atomnum, atom, &bondnum, bond, cinfo, minfo);
        }
        atomicnum_flag = 1;
        atomname_flag = 1;
        default_flag = 2;
        bondtype_flag = 1;
        checkbybondtype = 1;
    } else if (strcmp("alc", cinfo.intype) == 0 || strcmp("13", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag = ralc(ifilename, &atomnum, atom, &bondnum, bond, cinfo, minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = ralc(ifilename, &atomnum, atom, &bondnum, bond, cinfo, minfo);
        }
        atomicnum_flag = 1;
        atomname_flag = 1;
        bondtype_flag = 1;
        default_flag = 2;
        checkbybondtype = 1;
    } else if (strcmp("hin", cinfo.intype) == 0 || strcmp("16", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag = rhin(ifilename, &atomnum, atom, &bondnum, bond, cinfo, minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rhin(ifilename, &atomnum, atom, &bondnum, bond, cinfo, minfo);
        }
        atomicnum_flag = 1;
        atomname_flag = 1;
        bondtype_flag = 1;
        default_flag = 2;
        checkbybondtype = 1;
    } else if (strcmp("prepi", cinfo.intype) == 0 || strcmp("5", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag = rprepi(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag =
                rprepi(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
        }
        atomicnum_flag = 1;
        bondtype_flag = 2;
        connect_flag = 0;
        checkbyatomtype = 1;
    } else if (strcmp("prepc", cinfo.intype) == 0 || strcmp("6", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag = rprepc(ifilename, &atomnum, atom, &cinfo, &minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rprepc(ifilename, &atomnum, atom, &cinfo, &minfo);
        }
        atomicnum_flag = 1;
        bondtype_flag = 2;
        connect_flag = 1;
        checkbyatomtype = 1;
    } else if (strcmp("rst", cinfo.intype) == 0 || strcmp("17", cinfo.intype) == 0) {
        eprintf("File format (%s) can only be an additional "
                "file because it only has coordinate info.", cinfo.intype);
    } else if (strcmp("divcrt", cinfo.intype) == 0 || strcmp("21", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag = rdivcrt(ifilename, &atomnum, atom, cinfo, minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rdivcrt(ifilename, &atomnum, atom, cinfo, minfo);
        }
        atomicnum_flag = 1;
        atomname_flag = 1;
        default_flag = 2;
        connect_flag = 1;
        bondtype_flag = 2;
    } else if (strcmp("divout", cinfo.intype) == 0 || strcmp("22", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag = rdivout(ifilename, &atomnum, atom, &cinfo, &minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rdivout(ifilename, &atomnum, atom, &cinfo, &minfo);
        }
        atomicnum_flag = 1;
        atomname_flag = 1;
        default_flag = 2;
        connect_flag = 1;
        bondtype_flag = 2;
    } else if (strcmp("sqmcrt", cinfo.intype) == 0 || strcmp("23", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag = rsqmcrt(ifilename, &atomnum, atom, cinfo, minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rsqmcrt(ifilename, &atomnum, atom, cinfo, minfo);
        }
        atomicnum_flag = 1;
        atomname_flag = 1;
        default_flag = 2;
        connect_flag = 1;
        bondtype_flag = 2;
    } else if (strcmp("sqmout", cinfo.intype) == 0 || strcmp("24", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag = rsqmout(ifilename, &atomnum, atom, &cinfo, &minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rsqmout(ifilename, &atomnum, atom, &cinfo, &minfo);
        }
        atomicnum_flag = 1;
        atomname_flag = 1;
        default_flag = 2;
        connect_flag = 1;
        bondtype_flag = 2;
    } else if (strcmp("charmm", cinfo.intype) == 0 || strcmp("25", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag =
            rcharmm(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag =
                rcharmm(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
        }
        default_flag = 2;
        atomicnum_flag = 1;
        adjustatomname_flag = 1;
        if (usr_aan_flag == 0)
            adjustatomname_flag = 0;
    } else if (strcmp("gesp", cinfo.intype) == 0 || strcmp("26", cinfo.intype) == 0) {
        check_input_file_format(ifilename, cinfo.intype);
        overflow_flag = rgesp(ifilename, &atomnum, atom, cinfo, minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rgesp(ifilename, &atomnum, atom, cinfo, minfo);
        }
        atomicnum_flag = 1;
        atomname_flag = 1;
        default_flag = 2;
        connect_flag = 1;
        bondtype_flag = 2;
    } else if (strcmp("pdbqt", cinfo.intype) == 0 || strcmp("31", cinfo.intype) == 0) {
        overflow_flag = rpdbqt(ifilename, &atomnum, atom, cinfo, minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rpdbqt(ifilename, &atomnum, atom, cinfo, minfo);
        }
        adjustatomname_flag = 1;
        if (usr_aan_flag == 0)
            adjustatomname_flag = 0;
        atomicnum_flag = 1;
        default_flag = 2;
        bondtype_flag = 2;
        connect_flag = 1;
    } else {
        eprintf("Unknown input file format (%s).", cinfo.intype);
    }

/*****************************************************************************/
/* End of big block over input file formats */
/*****************************************************************************/

}

#endif                          // FILEFORMAT_C__YES_THAT_IS_C
