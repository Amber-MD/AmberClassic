// Molecule input structure validation.

// This file doubles as a header file and as an included (!) source file.
// Yes, we are on the road to multiple header files and separate compilation.

#ifndef CHECKMOLECULE_H
#define CHECKMOLECULE_H

extern int checkmolecule;       // check molecule in the manner of acdoctor
extern int checkbyatomtype;     // atom type based open valence checking
extern int checkbybondtype;     // bond type based weird bond checking

#endif                          // CHECKMOLECULE_H


#ifndef CHECKMOLECULE_C__YES_THAT_IS_C
#define CHECKMOLECULE_C__YES_THAT_IS_C

# include <assert.h>
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include "define.h"
# include "atom.h"
# include "eprintf.h"
# define MAXCORR 256
# define debug 0

int checkmolecule = 0;          // check molecule in the manner of acdoctor
int checkbyatomtype = 0;        // atom type based open valence checking
int checkbybondtype = 0;        // bond type based weird bond checking

typedef struct {
    char type[10];
    char element[10];
    int connum;
} CORR;

typedef struct {
    double low;
    double high;
    double bond;
    double offset;
    int valid;
} BONDLENGTH;

double bondradius[120];
BONDLENGTH bl;


void bondrange(int atomicnum1, int atomicnum2)
{
    double bond;
    double offset;
    bl.low = 0;
    bl.high = 0;
    if (bondradius[atomicnum1] > 0 && bondradius[atomicnum2] > 0) {
        bond = bondradius[atomicnum1] + bondradius[atomicnum2];
        if (bond <= 1.5)
            offset = bond * 0.15;
        if (bond > 1.5 && bond <= 1.90)
            offset = bond * 0.11;
        if (bond > 1.90 && bond <= 2.05)
            offset = bond * 0.09;
        if (bond > 2.05)
            offset = bond * 0.08;
        bl.offset = offset;
        bl.bond = bond;
/* definition in connect() */
/*
	bl.low = bond * 0.5;
	bl.high = bond + offset;
*/
/* now apply more stringent definition */
        bl.low = bond * 0.65;
        bl.high = bond + offset * 0.75;
        bl.valid = 1;
    } else
        bl.valid = 0;
    fflush(stdout);
}


void check_input_molecule()
{
    extern int atomnum;         // 1 more than the highest valid atom array index.
    extern int bondnum;
    extern ATOM *atom;          // the atom array.
    extern BOND *bond;
    extern MOLINFO minfo;
    int i, j, k;
    int numh = 0;
    int numn = 0;
    int numo = 0;
    int numc = 0;
    int nums = 0;
    int nump = 0;
    int numx = 0;               /* number of halogens */
    FILE *fp;
    char filename[MAXCHAR];
    char line[MAXCHAR];
    CORR *corr;
    int corrnum = 0;
    int tvalence;
    int val_ar;
    int connum;
    int tmpint;
    int selectflag;             /*used by unit check */
    int errorflag;
    int warnings;
    int numselect;
    int bondflag;
    double dist;
    double radius;

/*part1: element */
    printf("-- Check Unusual Elements --\n");
    for (i = 0; i < atomnum; i++) {
        if (atom[i].atomicnum == 1)
            numh++;
        else if (atom[i].atomicnum == 9 || atom[i].atomicnum == 17
                 || atom[i].atomicnum == 35 || atom[i].atomicnum == 53)
            numx++;
        else if (atom[i].atomicnum == 6)
            numc++;
        else if (atom[i].atomicnum == 7)
            numn++;
        else if (atom[i].atomicnum == 8)
            numo++;
        else if (atom[i].atomicnum == 15)
            nump++;
        else if (atom[i].atomicnum == 16)
            nums++;
        else
            fprintf(stdout, "Warning: Unusual element (%s)"
                            " for atom (ID: %d, Name: %s).\n",
                    atom[i].element, i + 1, atom[i].name);
    }
    if ((atomnum - numh - numn - numo - numc - nums - nump - numx) > 0) {
        eprintf("GAFF does not have sufficient parameters for molecules having unusual\n"
                "       elements (those other than H,C,N,O,S,P and halogens).\n"
                "       To ensure antechamber works properly, one may need to designate\n"
                "       bond types for bonds involved with unusual elements.\n"
                "       To do so, simply freeze the bond types by appending \"F\" or \"f\" \n"
                "       to the corresponding bond types in ac or mol2 files\n"
                "       and rerun antechamber without unusual element checking via:\n"
                "       antechamber -dr no \n"
                "       Alternatively for metals, see metalpdb2mol2.py in MCPB. \n");
    } else if (atomnum == numh + numn + numo + numc + nums + nump + numx) {
        printf("   Status: pass\n");
    }

/*part2: unfilled valences */
    printf("-- Check Open Valences --\n");
    errorflag = 0;
    warnings = 0;
// how many hydrogen and halogen atoms.
    if (numh + numx == 0) {
        fprintf(stdout,
                "Warning: This molecule has no hydrogens nor halogens.\n"
                "         It is quite possible that there are unfilled valences.\n");
        ++warnings;
    }

/* check if atom[].connum is consistent with that in CORR_NAME_TYPE.DAT */
    if (checkbyatomtype == 1) {
        corr = (CORR *) emalloc(sizeof(CORR) * MAXCORR);
        build_dat_path(filename, "CORR_NAME_TYPE.DAT", sizeof filename, 0);
        if ((fp = fopen(filename, "r")) == NULL) {
            eprintf("Cannot open the CORR_NAME_TYPE.DAT file (%s).",
                    filename);
        }
        for (;;) {
            if (fgets(line, MAXCHAR, fp) == NULL)
                break;
            sscanf(line, "%s%s%d", corr[corrnum].type, corr[corrnum].element,
                   &corr[corrnum].connum);
            corrnum++;
            if (corrnum == MAXCORR) {
                corr = (CORR *) erealloc(corr, sizeof(CORR) * (MAXCORR + corrnum));
            }
        }
        fclose(fp);

        if (debug)
            for (i = 0; i < atomnum; i++)
                printf("\nATOM %d %s %8.3lf %8.3lf %8.3lf %8.3lf %s\n", i + 1,
                       atom[i].name, atom[i].x, atom[i].y, atom[i].z, atom[i].charge,
                       atom[i].ambername);
        for (i = 0; i < atomnum; i++)
            for (j = 0; j < corrnum; j++)
                if (strcmp(atom[i].ambername, corr[j].type) == 0)
                    if (corr[j].connum != -1 && corr[j].connum != atom[i].connum) {
                        fprintf(stdout,
                                "Warning: The number of bonds (%d) for"
                                " atom (ID: %d, Name: %s) does not match\n",
                                atom[i].connum, i + 1, atom[i].name);
                        fprintf(stdout,
                                "         the connectivity (%d) for"
                                " atom type (%s) defined in CORR_NAME_TYPE.DAT.\n",
                                corr[j].connum, corr[j].type);
                        errorflag = 1;
/*		        exit(1); */
// Although error happens, it is ideal not to exit since many atom names (N2)
// same to AMBER atom types, such as N2.
                    }
        if (errorflag == 1)
                // coloron = "[32m";
                // coloroff = "[0m";
            fprintf(stdout,
                    "But, you may safely ignore the warnings if your molecule\n"
                    "         uses atom names or element names as atom types.\n");
        free(corr);
    }
    if (errorflag == 0 && warnings == 0) {
        printf("   Status: pass\n");
    }

/*part3: geometry */
    printf("-- Check Geometry --\n");
    errorflag = 0;
    warnings = 0;
    printf("      for those bonded   \n");

    if ((fp = fopen(minfo.connect_file, "r")) == NULL) {
        eprintf("Cannot open the minfo.connect_file (%s).", minfo.connect_file);
    }
    for (i = 0; i < 120; i++)
        bondradius[i] = -1;
    for (;;) {
        if (fgets(line, MAXCHAR, fp) == NULL)
            break;
        if (line[10] == '.') {
            sscanf(&line[3], "%d%lf", &tmpint, &radius);
            bondradius[tmpint] = radius;
        }
    }

    for (i = 0; i < bondnum; i++) {
        bondrange(atom[bond[i].bondi].atomicnum, atom[bond[i].bondj].atomicnum);
        if (bl.valid == 0)
            continue;
        dist =  (atom[bond[i].bondi].x - atom[bond[i].bondj].x) *
                (atom[bond[i].bondi].x - atom[bond[i].bondj].x);
        dist += (atom[bond[i].bondi].y - atom[bond[i].bondj].y) *
                (atom[bond[i].bondi].y - atom[bond[i].bondj].y);
        dist += (atom[bond[i].bondi].z - atom[bond[i].bondj].z) *
                (atom[bond[i].bondi].z - atom[bond[i].bondj].z);
        dist = sqrt(dist);
        if (dist < bl.low) {
            fprintf(stdout,
                    "Warning: Small distance for BOND\t%d\t%s\t%s\t%d\t%9.2lf  [%-5.2lf-%5.2lf]\n",
                    i + 1, atom[bond[i].bondi].name, atom[bond[i].bondj].name,
                    bond[i].type, dist, bl.low, bl.high);
            ++warnings;
        }
        if (dist > bl.high) {
            fprintf(stdout,
                    "Warning: Large distance for BOND\t%d\t%s\t%s\t%d\t%9.2lf  [%-5.2lf-%5.2lf]\n",
                    i + 1, atom[bond[i].bondi].name, atom[bond[i].bondj].name,
                    bond[i].type, dist, bl.low, bl.high);
            ++warnings;
        }
        if (debug)
            fprintf(stdout, "[33m\nBOND\t%d\t%s\t%s\t%d\t%9.2lf[0m", i + 1,
                    atom[bond[i].bondi].name, atom[bond[i].bondj].name, bond[i].type,
                    dist);
    }

/* now check if those unbonded atoms come too close */
    printf("      for those not bonded   \n");
    for (i = 0; i < atomnum - 1; i++)
        for (j = i + 1; j < atomnum; j++) {
            bondflag = 0;
            for (k = 0; k < bondnum; k++) {
                if (bond[k].bondi == i && bond[k].bondj == j) {
                    bondflag = 1;
                    break;
                }
                if (bond[k].bondi == j && bond[k].bondj == i) {
                    bondflag = 1;
                    break;
                }
            }
            if (bondflag == 1)
                continue;
            bondrange(atom[i].atomicnum, atom[j].atomicnum);
            if (bl.valid == 0)
                continue;
            dist =  (atom[i].x - atom[j].x) * (atom[i].x - atom[j].x);
            dist += (atom[i].y - atom[j].y) * (atom[i].y - atom[j].y);
            dist += (atom[i].z - atom[j].z) * (atom[i].z - atom[j].z);
            dist = sqrt(dist);
            if (dist < (bl.bond + bl.offset * 1.25)) {
                // coloron = "[33m";
                // coloroff = "[0m";
                fprintf(stdout,
                        "Warning: Close (%5.2lf) nonbonded atoms"
                        " (ID: %d, Name: %s) and (ID: %d, Name: %s).\n",
                        dist, i + 1, atom[i].name, j + 1, atom[j].name);
                ++warnings;
            }
        }
    if (errorflag == 0 && warnings == 0) {
        printf("   Status: pass\n");
    }

/*part4  now check bonds */
    printf("-- Check Weird Bonds --\n");
    if (checkbybondtype == 1) {
        if (debug)
            for (i = 0; i < bondnum; i++)
                printf("\nBOND %d  %d  %s  %d  %s  %d\n", i + 1, bond[i].bondi + 1,
                       atom[bond[i].bondi].name, bond[i].bondj + 1,
                       atom[bond[i].bondj].name, bond[i].type);
        for (i = 0; i < atomnum; i++) {
            tvalence = 0;
            val_ar = 1;
            for (j = 0; j < bondnum; j++)
                if (bond[j].bondi == i || bond[j].bondj == i) {
                    if (bond[j].type <= 3)
                        tvalence += bond[j].type;
                    else if (bond[j].type == 7)
                        tvalence += 1;
                    else if (bond[j].type == 8)
                        tvalence += 2;
                    else if (bond[j].type == 9 || bond[j].type == 10) {
                        tvalence += val_ar;
                        if (val_ar == 1)
                            val_ar = 2;
                        else
                            val_ar = 1;
                    }
                }
            if (atom[i].atomicnum == 6) {
                if (4 < tvalence || tvalence < 4) {
                    eprintf("Weird atomic valence (%d)"
                            " for atom (ID: %d, Name: %s).\n"
                            "       Possible open valence.",
                            tvalence, i + 1, atom[i].name);
                }
            }

            if (atom[i].atomicnum == 1 || atom[i].atomicnum == 9
                || atom[i].atomicnum == 17 || atom[i].atomicnum == 35
                || atom[i].atomicnum == 53)
                if (tvalence != 1) {
                    eprintf("Weird atomic valence (%d)"
                            " for atom (ID: %d, Name: %s).\n"
                            "       Please check atomic connectivity.",
                            tvalence, i + 1, atom[i].name);
                }

            if (atom[i].atomicnum == 7) {
                if (tvalence > 4) {
                    eprintf("Weird atomic valence (%d)"
                            " for atom (ID: %d, Name: %s).\n"
                            "       Please check atomic connectivity.",
                            tvalence, i + 1, atom[i].name);
                }
                if (tvalence < 3) {
                    eprintf("Weird atomic valence (%d)"
                            " for atom (ID: %d, Name: %s).\n"
                            "       Possible open valence.",
                            tvalence, i + 1, atom[i].name);
                }
            }
            if (atom[i].atomicnum == 6) {
                if (tvalence > 4) {
                    eprintf("Weird atomic valence (%d)"
                            " for atom (ID: %d, Name: %s).\n"
                            "       Please check atomic connectivity.",
                            tvalence, i + 1, atom[i].name);
                }
                if (tvalence <= 3) {
                    eprintf("Weird atomic valence (%d)"
                            " for atom (ID: %d, Name: %s).\n"
                            "       Possible open valence.",
                            tvalence, i + 1, atom[i].name);
                }
            }
            if (atom[i].atomicnum == 8) {
                if (tvalence > 2) {
                    eprintf("Weird atomic valence (%d)"
                            " for atom (ID: %d, Name: %s).\n"
                            "       Please check atomic connectivity.",
                            tvalence, i + 1, atom[i].name);
                }
                if (tvalence < 1) {
                    eprintf("Weird atomic valence (%d)"
                            " for atom (ID: %d, Name: %s).\n"
                            "       Possible open valence.",
                            tvalence, i + 1, atom[i].name);
                }
            }
        }
    }
/* check if bonded atoms exceeds 6 */
    if (checkbybondtype == 0)
        for (i = 0; i < atomnum; i++) {
            connum = 0;
            for (j = 0; j < bondnum; j++) {
                if (bond[j].bondi == i || bond[j].bondj == i)
                    connum++;
                if (connum > 6)
                    break;
            }
            if (connum > 6) {
                eprintf("Number of bonded atoms (%d) for atom"
                        " (ID: %d, Name: %s) exceeds 6.\n"
                        "       antechamber cannot handle such kinds of molecules.",
                        connum, i + 1, atom[i].name);
            }
        }

    printf("   Status: pass\n");

/*part5: check if all atoms are linked together through paths*/
    printf("-- Check Number of Units --\n");
    for (i = 0; i < atomnum; i++)
        atom[i].select = 0;
    atom[0].select = 1;
    selectflag = 1;
    while (selectflag) {
        selectflag = 0;
        for (i = 0; i < atomnum; i++) {
            if (atom[i].select == 0)
                continue;
            for (j = 0; j < atom[i].connum; j++)
                if (atom[atom[i].con[j]].select == 0) {
                    atom[atom[i].con[j]].select = 1;
                    selectflag = 1;
                }
        }
    }
    numselect = 0;
    for (i = 0; i < atomnum; i++)
        numselect += atom[i].select;
    if (numselect < atomnum) {
        // coloron = "[31m";
        // coloroff = "[0m";
        eprintf("This molecule may have more than one unit.\n"
                "       antechamber can only handle one unit.  If the input is a single unit\n"
                "       then the connectivity is wrong and the geometry may be bad.\n"
                "       Please convert your molecule to a mol2 file via:\n"
                "       antechamber -j 5 -at sybyl -dr no \n"
                "       And then check your molecule with a visualization program;\n"
                "       manually add missing bonds or delete unwanted bonds as appropriate.");
    }

    printf("   Status: pass\n");
}


#endif                          // CHECKMOLECULE_C__YES_THAT_IS_C
