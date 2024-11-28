/* PDB */
/*
1 - 6 Record name "ATOM "
7 - 11 Integer serial Atom serial number.
13 - 16 Atom name Atom name.
17 Character altLoc Alternate location indicator.
18 - 20 Residue name resName Residue name.
22 Character chainID Chain identifier.
23 - 26 Integer resSeq Residue sequence number.
27 AChar iCode Code for insertion of residues.
31 - 38 Real(8.3) x Orthogonal coordinates for X in Angstroms
39 - 46 Real(8.3) y Orthogonal coordinates for Y in Angstroms
47 - 54 Real(8.3) z Orthogonal coordinates for Z in Angstroms
55 - 60 Real(6.2) occupancy Occupancy.
61 - 66 Real(6.2) tempFactor Temperature factor.
77 - 78 LString(2) element Element symbol, right-justified.
79 - 80 LString(2) charge Charge on the atom.
12345678901234567890123456789012345678901234567890123456789012345678901234567890
ATOM    145  N   VAL A  25      32.433  16.336  57.540  1.00 11.92           N
ATOM    146  CA  VAL A  25      31.132  16.439  58.160  1.00 11.85           C  
*/

#include <ctype.h>
#include "eprintf.h"


int rpdbqt(char *filename, int *atomnum, ATOM * atom, CONTROLINFO cinfo, MOLINFO minfo)
{
    int numatom;
    int terindex;
    int overflow_flag = 0;
    int warning_flag = 1;
    int resno = -999;             /* used to check for a multiple residue PDB file */
    int resno_bak = -998;         /* used to check for a multiple residue PDB file */
    int resno_counter = 0;      /* used to check for a multiple residue PDB file */
    const int MULTIPLE_RESNO_WARNING_TOLERANCE = 10;
    char tmpchar[MAXCHAR];
    char line[MAXCHAR];
    char elem[MAXCHAR];
    char atomname[MAXCHAR];
    char resname[MAXCHAR];
    double tmpfloat1, tmpfloat2;
    double x, y, z, charge;
    int id;
    int tmpint;
    int i;
    int len;
    FILE *fpin;

    fpin = efopen(filename, "r");
    numatom = 0;
    initial(cinfo.maxatom, atom, minfo.resname);
    terindex = -1;
    for (;;) {
        if (fgets(line, MAXCHAR, fpin) == NULL) {
            /*  printf("\nFinished reading %s file.", filename); */
            break;
        }
        if (strncmp("TER", line, 3) == 0) {
            terindex = 1;
            continue;
        }
        if (strncmp("ATOM", line, 4) == 0) {
            if (overflow_flag == 0) {
                line[26] = ' ';
		sscanf(&line[11], "%s", atomname); 
                strcpy(atom[numatom].name, atomname);
		sscanf(&line[17], "%s", resname); 
		strcpy(atom[numatom].aa, resname);
                atom[numatom].chain[0] = line[21];
		if(line[25] == ' ') 
			resno = 1;
		else
			sscanf(&line[22], "%d", &resno);	

                atom[numatom].ter = terindex;
                atom[numatom].resno = resno;
                sscanf(&line[30], "%lf%lf%lf%lf%lf%lf", &x, &y, &z, &tmpfloat1, &tmpfloat2, &charge);
                atom[numatom].x = x;
                atom[numatom].y = y;
                atom[numatom].z = z;
                atom[numatom].charge = charge;

                if (atom[numatom].resno != resno_bak) {
                    resno_bak = atom[numatom].resno;
                    ++resno_counter;
                }
                if (terindex == 1) {
                    atom[numatom].ter = terindex;
                    terindex = -1;
                }
                if (strcmp(atom[numatom].name, "dumm") == 0)
                    continue;
                if (strcmp(atom[numatom].name, "Du") == 0)
                    continue;
                if (strcmp(atom[numatom].name, "DUMM") == 0)
                    continue;
            }
            numatom++;
            if (numatom >= cinfo.maxatom && overflow_flag == 0) {
                printf("Info: The number of atoms (%d) exceeded MAXATOM.\n", numatom);
                overflow_flag = 1;
            }
        }
    }
    *atomnum = numatom;
    fclose(fpin);
    return overflow_flag;

}

void wpdbqt()
{
    printf("Warning: Cannot write a pdbqt output file.\n"
           "         You must run a program provided by autodock developers.\n");
}

