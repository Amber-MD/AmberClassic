char *amberhome;
# define debug 0
# include "define.h"
# include "atom.h"
# include "eprintf.h"

# include "utility.c"
# include "common.c"
# include "rotate.c"
# include "ac.c"
# include "charmm.c"
# include "mol2.c"
# include "mmcif.c"
# include "mopcrt.c"
# include "divcrt.c"
# include "mopint.c"
# include "mopout.c"
# include "divout.c"
# include "sqmcrt.c"
# include "sqmout.c"
# include "gesp.c"
# include "gcrt.c"
# include "orca.c"
# include "gzmat.c"
# include "gout.c"
# include "orcout.c"
# include "pdb.c"
# include "csd.c"
# include "mdl.c"
# include "alc.c"
# include "hin.c"
# include "prep.c"
# include "rst.c"
# include "jzmat.c"
# include "jcrt.c"
# include "jout.c"
MOLINFO minfo2;

static int overflow_flag = 0;   /*if overflow_flag ==1, reallocate memory */
int atomtype_flag;              /*judge atom type? */
int bondtype_flag;              /*judge bond type? */
int default_flag = 0;           /*assign default information? */
int atomname_flag = 0;          /*assign atom name? */
int atomicnum_flag = 0;         /*judge atomic number according to atom name ? */
int adjustatomname_flag = 0;    /*adjust atom name? */
int usr_aan_flag = -99;         /*the user's input for adjusting atom names */
int duplicatedname_flag = 0;    /*check atom name duplication? */
int cartcoord_flag = 0;         /*generate coordinate from internal coordinate ? */
int connect_flag = 0;           /*judge atom connectivity and generate bond ? */

int atomnum = 0;                // 1 more than the highest valid atom array index.
int bondnum = 0;
int ringnum = 0;
ATOM *atom;                     // the atom array.
BOND *bond;
RING *ring;
AROM *arom;
MOLINFO minfo;
CONTROLINFO cinfo;

char line[MAXCHAR];
char ifilename[MAXCHAR];

# include "checkmolecule.c"
# include "fileformat.c"


void usage()
{
    if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
        printf("[31mUsage: acdoctor -i [0m input file name\n"
               "[31m                -f, -fi [0m input file format\n");
        printf("\n               [31m List of the File Formats [0m \n");
        printf
            ("\n	 	file format type  abbre. index | file format type abbre. index");
        printf
            ("\n		--------------------------------------------------------------- ");
        printf
            ("\n		Antechamber        ac       1  | Sybyl Mol2         mol2    2 ");
        printf
            ("\n		PDB                pdb      3  | Modified PDB       mpdb    4 ");
        printf
            ("\n		AMBER PREP (int)   prepi    5  | AMBER PREP (car)   prepc   6 ");
        printf
            ("\n		Gaussian Z-Matrix  gzmat    7  | Gaussian Cartesian gcrt    8 ");
        printf
            ("\n		Mopac Internal     mopint   9  | Mopac Cartesian    mopcrt 10 ");
        printf
            ("\n		Gaussian Output    gout    11  | Mopac Output       mopout 12 ");
        printf
            ("\n		Alchemy            alc     13  | CSD                csd    14 ");
        printf
            ("\n		MDL                mdl     15  | Hyper              hin    16 ");
        printf
            ("\n		AMBER Restart      rst     17  | Jaguar Cartesian   jcrt   18 ");
        printf
            ("\n		Jaguar Z-Matrix    jzmat   19  | Jaguar Output      jout   20 ");
        printf
            ("\n		Divcon Input       divcrt  21  | Divcon Output      divout 22 ");
        printf
            ("\n		SQM Input          sqmcrt  23  | SQM Output         sqmout 24 ");
        printf
            ("\n		Charmm             charmm  25  | Gaussian ESP       gesp   26 ");
        printf
            ("\n		geostd cif         ccif    27  |                              ");

        printf
            ("\n		--------------------------------------------------------------\n");
    } else {
        printf("Usage: acdoctor -i input file name\n"
               "                -f, -fi input file format\n");
        printf("\n                List of the File Formats \n");
        printf
            ("\n            file format type  abbre. index | file format type abbre. index");
        printf
            ("\n            --------------------------------------------------------------- ");
        printf
            ("\n            Antechamber        ac       1  | Sybyl Mol2         mol2    2 ");
        printf
            ("\n            PDB                pdb      3  | Modified PDB       mpdb    4 ");
        printf
            ("\n            AMBER PREP (int)   prepi    5  | AMBER PREP (car)   prepc   6 ");
        printf
            ("\n            Gaussian Z-Matrix  gzmat    7  | Gaussian Cartesian gcrt    8 ");
        printf
            ("\n            Mopac Internal     mopint   9  | Mopac Cartesian    mopcrt 10 ");
        printf
            ("\n            Gaussian Output    gout    11  | Mopac Output       mopout 12 ");
        printf
            ("\n            Alchemy            alc     13  | CSD                csd    14 ");
        printf
            ("\n            MDL                mdl     15  | Hyper              hin    16 ");
        printf
            ("\n            AMBER Restart      rst     17  | Jaguar Cartesian   jcrt   18 ");
        printf
            ("\n            Jaguar Z-Matrix    jzmat   19  | Jaguar Output      jout   20 ");
        printf
            ("\n            Divcon Input       divcrt  21  | Divcon Output      divout 22 ");
        printf
            ("\n	        SQM Input          sqmcrt  23  | SQM Output         sqmout 24 ");
        printf
            ("\n	        Charmm             charmm  25  | Gaussian ESP       gesp   26 ");
        printf
            ("\n	        geostd cif         ccif    27  |                              ");
        printf
            ("\n            --------------------------------------------------------------\n");
    }
}


void judgebondtype(int atomnum, ATOM * atom, int bondnum, BOND * bond, CONTROLINFO cinfo,
                   MOLINFO minfo, int bondtype_flag)
{
    size_t copied_size;
    char *options;
    char tmpchar[MAXCHAR];
    wac("ACDOCTOR_INIT.ac", atomnum, atom, bondnum, bond, cinfo, minfo);

    copied_size = build_exe_path(tmpchar, "bondtype", sizeof tmpchar, 1);
    if (bondtype_flag == 1)
        options = " -j part -i ACDOCTOR_INIT.ac" " -o ACDOCTOR_BOND.ac -f ac";
    else
        options = " -j full -i ACDOCTOR_INIT.ac" " -o ACDOCTOR_BOND.ac -f ac";
    strncat(tmpchar, options, MAXCHAR - copied_size);

    fprintf(stdout, "Running: %s\n", tmpchar);
    esystem(tmpchar);
    minfo2 = minfo;
    rac("ACDOCTOR_BOND.ac", &atomnum, atom, &bondnum, bond, &cinfo, &minfo2);
}


void trybondandatomtypes(void)
{
    char tmpchar[MAXCHAR];

    fprintf(stdout, "-- Now try to judge bond type -- \n");
    if (bondnum > 0)
        judgebondtype(atomnum, atom, bondnum, bond, cinfo, minfo, bondtype_flag);

    fprintf(stdout, "-- Now try to assign atom type -- \n");
    build_exe_path(tmpchar, "atomtype", sizeof tmpchar, 1);
    strcat(tmpchar, " -i ACDOCTOR_BOND.ac -o ACDOCTOR_ATOM.ac -p ");
    strcat(tmpchar, minfo.atom_type_def);
    fprintf(stdout, "Running: %s\n", tmpchar);
    esystem(tmpchar);

    fprintf(stdout, "-- Now write out ACDOCTOR.mol2 -- \n");
    build_exe_path(tmpchar, "atomtype", sizeof tmpchar, 1);
    strcat(tmpchar, " -i ACDOCTOR_BOND.ac -o ACDOCTOR_SYBYL.ac -p sybyl");
    esystem(tmpchar);
    minfo2 = minfo;
    rac("ACDOCTOR_SYBYL.ac", &atomnum, atom, &bondnum, bond, &cinfo, &minfo2);
    wmol2("ACDOCTOR.mol2", atomnum, atom, bondnum, bond, arom, cinfo, minfo);

    fprintf(stdout,
            "\n-- IF no error occurs in bond and atom type assignments,"
            "\n   %s should present no problems to the antechamber package. -- \n"
            "\nACDOCTOR_INIT.ac  : ac file before bond type and atom type assignment"
            "\nACDOCTOR_BT.ac    : ac file after bond type assignment"
            "\nACDOCTOR_AT.ac    : ac file after gaff atom type assignment"
            "\nACDOCTOR_SYBYL.ac : ac file after sybyl atom type assignment"
            "\nACDOCTOR.mol2     : mol2 file after sybyl atom type assignment\n",
            ifilename);
}

int main(int argc, char *argv[])
{
    int i, j, k;
    int index;

    fprintf(stdout, "\nWelcome to acdoctor %s: check and diagnose problems"
            " in molecular input files.\n\n", ANTECHAMBER_VERSION);
    esetprogramname(argv[0]);
    amberhome = (char *) getenv("AMBERCLASSICHOME");
    if (amberhome == NULL) {
        eprintf("AMBERCLASSICHOME is not set.");
    }
    if (argc == 2)
        if (strncmp(argv[1], "-h", 2) == 0 || strncmp(argv[1], "-H", 2) == 0) {
            usage();
            exit(0);
        }
    if (argc == 1) {
        usage();
        exit(1);
    }

/* 	set defaults information */
    default_cinfo(&cinfo);
    default_minfo(&minfo);
    atomtype_flag = 0;
    bondtype_flag = 0;
/* 	
	0 <-> not set 
	1 <-> set atom[].chain, atom[].ter and atom[].ambername	
	2 <-> set atom[].chain, atom[].ter
*/
    default_flag = 0;
    atomname_flag = 0;
    atomicnum_flag = 0;
    adjustatomname_flag = 0;
    duplicatedname_flag = 1;
    cartcoord_flag = 0;
    connect_flag = 0;


    for (i = 1; i < argc - 1; i += 2) {
        if (strcmp(argv[i], "-i") == 0) {
            strcpy(ifilename, argv[i + 1]);
            continue;
        } else if (strcmp(argv[i], "-f") == 0) {
            strcpy(cinfo.intype, argv[i + 1]);
            continue;
        } else if (strcmp(argv[i], "-fi") == 0) {
            strcpy(cinfo.intype, argv[i + 1]);
            continue;
        }
    }

/* 	for connect.tpl and radius parameter files */
    build_dat_path(minfo.connect_file, "CONNECT.TPL", sizeof minfo.connect_file, 0);
    build_dat_path(minfo.radius_file, "RADIUS.DAT", sizeof minfo.radius_file, 0);


/*      allocate memory using default parameters MAXATOM and MAXBOND */
    memory(0, MAXATOM, MAXBOND, MAXRING);

/******************************************/
/* 	The following codes readin input file */
/******************************************/

    if (strcmp("ac", cinfo.intype) == 0 || strcmp("1", cinfo.intype) == 0) {
        check_input_file_format(ifilename, "ac");
        overflow_flag = rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag =
                rac(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
        }
        adjustatomname_flag = 1;
        atomicnum_flag = 1;
        checkbyatomtype = 1;
        checkbybondtype = 1;
    }

    if (strcmp("mol2", cinfo.intype) == 0 || strcmp("2", cinfo.intype) == 0) {

        check_input_file_format(ifilename, "mol2");
        overflow_flag =
            rmol2(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo, 1);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag =
                rmol2(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo, 1);
        }
        default_flag = 2;
        atomicnum_flag = 1;
        adjustatomname_flag = 1;
        checkbyatomtype = 1;
        checkbybondtype = 1;
    }

    if (strcmp("ccif", cinfo.intype) == 0 || strcmp("27", cinfo.intype) == 0) {

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
        atomtype_flag = 1;
        bondtype_flag = 2;
        connect_flag = 1;
    }

    if (strcmp("mopint", cinfo.intype) == 0 || strcmp("9", cinfo.intype) == 0) {

        check_input_file_format(ifilename, "mopint");
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
        default_flag = 1;
        connect_flag = 1;
        bondtype_flag = 2;
    }

    if (strcmp("mopcrt", cinfo.intype) == 0 || strcmp("10", cinfo.intype) == 0) {

        check_input_file_format(ifilename, "mopcrt");
        overflow_flag = rmopcrt(ifilename, &atomnum, atom, cinfo, minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rmopcrt(ifilename, &atomnum, atom, cinfo, minfo);
        }
        atomicnum_flag = 1;
        atomname_flag = 1;
        default_flag = 1;
        connect_flag = 1;
        bondtype_flag = 2;
    }

    if (strcmp("mopout", cinfo.intype) == 0 || strcmp("12", cinfo.intype) == 0) {

        check_input_file_format(ifilename, "mopout");
        overflow_flag = rmopout(ifilename, &atomnum, atom, &cinfo, &minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rmopout(ifilename, &atomnum, atom, &cinfo, &minfo);
        }
        atomicnum_flag = 1;
        atomname_flag = 1;
        default_flag = 1;
        connect_flag = 1;
        bondtype_flag = 2;
    }

    if (strcmp("orcinp", cinfo.intype) == 0 || strcmp("29", cinfo.intype) == 0) {

        check_input_file_format(ifilename, "orcinp");
        overflow_flag = rorca(ifilename, &atomnum, atom, cinfo, minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rorca(ifilename, &atomnum, atom, cinfo, minfo);
        }
        atomname_flag = 1;
        default_flag = 1;
        connect_flag = 1;
        bondtype_flag = 2;
    }

    if (strcmp("gcrt", cinfo.intype) == 0 || strcmp("8", cinfo.intype) == 0) {

        check_input_file_format(ifilename, "gcrt");
        overflow_flag = rgcrt(ifilename, &atomnum, atom, cinfo, minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rgcrt(ifilename, &atomnum, atom, cinfo, minfo);
        }
        atomname_flag = 1;
        default_flag = 1;
        connect_flag = 1;
        bondtype_flag = 2;
    }

    if (strcmp("gzmat", cinfo.intype) == 0 || strcmp("7", cinfo.intype) == 0) {

        check_input_file_format(ifilename, "gzmat");
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
        default_flag = 1;
        connect_flag = 1;
        bondtype_flag = 2;
    }

    if (strcmp("orcout", cinfo.intype) == 0 || strcmp("30", cinfo.intype) == 0) {

        check_input_file_format(ifilename, "orcout");
        overflow_flag = rorcout(ifilename, &atomnum, atom, cinfo, &minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rorcout(ifilename, &atomnum, atom, cinfo, &minfo);
        }
        atomname_flag = 1;
        default_flag = 1;
        connect_flag = 1;
        bondtype_flag = 2;
    }


    if (strcmp("gout", cinfo.intype) == 0 || strcmp("11", cinfo.intype) == 0) {

        check_input_file_format(ifilename, "gout");
        overflow_flag = rgout(ifilename, &atomnum, atom, cinfo, &minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rgout(ifilename, &atomnum, atom, cinfo, &minfo);
        }
        atomname_flag = 1;
        default_flag = 1;
        connect_flag = 1;
        bondtype_flag = 2;
    }

    if (strcmp("jcrt", cinfo.intype) == 0 || strcmp("18", cinfo.intype) == 0) {

        check_input_file_format(ifilename, "jcrt");
        overflow_flag = rjcrt(ifilename, &atomnum, atom, cinfo, minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rjcrt(ifilename, &atomnum, atom, cinfo, minfo);
        }
        default_flag = 1;
        connect_flag = 1;
        bondtype_flag = 2;
    }

    if (strcmp("jzmat", cinfo.intype) == 0 || strcmp("19", cinfo.intype) == 0) {

        check_input_file_format(ifilename, "jzmat");
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
        default_flag = 1;
        connect_flag = 1;
        bondtype_flag = 2;
    }

    if (strcmp("jout", cinfo.intype) == 0 || strcmp("20", cinfo.intype) == 0) {

        check_input_file_format(ifilename, "jout");
        overflow_flag = rjout(ifilename, &atomnum, atom, cinfo, &minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rjout(ifilename, &atomnum, atom, cinfo, &minfo);
        }
        default_flag = 1;
        connect_flag = 1;
        bondtype_flag = 2;
    }

    if (strcmp("pdb", cinfo.intype) == 0 || strcmp("3", cinfo.intype) == 0) {
        check_input_file_format(ifilename, "pdb");
        overflow_flag = rpdb(ifilename, &atomnum, atom, cinfo, minfo, 0);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rpdb(ifilename, &atomnum, atom, cinfo, minfo, 0);
        }
        adjustatomname_flag = 1;
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
    }

    if (strcmp("mpdb", cinfo.intype) == 0 || strcmp("4", cinfo.intype) == 0) {

        check_input_file_format(ifilename, "mpdb");
        overflow_flag = rpdb(ifilename, &atomnum, atom, cinfo, minfo, 1);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rpdb(ifilename, &atomnum, atom, cinfo, minfo, 1);
        }
        adjustatomname_flag = 1;
        atomicnum_flag = 1;
        default_flag = 0;
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
    }
    if (strcmp("csd", cinfo.intype) == 0 || strcmp("14", cinfo.intype) == 0) {

        check_input_file_format(ifilename, "csd");
        overflow_flag = rcsd(ifilename, &atomnum, atom, cinfo, minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rcsd(ifilename, &atomnum, atom, cinfo, minfo);
        }
        atomicnum_flag = 1;
        default_flag = 1;
        connect_flag = 1;
        bondtype_flag = 2;
    }

    if (strcmp("mdl", cinfo.intype) == 0 || strcmp("sd", cinfo.intype) == 0
        || strcmp("sdf", cinfo.intype) == 0 || strcmp("15", cinfo.intype) == 0) {

        check_input_file_format(ifilename, "mdl");
        overflow_flag = rmdl(ifilename, &atomnum, atom, &bondnum, bond, cinfo, minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rmdl(ifilename, &atomnum, atom, &bondnum, bond, cinfo, minfo);
        }
        atomicnum_flag = 1;
        atomname_flag = 1;
        default_flag = 1;
        bondtype_flag = 1;
        checkbybondtype = 1;
    }

    if (strcmp("alc", cinfo.intype) == 0 || strcmp("13", cinfo.intype) == 0) {

        check_input_file_format(ifilename, "alc");
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
        default_flag = 1;
        checkbybondtype = 1;
    }

    if (strcmp("hin", cinfo.intype) == 0 || strcmp("16", cinfo.intype) == 0) {

        check_input_file_format(ifilename, "hin");
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
        default_flag = 1;
        checkbybondtype = 1;
    }

    if (strcmp("prepi", cinfo.intype) == 0 || strcmp("5", cinfo.intype) == 0) {

        check_input_file_format(ifilename, "prepi");
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
    }

    if (strcmp("prepc", cinfo.intype) == 0 || strcmp("6", cinfo.intype) == 0) {

        check_input_file_format(ifilename, "prepc");
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
    }

    if (strcmp("rst", cinfo.intype) == 0 || strcmp("17", cinfo.intype) == 0) {
        eprintf("RST (17) file format can only be an additional file\n"
                "because it only has coordinate information.");
    }

    if (strcmp("divcrt", cinfo.intype) == 0 || strcmp("21", cinfo.intype) == 0) {

        check_input_file_format(ifilename, "divcrt");
        overflow_flag = rdivcrt(ifilename, &atomnum, atom, cinfo, minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rdivcrt(ifilename, &atomnum, atom, cinfo, minfo);
        }
        atomicnum_flag = 1;
        atomname_flag = 1;
        default_flag = 1;
        connect_flag = 1;
        bondtype_flag = 2;
    }

    if (strcmp("divout", cinfo.intype) == 0 || strcmp("22", cinfo.intype) == 0) {

        check_input_file_format(ifilename, "divout");
        overflow_flag = rdivout(ifilename, &atomnum, atom, &cinfo, &minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rdivout(ifilename, &atomnum, atom, &cinfo, &minfo);
        }
        atomicnum_flag = 1;
        atomname_flag = 1;
        default_flag = 1;
        connect_flag = 1;
        bondtype_flag = 2;
    }

    if (strcmp("sqmcrt", cinfo.intype) == 0 || strcmp("23", cinfo.intype) == 0) {

        check_input_file_format(ifilename, "sqmcrt");
        overflow_flag = rsqmcrt(ifilename, &atomnum, atom, cinfo, minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rsqmcrt(ifilename, &atomnum, atom, cinfo, minfo);
        }
        atomicnum_flag = 1;
        atomname_flag = 1;
        default_flag = 1;
        connect_flag = 1;
        bondtype_flag = 2;
    }

    if (strcmp("sqmout", cinfo.intype) == 0 || strcmp("24", cinfo.intype) == 0) {

        check_input_file_format(ifilename, "sqmout");
        overflow_flag = rsqmout(ifilename, &atomnum, atom, &cinfo, &minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rsqmout(ifilename, &atomnum, atom, &cinfo, &minfo);
        }
        atomicnum_flag = 1;
        atomname_flag = 1;
        default_flag = 1;
        connect_flag = 1;
        bondtype_flag = 2;
    }

    if (strcmp("charmm", cinfo.intype) == 0 || strcmp("25", cinfo.intype) == 0) {

        overflow_flag =
            rcharmm(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag =
                rcharmm(ifilename, &atomnum, atom, &bondnum, bond, &cinfo, &minfo);
        }
        default_flag = 1;
        atomicnum_flag = 1;
        adjustatomname_flag = 1;
    }

    if (strcmp("gesp", cinfo.intype) == 0 || strcmp("26", cinfo.intype) == 0) {

        overflow_flag = rgesp(ifilename, &atomnum, atom, cinfo, minfo);
        if (overflow_flag) {
            cinfo.maxatom = atomnum + 10;
            cinfo.maxbond = bondnum + 10;
            memory(7, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag = rgesp(ifilename, &atomnum, atom, cinfo, minfo);
        }
        atomicnum_flag = 1;
        atomname_flag = 1;
        default_flag = 1;
        connect_flag = 1;
        bondtype_flag = 2;
    }

    if (adjustatomname_flag) {
        if (strcmp(cinfo.intype, "mol2") == 0 || strcmp(cinfo.intype, "2") == 0
            || strcmp(cinfo.intype, "ac") == 0 || strcmp(cinfo.intype, "1") == 0)
            adjustatomname(atomnum, atom, 1);
        else
            adjustatomname(atomnum, atom, 0);
    }
    if (atomicnum_flag) {
        if (cinfo.intstatus == 2)
            printf("Info: Determining atomic numbers from atomic symbols which "
                   "are case sensitive.\n");
        atomicnum(atomnum, atom);
    }
    if (atomname_flag)
        atomname(atomnum, atom);
    if (default_flag)
        default_inf(atomnum, atom, default_flag);
    if (cartcoord_flag)
        cartcoord(atomnum, atom);
    if (connect_flag) {
        overflow_flag =
            connect(minfo.connect_file, atomnum, atom, &bondnum, bond, cinfo.maxbond);
        if (overflow_flag) {
            cinfo.maxbond = bondnum + 10;
            memory(2, cinfo.maxatom, cinfo.maxbond, cinfo.maxring);
            overflow_flag =
                connect(minfo.connect_file, atomnum, atom, &bondnum, bond, cinfo.maxbond);
        }
    }
    if (duplicatedname_flag)
        duplicatedname(atomnum, atom);

/* 	print out information */
    if (debug) {
        for (i = 0; i < atomnum; i++)
            printf
                ("ATOM\t%d\t%s\t%d\t%s\t%8.3lf\t%8.3lf\t%8.3lf\t%8.3lf\t%d\t%d\t%d\t%d\t%d\t%d\n",
                 i + 1, atom[i].name, atom[i].resno, atom[i].aa, atom[i].x, atom[i].y,
                 atom[i].z, atom[i].charge, atom[i].con[0] + 1, atom[i].con[1] + 1,
                 atom[i].con[2] + 1, atom[i].con[3] + 1, atom[i].con[4] + 1,
                 atom[i].con[5] + 1);
        for (i = 0; i < bondnum; i++)
            printf("BOND\t%d\t%d\t%s\t%d\t%s\n", i + 1, bond[i].bondi + 1,
                   atom[bond[i].bondi].name, bond[i].bondj + 1, atom[bond[i].bondj].name);
    }

    check_input_molecule();
    printf("\n");
    trybondandatomtypes();

    printf("\n");
    return (0);
}
