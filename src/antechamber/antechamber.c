char *amberhome;
# include <assert.h>
# include <ctype.h>
# include <errno.h>
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "define.h"
# include "atom.h"
# include "eprintf.h"

# include "utility.c"
# include "common.c"
# include "equatom.c"
# include "ac.c"
# include "charmm.c"
# include "mol2.c"
# include "mmcif.c"
# include "mopcrt.c"
# include "divcrt.c"
# include "sqmcrt.c"
# include "sqmout.c"
# include "mopint.c"
# include "mopout.c"
# include "divout.c"
# include "gcrt.c"
# include "orca.c"
# include "gzmat.c"
# include "gout.c"
# include "orcout.c"
# include "gamess.c"
# include "pdb.c"
# include "pdbqt.c"
# include "csd.c"
# include "mdl.c"
# include "alc.c"
# include "hin.c"
# include "prep.c"
# include "rst.c"
# include "jzmat.c"
# include "jcrt.c"
# include "jout.c"
# include "gesp.c"
# include "charge.c"
# include "esp_quick.c"

int ao_flag = 0;
/*For addtional file 
 1, only readin coordinates
 2, only readin charge
 3, only readin atom names
 4, only readin atom types
 5, only readin bond types 
*/
int atomtype_flag = 0;          /*judge atom type? */
int bondtype_flag = 0;          /*judge bond type? */
int default_flag = 0;           /*assign default information? */
int atomname_flag = 0;          /*assign atom name? */
int atomicnum_flag = 0;         /*judge atomic number according to atom name ? */
int adjustatomname_flag = 0;    /*adjust atom name? */
int usr_aan_flag = -99;         /*the user's input for adjusting atom names */
int duplicatedname_flag = 0;    /*check atom name duplication? */
int cartcoord_flag = 0;         /*generate coordinate from internal coordinate ? */
int connect_flag = 0;           /*judge atom connectivity and generate bond ? */
int chargemethod_flag = 0;      /* is there a "-c" flag on the input line? */
int divcon_flag = 2;
int ek_flag = 0;                /* read empircal calculation keyword or not */
int max_path_length = -1;

int atomnum = 0;                // 1 more than the highest valid atom array index.
int bondnum = 0;
int ringnum = 0;
ATOM *atom;                     // the atom array.
BOND *bond;
RING *ring;
AROM *arom;
MOLINFO minfo;
MOLINFO minfo2;
CONTROLINFO cinfo;

int atomnum_tmp = 0;
int bondnum_tmp = 0;
ATOM *atom_tmp;
BOND *bond_tmp;
char ifilename[MAXCHAR];
char ofilename[MAXCHAR];
char afilename[MAXCHAR] = "";
char cfilename[MAXCHAR];

# include "fileformat.c"
# include "checkmolecule.c"

/*The following four functions, read_at(), write_at(), read_bt() and write_bt() are
used in amber sybyl interface development and are unreachable to the users
*/
int ra_flag = 0;
int rb_flag = 0;
int wa_flag = 0;
int wb_flag = 0;
char at_filename[MAXCHAR];
char bt_filename[MAXCHAR];

size_t copied_size;

void usage()
{
    if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
        printf("[31mUsage: antechamber -i    [0m input file name\n"
               "[31m                   -fi   [0m input file format\n"
               "[31m                   -o    [0m output file name\n"
               "[31m                   -fo   [0m output file format\n"
               "[31m                   -c    [0m charge method\n"
               "[31m                   -cf   [0m charge file name\n"
               "[31m                   -nc   [0m net molecular charge (int)\n"
               "[31m                   -a    [0m additional file name\n"
               "[31m                   -fa   [0m additional file format\n"
               "[31m                   -ao   [0m additional file operation\n"
               "[34m                          crd   [0m: only read in coordinate\n"
               "[34m                          crg   [0m: only read in charge\n"
               "[34m                          radius[0m: only read in radius\n"
               "[34m                          name  [0m: only read in atom name\n"
               "[34m                          type  [0m: only read in atom type\n"
               "[34m                          bond  [0m: only read in bond type \n"
               "[31m                   -m    [0m multiplicity (2S+1), default is 1\n"
               "[31m                   -rn   [0m residue name, overrides input file, default is MOL\n"
               "[31m                   -rf   [0m residue toplogy file name in prep input file,\n"
               "                          default is molecule.res\n"
               "[31m                   -ch   [0m check file name for gaussian, default is 'molecule'\n"
               "[31m                   -ek   [0m mopac or sqm keyword, inside quotes; overwrites previous ones\n"
               "[31m                   -gk   [0m gaussian job keyword, inside quotes, is ignored when both -gopt and -gsp are used\n"
               "[31m                   -gopt [0m gaussian job keyword for optimization, inside quotes\n"
               "[31m                   -gsp  [0m gaussian job keyword for single point calculation, inside quotes\n"
               "[31m                   -gm   [0m gaussian memory keyword, inside quotes, such as \"%%mem=1000MB\"\n"
               "[31m                   -gn   [0m gaussian number of processors keyword, inside quotes, such as \"%%nproc=8\"\n"
               "[31m                   -gdsk [0m gaussian maximum disk usage keyword, inside quotes, such as \"%%maxdisk=50GB\"\n"
               "[31m                   -gv   [0m add keyword to generate gesp file (for Gaussian 09 only)\n"
               "[34m                          1    [0m: yes\n"
               "[34m                          0    [0m: no, the default\n"
               "[31m                   -ge   [0m gaussian esp file generated by iop(6/50=1), default is g09.gesp\n"
               "[31m                   -tor  [0m torsional angle list, inside a pair of quotes, such as \"1-2-3-4:0,5-6-7-8\"\n"
               "[34m                         [0m ':1' or ':0' indicates the torsional angle is frozen or not\n"
               "[31m                   -df   [0m am1-bcc precharge flag, 2 - use sqm(default); 0 - use mopac\n"
               "[31m                   -at   [0m atom type\n"
               "[34m                          gaff [0m: the default\n"
               "[34m                          gaff2[0m: for gaff2 (beta-version)\n"
               "[34m                          amber[0m: for PARM94/99/99SB\n"
               "[34m                          bcc  [0m: bcc \n"
               "[34m                          abcg2[0m: abcg2 \n"
               "[34m                          sybyl[0m: sybyl \n"
               "[31m                   -du   [0m fix duplicate atom names: yes(y)[default] or no(n)\n"
               "[31m                   -bk   [0m component/block Id, for ccif\n"
               "[31m                   -an   [0m adjust atom names: yes(y) or no(n)\n"
               "                          the default is 'y' for 'mol2' and 'ac' and 'n' for the other formats \n"
               "[31m                   -j    [0m atom type and bond type prediction index, default is 4 \n"
               "[34m                          0    [0m: no assignment\n"
               "[34m                          1    [0m: atom type \n"
               "[34m                          2    [0m: full  bond types \n"
               "[34m                          3    [0m: part  bond types \n"
               "[34m                          4    [0m: atom and full bond type \n"
               "[34m                          5    [0m: atom and part bond type \n"
               "[31m                   -s    [0m status information: 0(brief), 1(default) or 2(verbose)\n"
               "[31m                   -eq   [0m equalizing atomic charge, default is 1 for '-c resp', '-c bcc', '-c abcg2' and 0 for the other charge methods \n"
               "[34m                          0    [0m: no use\n"
               "[34m                          1    [0m: by atomic paths \n"
               "[34m                          2    [0m: by atomic paths and structural information, i.e. E/Z configurations \n"
               "[31m                   -pf   [0m remove intermediate files: yes(y) or no(n)[default]\n"
               "[31m                   -pl   [0m maximum path length to determin equivalence of atomic charges for resp, bcc and abcg2,\n"
               "[31m                         [0m the smaller the value, the faster the algorithm, default is -1 (use full length),\n"
               "[31m                         [0m set this parameter to 10 to 30 if your molecule is big (# atoms >= 100)\n"
               "[31m                   -seq  [0m atomic sequence order changable: yes(y)[default] or no(n)\n"
               "[31m                   -dr   [0m acdoctor mode: yes(y)[default] or no(n)\n"
               "^[[31m                   -f    ^[[0m The reweighting factor used to scale RESP charges obtained using QUICK, default is 17.00\n"
               "[31m                   -i -o -fi and -fo must appear; others are optional[0m\n"
               "[32m                   Use 'antechamber -L' to list the supported file formats and charge methods[0m\n");
    } else {
        printf("Usage: antechamber -i     input file name\n"
               "                   -fi    input file format\n"
               "                   -o     output file name\n"
               "                   -fo    output file format\n"
               "                   -c     charge method\n"
               "                   -cf    charge file name\n"
               "                   -nc    net molecular charge (int)\n"
               "                   -a     additional file name\n"
               "                   -fa    additional file format\n"
               "                   -ao    additional file operation\n"
               "                          crd   : only read in coordinate\n"
               "                          crg   : only read in charge\n"
               "                          radius: only read in radius\n"
               "                          name  : only read in atom name\n"
               "                          type  : only read in atom type\n"
               "                          bond  : only read in bond type \n"
               "                   -m     multiplicity (2S+1), default is 1\n"
               "                   -rn    residue name, overrides input file, default is MOL\n"
               "                   -rf    residue toplogy file name in prep input file,\n"
               "                          default is molecule.res\n"
               "                   -ch    check file name for gaussian, default is 'molecule'\n"
               "                   -ek    mopac or sqm keyword, inside quotes; overwrites previous ones\n"
               "                   -gk    gaussian job keyword, inside quotes, is ignored when both -gopt and -gsp are used\n"
               "                   -gopt  gaussian job keyword for optimization, inside quotes\n"
               "                   -gsp   gaussian job keyword for single point calculation, inside quotes\n"
               "                   -gm    gaussian memory keyword, inside quotes, such as \"%%mem=1000MB\"\n"
               "                   -gn    gaussian number of processors keyword, inside quotes, such as \"%%nproc=8\"\n"
               "                   -gdsk  gaussian maximum disk usage keyword, inside quotes, such as \"%%maxdisk=50GB\"\n"
               "                   -gv    add keyword to generate gesp file (for Gaussian 09 only)\n"
               "                          1    : yes\n"
               "                          0    : no, the default\n"
               "                   -ge    gaussian esp file generated by iop(6/50=1), default is g09.gesp\n"
               "                   -tor   torsional angle list, inside a pair of quotes, such as \"1-2-3-4:0,5-6-7-8\"\n"
               "                          ':1' or ':0' indicates the torsional angle is frozen or not\n"
               "                   -df    am1-bcc precharge flag, 2 - use sqm(default); 0 - use mopac\n"
               "                   -at    atom type\n"
               "                          gaff : the default\n"
               "                          gaff2: for gaff2 (beta-version)\n"
               "                          amber: for PARM94/99/99SB\n"
               "                          bcc  : bcc \n"
               "                          abcg2: abcg2 \n"
               "                          sybyl: sybyl \n"
               "                   -du    fix duplicate atom names: yes(y)[default] or no(n)\n"
               "                   -bk    component/block Id, for ccif\n"
               "                   -an    adjust atom names: yes(y) or no(n)\n"
               "                          the default is 'y' for 'mol2' and 'ac' and 'n' for the other formats \n"
               "                   -j     atom type and bond type prediction index, default is 4 \n"
               "                          0    : no assignment\n"
               "                          1    : atom type \n"
               "                          2    : full  bond types \n"
               "                          3    : part  bond types \n"
               "                          4    : atom and full bond type \n"
               "                          5    : atom and part bond type \n"
               "                   -s     status information: 0(brief), 1(default) or 2(verbose)\n"
               "                   -eq    equalizing atomic charge, default is 1 for '-c resp', '-c bcc', '-c abcg2' and 0 for the other charge methods \n"
               "                          0    : no use\n"
               "                          1    : by atomic paths \n"
               "                          2    : by atomic paths and structural information, i.e. E/Z configurations \n"
               "                   -pf    remove intermediate files: yes(y) or no(n)[default]\n"
               "                   -pl    maximum path length to determin equivalence of atomic charges for resp, bcc and abcg2\n"
               "                          the smaller the value, the faster the algorithm, default is -1 (use full length),\n"
               "                          set this parameter to 10 to 30 if your molecule is big (# atoms >= 100)\n"
               "                   -seq   atomic sequence order changable: yes(y)[default] or no(n)\n"
               "                   -dr    acdoctor mode: yes(y)[default] or no(n)\n"
               "                   -f     The reweighting factor used to scale RESP charges obtained using QUICK, default is 17.00\n"
               "                   -i -o -fi and -fo must appear; others are optional[0m\n"
               "                   Use 'antechamber -L' to list the supported file formats and charge methods[0m\n");
    }
}

void list()
{

    if (strcmp(COLORTEXT, "YES") == 0 || strcmp(COLORTEXT, "yes") == 0) {
        printf("	         	    [31m List of the File Formats [0m \n"
             "\n	 	file format type  abbre. index | file format type abbre. index"
             "\n		--------------------------------------------------------------"
             "\n		Antechamber        ac       1  | Sybyl Mol2         mol2    2 "
             "\n		PDB                pdb      3  | Modified PDB       mpdb    4 "
             "\n		AMBER PREP (int)   prepi    5  | AMBER PREP (car)   prepc   6 "
             "\n		Gaussian Z-Matrix  gzmat    7  | Gaussian Cartesian gcrt    8 "
             "\n		Mopac Internal     mopint   9  | Mopac Cartesian    mopcrt 10 "
             "\n		Gaussian Output    gout    11  | Mopac Output       mopout 12 "
             "\n		Alchemy            alc     13  | CSD                csd    14 "
             "\n		MDL                mdl     15  | Hyper              hin    16 "
             "\n		AMBER Restart      rst     17  | Jaguar Cartesian   jcrt   18 "
             "\n		Jaguar Z-Matrix    jzmat   19  | Jaguar Output      jout   20 "
             "\n		Divcon Input       divcrt  21  | Divcon Output      divout 22 "
             "\n		SQM Input          sqmcrt  23  | SQM Output         sqmout 24 "
             "\n		Charmm             charmm  25  | Gaussian ESP       gesp   26 "
             "\n		geostd cif         ccif    27  | GAMESS dat         gamess 28 "
             "\n		Orca input         orcinp  29  | Orca output        orcout 30 "
             "\n		pdbqt              pdbqt   31  |"
             "\n		--------------------------------------------------------------\n"
               "                AMBER restart file can only be read in as additional file.\n"
             "\n	         	    [31m List of the Charge Methods [0m \n"
             "\n		charge method     abbre. index | charge method    abbre. index"
             "\n		--------------------------------------------------------------"
             "\n		RESP               resp     1  |  AM1-BCC           bcc     2"
             "\n		CM1                cm1      3  |  CM2               cm2     4"
             "\n		ESP (Kollman)      esp      5  |  Mulliken          mul     6"
             "\n		Gasteiger          gas      7  |  ABCG2             abcg2   8"
             "\n		Read in charge     rc       9  |  Write out charge  wc     10"
             "\n		Delete Charge      dc      11  |"
             "\n		--------------------------------------------------------------\n");
    } else {
        printf("	         	     List of the File Formats \n"
             "\n	 	file format type  abbre. index | file format type abbre. index"
             "\n		--------------------------------------------------------------"
             "\n		Antechamber        ac       1  | Sybyl Mol2         mol2    2 "
             "\n		PDB                pdb      3  | Modified PDB       mpdb    4 "
             "\n		AMBER PREP (int)   prepi    5  | AMBER PREP (car)   prepc   6 "
             "\n		Gaussian Z-Matrix  gzmat    7  | Gaussian Cartesian gcrt    8 "
             "\n		Mopac Internal     mopint   9  | Mopac Cartesian    mopcrt 10 "
             "\n		Gaussian Output    gout    11  | Mopac Output       mopout 12 "
             "\n		Alchemy            alc     13  | CSD                csd    14 "
             "\n		MDL                mdl     15  | Hyper              hin    16 "
             "\n		AMBER Restart      rst     17  | Jaguar Cartesian   jcrt   18 "
             "\n		Jaguar Z-Matrix    jzmat   19  | Jaguar Output      jout   20 "
             "\n		Divcon Input       divcrt  21  | Divcon Output      divout 22 "
             "\n		SQM Input          sqmcrt  23  | SQM Output         sqmout 24 "
             "\n		Charmm             charmm  25  | Gaussian ESP       gesp   26 "
             "\n		geostd cif         ccif    27  | GAMESS dat         gamess 28 "
             "\n		Orca input         orcinp  29  | Orca output        orcout 30 "
             "\n		pdbqt              pdbqt   31  | QUICK output       quick  32 "
             "\n		--------------------------------------------------------------\n"
               "                AMBER restart file can only be read in as additional file.\n"
             "\n	         	     List of the Charge Methods \n"
             "\n		charge method     abbre. index | charge method    abbre. index"
             "\n		--------------------------------------------------------------  "
             "\n		RESP               resp     1  |  AM1-BCC            bcc    2"
             "\n		CM1                cm1      3  |  CM2                cm2    4"
             "\n		ESP (Kollman)      esp      5  |  Mulliken           mul    6"
             "\n		Gasteiger          gas      7  |  ABCG2             abcg2   8"
             "\n		Read in charge     rc       9  |  Write out charge  wc     10"
             "\n		Delete Charge      dc      11  |"
             "\n		--------------------------------------------------------------\n");
    }
}


void judgebondtype(int atomnum, ATOM * atom, int bondnum, BOND * bond, CONTROLINFO cinfo,
                   MOLINFO minfo, int bondtype_flag, int adjustatomname_flag)
{
    size_t copied_size;
    char *options;
    char tmpchar[MAXCHAR];
    wac("ANTECHAMBER_BOND_TYPE.AC0", atomnum, atom, bondnum, bond, cinfo, minfo);

    copied_size = build_exe_path(tmpchar, "bondtype", sizeof tmpchar, 1);
    if (bondtype_flag == 1)
        options =
            " -j part -i ANTECHAMBER_BOND_TYPE.AC0" " -o ANTECHAMBER_BOND_TYPE.AC -f ac";
    else
        options =
            " -j full -i ANTECHAMBER_BOND_TYPE.AC0" " -o ANTECHAMBER_BOND_TYPE.AC -f ac";

    strncat(tmpchar, options, MAXCHAR - copied_size);
    if (adjustatomname_flag == 0)
        strncat(tmpchar, " -an no", 7);

    if (cinfo.intstatus == 2)
        fprintf(stdout, "Running: %s\n", tmpchar);
    esystem(tmpchar);
    minfo2 = minfo;
    rac("ANTECHAMBER_BOND_TYPE.AC", &atomnum, atom, &bondnum, bond, &cinfo, &minfo2);
}

int read_at(char *filename)
{
    FILE *fpin;
    int num = 0;
    char tmpchar1[MAXCHAR];
    char tmpchar2[MAXCHAR];
    char tmpchar3[MAXCHAR];
    char line[MAXCHAR];

    if ((fpin = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "Cannot open file %s in read_at(), exit\n", filename);
        return 0;
    }
    for (;;) {
        if (fgets(line, MAXCHAR, fpin) == NULL)
            break;
        sscanf(line, "%s%s%s", tmpchar1, tmpchar2, tmpchar3);
        strcpy(atom[num++].ambername, tmpchar3);
    }
    fclose(fpin);
    return 0;
}

int read_bt(char *filename)
{
    FILE *fpin;
    int num = 0;
    char tmpchar1[MAXCHAR];
    char tmpchar2[MAXCHAR];
    char tmpchar3[MAXCHAR];
    char line[MAXCHAR];
    if ((fpin = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "Cannot open file %s in read_bt().\n", filename);
        return 0;
    }
    for (;;) {
        if (fgets(line, MAXCHAR, fpin) == NULL)
            break;
        sscanf(line, "%s%s%s%d", tmpchar1, tmpchar2, tmpchar3, &bond[num++].type);
    }
    fclose(fpin);
    return 0;

}

int write_at(char *filename)
{
    FILE *fpout;
    int i;
    if ((fpout = fopen(filename, "w")) == NULL) {
        fprintf(stderr, "Cannot open file %s to write in write_at().\n", filename);
        return 0;
    }
    for (i = 0; i < atomnum; i++)
        fprintf(fpout, "%5d %5s %5s\n", i + 1, atom[i].name, atom[i].ambername);
    fclose(fpout);
    return 0;

}

int write_bt(char *filename)
{
    FILE *fpout;
    int i;
    if ((fpout = fopen(filename, "w")) == NULL) {
        fprintf(stderr, "Cannot open file %s to write in write_bt().\n", filename);
        return 0;
    }
    for (i = 0; i < bondnum; i++)
        fprintf(fpout, "%5d %5d %5d %5d\n", i + 1, bond[i].bondi + 1, bond[i].bondj,
                bond[i].type);
    fclose(fpout);
    return 0;

}

int main(int argc, char *argv[])
{
    int i;
    int index;
    char tmpchar[MAXCHAR];
    double fraction;
    double tmpf;
    int overflow_flag = 0;      /* if overflow_flag ==1, reallocate memory */
    char *atomtype_args;
    char ftype[3];

    if( cinfo.intstatus) 
    fprintf(stdout, "\nWelcome to antechamber %s: molecular input file processor.\n\n",
            ANTECHAMBER_VERSION);
    esetprogramname(argv[0]);
    amberhome = (char *) getenv("AMBERCLASSICHOME");
    if (amberhome == NULL) {
        eprintf("AMBERCLASSICHOME is not set.");
    }
    if (argc == 2) {
        if (strncmp(argv[1], "-h", 2) == 0 || strncmp(argv[1], "-H", 2) == 0) {
            usage();
            exit(0);
        }
        if (strncmp(argv[1], "-l", 2) == 0 || strncmp(argv[1], "-L", 2) == 0
            || strncmp(argv[1], "-list", 4) == 0 || strncmp(argv[1], "-List", 4) == 0
            || strncmp(argv[1], "-LIST", 4) == 0) {
            list();
            exit(0);
        }
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
    default_flag = 0;
    atomname_flag = 0;
    atomicnum_flag = 0;
    adjustatomname_flag = 0;
    duplicatedname_flag = 1;
    cartcoord_flag = 0;
    connect_flag = 0;
    chargemethod_flag = 0;

    index = 0;
    for (i = 1; i < argc - 1; i += 2) {
        if (strcmp(argv[i], "-i") == 0) {
            index++;
            strcpy(ifilename, argv[i + 1]);
            continue;
        } else if (strcmp(argv[i], "-o") == 0) {
            index++;
            strcpy(ofilename, argv[i + 1]);
            continue;
        } else if (strcmp(argv[i], "-fi") == 0) {
            index++;
            strcpy(cinfo.intype, argv[i + 1]);
            continue;
        } else if (strcmp(argv[i], "-fo") == 0) {
            index++;
            strcpy(cinfo.outtype, argv[i + 1]);
            continue;
        } else if (strcmp(argv[i], "-j") == 0) {
            // TODO make this a function and apply to -df, -s, etc.
            errno = 0;
            char *endp;
            cinfo.prediction_index = (int) strtol(argv[i + 1], & endp, 10);
            if (errno != 0) {
                perror("Error: errno is nonzero");
            }
            if (endp == argv[i + 1] || *endp != '\0' || errno != 0) {
                eprintf("Invalid argument (%s) to -j option.", argv[i + 1]);
            }
            continue;
        } else if (strcmp(argv[i], "-a") == 0) {
            strcpy(afilename, argv[i + 1]);
            continue;
        } else if (strcmp(argv[i], "-fa") == 0) {
            strcpy(cinfo.atype, argv[i + 1]);
            continue;
        } else if (strcmp(argv[i], "-c") == 0) {
            strcpy(cinfo.chargetype, argv[i + 1]);
            chargemethod_flag = 1;
            continue;
        } else if (strcmp(argv[i], "-m") == 0) {
            minfo.multiplicity = atoi(argv[i + 1]);
            continue;
        } else if (strcmp(argv[i], "-s") == 0) {
            cinfo.intstatus = atoi(argv[i + 1]);
            continue;
        } else if (strcmp(argv[i], "-f") == 0) {
            strcpy(cinfo.sfac, argv[i + 1]);
            continue;
        } else if (strcmp(argv[i], "-ra") == 0) {
            strcpy(at_filename, argv[i + 1]);
            ra_flag = 1;
            continue;
        } else if (strcmp(argv[i], "-rb") == 0) {
            strcpy(bt_filename, argv[i + 1]);
            rb_flag = 1;
            continue;
        } else if (strcmp(argv[i], "-wa") == 0) {
            strcpy(at_filename, argv[i + 1]);
            wa_flag = 1;
            continue;
        } else if (strcmp(argv[i], "-wb") == 0) {
            strcpy(bt_filename, argv[i + 1]);
            wb_flag = 1;
            continue;
        } else if (strcmp(argv[i], "-bk") == 0) {
            strncpy(blockId, argv[i + 1], MAX_CIF_BLOCK_CODE_LENGTH);
            blockId[MAX_CIF_BLOCK_CODE_LENGTH] = '\0';
            continue;
        } else if (strcmp(argv[i], "-ao") == 0) {
            if (strcmp(argv[i + 1], "crd") == 0 || strcmp(argv[i + 1], "CRD") == 0) {
                ao_flag = 1;
                continue;
            }
            if (strcmp(argv[i + 1], "crg") == 0 || strcmp(argv[i + 1], "CRG") == 0) {
                ao_flag = 2;
                continue;
            }
            if (strcmp(argv[i + 1], "name") == 0 || strcmp(argv[i + 1], "NAME") == 0) {
                ao_flag = 3;
                continue;
            }
            if (strcmp(argv[i + 1], "type") == 0 || strcmp(argv[i + 1], "NAME") == 0) {
                ao_flag = 4;
                continue;
            }
            if (strcmp(argv[i + 1], "bond") == 0 || strcmp(argv[i + 1], "BOND") == 0) {
                ao_flag = 5;
                continue;
            }
            if (strcmp(argv[i + 1], "radius") == 0 || strcmp(argv[i + 1], "RADIUS") == 0) {
                ao_flag = 6;
                continue;
            }
        } else if (strcmp(argv[i], "-cf") == 0) {
            strcpy(cfilename, argv[i + 1]);
            continue;
        } else if (strcmp(argv[i], "-nc") == 0) {
            minfo.usercharge = atoi(argv[i + 1]);
            continue;
        } else if (strcmp(argv[i], "-rn") == 0) {
            strcpy(minfo.longresname, argv[i + 1]);
/*			strncpy(minfo.resname, argv[i + 1], 3); */
            strcpy(minfo.resname, argv[i + 1]);
            cinfo.rnindex = 1;
            continue;
        } else if (strcmp(argv[i], "-ch") == 0) {
            strcpy(minfo.chkfile, argv[i + 1]);
            continue;
        } else if (strcmp(argv[i], "-at") == 0) {
            strcpy(minfo.atom_type_def, argv[i + 1]);
            continue;
        } else if (strcmp(argv[i], "-du") == 0) {
            if (strcmp(argv[i + 1], "Yes") == 0 || strcmp(argv[i + 1], "Y") == 0
                || strcmp(argv[i + 1], "yes") == 0 || strcmp(argv[i + 1], "y") == 0)
                duplicatedname_flag = 1;
            if (strcmp(argv[i + 1], "NO") == 0 || strcmp(argv[i + 1], "No") == 0
                || strcmp(argv[i + 1], "N") == 0 || strcmp(argv[i + 1], "no") == 0
                || strcmp(argv[i + 1], "n") == 0)
                duplicatedname_flag = 0;
            continue;
        } else if (strcmp(argv[i], "-an") == 0) {
            // TODO make this a function and apply to other y/n options
            if (strcmp(argv[i + 1], "Yes") == 0 || strcmp(argv[i + 1], "Y") == 0
                || strcmp(argv[i + 1], "yes") == 0 || strcmp(argv[i + 1], "y") == 0)
                usr_aan_flag = 1;
            else if (strcmp(argv[i + 1], "NO") == 0 || strcmp(argv[i + 1], "No") == 0
                     || strcmp(argv[i + 1], "N") == 0 || strcmp(argv[i + 1], "no") == 0
                     || strcmp(argv[i + 1], "n") == 0)
                usr_aan_flag = 0;
            else {
                eprintf("Invalid argument (%s) to -an option.", argv[i + 1]);
            }
            continue;
        } else if (strcmp(argv[i], "-rf") == 0) {
            strcpy(minfo.resfilename, argv[i + 1]);
            continue;
        } else if (strcmp(argv[i], "-pf") == 0) {
            if (strcmp(argv[i + 1], "yes") == 0 || strcmp(argv[i + 1], "YES") == 0
                || strcmp(argv[i + 1], "Y") == 0 || strcmp(argv[i + 1], "y") == 0
                || strcmp(argv[i + 1], "Yes") == 0)
                cinfo.pfindex = 1;
            if (strcmp(argv[i + 1], "no") == 0 || strcmp(argv[i + 1], "NO") == 0
                || strcmp(argv[i + 1], "N") == 0 || strcmp(argv[i + 1], "n") == 0
                || strcmp(argv[i + 1], "No") == 0)
                cinfo.pfindex = 0;
            continue;
        } else if (strcmp(argv[i], "-ek") == 0) {
            strcpy(minfo.ekeyword, argv[i + 1]);
            ek_flag = 1;
            continue;
        } else if (strcmp(argv[i], "-gk") == 0) {
            strcpy(minfo.gkeyword, argv[i + 1]);
            minfo.igkeyword = 1;
            continue;
        } else if (strcmp(argv[i], "-gopt") == 0) {
            strcpy(minfo.gopt, argv[i + 1]);
            minfo.igopt = 1;
            continue;
        } else if (strcmp(argv[i], "-gsp") == 0) {
            strcpy(minfo.gsp, argv[i + 1]);
            minfo.igsp = 1;
            continue;
        } else if (strcmp(argv[i], "-gdsk") == 0) {
            strcpy(minfo.gdsk, argv[i + 1]);
            minfo.igdsk = 1;
            continue;
        } else if (strcmp(argv[i], "-tor") == 0) {
            strcpy(minfo.tor, argv[i + 1]);
            minfo.itor = 1;
            continue;
        } else if (strcmp(argv[i], "-gm") == 0) {
            strcpy(minfo.gm, argv[i + 1]);
            continue;
        } else if (strcmp(argv[i], "-gn") == 0) {
            strcpy(minfo.gn, argv[i + 1]);
            continue;
        } else if (strcmp(argv[i], "-ge") == 0) {
            strcpy(minfo.gesp, argv[i + 1]);
            continue;
        } else if (strcmp(argv[i], "-gv") == 0) {
            minfo.gv = atoi(argv[i + 1]);
            continue;
        } else if (strcmp(argv[i], "-eq") == 0) {
            minfo.eqcharge = atoi(argv[i + 1]);
            continue;
        } else if (strcmp(argv[i], "-df") == 0) {
            divcon_flag = atoi(argv[i + 1]);
            if (divcon_flag != 0 && divcon_flag != 2) {
                eprintf("Invalid argument (%s) to -df option.", argv[i + 1]);
            }
            minfo.divcon = divcon_flag;
            continue;
        } else if (strcmp(argv[i], "-pl") == 0) {
            max_path_length = atoi(argv[i + 1]);
            continue;
        } else if (strcmp(argv[i], "-seq") == 0) {
            cinfo.atseq = 1;
            if(strcmp(argv[i + 1], "no") == 0 || strcmp(argv[i + 1], "NO") == 0 ||
               strcmp(argv[i + 1], "No") == 0 || strcmp(argv[i + 1], "n") == 0  || strcmp(argv[i + 1], "N") == 0)
                    cinfo.atseq = 0;
       	    continue;
        } else if (strcmp(argv[i], "-dr") == 0) {
            if (strcmp(argv[i + 1], "yes") == 0 || strcmp(argv[i + 1], "YES") == 0
                || strcmp(argv[i + 1], "Y") == 0 || strcmp(argv[i + 1], "y") == 0
                || strcmp(argv[i + 1], "Yes") == 0) {
                checkformat = 1;
                checkmolecule = 1;
            }
            if (strcmp(argv[i + 1], "no") == 0 || strcmp(argv[i + 1], "NO") == 0
                || strcmp(argv[i + 1], "N") == 0 || strcmp(argv[i + 1], "n") == 0
                || strcmp(argv[i + 1], "No") == 0) {
                checkformat = 0;
                checkmolecule = 0;
            }
            continue;
        } else {
            eprintf("Unrecognized option (%s).\n"
                    "Use antechamber -h for command-line syntax.", argv[i]);
        }
    }

    if (index != 4 && (strcmp(cinfo.chargetype, "wc") != 0)) {
        eprintf("Need both input and output files and formats.\n"
                "Use antechamber -h for command-line syntax.");
    }

    if (argc != i) {
        eprintf("Number of arguments is odd - arguments must come in pairs.\n"
                "Use antechamber -h for command-line syntax.");
    }

/*	reset charge equlibration option*/
    if (minfo.eqcharge == -1) {
        minfo.eqcharge = 0;
        if (strcmp(cinfo.chargetype, "resp") == 0 || strcmp(cinfo.chargetype, "RESP") == 0
            || strcmp(cinfo.chargetype, "1") == 0)
            minfo.eqcharge = 1;
        if (strcmp(cinfo.chargetype, "bcc") == 0 || strcmp(cinfo.chargetype, "BCC") == 0
            || strcmp(cinfo.chargetype, "2") == 0)
            minfo.eqcharge = 1;
        if (strcmp(cinfo.chargetype, "abcg2") == 0 || strcmp(cinfo.chargetype, "ABCG2") == 0
            || strcmp(cinfo.chargetype, "8") == 0)
            minfo.eqcharge = 1;
    }
    if (minfo.eqcharge == 1)
        cinfo.max_path_length = max_path_length;
/* 	set ekeyword if it is not read in*/
    if (ek_flag == 0) {
        if (divcon_flag == 0)
            strcpy(minfo.ekeyword, minfo.mkeyword);
        else if (divcon_flag == 1)
            strcpy(minfo.ekeyword, minfo.dkeyword);
        else if (divcon_flag == 2)
            strcpy(minfo.ekeyword, minfo.skeyword);
    }

    create_output_file_comment(argc, argv);
/* 	for connect.tpl and radius parameter files */
    build_dat_path(minfo.connect_file, "CONNECT.TPL", sizeof minfo.connect_file, 0);
    build_dat_path(minfo.radius_file, "RADIUS.DAT", sizeof minfo.radius_file, 0);
/*      allocate memory using default parameters MAXATOM and MAXBOND */
    memory(0, MAXATOM, MAXBOND, MAXRING);

    if ( strcmp("quick", cinfo.outtype) == 0 ) {

      // quick  input file has to be a .in file
      strcpy(ftype,ofilename+strlen(ofilename)-3);
      if(strcmp(ftype,".in") != 0) eprintf("output file needs to be in .in format\n");

      if ( strcmp("pdb", cinfo.intype) == 0 ) {

        atom_tmp = (ATOM *) emalloc(sizeof(ATOM) * cinfo.maxatom);
        pdb_to_quick_input(ifilename, ofilename, atom_tmp, &minfo);

        return(0);

      }else if ( strcmp("xyz", cinfo.intype) == 0 ) {

        atom_tmp = (ATOM *) emalloc(sizeof(ATOM) * cinfo.maxatom);
        xyz_to_quick_input(ifilename, ofilename, atom_tmp, &minfo);

        return(0);

      };

    };

    if (checkformat) {
        printf("Info: acdoctor mode is on: check and diagnose problems in the input file.\n");
    }
    if( cinfo.intstatus )
    printf("Info: The atom type is set to %s; the options available to the -at flag are\n"
           "      gaff, gaff2, amber, bcc, abcg2, and sybyl.\n\n", minfo.atom_type_def);
    read_and_validate_input_file();

/* when read in a mopout, divout, or sqmout file, always output a pdb file
   that has the optimized coordinates */

    if (strcmp("mopout", cinfo.intype) == 0 || strcmp("12", cinfo.intype) == 0) {
        atom_tmp = (ATOM *) emalloc(sizeof(ATOM) * cinfo.maxatom);
        for (i = 0; i < atomnum; i++) {
            atom_tmp[i] = atom[i];
            atom_tmp[i].connum = 0;
        }
        rmopout_coord(ifilename, atom_tmp);
        wpdb("mopac.pdb", atomnum, atom_tmp);
        free(atom_tmp);
    } else if (strcmp("divout", cinfo.intype) == 0 || strcmp("22", cinfo.intype) == 0) {
        atom_tmp = (ATOM *) emalloc(sizeof(ATOM) * cinfo.maxatom);
        for (i = 0; i < atomnum; i++) {
            atom_tmp[i] = atom[i];
            atom_tmp[i].connum = 0;
        }
        rdivout_coord(ifilename, atom_tmp);
        wpdb("divcon.pdb", atomnum, atom_tmp);
        free(atom_tmp);
    } else if (strcmp("sqmout", cinfo.intype) == 0 || strcmp("24", cinfo.intype) == 0) {
        atom_tmp = (ATOM *) emalloc(sizeof(ATOM) * cinfo.maxatom);
        for (i = 0; i < atomnum; i++) {
            atom_tmp[i] = atom[i];
            atom_tmp[i].connum = 0;
        }
        rsqmout_coord(ifilename, atom_tmp);
        wpdb("sqm.pdb", atomnum, atom_tmp);
        free(atom_tmp);
    } else {
        // empty else; this is not an error condition.
    }

/*****************************************************************************/
/* Begin big if over additional input file formats */
/* This very big if reads in the additional files */
/*****************************************************************************/

    if (strlen(afilename) > 1 && strlen(cinfo.atype) > 1) {
        // expand the array size a little bit, since atom file type such as 
        // prepi may have additional atom records.
        cinfo.maxatom = atomnum + 10;
        cinfo.maxbond = bondnum + 10;
        atom_tmp = (ATOM *) emalloc(sizeof(ATOM) * cinfo.maxatom);
        bond_tmp = (BOND *) emalloc(sizeof(BOND) * cinfo.maxbond);
        for (i = 0; i < cinfo.maxbond; ++i) {
            bond[i].jflag = -1;  // bond type has not been assigned
        }

        if (strcmp("ac", cinfo.atype) == 0 || strcmp("1", cinfo.atype) == 0) {
            minfo2 = minfo;
            overflow_flag =
                rac(afilename, &atomnum_tmp, atom_tmp, &bondnum_tmp, bond_tmp, &cinfo,
                    &minfo2);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }
        if (strcmp("mol2", cinfo.atype) == 0 || strcmp("2", cinfo.atype) == 0) {
            overflow_flag =
                rmol2(afilename, &atomnum_tmp, atom_tmp, &bondnum_tmp, bond_tmp, &cinfo,
                      &minfo, 0);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }
        if (strcmp("mopint", cinfo.atype) == 0 || strcmp("9", cinfo.atype) == 0) {
            overflow_flag = rmopint(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }
        if (strcmp("mopcrt", cinfo.atype) == 0 || strcmp("10", cinfo.atype) == 0) {
            overflow_flag = rmopcrt(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }

        if (strcmp("mopout", cinfo.atype) == 0 || strcmp("12", cinfo.atype) == 0) {
            overflow_flag = rmopout(afilename, &atomnum_tmp, atom_tmp, &cinfo, &minfo);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }

        if (strcmp("orcinp", cinfo.atype) == 0 || strcmp("29", cinfo.atype) == 0) {
            overflow_flag = rorca(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }

        if (strcmp("gcrt", cinfo.atype) == 0 || strcmp("8", cinfo.atype) == 0) {
            overflow_flag = rgcrt(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }
        if (strcmp("gzmat", cinfo.atype) == 0 || strcmp("7", cinfo.atype) == 0) {
            overflow_flag = rgzmat(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }

        if (strcmp("orcout", cinfo.atype) == 0 || strcmp("30", cinfo.atype) == 0) {
            overflow_flag = rorcout(afilename, &atomnum_tmp, atom_tmp, cinfo, &minfo);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }

        if (strcmp("gout", cinfo.atype) == 0 || strcmp("11", cinfo.atype) == 0) {
            overflow_flag = rgout(afilename, &atomnum_tmp, atom_tmp, cinfo, &minfo);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }

        if (strcmp("quick", cinfo.atype) == 0 || strcmp("32", cinfo.atype) == 0) {
            overflow_flag = rqout(afilename, &atomnum_tmp, atom_tmp, cinfo, &minfo);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }

        if (strcmp("gamess", cinfo.atype) == 0 || strcmp("27", cinfo.atype) == 0) {
            overflow_flag = rgamess(afilename, &atomnum_tmp, atom_tmp, cinfo, &minfo);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }

        if (strcmp("jcrt", cinfo.atype) == 0 || strcmp("18", cinfo.atype) == 0) {
            overflow_flag = rjcrt(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }

        if (strcmp("jzmat", cinfo.atype) == 0 || strcmp("19", cinfo.atype) == 0) {
            overflow_flag = rjzmat(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }

        if (strcmp("jout", cinfo.atype) == 0 || strcmp("20", cinfo.atype) == 0) {
            overflow_flag = rjout(afilename, &atomnum_tmp, atom_tmp, cinfo, &minfo);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }

        if (strcmp("pdb", cinfo.atype) == 0 || strcmp("3", cinfo.atype) == 0) {
            overflow_flag = rpdb(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo, 0);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }

        if (strcmp("pdbqt", cinfo.atype) == 0 || strcmp("31", cinfo.atype) == 0) {
            overflow_flag = rpdbqt(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }

        if (strcmp("mpdb", cinfo.atype) == 0 || strcmp("4", cinfo.atype) == 0) {
            overflow_flag = rpdb(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo, 1);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }

        if (strcmp("csd", cinfo.atype) == 0 || strcmp("14", cinfo.atype) == 0) {
            overflow_flag = rcsd(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }

        if (strcmp("mdl", cinfo.atype) == 0 || strcmp("sd", cinfo.atype) == 0
            || strcmp("sdf", cinfo.atype) == 0 || strcmp("15", cinfo.atype) == 0) {
            overflow_flag =
                rmdl(afilename, &atomnum_tmp, atom_tmp, &bondnum_tmp, bond_tmp, cinfo,
                     minfo);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }

        if (strcmp("alc", cinfo.atype) == 0 || strcmp("13", cinfo.atype) == 0) {
            overflow_flag =
                ralc(afilename, &atomnum_tmp, atom_tmp, &bondnum_tmp, bond_tmp, cinfo,
                     minfo);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }

        if (strcmp("hin", cinfo.atype) == 0 || strcmp("16", cinfo.atype) == 0) {
            overflow_flag =
                rhin(afilename, &atomnum_tmp, atom_tmp, &bondnum_tmp, bond_tmp, cinfo,
                     minfo);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }

        if (strcmp("prepi", cinfo.atype) == 0 || strcmp("5", cinfo.atype) == 0) {
            overflow_flag =
                rprepi(afilename, &atomnum_tmp, atom_tmp, &bondnum_tmp, bond_tmp, &cinfo,
                       &minfo);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }

        if (strcmp("prepc", cinfo.atype) == 0 || strcmp("6", cinfo.atype) == 0) {
            overflow_flag = rprepc(afilename, &atomnum_tmp, atom_tmp, &cinfo, &minfo);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }
        if (strcmp("rst", cinfo.atype) == 0 || strcmp("17", cinfo.atype) == 0) {
            overflow_flag = rrst(afilename, &atomnum_tmp, atom_tmp, cinfo);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }
        if (strcmp("divcrt", cinfo.atype) == 0 || strcmp("21", cinfo.atype) == 0) {
            overflow_flag = rdivcrt(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }

        if (strcmp("divout", cinfo.atype) == 0 || strcmp("22", cinfo.atype) == 0) {
            overflow_flag = rdivout(afilename, &atomnum_tmp, atom_tmp, &cinfo, &minfo);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }
        if (strcmp("sqmcrt", cinfo.atype) == 0 || strcmp("23", cinfo.atype) == 0) {
            overflow_flag = rsqmcrt(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }

        if (strcmp("sqmout", cinfo.atype) == 0 || strcmp("24", cinfo.atype) == 0) {
            overflow_flag = rsqmout(afilename, &atomnum_tmp, atom_tmp, &cinfo, &minfo);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }

        if (strcmp("charmm", cinfo.intype) == 0 || strcmp("25", cinfo.intype) == 0) {
            overflow_flag =
                rcharmm(afilename, &atomnum_tmp, atom_tmp, &bondnum, bond, &cinfo,
                        &minfo);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }

        if (strcmp("gesp", cinfo.intype) == 0 || strcmp("26", cinfo.intype) == 0) {
            overflow_flag = rgesp(afilename, &atomnum_tmp, atom_tmp, cinfo, minfo);
            if (overflow_flag) {
                fprintf(stderr, "Overflow happens for additional files, exit");
                exit(1);
            }
        }
        if (ao_flag == 1)
            for (i = 0; i < atomnum; i++) {
                atom[i].x = atom_tmp[i].x;
                atom[i].y = atom_tmp[i].y;
                atom[i].z = atom_tmp[i].z;
            }

        if (ao_flag == 2)
            for (i = 0; i < atomnum; i++)
                atom[i].charge = atom_tmp[i].charge;

        if (ao_flag == 3)
            for (i = 0; i < atomnum; i++)
                strcpy(atom[i].name, atom_tmp[i].name);
        if (ao_flag == 4)
            for (i = 0; i < atomnum; i++)
                strcpy(atom[i].ambername, atom_tmp[i].ambername);
        if (ao_flag == 5)
            for (i = 0; i < bondnum; i++)
                bond[i].type = bond_tmp[i].type;
        if (ao_flag == 6)
            for (i = 0; i < atomnum; i++)
                atom[i].radius = atom_tmp[i].radius;
        free(atom_tmp);
        free(bond_tmp);
    }

/*****************************************************************************/
/*   End of block over additional input file formats    */
/*****************************************************************************/


/* if charge type is bcc, the prediction_index must be '4'
   (since we need to assign am1-bcc bond types and atom types) */
    if ((strcmp("bcc", cinfo.chargetype) == 0 || strcmp("2", cinfo.chargetype) == 0 ||
         strcmp("abcg2", cinfo.chargetype) == 0 || strcmp("8", cinfo.chargetype) == 0)
        && cinfo.prediction_index < 4) {
        if (cinfo.prediction_index <= 2) {
            fprintf(stdout,
                    "Warning: the antechamber program automatically changes the prediction type from '%d' to '4' for bcc/abcg2 calculations !\n",
                    cinfo.prediction_index);
            cinfo.prediction_index = 4;
        } else if (cinfo.prediction_index == 3) {
            fprintf(stdout,
                    "Warning: the antechamber program automatically changes the prediction type from '3' to '4' for bcc/abcg2 calculations !\n");
            cinfo.prediction_index = 5;
            // TODO should it be 4  ^^^
        }
    }

/*****************************************************************************/
/* Assign atomtype_flag, bondtype_flag and charge_flag according to -j flag */
/*****************************************************************************/

    if (cinfo.prediction_index == 0) {
        atomtype_flag = 0;
        bondtype_flag = 0;
    }
    else if (cinfo.prediction_index == 1) {
        atomtype_flag = 1;
        bondtype_flag = 0;
    }
    else if (cinfo.prediction_index == 2) {
        atomtype_flag = 0;
        bondtype_flag = 2;
    }
    else if (cinfo.prediction_index == 3) {
        atomtype_flag = 0;
        bondtype_flag = 1;
    }
    else if (cinfo.prediction_index == 4) {
        atomtype_flag = 1;
        bondtype_flag = 2;
    }
    else if (cinfo.prediction_index == 5) {
        atomtype_flag = 1;
        bondtype_flag = 1;
    }
    else {
        eprintf("Out of range argument (%i) to -j option.",
            cinfo.prediction_index);
    }

/* the following code judge or assign atom name, atom type, bond type etc according to flags */
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

    if (checkmolecule) {
        check_input_molecule();
        printf("acdoctor mode has completed checking the input file.\n\n");
    }

/* Reassign the connect_flag according to the output types */
    if (strcmp("mopcrt", cinfo.outtype) == 0 || strcmp("mopout", cinfo.outtype) == 0
        || strcmp("gcrt", cinfo.outtype) == 0 || strcmp("gout", cinfo.outtype) == 0
        || strcmp("jcrt", cinfo.outtype) == 0 || strcmp("jout", cinfo.outtype) == 0
        || strcmp("pdb", cinfo.outtype) == 0 || strcmp("rst", cinfo.outtype) == 0
        || strcmp("3", cinfo.outtype) == 0 || strcmp("8", cinfo.outtype) == 0
        || strcmp("10", cinfo.outtype) == 0 || strcmp("11", cinfo.outtype) == 0
        || strcmp("12", cinfo.outtype) == 0 || strcmp("17", cinfo.outtype) == 0
        || strcmp("18", cinfo.outtype) == 0 || strcmp("20", cinfo.outtype) == 0
        || strcmp("29", cinfo.outtype) == 0 || strcmp("30", cinfo.outtype) == 0) {
        connect_flag = 0;
        bondtype_flag = 0;
        atomtype_flag = 0;
    }
    if (strcmp("pdb", cinfo.outtype) == 0 || strcmp("3", cinfo.outtype) == 0) {
        /* Effectively undo connect settings of input molecule */
        for (i = 0; i < atomnum; i++)
            atom[i].connum = 0;
    }

    if (strcmp("mopout", cinfo.outtype) == 0 || strcmp("gout", cinfo.outtype) == 0
        || strcmp("jout", cinfo.outtype) == 0 || strcmp("rst", cinfo.outtype) == 0
        || strcmp("11", cinfo.outtype) == 0 || strcmp("12", cinfo.outtype) == 0
        || strcmp("17", cinfo.outtype) == 0 || strcmp("20", cinfo.outtype) == 0) {
        duplicatedname_flag = 0;
    }

    if (strcmp("mopint", cinfo.outtype) == 0 || strcmp("gzmat", cinfo.outtype) == 0
        || strcmp("jzmat", cinfo.outtype) == 0 || strcmp("7", cinfo.outtype) == 0
        || strcmp("9", cinfo.outtype) == 0 || strcmp("19", cinfo.outtype) == 0) {
        bondtype_flag = 0;
        atomtype_flag = 0;
    }

    if (bondtype_flag && bondnum > 0) {
        judgebondtype(atomnum, atom, bondnum, bond, cinfo, minfo, bondtype_flag, adjustatomname_flag);
        if (cinfo.prediction_index == 2 || cinfo.prediction_index == 3) {
            cinfo.prediction_index = 0;
            atomtype_flag = 0;
            bondtype_flag = 0;
        }
        if (cinfo.prediction_index == 4 || cinfo.prediction_index == 5) {
            cinfo.prediction_index = 1;
            atomtype_flag = 1;
            bondtype_flag = 0;
        }
    }

    if (atomtype_flag) {
        wac("ANTECHAMBER_AC.AC0", atomnum, atom, bondnum, bond, cinfo, minfo);
        atomtype_args = " -i ANTECHAMBER_AC.AC0 -o ANTECHAMBER_AC.AC -p ";
        copied_size = build_exe_path(tmpchar, "atomtype", sizeof tmpchar, 1);
        strncat(tmpchar, atomtype_args, MAXCHAR - copied_size);
        strncat(tmpchar, minfo.atom_type_def,
                MAXCHAR - copied_size - strlen(atomtype_args));
        if (adjustatomname_flag == 0)
            strncat(tmpchar, " -an no", 7);
        if (cinfo.intstatus == 2)
            fprintf(stdout, "\nRunning: %s\n", tmpchar);
        esystem(tmpchar);
        minfo2 = minfo;
        rac("ANTECHAMBER_AC.AC", &atomnum, atom, &bondnum, bond, &cinfo, &minfo2);
    }

/*     the following code readin or calculate charges */
/* 	usercharge info*/
    if (minfo.usercharge > -9999) {     /*charge read in with -nc flag, it is unlikely users input a charge smaller than -9999 */
        minfo.icharge = minfo.usercharge;
        minfo.dcharge = minfo.usercharge;
    } else {
        if (minfo.dcharge < -9990) {
            minfo.dcharge = 0.0;
            for (i = 0; i < atomnum; i++)
                minfo.dcharge += atom[i].charge;
            fraction = modf(minfo.dcharge, &tmpf);
            minfo.icharge = (int) tmpf;
            if (fabs(fraction) >= 0.50) {
                if (minfo.dcharge < 0)
                    minfo.icharge--;
                if (minfo.dcharge > 0)
                    minfo.icharge++;
            }
        }
    }
    if (minfo.multiplicity < -9990)
        minfo.multiplicity = 1;
/*zero weird charges */
    if (minfo.usercharge < -9990 && (minfo.icharge <= -10 || minfo.icharge >= 10)) {
        fprintf(stdout,
                "Warning: Weird total charge (%d)!\n"
                "         Assuming the net charge is 0.\n"
                "         If the weird charge was correct, "
                "specify it via the -nc net charge flag.\n", minfo.icharge);
        minfo.icharge = 0;
    }
    if (strcmp("resp", cinfo.chargetype) == 0 || strcmp("1", cinfo.chargetype) == 0)
        resp(ifilename, atomnum, atom, bondnum, bond, cinfo, minfo);
    else if (strcmp("bcc", cinfo.chargetype) == 0 || strcmp("2", cinfo.chargetype) == 0)
        bcc(ifilename, atomnum, atom, bondnum, bond, arom, &cinfo, &minfo);
    else if (strcmp("cm1", cinfo.chargetype) == 0 || strcmp("3", cinfo.chargetype) == 0)
        cm1(atomnum, atom, &cinfo, &minfo);
    else if (strcmp("cm2", cinfo.chargetype) == 0 || strcmp("4", cinfo.chargetype) == 0)
        cm2(atomnum, atom, &cinfo, &minfo);
    else if (strcmp("esp", cinfo.chargetype) == 0 || strcmp("5", cinfo.chargetype) == 0)
        esp(ifilename, atomnum, atom, cinfo, minfo);
    else if (strcmp("mul", cinfo.chargetype) == 0 || strcmp("6", cinfo.chargetype) == 0)
        mul(ifilename, atomnum, atom, &cinfo, &minfo);
    else if (strcmp("gas", cinfo.chargetype) == 0 || strcmp("7", cinfo.chargetype) == 0)
        gascharge(atomnum, atom, bondnum, bond, cinfo, &minfo);
    else if (strcmp("abcg2", cinfo.chargetype) == 0 || strcmp("8", cinfo.chargetype) == 0)
        bcc(ifilename, atomnum, atom, bondnum, bond, arom, &cinfo, &minfo);
    else if (strcmp("rc", cinfo.chargetype) == 0 || strcmp("9", cinfo.chargetype) == 0)
        rcharge(cfilename, atomnum, atom, cinfo, &minfo);
    else if (strcmp("wc", cinfo.chargetype) == 0 || strcmp("10", cinfo.chargetype) == 0)
        wcharge(cfilename, atomnum, atom, cinfo, minfo);
    else if (strcmp("dc", cinfo.chargetype) == 0 || strcmp("11", cinfo.chargetype) == 0) {
        for (i = 0; i < atomnum; i++)
            atom[i].charge = 0.0;
    } else if (chargemethod_flag == 1) {
        eprintf("Unknown charge method (%s).", cinfo.chargetype);
    }

/*	read the radii*/
    if (strcmp("mpdb", cinfo.intype) != 0 && strcmp("4", cinfo.intype) != 0
        && (strcmp("mpdb", cinfo.outtype) == 0 || strcmp("4", cinfo.outtype) == 0))
        read_radius(minfo.radius_file, atomnum, atom);

    if (ra_flag == 1)
        read_at(at_filename);
    if (rb_flag == 1)
        read_bt(bt_filename);

/*      write out files */
    if (strcmp("ac", cinfo.outtype) == 0 || strcmp("1", cinfo.outtype) == 0) {
        wac(ofilename, atomnum, atom, bondnum, bond, cinfo, minfo);
    } else if (strcmp("charmm", cinfo.outtype) == 0 || strcmp("25", cinfo.outtype) == 0) {
        wcharmm(ofilename, ifilename, atomnum, atom, bondnum, bond, cinfo, minfo);
    } else if (strcmp("mol2", cinfo.outtype) == 0 || strcmp("2", cinfo.outtype) == 0) {
        wmol2(ofilename, atomnum, atom, bondnum, bond, arom, cinfo, minfo);
    } else if (strcmp("mopint", cinfo.outtype) == 0 || strcmp("9", cinfo.outtype) == 0) {
        adjust_sequence_order(atomnum, atom, bondnum, bond);
        wmopint(ofilename, atomnum, atom, minfo);
    } else if (strcmp("mopcrt", cinfo.outtype) == 0 || strcmp("10", cinfo.outtype) == 0) {
        wmopcrt(ofilename, atomnum, atom, minfo);
    } else if (strcmp("mopout", cinfo.outtype) == 0 || strcmp("12", cinfo.outtype) == 0) {
        wmopout(ofilename, atomnum, atom, cinfo, minfo);
    } else if (strcmp("orcinp", cinfo.outtype) == 0 || strcmp("29", cinfo.outtype) == 0) {
        worca(ofilename, atomnum, atom, minfo);
    } else if (strcmp("gcrt", cinfo.outtype) == 0 || strcmp("8", cinfo.outtype) == 0) {
        wgcrt(ofilename, atomnum, atom, minfo);
    } else if (strcmp("gzmat", cinfo.outtype) == 0 || strcmp("7", cinfo.outtype) == 0) {
        adjust_sequence_order(atomnum, atom, bondnum, bond);
        wgzmat(ofilename, atomnum, atom, minfo);
    } else if (strcmp("gout", cinfo.outtype) == 0 || strcmp("11", cinfo.outtype) == 0) {
        wgout();
    } else if (strcmp("gamess", cinfo.outtype) == 0 || strcmp("28", cinfo.outtype) == 0) {
        wgamess();
    } else if (strcmp("jcrt", cinfo.outtype) == 0 || strcmp("18", cinfo.outtype) == 0) {
        wjcrt(ofilename, atomnum, atom, minfo);
    } else if (strcmp("jzmat", cinfo.outtype) == 0 || strcmp("19", cinfo.outtype) == 0) {
        adjust_sequence_order(atomnum, atom, bondnum, bond);
        wjzmat(ofilename, atomnum, atom, minfo);
    } else if (strcmp("jout", cinfo.outtype) == 0 || strcmp("20", cinfo.outtype) == 0) {
        wjout();
    } else if (strcmp("pdb", cinfo.outtype) == 0 || strcmp("3", cinfo.outtype) == 0) {
        wpdb(ofilename, atomnum, atom);
    } else if (strcmp("mpdb", cinfo.outtype) == 0 || strcmp("4", cinfo.outtype) == 0) {
        wmpdb(ofilename, atomnum, atom);
    } else if (strcmp("csd", cinfo.outtype) == 0 || strcmp("14", cinfo.outtype) == 0) {
        wcsd(ofilename, atomnum, atom);
    } else if (strcmp("mdl", cinfo.outtype) == 0 || strcmp("sdf", cinfo.outtype) == 0
               || strcmp("sd", cinfo.outtype) == 0 || strcmp("15", cinfo.outtype) == 0) {
        wmdl(ofilename, atomnum, atom, bondnum, bond, cinfo);
    } else if (strcmp("alc", cinfo.outtype) == 0 || strcmp("13", cinfo.outtype) == 0) {
        walc(ofilename, atomnum, atom, bondnum, bond);
    } else if (strcmp("hin", cinfo.outtype) == 0 || strcmp("16", cinfo.outtype) == 0) {
        whin(ofilename, atomnum, atom, bondnum, bond);
    } else if (strcmp("prepi", cinfo.outtype) == 0 || strcmp("5", cinfo.outtype) == 0) {
        wprep(ofilename, ifilename, atomnum, atom, bondnum, bond, cinfo, &minfo, 1);
    } else if (strcmp("prepc", cinfo.outtype) == 0 || strcmp("6", cinfo.outtype) == 0) {
        wprep(ofilename, ifilename, atomnum, atom, bondnum, bond, cinfo, &minfo, 0);
    } else if (strcmp("rst", cinfo.outtype) == 0 || strcmp("17", cinfo.outtype) == 0) {
        wrst(ofilename, atomnum, atom);
    } else if (strcmp("divcrt", cinfo.outtype) == 0 || strcmp("21", cinfo.outtype) == 0) {
        wdivcrt(ofilename, atomnum, atom, minfo);
    } else if (strcmp("divout", cinfo.outtype) == 0 || strcmp("22", cinfo.outtype) == 0) {
        wdivout(ofilename, atomnum, atom, cinfo, minfo);
    } else if (strcmp("gesp", cinfo.outtype) == 0 || strcmp("26", cinfo.outtype) == 0) {
        wgesp();
    } else if (strcmp("sqmcrt", cinfo.outtype) == 0 || strcmp("23", cinfo.outtype) == 0) {
        wsqmcrt(ofilename, atomnum, atom, minfo);
    } else if (strcmp("sqmout", cinfo.outtype) == 0 || strcmp("24", cinfo.outtype) == 0) {
        wsqmout(ofilename, atomnum, atom, cinfo, minfo);
/*	info(atomnum, atom, bondnum, bond, arom, cinfo, minfo); */
    } else if (strcmp("pdbqt", cinfo.outtype) == 0 || strcmp("31", cinfo.outtype) == 0) {
        wpdbqt();
    } else {
        eprintf("Unknown output file format (%s).", cinfo.outtype);
    }
/*
        free(atom);
        free(bond);
        free(arom);
*/
    if (wa_flag == 1)
        write_at(at_filename);
    if (wb_flag == 1)
        write_bt(bt_filename);
    if (cinfo.pfindex == 1) {
        esystem("rm -f ANTECHAMBER* ATOMTYPE.INF BCCTYPE.INF NEWPDB.PDB PREP.INF");
    }
    printf("\n");
    return (0);
}
