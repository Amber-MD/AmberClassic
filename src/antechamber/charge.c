# include "eprintf.h"
# define CONVERG 0.00001
# define GASMAXITER 500
# define DAMPFACTOR 0.5
GASTEIGER gas[MAXGAS];
MOLINFO minfo2;

void wpdb_optimized(char *filename, int atomnum, ATOM atom[], int flag)
{
    ATOM *atom_tmp;
    int i;
/* when read in a mopout file, always output a pdb file that has the optimized coordinates */
    atom_tmp = (ATOM *) emalloc(sizeof(ATOM) * atomnum);
    for (i = 0; i < atomnum; i++) {
        atom_tmp[i] = atom[i];
        atom_tmp[i].connum = 0;
    }
    if (flag == 0) {
        rmopout_coord(filename, atom_tmp);
        wpdb("mopac.pdb", atomnum, atom_tmp);
    }
    if (flag == 1) {
        rdivout_coord(filename, atom_tmp);
        wpdb("divcon.pdb", atomnum, atom_tmp);
    }
    if (flag == 2) {
        rsqmout_coord(filename, atom_tmp);
        wpdb("sqm.pdb", atomnum, atom_tmp);
    }
    free(atom_tmp);
}

void rsqmcharge(char *filename, int atomnum, ATOM atom[], MOLINFO * minfo)
{

    int index;
    int number = 0;
    int tmpint;
    FILE *fpout;
    char tmpchar[MAXCHAR];
    char tmpchar1[MAXCHAR];
    char tmpchar2[MAXCHAR];
    char tmpchar3[MAXCHAR];
    char tmpchar4[MAXCHAR];
    char tmpchar5[MAXCHAR];
    char tmpchar6[MAXCHAR];
    char line[MAXCHAR];
    double tmpfloat1;

    fpout = efopen(filename, "r");
    index = 0;
    for (;;) {
        if (fgets(line, MAXCHAR, fpout) == NULL)
            break;

        sscanf(line, "%s%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3, tmpchar4, tmpchar5,
               tmpchar6);
        if (strcmp("Mulliken", tmpchar3) == 0) {
            index = 1;
            continue;
        }
        if (index == 1) {
            sscanf(line, "%d%s%lf", &tmpint, tmpchar, &tmpfloat1);
            atom[number].charge = tmpfloat1;
            number++;
        }
    }
    fclose(fpout);
    if (number == 0) {
        eprintf("Unable to find sqm charges in file (%s).\n"
                "Verify the filename and the file contents.", filename);
    }
//      if ((*minfo).usercharge < -9990)
    if ((*minfo).icharge < -9990)
        (*minfo).icharge = intcharge(atomnum, atom);
}

void rmopcharge(char *filename, int atomnum, ATOM atom[], MOLINFO * minfo)
{

    int index;
    int number = 0;
    int tmpint;
    FILE *fpout;
    char tmpchar[MAXCHAR];
    char tmpchar1[MAXCHAR];
    char tmpchar2[MAXCHAR];
    char tmpchar3[MAXCHAR];
    char tmpchar4[MAXCHAR];
    char tmpchar5[MAXCHAR];
    char tmpchar6[MAXCHAR];
    char line[MAXCHAR];
    double tmpfloat1;
    double tmpfloat2;

    fpout = efopen(filename, "r");
    index = 0;
    for (;;) {
        if (fgets(line, MAXCHAR, fpout) == NULL)
            break;

        sscanf(line, "%s%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3, tmpchar4, tmpchar5,
               tmpchar6);
        if (strcmp("NET", tmpchar1) == 0 && strcmp("ATOMIC", tmpchar2) == 0
            && strcmp("CHARGES", tmpchar3) == 0 && strcmp("AND", tmpchar4) == 0
            && strcmp("DIPOLE", tmpchar5) == 0
            && strcmp("CONTRIBUTIONS", tmpchar6) == 0) {
            index = 1;
            continue;
        }
        if (index == 1 && strcmp("ATOM", tmpchar1) == 0 && strcmp("NO.", tmpchar2) == 0) {
            index = 2;
            continue;
        }
        if (index == 2 && strcmp("DIPOLE", tmpchar1) == 0 && strcmp("X", tmpchar2) == 0
            && strcmp("Y", tmpchar3) == 0 && strcmp("Z", tmpchar4) == 0
            && strcmp("TOTAL", tmpchar5) == 0) {
            index = 3;
            break;
        }
        if (index == 2) {
            sscanf(line, "%d%s%lf%lf", &tmpint, tmpchar, &tmpfloat1, &tmpfloat2);
            atom[number].charge = tmpfloat1;
            number++;
        }
    }
    fclose(fpout);
    if (number == 0) {
        eprintf("Unable to find mopac charges in file (%s).\n"
                "Verify the filename and the file contents.", filename);
    }
//      if ((*minfo).usercharge < -9990)
    if ((*minfo).icharge < -9990)
        (*minfo).icharge = intcharge(atomnum, atom);
}


void rdivcharge(char *filename, int atomnum, ATOM atom[], MOLINFO * minfo, int flag)
{

    /* now modified for divcon.out output */

    int index;
    int number = 0;
    int tmpint;
    FILE *fpout;
    char tmpchar[MAXCHAR];
    char tmpchar1[MAXCHAR];
    char tmpchar2[MAXCHAR];
    char tmpchar3[MAXCHAR];
    char tmpchar4[MAXCHAR];
    char tmpchar5[MAXCHAR];
    char tmpchar6[MAXCHAR];
    char line[MAXCHAR];
    double tmpfloat1;
    double tmpfloat2;
    double tmpfloat3;

    fpout = efopen(filename, "r");
    index = 0;
    for (;;) {
        if (fgets(line, MAXCHAR, fpout) == NULL)
            break;

        sscanf(line, "%s%s%s%s%s%s", tmpchar1, tmpchar2, tmpchar3, tmpchar4, tmpchar5,
               tmpchar6);
        /* parse DIVCON output  */
        if (strcmp("NO.", tmpchar1) == 0 && strcmp("SYMBOL", tmpchar2) == 0
            && strcmp("MULLIKEN", tmpchar3) == 0 && strcmp("CM1", tmpchar4) == 0
            && strcmp("CM2", tmpchar5) == 0) {
            index = 1;
            continue;
        }
        if (index == 1 && strcmp("CHARGE", tmpchar1) == 0
            && strcmp("CHARGE", tmpchar2) == 0) {
            index = 2;
            fgets(line, MAXCHAR, fpout);
            continue;
        }
        if (index == 2 && strcmp("\n", line) == 0) {
            index = 3;
            break;
        }
        if (index == 2) {
            sscanf(line, "%d%s%lf%lf%lf", &tmpint, tmpchar, &tmpfloat1, &tmpfloat2,
                   &tmpfloat3);
            if (flag == 1)
                atom[number].charge = tmpfloat1;
            if (flag == 2)
                atom[number].charge = tmpfloat2;
            if (flag == 3)
                atom[number].charge = tmpfloat3;
            number++;
        }
    }
    fclose(fpout);
    if (number == 0) {
        eprintf("Unable to find divcon charges in file (%s).\n"
                "Verify the filename and the file contents.", filename);
    }
//      if ((*minfo).usercharge < -9990)
    if ((*minfo).icharge < -9990)
        (*minfo).icharge = intcharge(atomnum, atom);
}

void rorccharge(char *filename, int atomnum, ATOM atom[], MOLINFO * minfo)
{
    int numatom = 0;
    int read_flag = 0;
    char line[MAXCHAR];
    FILE *fpin;

    fpin = efopen(filename, "r");
    for (;;) {
        if (!fgets(line, MAXCHAR, fpin))
            eprintf("Premature end of file");
        if (read_flag == 0) {
            if ((strncmp(line,"                 *** FINAL ENERGY", 33) == 0) ||
                (strncmp(line,"                       * Single Point ", 38) == 0)) {
            read_flag = 1;
            continue;
            }
            }
        
        if (read_flag == 1) {
            if (strncmp(line, "MULLIKEN ATOMIC CHARGES", 23) == 0 ) {
 //               printf("READING CHARGES\n");
                read_flag = 2;
            }
            continue;
            }
        if (read_flag == 2) {
            if (line[0] == '-') {
                continue;
            }
            else if (line[0] == 'S') {
                sscanf(line, "%*s%*s%*s%*s%lf", &minfo->dcharge);
                minfo->icharge = (int) minfo->dcharge;
//               printf("%f\n",minfo->dcharge);
//                printf("%f\n",minfo->icharge);
                break;
            }
            sscanf(line, "%*d%*s%*s%lf", &atom[numatom].charge);
//            printf("%d %f\n",numatom, atom[numatom].charge);
            numatom++;
            continue;
        }
    }
    fclose(fpin);
    if (numatom == 0) {
        eprintf("Unable to find Orca charges in file (%s).\n"
                "Verify the filename and the file contents.", filename);
    }
//      if ((*minfo).usercharge < -9990)
    if ((*minfo).icharge < -9990)
        (*minfo).icharge = intcharge(atomnum, atom);
}

void rgaucharge(char *filename, char *chargemethod, int atomnum, ATOM atom[],
                MOLINFO * minfo)
{
    int chargeindex;
    int number = 0;
    int Found_Stationary = 0;
    int index = 0;
    int suc = 0;
    char line[MAXCHAR];
    FILE *fpin;

    fpin = efopen(filename, "r");
    if (strcmp(chargemethod, "mul") == 0)
        chargeindex = 1;
    else
        chargeindex = 2;
    for (;;) {
        if (fgets(line, MAXCHAR, fpin) == NULL)
            break;
        if (strncmp("-- Stationary point found", &line[4], 25) == 0) {
            Found_Stationary = 1;
            continue;
        }
        if (chargeindex == 1 && Found_Stationary == 1 && index == 0) {
            if (strncmp("Total atomic charges:", &line[1], 21) == 0)
                index = 1;
            if (strncmp("Mulliken atomic charges:", &line[1], 24) == 0)
                index = 1;
            continue;
        }
        if (chargeindex == 2 && Found_Stationary == 1 && index == 0
            && strncmp("Fitting point charges", &line[1], 21) == 0)
            index = 1;
        if (index == 1
            && strncmp("Atomic charges with hydrogens summed into heavy atoms", &line[1],
                       52) == 0)
            index = 2;
        if (Found_Stationary == 1 && index == 1) {
            if (line[11] == '.') {
                sscanf(&line[8], "%lf", &atom[number].charge);
		suc = 1;
                number++;
            }
            if (line[14] == '.') {
                sscanf(&line[11], "%lf", &atom[number].charge);
		suc = 1;
                number++;
            }
        }
        if (number > atomnum)
            break;
    }
    if(suc == 0) {
	rewind(fpin);
	number = 0;
    	for (;;) {
        	if (fgets(line, MAXCHAR, fpin) == NULL)
            		break;
        	if (chargeindex == 1 && index == 0) {
            		if (strncmp("Total atomic charges:", &line[1], 21) == 0)
                		index = 1;
            		if (strncmp("Mulliken atomic charges:", &line[1], 24) == 0)
                		index = 1;
            		if (strncmp("Mulliken charges:", &line[1], 16) == 0)
                		index = 1;
			continue;
        	}
        	if (chargeindex == 2 && index == 0
            		&& strncmp("Fitting point charges", &line[1], 21) == 0)
            		index = 1;
        	if (index == 1 && strncmp("Atomic charges with hydrogens summed into heavy atoms", &line[1], 52) == 0)
            		index = 2;
        	if (index == 1 && strncmp("Sum of Mulliken charges", &line[1], 23) == 0)
            		index = 2;
        	if (index == 1) {
            		if (line[11] == '.') {
                		sscanf(&line[8], "%lf", &atom[number].charge);
                		number++;
            		}
            		if (line[14] == '.') {
                		sscanf(&line[11], "%lf", &atom[number].charge);
                		suc = 1;
                		number++;
            		}
        	}
        	if (number > atomnum)
            		break;
    	}
    }
    fclose(fpin);
    if (number == 0) {
        eprintf("Unable to find gaussian charges in file (%s).\n"
                "Verify the filename and the file contents.", filename);
    }
//      if ((*minfo).usercharge < -9990)
    if ((*minfo).icharge < -9990)
        (*minfo).icharge = intcharge(atomnum, atom);
}

void rjoutcharge(char *filename, char *chargemethod, int atomnum,
                                ATOM atom[], MOLINFO *minfo)
{
        int num;
	int pos = -1;
	int i;
        char line[MAXCHAR];
	char tmpc1[MAXCHAR], tmpc2[MAXCHAR];
	double c[5];
        FILE *fpin;

        if ((fpin = fopen(filename, "r")) == NULL) {
                fprintf(stdout, "Cannot open the jaguar output file: %s in rjoutcharge(), exit\n", filename);
                return;
        }
/*
        int chargeindex;
        if (strcmp(chargemethod, "mul") == 0)
                chargeindex = 1;
        else
                chargeindex = 2;
*/
        for (;;) {
                if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
		sscanf(line, "%s%s", tmpc1, tmpc2);
		if(strcmp(tmpc1, "Atomic") == 0 && strcmp(tmpc2, "charges") == 0) 
			pos = ftell(fpin);
	}		
	if(pos < 0) {
        	fclose(fpin);
		return;
	} 
        num = 0;
        fseek(fpin, pos, 0);
        for (;;) {
                if (fgets(line, MAXCHAR, fpin) == NULL)
                        break;
		strcpy(tmpc1, "");
		strcpy(tmpc2, "");
		tmpc1[0]='\0';
		tmpc2[0]='\0';
		sscanf(line, "%s%s", tmpc1, tmpc2);
		if(strcmp(tmpc1, "Charge") == 0) {
			sscanf(&line[7],"%lf%lf%lf%lf%lf", &c[0], &c[1], &c[2], &c[3], &c[4]);
			for(i=0; i< 5; i++) {
				atom[num].charge = c[i];
				num++;
        			if (num > atomnum) 
					break;
			}
		}
		if(strcmp(tmpc1, "sum") == 0 && strcmp(tmpc2, "of") == 0) 
			break;
	}		
	
        fclose(fpin);
/*      if ((*minfo).usercharge < -9990) */
        if ((*minfo).icharge < -9990)
       	    (*minfo).icharge = intcharge(atomnum, atom);
}


/*CHARGE METHOD : READ CHARGE */
void rcharge(char *filename, int atomnum, ATOM atom[], CONTROLINFO cinfo, MOLINFO * minfo)
{
    FILE *fpcharge;
    char line[MAXCHAR];
    int i, j;
    int flag;
    int break_flag = 0;
    int ncol = 0;
    int number = 0;
    double chg[10];

    fpcharge = efopen(filename, "r");
    i = 0;
    flag = 0;
    for (;;) {
        if (fgets(line, MAXCHAR, fpcharge) == NULL) {
            if (cinfo.intstatus == 2)
                printf("\nInfo: Finished reading charge file (%s).", filename);
            break;
        }
        sscanf(line, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &chg[0], &chg[1], &chg[2], &chg[3],
               &chg[4], &chg[5], &chg[6], &chg[7], &chg[8], &chg[9]);
//      check how many columns per line
        if (flag == 0) {
            for (j = 0; j < strlen(line); j++)
                if (line[j] == '.')
                    ncol++;
            flag = 1;
        }
        for (j = 0; j < ncol; j++) {
            if (i == atomnum - j) {
                break_flag = 1;
                break;
            }
            atom[i + j].charge = chg[j];
            number++;
        }
        if (break_flag == 1)
            break;
        else
            i = i + ncol;
    }
    fclose(fpcharge);
    if (number == 0) {
        eprintf("\nUnable to read charges in the -cf file (%s).\n"
                "Using charges from the -i input file.", filename);
    }
//      if ((*minfo).usercharge < -9990.)
    if ((*minfo).icharge < -9990.)
        (*minfo).icharge = intcharge(atomnum, atom);
}

/*CHARGE METHOD : WRITE CHARGE */

void wcharge(char *filename, int atomnum, ATOM atom[], CONTROLINFO cinfo, MOLINFO minfo)
{
    FILE *fpcharge;
    int i;
    int num = 0;

    fpcharge = efopen(filename, "w");
    for (i = 0; i < atomnum; i++) {
        fprintf(fpcharge, "%10.6lf", atom[i].charge);
        num++;
        if (num == 8) {
            num = 0;
            fprintf(fpcharge, "\n");
        }
    }
    fclose(fpcharge);
}



/*CHARGE METHOD : RESP CHARGE */

void resp(char *filename, int atomnum, ATOM * atom, int bondnum, BOND * bond,
          CONTROLINFO cinfo, MOLINFO minfo)
{
    char tmpchar1[MAXCHAR];
    char tmpchar2[MAXCHAR];
    char tmpchar3[MAXCHAR];
    char tmpchar4[MAXCHAR];
    char tmpchar5[MAXCHAR];
    char str_eqcharge[20];
    char str_max_path_len[20];
    size_t copied_size;
    int flag;
    int use_max_path_len = 0;
    char vdwfile[MAXCHAR]="";

    flag = 0;
    if (strcmp(cinfo.intype, "gout") == 0 || strcmp(cinfo.intype, "11") == 0)
        flag = 1;
    if (strcmp(cinfo.intype, "gesp") == 0 || strcmp(cinfo.intype, "25") == 0)
        flag = 2;
    if (strcmp(cinfo.intype, "gamess") == 0 || strcmp(cinfo.intype, "27") == 0)
        flag = 3;
    if (strcmp(cinfo.intype, "quick") == 0 || strcmp(cinfo.intype, "32") == 0) {
        strncpy(vdwfile,filename,strlen(filename)-3);
        strcat(vdwfile,"vdw");
//        strncpy(filename,tmpfile,strlen(tmpfile));
        flag = 4;
    };

    if (flag == 0) {
        eprintf ("Cannot generate RESP charges.\n"
                  "Invalid input selections: The RESP charge method requires\n"
                  "a Gaussian output file, i.e. -fi gout -gv 0\n"
                  "or a Gaussian ESP file, i.e. -fi gesp -gv 1\n"
                  "or a GAMESS dat file, i.e. -fi gamess.\n"
                  "or a QUICK output file, i.e. -fi quick.\n");
        return;
    }
    wac("ANTECHAMBER_RESP.AC", atomnum, atom, bondnum, bond, cinfo, minfo);

    sprintf(str_eqcharge, "%d", minfo.eqcharge);
    if (cinfo.max_path_length > 0) {
        sprintf(str_max_path_len, "%d", cinfo.max_path_length);
        use_max_path_len = 1;
    }
    copied_size = build_exe_path(tmpchar1, "espgen", sizeof tmpchar1, 1);
    strncat(tmpchar1, " -o ANTECHAMBER.ESP -i ", sizeof tmpchar1 - copied_size);
    strncat(tmpchar1, filename, sizeof(tmpchar1) - strlen(tmpchar1) - 1);
    /* if quick output file is use flag is set to 4 */
    if(flag == 4) {
      strncat(tmpchar1, " -f 4 -vdw ", 11);
      strncat(tmpchar1, vdwfile, sizeof(tmpchar1) - strlen(tmpchar1) - 1);
    }

    copied_size = build_exe_path(tmpchar2, "respgen", sizeof tmpchar2, 1);
    strncat(tmpchar2, " -i ANTECHAMBER_RESP.AC -o ANTECHAMBER_RESP1.IN" " -f resp1 -e ",
            sizeof tmpchar2 - copied_size);
    strcat(tmpchar2, str_eqcharge);
    if (use_max_path_len == 1) {
        strcat(tmpchar2, " -l ");
        strcat(tmpchar2, str_max_path_len);
    }

    copied_size = build_exe_path(tmpchar3, "respgen", sizeof tmpchar3, 1);
    strncat(tmpchar3, " -i ANTECHAMBER_RESP.AC -o ANTECHAMBER_RESP2.IN" " -f resp2 -e ",
            sizeof tmpchar3 - copied_size);
    strcat(tmpchar3, str_eqcharge);
    if (use_max_path_len == 1) {
        strcat(tmpchar3, " -l ");
        strcat(tmpchar3, str_max_path_len);
    }

    copied_size = build_exe_path(tmpchar4, "resp", sizeof tmpchar4, 1);

    copied_size = build_exe_path(tmpchar5, "resp", sizeof tmpchar5, 1);

    /*  run espgen:  */
    if (cinfo.intstatus == 2)
        fprintf(stdout, "\nRunning: %s\n\n", tmpchar1);
    esystem(tmpchar1);

    /*  run  respgen to make ANTECHAMBER_RESP1.IN:   */
    if (cinfo.intstatus == 2)
        fprintf(stdout, "\nRunning: %s\n\n", tmpchar2);
    esystem(tmpchar2);

    /*  run  respgen to make ANTECHAMBER_RESP2.IN:   */
    if (cinfo.intstatus == 2)
        fprintf(stdout, "\nRunning: %s\n\n", tmpchar3);
    esystem(tmpchar3);

    /* run first stage of resp:  */
    if(flag == 4) {
        if(cinfo.sfac[0] == '\0'){
            strcat(tmpchar4,
               " -C -O -i ANTECHAMBER_RESP1.IN -o ANTECHAMBER_RESP1.OUT -e ANTECHAMBER.ESP -t qout");
        }else{
            strcat(tmpchar4,
               " -C -O -i ANTECHAMBER_RESP1.IN -o ANTECHAMBER_RESP1.OUT -e ANTECHAMBER.ESP -t qout -f ");
            strcat(tmpchar4,cinfo.sfac);
        }
    }else{
        strcat(tmpchar4,
           " -O -i ANTECHAMBER_RESP1.IN -o ANTECHAMBER_RESP1.OUT -e ANTECHAMBER.ESP -t qout");
    }
    if (cinfo.intstatus == 2)
        fprintf(stdout, "\nRunning: %s\n\n", tmpchar4);
    esystem(tmpchar4);

    /* run second stage of resp:  */
    if(flag == 4) {
        if(cinfo.sfac[0] == '\0'){
            strcat(tmpchar5,
               " -C -O -i ANTECHAMBER_RESP2.IN -o ANTECHAMBER_RESP2.OUT -e ANTECHAMBER.ESP -q qout -t QOUT");
        }else {
            strcat(tmpchar5,
               " -C -O -i ANTECHAMBER_RESP2.IN -o ANTECHAMBER_RESP2.OUT -e ANTECHAMBER.ESP -q qout -t QOUT -f ");
            strcat(tmpchar5,cinfo.sfac);
        }
    }else{
        strcat(tmpchar5,
           " -O -i ANTECHAMBER_RESP2.IN -o ANTECHAMBER_RESP2.OUT -e ANTECHAMBER.ESP -q qout -t QOUT");
    }
    if (cinfo.intstatus == 2)
        fprintf(stdout, "\nRunning: %s\n\n", tmpchar5);
    esystem(tmpchar5);
    rcharge("QOUT", atomnum, atom, cinfo, &minfo);
}

/* CHARGE METHOD : BCC-AM1 */
void bccharge(int atomnum, ATOM atom[], int bondnum, BOND bond[], AROM arom[],
              CONTROLINFO * cinfo, MOLINFO * minfo)
{
    int i, j;
    int atomnum_tmp;
    int bondnum_tmp;
    ATOM *atom_tmp;
    BOND *bond_tmp;
    int *equ_atom_id;
    char tmpchar[MAXCHAR];
    char tmpchar1[MAXCHAR];
    char tmpchar2[MAXCHAR];
    size_t copied_size;

    int nequatom;
    double tcharge;

//      find equal atoms and force their charge same
    if ((*minfo).eqcharge == 1 || (*minfo).eqcharge == 2) {
        equ_atom_id = (int *) emalloc(sizeof(int) * atomnum);
        if ((*minfo).eqcharge == 1)
            identify_equatom(atomnum, atom, equ_atom_id, (*cinfo).max_path_length,
                             bondnum, bond, 0);
        if ((*minfo).eqcharge == 2)
            identify_equatom(atomnum, atom, equ_atom_id, (*cinfo).max_path_length,
                             bondnum, bond, 1);
        for (i = 0; i < atomnum - 1; i++) {
            if (equ_atom_id[i] != -1)
                continue;
            nequatom = 1;
            tcharge = atom[i].charge;
            for (j = i + 1; j < atomnum; j++) {
                if (equ_atom_id[j] == i) {
                    nequatom++;
                    tcharge += atom[j].charge;
                }
            }
            if (nequatom == 1)
                continue;
            atom[i].charge = tcharge / nequatom;
        }
        for (i = 0; i < atomnum; i++) {
            if (equ_atom_id[i] == -1)
                continue;
            atom[i].charge = atom[equ_atom_id[i]].charge;
        }
    }
//      done
    wac("ANTECHAMBER_AM1BCC_PRE.AC", atomnum, atom, bondnum, bond, *cinfo, *minfo);
    copied_size = build_exe_path(tmpchar1, "am1bcc", sizeof tmpchar1, 1);
    if (strcmp((*cinfo).chargetype, "abcg2") == 0 || strcmp((*cinfo).chargetype, "ABCG2") == 0
    	|| strcmp((*cinfo).chargetype, "8") == 0) {
    	strncat(tmpchar1,
            " -i ANTECHAMBER_AM1BCC_PRE.AC -o ANTECHAMBER_AM1BCC.AC" " -f ac -t abcg2 -p ",
            sizeof tmpchar1 - copied_size);
    	build_dat_path(tmpchar2, "BCCPARM_ABCG2.DAT", sizeof tmpchar2, 1);
    }
    else {
    	strncat(tmpchar1,
            " -i ANTECHAMBER_AM1BCC_PRE.AC -o ANTECHAMBER_AM1BCC.AC" " -f ac -t bcc -p ",
            sizeof tmpchar1 - copied_size);
    	build_dat_path(tmpchar2, "BCCPARM.DAT", sizeof tmpchar2, 1);
    }
    strncat(tmpchar1, tmpchar2, sizeof(tmpchar1) - strlen(tmpchar1) - 1);
    strcat(tmpchar1, " -s ");
    sprintf(tmpchar2, "%d", (*cinfo).intstatus);
/*      newitoa((*cinfo).intstatus, tmpchar2); */
    strcat(tmpchar1, tmpchar2);
/*
	if (strcmp((*cinfo).intype, "ac") == 0 || strcmp((*cinfo).intype, "1") == 0
		|| strcmp((*cinfo).intype, "mol2") == 0
		|| strcmp((*cinfo).intype, "2") == 0) {
		if ((*cinfo).prediction_index == -1)
			(*cinfo).prediction_index = 0;
		sprintf(tmpchar, "%d", (*cinfo).prediction_index); 
		strcat(tmpchar1, " -j ");
		strcat(tmpchar1, tmpchar);
	} else
		strcat(tmpchar1, " -j 4");
*/
    if ((*cinfo).prediction_index == -1)
        (*cinfo).prediction_index = 0;
/*	newitoa((*cinfo).prediction_index, tmpchar); */
    sprintf(tmpchar, "%d", (*cinfo).prediction_index);
    strcat(tmpchar1, " -j ");
    strcat(tmpchar1, tmpchar);

    if ((*cinfo).intstatus == 2)
        fprintf(stdout, "\nRunning: %s\n", tmpchar1);
    esystem(tmpchar1);
    (*cinfo).maxatom = atomnum + 10;
    (*cinfo).maxbond = bondnum + 10;
    atom_tmp = (ATOM *) emalloc(sizeof(ATOM) * (*cinfo).maxatom);
    bond_tmp = (BOND *) emalloc(sizeof(BOND) * (*cinfo).maxbond);
    for (i = 0; i < (*cinfo).maxbond; ++i) {
        bond_tmp[i].jflag = -1; /* bond type has not been assigned */
    }
    minfo2 = *minfo;
    rac("ANTECHAMBER_AM1BCC.AC", &atomnum_tmp, atom_tmp, &bondnum_tmp, bond_tmp, cinfo,
        &minfo2);

    for (i = 0; i < atomnum; i++)
        atom[i].charge = atom_tmp[i].charge;
    free(atom_tmp);
    free(bond_tmp);

}

/* CHARGE METHOD : BCC-AM1 */
void bcc(char *filename, int atomnum, ATOM atom[], int bondnum, BOND bond[], AROM arom[],
         CONTROLINFO * cinfo, MOLINFO * minfo)
{
    char tmpchar[MAXCHAR];

    if (strcmp((*cinfo).intype, "mopout") == 0 || strcmp((*cinfo).intype, "12") == 0)
        rmopcharge(filename, atomnum, atom, minfo);
    else if (strcmp((*cinfo).intype, "divout") == 0 || strcmp((*cinfo).intype, "22") == 0)
        rdivcharge(filename, atomnum, atom, minfo, 1);
    else if (strcmp((*cinfo).intype, "sqmout") == 0 || strcmp((*cinfo).intype, "24") == 0)
        rsqmcharge(filename, atomnum, atom, minfo);
    else {
        if ((*minfo).divcon == 0) {
            wmopcrt("mopac.in", atomnum, atom, *minfo);
            build_exe_path(tmpchar, "mopac.sh", sizeof tmpchar, 1);
        }
        if ((*minfo).divcon == 1) {
            wdivcrt("divcon.in", atomnum, atom, *minfo);
            build_exe_path(tmpchar, "divcon", sizeof tmpchar, 1);
        }
        if ((*minfo).divcon == 2) {
            wsqmcrt("sqm.in", atomnum, atom, *minfo);
            build_exe_path(tmpchar, "sqm", sizeof tmpchar, 1);
            strcat(tmpchar, " -O -i sqm.in -o sqm.out");
        }
        if ((*cinfo).intstatus == 2 || (*cinfo).intstatus == 1)
            fprintf(stdout, "\nRunning: %s\n", tmpchar);
        esystem(tmpchar);
        if ((*minfo).divcon == 0) {
            rmopcharge("mopac.out", atomnum, atom, minfo);
            wpdb_optimized("mopac.out", atomnum, atom, 0);
        }
        if ((*minfo).divcon == 1) {
            rdivcharge("divcon.out", atomnum, atom, minfo, 1);
            wpdb_optimized("divcon.out", atomnum, atom, 1);
        }
        if ((*minfo).divcon == 2) {
            rsqmcharge("sqm.out", atomnum, atom, minfo);
            wpdb_optimized("sqm.out", atomnum, atom, 2);
        }
    }
    bccharge(atomnum, atom, bondnum, bond, arom, cinfo, minfo);
}

/* CHARGE METHOD : CM1 */
void cm1(int atomnum, ATOM atom[], CONTROLINFO * cinfo, MOLINFO * minfo)
{
    char tmpchar[MAXCHAR];

    if ((*minfo).divcon == 0) {
        eprintf("The mopac program cannot generate CM1 charges.");
        return;
    }
    if ((*minfo).divcon == 1) {
        wdivcrt("divcon.in", atomnum, atom, *minfo);
        build_exe_path(tmpchar, "divcon", sizeof tmpchar, 1);
    }
    if ((*minfo).divcon == 2) {
        eprintf("The sqm program cannot generate CM1 charges.");
        return;
    }
    if ((*cinfo).intstatus == 2 || (*cinfo).intstatus == 1)
        fprintf(stdout, "\nRunning: %s\n", tmpchar);
    esystem(tmpchar);
    rdivcharge("divcon.out", atomnum, atom, minfo, 2);
    wpdb_optimized("divcon.out", atomnum, atom, 1);
}

/* CHARGE METHOD : CM2 */
void cm2(int atomnum, ATOM atom[], CONTROLINFO * cinfo, MOLINFO * minfo)
{
    char tmpchar[MAXCHAR];

    if ((*minfo).divcon == 0) {
        eprintf("The mopac program cannot generate CM2 charges.");
        return;
    }
    if ((*minfo).divcon == 1) {
        wdivcrt("divcon.in", atomnum, atom, *minfo);
        build_exe_path(tmpchar, "divcon", sizeof tmpchar, 1);
    }
    if ((*minfo).divcon == 2) {
        eprintf("The sqm program cannot generate CM2 charges.");
        return;
    }
    if ((*cinfo).intstatus == 2 || (*cinfo).intstatus == 1)
        fprintf(stdout, "\nRunning: %s\n", tmpchar);
    esystem(tmpchar);
    rdivcharge("divcon.out", atomnum, atom, minfo, 3);
    wpdb_optimized("divcon.out", atomnum, atom, 1);
}

/* CHARGE METHOD : Mulliken  */
void mul(char *filename, int atomnum, ATOM atom[], CONTROLINFO * cinfo, MOLINFO * minfo)
{
    char tmpchar[MAXCHAR];

    if (strcmp((*cinfo).intype, "gout") == 0 || strcmp((*cinfo).intype, "11") == 0)
        rgaucharge(filename, "mul", atomnum, atom, minfo);
    else if (strcmp((*cinfo).intype, "orcout") == 0 || strcmp((*cinfo).intype, "29") == 0)
        rorccharge(filename, atomnum, atom, minfo);
    else if (strcmp((*cinfo).intype, "mopout") == 0 || strcmp((*cinfo).intype, "12") == 0)
        rmopcharge(filename, atomnum, atom, minfo);
    else if (strcmp((*cinfo).intype, "divout") == 0 || strcmp((*cinfo).intype, "22") == 0)
        rdivcharge(filename, atomnum, atom, minfo, 1);
    else if (strcmp((*cinfo).intype, "sqmout") == 0 || strcmp((*cinfo).intype, "24") == 0)
        rsqmcharge(filename, atomnum, atom, minfo);
    else {
        if ((*minfo).divcon == 0) {
            wmopcrt("mopac.in", atomnum, atom, *minfo);
            build_exe_path(tmpchar, "mopac.sh", sizeof tmpchar, 1);
        }
        if ((*minfo).divcon == 1) {
            wdivcrt("divcon.in", atomnum, atom, *minfo);
            build_exe_path(tmpchar, "divcon", sizeof tmpchar, 1);
        }
        if ((*minfo).divcon == 2) {
            wsqmcrt("sqm.in", atomnum, atom, *minfo);
            build_exe_path(tmpchar, "sqm", sizeof tmpchar, 1);
            strcat(tmpchar, " -O -i sqm.in -o sqm.out");
        }
        if ((*cinfo).intstatus == 2 || (*cinfo).intstatus == 1)
            fprintf(stdout, "\nRunning: %s\n", tmpchar);
        esystem(tmpchar);
        if ((*minfo).divcon == 0) {
            rmopcharge("mopac.out", atomnum, atom, minfo);
            wpdb_optimized("mopac.out", atomnum, atom, 0);
        }
        if ((*minfo).divcon == 1) {
            rdivcharge("divcon.out", atomnum, atom, minfo, 1);
            wpdb_optimized("divcon.out", atomnum, atom, 1);
        }
        if ((*minfo).divcon == 2) {
            rsqmcharge("sqm.out", atomnum, atom, minfo);
            wpdb_optimized("sqm.out", atomnum, atom, 2);
        }
    }
}

/* CHARGE METHOD : ESP(Kollman) */
void esp(char *filename, int atomnum, ATOM atom[], CONTROLINFO cinfo, MOLINFO minfo)
{
    if (strcmp(cinfo.intype, "gout") == 0)
        rgaucharge(filename, "esp", atomnum, atom, &minfo);
    else if (strcmp(cinfo.intype, "jout") == 0)
        rjoutcharge(filename, "esp", atomnum, atom, &minfo);
    else
        eprintf ("Cannot generate Kollman ESP charges.\n"
                  "Invalid input selections: The Kollman ESP charge method requires\n"
                  "a Gaussian output file, i.e. -fi gout.");
}

/* CHARGE METHOD : Gasteiger charge */
void rgasparm(char *filename, int *gasparmnum, GASTEIGER gas[])
{
    FILE *fp;
    int num = 0;
    char line[MAXCHAR];

    fp = efopen(filename, "r");
    for (;;) {
        if (fgets(line, MAXCHAR, fp) == NULL)
            break;
        if (strncmp(line, "GASPARM", 7) == 0) {
            gas[num].name[0] = '\0';
            sscanf(&line[8], "%s%lf%lf%lf%lf%lf", gas[num].name, &gas[num].a, &gas[num].b,
                   &gas[num].c, &gas[num].d, &gas[num].charge);
            num++;
        }
    }
    *gasparmnum = num;
    fclose(fp);
}

void assign(int atomnum, ATOM atom[], int gasparmnum, GASTEIGER gas[], int gasparmindex[])
{
    int i, j;
    int flag;
    for (i = 0; i < atomnum; i++) {
        flag = 0;
        for (j = 0; j < gasparmnum; j++)
            if (strcmp(gas[j].name, atom[i].ambername) == 0) {
                gasparmindex[i] = j;
                flag = 1;
                break;
            }
        if (flag == 0) {
            eprintf("No Gasteiger parameter for atom (ID: %d, Name: %s, Type: %s).",
                    i, atom[i].name, atom[i].ambername);
        }
    }
}

double rmscal(int atomnum, double gaschargep[], double gaschargea[])
{
    double total = 0.0;
    double rms;
    int i;

    for (i = 0; i < atomnum; i++) {
        total += (gaschargep[i] - gaschargea[i]) * (gaschargep[i] - gaschargea[i]);
        gaschargep[i] = gaschargea[i];
    }
    rms = total / atomnum;
    rms = sqrt(rms);
    return rms;
}

void gasiter(int atomnum, ATOM atom[], int gasparmnum, GASTEIGER gas[], int netcharge)
{
    int i, j;
    int iteration = 0;
    double q;
    double xx;
    double rmsd;
    int *gasparmindex;
    double *x;
    double *gaschargep;
    double *gaschargea;
    double TotalGasCharge = 0.0;

    gasparmindex = (int *) emalloc(sizeof(int) * atomnum);
    x = (double *) emalloc(sizeof(double) * atomnum);
    gaschargep = (double *) emalloc(sizeof(double) * atomnum);
    gaschargea = (double *) emalloc(sizeof(double) * atomnum);
    assign(atomnum, atom, gasparmnum, gas, gasparmindex);
    for (i = 0; i < atomnum; i++) {
        gaschargep[i] = gas[gasparmindex[i]].charge;
        gaschargea[i] = gas[gasparmindex[i]].charge;
        TotalGasCharge += gas[gasparmindex[i]].charge;
    }
    // TODO check this - fp vs int ?
    if (TotalGasCharge != netcharge) {
        printf("\nWarning: The net charge of the molecule (%-5d) does not equal the"
                 "         total formal charge (%-5.2lf) according to the"
                 " Gasteiger atom types.\n", netcharge, TotalGasCharge);
    }
    do {
        for (i = 0; i < atomnum; i++) {
            x[i] = gas[gasparmindex[i]].a + gas[gasparmindex[i]].b * gaschargep[i];
            x[i] += gas[gasparmindex[i]].c * gaschargep[i] * gaschargep[i];
            if (x[i] == 0.0)
                x[i] = 0.0000000001;
        }
        for (i = 0; i < atomnum; i++)
            for (j = i + 1; j < atomnum; j++)
                if (atom[i].con[0] == j || atom[i].con[1] == j || atom[i].con[2] == j
                    || atom[i].con[3] == j || atom[i].con[4] == j
                    || atom[i].con[5] == j) {
                    if (x[i] <= x[j]) {
                        xx = gas[gasparmindex[i]].d;
                        q = (x[j] - x[i]) / xx * pow(DAMPFACTOR, iteration + 1);
                        gaschargea[i] += q;
                        gaschargea[j] -= q;
                    }
                    if (x[i] > x[j]) {
                        xx = gas[gasparmindex[j]].d;
                        q = (x[i] - x[j]) / xx * pow(DAMPFACTOR, iteration + 1);
                        gaschargea[i] -= q;
                        gaschargea[j] += q;
                    }
                }
        iteration++;
/*
		printf("\nIteration %5d", iteration);
		for (i = 0; i < atomnum; i++) 
			printf("\n%5d %5s %8.4lf %8.4lf", i+1, atom[i].name, gaschargep[i], gaschargea[i]);
*/
        rmsd = rmscal(atomnum, gaschargep, gaschargea);
    } while (rmsd > CONVERG && iteration < GASMAXITER);

    for (i = 0; i < atomnum; i++)
        atom[i].charge = gaschargea[i];
    free(gasparmindex);
    free(x);
    free(gaschargep);
    free(gaschargea);
}

void gascharge(int atomnum, ATOM atom[], int bondnum, BOND bond[], CONTROLINFO cinfo,
               MOLINFO * minfo)
{
    int i;
    int atomnum_tmp;
    int bondnum_tmp;
    ATOM *atom_tmp;
    BOND *bond_tmp;

    int gasparmnum = 0;
    char tmpchar1[MAXCHAR];
    char tmpchar2[MAXCHAR];
    size_t copied_size;
    GASTEIGER *gas;

    gas = (GASTEIGER *) emalloc(sizeof(GASTEIGER) * MAXGAS);
    build_dat_path((*minfo).gfilename, "GASPARM.DAT", sizeof(*minfo).gfilename, 0);
    rgasparm((*minfo).gfilename, &gasparmnum, gas);     /*principle parameter file */
    wac("ANTECHAMBER_GAS.AC", atomnum, atom, bondnum, bond, cinfo, (*minfo));
    copied_size = build_exe_path(tmpchar1, "atomtype", sizeof tmpchar1, 1);
    strncat(tmpchar1, " -i ANTECHAMBER_GAS.AC -o ANTECHAMBER_GAS_AT.AC -d ",
            sizeof tmpchar1 - copied_size);
    build_dat_path(tmpchar2, "ATOMTYPE_GAS.DEF", sizeof tmpchar2, 1);
    strncat(tmpchar1, tmpchar2, sizeof(tmpchar1) - strlen(tmpchar1) - 1);

    if (cinfo.intstatus == 2)
        fprintf(stdout, "Running: %s\n", tmpchar1);
    esystem(tmpchar1);
    cinfo.maxatom = atomnum + 10;
    cinfo.maxbond = bondnum + 10;
    atom_tmp = (ATOM *) emalloc(sizeof(ATOM) * cinfo.maxatom);
    bond_tmp = (BOND *) emalloc(sizeof(BOND) * cinfo.maxbond);
    for (i = 0; i < cinfo.maxbond; ++i) {
        bond_tmp[i].jflag = -1; /* bond type has not been assigned */
    }

    minfo2 = *minfo;
    rac("ANTECHAMBER_GAS_AT.AC", &atomnum_tmp, atom_tmp, &bondnum_tmp, bond_tmp, &cinfo,
        &minfo2);
//      gasiter(atomnum_tmp, atom_tmp, gasparmnum, gas, (*minfo).usercharge);
    gasiter(atomnum_tmp, atom_tmp, gasparmnum, gas, (*minfo).icharge);
    for (i = 0; i < atomnum; i++)
        atom[i].charge = atom_tmp[i].charge;
    free(atom_tmp);
    free(bond_tmp);
    free(gas);
}



void write_sybyl_bat(char *str)
{
    FILE *fpout;
    char tmpchar[MAXCHAR];

    fpout = efopen("antechamber_sybyl.bat", "w");
    build_dat_path(tmpchar, "charge.spl", sizeof tmpchar, 0);

    fprintf(fpout, "%s", "#!/bin/csh");
    fprintf(fpout, "\n%s", "sybyl << @");
    fprintf(fpout, "\n%s", "take ");
    fprintf(fpout, "%s", tmpchar);
    if (strcmp(str, "gas1") == 0)
        fprintf(fpout, "\n%s", "CHARGE GASTEIGER");
    if (strcmp(str, "gas2") == 0)
        fprintf(fpout, "\n%s", "CHARGE GAST_HUCK");
    if (strcmp(str, "del") == 0)
        fprintf(fpout, "\n%s", "CHARGE DELRE");
    if (strcmp(str, "pull") == 0)
        fprintf(fpout, "\n%s", "CHARGE PULLMAN ");
    if (strcmp(str, "huc") == 0)
        fprintf(fpout, "\n%s", "CHARGE HUCKEL");
    if (strcmp(str, "mmff") == 0)
        fprintf(fpout, "\n%s", "CHARGE MMFF94");
    fprintf(fpout, "\n%s", "QUIT YES");
    fprintf(fpout, "\n%s\n\n", "@");
    fclose(fpout);
    esystem("chmod +x antechamber_sybyl.bat");
    esystem("antechamber_sybyl.bat");
}


/* CHARGE METHOD : Gasteiger, using sybyl */
void gas1(int atomnum, ATOM atom[], int bondnum, BOND bond[], AROM arom[],
          CONTROLINFO cinfo, MOLINFO minfo)
{
    wmol2("ANTECHAMBER_INPUT.mol2", atomnum, atom, bondnum, bond, arom, cinfo, minfo);
    write_sybyl_bat("gas1");
    rmol2("ANTECHAMBER_INPUT.mol2", &atomnum, atom, &bondnum, bond, &cinfo, &minfo, 0);
}

/* CHARGE METHOD : Del Re, using sybyl*/
void del(int atomnum, ATOM atom[], int bondnum, BOND bond[], AROM arom[],
         CONTROLINFO cinfo, MOLINFO minfo)
{
    wmol2("ANTECHAMBER_INPUT.mol2", atomnum, atom, bondnum, bond, arom, cinfo, minfo);
    write_sybyl_bat("del");
    rmol2("ANTECHAMBER_INPUT.mol2", &atomnum, atom, &bondnum, bond, &cinfo, &minfo, 0);
}

/* CHARGE METHOD : Pullman, using sybyl */
void pull(int atomnum, ATOM atom[], int bondnum, BOND bond[], AROM arom[],
          CONTROLINFO cinfo, MOLINFO minfo)
{
    wmol2("ANTECHAMBER_INPUT.mol2", atomnum, atom, bondnum, bond, arom, cinfo, minfo);
    write_sybyl_bat("pull");
    rmol2("ANTECHAMBER_INPUT.mol2", &atomnum, atom, &bondnum, bond, &cinfo, &minfo, 0);
}

/* CHARGE METHOD :Gasteiger-Huckel, using sybyl*/
void gas2(int atomnum, ATOM atom[], int bondnum, BOND bond[], AROM arom[],
          CONTROLINFO cinfo, MOLINFO minfo)
{
    wmol2("ANTECHAMBER_INPUT.mol2", atomnum, atom, bondnum, bond, arom, cinfo, minfo);
    write_sybyl_bat("gas2");
    rmol2("ANTECHAMBER_INPUT.mol2", &atomnum, atom, &bondnum, bond, &cinfo, &minfo, 0);
}

/* CHARGE METHOD : Huckel, using sybyl */
void huc(int atomnum, ATOM atom[], int bondnum, BOND bond[], AROM arom[],
         CONTROLINFO cinfo, MOLINFO minfo)
{

    wmol2("ANTECHAMBER_INPUT.mol2", atomnum, atom, bondnum, bond, arom, cinfo, minfo);
    write_sybyl_bat("huc");
    rmol2("ANTECHAMBER_INPUT.mol2", &atomnum, atom, &bondnum, bond, &cinfo, &minfo, 0);
}

/* CHARGE METHOD :MMFF94, using sybyl*/
void mmff(int atomnum, ATOM atom[], int bondnum, BOND bond[], AROM arom[],
          CONTROLINFO cinfo, MOLINFO minfo)
{
    wmol2("ANTECHAMBER_INPUT.mol2", atomnum, atom, bondnum, bond, arom, cinfo, minfo);
    write_sybyl_bat("mmff");
    rmol2("ANTECHAMBER_INPUT.mol2", &atomnum, atom, &bondnum, bond, &cinfo, &minfo, 0);
}
