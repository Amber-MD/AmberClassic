#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "topRead.h"
#include "vector.h"
#include "matrix.h"
#include "grid.h"
#include "crdmanip.h"
#include "macros.h"
#define PI                     3.1415926535898

/***=======================================================================***/
/*** GetPrmTop: load an AMBER 7 prmtop topology file and return it as a    ***/
/***            prmtop data structure.                                     ***/
/***=======================================================================***/
void GetPrmTop(prmtop * tp, int verbosity, int adjbnd)
{
    int i, natm;
    FILE *inp;

  /*** Open the source file ***/
    if ((inp = fopen(tp->source, "r")) == NULL) {
        printf("GetPrmTop >> Error.  AMBER prmtop file %s not found!\n",
               tp->source);
        exit(1);
    }

  /*** Version stamp ***/
    fgets(tp->vstamp, 128, inp);

  /*** First, the pointer information ***/
    ReadToFlag("POINTERS", inp);
    fscanf(inp,
           "%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d",
           &tp->natom, &tp->ntypes, &tp->withH.nbond, &tp->woH.nbond,
           &tp->withH.nangl, &tp->woH.nangl, &tp->withH.ndihe,
           &tp->woH.ndihe, &tp->nhparm, &tp->nparm, &tp->tnexcl, &tp->nres,
           &tp->woHC.nbond, &tp->woHC.nangl, &tp->woHC.ndihe,
           &tp->nBAH.nbond, &tp->nBAH.nangl, &tp->nBAH.ndihe, &tp->natyp,
           &tp->nphb, &tp->ifpert, &tp->pert.nbond, &tp->pert.nangl,
           &tp->pert.ndihe, &tp->wpert.nbond, &tp->wpert.nangl,
           &tp->wpert.ndihe, &tp->ifbox, &tp->nmxrs, &tp->ifcap,
           &tp->blank);

  /*** Now, read each atom name, charge, mass, type, type name, and class ***/
    natm = tp->natom;
    ReadToFlag("ATOM_NAME", inp);
    tp->AtomNames = LoadChar(inp, natm, 4, 20);
    ReadToFlag("CHARGE", inp);
    tp->Charges = LoadDouble(inp, natm);
    ReadToFlag("MASS", inp);
    tp->Masses = LoadDouble(inp, natm);
    ReadToFlag("ATOM_TYPE_INDEX", inp);
    tp->LJIdx = LoadInteger(inp, natm);
    if (adjbnd == 1) {
        AddToIVec(tp->LJIdx, natm, -1);
    }
    ReadToFlag("AMBER_ATOM_TYPE", inp);
    tp->AtomTypes = LoadChar(inp, natm, 4, 20);
    ReadToFlag("TREE_CHAIN_CLASSIFICATION", inp);
    tp->TreeSymbols = LoadChar(inp, natm, 4, 20);

  /*** Read the number of excluded atoms and make lists ***/
    ReadToFlag("NUMBER_EXCLUDED_ATOMS", inp);
    tp->NExcl = LoadInteger(inp, natm);
    tp->ConExcl = (int *) calloc(natm, sizeof(int));
    for (i = 0; i < natm; i++) {
        tp->ConExcl[i] = ISum(tp->NExcl, i);
    }
    ReadToFlag("EXCLUDED_ATOMS_LIST", inp);
    tp->ExclList =
        LoadInteger(inp, tp->ConExcl[natm - 1] + tp->NExcl[natm - 1]);
    if (adjbnd == 1) {
        AddToIVec(tp->ExclList,
                  tp->ConExcl[natm - 1] + tp->NExcl[natm - 1], -1);
    }

  /*** Now get the force constants for all sorts of interactions ***/
    ReadToFlag("NONBONDED_PARM_INDEX", inp);
    tp->NBParmIdx = LoadInteger(inp, tp->ntypes * tp->ntypes);
    if (adjbnd == 1) {
        AddToIVec(tp->NBParmIdx, tp->ntypes * tp->ntypes, -1);
    }
    ReadToFlag("BOND_FORCE_CONSTANT", inp);
    tp->BondK = LoadDouble(inp, tp->nBAH.nbond);
    ReadToFlag("BOND_EQUIL_VALUE", inp);
    tp->BondEq = LoadDouble(inp, tp->nBAH.nbond);
    ReadToFlag("ANGLE_FORCE_CONSTANT", inp);
    tp->AnglK = LoadDouble(inp, tp->nBAH.nangl);
    ReadToFlag("ANGLE_EQUIL_VALUE", inp);
    tp->AnglEq = LoadDouble(inp, tp->nBAH.nangl);
    ReadToFlag("DIHEDRAL_FORCE_CONSTANT", inp);
    tp->DiheK = LoadDouble(inp, tp->nBAH.ndihe);
    ReadToFlag("DIHEDRAL_PERIODICITY", inp);
    tp->DiheN = LoadDouble(inp, tp->nBAH.ndihe);
    ReadToFlag("DIHEDRAL_PHASE", inp);
    tp->DihePhi = LoadDouble(inp, tp->nBAH.ndihe);
    ReadToFlag("LENNARD_JONES_ACOEF", inp);
    tp->LJA = LoadDouble(inp, tp->ntypes * (tp->ntypes + 1) / 2);
    ReadToFlag("LENNARD_JONES_BCOEF", inp);
    tp->LJB = LoadDouble(inp, tp->ntypes * (tp->ntypes + 1) / 2);
    ReadToFlag("HBOND_ACOEF", inp);
    tp->SolA = LoadDouble(inp, tp->nphb * (tp->nphb + 1) / 2);
    ReadToFlag("HBOND_BCOEF", inp);
    tp->SolB = LoadDouble(inp, tp->nphb * (tp->nphb + 1) / 2);
    ReadToFlag("HBCUT", inp);
    tp->HBCut = LoadDouble(inp, tp->nphb * (tp->nphb + 1) / 2);
    ReadToFlag("RADIUS_SET", inp);
    fgets(tp->RadSet, 128, inp);
    ReadToFlag("RADII", inp);
    tp->Radii = LoadDouble(inp, tp->natom);
    ReadToFlag("SCREEN", inp);
    tp->Screen = LoadDouble(inp, tp->natom);

  /*** Get residue limits and labels ***/
    ReadToFlag("RESIDUE_LABEL", inp);
    tp->ResNames = LoadChar(inp, tp->nres, 4, 20);
    ReadToFlag("RESIDUE_POINTER", inp);
    tp->ResLims = LoadInteger(inp, tp->nres);
    if (adjbnd == 1) {
        AddToIVec(tp->ResLims, tp->nres, -1);
    }
    tp->ResLims =
        (int *) realloc(tp->ResLims, (tp->nres + 1) * sizeof(int));
    tp->ResLims[tp->nres] = tp->natom;

  /*** Bond, angle, and dihedral lists ***/
    ReadToFlag("BONDS_INC_HYDROGEN", inp);
    tp->BIncH = (bond *) LoadInteger(inp, 3 * tp->withH.nbond);
    if (adjbnd == 1) {
        AdjustBondArray((int *) tp->BIncH, tp->withH.nbond, 3);
    }
    ReadToFlag("BONDS_WITHOUT_HYDROGEN", inp);
    tp->BNoH = (bond *) LoadInteger(inp, 3 * tp->woH.nbond);
    if (adjbnd == 1) {
        AdjustBondArray((int *) tp->BNoH, tp->woH.nbond, 3);
    }
    ReadToFlag("ANGLES_INC_HYDROGEN", inp);
    tp->AIncH = (angle *) LoadInteger(inp, 4 * tp->withH.nangl);
    if (adjbnd == 1) {
        AdjustBondArray((int *) tp->AIncH, tp->withH.nangl, 4);
    }
    ReadToFlag("ANGLES_WITHOUT_HYDROGEN", inp);
    tp->ANoH = (angle *) LoadInteger(inp, 4 * tp->woH.nangl);
    if (adjbnd == 1) {
        AdjustBondArray((int *) tp->ANoH, tp->woH.nangl, 4);
    }
    ReadToFlag("DIHEDRALS_INC_HYDROGEN", inp);
    tp->HIncH = (dihedral *) LoadInteger(inp, 5 * tp->withH.ndihe);
    if (adjbnd == 1) {
        AdjustBondArray((int *) tp->HIncH, tp->withH.ndihe, 5);
    }
    ReadToFlag("DIHEDRALS_WITHOUT_HYDROGEN", inp);
    tp->HNoH = (dihedral *) LoadInteger(inp, 5 * tp->woH.ndihe);
    if (adjbnd == 1) {
        AdjustBondArray((int *) tp->HNoH, tp->woH.ndihe, 5);
    }

  /*** Miscellaneous ***/
    ReadToFlag("JOIN_ARRAY", inp);
    tp->Join = LoadInteger(inp, natm);
    if (adjbnd == 1) {
        AddToIVec(tp->Join, natm, -1);
    }
    ReadToFlag("IROTAT", inp);
    tp->Rotat = LoadInteger(inp, natm);
    if (adjbnd == 1) {
        AddToIVec(tp->Rotat, natm, -1);
    }

  /*** For periodic boxes ***/
    if (tp->ifbox > 0) {
        ReadToFlag("SOLVENT_POINTERS", inp);
        fscanf(inp, "%d%d%d", &tp->iptres, &tp->nspm, &tp->nspsol);
        ReadToFlag("ATOMS_PER_MOLECULE", inp);
        tp->Nsp = LoadInteger(inp, tp->nspm);
        if (adjbnd == 1) {
            AddToIVec(tp->Nsp, tp->nspm, -1);
        }
        ReadToFlag("BOX_DIMENSIONS", inp);
        fscanf(inp, "%lf%lf%lf%lf", &tp->smbx.beta, &tp->smbx.x,
               &tp->smbx.y, &tp->smbx.z);
    }

  /*** Currently no support for IFCAP > 0, IFPERT > 0, and/or IPOL == 1 ***/
}

/***=======================================================================***/
/*** GetRst: reads a set of coordinates, neglecting whether it is a        ***/
/***         restart file or not.                                          ***/
/***=======================================================================***/
void GetRst(coord * tc, prmtop * tp)
{
    int i;
    char line[MAXLINE];
    FILE *inp;

    if ((inp = fopen(tc->source, "r")) == NULL) {
        printf("GetRst >> Error.  Coordinate file %s not found!\n",
               tc->source);
        exit(1);
    }
    fgets(line, MAXLINE, inp);
    fgets(line, MAXLINE, inp);
    if (tc->rst == 0) {
        sscanf(line, "%d", &tc->natom);
    } else {
        sscanf(line, "%d%lf", &tc->natom, &tc->t);
    }
    if (tc->natom != tp->natom) {
        printf
            ("GetRst >> Error.  %d atoms in topology, %d in coordinate file!\n",
             tp->natom, tc->natom);
        exit(1);
    }
    tc->loc = (double *) malloc(3 * tc->natom * sizeof(double));
    for (i = 0; i < 3 * tc->natom; i++) {
        fscanf(inp, "%lf", &tc->loc[i]);
    }
    if (tc->rst == 1) {
        tc->vel = (double *) malloc(3 * tc->natom * sizeof(double));
        for (i = 0; i < 3 * tc->natom; i++) {
            fscanf(inp, "%lf", &tc->vel[i]);
        }
    }
    if (tp->ifbox > 0) {
        for (i = 0; i < 6; i++) {
            fscanf(inp, "%lf", &tc->box[i]);
            if (i >= 3) {
                tc->box[i] *= PI / 180.0;
            }
        }
    }
}

/***=======================================================================***/
/*** ReadToFlag: reads a topology file from the beginning until a          ***/
/***             particular flag is found, then skips one additional line  ***/
/***             and checks that this line is also a comment line before   ***/
/***             exiting, leaving the file pointer inp set to the next     ***/
/***             line after the desired flag.                              ***/
/***=======================================================================***/
void ReadToFlag(char *fname, FILE * inp)
{
    int collect, flen;
    char line[MAXLINE];

    collect = 0;
    flen = strlen(fname);
    rewind(inp);
    while (collect == 0) {
        if (fgets(line, MAXLINE, inp) == NULL) {
            printf
                ("ReadToFlag >> Error.  End of File reached seeking %s.\n",
                 fname);
            exit(1);
        }
        if (line[0] == '%' && strncmp(&line[6], fname, flen) == 0) {
            fgets(line, MAXLINE, inp);
            if (line[0] != '%') {
                printf
                    ("ReadToFlag >> Error.  Comment line expected following flag.\n"
                     "ReadToFlag >> Found \"%s\"\n", line);
                exit(1);
            }
            collect = 1;
        }
    }
}

/***=======================================================================***/
/*** LoadChar: read in N sets of P characters each, at the rate of Q sets  ***/
/***           per line from file inp.                                     ***/
/***=======================================================================***/
char *LoadChar(FILE * inp, int N, int P, int Q)
{
    int i, j, k;
    char line[MAXLINE];
    char *s;

    i = 0;
    s = (char *) malloc(N * P * sizeof(char));
    while (i < N) {
        if (fgets(line, MAXLINE, inp) == NULL) {
            printf("LoadChar >> Error.  EOF reached.  Terminating.\n");
            exit(1);
        }
        j = 0;
        while (j < Q && i < N) {
            for (k = P * j; k < P * (j + 1); k++) {
                s[(i - j) * P + k] = line[k];
            }
            i++;
            j++;
        }
    }

    return s;
}

/***=======================================================================***/
/*** LoadInteger: read in N integers from the a file inp.                  ***/
/***=======================================================================***/
int *LoadInteger(FILE * inp, int N)
{
    int i;
    int *s;

    s = (int *) malloc(N * sizeof(int));
    for (i = 0; i < N; i++) {
        fscanf(inp, "%d", &s[i]);
    }

    return s;
}

/***=======================================================================***/
/*** LoadDouble: read in N double-precision reals from a file inp.         ***/
/***=======================================================================***/
double *LoadDouble(FILE * inp, int N)
{
    int i;
    double *s;

    s = (double *) malloc(N * sizeof(double));
    for (i = 0; i < N; i++) {
        fscanf(inp, "%lf", &s[i]);
    }

    return s;
}

/***=======================================================================***/
/*** MapResidues: find distinct residues within a topology and map them to ***/
/***              human-readable formats.                                  ***/
/***=======================================================================***/
void MapResidues(prmtop * tp, coord * tc, char *outname)
{
    int i, j, k, nres, rlen, ai, aj, match, ninst, finst, natm;
    int *cpylist;
    double qi, qj, lja, ljb, sig, eps, bfac;
    char *ctmp;
    FILE *outp;

  /*** Some conversion factors ***/
    bfac = 1.0 / sqrt(332.0636);

    cpylist = (int *) malloc(tp->nres * sizeof(int));
    SetIVec(cpylist, tp->nres, -1);
    nres = 0;
    for (i = 0; i < tp->nres; i++) {

    /*** Make sure we haven't done this already ***/
        if (cpylist[i] > -1) {
            continue;
        }

    /*** This residue is number "nres" ***/
        cpylist[i] = nres;

    /*** Now, look for all other instances ***/
        ctmp = &tp->ResNames[4 * i];
        rlen = tp->ResLims[i + 1] - tp->ResLims[i];
        for (j = i + 1; j < tp->nres; j++) {
            if (strncmp(ctmp, &tp->ResNames[4 * j], 4) == 0 &&
                rlen == tp->ResLims[j + 1] - tp->ResLims[j]) {
                match = 1;
                for (k = 0; k < rlen; k++) {
                    qi = tp->Charges[tp->ResLims[i] + k];
                    qj = tp->Charges[tp->ResLims[j] + k];
                    ai = tp->LJIdx[tp->ResLims[i] + k];
                    aj = tp->LJIdx[tp->ResLims[j] + k];
                    if (DNEQ(qi, qj) || ai != aj) {
                        match = 0;
                    }
                }
                if (match == 1) {
                    cpylist[j] = nres;
                }
            }
        }

    /*** Increment our counter ***/
        nres++;
    }

  /*** Now, write the map for each residue ***/
    outp = fopen(outname, "w");
    for (i = 0; i < nres; i++) {

    /*** Find the first instance of the residue and store it as "finst" ***/
        finst = -1;
        ninst = 0;
        for (j = 0; j < tp->nres; j++) {
            if (cpylist[j] == i) {
                if (finst < 0) {
                    finst = j;
                }
                ninst++;
            }
        }

    /*** Print the atomic information about this residue ***/
        natm = tp->ResLims[finst + 1] - tp->ResLims[finst];
        fprintf(outp,
                "%% RESIDUE  %.4s  [ %4d atoms, %8.4lf net charge, %6d "
                "occurences ]\n"
                " Atom Name      Charge   Type Name    LJ Index    LJ Sigma  "
                "LJ Epsilon  Tree Class\n----------  ----------  ----------  "
                "----------  ----------  ----------  ----------\n",
                &tp->ResNames[4 * finst], natm,
                bfac * DSum(&tp->Charges[tp->ResLims[finst]], natm),
                ninst);
        for (j = tp->ResLims[finst]; j < tp->ResLims[finst + 1]; j++) {
            lja = tp->LJA[tp->NBParmIdx[(tp->ntypes + 1) * tp->LJIdx[j]]];
            ljb = tp->LJB[tp->NBParmIdx[(tp->ntypes + 1) * tp->LJIdx[j]]];
            sig = (lja > 1.0e-08
                   && ljb > 1.0e-08) ? pow(lja / ljb, 1.0 / 6.0) : 0.0;
            eps = (sig > 1.0e-08) ? 0.25 * ljb / pow(sig, 6.0) : 0.0;
            fprintf(outp,
                    "      %.4s  %10.6lf        %.4s  %10d  %10.6lf  %10.6lf  "
                    "      %.4s\n", &tp->AtomNames[4 * j],
                    tp->Charges[j] * bfac, &tp->AtomTypes[4 * j],
                    tp->LJIdx[j], sig, eps, &tp->TreeSymbols[4 * j]);
        }
        fprintf(outp, "\n");

    /*** Now seek the bond information on this residue ***/
        fprintf(outp,
                "  Bond A:B   Stiffness  Eq. Length  Mean Value    St. Dev."
                "\n----------  ----------  ----------  ----------  ----------\n");
        PrintBondInfo(tp, (int *) tp->BIncH, tp->withH.nbond, 3, finst,
                      tp->BondK, tp->BondEq, tp->DihePhi, outp, cpylist,
                      tc);
        PrintBondInfo(tp, (int *) tp->BNoH, tp->woH.nbond, 3, finst,
                      tp->BondK, tp->BondEq, tp->DihePhi, outp, cpylist,
                      tc);
        fprintf(outp, "\n");

    /*** Now seek the angle information on this residue ***/
        fprintf(outp,
                "  Angle A:B:C     Stiffness  Eq. Angle  Mean Value    "
                "St. Dev.\n---------------  ----------  ----------  ----------  "
                "----------\n");
        PrintBondInfo(tp, (int *) tp->AIncH, tp->withH.nangl, 4, finst,
                      tp->AnglK, tp->AnglEq, tp->DihePhi, outp, cpylist,
                      tc);
        PrintBondInfo(tp, (int *) tp->ANoH, tp->woH.nangl, 4, finst,
                      tp->AnglK, tp->AnglEq, tp->DihePhi, outp, cpylist,
                      tc);
        fprintf(outp, "\n");

    /*** Now seek the dihedral information on this residue ***/
        fprintf(outp,
                "  Dihedral A:B:C:D       Penalty     Period       Phase  "
                "Mean Value    St. Dev.\n"
                "--------------------  ----------  ----------  ----------  "
                "----------  ----------\n");
        PrintBondInfo(tp, (int *) tp->HIncH, tp->withH.ndihe, 5, finst,
                      tp->DiheK, tp->DiheN, tp->DihePhi, outp, cpylist,
                      tc);
        PrintBondInfo(tp, (int *) tp->HNoH, tp->woH.ndihe, 5, finst,
                      tp->DiheK, tp->DiheN, tp->DihePhi, outp, cpylist,
                      tc);
        fprintf(outp, "\n");
    }
    fclose(outp);
}

/***=======================================================================***/
/*** AdjustBondArray: the bonded arrays are tricky to adjust, because the  ***/
/***                  atom numbers are given as pointer indices and the    ***/
/***                  bond constants are in need of the typical Fortran-C  ***/
/***                  array indexing adjustment.                           ***/
/***                  We assume array A to have N sets of P numbers each,  ***/
/***                  of which the first P-1 are atom numbers.             ***/
/***=======================================================================***/
void AdjustBondArray(int *A, int N, int P)
{
    int h, i, j;

    h = 0;
    for (i = 0; i < N; i++) {
        for (j = 0; j < P - 1; j++) {
            A[h] = abs(A[h]) / 3;
            h++;
        }
        A[h] -= 1;
        h++;
    }
}

/***=======================================================================***/
/*** PrintBondInfo: a general routine for printing bond information about  ***/
/***                a residue.                                             ***/
/***=======================================================================***/
void PrintBondInfo(prmtop * tp, int *B, int nB, int P, int finst,
                   double *K, double *Eq, double *Pd, FILE * outp,
                   int *cpylist, coord * tc)
{
    int i, j, rmin, rmax, okgo, iC, iD;
    double mlen, rmsval;

    rmin = tp->ResLims[finst];
    rmax = tp->ResLims[finst + 1];
    for (i = 0; i < nB; i++) {

    /*** Make sure this bonded term is solely within the residue ***/
        okgo = 1;
        for (j = 0; j < P - 1; j++) {
            if (B[i * P + j] < rmin || B[i * P + j] >= rmax) {
                okgo = 0;
                break;
            }
        }
        if (okgo == 0) {
            continue;
        }

    /*** Print ***/
        fprintf(outp, " ");
        for (j = 0; j < P - 1; j++) {
            fprintf(outp, "%.4s", &tp->AtomNames[4 * B[i * P + j]]);
            if (j < P - 2) {
                fprintf(outp, ":");
            }
        }
        fprintf(outp, "  %10.6lf  %10.6lf", K[B[(i + 1) * P - 1]],
                Eq[B[(i + 1) * P - 1]]);
        if (P == 5) {
            fprintf(outp, "  %10.6lf", Pd[B[(i + 1) * P - 1]]);
        }
        iC = (P > 3) ? B[i * P + 2] : -1;
        iD = (P > 4) ? B[i * P + 3] : -1;
        BondScan(tp, cpylist, finst, B[i * P], B[i * P + 1], iC, iD, tc,
                 &mlen, &rmsval, P - 1);
        fprintf(outp, "  %10.6lf  %10.6lf", mlen, rmsval);
        fprintf(outp, "\n");
    }
}

/***=======================================================================***/
/*** BondScan: computes the mean value of some bond length over all        ***/
/***           matching residues of a simulation.                          ***/
/***=======================================================================***/
void BondScan(prmtop * tp, int *cpylist, int finst, int A, int B, int C,
              int D, coord * tc, double *mval, double *rmsval,
              int bondorder)
{
    int i, j, tcon, aA, aB, aC, aD;
    double ax[3], dx[3];
    double *blist;
    char *cA, *cB, *cC, *cD;

  /*** Determine which atoms "A" and "B" are ***/
    cA = &tp->AtomNames[4 * A];
    cB = &tp->AtomNames[4 * B];
    cC = &tp->AtomNames[4 * C];
    cD = &tp->AtomNames[4 * D];
    if (bondorder > 2) {
        cC = &tp->AtomNames[4 * C];
    }
    if (bondorder > 3) {
        cD = &tp->AtomNames[4 * D];
    }

    tcon = 0;
    blist = (double *) malloc(tp->nres * sizeof(double));
    for (i = 0; i < tp->nres; i++) {
        if (cpylist[i] != cpylist[finst]) {
            continue;
        }

    /*** Find the bond length in this residue ***/
        aA = -1;
        aB = -1;
        aC = -1;
        aD = -1;
        for (j = tp->ResLims[i]; j < tp->ResLims[i + 1]; j++) {
            if (strncmp(&tp->AtomNames[4 * j], cA, 4) == 0) {
                aA = j;
            } else if (strncmp(&tp->AtomNames[4 * j], cB, 4) == 0) {
                aB = j;
            } else if (bondorder > 2
                       && strncmp(&tp->AtomNames[4 * j], cC, 4) == 0) {
                aC = j;
            } else if (bondorder > 3
                       && strncmp(&tp->AtomNames[4 * j], cD, 4) == 0) {
                aD = j;
            }
        }
        if (aA == -1) {
            printf
                ("BondScan >> Error.  Atom %.4s not found in residue %d!\n",
                 cA, i);
            exit(1);
        }
        if (aB == -1) {
            printf
                ("BondScan >> Error.  Atom %.4s not found in residue %d!\n",
                 cB, i);
            exit(1);
        }
        if (bondorder > 2 && aC == -1) {
            printf
                ("BondScan >> Error.  Atom %.4s not found in residue %d!\n",
                 cC, i);
            exit(1);
        }
        if (bondorder > 3 && aD == -1) {
            printf
                ("BondScan >> Error.  Atom %.4s not found in residue %d!\n",
                 cD, i);
            exit(1);
        }
        if (bondorder == 2) {
            for (j = 0; j < 3; j++) {
                dx[j] = tc->loc[3 * aB + j] - tc->loc[3 * aA + j];
            }
            blist[tcon] = DIST(dx[0], dx[1], dx[2]);
        } else if (bondorder == 3) {
            for (j = 0; j < 3; j++) {
                dx[j] = tc->loc[3 * aC + j] - tc->loc[3 * aB + j];
                ax[j] = tc->loc[3 * aA + j] - tc->loc[3 * aB + j];
            }
            blist[tcon] =
                acos(TrimAcos
                     (DotP(ax, dx, 3) / (DMag(dx, 3) * DMag(ax, 3))));
        } else if (bondorder == 4) {
            blist[tcon] =
                Dihedral(&tc->loc[3 * aA], &tc->loc[3 * aB],
                         &tc->loc[3 * aC], &tc->loc[3 * aD]);
        }
        tcon++;
    }
    *mval = DAverage(blist, tcon);
    *rmsval = DStDev(blist, tcon);

  /*** Free Allocated Memory ***/
    free(blist);
}

/***=======================================================================***/
/*** FindDisulfides: find sulfur atoms on CYS residues that are very close ***/
/***                 together and prompt the user if they are not bonded.  ***/
/***=======================================================================***/
void FindDisulfides(prmtop * tp, coord * tc)
{
    int i, j, k, iresno, jresno, bonded, ndss;
    double dx[3];

    iresno = jresno = 0;  /* avoid compiler warnings */
    ndss = 0;
    for (i = 0; i < tp->natom; i++) {
        if ((tp->Masses[i] <= 32.0 || tp->Masses[i] >= 33.0) &&
            tp->AtomNames[4 * i] != 'S') {
            continue;
        }

    /*** This looks like sulfur, but is it part of a cysteine? ***/
        for (j = 0; j < tp->nres; j++) {
            if (tp->ResLims[j] < i) {
                iresno = j;
            }
        }
        if (strncmp(&tp->ResNames[4 * iresno], "CY", 2) != 0 &&
            strncmp(&tp->ResNames[4 * iresno + 1], "CY", 2) != 0) {
            continue;
        }
        for (j = i + 1; j < tp->natom; j++) {
            if ((tp->Masses[j] <= 32.0 || tp->Masses[j] >= 33.0) &&
                tp->AtomNames[4 * j] != 'S') {
                continue;
            }

            for (k = 0; k < tp->nres; k++) {
                if (tp->ResLims[k] < j) {
                    jresno = k;
                }
            }
            if (strncmp(&tp->ResNames[4 * jresno], "CY", 2) != 0 &&
                strncmp(&tp->ResNames[4 * jresno + 1], "CY", 2) != 0) {
                continue;
            }

      /*** So, we've got two atoms that look interesting ***/
            for (k = 0; k < 3; k++) {
                dx[k] = tc->loc[3 * i + k] - tc->loc[3 * j + k];
            }
            if (SQ_DIST(dx[0], dx[1], dx[2]) < 16.0) {

    /*** They're close together... are they bonded? ***/
                bonded = 0;
                for (k = 0; k < tp->woH.nbond; k++) {
                    if ((tp->BNoH[k].a == i && tp->BNoH[k].b == j) ||
                        (tp->BNoH[k].b == i && tp->BNoH[k].a == j)) {
                        bonded = 1;
                        ndss++;
                        break;
                    }
                }
                if (bonded == 0) {
                    printf
                        ("FindDisulfides >> Look at Atoms %.4s %d %.4s and "
                         "%.4s %d %.4s.\n", &tp->ResNames[4 * iresno],
                         iresno + 1, &tp->AtomNames[4 * i],
                         &tp->ResNames[4 * jresno], jresno + 1,
                         &tp->ResNames[4 * jresno]);
                    printf
                        ("FindDisulfides >> There may be a missing disulfide bond.\n");
                }
            }
        }
    }
    if (ndss > 0) {
        printf("FindDisulfides >> There are %d disulfide bonds.\n", ndss);
    }
}

/***=======================================================================***/
/*** VacuumBubble: check for vacuum bubbles in a simulation.               ***/
/***=======================================================================***/
void VacuumBubble(prmtop * tp, coord * tc, double gspc, double bubspc)
{
    int i, j, k, natm, nbub;
    double lja, ljb, sig, eps;
    double ccm[3];
    double *tmpcrd;
    bgrid tg;
    dmat U, invU;

  /*** Transform into crystal lattice space ***/
    U = CreateDmat(3, 3);
    invU = CreateDmat(3, 3);
    CmpXfrm(tc->box, U, invU);
    RotateCrd(tc->loc, tc->natom, U);

  /*** Center the protein and move solvent accordingly ***/
    FindCrdCenter(tc->loc, tp->Masses, 1, tc->natom, ccm);
    TransCrd(tc->loc, tc->natom, ccm, -1.0);

  /*** Re-image the protein ***/
    ReImage(tc->loc, tc->natom);

  /*** Now, color the grid with each atom type ***/
    tg = CreateBoolGrid(tc->box[0] / gspc + 1, tc->box[1] / gspc + 1,
                        tc->box[2] / gspc + 1, -0.5, -0.5, -0.5, 0.5, 0.5,
                        0.5);
    tmpcrd = (double *) malloc(3 * tc->natom * sizeof(double));
    for (i = 0; i < tp->ntypes; i++) {
        fprintf(stderr,
                "\rVacuumBubble >> Coloring grid for atom type %d of %d",
                i, tp->ntypes);
        fflush(stderr);
        natm = 0;
        for (j = 0; j < tp->natom; j++) {
            if (tp->LJIdx[j] == i) {
                for (k = 0; k < 3; k++) {
                    tmpcrd[3 * natm + k] = tc->loc[3 * j + k];
                }
                natm++;
            }
        }
        lja = tp->LJA[tp->NBParmIdx[(tp->ntypes + 1) * i]];
        ljb = tp->LJB[tp->NBParmIdx[(tp->ntypes + 1) * i]];
        sig = (lja > 1.0e-04
               && ljb > 1.0e-04) ? pow(lja / ljb, 1.0 / 6.0) : 0.0;
        eps = (sig > 1.0e-04) ? 0.25 * ljb / pow(sig, 6.0) : 0.0;
        if (eps > 0.0001) {
            sig = 0.5 * (sig + 2.0 * bubspc);
            BiColorGrid(tg, tmpcrd, natm, U, invU, sig, 0);
        }
    }

  /*** Now, scan the grid for any possible vacuum bubbles ***/
    nbub = 0;
    for (i = 0; i < tg.lx; i++) {
        for (j = 0; j < tg.ly; j++) {
            for (k = 0; k < tg.lz; k++) {
                if (ReadBoolGrid(tg, i, j, k) == 0) {
                    tmpcrd[0] = tg.ox + tg.gx * i;
                    tmpcrd[1] = tg.oy + tg.gy * j;
                    tmpcrd[2] = tg.oz + tg.gz * k;
                    BiColorGrid(tg, tmpcrd, 1, U, invU, 2.0 * bubspc, 0);
                    nbub++;
                }
            }
        }
    }

  /*** Tell me how many bubbles I have ***/
    printf
        ("\nVacuumBubble >> There is vacuum space in the box for %d %8.3lfA\n"
         "VacuumBubble >> diameter spheres.  Significant free space implies "
         "the\nVacuumBubble >> existence of unnatural cavities that should "
         "collapse\nVacuumBubble >> under NPT dynamics.\n", nbub, bubspc);

  /*** Free Allocated Memory ***/
    free(tg.data);
    free(tmpcrd);
    DestroyDmat(&U);
    DestroyDmat(&invU);
}

/***=======================================================================***/
/*** FindAtom: finds an atom in a list of names, between limits il and ih. ***/
/***=======================================================================***/
int FindAtom(prmtop * tp, int il, int ih, char *aname)
{
    int i;

    for (i = il; i < ih; i++) {
        if (strncmp(aname, &tp->AtomNames[4 * i], 4) == 0) {
            return i;
        }
    }

  /*** Error message ***/
    printf("FindAtom >> Error.  Could not locate atom %.s\n", aname);
    exit(1);

    return -1;
}

/***=======================================================================***/
/*** Chirality: get the chirality of a center, given the number CA of the  ***/
/***            atom at the center and four atoms around it.  We expect A  ***/
/***            to be the first atom in the list, B the second, C the      ***/
/***            third, and D the fourth atom in order of greatest to least ***/
/***            chiral dominance.                                          ***/
/***=======================================================================***/
int Chirality(prmtop * tp, coord * tc, int CA, int A, int B, int C, int D)
{
    int i;
    double cax[3], ax[3], bx[3], cx[3], dx[3], ra[3], rb[3], rc[3], rd[3];
    double prja[3], prjb[3], prjc[3], acrb[3], bcrc[3];

  /*** Get coordinates ***/
    for (i = 0; i < 3; i++) {
        cax[i] = tc->loc[3 * CA + i];
        ax[i] = tc->loc[3 * A + i];
        bx[i] = tc->loc[3 * B + i];
        cx[i] = tc->loc[3 * C + i];
        dx[i] = tc->loc[3 * D + i];
        ra[i] = ax[i] - cax[i];
        rb[i] = bx[i] - cax[i];
        rc[i] = cx[i] - cax[i];
        rd[i] = dx[i] - cax[i];
    }

  /*** Remove projection ***/
    Project(ra, rd, prja, 3);
    Project(rb, rd, prjb, 3);
    Project(rc, rd, prjc, 3);
    for (i = 0; i < 3; i++) {
        ra[i] -= prja[i];
        rb[i] -= prjb[i];
        rc[i] -= prjc[i];
    }

  /*** Find chirality via cross products ***/
    CrossP(ra, rb, acrb);
    CrossP(rb, rc, bcrc);

  /*** L chirality is when ACRB and BCRC point opposite RD ***/
    if (DotP(acrb, rd, 3) < 0.0 && DotP(bcrc, rd, 3) < 0.0) {
        return 1;
    }

  /*** Otherwise, this has D chirality ***/
    return -1;
}

/***=======================================================================***/
/*** ProteinChiralityCheck: checks the chiral centers of amino acids.  Any ***/
/***                        residues with D chirality are reported.        ***/
/***=======================================================================***/
void ProteinChiralityCheck(prmtop * tp, coord * tc)
{
    int i, nca, nn, nc, ncb, nha;
    char *ctmp;

    for (i = 0; i < tp->nres; i++) {
        ctmp = &tp->ResNames[4 * i];
        if (!
            (strncmp(ctmp, "ALA", 3) == 0 || strncmp(ctmp, "ARG", 3) == 0
             || strncmp(ctmp, "ASN", 3) == 0
             || strncmp(ctmp, "ASP", 3) == 0
             || strncmp(ctmp, "ASH", 3) == 0
             || strncmp(ctmp, "GLH", 3) == 0
             || strncmp(ctmp, "CYS", 3) == 0
             || strncmp(ctmp, "CYX", 3) == 0
             || strncmp(ctmp, "GLN", 3) == 0
             || strncmp(ctmp, "GLU", 3) == 0
             || strncmp(ctmp, "HID", 3) == 0
             || strncmp(ctmp, "HIE", 3) == 0
             || strncmp(ctmp, "HIP", 3) == 0
             || strncmp(ctmp, "HIS", 3) == 0
             || strncmp(ctmp, "ILE", 3) == 0
             || strncmp(ctmp, "LEU", 3) == 0
             || strncmp(ctmp, "LYS", 3) == 0
             || strncmp(ctmp, "MET", 3) == 0
             || strncmp(ctmp, "PHE", 3) == 0
             || strncmp(ctmp, "PRO", 3) == 0
             || strncmp(ctmp, "SER", 3) == 0
             || strncmp(ctmp, "THR", 3) == 0
             || strncmp(ctmp, "TRP", 3) == 0
             || strncmp(ctmp, "TYR", 3) == 0
             || strncmp(ctmp, "VAL", 3) == 0)) {
            continue;
        }

    /*** If we're still here, we've got some amino acid ***/
        nca = FindAtom(tp, tp->ResLims[i], tp->ResLims[i + 1], "CA  ");
        nn = FindAtom(tp, tp->ResLims[i], tp->ResLims[i + 1], "N   ");
        nc = FindAtom(tp, tp->ResLims[i], tp->ResLims[i + 1], "C   ");
        ncb = FindAtom(tp, tp->ResLims[i], tp->ResLims[i + 1], "CB  ");
        nha = FindAtom(tp, tp->ResLims[i], tp->ResLims[i + 1], "HA  ");

    /*** Now, get the chirality of these four atoms ***/
        if (Chirality(tp, tc, nca, nn, nc, ncb, nha) == -1) {
            printf
                ("ChiralityCheck >> Warning.  Residue %.4s %d has D chirality "
                 "at its CA atom.\n", &tp->ResNames[4 * i], i);
        }
    }
}

/***=======================================================================***/
/*** Format5e16: print a chunk of topology file with this numerical        ***/
/***             format.                                                   ***/
/***=======================================================================***/
void Format5e16(char *cname, double *values, int N, FILE * outp)
{
    int h, i, j, slen;
    char nstr[MAXNAME];

    fprintf(outp, "%%FLAG %s\n%%FORMAT(5E16.8)\n", cname);
    h = 0;
    for (i = 0; i < N; i++) {
        sprintf(nstr, "%16.8e", values[i]);
        slen = strlen(nstr);
        for (j = 0; j < slen; j++) {
            if (nstr[j] == 'e') {
                nstr[j] = 'E';
            }
        }
        fprintf(outp, "%s", nstr);
        h++;
        if (h == 5) {
            h = 0;
            fprintf(outp, "\n");
        }
    }
    if (h != 0 || N == 0) {
        fprintf(outp, "\n");
    }
}

/***=======================================================================***/
/*** Format10i8: print a chunk of topology file with this numerical        ***/
/***             format.                                                   ***/
/***=======================================================================***/
void Format10i8(char *cname, int *values, int N, FILE * outp)
{
    int h, i;

    fprintf(outp, "%%FLAG %s\n%%FORMAT(10I8)\n", cname);
    h = 0;
    for (i = 0; i < N; i++) {
        fprintf(outp, "%8d", values[i]);
        h++;
        if (h == 10) {
            h = 0;
            fprintf(outp, "\n");
        }
    }
    if (h != 0 || N == 0) {
        fprintf(outp, "\n");
    }
}

/***=======================================================================***/
/*** Format20a4: print a chunk of topology file with this character        ***/
/***             format.                                                   ***/
/***=======================================================================***/
void Format20a4(char *cname, char *values, int N, FILE * outp)
{
    int h, i;

    fprintf(outp, "%%FLAG %s\n%%FORMAT(20a4)\n", cname);
    h = 0;
    for (i = 0; i < N; i++) {
        fprintf(outp, "%.4s", &values[4 * i]);
        h++;
        if (h == 20) {
            h = 0;
            fprintf(outp, "\n");
        }
    }
    if (h != 0 || N == 0) {
        fprintf(outp, "\n");
    }
}

/***=======================================================================***/
/*** PutPrmTop: put a topology into a file.                                ***/
/***=======================================================================***/
void PutPrmTop(prmtop * tp, char *fname, char *title)
{
    int *itmp;
    double *solty;
    FILE *outp = fopen(fname, "w");

    fprintf(outp, "%s%%FLAG TITLE\n%%FORMAT(20a4)\n%s\n", tp->vstamp,
            title);
    fprintf(outp, "%%FLAG POINTERS\n%%FORMAT(10I8)\n");
    fprintf(outp, "%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n"
            "%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n"
            "%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n%8d\n", tp->natom,
            tp->ntypes, tp->withH.nbond, tp->woH.nbond, tp->withH.nangl,
            tp->woH.nangl, tp->withH.ndihe, tp->woH.ndihe, tp->nhparm,
            tp->nparm, tp->tnexcl, tp->nres, tp->woHC.nbond,
            tp->woHC.nangl, tp->woHC.ndihe, tp->nBAH.nbond, tp->nBAH.nangl,
            tp->nBAH.ndihe, tp->natyp, tp->nphb, tp->ifpert,
            tp->pert.nbond, tp->pert.nangl, tp->pert.ndihe,
            tp->wpert.nbond, tp->wpert.nangl, tp->wpert.ndihe, tp->ifbox,
            tp->nmxrs, tp->ifcap, tp->blank);
    Format20a4("ATOM_NAME", tp->AtomNames, tp->natom, outp);
    Format5e16("CHARGE", tp->Charges, tp->natom, outp);
    Format5e16("MASS", tp->Masses, tp->natom, outp);
    Format10i8("ATOM_TYPE_INDEX", tp->LJIdx, tp->natom, outp);
    Format10i8("NUMBER_EXCLUDED_ATOMS", tp->NExcl, tp->natom, outp);
    Format10i8("NONBONDED_PARM_INDEX", tp->NBParmIdx,
               tp->ntypes * tp->ntypes, outp);
    Format20a4("RESIDUE_LABEL", tp->ResNames, tp->nres, outp);
    Format10i8("RESIDUE_POINTER", tp->ResLims, tp->nres, outp);
    Format5e16("BOND_FORCE_CONSTANT", tp->BondK, tp->nBAH.nbond, outp);
    Format5e16("BOND_EQUIL_VALUE", tp->BondEq, tp->nBAH.nbond, outp);
    Format5e16("ANGLE_FORCE_CONSTANT", tp->AnglK, tp->nBAH.nangl, outp);
    Format5e16("ANGLE_EQUIL_VALUE", tp->AnglEq, tp->nBAH.nangl, outp);
    Format5e16("DIHEDRAL_FORCE_CONSTANT", tp->DiheK, tp->nBAH.ndihe, outp);
    Format5e16("DIHEDRAL_PERIODICITY", tp->DiheN, tp->nBAH.ndihe, outp);
    Format5e16("DIHEDRAL_PHASE", tp->DihePhi, tp->nBAH.ndihe, outp);
    solty = (double *) calloc(tp->natyp, sizeof(double));
    Format5e16("SOLTY", solty, tp->natyp, outp);
    free(solty);
    Format5e16("LENNARD_JONES_ACOEF", tp->LJA,
               tp->ntypes * (tp->ntypes + 1) / 2, outp);
    Format5e16("LENNARD_JONES_BCOEF", tp->LJB,
               tp->ntypes * (tp->ntypes + 1) / 2, outp);
    Format10i8("BONDS_INC_HYDROGEN", (int *) tp->BIncH,
               3 * tp->withH.nbond, outp);
    Format10i8("BONDS_WITHOUT_HYDROGEN", (int *) tp->BNoH,
               3 * tp->woH.nbond, outp);
    Format10i8("ANGLES_INC_HYDROGEN", (int *) tp->AIncH,
               4 * tp->withH.nangl, outp);
    Format10i8("ANGLES_WITHOUT_HYDROGEN", (int *) tp->ANoH,
               4 * tp->woH.nangl, outp);
    Format10i8("DIHEDRALS_INC_HYDROGEN", (int *) tp->HIncH,
               5 * tp->withH.ndihe, outp);
    Format10i8("DIHEDRALS_WITHOUT_HYDROGEN", (int *) tp->HNoH,
               5 * tp->woH.ndihe, outp);
    Format10i8("EXCLUDED_ATOMS_LIST", tp->ExclList,
               tp->ConExcl[tp->natom - 1] + tp->NExcl[tp->natom - 1],
               outp);
    Format5e16("HBOND_ACOEF", tp->SolA, tp->nphb * (tp->nphb + 1) / 2,
               outp);
    Format5e16("HBOND_BCOEF", tp->SolB, tp->nphb * (tp->nphb + 1) / 2,
               outp);
    Format5e16("HBCUT", tp->HBCut, tp->nphb * (tp->nphb + 1) / 2, outp);
    Format20a4("AMBER_ATOM_TYPE", tp->AtomTypes, tp->natom, outp);
    Format20a4("TREE_CHAIN_CLASSIFICATION", tp->TreeSymbols, tp->natom,
               outp);
    Format10i8("JOIN_ARRAY", tp->Join, tp->natom, outp);
    Format10i8("IROTAT", tp->Rotat, tp->natom, outp);
    if (tp->ifbox > 0) {
        itmp = (int *) malloc(3 * sizeof(int));
        itmp[0] = tp->iptres;
        itmp[1] = tp->nspm;
        itmp[2] = tp->nspsol;
        Format10i8("SOLVENT_POINTERS", itmp, 3, outp);
        free(itmp);
        Format10i8("ATOMS_PER_MOLECULE", tp->Nsp, tp->nspm, outp);
        solty = (double *) malloc(4 * sizeof(double));
        solty[0] = tp->smbx.beta;
        solty[1] = tp->smbx.x;
        solty[2] = tp->smbx.y;
        solty[3] = tp->smbx.z;
        Format5e16("BOX_DIMENSIONS", solty, 4, outp);
    }
    fprintf(outp, "%%FLAG RADIUS_SET\n%%FORMAT(1a80)\n" "%s", tp->RadSet);
    Format5e16("RADII", tp->Radii, tp->natom, outp);
    Format5e16("SCREEN", tp->Screen, tp->natom, outp);
    fclose(outp);
}
