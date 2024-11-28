#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib>

#include "const.h"
#include "saxsDS.h"

///////////////////////////////////////////////////////////////////
std::vector<coeff_f> assign_F_table () {
    std::vector<coeff_f> result;
    coeff_f H, C, O, N, P, S, Fe, CL, BR, NA, K, RB, CS, MG, SR;

    // Coefficients for neutral atom taken from Su and Coppens, Acta. Cryst. 1997, A53, 749-762
    //              for ions         taken from Macchi and Coppens, Acta. Cryst. 2001, A57, 656-662
    // H
    H.atomtype = "H";
    H.a1 = 0.43028;         H.b1 = 23.02312;
    H.a2 = 0.28537;         H.b2 = 10.20138;
    H.a3 = 0.17134;         H.b3 = 51.25444;
    H.a4 = 0.09451;         H.b4 = 4.13511;
    H.a5 = 0.01725;         H.b5 = 1.35427;
    H.a6 = 0.00114;         H.b6 = 0.24269;
    result.push_back (H);
    // C
    C.atomtype = "C";
    C.a1 = 2.09921;         C.b1 = 13.18997;
    C.a2 = 1.80832;         C.b2 = 30.37956;
    C.a3 = 1.26159;         C.b3 = 0.69255;
    C.a4 = 0.56775;         C.b4 = 0.16381;
    C.a5 = 0.26303;         C.b5 = 68.42774;
    C.a6 = 0.00010;         C.b6 = 0.44083;
    result.push_back (C);
    // O
    O.atomtype = "O";
    O.a1 = 2.34752;         O.b1 = 9.69710;
    O.a2 = 1.83006;         O.b2 = 18.59876;
    O.a3 = 1.61538;         O.b3 = 5.19879;
    O.a4 = 1.52402;         O.b4 = 0.32408;
    O.a5 = 0.41423;         O.b5 = 39.79099;
    O.a6 = 0.26867;         O.b6 = 0.01150;
    result.push_back (O);
    // N
    N.atomtype = "N";
    N.a1 = 2.45424;         N.b1 = 18.66694;
    N.a2 = 2.15782;         N.b2 = 8.31271;
    N.a3 = 1.05782;         N.b3 = 0.46989;
    N.a4 = 0.57557;         N.b4 = 42.44646;
    N.a5 = 0.44959;         N.b5 = 0.08747;
    N.a6 = 0.30480;         N.b6 = 0.47126;
    result.push_back (N);
    // P, correction from Su and Coppens, Acta. Cryst. 1998, A54, 357
    P.atomtype = "P";
    P.a1 = 6.48197;         P.b1 = 1.89537;
    P.a2 = 4.31666;         P.b2 = 27.61455;
    P.a3 = 1.73759;         P.b3 = 0.50991;
    P.a4 = 1.35793;         P.b4 = 66.28296;
    P.a5 = 1.10559;         P.b5 = 0.00010;
    P.a6 = 0.00010;         P.b6 = 12.05652;
    result.push_back (P);
    // S
    S.atomtype = "S";
    S.a1 = 6.90565;         S.b1 = 1.46764;
    S.a2 = 5.24410;         S.b2 = 22.31576;
    S.a3 = 1.54516;         S.b3 = 56.06328;
    S.a4 = 1.42922;         S.b4 = 0.25588;
    S.a5 = 0.87564;         S.b5 = 0.00010;
    S.a6 = 0.00010;         S.b6 = 26.96892;
    result.push_back (S);
    // Fe
    Fe.atomtype = "Fe";
    Fe.a1 = 11.18858;       Fe.b1 = 4.64599;
    Fe.a2 = 7.37206;        Fe.b2 = 0.30327;
    Fe.a3 = 3.55141;        Fe.b3 = 12.07655;
    Fe.a4 = 1.68125;        Fe.b4 = 44.15316;
    Fe.a5 = 1.20893;        Fe.b5 = 104.11866;
    Fe.a6 = 0.99652;        Fe.b6 = 0.00010;
    result.push_back (Fe);
    // Cl-
    CL.atomtype = "Cl-";
    CL.a1 = 7.13932;        CL.b1 = 1.18073;
    CL.a2 = 6.34213;        CL.b2 = 19.52901;
    CL.a3 = 2.29801;        CL.b3 = 61.04850;
    CL.a4 = 1.97826;        CL.b4 = 0.08057;
    CL.a5 = 0.22854;        CL.b5 = 23.18225;
    CL.a6 = 0.00983;        CL.b6 = 0.09759;
    result.push_back (CL);
    // Br-
    BR.atomtype = "Br-";
    BR.a1 = 14.72809;       BR.b1 = 1.87781;
    BR.a2 = 7.73340;        BR.b2 = 0.11285;
    BR.a3 = 4.08153;        BR.b3 = 23.45650;
    BR.a4 = 3.89920;        BR.b4 = 3.65207;
    BR.a5 = 2.84995;        BR.b5 = 21.50646;
    BR.a6 = 2.70412;        BR.b6 = 68.50430;
    result.push_back (BR);
    // Na+
    NA.atomtype = "Na+";
    NA.a1 = 3.69529;        NA.b1 = 3.24183;
    NA.a2 = 3.30459;        NA.b2 = 7.07179;
    NA.a3 = 1.68333;        NA.b3 = 0.12279;
    NA.a4 = 0.69149;        NA.b4 = 15.45334;
    NA.a5 = 0.62431;        NA.b5 = 1.43664;
    NA.a6 = 0.00088;        NA.b6 = 35.26383;
    result.push_back (NA);
    // K+
    K.atomtype = "K+";
    K.a1 = 8.00372;         K.b1 = 12.70476;
    K.a2 = 7.44077;         K.b2 = 0.77473;
    K.a3 = 1.42217;         K.b3 = 0.00010;
    K.a4 = 1.13491;         K.b4 = 32.44270;
    K.a5 = 0.00010;         K.b5 = 199.99900;
    K.a6 = 0.00010;         K.b6 = 82.98298;
    result.push_back (K);
    // Rb+
    RB.atomtype = "Rb+";
    RB.a1 = 17.72736;       RB.b1 = 1.68258;
    RB.a2 = 7.70846;        RB.b2 = 0.09962;
    RB.a3 = 6.22707;        RB.b3 = 13.34713;
    RB.a4 = 4.23320;        RB.b4 = 25.64859;
    RB.a5 = 0.10456;        RB.b5 = 76.90928;
    RB.a6 = 0.00010;        RB.b6 = 199.99860;
    result.push_back (RB);
    // Mg2+
    MG.atomtype = "Mg2+";
    MG.a1 = 4.30385;        MG.b1 = 4.02045;
    MG.a2 = 2.58390;        MG.b2 = 1.85304;
    MG.a3 = 1.71397;        MG.b3 = 0.10693;
    MG.a4 = 1.39368;        MG.b4 = 8.78523;
    MG.a5 = 0.00470;        MG.b5 = 58.58712;
    MG.a6 = 0.00010;        MG.b6 = 125.50050;
    result.push_back (MG);
    // Sr2+
    SR.atomtype = "Sr2+";
    SR.a1 = 13.56253;       SR.b1 = 1.52639;
    SR.a2 = 9.15282;        SR.b2 = 13.37893;
    SR.a3 = 7.57461;        SR.b3 = 0.09009;
    SR.a4 = 4.23621;        SR.b4 = 1.50827;
    SR.a5 = 1.47524;        SR.b5 = 28.97999;
    SR.a6 = 0.00010;        SR.b6 = 162.86130;
    result.push_back (SR);

    return result;
}

//////////////////////////////////////////////
double atom_fact (const coeff_f &atom,
                  double q) {
    double q2 = q*q/(16*PI*PI);
    return atom.a1*exp(-atom.b1*q2) + atom.a2*exp(-atom.b2*q2) + atom.a3*exp(-atom.b3*q2) + atom.a4*exp(-atom.b4*q2) +\
         + atom.a5*exp(-atom.b5*q2) + atom.a6*exp(-atom.b6*q2);
}

//////////////////////////////////////
double count_e (const std::string &type,
                double anom_f) {
    if ((type == "Ow") or (type == "O"))
        return 8;
    else if ((type == "Hw") or (type == "H"))
        return 1;
    else if (type == "C")
        return 6;
    else if (type == "N")
        return 7;
    else if (type == "P")
        return 15;
    else if (type == "S")
        return 16;
    else if ((type == "Na+") or (type == "Mg2+") or (type == "F-"))
        return 10;
    else if ((type == "K+") or (type == "Cl-") or (type == "Ca2+"))
        return 18;
    else if ((type == "Rb+") or (type == "Br-") or (type == "Sr2+"))
        return 36 + anom_f;
    else if ((type == "Cs+") or (type == "Ba2+") or (type == "I-"))
        return 54;
    else {
        std::cout << "ERROR!!! Not able to recognize atom   " << type << std::endl;
        exit (0);
    }
}

////////////////////////////////////////////////
double atom_rho (const std::string &type,
                 double dist_sqr,
                 const std::vector<coeff_f> &F_table) {
    int index = -1;
    std::string type2 = type;
    if (type2 == "Ow") type2 = "O";
    else if (type2 == "Hw") type2 = "H";
    for (size_t i = 0; i < F_table.size(); i++)
        if (type2 == F_table[i].atomtype) {
            index = i;
            break;
        }
    if (index == -1) {
        std::cerr << "Unable to recognize atom " << type << std::endl;
        exit (0);
    }
    coeff_f atom = F_table[index];
    return 8*PI*sqrt(PI) * (atom.a1 / (atom.b1*sqrt(atom.b1)) * exp(-4*PI*PI*dist_sqr / atom.b1) +\
                            atom.a2 / (atom.b2*sqrt(atom.b2)) * exp(-4*PI*PI*dist_sqr / atom.b2) +\
                            atom.a3 / (atom.b3*sqrt(atom.b3)) * exp(-4*PI*PI*dist_sqr / atom.b3) +\
                            atom.a4 / (atom.b4*sqrt(atom.b4)) * exp(-4*PI*PI*dist_sqr / atom.b4) +\
                            atom.a5 / (atom.b5*sqrt(atom.b5)) * exp(-4*PI*PI*dist_sqr / atom.b5) +\
                            atom.a6 / (atom.b6*sqrt(atom.b6)) * exp(-4*PI*PI*dist_sqr / atom.b6));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Assign the atomic factor for atoms, summing up the atomic factors of hydrogen atoms into heavy atom
// Using the analytical approximation for scattering factors,
//////////////////////////////////////////////////////////////////////////////////////////////////////
double f_atm (const std::string &type,
              double q,
              size_t nHyd,
              bool expli,
              bool grid,
              const std::vector<coeff_f> &F_table) {

    double atomic_factor;
    if (grid) {
        if ((type == "ex_e") or ((type == "Hw") and (expli)))
            atomic_factor = 1;
        else if (type == "Ow") {
            if (expli) atomic_factor = 8;
            else atomic_factor = 10;
        } else if ((type == "F-") or (type == "Na+") or (type == "Mg2+"))
            atomic_factor = 10;
        else if ((type == "Cl-") or (type == "K+") or (type == "Ca2+"))
            atomic_factor = 18;
        else if ((type == "Br-") or (type == "Rb+") or (type == "Sr2+"))
            atomic_factor = 36;
        else if ((type == "I-") or (type == "Cs+") or (type == "Ba2+"))
            atomic_factor = 54;
        else {
            std::cerr << "Unable to recognize grid type " << type << std::endl;
            exit (0);
        }
    } else {
        int index = -1;
        std::string type_atom = type;
        if (type == "Ow")
            type_atom = "O";
        else if (type == "Hw")
            type_atom = "H";
        for (size_t i = 0; i < F_table.size(); i++)
            if (type_atom == F_table[i].atomtype) {
                index = i;
                break;
            }
        if (index == -1) {
            std::cerr << "Unable to recognize atom " << type_atom << std::endl;
            exit (0);
        }
        if (expli)
            atomic_factor = atom_fact (F_table[index], q);
        else
            atomic_factor = atom_fact (F_table[index], q) + nHyd*atom_fact (F_table[0], q);
    }
    return atomic_factor;
}
