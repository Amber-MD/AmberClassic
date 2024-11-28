#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <limits>

#include "saxsDS.h"
#include "const.h"
#include "IO.h"

void saxs_rism_usage () {
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    std::cout << "                           Computing Small Angle X-ray Scattering intensity from 3D-RISM\n";
    std::cout << "                                      Author - Hung Nguyen    tienhung@rutgers.edu\n";
    std::cout << "                                                   Casegroup 2013\n\n";
    std::cout << "Usage: saxs_rism  -g   --grid_dir     folder where all the rism3d output found (expecting *.dx files there)\n";
    std::cout << "                  -s   --solute       pdb file of the solute\n";
    std::cout << "                  -m   --conc_ion     bulk concentration of salt [M]\n";
    std::cout << "                  -w   --conc_wat     water concentration [default 55.34M]\n";
    std::cout << "                  -q   --qcut         momentum transfer q cutoff [default 0.5 A^-1]\n";
    std::cout << "                  -i   --dq           q spacing [default 0.01 A^-1]\n";
    std::cout << "                  -c   --cutoff       distance cutoff [default 20 A]\n";
    std::cout << "                  -a   --anom_f       f' of atomic scattering factor, used for ASAXS calculation,\n";
    std::cout << "                                      currently only applied to Rb+, Sr2+ and Br- [default 0: off-edge]\n";
    std::cout << "                  -x   --exper        the experimental data file to read q from, once specified this overrides dq and qcut\n";
    std::cout << "                                          expect the first column is q (A^-1)\n";
    std::cout << "                  -e   --expli        flag for accounting for explicit H atoms in pdb file\n";
    std::cout << "                  -v   --exclV        flag for merging those contribution of the grid points inside the excluded volume of the solute into the solute\n";
    std::cout << "                  -d   --decomp       flag for decomposing SAXS intensity into site contributions (lead to 2-5x computational time)\n";
    std::cout << "                  -p   --phase        output phase and error analysis instead of partial intensities\n";
    std::cout << "                  -t   --tight        flag for using tighter convergence criteria for Lebedev quadrature (expect more time)\n";
    std::cout << "                  -f   --off_cutoff   flag for turning off cutoff, using all grid points for the calculation\n";
    std::cout << "                  -h   --coion        output a hypothetical contribution from co-ion\n";
    std::cout << "                                      the co-ions are assumed to be completely depleted within a certain range from the solute\n";
    std::cout << "                  -b   --bfactor      using B-factor in the PDB file to account for solute flexibility\n";
    std::cout << "                  -o   --output       output file\n";
#ifdef _OPENMP
    std::cout << "                  -n   --ncpus        number of cpus used [default: 0 - using all available cpus]\n";
#endif
}

/////////////////////////////////////////////////
void saxs_md_usage () {
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    std::cout << "                        A program for computing Small angle X-ray scattering curves from MD simulation\n";
    std::cout << "                                         Author - Hung Nguyen    tienhung@rutgers.edu\n";
    std::cout << "                                                   Casegroup 2013\n\n";
    std::cout << "Usage:  saxs_md  -i   --system       pdb file of the solute\n";
    std::cout << "                 -w   --solvent      pdb file of the solvent\n";
    std::cout << "                 -q   --qcut         momentum transfer q cutoff [default 1.0 A^-1]\n";
    std::cout << "                 -d   --dq           q spacing [default 0.01 A^-1]\n";
    std::cout << "                 -c   --cutoff       generate a box with buffer=cutoff [default 10A]. Only keeping solvent molecules\n";
    std::cout << "                                          within this box for SAXS calculation\n";
    std::cout << "                 -t   --tight        use tighter convergence criteria for Lebedev quadrature\n";
    std::cout << "                 -a   --anom_f       f' for anomalous scattering, used for ASAXS calculation,\n";
    std::cout << "                                     currently only support Rb+, Sr2+ and Br- [default 0: off-edge]\n";
    std::cout << "                 -e   --expli        flag for accounting for explicit H atoms in pdb file\n";
    std::cout << "                 -k   --corrected    correction d for bulk density difference between the blank and sample simulation\n";
    std::cout << "                                          the excess density will be g = (1-d)rho_sample - rho_blank\n";
    std::cout << "                 -b   --bcutoff      minimum distance between the solvent and solute to apply the above correction\n";
    std::cout << "                 -x   --exper        experiment data file to read q from, once specified this overrides dq and qcut\n";
    std::cout << "                                          Expect the first column is q (A^-1)\n";
    std::cout << "                 -o   --output       output file\n";
#ifdef _OPENMP
    std::cout << "                 -n   --ncpus        number of cpus used [default: 0 - using all available cpus]\n";
#endif
}

////////////////////////////////////////////////////
void pddf_rism_usage () {
    std::cout << "Usage: pddf_rism -p   --pdb          pdb file of the solute\n";
    std::cout << "                 -d   --dxfile       dx input file\n";
//  std::  cout << "                 -f   --rdf          electron rdf around water file (if dx is for water)\n";
    std::cout << "                 -b   --bulk_rho     bulk electron density [e/A^3] [default is 0.3333 of pure water]\n";
    std::cout << "                 -w   --conc         bulk concentration of solvent grid [default 55.34M for water]\n";
    std::cout << "                 -m   --merge        how many voxels will be merged into one large voxel\n";
    std::cout << "                                     input format \"x y z\" for three axes, for eg. \"1 2 3\"\n";
    std::cout << "                 -c   --cutoff       distance cutoff in PDDF file output [default 100 A]\n";
    std::cout << "                 -r   --dr           dr in PDDF [default 1 A]\n";
    std::cout << "                 -e   --emap         flag to output excess electron density [default 0]\n";
    std::cout << "                 -t   --type         1 - compute solute  - solvent PDDF [default]\n";
    std::cout << "                                     2 - compute solvent - solvent PDDF\n";
    std::cout << "                                     3 - compute solute  - solute  PDDF\n";
    std::cout << "                 -s   --smear_cut    how far to smear electron out from the nucleus [default 1.0 A]\n";
    std::cout << "                 -o   --output       PDDF output file\n";
#ifdef _OPENMP
    std::cout << "                 -n   --ncpus        number of cpus used [default: 0 - using all available cpus]\n";
#endif
}

//////////////////////////////////////////////////
void pddf_md_usage () {
    std::cout << "Usage:  pddf_md  -p   --pdb          trajectory in pdb format\n";
    std::cout << "                 -c   --cutoff       distance cutoff in PDDF file output [default 100 A]\n";
    std::cout << "                 -r   --dr           dr in PDDF, also grid spacing [default 2 A]\n";
    std::cout << "                 -e   --emap         flag to output excess electron density [default 0]\n";
    std::cout << "                 -b   --bulk_cut     distance from the solute to the bulk [default 20 A]\n";
    std::cout << "                 -t   --type         1 - compute solute  - solvent PDDF [default]\n";
    std::cout << "                                     2 - compute solvent - solvent PDDF\n";
    std::cout << "                                     3 - compute solute  - solute  PDDF\n";
    std::cout << "                 -s   --smear_cut    how far to smear electron out from the nucleus [default 1.0 A]\n";
    std::cout << "                 -o   --output       PDDF output file\n";
#ifdef _OPENMP
    std::cout << "                 -n   --ncpus        number of cpus used [default: 0 - using all available cpus]\n";
#endif
}

///////////////////////////////////////////////////
void saxs_rism_output (const std::vector< std::vector<double> > &intensity,
                       const std::vector<dx_type> &dx,
                       const std::vector<double> &q,
                       std::string &exper,
                       std::string &dx_dir, std::string &pdb_file,
                       double conc_salt, double conc_wat, double cutoff,
                       double anom_f, bool isFlex, bool isOff_cutoff, bool isExpli,
                       bool isExclV, bool isTight, bool isDecomp, bool isPhase,
                       int index_coion,
                       const std::vector< std::vector<DoubPair> > &phase,
                       const std::vector< std::vector<DoubPair> > &dphase,
                       std::ofstream &OUTPUT) {

    OUTPUT << "# saxs_rism ---- Computing Small Angle X-ray Scattering Intensity from 3D-RISM\n";
    OUTPUT << "#                       Author -- Hung Nguyen, tienhung@rutgers.edu\n";
    OUTPUT << "#                                  Casegroup 2013\n\n\n";
    OUTPUT << "#  Program options:\n\n";
    OUTPUT << "#   + Grid                   " << dx_dir << "/*.dx\n";
    OUTPUT << "#   + Pdb                    " << pdb_file << std::endl;
    OUTPUT << "#   + B-factor               ";
    if (isFlex)
        OUTPUT << "ON\n";
    else
        OUTPUT << "OFF\n";
    OUTPUT << "#   + Salt  conc.            " << conc_salt << " mol/l\n";
    OUTPUT << "#   + Water conc.            " << conc_wat << " mol/l\n";
    OUTPUT << "#   + Use all grid points    ";
    if (isOff_cutoff)
        OUTPUT << "ON\n";
    else {
        OUTPUT << "OFF\n";
        OUTPUT << "#   + Space cutoff           " << cutoff << " Angstrom\n";
    }
    OUTPUT << "#   + Anomalous f'           " << anom_f << std::endl;
    OUTPUT << "#   + Explicit hydrogen      ";
    if (isExpli)
        OUTPUT << "ON\n";
    else
        OUTPUT << "OFF\n";
    OUTPUT << "#   + Tight convergence      ";
    if (isTight)
        OUTPUT << "ON\n";
    else
        OUTPUT << "OFF\n";
    if (not exper.empty())
        OUTPUT << "#   + Exp file               " << exper << std::endl;
    if (not isPhase)
        saxs_rism_output_intensity (intensity, dx, q, isExclV, isDecomp, index_coion, OUTPUT);
    else
        saxs_rism_output_phase (intensity, dx, q, isExclV, isDecomp, index_coion, phase, dphase, OUTPUT);
}

///////////////////////////////////////////////////////////////////
void saxs_rism_output_intensity (const std::vector< std::vector<double> > &intensity,
                                 const std::vector<dx_type> &dx,
                                 const std::vector<double> &q,
                                 bool isExclV, bool isDecomp, int index_coion,
                                 std::ofstream &OUTPUT) {
    OUTPUT << "#   + Output                INTENSITY\n";
    OUTPUT << "\n##  q    Solute(vac)  ";
    if (isExclV)
        OUTPUT << "Solute(solv)   ";
    if (isDecomp)
        for (size_t i = 0; i < dx.size(); i++)
            OUTPUT << "     " << dx[i].type << "      ";
    else
        OUTPUT << "     Hyd      ";
    if (index_coion > -1)
        OUTPUT << "  Excl" << dx[index_coion].type << "   ";
    OUTPUT << "  Total\n";

    for (size_t i = 0; i < q.size(); i++) {
        OUTPUT << std::setw(6) << std::setprecision(3) << std::fixed << q[i];
        for (size_t j = 0; j < intensity[i].size(); j++)
            OUTPUT << std::setw(14) << std::setprecision(6) << std::fixed << std::scientific << intensity[i][j];
        OUTPUT << std::endl;
    }
}

////////////////////////////////////////////////////////////////////////
void saxs_rism_output_phase (const std::vector< std::vector<double> > &intensity,
                             const std::vector<dx_type> &dx,
                             const std::vector<double> &q,
                             bool isExclV, bool isDecomp, int index_coion,
                             const std::vector< std::vector<DoubPair> > &phase,
                             const std::vector< std::vector<DoubPair> > &dphase,
                             std::ofstream &OUTPUT) {
    OUTPUT << "#   + Output                 PHASE\n";
    OUTPUT << "\n##  q  Solu(vac) dSolu(vac) ";
    if (isExclV)
        OUTPUT << "Solu(solv) dSolu(solv) ";
    if (isDecomp)
        for (size_t i = 0; i < dx.size(); i++)
            OUTPUT << "  " << dx[i].type << "      d" << dx[i].type << "    ";
    else
        OUTPUT << " Hyd     dHyd     ";
    if (index_coion > -1)
        OUTPUT << " Excl" << dx[index_coion].type << "   dExcl" << dx[index_coion].type << "  ";
    OUTPUT << " Total   dTotal   ";
    if (isDecomp)
        for (size_t i = 0; i < dx.size(); i++)
            OUTPUT << "Solu-" << dx[i].type << " dSolu-" << dx[i].type << "  ";
    else OUTPUT << "Solu-Hyd   dSolu-Hyd  ";
    OUTPUT << "Solu-Solv  dSolu-Solv  Intensity\n";

    for (size_t i = 0; i < q.size(); i++) {
        OUTPUT << std::setw(6) << std::setprecision(3) << std::fixed << q[i];
        for (size_t j = 0; j < phase[i].size(); j++) {
            OUTPUT << std::setw(9) << std::setprecision(1) << std::fixed << phase[i][j].first;
            OUTPUT << std::setw(9) << std::setprecision(1) << std::fixed << phase[i][j].second;
        }
        for (size_t j = 0; j < dphase[i].size(); j++) {
            OUTPUT << std::setw(9) << std::setprecision(1) << std::fixed << dphase[i][j].first;
            OUTPUT << std::setw(9) << std::setprecision(1) << std::fixed << dphase[i][j].second;
        }
        OUTPUT << std::setw(18) << std::setprecision(6) << std::fixed << std::scientific << intensity[i].back() << std::endl;
    }
}

///////////////////////////////////////////////////////////////////////
void saxs_md_output (const std::vector<double> &intensity,
                     double dq,
                     std::string pdb_solu, std::string pdb_solv,
                     double cutoff, double bulk_corr, double bulk_cutoff,
                     double anom_f, bool isExpli, bool isTight,
                     std::ofstream &OUTPUT) {

    OUTPUT << " saxs_md ---- Computing Small angle X-ray scattering curves from MD simulation\n";
    OUTPUT << "                        Author - Hung Nguyen    tienhung@rutgers.edu\n";
    OUTPUT << "                                     Casegroup 2013\n\n";
    OUTPUT << "# Program options:\n";
    OUTPUT << "#      + Solute          " << pdb_solu << std::endl;
    OUTPUT << "#      + Solvent         " << pdb_solv << std::endl;
    OUTPUT << "#      + Cutoff          " << cutoff << " Angstrom\n";
    if (std::abs(bulk_corr) > std::numeric_limits<double>::epsilon()) {
        OUTPUT << "#      + Bulk correction " << bulk_corr << std::endl;
        OUTPUT << "#      + Bulk cutoff     " << bulk_cutoff << " Angstrom\n";
    }
    OUTPUT << "#      + Anomalous f'    " << anom_f << std::endl;
    OUTPUT << "#      + Explicit H      ";
    if (isExpli)
        OUTPUT << "ON\n";
    else OUTPUT << "OFF\n";
    OUTPUT << "#      + Tight conv      ";
    if (isTight)
        OUTPUT << "ON\n";
    else OUTPUT << "OFF\n";
    OUTPUT << "\n\n##  q       Intensity\n";

    double q = 0;
    for (size_t i = 0; i < intensity.size(); i++) {
        OUTPUT << std::setw(8) << std::setprecision(3) << std::fixed << q;
        OUTPUT << std::setw(14) << std::setprecision(6) << std::fixed << std::scientific << intensity[i] << std::endl;
        q += dq;
    }
}

/////////////////////////////////////////////////////////////
void pddf_output (const std::vector<double> &pddf,
                  size_t type,
                  double dr,
                  std::ofstream &OUTPUT) {
    if (type == 1)
        OUTPUT << "###  Solute - solvent PDDF\n";
    else if (type == 2)
        OUTPUT << "###  Solvent - solvent PDDF\n";
    else if (type == 3)
        OUTPUT << "###  Solute - solute PDDF\n";
    OUTPUT << "# r(A)       p(r)\n  0.00  0.000000e+00\n";
    for (size_t i = 0; i < pddf.size(); i++)
        OUTPUT << std::setw(6)  << std::setprecision(2) << std::fixed << (0.5 + i)*dr\
            << std::setw(14) << std::setprecision(6) << std::fixed << std::scientific << pddf[i] << std::endl;
}
////////////////////////////////////////////////////////////
// Read RDF of electron around water (for smearing purpose)
///////////////////////////////////////////////////////////
void read_RDF (const std::string &filename,
               double &dr,
               std::vector<double> &rdf) {

    std::ifstream OPENFILE (filename.c_str());
    if (OPENFILE.is_open()) {
        std::string line;
        double diff_bk = 0;
        double tmp1 = 0; 
        while (getline(OPENFILE, line)) {
            std::istringstream iss(line);
            // Skip line starting with #
            size_t found = line.find_first_not_of(" \t");
            if (found != std::string::npos)
                if (line[found] == '#')
                    continue;
            double tmp2, value;
            while (iss >> tmp2 >> value)
                rdf.push_back(value);
            double diff = tmp2 - tmp1;
            if ((std::abs(diff_bk) >= 1e-10) and (diff - diff_bk > 1e-10)) {
                std::cout << "dr is not a constant!! Quit!!!!\n";
                std::cout << diff_bk << std::endl << diff << std::endl;
                exit(0);
            }
            diff_bk = diff;
            tmp1 = tmp2;
        }
        OPENFILE.close();
        dr = diff_bk;
    }
}

///////////////////////////////////////////////////////////////////////
// Read experimental data, expecting the first column to be q(A^-1)
// # q(A^-1)
//////////////////////////////////////////////////////////////////////
void read_exp (const std::string &exper,
               std::vector<double> &q) {

    std::ifstream EXP (exper.c_str());
    if (EXP.is_open()) {
        q.push_back(0);
        std::string line;
        while (getline(EXP, line)) {
            std::istringstream iss(line);
            double value;
            if (iss >> value) {
                q.push_back(value);
                std::string tmp;
                while (iss >> tmp) {        // Do nothing
                }
            }
        }
        EXP.close();
    } else {
        std::cerr << "Unable to open file " << exper << std::endl;
        exit (0);
    }
}
