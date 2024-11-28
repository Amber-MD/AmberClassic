/* A program to generate SAXS and ASAXS profiles from 3D-RISM grids 
   More details are in Hung Nguyen et al, JCP 141, 2014; DOI: 10.1063/1.4896220

   Currently only support ASAXS calculation for Rb+ Sr2+ and Br-, and single salt solution (will add mixed soon?)
 
   Written by Hung Nguyen, Case's group, 2013 */

#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include <vector>
#include <cstdlib>
#include <getopt.h>
#include <dirent.h>

#include "saxsDS.h"
#include "lebedev.h"
#include "dx.h"
#include "vol_decomp.h"
#include "saxs_rism.h"
#include "atm_f.h"
#include "f_calc.h"
#include "const.h"
#include "pdb.h"
#include "IO.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

struct comparatorStruct {
    bool operator() (const list_cutoff &a,
                     const list_cutoff &b) {
        return a.index != b.index;
    }
};

//////////////////////////////////////////////////////////////////////////
int main (int argc,
          char *argv[]) {

    /// Declare variables and their default values
    string dx_dir, pdb_file, outfile;
    double qcut = 0.5;      // Cutoff q
    double anom_f = 0;      // f' in anomalous scattering
    double conc_salt = 0;   // Bulk concentration of salt
    double conc_wat = 55.34;// Water concentration
    double cutoff = 20;
    bool isExpli = 0;       // Account for explicit H atoms in pdb file
    bool isTight = 0;       // Control Lebedev quadrature tight or loose convergence
    bool isFlex = 0;        // Account for flexibility using B-factor in the PDB file
    //bool corr = 0;        // Using corrected atomic factor for water as in J Chem Phys 2000, 113, 9149
    bool isDecomp = 0;      // Decompose intensity into site contributions
    bool isCoion = 0;       // Output hypothetical co-ion exclusion contribution
    bool isExclV = 0;       // Control the decomp output
    bool isPhase = 0;       // Control the phase and error output
    string exper;           // Read q from this file
    bool isOff_cutoff = 0;
    double dq = 0.01;
    unsigned ncpus = 0;

    int option_char;
    do {
        struct option long_options[] =
        {
            {"anom_f",      required_argument,  0,  'a'},
            {"grid_dir",    required_argument,  0,  'g'},
            {"solute",      required_argument,  0,  's'},
            {"bfactor",     no_argument,        0,  'b'},
            {"conc_salt",   required_argument,  0,  'm'},
            {"conc_wat",    required_argument,  0,  'w'},
            {"cutoff",      required_argument,  0,  'c'},
            {"qcut",        required_argument,  0,  'q'},
            {"dq",          required_argument,  0,  'i'},
            {"exper",       required_argument,  0,  'x'},
            {"expli",       no_argument,        0,  'e'},
            {"exclV",       no_argument,        0,  'v'},
            {"off_cutoff",  no_argument,        0,  'f'},
            {"coion",       no_argument,        0,  'h'},
            {"decomp",      no_argument,        0,  'd'},
            {"phase",       no_argument,        0,  'p'},
            {"tight",       no_argument,        0,  't'},
            {"output",      required_argument,  0,  'o'},
            {"ncpus",       required_argument,  0,  'n'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        option_char = getopt_long (argc, argv, "a:g:s:bm:w:c:q:i:x:evfhdpto:n:", long_options, &option_index);
        // followed by 1 colon - require an argument; no colon - not argument required

        switch (option_char) {
            case '?':
                saxs_rism_usage();
                exit (0);
            case 'a':
                anom_f = atof (optarg);
                break;
            case 'g':
                dx_dir = optarg;
                break;
            case 's':
                pdb_file = optarg;
                break;
            case 'm':
                conc_salt = atof (optarg);
                break;
            case 'w':
                conc_wat = atof (optarg);
                break;
            case 'c':
                cutoff = atof (optarg);
                break;
            case 'q':
                qcut = atof (optarg);
                break;
            case 'i':
                dq = atof (optarg);
                break;
            case 'x':
                exper = optarg;
                break;
            case 'e':
                isExpli = 1;
                break;
            case 'v':
                isExclV = 1;
                break;
            case 'f':
                isOff_cutoff = 1;
                break;
            case 'h':
                isCoion = 1;
                break;
            case 't':
                isTight = 1;
                break;
            case 'd':
                isDecomp = 1;
                break;
            case 'p':
                isPhase = 1;
                break;
            case 'b':
                isFlex = 1;
                break;
            case 'o':
                outfile = optarg;
                break;
            case 'n':
                ncpus = atoi (optarg);
                break;
        }
    } while (option_char != -1);

    if (not isDecomp)
        isCoion = 0;

    ofstream OUTPUT (outfile.c_str());
    if (OUTPUT.is_open()) {
        // cout << "Reading input files ...\n";
        vector<dx_type> dx;
        vector<coordinate> pdb_coord;
        read_dx_dir (dx_dir, dx, conc_wat, conc_salt);
        read_pdb (pdb_file, pdb_coord);
        vector<double> q;

        if (exper.empty())
            for (size_t i = 0; i <= floor (qcut/dq); i++)
                q.push_back(i*dq);
        else
            read_exp (exper, q);

        vector< vector<DoubPair> > phase, dphase;
        int index_coion = -1;
        vector< vector<double> > intensity = saxs_rism_calc (dx, pdb_coord, q, cutoff, isTight, isFlex, anom_f,
                                                             isExpli, isDecomp, isOff_cutoff, isExclV, isCoion, index_coion,
                                                             ncpus, phase, dphase);

        saxs_rism_output (intensity, dx, q, exper, dx_dir, pdb_file, conc_salt, conc_wat, cutoff,
                          anom_f, isFlex, isOff_cutoff, isExpli, isExclV, isTight, isDecomp, isPhase,
                          index_coion, phase, dphase, OUTPUT);
        OUTPUT.close();
    } else {
        cout << "Unable to write to file " << outfile << endl;
        exit (0);
    }
    return 0;
}

//////////////////// END MAIN ///////////////////////////

////////////////////////////////////////
// SAXS RISM calculator main subroutine
////////////////////////////////////////
vector< vector<double> > saxs_rism_calc (const vector<dx_type> &dx,
                                         vector<coordinate> &pdb_coord,
                                         const vector<double> &q,
                                         double cutoff,
                                         bool isTight, bool isFlex, double anom_f,
                                         bool isExpli, bool isDecomp, bool isOff_cutoff, bool isExclV,
                                         bool isCoion, int &index_coion, unsigned ncpus,
                                         vector< vector<DoubPair> > &phase,
                                         vector< vector<DoubPair> > &dphase) {
#ifdef _OPENMP
    omp_set_dynamic(0);
    if (ncpus == 0)
        ncpus = omp_get_max_threads();
    omp_set_num_threads(ncpus);
#endif

    if (not isExpli)
        mergeH (pdb_coord);
    // cout << "Making lists ...\n";
    vector< vector<list_cutoff> > list_exclV, list_hyd;
    vector<list_cutoff> list_excl_coion;
    generate_list (dx, pdb_coord, cutoff, isOff_cutoff, isDecomp, isExclV,
                   list_exclV, list_hyd, list_excl_coion, isCoion, index_coion);

    vector<coeff_f> F_table = assign_F_table ();
    vector<dx_type> excess_e;
    if (not isDecomp) {
        dx_type ex_e = ex_elec (dx, anom_f, isExpli, F_table);
        excess_e.push_back (ex_e);
    }

    cout << "Calculating SAXS ...\n";
    vector< vector<double> > result;
    phase.resize  (q.size());
    dphase.resize (q.size());

    for (size_t i = 0; i < q.size(); i++) {
        size_t rule;
        if (not isTight)
            rule = 4 + floor (q[i]/.04);
        else
            rule = 5 + floor (q[i]/.03);
        while ((available_table(rule) == 0) and (rule < 65))
            // Maximum rule is 65, rarely use up to this number though
            rule++;

        vector<double> intensity;
        if (not isDecomp)
            intensity = saxs_rism_I (excess_e, 1, pdb_coord, list_hyd, 
                 list_exclV, list_excl_coion,
                 index_coion, q[i], rule, isFlex, anom_f, isExpli, isExclV, 
                 F_table, dphase[i], phase[i]);
        else
            intensity = saxs_rism_I (dx,       0, pdb_coord, list_hyd, 
                 list_exclV, list_excl_coion,
                 index_coion, q[i], rule, isFlex, anom_f, isExpli, isExclV, 
                 F_table, dphase[i], phase[i]);
        if( i%10 == 1 ) cout << "   done with " << q[i] << endl;
        result.push_back (intensity);
    }
    return result;
}

/////////////////////////
// Read all the dx files
/////////////////////////
void read_dx_dir (const string &dx_dir,
                  vector<dx_type> &dx,
                  double conc_wat, double conc_salt) {
    DIR *dir = NULL;
    struct dirent *file = NULL;
    if ((dir = opendir (dx_dir.c_str())) == NULL) {
        cout << "Not able to open " << dx_dir << endl;
        exit (0);
    }
    size_t valence = 0;  // monovalent ion -> 1; divalent ion -> 2
    size_t flen;
    while ((file = readdir (dir)) != NULL) {
        string filename = file->d_name;
        flen = filename.length();
        
        if (flen > 3 && filename.substr(flen-3,3) == ".dx") {
            string path = dx_dir + "/" + filename;
            dx_type newdx;

            size_t pos = filename.find ("O");
            if (pos != std::string::npos) {
                newdx.conc = conc_wat;
                newdx.type = "Ow";
            }
            pos = filename.find ("H1");
            if (pos != std::string::npos) {
                newdx.conc = 2*conc_wat;
                newdx.type = "Hw";
            }

            pos = filename.find ("F-");
            if (pos != std::string::npos) newdx.type = "F-";

            pos = filename.find ("Cl-");
            if (pos != std::string::npos) newdx.type = "Cl-";

            pos = filename.find ("Br-");
            if (pos != std::string::npos) newdx.type = "Br-";

            pos = filename.find ("I-");
            if (pos != std::string::npos) newdx.type = "I-";

            pos = filename.find ("Na+");
            if (pos != std::string::npos) {
                newdx.type = "Na+";
                valence = 1;
            }
            pos = filename.find ("K+");
            if (pos != std::string::npos) {
                newdx.type = "K+";
                valence = 1;
            }

            pos = filename.find ("Rb+");
            if (pos != std::string::npos) {
                newdx.type = "Rb+";
                valence = 1;
            }
            pos = filename.find ("Cs+");
            if (pos != std::string::npos) {
                newdx.type = "Cs+";
                valence = 1;
            }
            pos = filename.find ("Li+");
            if (pos != std::string::npos) {
                newdx.type = "Li+";
                valence = 1;
            }
            pos = filename.find ("Mg2+");
            if (pos != std::string::npos) {
                newdx.type = "Mg2+";
                valence = 2;
            }
            pos = filename.find ("Ca2+");
            if (pos != std::string::npos) {
                newdx.type = "Ca2+";
                valence = 2;
            }
            pos = filename.find ("Sr2+");
            if (pos != std::string::npos) {
                newdx.type = "Sr2+";
                valence = 2;
            }
            pos = filename.find ("Ba2+");
            if (pos != std::string::npos) {
                newdx.type = "Ba2+";
                valence = 2;
            }

            read_dx (path, newdx);
            dx.push_back(newdx);
        }
    }
    closedir (dir);

    // Swap Ow and Hw to the beginning, for output reading convenience
    for (size_t i = 0; i < dx.size(); i++) {
        if ((dx[i].type == "Ow") and (i > 0)) {
            dx_type tmp = dx[0];
            dx[0] = dx[i];
            dx[i] = tmp;
        }
        if ((dx[i].type == "Hw") and (i != 1)) {
            dx_type tmp = dx[1];
            dx[1] = dx[i];
            dx[i] = tmp;
        }
    }

    // Set concentration depending on valence
    for (size_t i = 2; i < dx.size(); i++){
        if (valence == 1)
            dx[i].conc = conc_salt;
        else if (valence == 2) {
            if ((dx[i].type == "Sr2+") or (dx[i].type == "Mg2+") or
                (dx[i].type == "Ca2+") or (dx[i].type == "Ba2+"))
                dx[i].conc = conc_salt;
            else if ((dx[i].type == "Cl-") or (dx[i].type == "F-") or
                     (dx[i].type == "Br-") or (dx[i].type == "I-"))
                dx[i].conc = 2*conc_salt;
        }
    }
}

////////////////////////////////////////////////////////////////
// Group the grid indexes into excluded volume, hydration ....
////////////////////////////////////////////////////////////////
void generate_list (const vector<dx_type> &dx,
                    const vector<coordinate> &pdb_coord,
                    double cutoff,
                    bool isOff_cutoff, bool &isDecomp, bool isExclV,
                    vector< vector<list_cutoff> > &list_exclV,
                    vector< vector<list_cutoff> > &list_hyd,
                    vector<list_cutoff> &list_excl_coion,
                    bool isCoion, int &index_coion) {

    vector< vector<list_cutoff> > list;
    bool isSimilar = check_dx (dx);

    double ex_int = 0;
    if (isCoion) {
        for (size_t i = 0; i < dx.size(); i++)
        if ((dx[i].type != "Ow") and (dx[i].type != "Hw")) {
            double excess_integral = integral_dx(dx[i]) - dx[i].value.size();
            if (excess_integral < 0) {
                ex_int = excess_integral;
                index_coion = i;
                break;
            }
        }
    }

    // Making cutoff list
    // Heterogeneous grid
    if (isSimilar) {
        if (not isDecomp) {
            // Currently don't support merging heterogeneous grids into a single excess electron grid
            cout << "Using heterogeneous grids -> Auto switch to decomp = 1\n";
            isDecomp = 1;
        }
        list.resize(dx.size());
        if (isExclV) {
            list_exclV.resize(dx.size());
            for (size_t i = 0; i < dx.size(); i++) {
                bool flag = 0;
                if (i > 0)
                    for (int j = i-1; j >= 0; j--)
                        if ((dx[i].ngrid == dx[j].ngrid) and (dx[i].origin == dx[j].origin)) {
                            // These two grids are essential the same
                            flag = 1;
                            list_exclV[i] = list_exclV[j];
                        }
                if (not flag)
                    list_exclV[i] = excluded_volume_list (dx[i], pdb_coord);
//                    list_exclV[i] = dx_list_cutoff (dx[i], pdb_coord, max_pdb, min_pdb, -1);
            }
        }

        if (isOff_cutoff)
            for (size_t i = 0; i < dx.size(); i++)
                // Only compute distance for co-ion grid
                if (int(i) != index_coion)
                    for (size_t j = 0; j < dx[i].value.size(); j++) {
                        list_cutoff append_this;
                        append_this.index = j;
                        append_this.dist  = 0;
                        list[i].push_back(append_this);
                    }
                else list[i] = dx_list_cutoff (dx[i], pdb_coord, INF);
        else
            for (size_t i = 0; i < dx.size(); i++) {
                bool flag = 0;
                if (i > 0)
                    for (int j = i-1; j >=0; j--)
                        if ((dx[i].ngrid == dx[j].ngrid) and (dx[i].origin == dx[j].origin)) {
                            // These two grids are essential the same
                            flag = 1;
                            list[i] = list[j];
                        }
                if (not flag)
                    list[i] = dx_list_cutoff (dx[i], pdb_coord, cutoff);
            }
        if (index_coion > -1)
            list_excl_coion = list[index_coion];

    // Homogeneous grids
    } else {
        list.resize(1);
        if (isExclV)
//            list_exclV.push_back (dx_list_cutoff (dx[0], pdb_coord, max_pdb, min_pdb, -1));
            list_exclV.push_back (excluded_volume_list (dx[0], pdb_coord));
        if (isOff_cutoff)
            if (index_coion == -1)
                for (size_t i = 0; i < dx[0].value.size(); i++) {
                    list_cutoff append_this;
                    append_this.index = i;
                    append_this.dist  = 0;
                    list[0].push_back(append_this);
                }
            else list[0] = dx_list_cutoff (dx[0], pdb_coord, INF);
        else list[0] = dx_list_cutoff (dx[0], pdb_coord, cutoff);
        if (index_coion > -1)
            list_excl_coion = list[0];
    }

    // Generate list of indexes outside the exclV
    if (isExclV)
        for (size_t i = 0; i < list.size(); i++) {
            vector<list_cutoff> diff;
            set_difference (list[i].begin(), list[i].end(), list_exclV[i].begin(), list_exclV[i].end(),
                            back_inserter(diff), comparatorStruct());
            list_hyd.push_back (diff);
        }
    else list_hyd = list;

    // Generate list of excluded co-ion
    if (index_coion > -1) {
        sort (list_excl_coion.begin(), list_excl_coion.end(), &dist_sort);
        size_t newsize = round(-ex_int);
        if (newsize > list_excl_coion.size()) {
            cout << "Need to use a bigger box for co-ion!!!\n";
            exit (0);
        }
        list_excl_coion.resize (newsize);
    }
}

////////////////////////////////////
// Sort the list based on distance
///////////////////////////////////
bool dist_sort (const list_cutoff &left,
                const list_cutoff &right) {
    return left.dist < right.dist;
}

////////////////////////////////////////
// Calculate the total excess amplitude
////////////////////////////////////////
vector< complex<double> > saxs_rism_amplitude (const vector<dx_type> &dx,
                                               bool isExcess,
                                               const vector<coordinate> &pdb_coord,
                                               const vector< vector<list_cutoff> > &list_hyd,
                                               const vector< vector<list_cutoff> > &list_exclV,
                                               const vector<list_cutoff> &list_excl_coion,
                                               int index_coion, double q,
                                               const coordinate &Leb,
                                               bool isFlex, double anom_f, bool isExpli, bool isExclV,
                                               const vector<coeff_f> &F_table) {

    vector< complex<double> > result;
    size_t size_increase = 2;
    if (isExclV)
        size_increase++;
    if (index_coion > -1)
        size_increase++;
    result.resize (dx.size() + size_increase, complex<double> (0, 0));

    result[0] = form_factor (pdb_coord, q, Leb, anom_f, isFlex, isExpli, F_table);
    if (isExclV)
        result[1] = result[0];
    result.back() = result[0];

    for (size_t j = 0; j < dx.size(); j++) {
        double atomic_factor;
        if (isExcess)     // Excess electron map
            atomic_factor = 1;
        else if ((isExpli) or (dx[j].type != "Hw")) {    // Not account for Hw grid in the implicit case
            atomic_factor = f_atm (dx[j].type, q, 0, isExpli, 1, F_table);
            if ((dx[j].type == "Rb+") or (dx[j].type == "Sr2+") or (dx[j].type == "Br-"))
                atomic_factor += anom_f;
        }
        if ((isExpli) or (dx[j].type != "Hw")) {   // Not spend time for Hw grid in the implicit case
            coordinate q_vector;
            q_vector.x = Leb.x * q;
            q_vector.y = Leb.y * q;
            q_vector.z = Leb.z * q;

            complex<double> grid_hyd;
            if (list_hyd.size() == 1)
                grid_hyd = grid_factor (dx[j], int(isExcess), atomic_factor, list_hyd[0], q_vector);
            else
                grid_hyd = grid_factor (dx[j], int(isExcess), atomic_factor, list_hyd[j], q_vector);

            if (isExclV) {
                complex<double> grid_exclV;
                if (list_exclV.size() == 1)
                    grid_exclV = grid_factor (dx[j], int(isExcess), atomic_factor, list_exclV[0], q_vector);
                else
                    grid_exclV = grid_factor (dx[j], int(isExcess), atomic_factor, list_exclV[j], q_vector);

                result[1] = result[1] + grid_exclV;
                result[j+2] = grid_hyd;
                result.back() = result.back() + grid_hyd + grid_exclV;
            } else {
                result[j+1] = grid_hyd;
                result.back() = result.back() + grid_hyd;
            }
            if (int(j) == index_coion)
                result[result.size()-2] = grid_factor (dx[j], -1, atomic_factor, list_excl_coion, q_vector);
        }
    }
    return result;
}

/////////////////////////////////////////////////////////
// Integrate over the sphere using Lebedev quadrature
////////////////////////////////////////////////////////
vector<double> saxs_rism_I (const vector<dx_type> &dx,
                            bool isExcess,
                            const vector<coordinate> &pdb_coord,
                            const vector< vector<list_cutoff> > &list_hyd,
                            const vector< vector<list_cutoff> > &list_exclV,
                            const vector<list_cutoff> &list_excl_coion,
                            int index_coion, double q, size_t rule,
                            bool isFlex, double anom_f, bool isExpli, bool isExclV,
                            const vector<coeff_f> &F_table,
                            vector<DoubPair> &dphase,
                            vector<DoubPair> &averaged_phase) {

    vector<double> intensity;
    size_t size_increase = 2;
    if (isExclV)
        size_increase++;
    if (index_coion > -1)
        size_increase++;
    intensity.resize (dx.size() + size_increase, 0);

    if (q > 0) {
        size_t Npoint = order_table (rule);     // Number of points in the unit sphere, for integration using Lebedev quadrature

        // Generate points on the unit sphere; x,y,z coordinates; w weight
        double *w_leb = new double[Npoint];
        double *x_leb = new double[Npoint];
        double *y_leb = new double[Npoint];
        double *z_leb = new double[Npoint];
        ld_by_order (Npoint, x_leb, y_leb, z_leb, w_leb);

        averaged_phase.resize (intensity.size(), make_pair(0, 0));
        dphase.resize (dx.size() + 1, make_pair(0, 0));

        vector< vector<double> > phase_matrix, dphase_matrix;  // To save all phase and phase difference information
        phase_matrix.resize  (Npoint);
        dphase_matrix.resize (Npoint);

        #pragma omp parallel for schedule (dynamic) shared (isExcess, q, x_leb, y_leb, z_leb, w_leb, anom_f,\
                                                            index_coion, isExpli, isExclV, isFlex,\
                                                            averaged_phase, dphase, phase_matrix, dphase_matrix, intensity)
        for (size_t i = 0; i < Npoint; i++)
            // Only need to compute I for one hemisphere since I(q) = I(-q)
            if (z_leb[i] >= 0.) {
                unsigned weight = 1;
                if (z_leb[i] > 0.)
                    weight = 2;

                coordinate Leb;
                Leb.x = x_leb[i];   Leb.y = y_leb[i];   Leb.z = z_leb[i];
                vector< complex<double> > ampl = saxs_rism_amplitude (dx, isExcess, pdb_coord, list_hyd, list_exclV,
                                                                      list_excl_coion, index_coion, q, Leb, isFlex,
                                                                      anom_f, isExpli, isExclV, F_table);

                // Compute phase for each amplitude
                for (size_t j = 0; j < ampl.size(); j++) {
                    double p = arg (ampl[j]) / PI * 180;
                    phase_matrix[i].push_back (p);
                    p = p * weight * w_leb[i];

                    #pragma omp atomic
                        averaged_phase[j].first += p;
                }

                // Compute phase difference between Solu and each dx
                complex<double> sum_solv (0, 0);
                complex<double> solu_ampl;
                size_t c;
                if (isExclV) {
                    c = 2;
                    solu_ampl = ampl[1];
                } else {
                    c = 1;
                    solu_ampl = ampl[0];
                }
                for (size_t j = 0; j < dx.size(); j++) {
                    complex<double> ampl_grid = ampl[j+c];
                    sum_solv = sum_solv + ampl_grid;
                    double dot_product = real (conj(solu_ampl) * ampl_grid);
//                    double dot_product = real(solu_ampl)*real(ampl_grid) + imag(solu_ampl)*imag(ampl_grid);
                    dphase_matrix[i].push_back (acos (dot_product / (abs(solu_ampl) * abs(ampl_grid))) / PI * 180);
                }
                dphase_matrix[i].push_back (acos (real (conj(solu_ampl) * sum_solv) / (abs(solu_ampl) * abs(sum_solv))) / PI * 180);
//                phase_diff[i] = .5*(norm(ampl.back()) - norm(solu_ampl) - norm(sum_solv)) / (abs(solu_ampl) * abs(sum_solv));

                for (size_t j = 0; j < dphase_matrix[i].size(); j++) {
                    double gain = dphase_matrix[i][j] * weight * w_leb[i];

                    #pragma omp atomic
                    	dphase[j].first += gain;


                }

                // Compute intensity for every amplitude by spherical averaging
                for (size_t j = 0; j < ampl.size(); j++) {
                    double gain = w_leb[i] * weight * norm(ampl[j]);

                    #pragma omp atomic
                        intensity[j] += gain;
                }
            }

        // Now calculate s^2 for all the phase
        for (size_t i = 0; i < Npoint; i++)
            if (z_leb[i] >= 0.) {
                unsigned weight = 1;
                if (z_leb[i] > 0)
                    weight = 2;
                for (size_t j = 0; j < averaged_phase.size(); j++) {
                    double phase_diff = phase_matrix[i][j] - averaged_phase[j].first;
                    averaged_phase[j].second += phase_diff*phase_diff * weight * w_leb[i];
                }
                for (size_t j = 0; j < dphase.size(); j++) {
                    double dphase_diff = dphase_matrix[i][j] - dphase[j].first;
                    dphase[j].second += dphase_diff*dphase_diff * weight * w_leb[i];
                }
            }

        // Back to s
        for (size_t i = 0; i < averaged_phase.size(); i++)
            averaged_phase[i].second = sqrt (averaged_phase[i].second);
        for (size_t i = 0; i < dphase.size(); i++)
            dphase[i].second = sqrt (dphase[i].second);
    } else {        // q = 0 case
        coordinate Leb;
        Leb.x = 0;     Leb.y = 0;     Leb.z = 0;
        vector< complex<double> > ampl = saxs_rism_amplitude (dx, isExcess, pdb_coord, list_hyd, list_exclV,
                                                              list_excl_coion, index_coion, 0, Leb, isFlex,
                                                              anom_f, isExpli, isExclV, F_table);
        // Compute intensity
        for (size_t j = 0; j < ampl.size(); j++)
            intensity[j] = norm(ampl[j]);
        averaged_phase.resize (intensity.size(), make_pair (0, 0));
        dphase.resize         (dx.size() + 1,    make_pair (0, 0));
    }
    // Since this is averaging, not integrating
    // we do not multiply by 4*PI, aka the two 4PI cancel each other
    return intensity;
}
