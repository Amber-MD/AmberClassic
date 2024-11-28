/* A program to compute SAXS and ASAXS profiles from MD simulation
   More details are in Park et al JCP 2009; 130; 134114
                       Nguyen et al JCP 2014; 141; 22D508

   Written by Hung Nguyen, Case's group, 2013 */

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <cmath>
#include <complex>
#include <vector>
#include <cstdlib>
#include <limits>
#include <getopt.h>

#include "saxsDS.h"
#include "lebedev.h"
#include "pdb.h"
#include "atm_f.h"
#include "f_calc.h"
#include "saxs_md.h"
#include "IO.h"
#include "const.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

int main (int argc,
          char *argv[]) {

    // Declare variables and their default values
    double qcut = 0.5;        // Cutoff q
    double anom_f = 0;        // Variable controls whether on-edge or off-edge calculation
    double dcutoff = 10.;
    bool isExpli = 0;         // Account for isExplicit H atoms in pdb file
    bool isTight = 0;         // Use isTighter convergence for Lebedev quadrature
//    bool corr = 0;          // Using corrected atomic factor for water as in J Chem Phys 2000, 113, 9149
    double bulk_corr = 0.;    // Correction for the difference in bulk density of water between the "sample" and "blank" simulation
                              // The true excess density outside of bulk_cutoff will be computed as g = (1-bulk_corr)*rho_wat(sample) - rho_wat(blank)
    double bulk_cutoff = 10.;
    double dq = 0.01;
    string exper, pdb_solu, pdb_solv, outfile;
    unsigned ncpus = 0;

    int option_char;
    do {
        static struct option long_options[] =
        {
            {"help",        no_argument,        0,    'h'},
            {"anom_f",      required_argument,  0,    'a'},
            {"system",      required_argument,  0,    'i'},
            {"solvent",     required_argument,  0,    'w'},
            {"qcut",        required_argument,  0,    'q'},
            {"dq",          required_argument,  0,    'd'},
            {"cutoff",      required_argument,  0,    'c'},
            {"isExpli",     no_argument,        0,    'e'},
            {"corrected",   required_argument,  0,    'k'},
            {"bcutoff",     required_argument,  0,    'b'},
            {"exper",       required_argument,  0,    'x'},
            {"output",      required_argument,  0,    'o'},
            {"ncpus",       required_argument,  0,    'n'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        option_char = getopt_long (argc, argv, "ha:i:w:q:d:c:eb:k:o:n:", long_options, &option_index);
        // followed by 1 colon - required an argument; 2 colon - not required argument
        if (option_char == -1)
            saxs_md_usage();
        else switch (option_char) {
            case 'h':
            case '?':
                saxs_md_usage();
                exit (0);
            case 'a':
                anom_f = atof (optarg);
                break;
            case 'i':
                pdb_solu = optarg;
                break;
            case 'w':
                pdb_solv = optarg;
                break;
            case 'q':
                qcut = atof (optarg);
                break;
            case 'd':
                dq = atof (optarg);
                break;
            case 'c':
                dcutoff = atof (optarg);
                break;
            case 'e':
                isExpli = 1;
                break;
            case 'k':
                bulk_corr = atof (optarg);
                break;
            case 'b':
                bulk_cutoff = atof (optarg);
                break;
            case 'x':
                exper = optarg;
                break;
            case 'o':
                outfile = optarg;
                break;
            case 'n':
                ncpus = atoi (optarg);
                break;
        }
    } while (option_char != -1);

    ofstream OUTPUT (outfile.c_str());
    if (OUTPUT.is_open()) {
        cout << "Reading pdb ...\n";
        vector< vector<coordinate> > solu_box, solv_box;
        vector<unsigned> weight_solu, weight_solv;
        read_pdb_weight (pdb_solu, solu_box, weight_solu);
        read_pdb_weight (pdb_solv, solv_box, weight_solv);

        vector<double> q;
        if (exper.empty())
            for (size_t i = 0; i<= floor(qcut/dq); i++)
                q.push_back(i*dq);
        else
            read_exp (exper, q);

        vector<double> intensity = saxs_md_calc (solu_box, solv_box, weight_solu, weight_solv, q,
                                                 bulk_corr, bulk_cutoff, anom_f, qcut, dq, dcutoff, isExpli, isTight, ncpus);
        saxs_md_output (intensity, dq, pdb_solu, pdb_solv, dcutoff,
                        bulk_corr, bulk_cutoff, anom_f, isExpli, isTight, OUTPUT);
        OUTPUT.close();
    } else {
        cout << "Unable to write to file " << outfile << endl;
        exit (0);
    }
    return 0;
}

//////////////////////// END MAIN ///////////////////////////

vector<double> saxs_md_calc (vector< vector<coordinate> > &solu_box,
                             vector< vector<coordinate> > &solv_box,
                             const vector<unsigned> &weight_solu,
                             const vector<unsigned> &weight_solv,
                             const vector<double> &q,
                             double bulk_corr, double bulk_cutoff,
                             double anom_f, double qcut, double dq,
                             double dcutoff, bool isExpli, bool isTight,
                             unsigned ncpus) {
#ifdef _OPENMP
    omp_set_dynamic(0);
    if (ncpus == 0)
        ncpus = omp_get_max_threads();
    omp_set_num_threads(ncpus);
#endif

    /* Currently assume :
       + weight_solu and weight_solv are equal, this is needed for the statistics later
       + sizes of solu and solv must be equal*/

    if (solu_box.size() != solv_box.size()) {
        cout << "!!!!!!!!!!!  Total snapshots of solute and solvent are not equal   !!!!!!!!!!!!!\n";
        cout << "!!!!!!!!!!!  Truncation   ";
        if (solu_box.size() > solv_box.size()) {
            cout << "SOLUTE box\n";
            solu_box.resize (solv_box.size());
        } else {
            cout << "SOLVENT box\n";
            solv_box.resize (solu_box.size());
        }
    }
    //for (size_t i = 0; i < weight_solu.size(); i++)
    //    weight_solv[i] = weight_solu[i];

    if (not isExpli) {
        cout << "Merging H atoms ...\n";

        #pragma omp parallel for shared (solu_box, solv_box)
        for (size_t i = 0; i < solu_box.size(); i++) {
            mergeH (solu_box[i]);
            mergeH (solv_box[i]);
        }
    }
    cout << "Stripping ...\n";
    vector< vector<coordinate> > solu;
    solu.resize(solu_box.size());

    coordinate max_solu, min_solu;
    max_solu.x = -INF;  max_solu.y = -INF;  max_solu.z = -INF;
    min_solu.x =  INF;  min_solu.y =  INF;  min_solu.z =  INF;
    for (size_t i = 0; i < solu.size(); i++) {
        solu[i] = solute_coord (solu_box[i]);
        coordinate max_model, min_model;
        maxmin_coord_pdb (solu[i], max_model, min_model);
        max_solu.x = max (max_model.x, max_solu.x);
        max_solu.y = max (max_model.y, max_solu.y);
        max_solu.z = max (max_model.z, max_solu.z);

        min_solu.x = min (min_model.x, min_solu.x);
        min_solu.y = min (min_model.y, min_solu.y);
        min_solu.z = min (min_model.z, min_solu.z);
    }

    vector< vector<coordinate> > solu_strip, solv_strip;
    solu_strip.resize(solu_box.size());
    solv_strip.resize(solv_box.size());

    #pragma omp parallel for schedule(dynamic) shared (solu_box, solv_box, solu_strip, solv_strip, max_solu, min_solu, dcutoff)
    for (size_t i = 0; i < solu_box.size(); i++) {
        solu_strip[i] = strip_buffer_atom (solu_box[i], dcutoff, max_solu, min_solu);
        solv_strip[i] = strip_buffer_atom (solv_box[i], dcutoff, max_solu, min_solu);
    }
    solu_box.clear();   solv_box.clear();

    vector<double> intensity;
    vector<coeff_f> F_table = assign_F_table ();

    coordinate max_buffer_cut;
    coordinate min_buffer_cut;
    max_buffer_cut.x = max_solu.x + bulk_cutoff;
    max_buffer_cut.y = max_solu.y + bulk_cutoff;
    max_buffer_cut.z = max_solu.z + bulk_cutoff;
    min_buffer_cut.x = min_solu.x - bulk_cutoff;
    min_buffer_cut.y = min_solu.y - bulk_cutoff;
    min_buffer_cut.z = min_solu.z - bulk_cutoff;
    for (size_t i = 0; i < q.size(); i++) {
        size_t rule;
        if (not isTight)
            rule = 4 + floor (q[i]/.04);
        else
            rule = 5 + floor (q[i]/.03);

        // Maximum rule is 65, rarely use up to this number though
        while ((available_table(rule) == 0) and (rule < 65))
            rule++;

        double I = saxs_md_I (solu_strip, weight_solu, solv_strip, weight_solv, q[i],
                              bulk_corr, max_buffer_cut, min_buffer_cut, rule, anom_f, isExpli, F_table);
        intensity.push_back (I);
    }
    return intensity;
}

///////////////////////////////////////////////////////////////////////////////////////
// Compute the mean of a complex vector with weights, v and w must have the same sizes
///////////////////////////////////////////////////////////////////////////////////////
complex<double> mean_complex (const vector< complex<double> > &v,
                              const vector<unsigned> &w) {

    complex<double> mean (0, 0);
    unsigned total_weight = 0;
    for (size_t i = 0; i < v.size(); i++) {
        double weight = w[i];
        mean = mean + weight*v[i];
        total_weight += w[i];
    }
    double inv_size = 1.0/total_weight;
    return mean*inv_size;
}
/////////////////////////////////////////////////////////////////////////////
// Compute D11(q), as described in Park et al, J Chem Phys 2009, 130, 134114
/////////////////////////////////////////////////////////////////////////////
double saxs_md_D11 (const vector< vector <coordinate> > &solu,
                    const vector<unsigned> &weight_solu,
                    const vector< vector <coordinate> > &solv,
                    const vector<unsigned> &weight_solv,
                    double q, double bulk_corr,
                    const coordinate &max_buffer_cut,
                    const coordinate &min_buffer_cut,
                    const coordinate &Leb,
                    double anom_f, bool isExpli,
                    const vector<coeff_f> &F_table) {

    vector< complex<double> > solu_f, solv_f;
    for (size_t i = 0; i < solu.size(); i++) {
        complex<double> f_solu (0, 0);
        for (size_t j = 0; j < solu[i].size(); j++)
            if ((isExpli) or (solu[i][j].type != "H")) {
                double atomic_factor;
                atomic_factor = f_atm (solu[i][j].type, q, solu[i][j].nHyd, isExpli, 0, F_table);
                if ((solu[i][j].type == "Rb+") or (solu[i][j].type == "Cs+") or (solu[i][j].type == "Br-"))
                    atomic_factor += anom_f;
                double qr = q * (Leb.x*solu[i][j].x + Leb.y*solu[i][j].y + Leb.z*solu[i][j].z);
                complex<double> f = atomic_factor * exp(complex<double> (0,1) * qr);
                if ((abs(bulk_corr) > std::numeric_limits<double>::epsilon()) and \
                   ((solu[i][j].x < min_buffer_cut.x) or (solu[i][i].y < min_buffer_cut.y) or (solu[i][j].z < min_buffer_cut.z) or \
                    (solu[i][j].x > max_buffer_cut.x) or (solu[i][j].y > max_buffer_cut.y) or (solu[i][j].z > max_buffer_cut.z)))
                    f = f * (1. - bulk_corr);
                f_solu = f_solu + f;
            }
        solu_f.push_back (f_solu);
    }
    for (size_t i = 0; i < solv.size(); i++)
        solv_f.push_back (form_factor (solv[i], q, Leb, anom_f, 0, isExpli, F_table));

    complex<double> ensbl_solu = mean_complex (solu_f, weight_solu);
    complex<double> ensbl_solv = mean_complex (solv_f, weight_solv);

    unsigned totalw_solu = 0;
    unsigned totalw_solv = 0;
    for (size_t i = 0; i < weight_solu.size(); i++)
        totalw_solu += weight_solu[i];
    for (size_t i = 0; i < weight_solv.size(); i++)
        totalw_solv += weight_solv[i];

    double diff_solu = 0;
    for (size_t i = 0; i < solu.size(); i++)
        diff_solu += weight_solu[i] * norm (solu_f[i] - ensbl_solu);
    diff_solu /= totalw_solu;

    double diff_solv = 0;
    for (size_t i = 0; i < solv.size(); i++)
        diff_solv += weight_solv[i] * norm (solv_f[i] - ensbl_solv);
    if (totalw_solv > 1)
        diff_solv *= (totalw_solv + 1) / (totalw_solv*(totalw_solv - 1));

    return norm(ensbl_solu - ensbl_solv) + diff_solu - diff_solv;
}

/////////////////////////////////////////////////////////
// Integrate over the sphere using Lebedev quadrature
////////////////////////////////////////////////////////
double saxs_md_I (const vector< vector<coordinate> > &solu,
                  const vector<unsigned> &weight_solu,
                  const vector< vector<coordinate> > &solv,
                  const vector<unsigned> &weight_solv,
                  double q, double bulk_corr,
                  const coordinate &max_buffer_cut,
                  const coordinate &min_buffer_cut,
                  size_t rule, double anom_f, bool isExpli,
                  const vector<coeff_f> &F_table) {

    double intensity = 0;
    if (q > 0) {
        int Npoint = order_table (rule);    // Number of points in the unit sphere, for integration using Lebedev quadrature

        // Generate points on the unit sphere; x,y,z coordinates; w weight
        double *w_leb = new double[Npoint];
        double *x_leb = new double[Npoint];
        double *y_leb = new double[Npoint];
        double *z_leb = new double[Npoint];
        ld_by_order (Npoint, x_leb, y_leb, z_leb, w_leb);

        #pragma omp parallel for schedule (dynamic) shared (q, bulk_corr, rule, anom_f, isExpli, x_leb, y_leb, z_leb, w_leb) reduction (+ : intensity)
        for (size_t i = 0; i < Npoint; i++)
            // Only need to compute I for one hemisphere since I(q) = I(-q)
            if (z_leb[i] >= 0.) {
                size_t scale = 1;
                if (z_leb[i] > 0.)
                    scale = 2;
                coordinate Leb;
                Leb.x = x_leb[i];    Leb.y = y_leb[i];    Leb.z = z_leb[i];
                intensity += w_leb[i] * scale * saxs_md_D11 (solu, weight_solu, solv, weight_solv, q, bulk_corr,
                                                             max_buffer_cut, min_buffer_cut, Leb, anom_f, isExpli, F_table);
            }
    } else {        // q = 0 case
        coordinate Leb;
        Leb.x = 0;     Leb.y = 0;     Leb.z = 0;
        intensity = saxs_md_D11 (solu, weight_solu, solv, weight_solv, q, bulk_corr,
                                 max_buffer_cut, min_buffer_cut, Leb, anom_f, isExpli, F_table);
    }
    return intensity;   // Since this is averaging, not integrating, one does not multiply by 4*PI, aka the two 4PI cancel each other
}
