#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <getopt.h>

#include "saxsDS.h"
#include "dx.h"
#include "atm_f.h"
#include "pdb.h"
#include "pddf.h"
#include "IO.h"
#include "vol_decomp.h"
#include "const.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

dx_type solute_excess_e_map (const vector<coordinate> &pdb_coord,
                             const dx_type &dx,
                             const vector<list_cutoff> &list_exclV,
                             double rho_bulk,
                             double smear_cut,
                             double anom_f,
                             const vector<coeff_f> &F_table);

//////////////////////////////////////////////////////////////////////////
int main (int argc,
          char *argv[]) {
    string dx_infile, pdb_file, pddf_outfile;
    double bulk_rho = 0.3333;
    double conc = 55.34;
    double cutoff = 100.;
    vector<size_t> merge;
    merge.resize (3, 1);
    double dr = 1.0;
    double smear_cut = 1.0;
    size_t type = 1;
    bool emap = 0;
    double anom_f = 0;
    unsigned ncpus = 0;

    int option_char;
    do {
        static struct option long_options[] =
        {
            {"pdb",         required_argument,  0,  'p'},
            {"dxfile",      required_argument,  0,  'd'},
//            {"rdf",         required_argument,  0,  'f'},
            {"bulk_rho",    required_argument,  0,  'b'},
            {"conc",        required_argument,  0,  'w'},
            {"merge",       required_argument,  0,  'm'},
            {"cutoff",      required_argument,  0,  'c'},
            {"emap",        no_argument,        0,  'e'},
            {"smear_cut",   required_argument,  0,  's'},
            {"dr",          required_argument,  0,  'r'},
            {"type",        required_argument,  0,  't'},
            {"output",      required_argument,  0,  'o'},
            {"ncpus",       required_argument,  0,  'n'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        option_char = getopt_long (argc, argv, "p:d:b:w:m:c:es:r:t:o:n:", long_options, &option_index);	// followed by 1 colon - require an argument; no colon - not argument required
        if (option_char == -1)
            pddf_rism_usage();
        else switch (option_char) {
            char *end1, *end2, *end3;
            case '?':
                pddf_rism_usage();
                exit (0);
            case 'p':
                pdb_file = optarg;
                break;
            case 'd':
                dx_infile = optarg;
                break;
/*            case 'f':
                rdf_infile = optarg;
                break;*/
            case 'b':
                bulk_rho = atof (optarg);
                break;
            case 'w':
                conc = atof (optarg);
                break;
            case 'm':
                merge[0] = strtol (optarg, &end1, 10);
                merge[1] = strtol (optarg, &end2, 10);
                merge[2] = strtol (optarg, &end3, 10);
                break;
            case 'c':
                cutoff = atof (optarg);
                break;
            case 'r':
                dr = atof (optarg);
                break;
            case 't':
                type = atoi (optarg);
                break;
            case 'e':
                emap = 1;
                break;
            case 's':
                smear_cut = atof (optarg);
                break;
            case 'o':
                pddf_outfile = optarg;
                break;
            case 'n':
                ncpus = atoi (optarg);
                break;
        }
    } while (option_char != -1);

#ifdef _OPENMP
	omp_set_dynamic(0);
	if (ncpus == 0)
		ncpus = omp_get_max_threads();
	omp_set_num_threads(ncpus);
#endif

    /////////////// Sanity check ////////////////////
    if ((merge[0] < 1) or (merge[0] < 1) or (merge[0] < 1)) {
        cout << "merge minimum is 1!!!!!\n";
        exit (0);
    }
    if ((type < 1) or (type > 3)) {
        cout << "type must be 1, 2 or 3!!!!\n";
        exit (0);
    }

    ofstream OUTPUT (pddf_outfile.c_str());
    if (OUTPUT.is_open()) {
        cout << "Reading input file ...\n";
        dx_type dx;
        read_dx (dx_infile, dx);
        dx.conc = conc;
        string filename = dx_infile.substr (dx_infile.find_last_of ("/") + 1);
        double Z;
        if (filename.substr(0,3) == "guv") {
            size_t pos = filename.find ("O");
            if (pos != std::string::npos) {
                dx.type = "Ow";
                Z = 10;
            } else {
                pos = filename.find ("Cl-");
                if (pos != std::string::npos) {
                    dx.type = "Cl-";
                    Z = 18;
                } else {
                    pos = filename.find ("Br-");
                    if (pos != std::string::npos) {
                        dx.type = "Br-";
                        Z = 36;
                    } else {
                        pos = filename.find ("Na+");
                        if (pos != std::string::npos) {
                            dx.type = "Na+";
                            Z = 10;
                        } else {
                            pos = filename.find ("Rb+");
                            if (pos != std::string::npos) {
                                dx.type = "Rb+";
                                Z = 36;
                            } else {
                                pos = filename.find ("Mg2+");
                                if (pos != std::string::npos) {
                                    dx.type = "Mg2+";
                                    Z = 10;
                                } else {
                                    pos = filename.find ("Sr2+");
                                    if (pos != std::string::npos) {
                                        dx.type = "Sr2+";
                                        Z = 36;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        vector<coordinate> pdb_coord;
        read_pdb (pdb_file, pdb_coord);
/*        vector<double> e_rdf_wat;
        double rdf_dr;
        if (dx.type == "Ow")
            read_RDF (rdf_infile, rdf_dr, e_rdf_wat);*/

        cout << "Separating excluded volume and hydration lists ...\n";
        vector<list_cutoff> list_exclV_fine, list_hyd_fine;
        gen_exV_hyd_list (pdb_coord, dx, list_exclV_fine, list_hyd_fine);

        vector<coeff_f> F_table = assign_F_table();
        dx_type solute_map = solute_excess_e_map (pdb_coord, dx, list_exclV_fine, bulk_rho, smear_cut, anom_f, F_table);

        dx_type solvent_map = dx;
        // Generate excess density map
        for (size_t i = 0; i < dx.value.size(); i++)
            solvent_map.value[i] -= 1.;

        dx_type solv_coarse, solu_coarse;
        if ((merge[0] > 1) or (merge[1] > 1) or (merge[2] > 1)) {
            cout << "Coarse graining grids ...\n";
            solv_coarse = gen_coarse_dx (solvent_map, merge[0], merge[1], merge[2]);
            solu_coarse = solv_coarse;
            solu_coarse.type = "solu_coarse";
            coarse_grid (solute_map,  solu_coarse, merge[0], merge[1], merge[2]);
            coarse_grid (solvent_map, solv_coarse, merge[0], merge[1], merge[2]);
        }

        if (emap) {
            cout << "Output excess electron map ...\n";
            string solu_emap_file = "solu_emap.dx";
            string solv_emap_file = dx.type + "_emap.dx";
            if ((merge[0] > 1) or (merge[1] > 1) or (merge[2] > 1)) {
                write_dx (solu_coarse, solu_emap_file);
                write_dx (solv_coarse, solv_emap_file);
            } else
                write_dx (solute_map, solu_emap_file);
        }

        cout << "Computing the PDDF ...\n";
        vector<double> pddf;
        double scaling = 1;
        if ((merge[0] > 1) or (merge[1] > 1) or (merge[2] > 1)) {
            vector<list_cutoff> list_exclV_coarse, list_hyd_coarse;
            gen_exV_hyd_list (pdb_coord, solu_coarse, list_exclV_coarse, list_hyd_coarse);
            double dvol_solu = solu_coarse.delta[0] * solu_coarse.delta[1] * solu_coarse.delta[2];
            double dvol_solv = solv_coarse.delta[0] * solv_coarse.delta[1] * solv_coarse.delta[2];
            if (type == 1) {
                pddf = calc_pddf (solu_coarse, list_exclV_coarse, solv_coarse, list_hyd_coarse, cutoff, dr);
                scaling = Z * AVOGADRO * 1e-27 * solv_coarse.conc * dvol_solu * dvol_solv;
            } else if (type == 2) {
                pddf = calc_pddf (solv_coarse, list_hyd_coarse, solv_coarse, list_hyd_coarse, cutoff, dr);
                double tmp = Z * AVOGADRO * 1e-27 * solv_coarse.conc * dvol_solv;
                scaling = tmp * tmp;
            } else if (type == 3) {
                pddf = calc_pddf (solu_coarse, list_exclV_coarse, solu_coarse, list_exclV_coarse, cutoff, dr);
                scaling = dvol_solu * dvol_solu;
            }
        } else {
            double dvol_solu = solute_map.delta[0] * solute_map.delta[1] * solute_map.delta[2];
            double dvol_solv = solvent_map.delta[0] * solvent_map.delta[1] * solvent_map.delta[2];
            if (type == 1) {
                pddf = calc_pddf (solute_map, list_exclV_fine, solvent_map, list_hyd_fine, cutoff, dr);
                scaling = Z * AVOGADRO * 1e-27 * solvent_map.conc * dvol_solu * dvol_solv;
            } else if (type == 2) {
                pddf = calc_pddf (solvent_map, list_hyd_fine, solvent_map, list_hyd_fine, cutoff, dr);
                double tmp = Z * AVOGADRO * 1e-27 * solvent_map.conc * dvol_solv;
                scaling = tmp * tmp;
            } else if (type == 3) {
                pddf = calc_pddf (solute_map, list_exclV_fine, solute_map, list_exclV_fine, cutoff, dr);
                scaling = dvol_solu * dvol_solu;
            }
        }
        for (size_t i = 0; i < pddf.size(); i++)
            pddf[i] *= scaling;
        pddf_output (pddf, type, dr, OUTPUT);
    } else {
        cout << "Unable to write to file " << pddf_outfile << endl;
        exit (0);
    }

    return 0;
}
/////////////////////////////////////    END MAIN      /////////////////////////////////////////

/////////////////////////////////////////////////////
// Generate an excess electron density map of solute
/////////////////////////////////////////////////////
dx_type solute_excess_e_map (const vector<coordinate> &pdb_coord,
                             const dx_type &dx,
                             const vector<list_cutoff> &list_exclV,
                             double rho_bulk,
                             double smear_cut,
                             double anom_f,
                             const vector<coeff_f> &F_table) {
    dx_type result;
    result.type = "drho_solu";
    result.ngrid = dx.ngrid;    result.origin = dx.origin;  result.delta = dx.delta;
    result.value.resize (dx.value.size(), 0);

    #pragma omp parallel for shared (result, smear_cut, anom_f)
    for (size_t i = 0; i < pdb_coord.size(); i++)
        distributing_e (pdb_coord[i], result, smear_cut, anom_f, F_table);

    #pragma omp parallel for shared (result, rho_bulk)
    for (size_t i = 0; i < list_exclV.size(); i++)
        result.value[list_exclV[i].index] -= rho_bulk;

    return result;
}
