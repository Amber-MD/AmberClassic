#include <iostream>
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

void gen_excess_e_map (const vector< vector<coordinate> > &snapshot,
                       dx_type &dx,
                       const vector<list_cutoff> &list_exclV,
                       const vector<list_cutoff> &list_hyd,
                       const coordinate &max_bulk,
                       const coordinate &min_bulk,
                       double smear_cut,
                       double anom_f,
                       const vector<coeff_f> &F_table);

//////////////////////////////////////////////////////////////////////////
int main (int argc,
          char *argv[]) {
    string pdb_file, pddf_outfile;
    double cutoff = 100.;
    double dr = 2.0;
    double smear_cut = 1.0;
    double dcutoff = 20;
    size_t type = 1;
    bool emap = 0;
    double anom_f = 0;
    unsigned ncpus = 0;

    int option_char;
    do {
        static struct option long_options[] =
        {
            {"pdb",         required_argument,  0,  'p'},
            {"cutoff",      required_argument,  0,  'c'},
            {"emap",        no_argument,        0,  'e'},
            {"smear_cut",   required_argument,  0,  's'},
            {"bulk_cut",    required_argument,  0,  'b'},
            {"dr",          required_argument,  0,  'r'},
            {"type",        required_argument,  0,  't'},
            {"output",      required_argument,  0,  'o'},
            {"ncpus",       required_argument,  0,  'n'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        option_char = getopt_long (argc, argv, "p:b:c:es:r:t:o:n:", long_options, &option_index);	// followed by 1 colon - require an argument; no colon - not argument required
        if (option_char == -1)
            pddf_md_usage();
        else switch (option_char) {
            case '?':
                pddf_md_usage();
                exit (0);
            case 'p':
                pdb_file = optarg;
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
            case 'b':
                dcutoff = atof (optarg);
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
    if ((type < 1) or (type > 3)) {
        cout << "type must be 1, 2 or 3!!!!\n";
        exit (0);
    }

    ofstream OUTPUT (pddf_outfile.c_str());
    if (OUTPUT.is_open()) {
        cout << "Reading input file ...\n";

        vector< vector<coordinate> > snapshot;
        vector<unsigned> weight;
        read_pdb_weight (pdb_file, snapshot, weight); // reading weight, but not using them

        vector< vector<coordinate> > solu;
        solu.resize(snapshot.size());

        coordinate max_solu, min_solu;
        max_solu.x = -INF;  max_solu.y = -INF;  max_solu.z = -INF;
        min_solu.x =  INF;  min_solu.y =  INF;  min_solu.z =  INF;

        for (size_t i = 0; i < snapshot.size(); i++) {
            solu[i] = solute_coord (snapshot[i]);
            coordinate max_model, min_model;
            maxmin_coord_pdb (solu[i], max_model, min_model);
            max_solu.x = max (max_model.x, max_solu.x);
            max_solu.y = max (max_model.y, max_solu.y);
            max_solu.z = max (max_model.z, max_solu.z);

            min_solu.x = min (min_model.x, min_solu.x);
            min_solu.y = min (min_model.y, min_solu.y);
            min_solu.z = min (min_model.z, min_solu.z);
        }

        coordinate min_bulk, max_bulk;
        min_bulk.x = min_solu.x - dcutoff;
        min_bulk.y = min_solu.y - dcutoff;
        min_bulk.z = min_solu.z - dcutoff;

        max_bulk.x = max_solu.x + dcutoff;
        max_bulk.y = max_solu.y + dcutoff;
        max_bulk.z = max_solu.z + dcutoff;

        cout << "Separating excluded volume and hydration lists ...\n";
        vector<list_cutoff> list_exclV, list_hyd;
        coordinate max_snapshot, min_snapshot;
        maxmin_coord_traj (snapshot, max_snapshot, min_snapshot);
        coordinate res;
        res.x = dr;    res.y = dr;    res.z = dr;
        dx_type dx = gen_grid (max_snapshot, min_snapshot, res);
        gen_exV_hyd_list (solu[0], dx, list_exclV, list_hyd);

        // Generate solute and solvent excess electron density grid
        vector<coeff_f> F_table = assign_F_table ();
        gen_excess_e_map (snapshot, dx, list_exclV, list_hyd, max_bulk, min_bulk, smear_cut, anom_f, F_table); 
        snapshot.clear();

        if (emap) {
            cout << "Outputting excess electron map ...\n";
            string emap_file = "emap.dx";
            write_dx (dx, emap_file);
        }

        cout << "Computing the PDDF ...\n";
        vector<double> pddf;
        double dvol = dx.delta[0] * dx.delta[1] * dx.delta[2];
        double scaling = dvol * dvol;
        if (type == 1)
            pddf = calc_pddf (dx, list_exclV, dx, list_hyd, cutoff, dr);
        else if (type == 2)
            pddf = calc_pddf (dx, list_hyd, dx, list_hyd, cutoff, dr);
        else if (type == 3)
            pddf = calc_pddf (dx, list_exclV, dx, list_exclV, cutoff, dr);

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

////////////////////////////////////////////////////////////
// Generate an excess electron density map from coordinates
////////////////////////////////////////////////////////////
void gen_excess_e_map (const vector< vector<coordinate> > &snapshot,
                       dx_type &dx,
                       const vector<list_cutoff> &list_exclV,
                       const vector<list_cutoff> &list_hyd,
                       const coordinate &max_bulk,
                       const coordinate &min_bulk,
                       double smear_cut,
                       double anom_f,
                       const vector<coeff_f> &F_table) {

    #pragma omp parallel for schedule (dynamic, 1000) shared (dx, smear_cut, anom_f) collapse (2)
    for (size_t i = 0; i < snapshot.size(); i++)
        for (size_t j = 0; j < snapshot[i].size(); j++)
            if (isSmear (snapshot[i][j], max_bulk, min_bulk))
                distributing_e (snapshot[i][j], dx, smear_cut, anom_f, F_table);
            else {
                int x_atom = floor((snapshot[i][j].x - dx.origin[0]) / dx.delta[0]);
                int y_atom = floor((snapshot[i][j].y - dx.origin[1]) / dx.delta[1]);
                int z_atom = floor((snapshot[i][j].z - dx.origin[2]) / dx.delta[2]);
                int index = dx_3Dindex_to_1Dindex (dx, x_atom, y_atom, z_atom);
                double rho_e = count_e (snapshot[i][j].type, anom_f) / (dx.delta[0] * dx.delta[1] * dx.delta[2]);
                #pragma omp atomic
                    dx.value[index] += rho_e;
            }

    double rho_bulk = calc_bulk_rho_e (dx, max_bulk, min_bulk) / snapshot.size();
    cout << "    rho_bulk   " << rho_bulk << "   e/A^3\n";

    #pragma omp parallel for shared (dx, rho_bulk) collapse (3)
    for (size_t x = 0; x < dx.ngrid[0]; x++)
        for (size_t y = 0; y < dx.ngrid[1]; y++)
            for (size_t z = 0; z < dx.ngrid[2]; z++) {
                size_t i = dx_3Dindex_to_1Dindex (dx, x, y, z);
                if ((x == 0) or (x == dx.ngrid[0] - 1) or\
                    (y == 0) or (y == dx.ngrid[1] - 1) or\
                    (z == 0) or (z == dx.ngrid[2] - 1))
                    dx.value[i] = 0.;           // Get rid of boundary points
                else dx.value[i] = dx.value[i] / snapshot.size() - rho_bulk;
            }
}
