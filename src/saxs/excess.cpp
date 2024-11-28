#include <iostream>
#include <string>
#include <vector>
#include <getopt.h>
#include <cstdlib>

#include "saxsDS.h"
#include "const.h"
#include "dx.h"
#include "pdb.h"
#include "vol_decomp.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

string dx_file, pdbfile;

double conc = 55.34;
double cutoff = 20;

// Declare functions to use
void usage ();
double excess_number (const dx_type &dx,
                      const vector<list_cutoff> &list,
                      double conc);

/////////////////////
int main (int argc, char *argv[]) {
    int option_char;
    do {
        static struct option long_options[] =
        {
            {"input",       required_argument,  0,  'i'},
            {"pdb",         required_argument,  0,  'p'},
            {"conc",        required_argument,  0,  'c'},
            {"cutoff",      required_argument,  0,  'r'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        option_char = getopt_long (argc, argv, "i:p:c:r:", long_options, &option_index);		// followed by 1 colon - required an argument; 2 colon - not required argument
        if (option_char == -1)
            usage();
        else switch (option_char) {
            case '?':
                usage();
                exit (0);
            case 'i':
                dx_file = optarg;
                break;
            case 'p':
                pdbfile = optarg;
                break;
            case 'c':
                conc = atof (optarg);
                break;
            case 'r':
                cutoff = atof (optarg);
                break;
        }
    } while (option_char != -1);

#ifdef _OPENMP
    omp_set_dynamic(0);
    unsigned ncpus = omp_get_max_threads();
    omp_set_num_threads(ncpus);
#endif

    dx_type dx;
    read_dx (dx_file, dx);
    vector<coordinate> pdb_coord;
    read_pdb (pdbfile, pdb_coord);

    vector<list_cutoff> exclV = excluded_volume_list (dx, pdb_coord);
    cout << "Excluded volume (A^3) = " << exclV.size() * dx.delta[0] * dx.delta[1] * dx.delta[2] << endl;
    vector<list_cutoff> list = dx_list_cutoff (dx, pdb_coord, cutoff);

    double N_excl = excess_number (dx, exclV, conc);
    cout << "N_exclV      = " << N_excl << endl;
    double N_excess = excess_number (dx, list, conc);
    cout << "N_solv       = " << N_excess - N_excl << endl;
    cout << "N_excess     = " << N_excess << endl;

    return 0;
}

/////////////////////
void usage () {
    cout << "Usage:    excess    -h   --help       display help\n";
    cout << "                    -i   --input      input dx file\n";
    cout << "                    -p   --pdb        pdb file\n";
    cout << "                    -c   --conc       concentration\n";
    cout << "                    -r   --cutoff     distance cutoff [default 20A]\n";
}

///////////////////////////////////////////////////////
double excess_number (const dx_type &dx,
                      const vector<list_cutoff> &list,
                      double conc) {
    double Nex = 0;

    #pragma omp parallel for shared (dx, list) reduction(+:Nex)
    for (size_t i = 0; i < list.size(); i++)
        Nex += dx.value[list[i].index] - 1;
    Nex *= AVOGADRO * 1e-27 * conc * dx.delta[0] * dx.delta[1] * dx.delta[2];
    return Nex;
}
