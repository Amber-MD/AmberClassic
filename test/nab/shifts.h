molecule	m;
atom		a;
residue		r;
int			i,k;                      // general loop indices
int			l_res;

#define MAXRINGS 500
#define MAXRINGPOS 5000
#define MAXRES 1000

int			fflush() c;
int			system() c;
int			bonded_atoms() c;         // variables used in reading in
int			sbcoil();
int			no_coil[ hashed ];

int			nswap;                    // for swapping pro-chiral pairs
int			swap_shifts();

molecule	linkprot();               //   pdb-file, etc.
int			getxyz_from_pdb();
int			getseq_from_pdb();
int			getresid_from_pdb();
float		dum;
string		matched, notmatched, seq[1];
int			numstrand, n_matched;
string		rlb, pdbfile, obsfile, rdbfile, sanderfile;
string		moltype[1];
string		strandname[1];
string		rid[ hashed ], id_res;

int			read_obs_shifts();
int			write_sander_inp();
int			obs_found;
file		temp_file, rdbf;
string		idp, key;
float		observed[ hashed ];
float		calculated[ hashed ];

                                      //  variables for ring currents:
int			get_ring_info();
string		atskip[ 26 ];
string		ringskip[ hashed ];       // protons to be skipped
int			ringsize[ MAXRINGS ];
string		ringname[ MAXRINGS ];
string		ringaname[ MAXRINGPOS ];
float		ringstr[ MAXRINGS ];
point		ringpos[ MAXRINGPOS ];
int			ringnum[ MAXRINGS ];
int			nrings, abs_res, ring, atring; // indices to get ring info
point		r1, r2;
point		rn[ MAXRINGS ];
float		b, s12, r1sq, d1, r2sq, d2;
int			nskip, skip, kp1;
float		hm[ hashed ];
float		hm_pA[ hashed ];
float		hm_pA6[ hashed ];
float		hm_pG[ hashed ];
float		hm_pG6[ hashed ];
float		hm_pC[ hashed ];
float		hm_pT[ hashed ];

atom		ae;                       // variables for electrostatics:
residue		re;
point		hx, ch;
float		rhx2, rhx, rch, helst;
int			ipep, npep, nbonds, e_res;
atom		neighbors[ 8 ];
float		hFlygare;

float		delta, rpep, rpep2;        // variables for anisotropy:
point		pc[ MAXRES ], pn[ MAXRES ];
float		pb[ MAXRES ];
int			pepres[ MAXRES ];
int			get_pep_info();

float		shhm, shhm_r, she, shp, shp_r;  // partial and total shifts
                                            //  for individual protons
float		shhm_t, she_t, shp_t;           // total shifts for equilvaent
                                            //  protons
float		shhm_p[ hashed ], she_p[ hashed ], shp_p[ hashed ],
			sconst_p[ hashed ];             // values for printing
string		fields_p[ 3 ];
float		shhm_save_D, she_save_D, shp_save_D;  // for averaging aromatics
float		shhm_save_E, she_save_E, shp_save_E;  // for averaging aromatics
int         nprot;                          // no. of protons in current
                                            //  equivalent set.
float		sconst;

int			get_sugar_info();

float		h1p[ MAXRES ], h2p1[ MAXRES ], h2p2[ MAXRES ], 
			h3p[ MAXRES ], h4p[ MAXRES ];

float		cH2[ hashed ], cH5[ hashed ], cH6[ hashed ], cH7[ hashed ],
			cH8[ hashed ];
