/* Warning:  these function declarations do not specify the arguments */
/*           this should be corrected. */

#define ERROR (-1)
#define PI   3.14159265358979323846
#define TWOPI 6.28318530717958647692
#define Rad2Deg      (180.0/PI)
#define Deg2Rad      (PI/180.0)
#define MAXAT_EDGE        30
#define MAXAT_CYCLES      10
#define MAXTOR_PROBE      10
#define MAXTOR_EDGE       20
#define MAX_FACE_EDGE     20
#define FREE_TORUS        -1
#define FREE_EDGE          1
#define BURIED_TORUS      -2
#define MAX_FACE_CYCLES    4
#define MAXTMP            20

REAL_T __dx__, __dy__, __dz__;

#define DOT(pi,pj) ( (pi)[0]*(pj)[0] + (pi)[1]*(pj)[1] + (pi)[2]*(pj)[2] )

#define DIST(pi,pj) \
        (sqrt( ( (pi)[0] - (pj)[0] )* ( (pi)[0] - (pj)[0] ) + \
               ( (pi)[1] - (pj)[1] )* ( (pi)[1] - (pj)[1] ) + \
	           ( (pi)[2] - (pj)[2] )* ( (pi)[2] - (pj)[2] )  ) )

typedef REAL_T POINT[3];

#define NAME_SIZE       8

/************* atom data structures ************/
typedef struct res {
#ifdef DEBUG
	char nam[NAME_SIZE];
#endif
	int num;
} RES;

typedef struct atom {
	POINT pos;
	REAL_T q, rad;
#ifdef DEBUG
	char anam[NAME_SIZE];
	int anum;
	int res;
#endif
	int buried;	/* 1 = buried */
	int neighbor_start, n_neighbors;
	int upper_start, n_upper;   /* points to neighbor_torus struct */
	int torus_start;	    /* points to torus struct */
    int ntorus; /* note that the toruslist[] only contains pairs of atoms 
			 in ascending order if you want all tori that contain an 
			 atom you have to look through the whole toruslist  */
	int n_convex_edges;	/*  convex edges associated with atom */
	int convex_edges[MAXAT_EDGE]; /*  convex edges associated with atom */
	int n_cycles;       /*  cycles of edges associated with atom */
	int cycle_start;    /*  points to start of cyclelist for the atom */
	REAL_T area;    /* accessible surface area  associated with the atom */
} ATOM;

typedef struct neighbor_torus {
	int iatom;
	int nprobes;	/* -1 = buried, 0 = free, +1 = partially free */
} NEIGHBOR_TORUS;

typedef struct neighbor{
	int iatom;
} NEIGHBOR;

/************* torus and probe data structures *********/
typedef struct {
        POINT center, uv;       /* torus center and axis unit vector */
        REAL_T rad;
        int a1, a2;
	int concave_edges[MAXTOR_EDGE], n_concave_edges;
	int convex_edges[MAXTOR_EDGE], n_convex_edges;
	int circle1, circle2;
	int low;		/* >1 = low, 0 = normal */
} TORUS;

typedef struct probe {
    POINT pos; 
	int a1, a2, a3;        /* atoms associated with probe */
	int c1, c2, c3;		/* circles associated with probe */
	REAL_T height;          /* height of probe above base plane */
	int low;		/* >1 = low, 0 = normal */
} PROBE; 

/************* edge data structures ************/

typedef struct vertex { 
	POINT pos; 
	int iatom;      /* atom associated with vertex */
    int iprobe;   /* probe associated with vertex */
	REAL_T beta;   /* concave triangle angle */
} VERTEX;

typedef struct edge {
    int vert1, vert2;
	int circle;
	int alive;      /* 0 = dead (not part of surface) */
} EDGE;

typedef struct circle {
    int torus;
	int atom_or_probe_num;		/* concave = probe  convex = atom */
	REAL_T rad;
	POINT center;
	POINT axis;
} CIRCLE;

/************* face data structures ************/

typedef struct concave_face {  /* each concave face has 3 edges, */
    int e1, e2, e3;	       /* edges                              */
	int probe;	       /* probe associated with concave face */
	int alive;		/* 1 = active  0 = dead */
	REAL_T area;
} CONCAVE_FACE;

typedef struct saddle_face {   /* edges oriented clockwise when view from "above" torus */
			       /* starting with concave edge */
    int e1_concave, e2_convex, e3_concave, e4_convex;
	int torus;		/* torus associated with saddle face */
	REAL_T area;
	int alive;
} SADDLE_FACE;

typedef struct cycle {
    int nedges;
	int edge[MAX_FACE_EDGE];
	int atom;
	REAL_T area;
} CYCLE;

typedef struct convex_face {
      int n_cycles;           /* 0 or more cycles border a convex face */
      int cycle[MAX_FACE_CYCLES];
      int atom;       /* atom associated with convex face */
      REAL_T area;
} CONVEX_FACE;


/************** Cusp trimming data structures *********************/

typedef struct low_torus {	/* torus that intersects itself */
      int itorus;		/* index in  TORUS[] array */
      int vert1, vert2;/* 2 vertices where the torus intersects itself */
      int nfaces;
      int face[MAXTOR_PROBE];  /* broken_concave_faces associated with low_torus */
      int ncones;             
      int cone[MAXTOR_PROBE];  /* cone faces associated with low torus*/
} LOW_TORUS;

typedef struct cone_face {	/* cones are remaining part of saddles from
				   a self-intersecting torus */
	int e1_convex;		/* only one convex edge	*/
	int e2_concave, e3_concave;	/* two edges joining convex edge to cusp point */
	int itorus;		/* index in TORUS array 	*/
	int cusp_vertex;	/* vertex that is cone apex */
	REAL_T area;
} CONE_FACE;

typedef struct broken_concave_face {	
	int itorus[3];			
	int probe;
	int n_cycles;			
	int concave_cycle[MAX_FACE_CYCLES];
	int alive;
	REAL_T area;
} BROKEN_CONCAVE_FACE;

/* concave face resulting from CONCAVE_FACE intersections  */
/* it's possible for face to be split up */
/* edges that form a cycle */

typedef struct concave_cycle {
    int nedges;
	int edge[MAX_FACE_EDGE];
	int edge_direction[MAX_FACE_EDGE];	/* 1 right hand rule; -1 left hand */
	int cusp_edge[MAX_FACE_EDGE];		/* cusp edge index, -1 = ordinary edge */
	int iprobe;
	int iface;
	int intersects_self;			/* 1 = 2 cusps intersect 0 = no intersecting cusps */
	REAL_T area;
} CONCAVE_CYCLE;

typedef struct cusp_edge {
	int cycle1;		/* first  cycle */
	int cycle2;		/* second cycle */
	int edge;				/* concave edge index */
	int probe1, probe2;			/* probes that form cusp edge */
	int alive;				/* 1 = alive; 0 = dead */
	int concentric_pair;			/* 1 = intersects w/ a conc. cusp -> */
	                                        /* can't intersect w/ anyone else    */
}CUSP_EDGE;

typedef struct cusp_group {
	int n_pairs;			/* number of cusps in group */
	int cusp_pair[MAX_FACE_EDGE];	/* new_cusps in group */
} CUSP_GROUP;

/* new_cusp is for storing non_axial cusps, which may or may not
   become real cusps, depending on the non_axial cusps that it intersects */

typedef struct cusp_pair {
	POINT circle_center, circle_axis, vert1, vert2;
	REAL_T circle_rad;
	int cycle1, cycle2, cycle3;
	int cusp1, cusp2;		/*  cusps that intersect at new_cusp */
	int group;			/* group to which new_cusp belogs */
} CUSP_PAIR;

typedef struct extreme_vertex {
  int cusp_pair;		/* index in cusp_pair[] array */
  int vert;			/* 1 = cusp_pair[].vert1  2 = cusp_pair[].vert2 */
  int vert_index;		/* index in vertex[] array */
} EXTREME_VERTEX;

#if 0
/*  routine declarations   */

static void check_broken_faces ();
static void add_edge ();
static void add_free_edge ();
static void add_probe ();
static void add_saddle_face ();
static void cross ();
static void memory_usage ();
static void vnorm ();
static void write_verts ();
static void atom_vertex_match ();
static int getneighbors ();
static int get_probes ();
static int probe_pos ();
static int no_bump ();
static int t_buried ();
static int bury_check ();
static int is_buried ();
static int get_torus ();
static void torus_data ();
static void sort_neighbors ();
static int convex_circles ();
static int concave_circles ();
static void concave_edges ();
static void addvert ();
static void convex_edges ();
static void check_data ();
static void face_info ();
static void check_convex_edges ();
static void add_convex_edge_to_torus ();
static void add_convex_edge_to_atom ();
static void convex_faces ();
static void sort_edges ();
static void cycles ();
static int new_edge ();
static int next_edge ();
static void dump_atom_edges ();
static void draw_arc ();
static void draw_circle ();
static void draw_edges ();
static int is_cycle_inside ();
static REAL_T get_angle ();
static REAL_T convex_area ();
static REAL_T cycle_piece ();
static void write_info ();
static int concave_area ();
static int saddle_area ();
static int id_torus ();
static void make_cones ();
static int one_sided_torus ();
static void add_2_verts ();
static int get_low_torus_index ();
static void make_2_cone_faces ();
static void kill_saddle_face ();
static void cone_init ();
static void write_cone_info ();
static void make_broken_faces ();
static void add_concave_cycle ();
static void broken_face_info ();
static REAL_T get_cone_area ();
static void axial_trim ();
static int get_cycle_edge ();
static void sort_faces ();
static void add_circle ();
static int cone_edge ();
static void add_edges_2_cycle ();
static void check_cycle ();
static void concentric_axial_cusps ();
static void reroute ();
static void split_cycle ();
static int next_cycle_edge ();
static void non_axial_trim ();
static int new_cusp_in_group ();
static void cusp_intersect ();
static int get_cycle_id ();
static void get_probe_id ();
static void add_new_cusp ();
static int center_cycle ();
static int number_of_cusps ();
static int new_cusp ();
static void make_circle ();
static void trim_2_cusps ();
static void unique_cycles ();
static void trim_3_cusps ();
static void split_3_cusps ();
static void trim_4_cusps ();
static void get_xvertex ();
static REAL_T dist2 ();
static void copy_vert ();
static void get_3_xvertex ();
static void copyvec ();
static int get_cusps ();
static int span_pairs ();
static void dump_cycle ();
static void add_2_cusps ();
static void add_2cusp_verts ();
static void split_old_cusps ();
static void make_new_cusp ();
static void get_faces ();
static int is_new_face ();
static void split_face ();
static void get_starters ();
static void copy_cycle ();
static int cusp_match ();
static int normal_match ();
static void add_non_axial_cusp ();
static void add_1_vert ();
static void add_cusp_verts ();
static void add_cusp_circle ();
static void get_2_cycles ();
static void add_1_cusp ();
static void cusp_intersect_2 ();
static int broken_concave_area ();
static REAL_T conc_cycle_piece ();
static REAL_T interior_angle ();
static void allocate_memory ();
static void free_memory ();
#endif
