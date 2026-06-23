
/* ------------------------------------------------------------------ */
/*  Data structures                                                    */
/* ------------------------------------------------------------------ */

#define MAX_SYMOPS 192   /* plenty for any space group */
#define XHM_LEN    32

typedef struct {
    int rot[3][3];   /* rotation matrix, entries -1/0/1            */
    int trans[3];    /* translation in 1/12 units                  */
} SYMOPt;

typedef struct {
    int    number;
    char   xHM[XHM_LEN];
    char   old[XHM_LEN];
    SYMOPt symops[MAX_SYMOPS];
    int    n_symops;
} SPACEGROUPt;

extern int parse_symop(const char *s, int rot[3][3], int trans[3]);
extern int parse_spacegroup_file(int number, const char *name, SPACEGROUPt *sg_return);
