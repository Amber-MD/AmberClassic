
/* ------------------------------------------------------------------ */
/*  Data structures                                                    */
/* ------------------------------------------------------------------ */

#define MAX_SYMOPS 192   /* plenty for any space group */

typedef struct {
    int rot[3][3];   /* rotation matrix, entries -1/0/1            */
    int trans[3];    /* translation in 1/12 units                  */
} SYMOPt;

typedef struct {
    int    number;
    char   xHM[24]; // max = 13
    char   old[24]; // max = 16
    char   pgrp[16]; // max is about 8
    SYMOPt symops[MAX_SYMOPS];
    int    n_symops;
} SPACEGROUPt;

extern int parse_symop(const char *s, int rot[3][3], int trans[3]);
extern int parse_spacegroup_file(int number, const char *name, SPACEGROUPt *sg_return);
extern void BuildSymopMatrices( UNIT uUnit,
                     SYMOPt *symmops, int nSymops,
                     MATRIX **maPSymops );
extern void BuildFractionalTransforms( UNIT uUnit, MATRIX M, MATRIX Mi);
