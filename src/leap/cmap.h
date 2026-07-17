#ifndef CMAP_H
#define CMAP_H

typedef char WRD[8];

// There are CMAP_TYPES number of CMAPs in a PRMTOP 
typedef struct {
    char title[256];
    WRD *reslist; // [nres]
    int nres;
    // linked list CMNT struct replaced by a single \n joined multiline string (up to 1068 chars)
    STRING comments; // accumulated comments stripped of trailing space, '\n'-joined
    int resolution; // Written to CMAP_RESOLUTION (size=CMAP_TYPES, var 'maptypes' in old code), entry number iParmIndex
    double *map; // [resolution * resolution]
                   // Written to "CMAP_PARAMETER_%02d",iParmIndex  -- includes %COMMENT comments
    WRD atmname[5];
    int residx[5]; // relative residue index typcially 1 through 3 are zero
                   // indices with 0 mark central residue and SHOULD all have the same resname,
                   // (but we do not verify this). The first entry with 0 marks the atom from
                   // which we compare the residue name to reslist.
                   // In ParmSetFindCMAP, even though only 1 resName is matched per CMAP lookup,
                   // which resName matches varies by CMAP so all 5 must be sent
} CMAPt;
typedef CMAPt *CMAP;

#endif // CMAP_H
