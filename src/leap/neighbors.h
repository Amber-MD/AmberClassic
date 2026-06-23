#ifndef NEIGHBORS_H
#define NEIGHBORS_H
#include <stdint.h>
#include <stdlib.h>
/* ---------------- Types ---------------- */

_Static_assert(sizeof(int)==sizeof(float), "Float and int sizes must match");
typedef struct {
    float x,y,z;
    int group;        // residue number, can be groups by whole molecule
    union {
        int member;   // Atom numer
        float r;      // or covalent radius
    };
} Point;

typedef struct {
    int from_idx, from_group;
    union { int from_member; float from_r; };
    int to_idx,   to_group;
    union { int to_member; float to_r; };
    double d2;
} Pair;


typedef struct NeighborGrid NeighborGrid;  // opaque to caller

/* ---------------- Public API: find_neighbors & free_neighbors ---------------- */

/*
 find_neighbors:
 - points: array of Point[total_points]
 - total_points: number of points
 - num_groups: number of groups (groups numbered 0..num_groups-1)
 - group_start: array length num_groups, start index in points[] for each group
 - r_cut: neighbor radius
 Never returns on failure;
*/
NeighborGrid *neighbor_grid_setup(const Point *points,
                                  unsigned int total_points,
                                  int num_groups,
                                  const unsigned int *group_start,
                                  float r_cut);

int neighbor_grid_query_group(NeighborGrid *grid,
                              int query_group,
                              const Pair **pairs_out,
                              unsigned int *count_out);

int neighbor_grid_query_point(NeighborGrid *grid,
                              float x, float y, float z,
                              int query_group,
                              int query_member,
                              const Pair **pairs_out,
                              size_t *count_out);

void neighbor_grid_free(NeighborGrid *grid);

#endif //NEIGHBORS_H
