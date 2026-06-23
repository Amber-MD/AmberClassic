/*
 * neighbors.c
 *
 * Sparse-grid neighbor search with reusable setup and single-point/group queries.
 *
 * Overview
 * --------
 * This module builds a spatial lookup structure over a set of 3D points and
 * supports efficient radius-based neighbor queries. The grid is constructed
 * once (setup phase) and reused for subsequent queries.
 *
 * The implementation uses a sparse cell grid:
 *   - Points are assigned to grid cells based on (x,y,z) coordinates.
 *   - Cells are keyed and sorted, then compressed into contiguous ranges.
 *   - Neighbor queries examine the 27 surrounding cells.
 *
 * API Summary
 * -----------
 *   NeighborGrid *neighbor_grid_setup(points, total_points, num_groups,
 *                                     group_start, r_cut);
 *
 *   int neighbor_grid_query_group(grid, group_id,
 *                                 &pairs, &count);
 *
 *   int neighbor_grid_query_point(grid, x, y, z,
 *                                 query_group, query_member,
 *                                 &pairs, &count);
 *
 *   void neighbor_grid_free(grid);
 *
 *
 * Data Model
 * ----------
 *   Point:
 *     - x,y,z      : coordinates (float)
 *     - r          : point radius
 *     - group      : group index (must match group_start layout)
 *     - member     : user-defined identifier (often index or per-group offset)
 *                    Typically for molecules, group is a residue and member is an atom.
 *                    The atom identifer can be the index within the residue or a global
 *                    index. Radius is typically a covelent or VdW radius for type
 *                    dependent anlysis. (Only present for speed up to avoid double
 *                    indexing into the atom array with member identifier.)
 *
 *   group_start:
 *     - length = num_groups + 1
 *     - group g occupies points[group_start[g] .. group_start[g+1]-1]
 *     - group_start[num_groups] == total_points (sentinel)
 *
 *
 * Query Behavior
 * --------------
 *   neighbor_grid_query_group():
 *     - finds all neighbors within radius r_cut of all points in a group
 *     - excludes same-group matches
 *
 *   neighbor_grid_query_point():
 *     - finds neighbors of a synthetic (non-dataset) point
 *     - if query_group >= 0, excludes that group
 *     - if query_group == -1, no exclusion (returns all neighbors)
 *
 *
 * Output
 * ------
 *   Returned pairs are stored internally and exposed as:
 *
 *     const Pair *pairs;
 *     size_t count;
 *
 *   Each Pair contains:
 *     - from_idx / to_idx     : indices into original points[]
 *     - from_group / to_group
 *     - from_member / to_member
 *     - d2 (squared distance)
 *
 *   NOTE:
 *     The returned pointer is owned by the grid and is overwritten on the next
 *     query call. Consume results before calling another query.
 *
 * Performance Characteristics
 * --------------------------
 *   - Setup: O(N log N) due to sorting
 *   - Query: O(k) per point, where k is number of nearby candidates
 *   - Memory: O(N)
 *
 *   Design choices:
 *     - Sparse grid (no full 3D allocation)
 *     - Float coordinates for reduced bandwidth
 *     - OpenMP parallelization in setup phase only
 *     - Query phase is single-threaded for small group sizes
 *
 * Assumptions / Requirements
 * -------------------------
 *   - group_start[] is valid and consistent with points[]
 *   - points are grouped contiguously by group
 *   - group indices are dense: 0..num_groups-1  (All groups must exist, but size zero group is allowed)
 *   - group_start[num_groups] == total_points
 *
 *
 * Limitations
 * -----------
 *   - Binary search per cell lookup (can be replaced with hash/dense grid)
 *   - Returned results are not persistent across queries
 *   - No internal validation of input beyond basic checks
 *
 *
 * Notes
 * -----
 *   - The epsilon adjustment on xmin/ymin/zmin avoids negative cell indices
 *     due to floating point rounding.
 *   - Synthetic queries allow reuse of the grid for arbitrary probe points.
 *
 */

#ifndef NEIGHBORS_H
#define NEIGHBORS_H
#include <stdint.h>
#include <stdlib.h>
/* ---------------- Types ---------------- */

typedef struct {
    float x,y,z,r;  // r = radius for covalent or vdw radius checking
    int group;   // generally residue group
    int member;  // any member identifier (e.g. global atom index, or index wihin residue)
} Point;

typedef struct {
    int from_idx, from_group, from_member;
    int to_idx,   to_group,   to_member;
    double d2;   // distance squared
} Pair;


typedef struct NeighborGrid NeighborGrid;  // opaque to caller

/* ---------------- Public API: find_neighbors & free_neighbors ---------------- */


/*
 :neighbor_grid_query_group
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
