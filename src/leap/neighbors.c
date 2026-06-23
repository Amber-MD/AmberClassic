#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include "neighbors.h"
#include "basics.h"

typedef struct {
    uint64_t key;
    int idx;
} KeyedIdx;

typedef struct {
    unsigned int start;
    unsigned int count;
    uint64_t key;
} CellRange;


struct NeighborGrid {
    const Point *points;
    unsigned int total_points;
    int num_groups;
    const unsigned int *group_start;   // length num_groups + 1

    float r_cut;
    float cell;
    float xmin, ymin, zmin;
    float xmax, ymax, zmax;

    KeyedIdx *kidx;
    CellRange *ranges;
    unsigned int ranges_count;

    Pair *results;
    unsigned int results_count;
    unsigned int results_cap;
};

static inline uint64_t pack_cell(int32_t ci, int32_t cj, int32_t ck)
{
    const uint64_t BIAS = (1ULL << 20);
    const uint64_t MASK21 = (1ULL << 21) - 1ULL;
    uint64_t a = ((uint64_t)(ci + (int32_t)BIAS)) & MASK21;
    uint64_t b = ((uint64_t)(cj + (int32_t)BIAS)) & MASK21;
    uint64_t c = ((uint64_t)(ck + (int32_t)BIAS)) & MASK21;
    return (a << 42) | (b << 21) | c;
}

static int cmp_keyedidx(const void *pa, const void *pb)
{
    const KeyedIdx *a = (const KeyedIdx *)pa;
    const KeyedIdx *b = (const KeyedIdx *)pb;
    return (a->key > b->key) - (a->key < b->key);
}

static CellRange *build_ranges(const KeyedIdx *kidx, unsigned int n, unsigned int *out_count)
{
    if (n == 0) {
        *out_count = 0;
        return NULL;
    }

    unsigned int est = n / 4 + 4;
    CellRange *ranges;
    MALLOC(ranges, CellRange *, est * sizeof(CellRange));

    unsigned int rcount = 0;
    unsigned int i = 0;
    while (i < n) {
        uint64_t key = kidx[i].key;
        unsigned int start = i;
        unsigned int cnt = 1;
        i++;
        while (i < n && kidx[i].key == key) {
            cnt++;
            i++;
        }

        if (rcount >= est) {
            est *= 2;
            CellRange *tmp = realloc(ranges, est * sizeof(CellRange));
            if (!tmp) {
                perror("realloc ranges");
                exit(1);
            }
            ranges = tmp;
        }

        ranges[rcount].key = key;
        ranges[rcount].start = start;
        ranges[rcount].count = cnt;
        rcount++;
    }

    *out_count = rcount;
    return ranges;
}

static CellRange *find_range(CellRange *ranges, unsigned int rcount, uint64_t key)
{
    unsigned int lo = 0, hi = rcount;
    while (lo < hi) {
        unsigned int m = (lo + hi) >> 1;
        if (ranges[m].key < key) lo = m + 1;
        else hi = m;
    }
    if (lo < rcount && ranges[lo].key == key) return &ranges[lo];
    return NULL;
}

NeighborGrid *neighbor_grid_setup(const Point *points,
                                  unsigned int total_points,
                                  int num_groups,
                                  const unsigned int *group_start,
                                  float r_cut)
{
    if (!points || !group_start || num_groups <= 0 || total_points == 0) return NULL;
    if (group_start[0] != 0) {
        fprintf(stderr, "group_start[0] must be 0\n");
        return NULL;
    }
    if (group_start[num_groups] != total_points) {
        fprintf(stderr, "group_start[num_groups] must equal total_points\n");
        return NULL;
    }
    for (int g = 0; g < num_groups; ++g) {
        if (group_start[g] > group_start[g + 1]) {
            fprintf(stderr, "group_start must be nondecreasing\n");
            return NULL;
        }
    }

    float xmin = 1e30f, ymin = 1e30f, zmin = 1e30f;
    float xmax = -1e30f, ymax = -1e30f, zmax = -1e30f;
    for (unsigned int i = 0; i < total_points; ++i) {
        if (points[i].x < xmin) xmin = points[i].x;
        if (points[i].y < ymin) ymin = points[i].y;
        if (points[i].z < zmin) zmin = points[i].z;
        if (points[i].x > xmax) xmax = points[i].x;
        if (points[i].y > ymax) ymax = points[i].y;
        if (points[i].z > zmax) zmax = points[i].z;
    }

    float eps = 1e-6f * (r_cut > 1.0f ? r_cut : 1.0f);
    xmin -= eps;
    ymin -= eps;
    zmin -= eps;

    KeyedIdx *kidx;
    MALLOC(kidx, KeyedIdx *, total_points * sizeof(KeyedIdx));

    const float cell = r_cut;

#if defined(_OPENMP)
    #pragma omp parallel for schedule(static)
#endif
    for (unsigned int i = 0; i < total_points; ++i) {
        int32_t ci = (int32_t)floorf((points[i].x - xmin) / cell);
        int32_t cj = (int32_t)floorf((points[i].y - ymin) / cell);
        int32_t ck = (int32_t)floorf((points[i].z - zmin) / cell);
        kidx[i].key = pack_cell(ci, cj, ck);
        kidx[i].idx = (int)i;
    }

    qsort(kidx, total_points, sizeof(KeyedIdx), cmp_keyedidx);

    unsigned int ranges_count = 0;
    CellRange *ranges = build_ranges(kidx, total_points, &ranges_count);

    NeighborGrid *grid;
    MALLOC(grid, NeighborGrid *, sizeof(*grid));

    grid->points = points;
    grid->total_points = total_points;
    grid->num_groups = num_groups;
    grid->group_start = group_start;
    grid->r_cut = r_cut;
    grid->cell = cell;
    grid->xmin = xmin;
    grid->ymin = ymin;
    grid->zmin = zmin;
    grid->xmax = xmax;
    grid->ymax = ymax;
    grid->zmax = zmax;
    grid->kidx = kidx;
    grid->ranges = ranges;
    grid->ranges_count = ranges_count;

    grid->results_cap = 1024;
    grid->results_count = 0;
    MALLOC(grid->results, Pair *, grid->results_cap * sizeof(Pair));

    return grid;
}

static void append_result(NeighborGrid *grid, const Pair *p)
{
    if (grid->results_count == grid->results_cap) {
        unsigned int newcap = grid->results_cap * 2;
        Pair *tmp;
        REALLOC(tmp, Pair *, grid->results, newcap * sizeof(Pair));
        grid->results = tmp;
        grid->results_cap = newcap;
    }
    grid->results[grid->results_count++] = *p;
}

int neighbor_grid_query_group(NeighborGrid *grid,
                              int query_group,
                              const Pair **pairs_out,
                              unsigned int *count_out)
{
    if (!grid || !pairs_out || !count_out) return -1;
    if (query_group < 0 || query_group >= grid->num_groups) return -1;

    const Point *points = grid->points;
    const unsigned int *group_start = grid->group_start;
    const float r2 = grid->r_cut * grid->r_cut;
    const float cell = grid->cell;

    grid->results_count = 0;

    unsigned int gs = group_start[query_group];
    unsigned int ge = group_start[query_group + 1];
    if (ge > grid->total_points) ge = grid->total_points;

    for (unsigned int qi = gs; qi < ge; ++qi) {
        Point Q = points[qi];

        int32_t ci = (int32_t)floorf((Q.x - grid->xmin) / cell);
        int32_t cj = (int32_t)floorf((Q.y - grid->ymin) / cell);
        int32_t ck = (int32_t)floorf((Q.z - grid->zmin) / cell);

        for (int di = -1; di <= 1; ++di)
        for (int dj = -1; dj <= 1; ++dj)
        for (int dk = -1; dk <= 1; ++dk) {
            uint64_t key = pack_cell(ci + di, cj + dj, ck + dk);
            CellRange *cr = find_range(grid->ranges, grid->ranges_count, key);
            if (!cr) continue;

            for (unsigned int t = 0; t < cr->count; ++t) {
                int pidx = grid->kidx[cr->start + t].idx;
                if (pidx == (int)qi) continue;
                if (points[pidx].group == query_group) continue;

                float dx = points[pidx].x - Q.x;
                float dy = points[pidx].y - Q.y;
                float dz = points[pidx].z - Q.z;
                float d2 = dx*dx + dy*dy + dz*dz;
                if (d2 > r2) continue;

                Pair p;
                p.from_idx = (int)qi;
                p.from_group = points[qi].group;
                p.from_member = points[qi].member;
                p.to_idx = pidx;
                p.to_group = points[pidx].group;
                p.to_member = points[pidx].member;
                p.d2 = d2;

                append_result(grid, &p);
            }
        }
    }

    *pairs_out = grid->results;
    *count_out = grid->results_count;
    return 0;
}

int neighbor_grid_query_point(NeighborGrid *grid,
                              float x, float y, float z,
                              int query_group,
                              int query_member,
                              const Pair **pairs_out,
                              size_t *count_out)
{
    if (!grid || !pairs_out || !count_out) return -1;

    grid->results_count = 0;

    // quick reject if query is definitely too far from the model box
    float r = grid->r_cut;
    if (x < grid->xmin - r || x > grid->xmax + r ||
        y < grid->ymin - r || y > grid->ymax + r ||
        z < grid->zmin - r || z > grid->zmax + r) {
        *pairs_out = grid->results;
        *count_out = 0;
        return 0;
    }

    const float cell = grid->cell;
    const float r2 = grid->r_cut * grid->r_cut;

    int32_t ci = (int32_t)floorf((x - grid->xmin) / cell);
    int32_t cj = (int32_t)floorf((y - grid->ymin) / cell);
    int32_t ck = (int32_t)floorf((z - grid->zmin) / cell);

    for (int di = -1; di <= 1; ++di)
    for (int dj = -1; dj <= 1; ++dj)
    for (int dk = -1; dk <= 1; ++dk) {
        uint64_t key = pack_cell(ci + di, cj + dj, ck + dk);
        CellRange *cr = find_range(grid->ranges, grid->ranges_count, key);
        if (!cr) continue;

        for (size_t t = 0; t < cr->count; ++t) {
            int pidx = grid->kidx[cr->start + t].idx;
            const Point *P = &grid->points[pidx];

            if (query_group >= 0 && P->group == query_group) continue;

            float dx = P->x - x;
            float dy = P->y - y;
            float dz = P->z - z;
            float d2 = dx*dx + dy*dy + dz*dz;
            if (d2 > r2) continue;

            Pair p;
            p.from_idx = -1;          // synthetic point
            p.from_group = query_group;
            p.from_member = query_member;
            p.to_idx = pidx;
            p.to_group = P->group;
            p.to_member = P->member;
            p.d2 = d2;

            append_result(grid, &p);
        }
    }

    *pairs_out = grid->results;
    *count_out = grid->results_count;
    return 0;
}

void neighbor_grid_free(NeighborGrid *grid)
{
    if (!grid) return;
    FREE(grid->kidx);
    FREE(grid->ranges);
    FREE(grid->results);
    FREE(grid);
}

