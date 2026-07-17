#ifndef SELECT_MASK_H
#define SELECT_MASK_H
#include <stddef.h>
#include "basics.h"
#include "unit.h"
#include "atom.h"
#include "varArray.h"
#include "neighbors.h"

typedef enum {
    SEL_NODE_AND,
    SEL_NODE_OR,
    SEL_NODE_NOT,
    SEL_NODE_ALL,

    SEL_NODE_RESNUM,
    SEL_NODE_RES_PDBSEQ,
    SEL_NODE_RES_INDEX,
    SEL_NODE_RESNAME,
    SEL_NODE_RESTYPE,
    SEL_NODE_CHAINID,
    SEL_NODE_ATOMNAME,
    SEL_NODE_ATOMTYPE,
    SEL_NODE_ATOM_INDEX,
    SEL_NODE_ATOM_RESIDX,
    SEL_NODE_ELEMENT,
    SEL_NODE_ELEMENT_NUM,
    SEL_NODE_INDEX,
    SEL_NODE_RANGE_ATOM,
    SEL_NODE_RANGE_ATOM_INDEX,
    SEL_NODE_RANGE_ATOM_RESIDX,
    SEL_NODE_RANGE_MOL,
    SEL_NODE_RANGE_RES_INDEX,
    SEL_NODE_RANGE_RES_PDBSEQ,
    SEL_NODE_RANGE_RESNUM,
    SEL_NODE_RANGE_ELEM,
    SEL_NODE_MOLNUM,
    SEL_NODE_RES_CONTAINS,
    SEL_NODE_MOL_CONTAINS,
    SEL_NODE_DIST_WITHIN_ATOM,
    SEL_NODE_DIST_WITHIN_RES,
    SEL_NODE_DIST_WITHIN_RESCEN,
    SEL_NODE_DIST_WITHIN_MOL,
    SEL_NODE_DIST_BEYOND_ATOM,
    SEL_NODE_DIST_BEYOND_RES,
    SEL_NODE_DIST_BEYOND_RESCEN,
    SEL_NODE_DIST_BEYOND_MOL,
} SELNODEKINDt;

typedef struct SELNODEt SELNODEt;
typedef SELNODEt *SELNODE;

  /* invalid sentinel, distinct from NULL */
#define SEL_CACHE_INVALID ((void*)(intptr_t)(-1))

struct SELNODEt {
    SELNODEKINDt kind;
    struct SELNODEt *left;   /* left child (binary), single child (unary), or ref selection (distance) */
    struct SELNODEt *right;  /* right child (binary), NULL otherwise */
    char *text;    /* for name nodes */
    long a, b;     /* for numeric, numeric-range nodes */
    double dist;  /* distance, distance nodes only */
    NeighborGrid *ngGrid; /* nonbond grid cache, distance nodes only */
    Point *pPGridPoints; /* point array for neighbor grid */
    void *cache_object;
    bool cache_result;
    bool forced_string;
    bool has_glob;
};

extern SELNODE selParseAtomMask(const char *mask);
extern void SelFree(SELNODE node);
extern LIST lAtomMaskSelect(UNIT uUnit, char *sMask );
extern VARARRAY vaAtomMaskSelect(UNIT uUnit, char *sMask );
extern LIST lUnitEvalSelection(const SELNODE node, const UNIT uUnit);
extern VARARRAY vaUnitEvalSelection(const SELNODE node, const UNIT uUnit);
extern bool bAtomEvalSelection(const SELNODE node, const ATOM atom);

#endif // SELECT_MASK_H
