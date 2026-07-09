// file: select_atoms.c
// requires processed selection.l, selection.y
#include <fnmatch.h>
#include "basics.h"
#include "container.h"
#include "selection.h"
#include "elements.h"
#include "residue.h"
#include "loop.h"
#include "molecule.h" // why does _Generic() macero expander need this here only???

bool cpptraj_compatible = TRUE;

void SelFree(SELNODE node)
{
    if (!node) return;
    SelFree(node->left);
    SelFree(node->right);
    free(node->text);
    free(node->ngGrid);
    free(node);
}

VARARRAY
vaUnitEvalSelection(const SELNODE node, const UNIT uUnit)
{
    VARARRAY vaResult = vaVarArrayCreate(sizeof(ATOM));
    LOOP lAtoms = lLoop((OBJEKT) uUnit, ATOMS);
    ATOM aAtom;
    while ( (aAtom = (ATOM) oNext(&lAtoms)) ) {
        if (bAtomEvalSelection(node, aAtom)) {
            VarArrayAdd(vaResult, (GENP)&aAtom);
        }
    }
    return vaResult;
}

static void NodeEnsureGrid(const SELNODE node, const ATOM atom) {
    if (node->ngGrid) return;
    UNIT uUnit = (UNIT) cContainerWithin(cContainerWithin(atom));
    VARARRAY vaAtoms = vaUnitEvalSelection(node->left, uUnit);
    ATOM *aPAtom = PVAI(vaAtoms, ATOM, 0);
    unsigned int npoints = iVarArrayElementCount(vaAtoms);
    Point *points = malloc(npoints*sizeof(Point));
    for(unsigned int i=0; i< npoints; i++,aPAtom++) {
        points[i].x = dVX(&vAtomPosition(*aPAtom));
        points[i].y = dVY(&vAtomPosition(*aPAtom));
        points[i].z = dVZ(&vAtomPosition(*aPAtom));
        points[i].group = 0;
        points[i].p = (void*)*aPAtom;
    }
    unsigned int group_start[2]={0,npoints};
    node->ngGrid = neighbor_grid_setup(points,npoints,1, group_start,node->dist);
    VarArrayDestroy(&vaAtoms);
}

static inline bool bSelMatchText(const SELNODE node, const char *value)
{
    if (!node || !value) return false;
    /* forced_string or no glob chars: exact compare */
    if (node->forced_string || !node->has_glob)
        return strcmp(node->text, value) == 0;
    /*
     * Globbing via fnmatch. The lexer preserves the literal pattern;
     * use ~ prefix for strings that should not be globbed.
     */
    return fnmatch(node->text, value, 0) == 0;
}

static bool bNodeSelectionWithin(const SELNODE node, const ATOM atom) {
    NodeEnsureGrid(node, atom);
    return neighbor_grid_query_point_bool(node->ngGrid,
            dVX(&vAtomPosition(atom)),
            dVY(&vAtomPosition(atom)),
            dVZ(&vAtomPosition(atom)), -1);
}

static bool bSelWithinResidue(SELNODE node, ATOM atom) {
    RESIDUE rRes = (RESIDUE)cContainerWithin(atom);
    if (node->cache_object == (void*)rRes) return node->cache_result;
    node->cache_object = (void*)rRes;
    /* any atom in this residue within cutoff? */
    node->cache_result = false;
    LOOP lAtoms = lLoop((OBJEKT)rRes, ATOMS);
    ATOM aAtom;
    while ((aAtom = (ATOM)oNext(&lAtoms))) {
        if (bNodeSelectionWithin(node, aAtom)) {
            node->cache_result = true;
            break;
        }
    }
    return node->cache_result;
}

static bool bSelWithinResCenter(SELNODE node, ATOM atom) {
    RESIDUE rRes = (RESIDUE)cContainerWithin(atom);
    if (node->cache_object == (void*)rRes) return node->cache_result;
    node->cache_object = (void*)rRes;
    /* compute residue centroid */
    double cx=0, cy=0, cz=0;
    int n=0;
    LOOP lAtoms = lLoop((OBJEKT)rRes, ATOMS);
    ATOM aAtom;
    while ((aAtom = (ATOM)oNext(&lAtoms))) {
        cx += dVX(&vAtomPosition(aAtom));
        cy += dVY(&vAtomPosition(aAtom));
        cz += dVZ(&vAtomPosition(aAtom));
        n++;
    }
    if (n > 0) { cx /= n; cy /= n; cz /= n; }
    NodeEnsureGrid(node, atom);
    node->cache_result = neighbor_grid_query_point_bool(node->ngGrid, cx, cy, cz, -1);
    return node->cache_result;
}

static bool bSelWithinMolecule(SELNODE node, ATOM atom) {
    long iMol = ((RESIDUE)cContainerWithin(atom))->iMolecule;
    if ((long)node->cache_object == iMol && node->cache_object != NULL)
        return node->cache_result;
    node->cache_object = (void*)iMol;
    node->cache_result = false;
    /* ensure grid is built before querying */
    bNodeSelectionWithin(node, atom);
    /* walk all residues in unit, check those with matching molecule tag */
    UNIT uUnit = (UNIT)cContainerWithin(cContainerWithin(atom));
    LOOP lResidues = lLoop((OBJEKT)uUnit, RESIDUES);
    RESIDUE rRes;
    while ((rRes = (RESIDUE)oNext(&lResidues))) {
        if (rRes->iMolecule != iMol) continue;
        LOOP lAtoms = lLoop((OBJEKT)rRes, ATOMS);
        ATOM aAtom;
        while ((aAtom = (ATOM)oNext(&lAtoms))) {
            if (neighbor_grid_query_point_bool(node->ngGrid,
                    dVX(&vAtomPosition(aAtom)),
                    dVY(&vAtomPosition(aAtom)),
                    dVZ(&vAtomPosition(aAtom)), -1)) {
                node->cache_result = true;
                goto done;
            }
        }
    }
done:
    return node->cache_result;
}

static bool bSelMolContains(const SELNODE node, const UNIT uUnit, long iMol)
{
    LOOP lResidues = lLoop((OBJEKT)uUnit, RESIDUES);
    RESIDUE rRes;
    while ((rRes = (RESIDUE)oNext(&lResidues))) {
        if (rRes->iMolecule != iMol) continue;
        LOOP lAtoms = lLoop((OBJEKT)rRes, ATOMS);
        ATOM aAtom;
        while ((aAtom = (ATOM)oNext(&lAtoms)))
            if (bAtomEvalSelection(node->left, aAtom))
                return true;
    }
    return false;
}

#define RANGE_CHECK(v) \
    { long a=node->a, b=node->b; \
      return (a<=b) ? (v>=a && v<=b) : (v>=b && v<=a); }

bool bAtomEvalSelection(const SELNODE node, const ATOM atom)
{
    if (!node) return false;

    switch (node->kind) {
    case SEL_NODE_AND:
        return bAtomEvalSelection(node->left,  atom) &&
               bAtomEvalSelection(node->right, atom);
    case SEL_NODE_OR:
        return bAtomEvalSelection(node->left,  atom) ||
               bAtomEvalSelection(node->right, atom);
    case SEL_NODE_NOT:
        return !bAtomEvalSelection(node->left, atom);
    case SEL_NODE_ALL:
        return true;

    case SEL_NODE_RESNAME:
        return bSelMatchText(node, sContainerName(cContainerWithin(atom)));
    case SEL_NODE_CHAINID:
        return bSelMatchText(node, sResidueChainId(cContainerWithin(atom)));
    case SEL_NODE_ATOMNAME:
        return bSelMatchText(node, sContainerName(atom));
    case SEL_NODE_ATOMTYPE:
        return bSelMatchText(node, sAtomType(atom));
    case SEL_NODE_ELEMENT:
        return bSelMatchText(node, sElementName(iAtomElement(atom), NULL));

    case SEL_NODE_MOL_CONTAINS:
        {
            long iMol = ((RESIDUE)cContainerWithin(atom))->iMolecule;
            if ((long)node->cache_object == iMol && node->cache_object != NULL)
                return node->cache_result;
            node->cache_object = (void*)iMol;
            UNIT uUnit = (UNIT)cContainerWithin(cContainerWithin(atom));
            node->cache_result = bSelMolContains(node, uUnit, iMol);
            return node->cache_result;
        }

    case SEL_NODE_RES_CONTAINS:
        {
            RESIDUE rRes = ((RESIDUE)cContainerWithin(atom));
            if (node->cache_object == (void*)rRes) return node->cache_result;
            node->cache_object = (void*)rRes;
            node->cache_result = false;
            LOOP lAtoms = lLoop((OBJEKT)rRes, ATOMS);
            ATOM aAtom;
            while ( (aAtom = (ATOM) oNext(&lAtoms)) ) {
                if (bAtomEvalSelection(node->left, aAtom)) {
                    node->cache_result = true;
                    break;
                }
            }
            return node->cache_result;
        }
    case SEL_NODE_INDEX:
        return iContainerSequence(atom) == node->a;

    case SEL_NODE_RESNUM:
        if (cpptraj_compatible) goto SEL_RES_INDEX;
        /* fallthrough */
    case SEL_NODE_RES_PDBSEQ:
        return iResiduePdbSequence(cContainerWithin(atom)) == node->a;
    case SEL_NODE_RES_INDEX:
    SEL_RES_INDEX:
        return iContainerSequence(cContainerWithin(atom)) == node->a;

    case SEL_NODE_MOLNUM:
        return ((RESIDUE)cContainerWithin(atom))->iMolecule == node->a;

    case SEL_NODE_RANGE_MOL:
        RANGE_CHECK(((RESIDUE)cContainerWithin(atom))->iMolecule);

    case SEL_NODE_RANGE_ATOM:
        RANGE_CHECK(iContainerSequence(atom));

    case SEL_NODE_RANGE_RESNUM:
        if (cpptraj_compatible) goto RANGE_RES_INDEX;
        /* fallthrough */
    case SEL_NODE_RANGE_RES_PDBSEQ:
        RANGE_CHECK(iResiduePdbSequence(cContainerWithin(atom)));
    case SEL_NODE_RANGE_RES_INDEX:
    RANGE_RES_INDEX:
        RANGE_CHECK(iContainerSequence(cContainerWithin(atom)));

    case SEL_NODE_RESTYPE:
        {
             char cResType = cResidueType(cContainerWithin(atom));
             if (node->text[1]==0) return (cLower(node->text[0])==cResType);
             return !strcmp(node->text,sResidueTypeNameFromChar(cResType));
        }
    case SEL_NODE_DIST_WITHIN_ATOM:
        return bNodeSelectionWithin(node, atom);
    case SEL_NODE_DIST_BEYOND_ATOM:
        return !bNodeSelectionWithin(node, atom);

    case SEL_NODE_DIST_WITHIN_RES:
        return bSelWithinResidue(node, atom);
    case SEL_NODE_DIST_BEYOND_RES:
        return !bSelWithinResidue(node, atom);
    case SEL_NODE_DIST_WITHIN_RESCEN:
        return bSelWithinResCenter(node, atom);
    case SEL_NODE_DIST_BEYOND_RESCEN:
        return !bSelWithinResCenter(node, atom);
    case SEL_NODE_DIST_WITHIN_MOL:
        return bSelWithinMolecule(node, atom);
    case SEL_NODE_DIST_BEYOND_MOL:
        return !bSelWithinMolecule(node, atom);

    }
    // not implmenented error
    return false;
}
