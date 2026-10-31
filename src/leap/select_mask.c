// file: select_mask.c
// requires processed selection.l, selection.y
#include <fnmatch.h>
#include "basics.h"
#include "container.h"
#include "select_mask.h"
#include "elements.h"
#include "residue.h"
#include "molecule.h"
#include "internal.h" // only due to loop.h _Generic() dispatch
#include "loop.h"
#include "unitio.h"
#include "defaults.h"

/* AtomMaskPrepare: stamp sequential TempInt on residues and atoms in
 * output order (respecting bReorderResidues flag). Must be called before
 * any selection evaluation on a unit. Always re-runs since LEAP has no
 * per-atom dirty tracking -- cost is O(N) traversal which is acceptable
 * even for large systems. Molecule numbers are also refreshed here via
 * iUnitLabelMolecules().
 *
 * Final index numbering must order atoms exactly as unitio.c routines:
 *   UnitDoAtom()
 *   UnitIOBuildTables()
 *   iUnitIOAmberOrderResiduesUnitIOBuildTables()
 *
 * The actual sort happens with DIRECTCONTENTSPARMORDER or
 * DIRECTCONTENTSBYSEQNUM (old LEAP used convoluted ad-hoc ordering chaos)
 */
void
AtomMaskPrepare(UNIT uUnit)
{
    /* marks molecule number in RESIDUE->iMolecule */
    int iResidueCount = iUnitLabelMolecules(uUnit);
    if (iResidueCount == 0) return;

    if ( !GDefaults.bReorderResidues )
        VP2("\"order_residues\" off: keep input residue order.\n");

    int iResIndex = 0, iAtomIndex = 0;
    OBJEKT oObj;
    LOOPOVERALL(uUnit, GDefaults.bReorderResidues ?
            DIRECTCONTENTSPARMORDER : DIRECTCONTENTSBYSEQNUM,
            oObj, lResidues) {
        if (iObjectType(oObj) != RESIDUEid) continue;
        RESIDUE rRes = RESIDUE_from(oObj);
        ContainerSetTempInt(rRes, ++iResIndex);
        int iAtomResIndex = 0;
        ATOM aAtom;
        LOOPOVERALL(rRes, GDefaults.bReorderResidues ?
                DIRECTCONTENTSPARMORDER : DIRECTCONTENTSBYSEQNUM,
                aAtom, lAtoms) {
            ContainerSetTempInt(aAtom, ++iAtomIndex);
            AtomSetIndex(aAtom,++iAtomResIndex);
        }
    }
}

LIST lAtomMaskSelect(UNIT uUnit, char *sMask)
{
    SELNODE selNode = selParseAtomMask(sMask);
    AtomMaskPrepare(uUnit);
    LIST lAtoms = lUnitEvalSelection(selNode, uUnit);
    SelFree(selNode);
    return lAtoms;
}

VARARRAY vaAtomMaskSelect(UNIT uUnit, char *sMask)
{
    SELNODE selNode = selParseAtomMask(sMask);
    AtomMaskPrepare(uUnit);
    VARARRAY vaAtoms = vaUnitEvalSelection(selNode, uUnit);
    SelFree(selNode);
    return vaAtoms;
}

void SelFree(SELNODE node)
{
    if (!node) return;
    SelFree(node->left);
    SelFree(node->right);
    free(node->text);
    neighbor_grid_free(node->ngGrid);
    free(node->pPGridPoints);
    free(node);
}

VARARRAY
vaUnitEvalSelection(const SELNODE node, const UNIT uUnit)
{
    VARARRAY vaResult = vaVarArrayCreate(sizeof(ATOM));
    ATOM aAtom;
    LOOPOVERALL(uUnit, ATOMS, aAtom, lAtoms) {
        if (bAtomEvalSelection(node, aAtom))
            VarArrayAdd(vaResult, (GENP)&aAtom);
    }
    return vaResult;
}

LIST
lUnitEvalSelection(const SELNODE node, const UNIT uUnit)
{
    LIST lResult = LIST_from(oCreate(LISTid));
    ATOM aAtom;
    LOOPOVERALL(uUnit, ATOMS, aAtom, lAtoms) {
        if (bAtomEvalSelection(node, aAtom))
            ListAddToEnd(lResult, OBJEKT_from(aAtom));
    }
    return lResult;
}

/* Build the nonbond grid from the reference selection (node->left).
 * Called lazily on first distance query; cached on node for reuse.
 * Points array is also cached on node (pPGridPoints) and freed in SelFree.
 * Grid must be invalidated via SelInvalidateCoords() if coordinates change. */
static void NodeEnsureGrid(const SELNODE node, const ATOM atom)
{
    if (node->ngGrid) return;
    UNIT uUnit = UNIT_from(cContainerWithin(cContainerWithin(atom)));
    VARARRAY vaAtoms = vaUnitEvalSelection(node->left, uUnit);
    ATOM *aPAtom = PVAI(vaAtoms, ATOM, 0);
    unsigned int npoints = iVarArrayElementCount(vaAtoms);
    Point *points = node->pPGridPoints = malloc(npoints * sizeof(Point));
    for (unsigned int i = 0; i < npoints; i++, aPAtom++) {
        points[i].x     = dVX(&vAtomPosition(*aPAtom));
        points[i].y     = dVY(&vAtomPosition(*aPAtom));
        points[i].z     = dVZ(&vAtomPosition(*aPAtom));
        points[i].group = 0;
        points[i].p     = (void*)*aPAtom;
    }
    /* single group -- group_start has n_groups+1 elements, last is sentinel */
    unsigned int group_start[2] = {0, npoints};
    node->ngGrid = neighbor_grid_setup(points, npoints, 1, group_start, node->dist);
    VarArrayDestroy(&vaAtoms);
}

// Invalidate cached grids after coordinate change. Walk the full AST.
// Call between trajectory frames or after any AtomSetPosition() calls.
// NOTE: Only needed for a persistent atomMask AST parse, but currently all are transient.
void SelInvalidateCoords(SELNODE node)
{
    if (!node) return;
    if (node->ngGrid) {
        neighbor_grid_free(node->ngGrid);
        node->ngGrid = NULL;
        free(node->pPGridPoints);
        node->pPGridPoints = NULL;
        node->cache_object = SEL_CACHE_INVALID;
    }
    SelInvalidateCoords(node->left);
    SelInvalidateCoords(node->right);
}

static inline bool bSelMatchText(const SELNODE node, const char *value)
{
    if (!node || !value) return false;
    /* ~ forced_string or no glob chars: exact compare */
    if (node->forced_string || !node->has_glob)
        return strcmp(node->text, value) == 0;
    /* globbing via fnmatch() */
    return fnmatch(node->text, value, 0) == 0;
}

static bool bNodeSelectionWithin(const SELNODE node, const ATOM atom)
{
    NodeEnsureGrid(node, atom);
    return neighbor_grid_query_point_bool(node->ngGrid,
            dVX(&vAtomPosition(atom)),
            dVY(&vAtomPosition(atom)),
            dVZ(&vAtomPosition(atom)), -1);
}

static bool bSelWithinResidue(SELNODE node, ATOM atom)
{
    RESIDUE rRes = RESIDUE_from(cContainerWithin(atom));
    if (node->cache_object == (void*)rRes) return node->cache_result;
    node->cache_object = (void*)rRes;
    node->cache_result = false;
    ATOM aAtom;
    FOR_ATOMS_IN_RESIDUE(aAtom, rRes) {
        if (bNodeSelectionWithin(node, aAtom)) {
            node->cache_result = true;
            break;
        }
    }
    return node->cache_result;
}

static bool bSelWithinResCenter(SELNODE node, ATOM atom)
{
    RESIDUE rRes = RESIDUE_from(cContainerWithin(atom));
    if (node->cache_object == (void*)rRes) return node->cache_result;
    node->cache_object = (void*)rRes;
    double cx = 0, cy = 0, cz = 0;
    int n = 0;
    ATOM aAtom;
    FOR_ATOMS_IN_RESIDUE(aAtom, rRes) {
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

static bool bSelWithinMolecule(SELNODE node, ATOM atom)
{
    long iMol = RESIDUE_from(cContainerWithin(atom))->iMolecule;
    if ((long)node->cache_object == iMol && node->cache_object != NULL)
        return node->cache_result;
    node->cache_object = (void*)iMol;
    node->cache_result = false;
    NodeEnsureGrid(node, atom);
    UNIT uUnit = UNIT_from(cContainerWithin(cContainerWithin(atom)));
    LOOP lResidues = lLoop(OBJEKT_from(uUnit), RESIDUES);
    RESIDUE rRes;
    while ( (rRes = RESIDUE_from(oNext(&lResidues))) ) {
        if (rRes->iMolecule != iMol) continue;
        ATOM aAtom;
        FOR_ATOMS_IN_RESIDUE(aAtom, rRes) {
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
    LOOP lResidues = lLoop(OBJEKT_from(uUnit), RESIDUES);
    RESIDUE rRes;
    while ( (rRes = RESIDUE_from(oNext(&lResidues))) ) {
        if (rRes->iMolecule != iMol) continue;
        ATOM aAtom;
        FOR_ATOMS_IN_RESIDUE(aAtom, rRes) {
            if (bAtomEvalSelection(node->left, aAtom))
                return true;
        }
    }
    return false;
}

static bool bCompare(double a, int op, double b) {
   switch(op) {
   case OP_LT: return a < b;
   case OP_LE: return a <= b;
   case OP_EQ: return a == b;
   case OP_NE: return a != b;
   case OP_GT: return a > b;
   case OP_GE: return a >= b;
   default: return false;
   }
}

#define RANGE_CHECK(v) \
    { long _a = node->a, _b = node->b; \
      return (_a<=_b) ? ((v)>=_a && (v)<=_b) : ((v)>=_b && (v)<=_a); }

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
    case SEL_NODE_ELEMENT_NUM:
        return iAtomElement(atom) == node->a;
    case SEL_NODE_RANGE_ELEM:
        RANGE_CHECK(iAtomElement(atom));

    case SEL_NODE_MOL_CONTAINS:
        {
            long iMol = RESIDUE_from(cContainerWithin(atom))->iMolecule;
            if ((long)node->cache_object == iMol && node->cache_object != NULL)
                return node->cache_result;
            node->cache_object = (void*)iMol;
            UNIT uUnit = UNIT_from(cContainerWithin(cContainerWithin(atom)));
            node->cache_result = bSelMolContains(node, uUnit, iMol);
            return node->cache_result;
        }

    case SEL_NODE_RES_CONTAINS:
        {
            RESIDUE rRes = RESIDUE_from(cContainerWithin(atom));
            if (node->cache_object == (void*)rRes) return node->cache_result;
            node->cache_object = (void*)rRes;
            node->cache_result = false;
            ATOM aAtom;
            FOR_ATOMS_IN_RESIDUE(aAtom, rRes) {
                if (bAtomEvalSelection(node->left, aAtom)) {
                    node->cache_result = true;
                    break;
                }
            }
            return node->cache_result;
        }

    case SEL_NODE_INDEX:       // @<number>
        return iContainerTempInt(atom) == node->a; // PDB global index
    case SEL_NODE_ATOM_RESIDX: // @;
        return iAtomIndex(atom) == node->a; // PDB index within residue
    case SEL_NODE_ATOM_INDEX:  // @#
        return iContainerSequence(atom) == node->a; // LEAP index within residue
    case SEL_NODE_RANGE_ATOM:
        RANGE_CHECK(iContainerTempInt(atom));
    case SEL_NODE_RANGE_ATOM_RESIDX:
        RANGE_CHECK(iAtomIndex(atom));
    case SEL_NODE_RANGE_ATOM_INDEX:
        RANGE_CHECK(iContainerSequence(atom));

    case SEL_NODE_RESNUM:      // :<number>
        return iContainerTempInt(cContainerWithin(atom)) == node->a;
    case SEL_NODE_RES_PDBSEQ:  // :;
        return iResiduePdbSequence(cContainerWithin(atom)) == node->a;
    case SEL_NODE_RES_INDEX:   // :#
        return iContainerSequence(cContainerWithin(atom)) == node->a;
    case SEL_NODE_RANGE_RESNUM:
        RANGE_CHECK(iContainerTempInt(cContainerWithin(atom)));
    case SEL_NODE_RANGE_RES_PDBSEQ:
        RANGE_CHECK(iResiduePdbSequence(cContainerWithin(atom)));
    case SEL_NODE_RANGE_RES_INDEX:
        RANGE_CHECK(iContainerTempInt(cContainerWithin(atom)));

    case SEL_NODE_MOLNUM:
        return RESIDUE_from(cContainerWithin(atom))->iMolecule == node->a;
    case SEL_NODE_RANGE_MOL:
        RANGE_CHECK(RESIDUE_from(cContainerWithin(atom))->iMolecule);

    case SEL_NODE_RESTYPE:
        {
            char cResType = cResidueType(cContainerWithin(atom));
            if (node->text[1] == '\0') return (cLower(node->text[0]) == cResType);
            return !strcmp(node->text, sResidueTypeNameFromChar(cResType));
        }
    case SEL_NODE_RES_FLAG:
        return (fResidueFlags(cContainerWithin(atom)) & node->a) != 0;
    case SEL_NODE_RES_FLAG_NOT:
        return (fResidueFlags(cContainerWithin(atom)) & node->a) == 0;
    case SEL_NODE_RES_ATTR_CMP:
        {
            double dValue, dTemp;
            RESIDUE rRes = RESIDUE_from(cContainerWithin(atom));
            if (node->cache_object == rRes) return node->cache_result;
            switch(node->a) {
            case RES_CHARGE: ContainerTotalCharge((CONTAINER)cContainerWithin(atom), &dValue, &dTemp );
            default: return false; // error
            }
            bool result = bCompare(dValue, node->b, node->dist);
            node->cache_object = rRes;
            node->cache_result = result;
            return result;
        }
    case SEL_NODE_ATOM_FLAG:
        return (fAtomFlags(atom) & node->a) != 0;
    case SEL_NODE_ATOM_FLAG_NOT:
        return (fAtomFlags(atom) & node->a) == 0;
    case SEL_NODE_ATOM_ATTR_CMP_INT:
        return bCompare( *(const int *)((const char *)atom + node->a), node->b, node->dist);
    case SEL_NODE_ATOM_ATTR_CMP_REAL:
        return bCompare( *(const double *)((const char *)atom + node->a), node->b, node->dist);

    case SEL_NODE_DIST_WITHIN_ATOM:   return  bNodeSelectionWithin(node, atom);
    case SEL_NODE_DIST_BEYOND_ATOM:   return !bNodeSelectionWithin(node, atom);
    case SEL_NODE_DIST_WITHIN_RES:    return  bSelWithinResidue(node, atom);
    case SEL_NODE_DIST_BEYOND_RES:    return !bSelWithinResidue(node, atom);
    case SEL_NODE_DIST_WITHIN_RESCEN: return  bSelWithinResCenter(node, atom);
    case SEL_NODE_DIST_BEYOND_RESCEN: return !bSelWithinResCenter(node, atom);
    case SEL_NODE_DIST_WITHIN_MOL:    return  bSelWithinMolecule(node, atom);
    case SEL_NODE_DIST_BEYOND_MOL:    return !bSelWithinMolecule(node, atom);
    }
    VPFATAL("bAtomEvalSelection: unhandled node kind %d\n", node->kind);
    return false;
}

FLAGS fAtomFlagNameToBit(const char *name) {
    struct { const char *name; FLAGS flag; } flag_names[12] = {
        { "posfxd", ATOMPOSITIONFIXED },
        { "posblt", ATOMPOSITIONBUILT },
        { "posdrwn", ATOMPOSITIONDRAWN },
        { "selected", ATOMSELECTED },
        { "pert", ATOMPERTURB },
        { "notdisp", ATOMNOTDISPLAYED },
        { "touched", ATOMTOUCHED },
        { "posknwn", ATOMPOSITIONKNOWN },
        { "internal", ATOMPOSITIONINTERNAL },
        { "needsmin", ATOMNEEDSMINIMIZER },
        { "needsbuild", ATOMNEEDSBUILD },
        { NULL, 0 }
    };
    for (int i=0; flag_names[i].name; i++) {
        if (!strcasecmp(flag_names[i].name,name)) {
            return flag_names[i].flag;
        }
    }
    VPWARN("AtomMask: '%s' is not a valid atom flag name. Known flags:\n",name);
    for (int i=0; flag_names[i].name; i++)
        VP0("   %s\n", flag_names[i].name);
    return 0x00008000; // return a bit that cannot be matched
}

FLAGS fResFlagNameToBit(const char *name) {
    struct { const char *name; FLAGS flag; } flag_names[7] = {
        { "unknown",RESIDUEUNKNOWN },
        { "bulksolv",RESIDUEBULKSOLVENT },
        { "incap", RESIDUEINCAP },
        { "noend",RESIDUENOEND },
        { "firstend",RESIDUEFIRSTEND },
        { "lastend",RESIDUELASTEND },
        { NULL, 0 }
    };
    for (int i=0; flag_names[i].name; i++) {
        if (!strcasecmp(flag_names[i].name,name)) {
            return flag_names[i].flag;
        }
    }
    VPWARN("AtomMask: '%s' is not a valid residue flag name. Known flags:\n",name);
    for (int i=0; flag_names[i].name; i++)
        VP0("   %s\n", flag_names[i].name);
    return 0x00008000; // return a bit that cannot be matched
}

bool bAtomAttrNameToOffset(const char *name, size_t *offset, bool *is_real) {
    struct { const char *name; size_t offset; bool is_real; } attr_names[] = {
        { "element", offsetof(ATOMt, iAtomicNumber), false },
        { "atomserial", offsetof(ATOMt, iAtomSerial), false },
        { "bfactor", offsetof(ATOMt, dBfactor), true },
        { "occupancy", offsetof(ATOMt, dOccupancy), true },
        { "charge", offsetof(ATOMt, dCharge), true },
        { "pertcharge", offsetof(ATOMt, dPertCharge), true },
        { "position.x", offsetof(ATOMt, vPosition)+offsetof(VECTOR,dX), true },
        { "position.y", offsetof(ATOMt, vPosition)+offsetof(VECTOR,dY), true },
        { "position.z", offsetof(ATOMt, vPosition)+offsetof(VECTOR,dZ), true }
    }; 
    for (int i=0; attr_names[i].name; i++) {
        if (!strcasecmp(attr_names[i].name,name)) {
            *offset = attr_names[i].offset;
            *is_real = attr_names[i].is_real;
            return true;
        }
    }
    VPWARN("AtomMask: '%s' is not a valid atom attribute. Supported attributes:\n",name);
    for (int i=0; attr_names[i].name; i++)
        VP0("   %s\n", attr_names[i].name);
    *is_real=false;
    return false;
}

int iResAttrNameToEnum(const char *name) {
    const char *attr_names[] = { "charge", NULL };
    for (int i=0; attr_names[i]; i++) {
        if (!strcasecmp(attr_names[i],name)) return i+1;
    }
    VPWARN("AtomMask: '%s' is not a valid residue attribute. Supported attributes:\n",name);
    for (int i=0; attr_names[i]; i++)
        VP0("   %s\n", attr_names[i]);
    return 0;
}
