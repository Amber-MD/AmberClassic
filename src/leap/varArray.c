/*
 *      File: varArray.c
 *
 ************************************************************************
 *                            LEAP                                      *
 *                                                                      *
 *                   Copyright (c) 1992, 1995                           *
 *           Regents of the University of California                    *
 *                     All Rights Reserved.                             *
 *                                                                      *
 *  This software provided pursuant to a license agreement containing   *
 *  restrictions on its disclosure, duplication, and use. This software *
 *  contains confidential and proprietary information, and may not be   *
 *  extracted or distributed, in whole or in part, for any purpose      *
 *  whatsoever, without the express written permission of the authors.  *
 *  This notice, and the associated author list, must be attached to    *
 *  all copies, or extracts, of this software. Any additional           *
 *  restrictions set forth in the license agreement also apply to this  *
 *  software.                                                           *
 ************************************************************************
 *                                                                      *
 *     Principal Investigator: Peter A. Kollman                         *
 *                                                                      *
 ************************************************************************
 *
 *      Description:
 *              A VARARRAY is an array that can grow in size.
 *              It is represented internally as a single block of
 *              memory that is REALLOCed periodically to allow growth.
 *
 *              VARARRAY allows the caller to:
 *                      create new VARARRAY
 *                      destroy VARARRAYs
 *                      define the size of elements
 *                      obtain pointers to the start of indexed elements
 *                      add objects to the array
 *
 *
 *      Body of varArray.[ch] modified by Vladimir Romanovski (1994).
 *
 *      SSchott (2026): Added Claude suggested modifications to allow
 *      parametrization of systems that require more than 32 bits:
 *      CRITICAL FIX FOR LARGE SYSTEMS:
 *      1. Removed hard memory limits (only check for size_t overflow)
 *      2. Fixed integer overflow in element() macro - was causing seg faults
 *         when accessing arrays with index*size > INT_MAX (~2.1GB)
 *      For systems with millions of residues/atoms, the calculation
 *      (index * element_size) can exceed INT_MAX, causing pointer
 *      arithmetic overflow. This has been fixed by casting to size_t.
 *
 *      Juno Krahn (2026):
 *      int -> size_t conversion pass, char* -> void*
 *      conversion pass (removing a 36-year-old K&R holdover). Renamed
 *      element() macro to ELEMENT(). The size_t casts no longer needed.
 *
 *      Unified sizing policy. The old SLOT_FOR_COUNT()
 *      macro was originally meant to compute growth overhead, but was
 *      too tight for VarArrayAdd's hot path, so a prior edit patched
 *      around it with separate inline growth math in VarArrayAdd
 *      instead of fixing SLOT_FOR_COUNT itself. That left two
 *      inconsistent notions of "how much padding" in the same file.
 *      SLOT_FOR_COUNT is retired. There are now exactly two sizing
 *      modes:
 *        - "grow mode": ziVarArrayDesiredSlot(count) - count plus a
 *          clamped percentage of headroom, for anything expected to
 *          keep growing (create, add, insert, delete-then-shrink,
 *          and copies of arrays that already carry headroom).
 *        - "fixed mode": slot set to exactly count, zero padding,
 *          for VarArraySetSize (an explicit "this is final" signal)
 *          and for copies of arrays that are already tight (their
 *          slot == count), on the assumption a tight source was
 *          probably finalized deliberately and copies should
 *          preserve that intent rather than invent new overhead.
 *      ziVarArrayAllocSize() DFATALs internally on overflow, and
 *      MALLOC/REALLOC DFATAL internally on allocation failure, so
 *      callers never check for NULL, and no FREE() is needed right
 *      before a DFATAL since DFATAL exits the process.
 *
 ************************************************************************/


#include        "basics.h"
#include        "varArray.h"
#include        <stdint.h>

// Pointers are all now void*, but strict C we need char* cast for pointer math
#define ELEMENT(a,i) ((void*)((unsigned char*)(a)->data + (i)*(a)->size))
#define PORTION  10

#define VARARRAY_GROWTH_NUM     3   /* grow by 50%: new = old * 3/2 */
#define VARARRAY_GROWTH_DEN     2
#define VARARRAY_MIN_GROWTH_BYTES   (PORTION * 8)      /* smallest growth chunk, in bytes */
#define VARARRAY_MIN_GROWTH_OBJECTS PORTION             /* smallest growth chunk, in elements */
#define VARARRAY_MAX_GROWTH_BYTES   (8 * 1024 * 1024)   /* largest growth chunk per grow, in bytes */

/*-----------------------------------------------------
 *      ziVarArrayAllocSize
 *
 *      Compute size*count in bytes for an allocation request.
 *      DFATALs (and therefore never returns) if the multiplication
 *      would overflow size_t. Callers can treat the return value
 *      as always valid.
 *
 *      'caller' should be the name of the calling function, for
 *      use in the fatal error message; use __func__ at call sites.
 */
static inline size_t ziVarArrayAllocSize(const char *caller, size_t size, size_t count)
{
    if (size == 0 || count == 0) return 0;

    if (count > SIZE_MAX / size) {
        DFATAL(" %s: allocation size overflow (would need %zu elements of %zu bytes)",
                caller, count, size);
    }

    return size * count;
}

/*-----------------------------------------------------
 *      ziVarArrayDesiredSlot
 *
 *      "Grow mode" sizing: given an element count and element size,
 *      return the slot count that represents a comfortable, stable
 *      capacity for that many elements - count plus a clamped
 *      percentage of headroom. The minimum growth is the larger of
 *      a byte floor (VARARRAY_MIN_GROWTH_BYTES, converted to
 *      elements via size) and a flat object-count floor
 *      (VARARRAY_MIN_GROWTH_OBJECTS), so tiny elements still get a
 *      meaningful byte cushion and huge elements still get at least
 *      a handful of extra slots rather than rounding down to 1. The
 *      maximum growth is a pure byte ceiling
 *      (VARARRAY_MAX_GROWTH_BYTES, converted to elements via size),
 *      so growth is bounded by actual memory footprint regardless
 *      of what's stored. This is the single source of truth for
 *      growth policy; every function that expects further growth
 *      (create, add, insert, and shrink-with-possible-future-growth)
 *      sizes through this, so they can never disagree about what
 *      "the right size" is for a given count.
 */
static inline size_t ziVarArrayDesiredSlot(size_t count, size_t size)
{
    size_t growth;
    size_t min_growth, max_growth;
    size_t min_growth_from_bytes;

    if (count > SIZE_MAX / VARARRAY_GROWTH_NUM) {
        /* multiply would overflow; fall straight to max clamp */
        growth = SIZE_MAX;
    } else {
        growth = count * VARARRAY_GROWTH_NUM / VARARRAY_GROWTH_DEN - count;
    }

    /* convert byte-based bounds to element counts for this element size;
       size == 0 is degenerate (nothing to bound), treat bounds as elements */
    min_growth_from_bytes = (size > 0) ? (VARARRAY_MIN_GROWTH_BYTES / size) : VARARRAY_MIN_GROWTH_BYTES;
    max_growth = (size > 0) ? (VARARRAY_MAX_GROWTH_BYTES / size) : VARARRAY_MAX_GROWTH_BYTES;

    /* minimum growth is whichever floor is larger in element terms;
       VARARRAY_MIN_GROWTH_OBJECTS is a fixed positive constant, so
       min_growth can never be 0 here regardless of size */
    min_growth = (min_growth_from_bytes > VARARRAY_MIN_GROWTH_OBJECTS)
                 ? min_growth_from_bytes : VARARRAY_MIN_GROWTH_OBJECTS;

    if (max_growth < min_growth) max_growth = min_growth;

    if (growth < min_growth) {
        growth = min_growth;
    } else if (growth > max_growth) {
        growth = max_growth;
    }

    if (count > SIZE_MAX - growth) DFATAL(" ziVarArrayDesiredSlot: slot count overflow");

    return count + growth;
}

/*-----------------------------------------------------
 *               iVarArrayPointerToIndex
 *
 *
 */
size_t iVarArrayPointerToIndex(VARARRAY header, void *data)
{
    if (header == NULL || data == NULL) DFATAL(" iVarArrayPointerToIndex: VARARRAY or Data is NULL");
    return ((data - header->data) / header->size);
}

#ifdef DEBUG // otherwise inline form header
/*-----------------------------------------------------
 *                iVarArrayElementSize
 *
 *
 */
size_t iVarArrayElementSize(VARARRAY header)
{
    if (header == NULL) DFATAL(" iVarArrayElementSize: VARARRAY is NULL");
    return header->size;
}

/*------------------------------------------------------
 *                iVarArrayElementCount
 *
 *
 */
size_t iVarArrayElementCount(VARARRAY header)
{
    if (header == NULL)
        return 0;
    return header->count;
}

/*-----------------------------------------------------
 *               PVarArrayIndex
 */
void *PVarArrayIndex(VARARRAY header, size_t pos)
{
    if (header == NULL) DFATAL(" PVarArrayIndex: VARARRAY is NULL");
    return (ELEMENT(header, pos));

}
#endif

/*-----------------------------------------------------
 *      vaVarArrayCreate
 *
 *
 *      Create a new VARARRAY and initialize it.
 *      The caller must initialize the size of the elements
 *      of the VARARRAY when they create it.
 *
 *      Grow mode: a fresh array is always expected to grow.
 */

VARARRAY vaVarArrayCreate(size_t size)
{
    VARARRAY new;
    size_t init_slot;
    size_t alloc_size;

    if (size == 0) DFATAL("vaVarArrayCreate: element size must be nonzero");

    new = (VARARRAY)MALLOC(sizeof(HeaderStruct));

    if (new == NULL) DFATAL("vaVarArrayCreate: not enough memory ");

    init_slot = ziVarArrayDesiredSlot(0, size);
    alloc_size = ziVarArrayAllocSize(__func__, size, init_slot);

    new->data = MALLOC(alloc_size);

    if (new->data == NULL) DFATAL("vaVarArrayCreate: not enough memory ");
    new->count = 0;
    new->size = size;
    new->slot = init_slot;

    MESSAGE("Created VARARRAY: element_size=%zu, initial_slots=%zu, initial_bytes=%zu",
                 size, init_slot, alloc_size);

    return new;
}

/*-----------------------------------------------------
 *      VarArrayDestroy
 *
 *      Destroy the VarArray, after this call the pointer will
 *      be undefined.
 */


void VarArrayDestroy(VARARRAY *header)
{

    if (*header == NULL) DFATAL(" VarArrayDestroy: VARARRAY is NULL");

    MESSAGE("Destroying VARARRAY: count=%zu, slots=%zu",
                 (*header)->count, (*header)->slot);

    if ((*header)->data != NULL)
        FREE((*header)->data);

    FREE(*header);
    *header = NULL;
}



/*-----------------------------------------------------
 *      VarArrayAdd
 *
 *      Add one element to the VARARRAY and copy the data into it.
 *      This will require REALLOCing the array and probably changing
 *      its address.
 *
 *      Grow mode: sizes through ziVarArrayDesiredSlot, the single
 *      shared growth policy.
 */

void VarArrayAdd(VARARRAY header, GENP data)
{
    size_t new_slot;
    size_t new_alloc_size;

    if (header == NULL) DFATAL(" VarArrayAdd: VARARRAY is NULL");

    if (header->count == header->slot) {
        new_slot = ziVarArrayDesiredSlot(header->slot, header->size);

        if (new_slot < header->slot) DFATAL(" VarArrayAdd: slot count overflow");

        new_alloc_size = ziVarArrayAllocSize(__func__, header->size, new_slot);

        MESSAGE("Growing VARARRAY: %zu -> %zu slots, %zu bytes",
                     header->slot, new_slot, new_alloc_size);

        header->data = REALLOC(header->data, new_alloc_size);

        header->slot = new_slot;
    }
    memcpy(ELEMENT(header, header->count), data, header->size);
    header->count++;
}

/*-----------------------------------------------------
 *      vaVarArrayCopy
 *        Copy the VARARRAY.
 *
 *      Sizing: if the source is already tight (slot == count), it
 *      was probably finalized deliberately (e.g. via SetSize), so
 *      the copy stays tight too. Otherwise the source still carries
 *      growth headroom, so the copy gets fresh grow-mode headroom.
 *
 */
VARARRAY vaVarArrayCopy(VARARRAY header)
{
    VARARRAY new;
    size_t alloc_size;

    if (header == NULL) DFATAL(" vaVarArrayCopy: VARARRAY is NULL");

    new = (VARARRAY)MALLOC(sizeof(HeaderStruct));
    if (new == NULL) DFATAL(" vaVarArrayCopy: cannot allocate header");

    new->size = header->size;
    new->count = header->count;

    if (header->slot == header->count) {
        new->slot = new->count;               /* source was tight: stay tight */
    } else {
        new->slot = ziVarArrayDesiredSlot(new->count, new->size);  /* source had headroom: refresh it */
    }

    alloc_size = ziVarArrayAllocSize(__func__, new->size, new->slot);

    new->data = MALLOC(alloc_size);
    if (new->data == NULL) DFATAL(" vaVarArrayCopy: cannot allocate %zu bytes", alloc_size);

    /* only header->count elements actually exist; copy that many bytes,
       not the full (possibly larger) new->slot allocation. Guard the
       zero-count case explicitly: memcpy with a NULL source pointer is
       technically UB even when the count is 0, and header->data can
       legitimately be NULL when header->count == 0. */
    if (header->count > 0) {
        memcpy(new->data, header->data,
               ziVarArrayAllocSize(__func__, header->size, header->count));
    }

    MESSAGE("Copied VARARRAY: %zu elements, %zu slots, %zu bytes",
                 new->count, new->slot, alloc_size);

    return new;
}

/*-----------------------------------------------------
 *      vaVarArrayCopy2
 *
 *      Sizing: if either source is already tight (slot == count),
 *      treat the merge as "probably finalized" and keep the result
 *      tight. Only if both sources still carry headroom does the
 *      merged result get fresh grow-mode headroom.
 */
VARARRAY vaVarArrayCopy2(VARARRAY header1, VARARRAY header2)
{
    VARARRAY new;
    size_t copysize, alloc_size;
    size_t total_count;

    if (header1 == NULL || header2 == NULL) DFATAL(" vaVarArrayCopy2: VARARRAY is NULL");
    if (header1->size != header2->size)
        DFATAL(" vaVarArrayCopy2: header sizes different\n");

    if (header1->count > SIZE_MAX - header2->count) DFATAL(" vaVarArrayCopy2: combined count overflow");
    total_count = header1->count + header2->count;

    new = (VARARRAY)MALLOC(sizeof(HeaderStruct));
    if (new == NULL) DFATAL(" vaVarArrayCopy2: cannot allocate header");

    new->size = header1->size;
    new->count = total_count;

    if (header1->slot == header1->count || header2->slot == header2->count) {
        new->slot = total_count;              /* either source tight: stay tight */
    } else {
        new->slot = ziVarArrayDesiredSlot(total_count, new->size);  /* both had headroom */
    }

    alloc_size = ziVarArrayAllocSize(__func__, new->size, new->slot);

    new->data = MALLOC(alloc_size);
    if (new->data == NULL) DFATAL(" vaVarArrayCopy2: cannot allocate %zu bytes", alloc_size);

    copysize = (size_t)new->size * header1->count;
    /* guard each copy independently: either source could be empty
       (count == 0, data possibly NULL), and memcpy with a NULL
       source is technically UB even at count 0 */
    if (header1->count > 0) {
        memcpy(new->data, header1->data, copysize);
    }
    if (header2->count > 0) {
        memcpy(new->data + copysize, header2->data, (size_t)new->size * header2->count);
    }

    MESSAGE("Merged 2 VARRAYs: %zu + %zu = %zu elements, %zu slots, %zu bytes",
                 header1->count, header2->count, total_count, new->slot, alloc_size);

    return new;
}

/*-----------------------------------------------------
 *      VarArrayInsertBeforeMore
 *
 *      Grow mode: an insert, like an add, implies the array is
 *      actively being built up.
 */
void VarArrayInsertBeforeMore(VARARRAY header, size_t pos, size_t num)
{
    size_t shift, nslot;
    size_t alloc_size;
    unsigned char *h;

    if (header == NULL)
        DFATAL(" VarArrayInsertBeforeMore: VARARRAY is NULL");

    if (pos > header->count)
        DFATAL(" VarArrayInsertBeforeMore: position=%zu", pos);

    if (header->count > SIZE_MAX - num)
        DFATAL(" VarArrayInsertBeforeMore: count overflow");

    if (header->count + num > header->slot) {
        nslot = ziVarArrayDesiredSlot(header->count + num, header->size);

        alloc_size = ziVarArrayAllocSize(__func__, header->size, nslot);

        MESSAGE("Growing for insert: %zu -> %zu slots", header->slot, nslot);

        header->data = REALLOC(header->data, alloc_size);
        header->slot = nslot;
    }

    /*
     *  update item count
     */
    header->count += num;

    /*
     *  open up insert space by shuffling remainder down
     */
    shift = header->size * num;
    h = ELEMENT(header, pos);
    memmove(h + shift, h, (header->count - num - pos) * header->size);
}

/*-----------------------------------------------------
 *      VarArrayInsertBefore
 *
 *        Add one element to the VARARRAY and move all of the data
 *        at index iPos and beyond up one element.
 *        Copy the data at data into the new element that
 *        has been opened up.
 *
 */
void VarArrayInsertBefore(VARARRAY header, size_t pos, GENP data)
{
    VarArrayInsertBeforeMore(header, pos, 1);
    memcpy(ELEMENT(header, pos), data, header->size);
}



/*-----------------------------------------------------
 *        VarArrayDelete
 *
 *        Remove an element in a VARARRAY.
 *        Move the data below the one to be removed, up one
 *
 *        Shrink policy: still grow mode. A delete is not the same
 *        signal as an explicit SetSize, so growth might resume -
 *        the shrink target is ziVarArrayDesiredSlot(count), the
 *        same "comfortable capacity" VarArrayAdd would itself pick,
 *        with an extra ~10% margin before shrinking at all so a
 *        shrink and the next append can't immediately re-trigger
 *        each other at the same boundary.
 */
void VarArrayDeleteMore(VARARRAY header, size_t pos, size_t num)
{
    size_t shift, desired;
    size_t alloc_size;
    unsigned char *h;

    if (header == NULL) {
        DFATAL(" VarArrayDelete: VARARRAY is NULL");
    }
    if (pos + num > header->count || num < 1) {
        DFATAL(" VarArrayDelete: position=%5zu num=%5zu count=%5zu", pos, num,
                header->count);
    }
    header->count -= num;

    shift = num * header->size;

    h = ELEMENT(header, pos);

    memmove(h, h + shift, (header->count - pos) * header->size);

    desired = ziVarArrayDesiredSlot(header->count, header->size);

    /* only shrink if current slot exceeds desired by more than ~10% */
    if (header->slot > desired + desired / 10) {
        alloc_size = ziVarArrayAllocSize(__func__, header->size, desired);

        MESSAGE("Shrinking VARARRAY: %zu -> %zu slots", header->slot, desired);

        header->data = REALLOC(header->data, alloc_size);
        header->slot = desired;
    }
}

/*-----------------------------------------------------
 *      VarArraySetSize
 *
 *      Change the size of the array in terms of elements.
 *      The size of the array will be adjusted so that it
 *      can contain iElements elements.
 *      All previous contents of the VARARRAY are still there,
 *      unless the VARARRAY was made smaller, then the tail is lost.
 *
 *      Fixed mode: this is the one explicit "this is final" signal
 *      in the API, so slot is set to exactly ncount - zero padding,
 *      not even rounded up. A caller that later adds more will pay
 *      for a fresh grow-mode realloc at that point, which is the
 *      right tradeoff since finalization is meant to reclaim any
 *      memory that growth headroom was holding onto.
 */
void VarArraySetSize(VARARRAY header, size_t ncount)
{
    size_t alloc_size;

    if (header == NULL) {
        DFATAL(" VarArraySetSize: VARARRAY is NULL");
    }

    header->count = ncount;

    if (ncount != header->slot) {
        alloc_size = ziVarArrayAllocSize(__func__, header->size, ncount);

        MESSAGE("Resizing VARARRAY: %zu -> %zu slots, %zu bytes",
                     header->slot, ncount, alloc_size);

        header->data = REALLOC(header->data, alloc_size);
        header->slot = ncount;
    }
}

/*-----------------------------------------------------
 *      VarArrayTrim
 *
 *      Convenience wrapper: drop any growth headroom and shrink
 *      the array to exactly its current element count. Equivalent
 *      to VarArraySetSize(header, header->count) - see that
 *      function for the fixed-mode rationale.
 */
void VarArrayTrim(VARARRAY header)
{
    if (header == NULL) DFATAL(" VarArrayTrim: VARARRAY is NULL");
    VarArraySetSize(header, header->count);
}


/*-----------------------------------------------------
 *      PVarArrayDebugIndex
 *
 *      Return a pointer to the element within the VARARRAY, but
 *      first check the bounds.  Report an  if there
 *      is an out of bound access.
 */

GENP PVarArrayDebugIndex(VARARRAY header, size_t pos, const char *file, int line)
{
    if (header == NULL) {
        DFATAL("Attempting to access an invalid VARARRAY (%s line %d).", file,
                line);
    }
    if (header->count == 0) {
        DFATAL("Attempting to access a no-data VARARRAY (%s line %d).", file,
                line);
    }
    if (pos >= header->count) {
        DFATAL("Attempted to access element: %zu in a VARARRAY of size: %zu",
                pos, header->count);
    }
    return (ELEMENT(header, pos));
}
