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
 *              The LINKEDLIST allows the caller to:
 *                      create new VARARRAY 
 *                      destroy VARARRAYs 
 *                      define the size of elements 
 *                      obtain pointers to the start of indexed elements
 *                      add objects to the . 
 *
 *
 *       Body of varArray.[ch] modified by Vladimir Romanovski (1994).  
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
 ************************************************************************/



#include        "basics.h"
#include        "varArray.h"
#include        <stdint.h>
#include        <limits.h>

#define element(a,i) ((char*)((a)->data + ((size_t)(i))*((size_t)(a)->size)))
#define PORTION  10

#define SLOT_FOR_COUNT(count) \
  ( count > 0 ? (((count + PORTION -1) / PORTION ) * PORTION) : PORTION)

/* NO hard memory limit - only check for overflow */
#define VARARRAY_DEBUG_MEMORY 0

#if VARARRAY_DEBUG_MEMORY
#define VARARRAY_LOG(fmt, ...) \
    fprintf(stderr, "VARARRAY: " fmt "\n", ##__VA_ARGS__)
#else
#define VARARRAY_LOG(fmt, ...)
#endif

/* Helper function to check if multiplication would overflow */
static inline int check_size_overflow(size_t size, int count, size_t *result)
{
    if (count < 0) return 1;  /* Invalid count */
    if (size == 0) {
        *result = 0;
        return 0;
    }
    
    size_t c = (size_t)count;
    
    /* Check if multiplication would overflow */
    if (c > SIZE_MAX / size) {
        return 1;  /* Would overflow */
    }
    
    *result = size * c;
    return 0;  /* OK */
}

/*-----------------------------------------------------
 *               iVarArrayPointerToIndex
 *
 *
 */
int iVarArrayPointerToIndex(VARARRAY header, char *data)
{
    if (header == NULL || data == NULL) {
        DFATAL((" iVarArrayPointerToIndex: VARARRAY or Data is NULL"));
    }
    return ((data - header->data) / header->size);
}

/*-----------------------------------------------------
 *                iVarArrayElementSize
 *
 *
 */
int iVarArrayElementSize(VARARRAY header)
{
    if (header == NULL) {
        DFATAL((" iVarArrayElementSize: VARARRAY is NULL"));
    }
    return (header->size);
}

/*------------------------------------------------------
 *                iVarArrayElementCount
 *
 *
 */
int iVarArrayElementCount(VARARRAY header)
{
    if (header == NULL)
        return (0);
    return (header->count);
}

/*-----------------------------------------------------
 *               PVarArrayIndex
 */
char *PVarArrayIndex(VARARRAY header, int pos)
{
    if (header == NULL) {
        DFATAL((" PVarArrayIndex: VARARRAY is NULL"));
    }
    return (element(header, pos));

}

/*-----------------------------------------------------
 *      vaVarArrayCreate
 *
 *
 *      Create a new VARARRAY and initialize it.
 *      The caller must initialize the size of the elements
 *      of the VARARRAY when they create it.
 */

VARARRAY vaVarArrayCreate(int size)
{
    VARARRAY new;
    size_t alloc_size;

    MALLOC(new, VARARRAY, sizeof(HeaderStruct));

    if (new == NULL) {
        DFATAL(("vaVarArrayCreate: not enough memory "));
    }
    
    if (check_size_overflow(size, PORTION, &alloc_size)) {
        DFATAL(("vaVarArrayCreate: element size too large (%d bytes)", size));
    }
    
    MALLOC(new->data, char *, alloc_size);

    if (new->data == NULL) {
        FREE(new);
        DFATAL(("vaVarArrayCreate: not enough memory "));
    }
    new->count = 0;
    new->size = size;
    new->slot = PORTION;

    VARARRAY_LOG("Created VARARRAY: element_size=%d, initial_slots=%d, initial_bytes=%zu",
                 size, PORTION, alloc_size);

    return (new);
}

/*-----------------------------------------------------
 *      VarArrayDestroy
 *
 *      Destroy the VarArray, after this call the pointer will
 *      be undefined.
 */


void VarArrayDestroy(VARARRAY *header)
{

    if (*header == NULL) {
        DFATAL((" VarArrayDestroy: VARARRAY is NULL"));
    }
    
    VARARRAY_LOG("Destroying VARARRAY: count=%d, slots=%d",
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
 *      SSchott Code modified with Claude to improve memory allocations
 *
 */
void VarArrayAdd(VARARRAY header, GENP data)
{
    int new_slot;
    size_t new_alloc_size;
    
    if (header == NULL) {
        DFATAL((" VarArrayAdd: VARARRAY is NULL"));
    }
    
    if (header->count == header->slot) {
        if (header->slot < 100) {
            new_slot = header->slot + PORTION;
        } else if (header->slot < 10000) {
            new_slot = header->slot + header->slot / 2;
        } else if (header->slot < 1000000) {
            new_slot = header->slot + header->slot / 4;
        } else {
            new_slot = header->slot + header->slot / 10;
            if (new_slot - header->slot < 10000) {
                new_slot = header->slot + 10000;
            }
        }
        
        if (new_slot < header->slot) {
            DFATAL((" VarArrayAdd: slot count overflow"));
        }
        
        if (check_size_overflow(header->size, new_slot, &new_alloc_size)) {
            DFATAL((" VarArrayAdd: allocation size overflow (would need %d elements of %d bytes)",
                    new_slot, header->size));
        }
        
        VARARRAY_LOG("Growing VARARRAY: %d -> %d slots, %zu bytes",
                     header->slot, new_slot, new_alloc_size);
        
        REALLOC(header->data, char *, header->data, new_alloc_size);
        
        if (header->data == NULL) {
            DFATAL((" VarArrayAdd: realloc failed for %zu bytes", new_alloc_size));
        }
        
        header->slot = new_slot;
    }

    memcpy(element(header, header->count), (char *) data, header->size);
    header->count++;
}

/*-----------------------------------------------------
 *      vaVarArrayCopy
 *        Copy the VARARRAY.
 *
 *
 */
VARARRAY vaVarArrayCopy(VARARRAY header)
{
    VARARRAY new;
    size_t alloc_size;

    if (header == NULL) {
        DFATAL((" vaVarArrayCopy: VARARRAY is NULL"));
    }

    MALLOC(new, VARARRAY, sizeof(HeaderStruct));
    if (new == NULL) {
        DFATAL((" vaVarArrayCopy: cannot allocate header"));
    }

    new->size = header->size;
    new->count = header->count;
    new->slot = header->slot;

    if (check_size_overflow(new->size, header->slot, &alloc_size)) {
        FREE(new);
        DFATAL((" vaVarArrayCopy: allocation size overflow"));
    }

    MALLOC(new->data, char *, alloc_size);
    if (new->data == NULL) {
        FREE(new);
        DFATAL((" vaVarArrayCopy: cannot allocate %zu bytes", alloc_size));
    }

    memcpy(new->data, header->data, alloc_size);

    VARARRAY_LOG("Copied VARARRAY: %d elements, %zu bytes",
                 new->count, alloc_size);

    return (new);
}

/*-----------------------------------------------------
 *      vaVarArrayCopy2
 */
VARARRAY vaVarArrayCopy2(VARARRAY header1, VARARRAY header2)
{
    VARARRAY new;
    size_t copysize, alloc_size;
    int total_count;

    if (header1 == NULL || header2 == NULL) {
        DFATAL((" vaVarArrayCopy2: VARARRAY is NULL"));
    }
    if (header1->size != header2->size)
        DFATAL((" vaVarArrayCopy2: header sizes different\n"));

    if (header1->count > INT_MAX - header2->count) {
        DFATAL((" vaVarArrayCopy2: combined count overflow"));
    }
    total_count = header1->count + header2->count;

    MALLOC(new, VARARRAY, sizeof(HeaderStruct));
    if (new == NULL) {
        DFATAL((" vaVarArrayCopy2: cannot allocate header"));
    }

    new->size = header1->size;
    new->count = total_count;
    new->slot = SLOT_FOR_COUNT(new->count);

    if (check_size_overflow(new->size, new->slot, &alloc_size)) {
        FREE(new);
        DFATAL((" vaVarArrayCopy2: allocation size overflow"));
    }

    MALLOC(new->data, char *, alloc_size);
    if (new->data == NULL) {
        FREE(new);
        DFATAL((" vaVarArrayCopy2: cannot allocate %zu bytes", alloc_size));
    }

    copysize = (size_t)new->size * header1->count;
    memcpy(new->data, header1->data, copysize);
    memcpy(new->data + copysize, header2->data, (size_t)new->size * header2->count);

    VARARRAY_LOG("Merged 2 VARRAYs: %d + %d = %d elements, %zu bytes",
                 header1->count, header2->count, total_count, alloc_size);

    return (new);
}

void VarArrayInsertBeforeMore(VARARRAY header, int pos, int num)
{
    int shift, nslot;
    size_t alloc_size;
    char *h;

    if (header == NULL)
        DFATAL((" VarArrayInsertBeforeMore: VARARRAY is NULL"));

    if ((pos >= header->count) && (pos < 0))
        DFATAL((" VarArrayInsertBeforeMore: position=%d", pos));

    if (header->count > INT_MAX - num)
        DFATAL((" VarArrayInsertBeforeMore: count overflow"));

    nslot = SLOT_FOR_COUNT(header->count + num);

    if (header->slot != nslot) {
        if (check_size_overflow(header->size, nslot, &alloc_size)) {
            DFATAL((" VarArrayInsertBeforeMore: allocation size overflow"));
        }
        
        VARARRAY_LOG("Growing for insert: %d -> %d slots", header->slot, nslot);
        
        REALLOC(header->data, char *, header->data, alloc_size);
        if (header->data == NULL) {
            DFATAL((" VarArrayInsertBeforeMore: realloc failed"));
        }
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
    h = element(header, pos);
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
void VarArrayInsertBefore(VARARRAY header, int pos, GENP data)
{
    VarArrayInsertBeforeMore(header, pos, 1);
    memcpy(element(header, pos), (char *) data, header->size);
}



/*-----------------------------------------------------
 *        VarArrayDelete
 *
 *        Remove an element in a VARARRAY.
 *        Move the data below the one to be removed, up one
 *
 *        SSchott Code modified with Claude to improve memory cleanup
 */
void VarArrayDeleteMore(VARARRAY header, int pos, int num)
{
    int shift, nslot;
    size_t alloc_size;
    char *h;

    if (header == NULL) {
        DFATAL((" VarArrayDelete: VARARRAY is NULL"));
    }
    if (((pos + num) > header->count) || (pos < 0) || (num < 1)) {
        DFATAL((" VarArrayDelete: position=%-5d num=%-5d count=%-5d", pos, num,
                header->count));
    }
    header->count -= num;

    shift = num * header->size;

    h = element(header, pos);

    memmove(h, h + shift, (header->count - pos) * header->size);

    nslot = SLOT_FOR_COUNT(header->count);

    if (header->slot != nslot) {
        if (header->slot - nslot > header->slot / 10) {
            if (check_size_overflow(header->size, nslot, &alloc_size)) {
                return;  /* Can't shrink, but that's OK */
            }
            
            VARARRAY_LOG("Shrinking VARARRAY: %d -> %d slots", header->slot, nslot);
            
            REALLOC(header->data, char *, header->data, alloc_size);
            if (header->data == NULL) {
                DFATAL((" VarArrayDelete: realloc failed during shrink"));
            }
            header->slot = nslot;
        }
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
 */
void VarArraySetSize(VARARRAY header, int ncount)
{
    int nslot;
    size_t alloc_size;

    if (header == NULL) {
        DFATAL((" VarArraySetSize: VARARRAY is NULL"));
    }

    if (ncount < 0) {
        DFATAL((" VarArraySetSize: elements=%5d", ncount));
    }
    
    nslot = SLOT_FOR_COUNT(ncount);
    header->count = ncount;

    if (nslot != header->slot) {
        if (check_size_overflow(header->size, nslot, &alloc_size)) {
            DFATAL((" VarArraySetSize: allocation size overflow (need %d elements of %d bytes each)",
                    nslot, header->size));
        }
        
        VARARRAY_LOG("Resizing VARARRAY: %d -> %d elements, %d -> %d slots, %zu bytes",
                     header->count, ncount, header->slot, nslot, alloc_size);
        
        REALLOC(header->data, char *, header->data, alloc_size);
        if (header->data == NULL) {
            DFATAL((" VarArraySetSize: realloc failed for %zu bytes", alloc_size));
        }
        header->slot = nslot;
    }
}


/*-----------------------------------------------------
 *      PVarArrayDebugIndex
 *
 *      Return a pointer to the element within the VARARRAY, but
 *      first check the bounds.  Report an  if there
 *      is an out of bound access.
 */

GENP PVarArrayDebugIndex(VARARRAY header, int pos, char *file, int line)
{
    if (header == NULL) {
        DFATAL(("Attempting to access an invalid VARARRAY (%s line %d).", file,
                line));
    }
    if (header->count == 0) {
        DFATAL(("Attempting to access a no-data VARARRAY (%s line %d).", file,
                line));
    }
    if (pos < 0 || pos >= header->count) {
        DFATAL(("Attempted to access element: %d in a VARARRAY of size: %d",
                pos, header->count));
    }
    return (element(header, pos));
}
