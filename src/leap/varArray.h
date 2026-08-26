#ifndef  VARARRAY_H
#define  VARARRAY_H
 
/*
 *      File: varArray.h
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
 *     Designed by:    Christian Schafmeister                           *
 *     Author:         Christian Schafmeister                           *
 *                                                                      *
 *     VERSION: 1.0                                                     *
 *     Programmers:                                                     *
 *             Christian Schafmeister                                   *
 *             David Rivkin                                             *
 *                                                                      *
 *     Principal Investigator: Peter A. Kollman                         *
 *                                                                      *
 ************************************************************************
 *
 *      Description:
 *              A VARARRAY is a regular array whose size can increase.
 *
 *
 *       Body of varArray.[ch] modified by Vladimir Romanovski (1994).
 *
 */

#include "basics.h"

typedef struct {
        size_t  count;  /* real count of element in array  */ 
        size_t  size;   /* all elements have the same size */
        size_t  slot;   /* max available count of element before new realloc*/
        void    *data;
} HeaderStruct, *VARARRAY;

extern size_t      iVarArrayPointerToIndex(VARARRAY header, void *data);
#ifdef DEBUG
extern size_t      iVarArrayElementSize(VARARRAY header);
extern size_t      iVarArrayElementCount(VARARRAY header);
extern void     *PVarArrayIndex( VARARRAY header, size_t pos);
#else
static inline size_t iVarArrayElementSize(VARARRAY header) {
    return (header->size);
}
static inline size_t iVarArrayElementCount(VARARRAY header) {
    return header ? header->count : 0;
}
static inline void * PVarArrayIndex( VARARRAY header, size_t pos ) {
    return header->data + pos*header->size;
}
#endif

extern VARARRAY vaVarArrayCreate(size_t size);
extern void     VarArrayDestroy(VARARRAY *header);
extern void     VarArrayAdd(VARARRAY header, GENP data);
extern VARARRAY vaVarArrayCopy(VARARRAY header);
extern VARARRAY vaVarArrayCopy2(VARARRAY header1, VARARRAY header2);
extern void     VarArraySetSize(VARARRAY header, size_t ncount);
extern GENP     PVarArrayDebugIndex(VARARRAY header, size_t pos, 
                        const char *file, int line);

extern void     VarArrayInsertBefore(VARARRAY header, size_t pos, GENP data);
extern void     VarArrayInsertBeforeMore(VARARRAY header, size_t pos, size_t num);
     
#define VarArrayDelete(t,i) VarArrayDeleteMore(t,i,1)

#define PVAI(va,tc,i) ((tc*)PVarArrayIndex(va,(i)))

#endif  /* VARARRAY_H */




































