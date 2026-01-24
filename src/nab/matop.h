#ifndef MATOP_H
#define MATOP_H 
#include "nabtypes.h"

/* Functions defined in matop.c */

int	MAT_fprint( FILE *, int, MATRIX_T [] );
int	MAT_sprint( char [], int, MATRIX_T [] );
int	MAT_fscan( FILE *, int, MATRIX_T [] );
int	MAT_sscan( char [], int, MATRIX_T [] );
REF_MATRIX_T	MAT_concat( MATRIX_T, MATRIX_T );
int	MAT_count( char [] );
char	*MAT_getsyminfo( void );
int	MAT_istrue( MATRIX_T );

#endif
