#ifndef  CGEN_H
#define  CGEN_H

#include "nab.h"

int   CG_init( char [], int );
void  CG_exit( int );
void  CG_genend( void );
void  CG_genedefs( int );
void  CG_genestmts( int );
void  CG_gennl( void );
void  CG_genmain( void );
void  CG_genvardecl( NODE_T *, int, int, int );
void  CG_genassert( NODE_T * );
void  CG_gendebug( NODE_T * );
void  CG_genexpr( NODE_T * );
char  *CG_genop( char *, int );
void  CG_genrword( int );
void  CG_genfhdr( NODE_T *, NODE_T * );
void  CG_genplist( NODE_T * );
void  CG_genpdecls( void );
void  CG_genfstart( void );
void  CG_genfend( NODE_T * );
void  CG_genpdecl( NODE_T * );
char  *CG_gentemp( int );

void  fixexpr( NODE_T *, int, int );

#endif
