#ifndef  SYMBOL_H
#define  SYMBOL_H

#include <stdio.h>
#include "nab.h"

	/* function states:	*/
#define	FS_NEW	0
#define	FS_DECL	1
#define	FS_DEF	2

SYMREC_T	*entersym( int, char [], int, int, int, int );
SYMREC_T	*findsym( char [] );
int		openusyms( SYMREC_T * );
void		closeusyms( SYMREC_T * );
void	freelsyms( void );
void    dumpexpr( FILE *, struct node_t *, int );
void	dumpsyms( FILE *, int, int, int );
void	dsym( FILE *, char [], SYMREC_T * );

#endif

