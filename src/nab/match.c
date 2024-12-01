#include <stdio.h>

#include "nab.h"

#define	EXPBUF_SIZE	256
static	char	expbuf[ EXPBUF_SIZE ];

int	NAB_aematch( ATOM_T *ap, char aex[] )
{
	STRAND_T	*sp;
	RESIDUE_T	*rp;
	MOLECULE_T	*mp;

	rp = ap->a_residue;
	sp = rp->r_strand;
	mp = sp->s_molecule;
	return( atom_in_aexpr( ap, aex ) );
/*
	select_atoms( mp, aex );
	return( AT_SELECT & ap->a_attr );
*/
}

int	NAB_rematch( char str[], char pat[] )
{

	compile( pat, expbuf, &expbuf[ EXPBUF_SIZE ], '\0' );
	return( step( str, expbuf ) );
}
