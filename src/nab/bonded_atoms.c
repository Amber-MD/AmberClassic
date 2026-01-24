#include <stdio.h>

#include "nab.h"

int	bonded_atoms( ATOM_T *a, ATOM_T *neighbors[] )
{
	STRAND_T	*sp;
	RESIDUE_T	*res, *res1;
	EXTBOND_T	*ebp;
	int		nb, ai, aj;

	res = a->a_residue;
	sp = res->r_strand;
	ai = a - res->r_atoms;
	for( nb = 0; nb < a->a_nconnect; nb++ ){
		aj = a->a_connect[ nb ];
		neighbors[ nb ] = &res->r_atoms[ aj ];
	}
	for( ebp = res->r_extbonds; ebp; ebp = ebp->eb_next ){
		if( ebp->eb_anum - 1 == ai ){
			res1 = sp->s_residues[ ebp->eb_rnum - 1 ];
			neighbors[ nb ] = &res1->r_atoms[ ebp->eb_ranum - 1 ];
			nb++;
		}
	}
	return( nb );
}
