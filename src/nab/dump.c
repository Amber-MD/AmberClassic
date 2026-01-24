#include <stdio.h>
#include "nab.h"

void	dumpstrand( FILE *fp, STRAND_T *sp, int dres, int datom, int dbond );
void	dumpresidue( FILE *fp, RESIDUE_T *res, int datom, int dbond );
void	dumpatom( FILE *fp, RESIDUE_T *res, int anum, int dbond );

static	char	*attr2str( int );

void	dumpmolecule( FILE *fp, MOLECULE_T *mp, int dres, int datom, int dbond )
{
	STRAND_T	*sp;

	fprintf( fp, "Molecule %p {\n", mp );
	fprintf( fp, "\tframe\n" );
	fprintf( fp, "\t\torigin %8.3f %8.3f %8.3f\n",
		mp->m_frame[0][0], mp->m_frame[0][1], mp->m_frame[0][2] );
	fprintf( fp, "\t\tX-axis %8.3f %8.3f %8.3f\n",
		mp->m_frame[1][0], mp->m_frame[1][1], mp->m_frame[1][2] );
	fprintf( fp, "\t\tY-axis %8.3f %8.3f %8.3f\n",
		mp->m_frame[2][0], mp->m_frame[2][1], mp->m_frame[2][2] );
	fprintf( fp, "\t\tZ-axis %8.3f %8.3f %8.3f\n",
		mp->m_frame[3][0], mp->m_frame[3][1], mp->m_frame[3][2] );
	fprintf( fp, "\t%d strands\n", mp->m_nstrands );
	for( sp = mp->m_strands; sp; sp = sp->s_next ){
		dumpstrand( fp, sp, dres, datom, dbond );
	}
	fprintf( fp, "}\n" );
}

void	dumpstrand( FILE *fp, STRAND_T *sp, int dres, int datom, int dbond )
{
	int	r;

	fprintf( fp, "\tStrand %p {\n", sp );
	fprintf( fp, "\t\tname %s\n", sp->s_strandname );
	fprintf( fp, "\t\tnum  %d\n", sp->s_strandnum );
	fprintf( fp, "\t\tslct %s\n", attr2str( sp->s_attr ) );
	fprintf( fp, "\t\tmol  %p\n", sp->s_molecule );
	fprintf( fp, "\t\tnext %p\n", sp->s_next );
	fprintf( fp, "\t\tnres   %4d\n", sp->s_nresidues );
	fprintf( fp, "\t\trsiz   %4d\n", sp->s_res_size );
	fprintf( fp, "\t\tres    %p\n", sp->s_residues );
	if( dres ){
		for( r = 0; r < sp->s_nresidues; r++ )
			dumpresidue( fp, sp->s_residues[ r ], datom, dbond );
	}
	fprintf( fp, "\t}\n" );
}

void	dumpresidue( FILE *fp, RESIDUE_T *res, int datom, int dbond )
{
	int	a;
	EXTBOND_T	*ep;
	char	*wp;
	char	name[ 20 ];

	fprintf( fp, "\tResidue %p {\n", res );
	fprintf( fp, "\t\tname  %8s\n", res->r_resname );
	fprintf( fp, "\t\tresid %8s\n", res->r_resid );
	fprintf( fp, "\t\tnumid %8d\n", res->r_num );
	fprintf( fp, "\t\tnum   %8d\n", res->r_resnum );
	fprintf( fp, "\t\ttnum  %8d\n", res->r_tresnum );
	if( res->r_kind == RT_UNDEF )
		wp = "Undef";
	else if( res->r_kind == RT_DNA )
		wp = "dna";
	else if( res->r_kind == RT_RNA )
		wp = "rna";
	else if( res->r_kind == RT_AA )
		wp = "aa";
	else{
		sprintf( name, "rk(%d)", res->r_kind );
		wp = name;
	}
	fprintf( fp, "\t\tkind  %8s\n", wp );
	if( res->r_atomkind == RAT_UNDEF )
		wp = "Undef";
	else if( res->r_atomkind == RAT_UNITED )
		wp = "united";
	else if( res->r_atomkind == RAT_ALLATOM )
		wp = "all";
	else{
		sprintf( name, "ak(%d)", res->r_atomkind );
		wp = name;
	}
	fprintf( fp, "\t\takind %8s\n", wp );
	fprintf( fp, "\t\tslct %s\n", attr2str( res->r_attr ) );
	fprintf( fp, "\t\tnext  %p\n", res->r_next );
	fprintf( fp, "\t\tstr   %p\n", res->r_strand );
	fprintf( fp, "\t\tnatm  %4d\n", res->r_natoms );
	fprintf( fp, "\t\tatoms %p\n", res->r_atoms );
	if( res->r_atoms != NULL && datom ){
		for( a = 0; a < res->r_natoms; a++ )
			dumpatom( fp, res, a, dbond );
	}
	if( dbond ){
		if( (ep = res->r_extbonds) ){
			fprintf( fp, "\tExternal bonds:\n" );
			for( ; ep; ep = ep->eb_next )
				fprintf( fp, "\t\t%d->%d:%d\n",
					ep->eb_anum,
					ep->eb_rnum, ep->eb_ranum );
		}else{
			fprintf( fp, "\tNo external bonds\n" );
		}
	}
	fprintf( fp, "\t}\n" );
}

void	dumpatom( FILE *fp, RESIDUE_T *res, int anum, int dbond )
{
	ATOM_T	*ap;
	int	a, ai;

	ap = &res->r_atoms[ anum ];
	ai = ( res->r_aindex ) ? res->r_aindex[ anum ] : -1;
	fprintf( fp, "ATM[%3d,%3d,%p,%p] %-4s %s %8.3f%8.3f%8.3f",
		anum, ai, ap, ap->a_residue, ap->a_atomname,
		attr2str( ap->a_attr ),
		ap->a_pos[ 0 ], ap->a_pos[ 1 ], ap->a_pos[ 2 ] );
	if( dbond ){
		fprintf( fp, " [" );
		for( a = 0; a < ap->a_nconnect; a++ ){
			fprintf( fp, " %2d", ap->a_connect[ a ] );
		}
		fprintf( fp, "]" );
	}
	putc( '\n', fp );
}

void	dumpmatrix( FILE *fp, MATRIX_T mat )
{
	int	i, j;

	
	for( i = 0; i < 4; i++ ){
		for( j = 0; j < 4; j++ )
			fprintf( fp, " %8.3f", mat[ i ][ j ] );
		putc( '\n', fp );
	}
}

static	char	*attr2str( int attr )
{
	char	*sp;
	static	char	str[ 10 ];

	sp = str;
	*sp++ = ( attr & AT_WORK ) ? 'W' : '-';
	*sp++ = ( attr & AT_SELECTED ) ? 'S' : '-';
	*sp++ = ( attr & AT_SELECT ) ? 's' : '-';
	*sp = '\0';
	return( str );
}
