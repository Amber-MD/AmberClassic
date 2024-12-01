
#ifndef _GNU_SOURCE
#define	_GNU_SOURCE
#endif

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#include "nab.h"
#include "select_atoms.h"
#include "molutil.h"

#define	LAST	(-1)

#define	AEXPR_SIZE	1000
#define	REXPR_SIZE	1000
#define	EXPBUF_SIZE	1000

static	int	eval_1_aexpr( MOLECULE_T *, char [] );
static	int	atom_in_1_aexpr( ATOM_T *, char [] );
static	void	get_aex_parts(const char *, char **, char **, char **);
static	int	is_pattern(char [], int *, int *, int *, int *);
static	void	select_all( MOLECULE_T * );
static	void	clear_select( MOLECULE_T * );
static	void	clear_work( MOLECULE_T * );
static	void	or_select( MOLECULE_T * );
static	void	set_select( MOLECULE_T * );
static	void	match_str_pat( MOLECULE_T *, char [] );
static	int	atom_in_str_pat( ATOM_T *, char [] );
static	void	match_str_range(MOLECULE_T *, int, int);
static	int	atom_in_str_range(ATOM_T *, int, int);
static	void	match_res_pat( MOLECULE_T *, char [] );
static	int	atom_in_res_pat( ATOM_T *, char [] );
static	void	match_res_range(MOLECULE_T *, int, int, int, int);
static	int	atom_in_res_range(ATOM_T *, int, int, int, int);
static	void	match_atom_pat( MOLECULE_T *, char [] );
static	int	atom_in_atom_pat( ATOM_T *, char [] );
static	void	aexpr2rexpr( char [], char [] );

/* strnlen.c - determine the length of a fixed-size string
 * Andrew Ho (andrew@zeuscat.com)
 * 
 * Implements strnlen(), which is defined on some BSD systems and Linux
 * but not on Solaris. This is just like strlen(), but you provide a
 * maximum number of bytes to count, which saves time if you know
 * what the maximum allowed size of the string is.
 *
 */

/* just like strlen(3), but cap the number of bytes we count */
static size_t strnlen2(const char *s, size_t max) {
    register const char *p;
    for(p = s; *p && max--; ++p);
    return(p - s);
}
/* A replacement function, for systems that lack strndup. */
static char *strndup2 (char const *s, size_t n)
{
  size_t len = strnlen2 (s, n);
  char *new = malloc (len + 1);

  if (new == NULL)
    return NULL;

  new[len] = '\0';
  return memcpy (new, s, len);
}

int	setpoint( MOLECULE_T *mol, char aexpr[], POINT_T point )
{
	int	r, a, ta;
	STRAND_T	*sp;
	RESIDUE_T	*res;
	ATOM_T		*ap;
	REAL_T	x, y, z;

	select_atoms( mol, aexpr );

	/* compute the point as the CM of the selected atoms	*/
	x = y = z = 0.0;
	for( ta = 0, sp = mol->m_strands; sp; sp = sp->s_next ){
		if( AT_SELECT & sp->s_attr ){
			for( r = 0; r < sp->s_nresidues; r++ ){
				res = sp->s_residues[ r ];
				if( AT_SELECT & res->r_attr ){
					for( a = 0; a < res->r_natoms; a++ ){
						ap = &res->r_atoms[ a ];
						if( AT_SELECT & ap->a_attr ){
							x += ap->a_pos[ 0 ];
							y += ap->a_pos[ 1 ];
							z += ap->a_pos[ 2 ];
							ta++;
						}
					}
				}
			}
		}
	}

	if( ta == 0 ){
		fprintf( stderr, "setpoint: %s: no atoms selected\n", aexpr );
		return( 1 );
	}else{
		point[ 0 ] = x / ta;
		point[ 1 ] = y / ta;
		point[ 2 ] = z / ta;
	}

	return( 0 );
}

int	select_atoms( MOLECULE_T *mol, char aex[] )
{
	char	aexpr[ AEXPR_SIZE ];
	char	*aep, *n_aep;
	int	ael;

	if( aex == NULL ){
		select_all( mol );
		return( 0 );
	}

	clear_work( mol );
	clear_select( mol );

	for( aep = aex, n_aep = strchr( aep, '|' ); aep; ){
		if( n_aep ){
			ael = n_aep - aep;
			n_aep++;
		}else
			ael = strlen( aep );
		if( ael >= AEXPR_SIZE ){
			fprintf( stderr,
				"select_atoms: atom-expr too complicated\n" );
			return( 1 );
		}
		strncpy( aexpr, aep, ael );
		aexpr[ ael ] = '\0';
		eval_1_aexpr( mol, aexpr );
		or_select( mol );
		aep = n_aep;
		if( aep )
			n_aep = strchr( aep, '|' );
		clear_select( mol );
	}
	set_select( mol );

	return( 0 );
}

int	atom_in_aexpr( ATOM_T *ap, char aex[] )
{
	char	aexpr[ AEXPR_SIZE ];
	char	*aep, *n_aep;
	int	ael;

	if( aex == NULL )
		return( 0 );

	for( aep = aex, n_aep = strchr( aep, '|' ); aep; ){
		if( n_aep ){
			ael = n_aep - aep;
			n_aep++;
		}else
			ael = strlen( aep );
		if( ael >= AEXPR_SIZE ){
			fprintf( stderr,
				"atom_in_aexpr: atom-expr too complicated\n" );
			return( 0 );
		}
		strncpy( aexpr, aep, ael );
		aexpr[ ael ] = '\0';
		if( atom_in_1_aexpr( ap, aexpr ) )
			return( 1 );
		aep = n_aep;
		if( aep )
			n_aep = strchr( aep, '|' );
	}
	return( 0 );
}

void	set_attr_if( MOLECULE_T *mol, int attr, int i_attr )
{
	int		a, r;
	STRAND_T	*sp;
	RESIDUE_T	*res;
	ATOM_T		*ap;

	for( sp = mol->m_strands; sp; sp = sp->s_next ){
		sp->s_attr |= ( sp->s_attr & i_attr ) ? attr : 0;
		for( r = 0; r < sp->s_nresidues; r++ ){
			res = sp->s_residues[ r ];
			res->r_attr |= ( res->r_attr & i_attr ) ?
				attr : 0;
			for( a = 0; a < res->r_natoms; a++ ){
				ap = &res->r_atoms[ a ];
				ap->a_attr |= ( ap->a_attr & i_attr ) ?
					attr : 0;
			}
		}
	}
}

void	clear_attr( MOLECULE_T *mol, int attr )
{
	int		a, r;
	STRAND_T	*sp;
	RESIDUE_T	*res;

	for( sp = mol->m_strands; sp; sp = sp->s_next ){
		sp->s_attr &= ~attr;
		for( r = 0; r < sp->s_nresidues; r++ ){
			res = sp->s_residues[ r ];
			res->r_attr &= ~attr;
			for( a = 0; a < res->r_natoms; a++ )
				res->r_atoms[ a ].a_attr &= ~attr;
		}
	}
}

static	int	eval_1_aexpr( MOLECULE_T *mol, char aex[] )
{
	char	*spart = NULL;
	char	*rpart = NULL;
	char	*apart = NULL;
	char	*pp, *e_pp, *wp;
	int	a_lo, a_hi;
	int	lo, hi;
	
	get_aex_parts(aex, &spart, &rpart, &apart);

	if( spart ){
		for(e_pp = pp = spart; e_pp; ){
			wp = NULL;
			if((e_pp = strchr(pp, ',')) != NULL){
				wp = strndup2(pp, e_pp - pp);
				pp = e_pp + 1;
			}else{
				wp = strdup(pp);
				e_pp = NULL;
			}
			if(is_pattern(wp, &a_lo, &lo, &a_hi, &hi))
				match_str_pat(mol, wp);
			else
				match_str_range(mol, lo, hi);
			if(wp != NULL)
				free(wp);
		}
	}else
		match_str_range(mol, 1, LAST);

	if( rpart ){
		for(e_pp = pp = rpart; e_pp; ){
			wp = NULL;
			if((e_pp = strchr(pp, ',')) != NULL){
				wp = strndup2(pp, e_pp - pp);
				pp = e_pp + 1;
			}else{
				wp = strdup(pp);
				e_pp = NULL;
			}
			if(is_pattern(wp, &a_lo, &lo, &a_hi, &hi))
				match_res_pat(mol, wp);
			else
				match_res_range(mol, a_lo, lo, a_hi, hi);
			if(wp != NULL)
				free(wp);
		}
	}else
		match_res_range(mol, FALSE, 1, FALSE, LAST);

	if( apart ){
		for(e_pp = pp = apart; e_pp; ){
			wp = NULL;
			if((e_pp = strchr(pp, ',')) != NULL){
				wp = strndup2(pp, e_pp - pp);
				pp = e_pp + 1;
			}else{
				wp = strdup(pp);
				e_pp = NULL;
			}
			if(is_pattern(wp, &a_lo, &lo, &a_hi, &hi))
				match_atom_pat(mol, wp);
			else
				fprintf(stderr, "atom range not allowed\n");
			if(wp != NULL)
				free(wp);
		}
	}else
		match_atom_pat(mol, "*");

	if(spart != NULL)
		free(spart);
	if(rpart != NULL)
		free(rpart);
	if(apart != NULL)
		free(apart);

	return(0);
}

static	int	atom_in_1_aexpr( ATOM_T *ap, char aex[] )
{
	char	*spart = NULL;
	char	*rpart = NULL;
	char	*apart = NULL;
	char	*pp, *e_pp, *wp;
	int	a_lo, a_hi;
	int	lo, hi;
	int	found;
	int	err = 0;
	int	rval = 0;

	get_aex_parts(aex, &spart, &rpart, &apart);

	found = 1;
	if( spart ){
		for(found = 0, e_pp = pp = spart; e_pp && !found; ){
			wp = NULL;
			if((e_pp = strchr(pp, ',')) != NULL){
				wp = strndup2(pp, e_pp - pp);
				pp = e_pp + 1;
			}else{
				wp = strdup(pp);
				e_pp = NULL;
			}
			if(is_pattern(wp, &a_lo, &lo, &a_hi, &hi)){
				if(atom_in_str_pat(ap, wp))
					found = 1;
			}else if(atom_in_str_range(ap, lo, hi))
				found = 1;
			if(wp != NULL)
				free(wp);
		}
	}

	 if(found && rpart){
		for(found = 0, e_pp = pp = rpart; e_pp && !found; ){
			wp = NULL;
			if((e_pp = strchr(pp, ',')) != NULL){
				wp = strndup2(pp, e_pp - pp);
				pp = e_pp + 1;
			}else{
				wp = strdup(pp);
				e_pp = NULL;
			}
			if(is_pattern(wp, &a_lo, &lo, &a_hi, &hi)){
				if(atom_in_res_pat(ap, wp))
					found = 1;
			}else if(atom_in_res_range(ap, a_lo, lo, a_hi, hi))
				found = 1;
			if(wp != NULL)
				free(wp);
		}
	}

	if(found && apart){
		for(err = 0, found = 0, e_pp = pp = apart; e_pp && !found && !err; ){
			wp = NULL;
			if((e_pp = strchr(pp, ',')) != NULL){
				wp = strndup2(pp, e_pp - pp);
				pp = e_pp + 1;
			}else{
				wp = strdup(pp);
				e_pp = NULL;
			}
			if(is_pattern(wp, &a_lo, &lo, &a_hi, &hi)){
				if(atom_in_atom_pat(ap, wp))
					found = 1;
			}else{
				fprintf( stderr, "atom range not allowed\n" );
				err = 1;
			}
			if(wp != NULL)
				free(wp);
		}
	}

	if(spart != NULL)
		free(spart);
	if(rpart != NULL)
		free(rpart);
	if(apart != NULL)
		free(apart);

	return found;
}

static	void	get_aex_parts(const char *aex, char **spart, char **rpart, char **apart)
{
	const char	*aep;
	const char	*colon;

	*spart = *rpart = *apart = NULL;

	aep = aex;
	if(*aep == ':'){
		*spart = NULL;
		aep++;
	}else{
		if((colon = strchr(aep, ':')) != NULL){
			*spart = strndup2(aep, colon - aep);
			aep = colon + 1;
		}else{
			*spart = strdup(aep);
			aep += strlen(aep);
		}
	}
	if(*aep == ':'){
		*rpart = NULL;
		aep++;
	}else{
		if((colon = strchr(aep, ':')) != NULL){
			*rpart = strndup2(aep, colon - aep);
			aep = colon + 1;
		}else{
			*rpart = strdup(aep);
			aep += strlen(aep);
		}
	}
	if(*aep)
		*apart = strdup(aep);
	else
		*apart = NULL;
}

static	int	is_pattern(char item[], int *a_lo, int *lo, int *a_hi, int *hi)
{
	int	val;
	char	*ip;

/*	working
	if(!isdigit(*item) && *item != '-' && *item != '#')
		return  TRUE;

	if(isdigit(*item)){
		for( val = 0, ip = item; isdigit( *ip ); ip++ )
			val = 10 * val + *ip - '0';
		*lo = val;
		if(!*ip){
			*hi = *lo;
			return FALSE;
		}else if(*ip == '-')
			ip++;
		if(!*ip){
			*hi = LAST;
			return  FALSE;
		}else if(!isdigit(*ip))
			return TRUE;
		for(val = 0; isdigit(*ip); ip++)
			val = 10 * val + *ip - '0'; 
		*hi = val;
		return  *ip;
	}else{
		*lo = 1;
		ip = &item[1];
	}
	if(!*ip){
		*hi = LAST;
		return FALSE;
	}else if(isdigit(*ip)){
		for(val = 0; isdigit(*ip); ip++)
			val = 10 * val + *ip - '0';
		*hi = val;
		return  *ip;
	}
	return FALSE;
*/
	// ranges start with 0-9, - or #, anything else is a pattern
	if(!isdigit(*item) && *item != '-' && *item != '#')
		return  TRUE;

	*a_lo = *a_hi = FALSE;
	if(*item == '#'){
		// # must be followed by a number, if not it's a pattern
		if(!isdigit(item[1]))
			return TRUE;
		*a_lo = TRUE;
		ip = &item[1];
	}else
		ip = item;

	if(isdigit(*ip)){	// must be n, n-, n-m
		for(val = 0; isdigit(*ip); ip++)
			val = 10 * val + *ip - '0';
		*lo = val;
		if(!*ip){	// n	-> n-n
			*a_hi = *a_lo;
			*hi = *lo;
			return FALSE;
		}else if(*ip == '-')
			ip++;
		if(!*ip){	// n-	-> n-LAST
			*a_hi = *a_lo;
			*hi = LAST;
			return FALSE;
		}
		if(*ip == '#'){
			// # must be followed by a number
			if(!isdigit(ip[1]))
				return TRUE;
			*a_hi = TRUE;
			ip++;
		}
		// n-m
		for(val = 0; isdigit(*ip); ip++)
			val = 10 * val + *ip - '0';
		*hi = val;
		return *ip;	// it's a pattern if anything follws the last digit
	}else{	// must be -, as in -, -m
		*lo = 1;
		ip++;
	}
	if(!*ip){	// must be - 	-> 1-LAST
		*hi = LAST;
		*a_hi = *a_lo;
		return FALSE;
	}
	if(*ip == '#'){
		// # must be followed by a number
		if(!isdigit(ip[1]))
			return TRUE;
		*a_hi = TRUE;
		ip++;
	}
	if(isdigit(*ip)){	// must be -m
		for(val = 0; isdigit(*ip); ip++)
			val = 10 * val + *ip - '0';
		*hi = val;
	}

	return *ip;	// it's a pattern if anything follows the last digit
}

static	void	select_all( MOLECULE_T *mol )
{
	int		a, r;
	STRAND_T	*sp;
	RESIDUE_T	*res;

	for( sp = mol->m_strands; sp; sp = sp->s_next ){
		sp->s_attr |= AT_SELECT;
		for( r = 0; r < sp->s_nresidues; r++ ){
			res = sp->s_residues[ r ];
			res->r_attr |= AT_SELECT;
			for( a = 0; a < res->r_natoms; a++ )
				res->r_atoms[ a ].a_attr |= AT_SELECT;
		}
	}
}

static	void	clear_select( MOLECULE_T *mol )
{
	int		a, r;
	STRAND_T	*sp;
	RESIDUE_T	*res;

	for( sp = mol->m_strands; sp; sp = sp->s_next ){
		sp->s_attr &= ~AT_SELECT;
		for( r = 0; r < sp->s_nresidues; r++ ){
			res = sp->s_residues[ r ];
			res->r_attr &= ~AT_SELECT;
			for( a = 0; a < res->r_natoms; a++ )
				res->r_atoms[ a ].a_attr &= ~AT_SELECT;
		}
	}
}

static	void	clear_work( MOLECULE_T *mol )
{
	int		a, r;
	STRAND_T	*sp;
	RESIDUE_T	*res;

	for( sp = mol->m_strands; sp; sp = sp->s_next ){
		sp->s_attr &= ~AT_WORK;
		for( r = 0; r < sp->s_nresidues; r++ ){
			res = sp->s_residues[ r ];
			res->r_attr &= ~AT_WORK;
			for( a = 0; a < res->r_natoms; a++ )
				res->r_atoms[ a ].a_attr &= ~AT_WORK;
		}
	}
}

static	void	or_select( MOLECULE_T *mol )
{
	int		a, r;
	STRAND_T	*sp;
	RESIDUE_T	*res;
	ATOM_T		*ap;

	for( sp = mol->m_strands; sp; sp = sp->s_next ){
		sp->s_attr |= ( sp->s_attr & AT_SELECT ) ? AT_WORK : 0;
		for( r = 0; r < sp->s_nresidues; r++ ){
			res = sp->s_residues[ r ];
			res->r_attr |= ( res->r_attr & AT_SELECT ) ?
				AT_WORK : 0;
			for( a = 0; a < res->r_natoms; a++ ){
				ap = &res->r_atoms[ a ];
				ap->a_attr |= ( ap->a_attr & AT_SELECT ) ?
					AT_WORK : 0;
			}
		}
	}
}

static	void	set_select( MOLECULE_T *mol )
{
	int		a, r;
	STRAND_T	*sp;
	RESIDUE_T	*res;
	ATOM_T		*ap;

	for( sp = mol->m_strands; sp; sp = sp->s_next ){
		sp->s_attr |= ( sp->s_attr & AT_WORK ) ? AT_SELECT : 0;
		for( r = 0; r < sp->s_nresidues; r++ ){
			res = sp->s_residues[ r ];
			res->r_attr |= ( res->r_attr & AT_WORK ) ?
				AT_SELECT : 0;
			for( a = 0; a < res->r_natoms; a++ ){
				ap = &res->r_atoms[ a ];
				ap->a_attr |= ( ap->a_attr & AT_WORK ) ?
					AT_SELECT : 0;
			}
		}
	}
}

static	void	match_str_pat( MOLECULE_T *mol, char pat[] )
{
	char	rexpr[ REXPR_SIZE ];
	char	expbuf[ EXPBUF_SIZE ];
	STRAND_T	*sp;

	aexpr2rexpr( pat, rexpr );
	compile( rexpr, expbuf, &expbuf[ EXPBUF_SIZE ], '\0' );
	for( sp = mol->m_strands; sp; sp = sp->s_next ){
		sp->s_attr |= step( sp->s_strandname, expbuf ) ? AT_SELECT : 0;
	}
}

static	int	atom_in_str_pat( ATOM_T *ap, char pat[] )
{
	RESIDUE_T	*res;
	STRAND_T	*sp;
	char	rexpr[ REXPR_SIZE ];
	char	expbuf[ EXPBUF_SIZE ];

	res = ap->a_residue;
	sp = res->r_strand;
	aexpr2rexpr( pat, rexpr );
	compile( rexpr, expbuf, &expbuf[ EXPBUF_SIZE ], '\0' );
	return(	step( sp->s_strandname, expbuf ) );
}

static	void	match_str_range(MOLECULE_T *mol, int lo, int hi)
{
	int		m;
	STRAND_T	*sp;

	if( hi == UNDEF )
		hi = mol->m_nstrands;
	for( m = 1, sp = mol->m_strands; m <= mol->m_nstrands;
		m++, sp = sp->s_next ){
		if( lo <= m && m <= hi )
			sp->s_attr |= AT_SELECT;
	}
}

static	int	atom_in_str_range(ATOM_T *ap, int lo, int hi)
{
	int		m;
	RESIDUE_T	*res;
	STRAND_T	*sp, *sp1;
	MOLECULE_T	*mol;

	res = ap->a_residue;
	sp = res->r_strand;
	mol = sp->s_molecule;
	if( hi == UNDEF )
		hi = mol->m_nstrands;
	for( m = 1, sp1 = mol->m_strands; m <= mol->m_nstrands; m++, sp1 = sp1->s_next ){
		if( sp == sp1 ){
			if( lo <= m && m <= hi )
				return( 1 );
		}
	}
	return( 0 );
}

static	void	match_res_pat( MOLECULE_T *mol, char pat[] )
{
	char	rexpr[ REXPR_SIZE ];
	char	expbuf[ EXPBUF_SIZE ];
	int		r;
	STRAND_T	*sp;
	RESIDUE_T	*res;

	aexpr2rexpr( pat, rexpr );
	compile( rexpr, expbuf, &expbuf[ EXPBUF_SIZE ], '\0' );
	for( sp = mol->m_strands; sp; sp = sp->s_next ){
		if( AT_SELECT & sp->s_attr ){
			for( r = 0; r < sp->s_nresidues; r++ ){
				res = sp->s_residues[ r ];
				res->r_attr |= step( res->r_resname, expbuf ) ?
					AT_SELECT : 0;
			}
		}
	}
}

static	int	atom_in_res_pat( ATOM_T *ap, char pat[] )
{
	RESIDUE_T	*res;
	char	rexpr[ REXPR_SIZE ];
	char	expbuf[ EXPBUF_SIZE ];

	res = ap->a_residue;
	aexpr2rexpr( pat, rexpr );
	compile( rexpr, expbuf, &expbuf[ EXPBUF_SIZE ], '\0' );
	return( step( res->r_resname, expbuf ) );
}

static	void	match_res_range(MOLECULE_T *mol, int a_lo, int lo, int a_hi, int hi)
{
	int		r, rlo, rhi;
	STRAND_T	*sp;
	RESIDUE_T	*res;

	if(a_lo && a_hi){	// both are absolute res. numbers
		if(!mol->m_nvalid)
			upd_molnumbers(mol);
		for(sp = mol->m_strands; sp; sp = sp->s_next){
			if( AT_SELECT & sp->s_attr ){
				if(sp->s_nresidues == 0)
					continue;
				else if(sp->s_residues[sp->s_nresidues-1]->r_tresnum < lo)
					continue;
				else if(sp->s_residues[0]->r_tresnum > hi)
					break;
				for(r = 0; r < sp->s_nresidues; r++){
					res = sp->s_residues[r];
					if(lo <= res->r_tresnum && res->r_tresnum <= hi)
						res->r_attr |= AT_SELECT;
				}
			}
		}
	}else if(!a_lo && !a_hi){	// both are relative res. numbers
		for( sp = mol->m_strands; sp; sp = sp->s_next ){
			if( AT_SELECT & sp->s_attr ){
				rhi = ( hi == UNDEF ) ? sp->s_nresidues : hi;
				for( r = 0; r < sp->s_nresidues; r++ ){
					res = sp->s_residues[ r ];
					if( lo <= r + 1 && r + 1 <= rhi )
						res->r_attr |= AT_SELECT;
				}
			}
		}
	}else	// mixed abs + rel res numbers is an error
		fprintf(stderr, "match_res_range: mixed absolute and relative res numbers not allowed\n");
}

static	int	atom_in_res_range(ATOM_T *ap, int a_lo, int lo, int a_hi, int hi)
{
	int		r, rhi;
	MOLECULE_T	*mp;
	STRAND_T	*sp;
	RESIDUE_T	*res, *res1;

	res = ap->a_residue;
	sp = res->r_strand;
	mp = sp->s_molecule;
	
	if(a_lo && a_hi){
		if(!mp->m_nvalid)
			upd_molnumbers(mp);
		if(res->r_tresnum >= lo && res->r_tresnum <= hi)
			return TRUE;
	}else if(!a_lo && !a_hi){
		rhi = ( hi == UNDEF ) ? sp->s_nresidues : hi;
		for( r = 0; r < sp->s_nresidues; r++ ){
			res1 = sp->s_residues[ r ];
			if( res == res1 ){
				if( lo <= r + 1 && r + 1 <= rhi )
					return TRUE;
			}
		}
	}else	// mixed abs + rel res numbers is an error
		fprintf(stderr, "atom_in_res_range: mixed absolute and relative res numbers not allowed\n");
	return FALSE;
}

static	void	match_atom_pat( MOLECULE_T *mol, char pat[] )
{
	char	rexpr[ REXPR_SIZE ];
	char	expbuf[ EXPBUF_SIZE ];
	int		r, a;
	STRAND_T	*sp;
	RESIDUE_T	*res;
	ATOM_T		*ap;

	aexpr2rexpr( pat, rexpr );
	compile( rexpr, expbuf, &expbuf[ EXPBUF_SIZE ], '\0' );
	for( sp = mol->m_strands; sp; sp = sp->s_next ){
		if( AT_SELECT & sp->s_attr ){
			for( r = 0; r < sp->s_nresidues; r++ ){
				res = sp->s_residues[ r ];
				if( AT_SELECT & res->r_attr ){
					for( a = 0; a < res->r_natoms; a++ ){
						ap = &res->r_atoms[ a ];
						ap->a_attr |= 
						    step(ap->a_atomname,expbuf)?
						    AT_SELECT : 0;
					}
				}
			}
		}
	}
}

static	int	atom_in_atom_pat( ATOM_T *ap, char pat[] )
{
	char	rexpr[ REXPR_SIZE ];
	char	expbuf[ EXPBUF_SIZE ];

	aexpr2rexpr(pat, rexpr);
	compile(rexpr, expbuf, &expbuf[ EXPBUF_SIZE ], '\0');
	return  step( ap->a_atomname, expbuf );
}

static	void	aexpr2rexpr(char aexpr[], char rexpr[])
{
	char	*aep, *rep;

	rep = rexpr;
	*rep++ = '^';
	for(aep = aexpr; *aep; aep++){
		if(*aep == '*'){
			*rep++ = '.';
			*rep++ = '*';
		}else if(*aep == '?')
			*rep++ = '.';
		else
			*rep++ = *aep;
	}
	*rep++ = '$';
	*rep = '\0';
}
