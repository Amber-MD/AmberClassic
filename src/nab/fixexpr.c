#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "nab.h"
#include "nabgrm.tab.h"
#include "cgen.h"
#include "errormsg.h"
#include "symbol.h"

	/* literals	*/
extern	VALUE_T	v_0;
extern	VALUE_T	v_1;
extern	VALUE_T	v_NULL;
extern	VALUE_T	v_ainit;
extern	VALUE_T	v_hinit;
extern	VALUE_T	v_hin;
extern	VALUE_T	v_hfirst;
extern	VALUE_T	v_hnext;
extern	VALUE_T	v_mfirst;
extern	VALUE_T	v_mnext;
extern	VALUE_T	v_rfirst;
extern	VALUE_T	v_rnext;
extern	VALUE_T	v_afirst;
extern	VALUE_T	v_anext;

extern	char	c_funcname[];
extern	int	CG_get_attr_access( NODE_T *, NODE_T *, NODE_T * );

extern	int	debug;

extern	SYMREC_T	*astk[];
extern	int		astkp;

NODE_T	*node();

static	NODE_T	*mk_iexpr( NODE_T * );
static	int	is_getline( NODE_T * );
static	int	is_trig( NODE_T * );
static	int	is_index( NODE_T * );
static	int	is_molio( NODE_T * );
static	int	is_substitute( NODE_T * );
static	void	fixparm( NODE_T *, int );
static	void	fixid( NODE_T *, int );
static	NODE_T	*in2pre( char [], NODE_T *, NODE_T *, NODE_T * );
static	void	fix_getline( NODE_T * );
static	void	fix_trigcall( NODE_T * );
static	void	fix_indexcall( NODE_T * );
static	void	fix_moliocall( NODE_T * );
static	void	fix_substitute( NODE_T * );
static	void	fix_iocall( NODE_T * );
static	void	fix_iosrc( NODE_T *, int );
static	void	fix_iodst( NODE_T *, int );
static	NODE_T	*fix_alloc( NODE_T * );
static	NODE_T	*fix_direct_attr( NODE_T * );
static	NODE_T	*fix_function_attr( NODE_T * );
static	NODE_T	*fix_struct_field( NODE_T * );

void	fixexpr( NODE_T *expr, int isactual, int cconv )
{
	NODE_T	*npl, *npr, *npf, *np1, *np2, *np3, *npt, *npv;
	SYMREC_T	*sp, *s_array, *stag;
	VALUE_T	*vf, *vn;
	VALUE_T v_flev, v_ptmp, v_1p, v_htype, v_stag;
	char	*pname;
	int	htype;
	int	acc, dim;
	int	l_isactual, l_cconv;

	if( expr ){
		l_isactual = expr->n_sym == SYM_PARM;
		if( expr->n_sym == SYM_CALL ){
			npl = expr->n_left;
			if( !( sp = findsym( npl->n_val.v_value.v_cval ) ) )
				l_cconv = CC_UNDEF;
			else
				l_cconv = sp->s_cconv;
		}else
			l_cconv = cconv;
		if( expr->n_sym == SYM_ALLOCATE ){
			for( npr = expr->n_right; npr; npr = npr->n_right ){
				if( npr->n_sym == SYM_LBRACK )
					break;
			}
			/* remove any parens	*/
			if( npr != expr->n_right )
				expr->n_right = npr;
			for( npr = npr->n_right; npr; npr = npr->n_right )
				npr->n_sym = SYM_LIST;
		}
		if( expr->n_sym == SYM_LBRACK ){
			s_array = findsym( expr->n_left->n_val.v_value.v_cval );
			astk[ astkp ] = s_array;
			astkp++;
			dim = 0;
			for( npr = expr->n_right; npr; npr = npr->n_right ){
				dim++;
				npr->n_val.v_type = T_INT;
				npr->n_val.v_value.v_ival = dim;
			}
		}
		npl = expr->n_left;
		npr = expr->n_right;
		fixexpr( npl, l_isactual, l_cconv );
		fixexpr( npr, l_isactual, l_cconv );
		switch( expr->n_sym ){

		/* operator classes */
		case SYM_ASSIGN :
			if( expr->n_type == T_STRING ){
				np1 = node( SYM_ADDRESS, 0, NULL, npl );
				npf = in2pre( "NAB_strcpy", np1, npr, NULL );
				expr->n_sym = SYM_CALL;
				expr->n_left = npf->n_left;
				expr->n_right = npf->n_right;
			}else if( expr->n_type == T_POINT ){
				npf =in2pre( "NAB_ptcpy", npl, npr, NULL );
				expr->n_sym = SYM_CALL;
				expr->n_left = npf->n_left;
				expr->n_right = npf->n_right;
			}else if( expr->n_type == T_MATRIX ){
				npf =in2pre( "NAB_matcpy", npl, npr, NULL );
				expr->n_sym = SYM_CALL;
				expr->n_left = npf->n_left;
				expr->n_right = npf->n_right;
			}
			break;

		case SYM_PLUS_ASSIGN :
			if( expr->n_type == T_STRING ){
				npf = in2pre( "NAB_strcat", npl, npr, NULL );
				np1 = node( SYM_ADDRESS, 0, NULL, npl );
				npf = in2pre( "NAB_strcpy", np1, npf, NULL );
				expr->n_sym = SYM_CALL;
				expr->n_left = npf->n_left;
				expr->n_right = npf->n_right;
			}else if( expr->n_type == T_POINT ){
				v_ptmp.v_type = T_STRING;
				v_ptmp.v_value.v_cval = CG_gentemp( T_POINT );
				np1 = node( SYM_IDENT, &v_ptmp, NULL, NULL );
				npf = in2pre( "NAB_ptadd", np1, npl, npr );
				npf =in2pre( "NAB_ptcpy", npl, npf, NULL );
				expr->n_sym = SYM_CALL;
				expr->n_left = npf->n_left;
				expr->n_right = npf->n_right;
			}
			break;

		case SYM_MINUS_ASSIGN :
			if( expr->n_type == T_POINT ){
				v_ptmp.v_type = T_STRING;
				v_ptmp.v_value.v_cval = CG_gentemp( T_POINT );
				np1 = node( SYM_IDENT, &v_ptmp, NULL, NULL );
				npf = in2pre( "NAB_ptsub", np1, npl, npr );
				npf =in2pre( "NAB_ptcpy", npl, npf, NULL );
				expr->n_sym = SYM_CALL;
				expr->n_left = npf->n_left;
				expr->n_right = npf->n_right;
			}
			break;

		case SYM_STAR_ASSIGN :
			if( expr->n_type == T_POINT ){
				v_ptmp.v_type = T_STRING;
				v_ptmp.v_value.v_cval = CG_gentemp( T_POINT );
				np1 = node( SYM_IDENT, &v_ptmp, NULL, NULL );
				if( npr->n_type == T_INT ||
					npr->n_type == T_SIZE_T )
				{
					npf = in2pre( "I2R", npr, NULL, NULL );
					npf = in2pre( "NAB_ptscl",
						np1, npf, npl );
				}else
					npf = in2pre( "NAB_ptscl",
						np1, npr, npl );
				npf =in2pre( "NAB_ptcpy", npl, npf, NULL );
				expr->n_sym = SYM_CALL;
				expr->n_left = npf->n_left;
				expr->n_right = npf->n_right;
			}
			break;

		case SYM_SLASH_ASSIGN :
			if( expr->n_type == T_POINT ){
				v_ptmp.v_type = T_STRING;
				v_ptmp.v_value.v_cval = CG_gentemp( T_POINT );
				np1 = node( SYM_IDENT, &v_ptmp, NULL, NULL );
				np2 = in2pre( "RECIP", npr, NULL, NULL );
				npf = in2pre( "NAB_ptscl", np1, np2, npl );
				npf =in2pre( "NAB_ptcpy", npl, npf, NULL );
				expr->n_sym = SYM_CALL;
				expr->n_left = npf->n_left;
				expr->n_right = npf->n_right;
			}
			break;

		case SYM_UPARROW_ASSIGN :
			if( expr->n_type == T_POINT ){
				v_ptmp.v_type = T_STRING;
				v_ptmp.v_value.v_cval = CG_gentemp( T_POINT );
				np1 = node( SYM_IDENT, &v_ptmp, NULL, NULL );
				npf = in2pre( "NAB_ptcrs", np1, npl, npr );
				npf =in2pre( "NAB_ptcpy", npl, npf, NULL );
				expr->n_sym = SYM_CALL;
				expr->n_left = npf->n_left;
				expr->n_right = npf->n_right;
			}
			break;

		case SYM_MODULUS_ASSIGN :
			break;

		case SYM_LESS :
		case SYM_LESS_EQUAL :
		case SYM_EQUAL :
		case SYM_NOT_EQUAL :
		case SYM_GREATER_EQUAL :
		case SYM_GREATER :
			if( npl->n_type == T_STRING ){
				switch( expr->n_sym ){
				case SYM_LESS :
					pname = "LT";
					break;
				case SYM_LESS_EQUAL :
					pname = "LE";
					break;
				case SYM_EQUAL :
					pname = "EQ";
					break;
				case SYM_NOT_EQUAL :
					pname = "NE";
					break;
				case SYM_GREATER_EQUAL :
					pname = "GE";
					break;
				case SYM_GREATER :
					pname = "GT";
					break;
				}
				npf = in2pre( pname, npl, npr, NULL );
				expr->n_sym = SYM_CALL;
				expr->n_left = npf->n_left;
				expr->n_right = npf->n_right;
			}else if( npl->n_type == T_POINT ){
				if( expr->n_sym == SYM_EQUAL )
					npf = in2pre( "PTEQ", npl, npr, NULL );
				else
					npf = in2pre( "PTNE", npl, npr, NULL );
				expr->n_sym = SYM_CALL;
				expr->n_left = npf->n_left;
				expr->n_right = npf->n_right;
			}else if( npl->n_type == T_MATRIX ){
				if( expr->n_sym == SYM_EQUAL )
					npf = in2pre("NAB_mateq",npl,npr,NULL);
				else
					npf = in2pre("NAB_matne",npl,npr,NULL);
				expr->n_sym = SYM_CALL;
				expr->n_left = npf->n_left;
				expr->n_right = npf->n_right;
			}
			break;

		case SYM_PERIOD :
			acc = CG_get_attr_access( expr, npl, npr );
			if( acc == A_DIRECT )
				expr = fix_direct_attr( expr );
			else if( acc == A_FUNCTION )
				expr = fix_function_attr( expr );
			else if( acc == A_STRUCT )
				expr = fix_struct_field( expr );
			
			break;

		/* field - old version:
		case SYM_PERIOD :
			npr->n_sym = SYM_STRING_LIT;
			if( npl->n_type == T_POINT ){
				fname = npr->n_val.v_value.v_cval;
				if( !strcmp( fname, "x" ) )
					pname = "PTX";
				else if( !strcmp( fname, "y" ) )
					pname = "PTY";
				else if( !strcmp( fname, "z" ) )
					pname = "PTZ";
				npf = in2pre( pname, npl, NULL, NULL );
				expr->n_sym = SYM_CALL;
				expr->n_left = npf->n_left;
				expr->n_right = npf->n_right;
			}else if( npl->n_type == T_ATOM ){
				if( expr->n_type == T_INT )
					pname = "NAB_ari";
				else if( expr->n_type == T_FLOAT ){
					pname = "NAB_arf";
					class = C_VAR;
				}else if( expr->n_type == T_STRING )
					pname = "NAB_arc";
				else if( expr->n_type == T_POINT ){
					pname = "NAB_arp";
					class = C_VAR;
				}
				npf = in2pre( pname, npl, npr, NULL );
				if( expr->n_type != T_POINT ){
					expr->n_sym = SYM_INDIRECT;
					expr->n_left = NULL;
					expr->n_right = npf;
				}else{
					*expr = *npf;
					expr->n_type = T_POINT;
					expr->n_class = C_VAR;
					expr->n_kind = K_SCALAR;
				}
			}else if( npl->n_type == T_RESIDUE ){
				if( expr->n_type == T_INT )
					pname = "NAB_rri";
				else if( expr->n_type == T_STRING )
					pname = "NAB_rrc";
				npf = in2pre( pname, npl, npr, NULL );
				expr->n_sym = SYM_INDIRECT;
				expr->n_class = C_EXPR;
				expr->n_left = NULL;
				expr->n_right = npf;
			}else if( npl->n_type == T_MOLECULE ){
				pname = "NAB_mri";
				npf = in2pre( pname, npl, npr, NULL );
				expr->n_sym = SYM_INDIRECT;
				expr->n_class = C_EXPR;
				expr->n_left = NULL;
				expr->n_right = npf;
			}
			break;
		   end of old field version: */

		/* incr op's ++/-- */ 
		case SYM_PLUS_PLUS :
		case SYM_MINUS_MINUS :
			break;

		/* int, float, string (concat), point add */
		case SYM_PLUS :
			if( npl->n_type == T_STRING ){
				npf = in2pre( "NAB_strcat", npl, npr, NULL );
				expr->n_sym = SYM_CALL;
				expr->n_left = npf->n_left;
				expr->n_right = npf->n_right;
			}else if( npl->n_type == T_POINT ){
				v_ptmp.v_type = T_STRING;
				v_ptmp.v_value.v_cval = CG_gentemp( T_POINT );
				np1 = node( SYM_IDENT, &v_ptmp, NULL, NULL );
				npf = in2pre( "NAB_ptadd", np1, npl, npr );
				expr->n_sym = SYM_CALL;
				expr->n_left = npf->n_left;
				expr->n_right = npf->n_right;
			}
			break;

		/* point sub */
		case SYM_MINUS :
			if( npl->n_type == T_POINT ){
				v_ptmp.v_type = T_STRING;
				v_ptmp.v_value.v_cval = CG_gentemp( T_POINT );
				np1 = node( SYM_IDENT, &v_ptmp, NULL, NULL );
				npf = in2pre( "NAB_ptsub", np1, npl, npr );
				expr->n_sym = SYM_CALL;
				expr->n_left = npf->n_left;
				expr->n_right = npf->n_right;
			}
			break;

		case SYM_STAR :
			if( npl->n_type == T_POINT ){
				v_ptmp.v_type = T_STRING;
				v_ptmp.v_value.v_cval = CG_gentemp( T_POINT );
				np1 = node( SYM_IDENT, &v_ptmp, NULL, NULL );
				if( npr->n_type == T_INT ||
					npr->n_type == T_SIZE_T )
				{
					npf = in2pre( "I2R", npr, NULL, NULL );
					npf = in2pre( "NAB_ptscl",
						np1, npf, npl );
				}else
					npf = in2pre( "NAB_ptscl",
						np1, npr, npl );
				expr->n_sym = SYM_CALL;
				expr->n_left = npf->n_left;
				expr->n_right = npf->n_right;
			}else if( npr->n_type == T_POINT ){
				v_ptmp.v_type = T_STRING;
				v_ptmp.v_value.v_cval = CG_gentemp( T_POINT );
				np1 = node( SYM_IDENT, &v_ptmp, NULL, NULL );
				if( npl->n_type == T_INT ||
					npl->n_type == T_SIZE_T )
				{
					npf = in2pre( "I2R", npl, NULL, NULL );
					npf = in2pre( "NAB_ptscl",
						np1, npf, npr );
				}else
					npf = in2pre( "NAB_ptscl",
						np1, npl, npr );
				expr->n_sym = SYM_CALL;
				expr->n_left = npf->n_left;
				expr->n_right = npf->n_right;
			}
			break;

		case SYM_SLASH :
			if( npl->n_type == T_POINT ){
				v_ptmp.v_type = T_STRING;
				v_ptmp.v_value.v_cval = CG_gentemp( T_POINT );
				np1 = node( SYM_IDENT, &v_ptmp, NULL, NULL );
				np2 = in2pre( "RECIP", npr, NULL, NULL );
				npf = in2pre( "NAB_ptscl", np1, np2, npl );
				expr->n_sym = SYM_CALL;
				expr->n_left = npf->n_left;
				expr->n_right = npf->n_right;
			}
			break;

		case SYM_NEGATE :
			if( npr->n_type == T_POINT ){
				v_ptmp.v_type = T_STRING;
				v_ptmp.v_value.v_cval = CG_gentemp( T_POINT );
				np1 = node( SYM_IDENT, &v_ptmp, NULL, NULL );
				v_1p.v_type = T_FLOAT;
				v_1p.v_value.v_fval = 1.0;
				np2 = node( SYM_FLOAT_LIT, &v_1p, NULL, NULL );
				np2 = node( SYM_NEGATE, 0, NULL, np2 );
				npf = in2pre( "NAB_ptscl", np1, np2, npr );
				expr->n_sym = SYM_CALL;
				expr->n_left = npf->n_left;
				expr->n_right = npf->n_right;
			}
			break;

		case SYM_FOREACH :
			if( npl->n_type == T_ATOM ){
				if( npr->n_type == T_MOLECULE ){
					vf = &v_mfirst;
					vn = &v_mnext;
				}else{
					vf = &v_afirst;
					vn = &v_anext;
				}
				np1 = node( SYM_LIST, 0, npl, NULL );
				np1 = node( SYM_LIST, 0, npr, np1 );
				np2 = node( SYM_IDENT, vn, NULL, NULL );
				np2 = node( SYM_CALL, 0, np2, np1 );
				np2 = node( SYM_ASSIGN, 0, npl, np2 );
				np3 = node( SYM_SEMICOLON, 0, np2, NULL );

				/* create: vf( npr, flev )	*/
				np1 = node( SYM_IDENT, &v_NULL, NULL, NULL );
				np2 = node( SYM_ASSIGN, 0, npl, np1 );
				np2 = node( SYM_SEMICOLON, 0, np2, np3 );
			}else if( npl->n_type == T_RESIDUE ){
				vf = &v_rfirst;
				vn = &v_rnext;
				np1 = node( SYM_LIST, 0, npl, NULL );
				np1 = node( SYM_LIST, 0, npr, np1 );
				np2 = node( SYM_IDENT, vn, NULL, NULL );
				np2 = node( SYM_CALL, 0, np2, np1 );
				np2 = node( SYM_ASSIGN, 0, npl, np2 );
				np3 = node( SYM_SEMICOLON, 0, np2, NULL );

				/* create: vf( npr, flev )	*/
				np1 = node( SYM_IDENT, &v_NULL, NULL, NULL );
				np2 = node( SYM_ASSIGN, 0, npl, np1 );
				np2 = node( SYM_SEMICOLON, 0, np2, np3 );
			}else{
				vf = &v_hfirst;
				vn = &v_hnext;

				/* create: npl = vn( npr, flev ) */
				v_flev.v_type = T_STRING;
				v_flev.v_value.v_cval = CG_gentemp( T_CURHASH );
				np1 = node( SYM_IDENT, &v_flev, NULL, NULL );
				np1 = node( SYM_ADDRESS, 0, NULL, np1 );
				np1 = node( SYM_LIST, 0, np1, NULL );
				sp = findsym( npr->n_val.v_value.v_cval );
				if( sp->s_isparm )
					npr = node( SYM_INDIRECT, 0, 0, npr );
				np1 = node( SYM_LIST, 0, npr, np1 );
				np2 = node( SYM_IDENT, vn, NULL, NULL );
				np2 = node( SYM_CALL, 0, np2, np1 );
				np1 = node( SYM_ADDRESS, 0, NULL, npl );
				np1 = in2pre( "NAB_strcpy", np1, np2, NULL );
				np3 = node( SYM_SEMICOLON, 0, np1, NULL );

				/* create: vf( npr, flev )	*/
				np1 = node( SYM_IDENT, &v_flev, NULL, NULL );
				np1 = node( SYM_ADDRESS, 0, NULL, np1 );
				np1 = node( SYM_LIST, 0, np1, NULL );
				np1 = node( SYM_LIST, 0, npr, np1 );
				np2 = node( SYM_IDENT, vf, NULL, NULL );
				np2 = node( SYM_CALL, 0, np2, np1 );
				np2 = node( SYM_SEMICOLON, 0, np2, np3 );
			}
			expr->n_sym = SYM_SEMICOLON;
			expr->n_left = np2->n_left;
			expr->n_right = np2->n_right;
			break;

		case SYM_MATCH :
		case SYM_DONT_MATCH :
			if( npl->n_type == T_ATOM )
				pname = "NAB_aematch";
			else
				pname = "NAB_rematch";
			npf = in2pre( pname, npl, npr, NULL );
			if( expr->n_sym == SYM_DONT_MATCH ){
				npf = node( SYM_NOT, 0, NULL, npf );
				expr->n_sym = SYM_NOT;
			}else
				expr->n_sym = SYM_CALL;
			expr->n_left = npf->n_left;
			expr->n_right = npf->n_right;
			break;

		case SYM_NOT :
			if( npr->n_type == T_POINT ){
				npf = in2pre( "PT_ISTRUE", npr, NULL, NULL );
				npf = node( SYM_NOT, 0, NULL, npf );
				expr->n_sym = SYM_NOT;
				expr->n_left = npf->n_left;
				expr->n_right = npf->n_right;
			}else if( npr->n_type == T_MATRIX ){
				npf = in2pre( "MAT_istrue", npr, NULL, NULL );
				npf = node( SYM_NOT, 0, NULL, npf );
				expr->n_sym = SYM_NOT;
				expr->n_left = npf->n_left;
				expr->n_right = npf->n_right;
			}
			break;

		case SYM_TEST :
			if( npr->n_type == T_POINT ){
				npf = in2pre( "PT_ISTRUE", npr, NULL, NULL );
				expr->n_sym = SYM_CALL;
				expr->n_left = npf->n_left;
				expr->n_right = npf->n_right;
			}else if( npr->n_type == T_MATRIX ){
				npf = in2pre( "MAT_istrue", npr, NULL, NULL );
				expr->n_sym = SYM_CALL;
				expr->n_left = npf->n_left;
				expr->n_right = npf->n_right;
			}
			break;

		case SYM_IN :
			sp = findsym( npr->n_val.v_value.v_cval );
			if( sp->s_isparm )
				npr = node( SYM_INDIRECT, 0, NULL, npr );
			npf = in2pre( "NAB_hin", npr, npl, NULL );
			expr->n_sym = SYM_CALL;
			expr->n_left = npf->n_left;
			expr->n_right = npf->n_right;
			break;

		case SYM_ALLOCATE :
			fix_alloc( expr->n_right );
			break;

		case SYM_DEALLOCATE :
			np1 = node( SYM_LIST, 0, expr->n_right, NULL );
			v_flev.v_type = T_STRING;
			v_flev.v_value.v_cval = "free";
			np2 = node( SYM_IDENT, &v_flev, NULL, NULL );
			expr->n_sym = SYM_CALL;
			expr->n_left = np2;
			expr->n_right = np1;
			break;

		case SYM_DELETE :
			if( expr->n_right->n_sym == SYM_CALL ){
				expr->n_right->n_left->n_val.v_value.v_cval =
					"NAB_hdelete";
			}else{
				np1 = node( SYM_IDENT, &v_NULL, NULL, NULL );
				np1 = node( SYM_LIST, 0, np1, NULL );
				np2 = node( SYM_LIST, 0, expr->n_right, np1 );
				v_flev.v_type = T_STRING;
				v_flev.v_value.v_cval = "NAB_hdelete";
				np1 = node( SYM_IDENT, &v_flev, NULL, NULL );
				expr->n_sym = SYM_CALL;
				expr->n_left = np1;
				expr->n_right = np2;
			}

			break;

		/* indexing */
		case SYM_LBRACK :
			if( npl->n_kind == K_HASHED ){
				pname = "";
				switch( npl->n_type ){
				case T_INT :
					pname = "HRI";
					htype = T_INT;
					break;
				case T_SIZE_T :
					pname = "HRSZ";
					htype = T_SIZE_T;
					break;
				case T_FLOAT :
					pname = "HRF";
					htype = T_FLOAT;
					break;
				case T_STRING :
					pname = "HRC";
					htype = T_STRING;
					break;
				case T_POINT :
					pname = "HRPT";
					htype = T_POINT;
					break;
				case T_MATRIX :
					pname = "HRMAT";
					htype = T_MATRIX;
					break;
				case T_FILE :
					pname = "HRFP";
					htype = T_FILE;
					break;
				case T_ATOM :
					pname = "HRATOM";
					htype = T_ATOM;
					break;
				case T_RESIDUE :
					pname = "HRRES";
					htype = T_RESIDUE;
					break;
				case T_MOLECULE :
					pname = "HRMOL";
					htype = T_MOLECULE;
					break;
				case T_BOUNDS :
					pname = "HRB";
					htype = T_BOUNDS;
					break;
				case T_USER :
					pname = "HRU";
					htype = T_USER;
					break;
				}
				v_htype.v_type = T_INT;
				v_htype.v_value.v_ival = htype;
				if( htype == T_USER ){
					stag = findsym( npl->n_val.v_value.v_cval );
					v_stag.v_type = T_STRING;
					v_stag.v_value.v_cval = stag->s_uname;
					np1 = node( SYM_IDENT, &v_stag, NULL, NULL );
					np1 = node( SYM_STRUCT, NULL, np1, NULL );
					np1 = in2pre( "sizeof", np1, NULL, NULL );
					np2 = node( SYM_INDIRECT, NULL, NULL, NULL );
					np2 = node( SYM_IDENT, &v_stag, NULL, np2 );
					np2 = node( SYM_STRUCT, NULL, np2, NULL );
					np2 = node( SYM_LIST, NULL, np2, NULL );
					np1 = node( SYM_LIST, NULL, np1, np2 );
				}
				npt = node( SYM_INT_LIT, &v_htype, NULL, NULL );
				sp = findsym( npl->n_val.v_value.v_cval );
				if( !sp->s_isparm )
					npl = node( SYM_ADDRESS,NULL,NULL,npl );
				expr->n_left = npl;
				npf = in2pre( pname, npl, npr, npt );
				expr->n_sym = SYM_CALL;
				expr->n_left = npf->n_left;
				if( htype == T_USER ){
					for( np2 = npf->n_right; np2->n_right; np2 = np2->n_right )
						;
					np2->n_right = np1;
				}
				expr->n_right = npf->n_right;
			}
			astkp--;
			break;

		case SYM_INDEX :
			if( expr->n_type == T_INT || expr->n_type == T_SIZE_T ){
				expr->n_left = mk_iexpr( expr );
				expr->n_right = NULL;
			}
			break;

		/* I/O stuff */
		case SYM_CALL :
			if( is_getline( expr ) )
				fix_getline( expr );
			else if( is_trig( expr ) )
				fix_trigcall( expr );
			else if( is_index( expr ) )
				fix_indexcall( expr );
			else if( is_molio( expr ) )
				fix_moliocall( expr );
			else if( is_substitute(expr) )
				fix_substitute( expr );
			else
				fix_iocall( expr );
			break;

		case SYM_PARM :
			fixparm( expr, cconv );
			break;

		/* from on, ignored */
		case SYM_ADDRESS :
		case SYM_INDIRECT :
		case SYM_AND :
		case SYM_OR :
		case SYM_MODULUS :
			break;

		case SYM_ATSIGN :
			npv = in2pre( "NAB_ptdot", npl, npr, NULL );
			expr->n_sym = SYM_CALL;
			expr->n_left = npv->n_left;
			expr->n_right = npv->n_right;
			break;

		case SYM_UPARROW :
			v_ptmp.v_type = T_STRING;
			v_ptmp.v_value.v_cval = CG_gentemp( T_POINT );
			np1 = node( SYM_IDENT, &v_ptmp, NULL, NULL );
			npv = in2pre( "NAB_ptcrs", np1, npl, npr );
			expr->n_sym = SYM_CALL;
			expr->n_left = npv->n_left;
			expr->n_right = npv->n_right;
			break;

		case SYM_IDENT :
			fixid( expr, isactual );
			break;

		case SYM_ATTRIBUTE :
			break;

		case SYM_INT_LIT :
		case SYM_FLOAT_LIT :
		case SYM_STRING_LIT :
			break;

		case SYM_BREAK :
		case SYM_ELSE :
		case SYM_FOR :
		case SYM_IF :
		case SYM_RETURN :
		case SYM_WHILE :
			break;

		case SYM_DECL :
		case SYM_TYPE :
		case SYM_INT :
		case SYM_FLOAT :
		case SYM_STRING :
		case SYM_POINT :
		case SYM_MATRIX :
		case SYM_FILE :
		case SYM_ATOM :
		case SYM_RESIDUE :
		case SYM_MOLECULE :
		case SYM_BOUNDS :
		case SYM_HASHED :
			break;

		case SYM_LPAREN :
		case SYM_RPAREN :
		case SYM_RBRACK :
		case SYM_LBRACE :
		case SYM_RBRACE :
		case SYM_COMMA :
		case SYM_SEMICOLON :
		case SYM_LIST :
			break;

		case SYM_STMTLIST :
			break;

		case SYM_ERROR :
			break;
		}
	}
}

static	NODE_T	*mk_iexpr( NODE_T *expr )
{
	int	dim, ter;
	NODE_T	*np1;

	dim = expr->n_val.v_value.v_ival;
	ter = expr->n_right == NULL;
	if( ter ){
		if( dim == 1 ){
			np1 = node( SYM_INT_LIT, &v_1, 0, 0 );
			np1 = node( SYM_MINUS, 0, expr->n_left, np1 );
		}else{
			np1 = node( SYM_INT_LIT, &v_1, 0, 0 );
			np1 = node( SYM_MINUS, 0, expr->n_left, np1 );
			np1 = node( SYM_STAR, 0,
				astk[astkp-1]->s_parts[dim-2], np1 );
		}
	}else if( dim == 1 ){
		np1 = node( SYM_INT_LIT, &v_1, 0, 0 );
		np1 = node( SYM_MINUS, 0, expr->n_left, np1 );
		np1 = node( SYM_PLUS, 0, np1, expr->n_right->n_left );
		expr->n_right->n_left = NULL;
	}else{
		np1 = node( SYM_INT_LIT, &v_1, 0, 0 );
		np1 = node( SYM_MINUS, 0, expr->n_left, np1 );
		np1 = node( SYM_PLUS, 0, np1, expr->n_right->n_left );
		expr->n_right->n_left = NULL;
		np1 = node( SYM_STAR, 0, astk[astkp-1]->s_parts[dim-2], np1 );
	}
	return( np1 );
}

static	int	is_getline( NODE_T *expr )
{
	NODE_T	*npl;

	npl = expr->n_left;
	if( !strcmp( npl->n_val.v_value.v_cval, "getline" ) )
		return( 1 );
	else
		return( 0 );
}

static	void	fix_getline( NODE_T *expr )
{
	NODE_T	*npl;

	npl = expr->n_left;
	free( npl->n_val.v_value.v_cval );
	npl->n_val.v_value.v_cval = strdup( "NAB_getline" );
}

static	int	is_trig( NODE_T *expr )
{
	NODE_T	*npl;
	char	*fname;

	npl = expr->n_left;
	fname = npl->n_val.v_value.v_cval;
	if(	!strcmp( fname, "acos" ) ||
		!strcmp( fname, "asin" ) ||
		!strcmp( fname, "atan" ) ||
		!strcmp( fname, "atan2" ) ||
		!strcmp( fname, "cos" ) ||
		!strcmp( fname, "sin" ) ||
		!strcmp( fname, "tan" ) )
		return( 1 );
	else
		return( 0 );
}

static	void	fix_trigcall( NODE_T *expr )
{
	NODE_T	*npl;
	char	*fnp;

	npl = expr->n_left;
	for( fnp = npl->n_val.v_value.v_cval; *fnp; fnp++ )
		*fnp = islower( *fnp ) ? toupper( *fnp ) : *fnp;
}

static	int	is_index( NODE_T *expr )
{
	NODE_T	*npl;
	char	*fname;

	npl = expr->n_left;
	fname = npl->n_val.v_value.v_cval;
	if( !strcmp( fname, "index" ) )
		return( 1 );
	else
		return( 0 );
}

static	void	fix_indexcall( NODE_T * expr )
{
	NODE_T	*npl;

	npl = expr->n_left;
	npl->n_val.v_value.v_cval = "NAB_index";
}

static	int	is_molio( NODE_T *expr )
{
	NODE_T	*npl;
	char	*fname;

	npl = expr->n_left;
	fname = npl->n_val.v_value.v_cval;
	if( !strcmp( fname, "getpdb" ) )
		return( 1 );
	else if( !strcmp( fname, "putpdb" ) )
		return( 1 );
	else
		return( 0 );
}

static	void	fix_moliocall( NODE_T *expr )
{
	NODE_T	*npl, *npr, *nplast;
	char	*fname;
	int	nparm;

	NODE_T	*np1, *np2;
	npl = expr->n_left;
	npr = expr->n_right;
	for( nplast = NULL, nparm = 0; npr; npr = npr->n_right ){
		nparm++;
		nplast = npr;
	}
	fname = npl->n_val.v_value.v_cval;
	if( !strcmp( fname, "getpdb" ) ){
		if( nparm == 2 )
			return;
	}else if( nparm == 3 )
		return;
	np1 = node( SYM_IDENT, &v_NULL, NULL, NULL );
	np2 = node( SYM_PARM, 0, NULL, np1 );
	np1 = node( SYM_LIST, 0, np2, NULL );
	nplast->n_right = np1;
}

static	int	is_substitute(NODE_T *expr)
{
	NODE_T	*npl;
	char	*fname;

	npl = expr->n_left;
	fname = npl->n_val.v_value.v_cval;
	if( !strcmp( fname, "sub" ) )
		return( 1 );
	else if( !strcmp( fname, "gsub" ) )
		return( 1 );
	else
		return( 0 );
}

static	void	fix_substitute(NODE_T *expr)
{
	NODE_T	*npl, *npr;
	NODE_T	*np1, *np2;

	npl = expr->n_left;
	npr = expr->n_right;

	// insert 1 or 0 as 1st parm for gsub() or sub()
	if(!strcmp(npl->n_val.v_value.v_cval, "gsub"))
		np1 = node(SYM_INT_LIT, &v_1, NULL, NULL);
	else
		np1 = node(SYM_INT_LIT, &v_0, NULL, NULL);
	np2 = node(SYM_PARM, 0, NULL, np1);
	expr->n_right = node(SYM_LIST, 0, np2, npr);
	npl->n_val.v_value.v_cval = "NAB_gsub";
}

static	void	fixparm( NODE_T * expr, int cconv )
{
	NODE_T	*npt, *npp, *npv;
	SYMREC_T	*sp;
	char	*pname;
	VALUE_T	val;

	if( cconv == CC_UNDEF )
		return;
	npp = expr->n_right;
	if( npp->n_class == C_NULL )
		return;
	if( cconv == CC_NAB || cconv == CC_FORTRAN  ){
		if( npp->n_kind == K_ARRAY || npp->n_kind == K_DARRAY )
			return;
		else if( npp->n_kind == K_HASHED ){
			sp = findsym( npp->n_val.v_value.v_cval );
			if( !sp->s_isparm ){
				npv = node(SYM_IDENT,&npp->n_val,
					NULL,NULL);
				npv->n_type = expr->n_type;
				npv->n_class = C_VAR;
				npv->n_kind = K_SCALAR;
				npp->n_sym = SYM_ADDRESS;
				npp->n_right = npv;
			}
			return;
		}else if( npp->n_class == C_FUNC )
			return;
		else if( npp->n_class == C_LIT ||
			npp->n_class == C_DEFINE ||
			npp->n_class == C_EXPR ){
			if( npp->n_type == T_INT ){
				pname = "ITEMP";
				val.v_value.v_cval = CG_gentemp( T_INT );
			}else if( npp->n_type == T_SIZE_T ){
				pname = "SZTEMP";
				val.v_value.v_cval = CG_gentemp( T_SIZE_T );
			}else if( npp->n_type == T_FLOAT ){
				pname = "FTEMP";
				val.v_value.v_cval = CG_gentemp( T_FLOAT );
			}else if( npp->n_type == T_STRING ){
				pname = "STEMP";
				val.v_value.v_cval = CG_gentemp( T_STRING );
			}else if( npp->n_type == T_FILE ){
				pname = "FPTEMP";
				val.v_value.v_cval = CG_gentemp( T_FILE );
			}
			val.v_type = T_STRING;
			npt = node( SYM_IDENT, &val, NULL, NULL );
			npv = in2pre( pname, npt, npp, NULL );
			if( cconv == CC_FORTRAN && npp->n_type == T_STRING ){
				expr->n_sym = SYM_INDIRECT;
				expr->n_left = NULL;
				expr->n_right = npv;
			}else{
				expr->n_sym = SYM_CALL;
				expr->n_left = npv->n_left;
				expr->n_right = npv->n_right;
			}
		}else if( npp->n_class == C_VAR ){
			if( npp->n_type == T_POINT || npp->n_type == T_MATRIX )
				return;
			if( npp->n_sym == SYM_IDENT ){
				sp = findsym( npp->n_val.v_value.v_cval );
				if( sp->s_isparm ){
					if( cconv == CC_FORTRAN && npp->n_type == T_STRING ){
						npv = node(SYM_IDENT,&npp->n_val,NULL,NULL);
						npv->n_type = expr->n_type;
						npv->n_class = C_VAR;
						npv->n_kind = K_SCALAR;
						npp->n_sym = SYM_INDIRECT;
						npp->n_right = npv;
					}
					return;
				}else if( cconv != CC_FORTRAN || npp->n_type != T_STRING ){
						npv = node(SYM_IDENT,&npp->n_val,
							NULL,NULL);
						npv->n_type = expr->n_type;
						npv->n_class = C_VAR;
						npv->n_kind = K_SCALAR;
						npp->n_sym = SYM_ADDRESS;
						npp->n_right = npv;
				}
			}else{
				npv = node( npp->n_sym, NULL,
					npp->n_left, npp->n_right );
				npv->n_type = expr->n_type;
				npv->n_class = C_VAR;
				npv->n_kind = K_SCALAR;
				npp->n_sym = SYM_ADDRESS;
				npp->n_left = NULL;
				npp->n_right = npv;
			}
		}
	}else if( cconv == CC_CC ){
		if( npp->n_kind == K_ARRAY || npp->n_kind == K_HASHED )
			return;
		else if( npp->n_class == C_FUNC )
			return;
		else if( npp->n_class == C_LIT ||
			npp->n_class == C_DEFINE ||
			npp->n_class == C_EXPR ){
			return;
		}else if( npp->n_class == C_VAR ){
			if( npp->n_type == T_POINT || npp->n_type == T_MATRIX )
				return;
			if( npp->n_type == T_USER ){
				npv = node(npp->n_sym,&npp->n_val,npp->n_left,npp->n_right);
				npv->n_type = npp->n_type;
				npv->n_class = npp->n_class;
				npv->n_kind = npp->n_kind;
				npp->n_sym = SYM_ADDRESS;
				npp->n_left = NULL;
				npp->n_right = npv;
				return;
			}
			if( npp->n_sym == SYM_IDENT ){
				sp = findsym( npp->n_val.v_value.v_cval );
				if( !sp->s_isparm )
					return;
				else{
					npv = node(SYM_IDENT,&npp->n_val,NULL,NULL);
					npv->n_type = expr->n_type;
					npv->n_class = C_VAR;
					npv->n_kind = K_SCALAR;
					npp->n_sym = SYM_INDIRECT;
					npp->n_right = npv;
				}
			}
		}
	}else if( cconv == CC_IO )
		return;	/* fixed up in fix_iocall() */
}

static	void	fixid( NODE_T * expr, int isactual )
{
	SYMREC_T	*sp;
	NODE_T		*np;

	if( isactual )
		return;	/* fixed by the PARM operator */
	else if( expr->n_kind == K_ARRAY || expr->n_kind == K_HASHED )
		return;
	else if( expr->n_class != C_VAR )
		return;
	else if( expr->n_type == T_POINT || expr->n_type == T_MATRIX )
		return;
	sp = findsym( expr->n_val.v_value.v_cval );
	if( !sp->s_isparm )
		return;
	np = node( SYM_IDENT, &expr->n_val, NULL, NULL );
	np->n_type = expr->n_type;
	np->n_class = C_VAR;
	np->n_kind = K_SCALAR;
	expr->n_left = NULL;
	expr->n_right = np;
	expr->n_sym = SYM_INDIRECT;
}

static	NODE_T	*in2pre( char pname[], NODE_T *p1, NODE_T *p2, NODE_T *p3 )
{
	NODE_T	*npl, *npc;
	VALUE_T	val;

	npl = p3 ? node( SYM_LIST, 0, p3, NULL ) : NULL;
	npl = p2 ? node( SYM_LIST, 0, p2, npl ) : NULL;
	npl = node( SYM_LIST, 0, p1, npl );
	
	val.v_value.v_cval = pname;
	npc = node( SYM_IDENT, &val, NULL, NULL );
	npc = node( SYM_CALL, 0, npc, npl );
	return( npc );
}

static	void	fix_iocall( NODE_T *expr )
{
	int	ch_dir;	/* ch parm direction src <-> dst w/this parm	*/
	int	sprf, file_ok;
	NODE_T	*npl, *npr, *npp;
	NODE_T	*np1, *np2, *npp1, *npp2;
	NODE_T	*npl1, *npl2, *npi;
	NODE_T	*npc_snp, *npc_scp, *npc_asrt;
	VALUE_T	val;
	char	*pname;
	int	pcnt;

	sprf = FALSE;
	file_ok = FALSE;
	npl = expr->n_left;
	npr = expr->n_right;
	pname = npl->n_val.v_value.v_cval;
	if( !strcmp( "scanf", pname ) ){
		ch_dir = 1;
	}else if( !strcmp( "fscanf", pname ) || !strcmp( "sscanf", pname ) ){
		file_ok = TRUE;
		ch_dir = 2;
	}else if( !strcmp( "sprintf", pname ) ){
		ch_dir = 2;
		sprf = TRUE;
	}else if( !strcmp( "printf", pname ) || !strcmp( "fprintf", pname ) ){
		/* act like CC funcs, no action required */
		return;
	}else
		return;	/* Not I/O call */

	if( sprf ){
		/* results to NAB_rsbuf via snprintf() inside of 	*/
		/* assert( snprintf(...) < NAB_RSBUF_SIZE ) which will	*/
		/* abort if buf too small.  If OK, copy to the temp	*/	
		val.v_type = T_STRING;
		val.v_value.v_ival = NAB_RSBUF_SIZE;
		np2 = node( SYM_INT_LIT, &val, NULL, NULL );
		np2->n_type = T_INT;
		np2->n_class = C_LIT;
		np2->n_kind = K_SCALAR;
		npp2 = node( SYM_PARM, 0, NULL, np2 );
		npp2->n_type = T_STRING;
		npp2->n_class = C_DEFINE;
		npp2->n_kind = K_SCALAR;
		npl2 = node( SYM_LIST, 0, npp2, npr );

		val.v_type = T_STRING;
		val.v_value.v_cval = "NAB_rsbuf";
		np1 = node( SYM_IDENT, &val, NULL, NULL );
		np1->n_type = T_STRING;
		np1->n_class = C_VAR;
		np1->n_kind = K_SCALAR;
		npp1 = node( SYM_PARM, 0, NULL, np1 );
		npp1->n_type = T_STRING;
		npp1->n_class = C_VAR;
		npp1->n_kind = K_SCALAR;
		npr = npl1 = node( SYM_LIST, 0, npp1, npl2 );

		for( pcnt = 0, npp = npr; npp; npp = npp->n_right ){
			pcnt++;
			if( pcnt > ch_dir )
				fix_iosrc( npp, FALSE );
		}

		val.v_type = T_STRING;
		val.v_value.v_cval = "snprintf";
		npi = node( SYM_IDENT, &val, NULL, NULL );
		npc_snp = node( SYM_CALL, 0, npi, npl1 );

		val.v_type = T_STRING;
		val.v_value.v_ival = NAB_RSBUF_SIZE;
		np2 = node( SYM_INT_LIT, &val, NULL, NULL );
		np2->n_type = T_INT;
		np2->n_class = C_LIT;
		np2->n_kind = K_SCALAR;
		np1 = node( SYM_LESS, 0, npc_snp, np2 );
		np1->n_type = T_INT;
		np1->n_class = C_EXPR;
		np1->n_kind = K_SCALAR;
		npp1 = node( SYM_PARM, 0, NULL, np1 );
		npp1->n_type = T_INT;
		npp1->n_class = C_EXPR;
		npp1->n_kind = K_SCALAR;
		npl2 = node( SYM_LIST, 0, npp1, NULL );
		val.v_type = T_STRING;
		val.v_value.v_cval = "assert";
		npi = node( SYM_IDENT, &val, NULL, NULL ); 
		npc_asrt = node( SYM_CALL, 0, npi, npl2 );

		val.v_type = T_STRING;
		val.v_value.v_cval = "NAB_rsbuf";
		np2 = node( SYM_IDENT, &val, NULL, NULL );
		np2->n_type = T_STRING;
		np2->n_class = C_VAR;
		np2->n_kind = K_SCALAR;
		npp2 = node( SYM_PARM, 0, NULL, np2 );
		npp2->n_type = T_STRING;
		npp2->n_class = C_VAR;
		npp2->n_kind = K_SCALAR;
		npl2 = node( SYM_LIST, 0, npp2, NULL );

		val.v_type = T_STRING;
		val.v_value.v_cval = CG_gentemp( T_STRING );
		np1 = node( SYM_IDENT, &val, NULL, NULL );
		np1->n_type = T_STRING;
		np1->n_class = C_VAR;
		np1->n_kind = K_SCALAR;
		np1 = node( SYM_ADDRESS, 0, NULL, np1 );
		npp1 = node( SYM_PARM, 0, NULL, np1 );
		npp1->n_type = T_STRING;
		npp1->n_class = C_VAR;
		npp1->n_kind = K_SCALAR;
		npl1 = node( SYM_LIST, 0, npp1, npl2 );

		val.v_type = T_STRING;
		val.v_value.v_cval = "NAB_strcpy";
		npi = node( SYM_IDENT, &val, NULL, NULL ); 
		npc_scp = node( SYM_CALL, 0, npi, npl1 );

		npp2 = node( SYM_PARM, 0, NULL, npc_scp );
		npl2 = node( SYM_LIST, 0, npp2, NULL );
		npp1 = node( SYM_PARM, 0, NULL, npc_asrt );
		npl1 = node( SYM_LIST, 0, npp1, npl2 );

		val.v_value.v_cval = "SPRINTF";
		npi = node( SYM_IDENT, &val, 0, 0 );
		expr->n_left = npi;
		expr->n_right = npl1;
	}else{
		for( pcnt = 0, npp = npr; npp; npp = npp->n_right ){
			pcnt++;
			if( pcnt <= ch_dir )
				fix_iosrc( npp, file_ok );
			else
				fix_iodst( npp, 0 );
		}
	}
}

static	void	fix_iodst( NODE_T *npp, int istemp )
{
	NODE_T	*npl, *npa, *np1, *np2;
	SYMREC_T	*sp;
	VALUE_T	rsv;
	
	npl = npp->n_left;
	if( npl->n_type == T_INT || npl->n_type == T_SIZE_T ||
		npl->n_type == T_FLOAT )
	{
		if( npl->n_sym == SYM_PARM ){
			npl = npl->n_right;
		}
		if( npl->n_sym != SYM_IDENT || istemp ){
			npa = node( SYM_ADDRESS, 0, NULL, npl );
			npp->n_left = npa;
		}else{
			sp = findsym( npl->n_val.v_value.v_cval );
			if( !sp->s_isparm ){
				npa = node( SYM_ADDRESS, 0, NULL, npl );
				npp->n_left = npa;
			}
		}
	}else if( npl->n_type == T_STRING ){
		if( npl->n_sym == SYM_PARM ){
			npl = npl->n_right;
		}
		if( npl->n_sym != SYM_IDENT || istemp )
			np1 = node( SYM_ADDRESS, 0, NULL, npl ); /* addr of parm */
		else{
			sp = findsym( npl->n_val.v_value.v_cval );
			if( !sp->s_isparm )
				np1 = node( SYM_ADDRESS, 0, NULL, npl );
			else
				np1 = npl;
		}
		np1 = node( SYM_LIST, 0, np1, NULL ); 
		rsv.v_type = T_STRING;
		rsv.v_value.v_cval = "NAB_readstring";
		np2 = node( SYM_IDENT, &rsv, NULL, NULL );
		np2 = node( SYM_CALL, 0, np2, np1 );
		npp->n_left = np2;
	}else{
		errormsg( FALSE, "can only printf() int/float/string.\n" );
	}
}

static	void	fix_iosrc( NODE_T *npp, int file_ok )
{
	NODE_T	*npl, *npi;
	SYMREC_T	*sp;
	
	npl = npp->n_left;
	if(npl->n_type==T_INT || npl->n_type == T_SIZE_T ||
		npl->n_type==T_FLOAT || npl->n_type==T_STRING)
	{
		if( npl->n_sym == SYM_PARM ){
			npl = npl->n_right;
		}
		if( npl->n_sym == SYM_IDENT ){
			sp = findsym( npl->n_val.v_value.v_cval );
			if( sp->s_isparm ){
				npi = node( SYM_INDIRECT, 0, NULL, npl );
				npp->n_left = npi;
			}
		}
	}else if( npl->n_type == T_FILE && file_ok ){
		if( npl->n_sym == SYM_PARM ){
			npl = npl->n_right;
		}
		if( npl->n_sym == SYM_IDENT ){
			sp = findsym( npl->n_val.v_value.v_cval );
			if( sp->s_isparm ){
				npi = node( SYM_INDIRECT, 0, NULL, npl );
				npp->n_left = npi;
			}
		}
	}else{
		errormsg( FALSE, "can only scanf() int/float/string.\n" );
	}
}

static	NODE_T	*fix_alloc( NODE_T *expr )
{
	NODE_T	*npv, *npi;
	NODE_T	*npa, *npsm, *npsm1;
	NODE_T	*npt, *npty;
	NODE_T	*npsz, *npanm, *npfnm;
	SYMREC_T	*s_array;
	int	nd;
	VALUE_T	val, v_type;
	
	while( expr->n_sym != SYM_LBRACK )
		expr = expr->n_right;
	npv = expr->n_left;
	npi = expr->n_right;
	npsm = npt = NULL;
	s_array = findsym( npv->n_val.v_value.v_cval );

	for( nd = 0; nd < s_array->s_pcount; nd++ ){
		npa = node( SYM_ASSIGN, 0, s_array->s_parts[nd], npi->n_left );
		if( npsm == NULL ){
			npsm = node( SYM_SEMICOLON, 0, npa, NULL );
			npsm1 = npsm;
		}else{
			npsm1->n_right = node( SYM_SEMICOLON, 0, npa, NULL );
			npsm1 = npsm1->n_right;
		}
		if( npt == NULL )
			npt = s_array->s_parts[ nd ];
		else
			npt = node( SYM_STAR, 0, npt, s_array->s_parts[ nd ] );
		npi = npi->n_right;
	}

	npty = NULL;
	if( s_array->s_type == T_USER ){
		v_type.v_type = T_STRING;
		v_type.v_value.v_cval = s_array->s_uname;
		npty = node( SYM_IDENT, &v_type, NULL, NULL );
		npty = node( SYM_STRUCT, NULL, npty, NULL );
	}else if(
		s_array->s_type != T_INT &&
		s_array->s_type != T_SIZE_T &&
		s_array->s_type != T_FLOAT &&
		s_array->s_type != T_POINT &&
		s_array->s_type != T_MATRIX &&
		s_array->s_type != T_USER )
	{
		npty = node( SYM_INDIRECT, 0, NULL, NULL );
	}
	v_type.v_type = T_INT;
	v_type.v_value.v_ival = s_array->s_type;
	npty = node( SYM_TYPE, &v_type, NULL, npty );
	npty->n_type = s_array->s_type;
	if( npty->n_type == T_STRING || s_array->s_init ){
		npsz = in2pre( "sizeof", npty, NULL, NULL );
		npsz = node( SYM_LIST, 0, npsz, NULL );
		npsz = node( SYM_LIST, 0, npt, npsz );
		npsz = in2pre( "calloc", npsz, NULL, NULL );
	}else{
		npsz = in2pre( "sizeof", npty, NULL, NULL );
		npsz = node( SYM_STAR, 0, npt, npsz );
		npsz = in2pre( "malloc", npsz, NULL, NULL );
	}

	npty = node( SYM_INDIRECT, 0, NULL, NULL );
	if( s_array->s_type == T_USER ){
		v_type.v_type = T_STRING;
		v_type.v_value.v_cval = s_array->s_uname;
		npty = node( SYM_IDENT, &v_type, NULL, npty );
		npty = node( SYM_STRUCT, NULL, npty, NULL );
	}else if(	
		s_array->s_type != T_INT &&
		s_array->s_type != T_SIZE_T &&
		s_array->s_type != T_FLOAT &&
		s_array->s_type != T_POINT &&
		s_array->s_type != T_MATRIX &&
		s_array->s_type != T_USER )
	{
		npty = node( SYM_INDIRECT, 0, NULL, npty );
	}
	v_type.v_type = T_INT;
	v_type.v_value.v_ival = s_array->s_type;
	npty = node( SYM_TYPE, &v_type, NULL, npty );
	npty->n_type = s_array->s_type;
	npty = in2pre( "", npty, NULL, NULL );
	npsz = node( SYM_STMTLIST, 0, npty, npsz );
	npsz = node( SYM_ASSIGN, 0, npv, npsz );
	expr->n_sym = SYM_STMTLIST;
	expr->n_left = npsm;
	val.v_type = T_STRING;
	val.v_value.v_cval = npv->n_val.v_value.v_cval;
	npanm = node( SYM_STRING_LIT, &val, NULL, NULL );
	val.v_value.v_cval = c_funcname;
	npfnm = node( SYM_STRING_LIT, &val, NULL, NULL );
	expr->n_right = in2pre( "DA_ALLOC", npsz, npfnm, npanm );
	return( NULL );
}

static	NODE_T	*fix_direct_attr( NODE_T *expr )
{
	NODE_T	*npl, *npr;
	char	*aname, *fname;
	char	field[ 40 ];
	NODE_T	*npf, *npv, *npp;
	VALUE_T	val;

	npl = expr->n_left;
	npr = expr->n_right;
	aname = npr->n_val.v_value.v_cval;
	if( npl->n_type == T_POINT ){
		if( !strcmp( aname, "x" ) )
			fname = "PTX";
		else if( !strcmp( aname, "y" ) )
			fname = "PTY";
		else if( !strcmp( aname, "z" ) )
			fname = "PTZ";
		else
			return( expr );
		npf = in2pre( fname, npl, NULL, NULL );
		expr->n_sym = SYM_CALL;
		expr->n_left = npf->n_left;
		expr->n_right = npf->n_right;
		return( expr );
	}

	fname = NULL; 
	if( npl->n_type == T_ATOM ){
		if( !strcmp( aname, "strandname" ) ||
			!strcmp( aname, "strandnum" ) )

			sprintf( field, "a_residue->r_strand->s_%s", aname );

		else if( !strcmp( aname, "resid" ) ||
			!strcmp( aname, "resname" ) ||
			!strcmp( aname, "resnum" ) ||
			!strcmp( aname, "tresnum" ) )

			sprintf( field, "a_residue->r_%s", aname );

		else if( !strcmp( aname, "x" ) ){
			strcpy( field, "a_pos" );
			fname = "PTX";
		}else if( !strcmp( aname, "y" ) ){
			strcpy( field, "a_pos" );
			fname = "PTY";
		}else if( !strcmp( aname, "z" ) ){
			strcpy( field, "a_pos" );
			fname = "PTZ";
		}else
			sprintf( field, "a_%s", aname );

	}else if( npl->n_type == T_RESIDUE ){
		if( !strcmp( aname, "strandname" ) ||
			!strcmp( aname, "strandnum" ) )

			sprintf( field, "r_strand->s_%s", aname );

		else
			sprintf( field, "r_%s", aname );

	}else if( npl->n_type == T_MOLECULE )
		sprintf( field, "m_%s", aname );

	val.v_type = T_STRING;
	val.v_value.v_cval = field;
	npv = node( SYM_ATTRIBUTE, &val, NULL, NULL );
	if( fname == NULL ){
		expr->n_sym = SYM_POINTS_TO;
		expr->n_right = npv;
	}else{
		npp = node( SYM_POINTS_TO, NULL, npl, npv );
		npf = in2pre( fname, npp, NULL, NULL );
		expr->n_sym = SYM_CALL;
		expr->n_left = npf->n_left;
		expr->n_right = npf->n_right;
	}

	return( expr );
}

static	NODE_T	*fix_function_attr( NODE_T *expr )
{
	NODE_T	*npl, *npr;
	char	*aname, *fname;
	NODE_T	*npf;

	npl = expr->n_left;
	npr = expr->n_right;
	aname = npr->n_val.v_value.v_cval;
	npr->n_sym = SYM_STRING_LIT;
	if( npl->n_type == T_ATOM ){
		if( expr->n_type == T_INT )
			fname = "NAB_ari";
		else if( expr->n_type == T_STRING )
			fname = "NAB_arc";
		npf = in2pre( fname, npl, npr, NULL );
		expr->n_sym = SYM_INDIRECT;
		expr->n_left = NULL;
		expr->n_right = npf;
	}else if( npl->n_type == T_RESIDUE ){
		if( expr->n_type == T_INT )
			fname = "NAB_rri";
		npf = in2pre( fname, npl, npr, NULL );
		expr->n_sym = SYM_INDIRECT;
		expr->n_class = C_EXPR;
		expr->n_left = NULL;
		expr->n_right = npf;
	}else if( npl->n_type == T_MOLECULE ){
		fname = "NAB_mri";
		npf = in2pre( fname, npl, npr, NULL );
		expr->n_sym = SYM_INDIRECT;
		expr->n_class = C_EXPR;
		expr->n_left = NULL;
		expr->n_right = npf;
	}
	return( expr );
}

static	NODE_T	*fix_struct_field( NODE_T *expr )
{
	NODE_T	*npl;
	SYMREC_T	*st;

	npl = expr->n_left;
	if( npl->n_sym == SYM_IDENT ){
		st = findsym( npl->n_val.v_value.v_cval );
		if( st->s_isparm ){
			expr->n_left = node(SYM_LPAREN, NULL, npl, NULL );
		}
	}
	return( expr );
}
