#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "nab.h"
#include "nabgrm.tab.h"
#include "errormsg.h"

extern	int	cg_lineno;

NODE_T	*node( int sym, VALUE_T *vp, NODE_T *left, NODE_T *right )
{
	NODE_T	*np;
	int	clen;
	char	*cp;
	char	msg[ 256 ];

	if( ( np = ( NODE_T * )malloc( sizeof( NODE_T ) ) ) == NULL ){
		sprintf( msg, "NODE %d", sym );
		errormsg_s( TRUE, E_NOMEM_FOR_S, msg );
	}
	np->n_lineno = cg_lineno;
	np->n_sym = sym;
	np->n_type = T_UNDEF;
	np->n_class = C_UNDEF;
	np->n_kind = K_UNDEF;
/*
	np->n_val.v_type = T_INT;
	np->n_val.v_value.v_ival = ( int )vp;
*/
	np->n_val.v_type = T_UNDEF;
	np->n_val.v_value.v_cval = NULL;
	np->n_left = left;
	np->n_right = right;
	if( sym == SYM_IDENT || sym == SYM_ATTRIBUTE ){
		clen = strlen( vp->v_value.v_cval );
		if( ( cp = ( char * )malloc( (clen + 1) * sizeof(char) ) ) == NULL ){
			sprintf( msg, "%s %s\n",
				sym == SYM_IDENT ? "IDENT" : "ATTRIBUTE",
				vp->v_value.v_cval );
			errormsg_s( TRUE, E_NOMEM_FOR_S, msg );
		}
		strcpy( cp, vp->v_value.v_cval );
		np->n_val.v_type = T_STRING;
		np->n_val.v_value.v_cval = cp;
	}else if( sym == SYM_STRING_LIT ){
		clen = strlen( vp->v_value.v_cval );
		if( ( cp = ( char * )malloc( (clen + 1)*sizeof(char) ) ) == NULL ){
			sprintf( msg, "STRING %s", vp->v_value.v_cval );
			errormsg_s( TRUE, E_NOMEM_FOR_S, msg );
		}
		strcpy( cp, vp->v_value.v_cval );
		/*strext2int( cp, vp->v_value.v_cval );*/
		np->n_type = T_STRING;
		np->n_class = C_LIT;
		np->n_kind = K_SCALAR;
		np->n_val.v_type = T_STRING;
		np->n_val.v_value.v_cval = cp;
	}else if( sym == SYM_INT_LIT ){
		np->n_type = T_INT;
		np->n_class = C_LIT;
		np->n_kind = K_SCALAR;
		np->n_val.v_type = T_INT;
		np->n_val.v_value.v_ival = vp->v_value.v_ival;
	}else if( sym == SYM_FLOAT_LIT ){
		np->n_type = T_FLOAT;
		np->n_class = C_LIT;
		np->n_kind = K_SCALAR;
		np->n_val.v_type = T_FLOAT;
		np->n_val.v_value.v_fval = vp->v_value.v_fval;
	}else if( sym == SYM_TYPE ){
		np->n_val.v_type = vp->v_type;
		np->n_val.v_value.v_ival = vp->v_value.v_ival;
	}
	return( np );
}

NODE_T	*copynode( NODE_T *np )
{
	NODE_T	*npl, *npr, *np1;

	if( np ){
		npl = copynode( np->n_left );
		npr = copynode( np->n_right );
/*
		if(	np->n_sym == SYM_IDENT ||
			np->n_sym == SYM_ATTRIBUTE ||
			np->n_sym == SYM_INT_LIT ||
			np->n_sym == SYM_FLOAT_LIT ||
			np->n_sym == SYM_STRING_LIT  )
			np1 = node( np->n_sym, &np->n_val, npl, npr );
		else
			np1 = node( np->n_sym, np->n_val.v_value.v_ival, npl, npr );
*/
		if(	np->n_sym == SYM_IDENT ||
			np->n_sym == SYM_ATTRIBUTE ||
			np->n_sym == SYM_INT_LIT ||
			np->n_sym == SYM_FLOAT_LIT ||
			np->n_sym == SYM_STRING_LIT ||
			np->n_sym == SYM_TYPE )
			np1 = node( np->n_sym, &np->n_val, npl, npr );
		else
			np1 = node( np->n_sym, &np->n_val, npl, npr );
		np1->n_lineno = np->n_lineno;
		np1->n_type = np->n_type;
		np1->n_class = np->n_class;
		np1->n_kind = np->n_kind;
		return( np1 );
	}else
		return( NULL );
}

#if 0
static	void	strext2int( ip, ep )
char	*ip;
char	*ep;
{
	int	ec;

	while( *ep ){
		if( ( *ip = *ep++ ) == '\\' ){
			ec = *ep++;
			if( ec == 'n' )
				*ip = '\n';
			else if( ec == 't' )
				*ip = '\t';
			else
				*ip = ec;
		}
		ip++;
	}
	*ip = '\0';
}
#endif
