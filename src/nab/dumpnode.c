#include <stdio.h>
#include "nab.h"
#include "nabgrm.tab.h"

#define	INCR	2
#define	DUMPFILE	"dump.file"

static	FILE	*dfp = NULL;

void	dumpexpr( FILE *, NODE_T *, int );
void	dumpnode( FILE *, NODE_T *, int );

void	dumpexpr( FILE *fp, NODE_T *np, int indent )
{

	if( !fp ){
		if( !dfp ){
			if( ( dfp = fopen( DUMPFILE, "w" ) ) == NULL )
			return;
		}
		fp = dfp;
	}
	if( np ){
		dumpnode( fp, np, indent );
		dumpexpr( fp, np->n_left, indent + INCR );
		dumpexpr( fp, np->n_right, indent + INCR );
	}
}

void	dumpnode( FILE *fp, NODE_T *np, int indent )
{

	fprintf( fp, "%*s", indent, "" );
	fprintf( fp, "%10p: lf = %10p rt = %10p tp,cl,kn = %2d,%2d,%2d ",
		np, np->n_left, np->n_right, np->n_type, np->n_class, np->n_kind );
	switch( np->n_sym ){
	case SYM_ADDRESS :
		fprintf( fp, "SYM_ADDRESS\n" );
		break;
	case SYM_ALLOCATE :
		fprintf( fp, "SYM_ALLOCATE\n" );
		break;
	case SYM_AND :
		fprintf( fp, "SYM_AND\n" );
		break;
	case SYM_ASSERT :
		fprintf( fp, "SYM_ASSERT\n" );
		break;
	case SYM_ASSIGN :
		fprintf( fp, "SYM_ASSIGN\n" );
		break;
	case SYM_ATOM :
		fprintf( fp, "SYM_ATOM\n" );
		break;
	case SYM_ATSIGN :
		fprintf( fp, "SYM_ATSIGN\n" );
		break;
	case SYM_ATTRIBUTE :
		fprintf( fp, "SYM_ATTRIBUTE = %s\n", np->n_val.v_value.v_cval );
		break;
	case SYM_BOUNDS :
		fprintf( fp, "SYM_BOUNDS\n" );
		break;
	case SYM_BREAK :
		fprintf( fp, "SYM_BREAK\n" );
		break;
	case SYM_CALL :
		fprintf( fp, "SYM_CALL\n" );
		break;
	case SYM_COMMA :
		fprintf( fp, "SYM_COMMA\n" );
		break;
	case SYM_CONTINUE :
		fprintf( fp, "SYM_CONTINUE\n" );
		break;
	case SYM_DEALLOCATE :
		fprintf( fp, "SYM_DEALLOCATE\n" );
		break;
	case SYM_DEBUG :
		fprintf( fp, "SYM_DEBUG\n" );
		break;
	case SYM_DECL :
		fprintf( fp, "SYM_DECL\n" );
		break;
	case SYM_DELETE :
		fprintf( fp, "SYM_DELETE\n" );
		break;
	case SYM_DONT_MATCH :
		fprintf( fp, "SYM_DONT_MATCH\n" );
		break;
	case SYM_DYNAMIC :
		fprintf( fp, "SYM_DYNAMIC\n" );
		break;
	case SYM_ELSE :
		fprintf( fp, "SYM_ELSE\n" );
		break;
	case SYM_EQUAL :
		fprintf( fp, "SYM_EQUAL\n" );
		break;
	case SYM_ERROR :
		fprintf( fp, "SYM_ERROR\n" );
		break;
	case SYM_FILE :
		fprintf( fp, "SYM_FILE\n" );
		break;
	case SYM_FLOAT :
		fprintf( fp, "SYM_FLOAT\n" );
		break;
	case SYM_FLOAT_LIT :
		fprintf( fp, "SYM_FLOAT_LIT = %8.3f\n", np->n_val.v_value.v_fval );
		break;
	case SYM_FOR :
		fprintf( fp, "SYM_FOR\n" );
		break;
	case SYM_FOREACH :
		fprintf( fp, "SYM_FOREACH\n" );
		break;
	case SYM_GREATER :
		fprintf( fp, "SYM_GREATER\n" );
		break;
	case SYM_GREATER_EQUAL :
		fprintf( fp, "SYM_GREATER_EQUAL\n" );
		break;
	case SYM_HASHED :
		fprintf( fp, "SYM_HASHED\n" );
		break;
	case SYM_IDENT :
		fprintf( fp, "SYM_IDENT = %s\n", np->n_val.v_value.v_cval );
		break;
	case SYM_IF :
		fprintf( fp, "SYM_IF\n" );
		break;
	case SYM_IN :
		fprintf( fp, "SYM_IN\n" );
		break;
	case SYM_INDEX :
		fprintf( fp, "SYM_INDEX = %d\n", np->n_val.v_value.v_ival );
		break;
	case SYM_INDIRECT :
		fprintf( fp, "SYM_INDIRECT\n" );
		break;
	case SYM_INT :
		fprintf( fp, "SYM_INT\n" );
		break;
	case SYM_INT_LIT :
		fprintf( fp, "SYM_INT_LIT = %d\n", np->n_val.v_value.v_ival );
		break;
	case SYM_LBRACE :
		fprintf( fp, "SYM_LBRACE\n" );
		break;
	case SYM_LBRACK :
		fprintf( fp, "SYM_LBRACK\n" );
		break;
	case SYM_LESS :
		fprintf( fp, "SYM_LESS\n" );
		break;
	case SYM_LESS_EQUAL :
		fprintf( fp, "SYM_LESS_EQUAL\n" );
		break;
	case SYM_LIST :
		fprintf( fp, "SYM_LIST\n" );
		break;
	case SYM_LPAREN :
		fprintf( fp, "SYM_LPAREN\n" );
		break;
	case SYM_MATCH :
		fprintf( fp, "SYM_MATCH\n" );
		break;
	case SYM_MATRIX :
		fprintf( fp, "SYM_MATRIX\n" );
		break;
	case SYM_MINUS :
		fprintf( fp, "SYM_MINUS\n" );
		break;
	case SYM_MINUS_ASSIGN :
		fprintf( fp, "SYM_MINUS_ASSIGN\n" );
		break;
	case SYM_MINUS_MINUS :
		fprintf( fp, "SYM_MINUS_MINUS\n" );
		break;
	case SYM_MODULUS :
		fprintf( fp, "SYM_MODULUS\n" );
		break;
	case SYM_MODULUS_ASSIGN :
		fprintf( fp, "SYM_MODULUS_ASSIGN\n" );
		break;
	case SYM_MOLECULE :
		fprintf( fp, "SYM_MOLECULE\n" );
		break;
	case SYM_NEGATE :
		fprintf( fp, "SYM_NEGATE\n" );
		break;
	case SYM_NOT :
		fprintf( fp, "SYM_NOT\n" );
		break;
	case SYM_NOT_EQUAL :
		fprintf( fp, "SYM_NOT_EQUAL\n" );
		break;
	case SYM_OR :
		fprintf( fp, "SYM_OR\n" );
		break;
	case SYM_PARM :
		fprintf( fp, "SYM_PARM\n" );
		break;
	case SYM_PERIOD :
		fprintf( fp, "SYM_PERIOD\n" );
		break;
	case SYM_PLUS :
		fprintf( fp, "SYM_PLUS\n" );
		break;
	case SYM_PLUS_ASSIGN :
		fprintf( fp, "SYM_PLUS_ASSIGN\n" );
		break;
	case SYM_PLUS_PLUS :
		fprintf( fp, "SYM_PLUS_PLUS\n" );
		break;
	case SYM_POINT :
		fprintf( fp, "SYM_POINT\n" );
		break;
	case SYM_POINTS_TO :
		fprintf( fp, "SYM_POINTS_TO\n" );
		break;
	case SYM_RBRACE :
		fprintf( fp, "SYM_RBRACE\n" );
		break;
	case SYM_RBRACK :
		fprintf( fp, "SYM_RBRACK\n" );
		break;
	case SYM_RESIDUE :
		fprintf( fp, "SYM_RESIDUE\n" );
		break;
	case SYM_RETURN :
		fprintf( fp, "SYM_RETURN\n" );
		break;
	case SYM_RPAREN :
		fprintf( fp, "SYM_RPAREN\n" );
		break;
	case SYM_SEMICOLON :
		fprintf( fp, "SYM_SEMICOLON\n" );
		break;
	case SYM_SIZE_T :
		fprintf( fp, "SYM_SIZE_T\n" );
		break;
	case SYM_SLASH :
		fprintf( fp, "SYM_SLASH\n" );
		break;
	case SYM_SLASH_ASSIGN :
		fprintf( fp, "SYM_SLASH_ASSIGN\n" );
		break;
	case SYM_STAR :
		fprintf( fp, "SYM_STAR\n" );
		break;
	case SYM_STAR_ASSIGN :
		fprintf( fp, "SYM_STAR_ASSIGN\n" );
		break;
	case SYM_STMTLIST :
		fprintf( fp, "SYM_STMTLIST\n" );
		break;
	case SYM_STRING :
		fprintf( fp, "SYM_STRING\n" );
		break;
	case SYM_STRING_LIT :
		fprintf( fp, "SYM_STRING_LIT = \"%s\"\n", np->n_val.v_value.v_cval );
		break;
	case SYM_STRUCT :
		fprintf( fp, "SYM_STRUCT\n" );
		break;
	case SYM_TEST :
		fprintf( fp, "SYM_TEST\n" );
		break;
	case SYM_TYPE :
		fprintf( fp, "SYM_TYPE = %d\n", np->n_val.v_value.v_ival );
		break;
	case SYM_UPARROW :
		fprintf( fp, "SYM_UPARROW\n" );
		break;
	case SYM_UPARROW_ASSIGN :
		fprintf( fp, "SYM_UPARROW_ASSIGN\n" );
		break;
	case SYM_WHILE :
		fprintf( fp, "SYM_WHILE\n" );
		break;
	default :
		fprintf( fp, "dumpnode: unknown symbol %d\n", np->n_sym );
		break;
	}
}
