%option noyywrap

%{
#include "nab.h"
#include "nabgrm.tab.h"
#include "cgen.h"
#include <string.h>
#include <stdlib.h>

VALUE_T	val;
extern	FILE	*cfp;

%}
%%
\/\/.*		;

[\ \t\f\r]+		;
\n		CG_gennl();


([0-9]+[eE][\+\-]?[0-9]+)|((([0-9]+\.[0-9]*)|(\.[0-9]+))([eE][\+\-]?[0-9]+)?) {
		  val.v_type = T_FLOAT;
		  val.v_value.v_fval = atof( yytext );
		  return( SYM_FLOAT_LIT );
		}
[0-9]+		{ val.v_type = T_INT;
		  val.v_value.v_ival = atoi( yytext );
		  return( SYM_INT_LIT );
		}
\"[^"\n]*	{ if( yytext[ yyleng - 1 ] == '\\' )
			yymore();
		  else{
			yytext[ yyleng ] = input();
			yytext[ yyleng ] = '\0';
			/* orig:  strcpy( yytext, &yytext[ 1 ] ); */
			memmove( yytext, &yytext[1], yyleng );
			yyleng = strlen( yytext );
			val.v_type = T_STRING;
			val.v_value.v_cval = yytext;
			return( SYM_STRING_LIT );
		  }
		}

allocate	return( SYM_ALLOCATE );
assert		return( SYM_ASSERT );
atom		return( SYM_ATOM );
bounds		return( SYM_BOUNDS );
break		return( SYM_BREAK );
continue	return( SYM_CONTINUE );
deallocate	return( SYM_DEALLOCATE );
debug		return( SYM_DEBUG );
delete		return( SYM_DELETE );
dynamic		return( SYM_DYNAMIC );
else		return( SYM_ELSE );
file		return( SYM_FILE );
float		return( SYM_FLOAT );
for		return( SYM_FOR );
hashed		return( SYM_HASHED );
if		return( SYM_IF );
in		return( SYM_IN );
int		return( SYM_INT );
matrix		return( SYM_MATRIX );
molecule	return( SYM_MOLECULE );
point		return( SYM_POINT );
real		return( SYM_FLOAT );
residue		return( SYM_RESIDUE );
return		return( SYM_RETURN );
size_t		return( SYM_SIZE_T );		
string		return( SYM_STRING );
struct		return( SYM_STRUCT );
while		return( SYM_WHILE );

[a-zA-Z][_a-zA-Z0-9]*	{ val.v_type = T_STRING;
			  val.v_value.v_cval = yytext;
			  return( SYM_IDENT );
			}

=		return( SYM_ASSIGN );
\+=		return( SYM_PLUS_ASSIGN );
\-=		return( SYM_MINUS_ASSIGN );
\*=		return( SYM_STAR_ASSIGN );
\/=		return( SYM_SLASH_ASSIGN );
\%=		return( SYM_MODULUS_ASSIGN );
\^=		return( SYM_UPARROW_ASSIGN );
&&		return( SYM_AND );
\|\|		return( SYM_OR );
!		return( SYM_NOT );
\<		return( SYM_LESS );
\<=		return( SYM_LESS_EQUAL );
==		return( SYM_EQUAL );
!=		return( SYM_NOT_EQUAL );
>=		return( SYM_GREATER_EQUAL );
>		return( SYM_GREATER );
=~		return( SYM_MATCH );
!~		return( SYM_DONT_MATCH );
\+		return( SYM_PLUS );
\+\+		return( SYM_PLUS_PLUS );
\-		return( SYM_MINUS );
\-\-		return( SYM_MINUS_MINUS );
\*		return( SYM_STAR );
\/		return( SYM_SLASH );
\%		return( SYM_MODULUS );
@		return( SYM_ATSIGN );
\^		return( SYM_UPARROW );
\(		return( SYM_LPAREN );
\)		return( SYM_RPAREN );
\[		return( SYM_LBRACK );
\]		return( SYM_RBRACK );
\{		return( SYM_LBRACE );
\}		return( SYM_RBRACE );
\.		return( SYM_PERIOD );
\,		return( SYM_COMMA );
;		return( SYM_SEMICOLON );

.		return( SYM_ERROR );

%%
