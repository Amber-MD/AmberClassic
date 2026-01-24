
%{
#include <stdio.h>
#include "nab.h"
#include "cgen.h"
#include "errormsg.h"

extern	VALUE_T	val;
static	VALUE_T	v_type;

int yyerror();

# define YYSTYPE_IS_DECLARED 1
# define YYDEBUG 1

%}

%token	<ival>	SYM_ADDRESS
%token	<ival>	SYM_ALLOCATE
%token	<ival>	SYM_AND
%token	<ival>	SYM_ASSERT
%token	<ival>	SYM_ASSIGN
%token	<ival>	SYM_ATOM
%token	<ival>	SYM_ATSIGN
%token	<ival>	SYM_ATTRIBUTE
%token	<ival>	SYM_BOUNDS
%token	<ival>	SYM_BREAK
%token	<ival>	SYM_CALL
%token	<ival>	SYM_COMMA
%token	<ival>	SYM_CONTINUE
%token	<ival>	SYM_DEALLOCATE
%token	<ival>	SYM_DEBUG
%token	<ival>	SYM_DECL
%token	<ival>	SYM_DELETE
%token	<ival>	SYM_DONT_MATCH
%token	<ival>	SYM_DYNAMIC
%token	<ival>	SYM_ELSE
%token	<ival>	SYM_EQUAL
%token	<ival>	SYM_ERROR
%token	<ival>	SYM_FILE
%token	<ival>	SYM_FLOAT
%token	<ival>	SYM_FLOAT_LIT
%token	<ival>	SYM_FOR
%token	<ival>	SYM_FOREACH
%token	<ival>	SYM_GREATER
%token	<ival>	SYM_GREATER_EQUAL
%token	<ival>	SYM_HASHED
%token	<ival>	SYM_IDENT
%token	<ival>	SYM_IF
%token	<ival>	SYM_IN
%token	<ival>	SYM_INDEX
%token	<ival>	SYM_INDIRECT
%token	<ival>	SYM_INT
%token	<ival>	SYM_INT_LIT
%token	<ival>	SYM_LBRACE
%token	<ival>	SYM_LBRACK
%token	<ival>	SYM_LESS
%token	<ival>	SYM_LESS_EQUAL
%token	<ival>	SYM_LIST
%token	<ival>	SYM_LPAREN
%token	<ival>	SYM_MATCH
%token	<ival>	SYM_MATRIX
%token	<ival>	SYM_MINUS
%token	<ival>	SYM_MINUS_ASSIGN
%token	<ival>	SYM_MINUS_MINUS
%token	<ival>	SYM_MODULUS
%token	<ival>	SYM_MODULUS_ASSIGN
%token	<ival>	SYM_MOLECULE
%token	<ival>	SYM_NEGATE
%token	<ival>	SYM_NOT
%token	<ival>	SYM_NOT_EQUAL
%token	<ival>	SYM_OR
%token	<ival>	SYM_PARM
%token	<ival>	SYM_PERIOD
%token	<ival>	SYM_PLUS
%token	<ival>	SYM_PLUS_ASSIGN
%token	<ival>	SYM_PLUS_PLUS
%token	<ival>	SYM_POINT
%token	<ival>	SYM_POINTS_TO
%token	<ival>	SYM_RBRACE
%token	<ival>	SYM_RBRACK
%token	<ival>	SYM_RESIDUE
%token	<ival>	SYM_RETURN
%token	<ival>	SYM_RPAREN
%token	<ival>	SYM_SEMICOLON
%token	<ival>	SYM_SIZE_T
%token	<ival>	SYM_SLASH
%token	<ival>	SYM_SLASH_ASSIGN
%token	<ival>	SYM_STAR
%token	<ival>	SYM_STAR_ASSIGN
%token	<ival>	SYM_STMTLIST
%token	<ival>	SYM_STRING
%token	<ival>	SYM_STRING_LIT
%token	<ival>	SYM_STRUCT
%token	<ival>	SYM_TEST
%token	<ival>	SYM_TYPE
%token	<ival>	SYM_UPARROW
%token	<ival>	SYM_UPARROW_ASSIGN
%token	<ival>	SYM_WHILE

%type	<npval>	actuals
%type	<npval>	addop
%type	<npval>	a_expr
%type	<npval>	alloc_stmt
%type	<npval>	a_parm
%type	<npval>	ap_list
%type	<npval>	ar_lval
%type	<npval>	asize
%type	<npval>	as_list
%type	<npval>	aspec
%type	<npval>	assert_stmt
%type	<npval>	assignop
%type	<npval>	at_lval
%type	<npval>	break_stmt
%type	<npval>	cmpd_stmt
%type	<npval>	conj
%type	<npval>	continue_stmt
%type	<npval>	dbg_list
%type	<npval>	dealloc_stmt
%type	<npval>	debug_stmt
%type	<npval>	def
%type	<npval>	defpart
%type	<npval>	defs
%type	<npval>	delete_stmt
%type	<npval>	disj
%type	<npval>	e_list
%type	<npval>	expr
%type	<npval>	expr_stmt
%type	<npval>	factor
%type	<npval>	field
%type	<npval>	field_list
%type	<npval>	f_defpart
%type	<npval>	for_count
%type	<npval>	for_ctrl
%type	<npval>	for_expr
%type	<npval>	for_hdr
%type	<npval>	for_in
%type	<npval>	formals
%type	<npval>	for_stmt
%type	<npval>	for_test_expr
%type	<npval>	f_parm
%type	<npval>	fp_list
%type	<npval>	f_stmtpart
%type	<npval>	func_body
%type	<npval>	func_decl
%type	<npval>	func_def
%type	<npval>	func_hdr
%type	<npval>	id
%type	<npval>	id_list
%type	<npval>	if_hdr
%type	<npval>	if_stmt
%type	<npval>	i_list
%type	<npval>	incr
%type	<npval>	incrop
%type	<npval>	lval
%type	<npval>	lv_decls
%type	<npval>	mulop
%type	<npval>	num
%type	<npval>	primary
%type	<npval>	program
%type	<npval>	relop
%type	<npval>	return_stmt
%type	<npval>	rval
%type	<npval> simple_type
%type	<npval>	stmt
%type	<npval>	stmtpart
%type	<npval>	stmts
%type	<npval>	string
%type	<npval> struct_type
%type	<npval>	term
%type	<npval>	type
%type	<npval>	type_decl
%type	<npval>	unop
%type	<npval>	var
%type	<npval>	var_decl
%type	<npval>	var_list
%type	<npval>	while_hdr
%type	<npval>	while_stmt

%%
program		: defpart stmtpart { CG_genend(); } ;
defpart		: defs
		| 		{ $$ = NULL; } ;
defs		: def	
		| def defs ;
stmtpart	: stmts
		| 		{ $$ = NULL; };
stmts		: stmt
		| stmts stmt ;

def		: type_decl
		| var_decl
		| func_decl
		| func_def ;

type_decl	: struct_type SYM_SEMICOLON
				{ $$ = node( SYM_DECL, 0, $1, 0 );
				  CG_genvardecl( $$, 0, 0, 0 ); } ;
var_decl	: type var_list SYM_SEMICOLON 
				{ $$ = node( SYM_DECL, 0, $1, $2 );
				  CG_genvardecl( $$, 0, 0, 0 ); } ;
type	: simple_type	{ $$ = $1; }
		| struct_type	{ $$ = $1; } ;
simple_type		: SYM_INT	{ 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_INT;
			$$ = node( SYM_TYPE, &v_type, 0, 0 ); }
		| SYM_SIZE_T	{
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_SIZE_T;
			$$ = node( SYM_TYPE, &v_type, 0, 0 ); }
		| SYM_FLOAT	{ 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_FLOAT;
			$$ = node( SYM_TYPE, &v_type, 0, 0 ); }
		| SYM_STRING	{ 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_STRING;
			$$ = node( SYM_TYPE, &v_type, 0, 0 ); }
		| SYM_FILE	{ 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_FILE;
			$$ = node( SYM_TYPE, &v_type, 0, 0 ); }
		| SYM_POINT	{ 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_POINT;
			$$ = node( SYM_TYPE, &v_type, 0, 0 ); }
		| SYM_BOUNDS	{ 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_BOUNDS;
			$$ = node( SYM_TYPE, &v_type, 0, 0 ); }
		| SYM_MATRIX	{ 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_MATRIX;
			$$ = node( SYM_TYPE, &v_type, 0, 0 ); }
		| SYM_ATOM	{ 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_ATOM;
			$$ = node( SYM_TYPE, &v_type, 0, 0 ); }
		| SYM_RESIDUE	{ 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_RESIDUE;
			$$ = node( SYM_TYPE, &v_type, 0, 0 ); }
		| SYM_MOLECULE	{ 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_MOLECULE;
			$$ = node( SYM_TYPE, &v_type, 0, 0 ); } ;
struct_type	: SYM_STRUCT id {
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_USER;
			$$ = node( SYM_TYPE, &v_type, node( SYM_STRUCT, 0, $2, 0 ), 0 ); }
		| SYM_STRUCT id SYM_LBRACE field_list SYM_RBRACE {
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_USER;
			$$ = node( SYM_TYPE, &v_type, node( SYM_STRUCT, 0, $2, $4 ), 0 ); } ;
field_list : field	{ $$ = node( SYM_LIST, 0, $1, 0 ); }
		| field field_list
					{ $$ = node( SYM_LIST, 0, $1, $2 ); } ;
field	: simple_type id_list SYM_SEMICOLON
					{ $$ = node( SYM_DECL, 0, $1, $2 ); } ;
id_list	: id		{ $$ = node( SYM_LIST, 0, $1, 0 ); }
		| id SYM_COMMA id_list
					{ $$ = node( SYM_LIST, 0, $1, $3 ); } ;
var_list	: var		{ $$ = node( SYM_LIST, 0, $1, 0 ); }
		| var SYM_COMMA var_list { $$ = node( SYM_LIST, 0, $1, $3 ); } ;
var		: id 		{ $$ = $1; }
		| id SYM_LBRACK aspec SYM_RBRACK
				{ $$ = node( SYM_LBRACK, 0, $1, $3 ); } ;
aspec		: as_list	{ $$ = $1; } 
		| SYM_HASHED	{ $$ = node( SYM_HASHED, 0, 0, 0 ); } ;
as_list		: asize		{ $$ = node( SYM_LIST, 0, $1, 0 ); }
		| asize SYM_COMMA as_list
				{ $$ = node( SYM_LIST, 0, $1, $3 ); } ;
asize		: expr		{ $$ = $1; }
		| SYM_DYNAMIC	{ $$ = node( SYM_DYNAMIC, 0, 0, 0 ); } ;

func_decl	: func_hdr SYM_SEMICOLON
				{ CG_genop( NULL, SYM_SEMICOLON);
				  CG_genfend( NULL ); } 
		| func_hdr id SYM_SEMICOLON 
				{ CG_genop( NULL, SYM_SEMICOLON );
				  CG_genfend( $2 ); } ;

func_def	: func_hdr 	{ CG_genfstart(); }
		  func_body 	{ CG_genfend( NULL ); } ;
func_hdr	: type id 	{ CG_genfhdr( $1, $2 ); }
		  SYM_LPAREN formals SYM_RPAREN
				{ CG_genplist( $5 ); } ;
formals		: fp_list	{ $$ = $1; }
		| 		{ $$ = NULL; } ;
fp_list		: f_parm	{ $$ = node( SYM_LIST, 0, $1, 0 ); }
		| f_parm SYM_COMMA fp_list
				{ $$ = node( SYM_LIST, 0, $1, $3 ); } ;
f_parm		: type var 
				{ $$ = node( SYM_LIST, 0, $2, NULL ); 
			  	  $$ = node( SYM_DECL, 0, $1, $$ ); } ;
func_body	: SYM_LBRACE	{ CG_genpdecls();
				  CG_genop( NULL, SYM_LBRACE ); }
		  f_defpart	{ CG_genedefs( TRUE ); }
		  f_stmtpart SYM_RBRACE
				{ CG_genestmts( TRUE );
				  CG_genop( NULL, SYM_RBRACE ); }
		  SYM_SEMICOLON { $$=NULL; } ;
f_defpart	: lv_decls 
		| 		{ $$ = NULL; } ;
lv_decls	: var_decl
		| lv_decls var_decl ;
f_stmtpart	: stmts
		| 		{ $$ = NULL; };

stmt		: expr_stmt
		| alloc_stmt
		| assert_stmt
		| break_stmt
		| continue_stmt
		| cmpd_stmt
		| dealloc_stmt
		| debug_stmt
		| delete_stmt
		| for_stmt
		| if_stmt
		| return_stmt
		| while_stmt ;

alloc_stmt	: SYM_ALLOCATE expr SYM_SEMICOLON 
				{ CG_genmain();
				  $$ = node( SYM_ALLOCATE, 0, 0, $2 );
				  CG_genexpr( $$ );
				  CG_genop( NULL, SYM_SEMICOLON ); } ;
assert_stmt	: SYM_ASSERT expr SYM_SEMICOLON
				{ CG_genmain();
				  $$ = node( SYM_ASSERT, 0, 0, $2 );
				  CG_genassert( $$ ); } ;
break_stmt	: SYM_BREAK SYM_SEMICOLON
				{ CG_genmain();
				  CG_genrword( SYM_BREAK );
				  CG_genop( NULL, SYM_SEMICOLON ); } ;
cmpd_stmt	: SYM_LBRACE	{ CG_genmain();
				  CG_genop( NULL, SYM_LBRACE ); }
	 	  stmts SYM_RBRACE
				{ CG_genop( NULL, SYM_RBRACE ); } ;
continue_stmt	: SYM_CONTINUE SYM_SEMICOLON
				{ CG_genmain();
				  CG_genrword( SYM_CONTINUE );
				  CG_genop( NULL, SYM_SEMICOLON ); } ;
dealloc_stmt	: SYM_DEALLOCATE expr SYM_SEMICOLON 
				{ CG_genmain();
				  $$ = node( SYM_DEALLOCATE, 0, 0, $2 );
			 	  CG_genexpr($$);
				  CG_genop( NULL, SYM_SEMICOLON ); } ;
debug_stmt	: SYM_DEBUG dbg_list SYM_SEMICOLON 
				{ CG_genmain();
				  $$ = node( SYM_DEBUG, 0, 0, $2 );
				  CG_gendebug( $$ ); } ;
delete_stmt	: SYM_DELETE expr SYM_SEMICOLON
				{ CG_genmain();
				  $$ = node( SYM_DELETE, 0, 0, $2 );
			 	  CG_genexpr($$);
				  CG_genop( NULL, SYM_SEMICOLON ); } ;
expr_stmt	: expr SYM_SEMICOLON
				{ CG_genmain(); CG_genexpr( $1 );
				  CG_genop( NULL, SYM_SEMICOLON ); } ;
if_stmt		: if_hdr stmt 
		| if_hdr stmt SYM_ELSE
				{ CG_genrword( SYM_ELSE ); } stmt ;
for_stmt	: for_hdr stmt ;
return_stmt	: SYM_RETURN	{ CG_genmain();
				  CG_genrword(SYM_RETURN);
				  CG_genop( NULL, SYM_LPAREN ); }
		  expr		{ CG_genexpr( $3 ); }
		  SYM_SEMICOLON	{ CG_genop( NULL, SYM_RPAREN );
				  CG_genop( NULL, SYM_SEMICOLON ); } ;
while_stmt	: while_hdr stmt ;

if_hdr		: SYM_IF	{ CG_genmain(); CG_genrword( SYM_IF ); }
		  SYM_LPAREN	{ CG_genop( NULL, SYM_LPAREN ); }
		  expr		{ $<npval>$ = node( SYM_TEST, 0, 0, $5 );
				  CG_genexpr( $<npval>$ ); }
		  SYM_RPAREN	{ CG_genop( NULL, SYM_RPAREN ); } ;

for_hdr		: SYM_FOR	{ CG_genmain(); CG_genrword( SYM_FOR ); }
		  SYM_LPAREN	{ CG_genop( NULL, SYM_LPAREN ); }
		  for_ctrl SYM_RPAREN
				{ CG_genop( NULL, SYM_RPAREN ); } ;
for_ctrl	: for_in
		| for_count ;
for_in		: id SYM_IN id	{ $$ = node( SYM_FOREACH, 0, $1, $3 );
				  CG_genexpr( $$ ); } ;
for_count	: for_expr SYM_SEMICOLON
				{ CG_genop( NULL, SYM_SEMICOLON ); }
		  for_test_expr SYM_SEMICOLON
				{ CG_genop( NULL, SYM_SEMICOLON ); } for_expr ;
for_expr	: expr		{ CG_genexpr( $1 ); }
		| 		{ $$ = NULL; } ;
for_test_expr	: expr		{ $$ = node( SYM_TEST, 0, 0, $1 );
				  CG_genexpr( $$ ); }
		|		{ $$ = NULL; } ;

while_hdr	: SYM_WHILE	{ CG_genmain(); CG_genrword( SYM_WHILE ); }
		  SYM_LPAREN	{ CG_genop( NULL, SYM_LPAREN ); }
		  expr 		{ $<npval>$ = node( SYM_TEST, 0, 0, $5 );
				  CG_genexpr( $<npval>$ ); }
		  SYM_RPAREN	{ CG_genop( NULL, SYM_RPAREN ); } ;

dbg_list	: SYM_LPAREN e_list SYM_RPAREN 
				{ $$ = $2; }
		| e_list	{ $$ = $1; } ;
e_list		: expr		{ $$ = node( SYM_LIST, 0, $1, 0 ); }
		| expr SYM_COMMA e_list
				{ $$ = node( SYM_LIST, 0, $1, $3 ); } ;

expr		: rval		{ $$ = $1; }
		| lval assignop expr
				{ $$ = node( $<ival>2, 0, $1, $3 ); } ;
lval		: id		{ $$ = $1; }
		| ar_lval	{ $$ = $1; }
		| at_lval	{ $$ = $1; } ;
ar_lval		: lval SYM_LBRACK i_list SYM_RBRACK
				{ $$ = node( SYM_LBRACK, 0, $1, $3 ); } ; 
at_lval		: lval SYM_PERIOD SYM_IDENT
				{ $$ = node( SYM_PERIOD, 0, $1,
				  node( SYM_ATTRIBUTE, &val, 0, 0 ) ); } ; 
rval		: disj		{ $$ = $1; }
		| disj SYM_OR rval
				{ $$ = node( SYM_OR, 0, $1, $3 ); } ;
disj		: conj		{ $$ = $1; }
		| conj SYM_AND disj
				{ $$ = node( SYM_AND, 0, $1, $3 ); } ;
conj		: a_expr	{ $$ = $1; } 
		| a_expr relop a_expr
				{ $$ = node( $<ival>2, 0, $1, $3 ); } ;
a_expr		: term		{ $$ = $1; }
		| term addop a_expr
				{ $$ = node( $<ival>2, 0, $1, $3 ); } ;
term		: factor	{ $$ = $1; }
		| factor mulop term
				{ $$ = node( $<ival>2, 0, $1, $3 ); } ;
factor		: primary 	{ $$ = $1; }
		| primary SYM_UPARROW factor
				{ $$ = node( SYM_UPARROW, 0, $1, $3 ); };
primary		: lval		{ $$ = $1; }
		| num		{ $$ = $1; }
		| string	{ $$ = $1; }
		| incr		{ $$ = $1; }
		| unop primary	{ $$ = node( $<ival>1, 0, NULL, $2 ); }
		| id SYM_LPAREN actuals SYM_RPAREN
				{ $$ = node( SYM_CALL, 0, $1, $3 ); }
		| SYM_LPAREN expr SYM_RPAREN
				{ $$ = node( SYM_LPAREN, 0, NULL, $2 ); } ;
incr		: incrop lval	{ $$ = node( $<ival>1, 0, 0, $2 ); }
		| lval incrop	{ $$ = node( $<ival>2, 0, $1, 0 ); } ;
actuals		: ap_list	{ $$ = $1; }
		| 		{ $$ = NULL; } ;
ap_list		: a_parm	{ $$ = node( SYM_LIST, 0, $1, 0 ); }
		| a_parm SYM_COMMA ap_list
				{ $$ = node( SYM_LIST, 0, $1, $3 ); } ;
a_parm		: expr		{ $$ = node( SYM_PARM, 0, 0, $1 ); } ;
i_list		: expr		{ $$ = node( SYM_INDEX, 0, $1, 0 ); }
		| expr SYM_COMMA i_list
				{ $$ = node( SYM_INDEX, 0, $1, $3 ); } ;
assignop	: SYM_ASSIGN	{ $<ival>$ = SYM_ASSIGN; }
		| SYM_PLUS_ASSIGN
				{ $<ival>$ = SYM_PLUS_ASSIGN; }
		| SYM_MINUS_ASSIGN
				{ $<ival>$ = SYM_MINUS_ASSIGN; }
		| SYM_STAR_ASSIGN
				{ $<ival>$ = SYM_STAR_ASSIGN; }
		| SYM_SLASH_ASSIGN
				{ $<ival>$ = SYM_SLASH_ASSIGN; }
		| SYM_MODULUS_ASSIGN
				{ $<ival>$ = SYM_MODULUS_ASSIGN; }
		| SYM_UPARROW_ASSIGN
				{ $<ival>$ = SYM_UPARROW_ASSIGN; } ;
relop		: SYM_LESS	{ $<ival>$ = SYM_LESS; }
		| SYM_LESS_EQUAL
				{ $<ival>$ = SYM_LESS_EQUAL; }
		| SYM_EQUAL	{ $<ival>$ = SYM_EQUAL; }
		| SYM_NOT_EQUAL	{ $<ival>$ = SYM_NOT_EQUAL; }
		| SYM_GREATER_EQUAL
				{ $<ival>$ = SYM_GREATER_EQUAL; }
		| SYM_GREATER	{ $<ival>$ = SYM_GREATER; }
		| SYM_MATCH	{ $<ival>$ = SYM_MATCH; }
		| SYM_DONT_MATCH
				{ $<ival>$ = SYM_DONT_MATCH; }
		| SYM_IN	{ $<ival>$ = SYM_IN; } ;
addop		: SYM_PLUS	{ $<ival>$ = SYM_PLUS; }
		| SYM_MINUS	{ $<ival>$ = SYM_MINUS; } ;
mulop		: SYM_STAR	{ $<ival>$ = SYM_STAR; }
		| SYM_SLASH	{ $<ival>$ = SYM_SLASH; }
		| SYM_MODULUS 	{ $<ival>$ = SYM_MODULUS; }
		| SYM_ATSIGN 	{ $<ival>$ = SYM_ATSIGN; } ;
incrop		: SYM_PLUS_PLUS	{ $<ival>$ = SYM_PLUS_PLUS; } 
		| SYM_MINUS_MINUS
				{ $<ival>$ = SYM_MINUS_MINUS; } ; 
unop		: SYM_MINUS	{ $<ival>$ = SYM_NEGATE; }
		| SYM_NOT 	{ $<ival>$ = SYM_NOT; } ;

id		: SYM_IDENT	{ $$ = node( SYM_IDENT, &val, 0, 0 ); } ;
num		: SYM_INT_LIT 	{ $$ = node( SYM_INT_LIT, &val, 0, 0 ); }
		| SYM_FLOAT_LIT	{ $$ = node( SYM_FLOAT_LIT, &val, 0, 0 ); } ;
string		: SYM_STRING_LIT
				{ $$ = node( SYM_STRING_LIT, &val, 0, 0 ); } ;
%%

#include "lex.yy.c"
