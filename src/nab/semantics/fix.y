%token	LEFT
%token	NULL
%token	RIGHT	
%token	TEMP
%token	IDENT
%token	STRING
%token	EQUAL
%token	LPAREN
%token	RPAREN
%token	COLON
%token	QUEST
%token	SEMI

%%

sentence	: rw_rule_list ;
rw_rule_list	: rw_rule
		| rw_rule rw_rule_list ;
rw_rule		: IDENT tree COLON case_list SEMI ;
tree		: LPAREN IDENT sub_tree sub_tree RPAREN
		| LPAREN RPAREN
		| NULL ;
sub_tree	: STRING
		| tree ;
case_list	: case
		| case case_list ;
case		: expr QUEST tree
		| tree ;
expr		: LEFT EQUAL IDENT
		| RIGHT EQUAL IDENT ;
