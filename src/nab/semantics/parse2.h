
// header file for parse2.c.
// for documentation, look at parse2.c

#include <stdio.h>

#ifndef PARSE2_H
#define PARSE2_H

// structs

typedef	struct	atrec_t	{
	char	*a_name;
	int	a_type;
	int	a_class;
	int	a_kind;
	int	a_in;
	int	a_access;
} ATREC_T;

#define	DEFNAME_SIZE	32
typedef	struct	def_t	{
	char	d_name[ DEFNAME_SIZE ];
	int	d_val;
	int	d_used;
} DEF_T;

typedef	struct	resword {
	char	*r_name;
	int	r_token;
} RESWORD;

typedef	struct	trow_t	{
	char	*t_name;
	int	t_nval;
	int	t_nrvals;
	int	*t_rvals;
} TROW_T;

typedef	struct	ttab_t	{
	int	t_index;
	int	t_nrows;
	TROW_T	**t_rows;
} TTAB_T;

typedef	struct	rule_t	{
	int	r_sym;
	int	r_prec;
	char	*r_print;
	int	r_lclass;
	int	r_rclass;
	int	r_octok;
	int	r_ocval;
	int	r_nakind;
	int	*r_akind;
	int	r_nlkind;
	int	*r_lkind;
	int	r_nrkind;
	int	*r_rkind;
	int	r_nokind;
	int	*r_oktok;
	int	*r_okval;
	int	*r_okqtok;
	int	*r_okqval;
	int	r_ntypes;
	TTAB_T	*r_types;
} RULE_T;


// function definitions
// these were automatically scraped from parse2.c using this tool: http://blog.olkie.com/2013/11/05/online-c-function-prototype-header-generator-tool/
int main(int argc, char *argv[] );
static char *str2hstr(char str[], char hstr[] );
static DEF_T *finddef_withval(int val, int n_tab, DEF_T tab[] );
static DEF_T *finddef(char name[], int n_tab, DEF_T tab[] );
static int check_rule(RULE_T *rp );
static int check_type(TTAB_T *ttp );
static int class(FILE *fp, RULE_T *rp );
static int countudefs(int n_tab, DEF_T tab[] );
static int defcmp(void const *v1, void const *v2 );
static int getovalue(FILE *fp, char ctxt[], int n_tab, DEF_T tab[], int *otok, int *oval, int *qtok, int *qval );
static int getsetvalue(FILE *fp, char ctxt[], int n_tab, DEF_T tab[], int *val, int *t_idx );
static int getsyms(FILE *fp, DEF_T symtab[] );
static int gettcka(FILE* fp, int *l_utype, int *v_error, int *n_tvtab, DEF_T tvtab[], int *n_cvtab, DEF_T cvtab[], int *n_kvtab, DEF_T kvtab[], int *n_avtab, DEF_T avtab[] );
static int gettoken(FILE *fp );
static int hprint_rule(char* rulePath, RULE_T * rp, char* outputFolder);
static int kind(FILE *fp, RULE_T *rp );
static int mk_checkexpr(char* output_file_path, int n_rules );
static int mk_controlfile(char* html_output_folder);
static int mk_frameset(char* html_output_folder);
static int operator(FILE *fp, RULE_T *rp );
static int parse_attrs(FILE *fp, int *n_attab, ATREC_T attab[] );
static int parse_rule(FILE *fp, RULE_T *rp );
static int split(char str[], char *fields[], char *fsep );
static int type(FILE *fp, RULE_T *rp );
static int getvecvalue( FILE *fp, TROW_T *trp );
static RULE_T *save_rule(RULE_T *rp );
static TROW_T *findtrow(char name[], int n_trows, TROW_T *trows[] );
static void clear_rule(RULE_T *rp );
static void fprint_ovalue(FILE *fp, char *oval, int n_tab, DEF_T tab[] );
static void fprint_rule(FILE *fp, RULE_T *rp );
static void fprint_set(FILE *fp, int set, int n_tab, DEF_T tab[] );
static void mk_checkattr(FILE *fp );
static void skiptonl(FILE *fp );
static void ungettoken();
void mk_rulecase(FILE *fp, RULE_T *rp );
void mk_typetabent(FILE *fp, int first, RULE_T *rp );

/* tokens from the rule files & attribute file	*/

#define	TOK_EOF		0

#define	TOK_AKIND	1
#define	TOK_CLASS	2
#define	TOK_CODE	3
#define	TOK_IF		4
#define	TOK_KIND	5
#define	TOK_LEFT	6
#define	TOK_OPERATOR	7
#define	TOK_OUTPUT	8
#define	TOK_PREC	9
#define	TOK_PRINT	10
#define	TOK_RIGHT	11
#define	TOK_SYM		12
#define	TOK_TYPE	13
#define	TOK_USE		14

#define	TOK_IDENT	15
#define	TOK_STRING	16
#define	TOK_INT		17

#define	TOK_COMMENT	18

#define	TOK_EQUAL	19
#define	TOK_LBRACE	20
#define	TOK_RBRACE	21
#define	TOK_LBRACK	22
#define	TOK_RBRACK	23
#define	TOK_COMMA	24
#define	TOK_NL		25
#define	TOK_ERROR	26

#endif
