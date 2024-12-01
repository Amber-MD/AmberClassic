/*
 * parse2.c is a modified version of parse.c which takes its arguments on the command line instead
 * of through preprocessor definitions.  This way, it can be used with the CMake build system and
 * Amber Host Tools.
 *
 * In addition to the command-line interface, I also made a few other changes:
 *  * converted old-style K&R function definitions to ISO C
 *  * created a header file and moved what I could to it
 *  * added functions' arguments to their prototypes, so that the compiler can check them for errors
 *  * made it open the rules files for rewriting with mode "wb", so it will still use LF newlines even on Windows (that's the bug I started this project to track down)
 *  * renamed some variables in and around void main() to make it clear what they do
 *  * added Windows compatibility
 *
 * -JS
 *
 * Original documentation follows:
 *
 *	This program converts the actions described in the tables
 *	SYM_*.rules into checkexpr(), the part of the nab compiler that 
 *	ensures that the operands of the various operators have the
 *	correct attributes and sets the attributes of the operator's 
 *	result.
 *
 *	Three operators require special handling.  These are the array
 *	operator LBRACK and the identifier operators IDENT and ATTRIBUTE.
 *
 *	When an array element is used in a nab program, say a[i], this
 *	is converted into this parse tree:
 *
 *		( LBRACK ( IDENT A ) ( INDEX ( IDENT i ) ) )
 *
 *	Because nab supports both ordinary (indexed by int valued expr's)
 *	and hashed (indexed by string valued expr's), the INDEX node needs
 *	to know what kind of array it is associated with in order to decided
 *	if its subexpression is valid.  However the parse tree contains no
 *	upward pointers and it is impossible for the INDEX node to determine
 *	its context.  This problem is solved by maintaining a stack of active
 *	array refs -- astk[] and astkp.  Each time an INDEX node is encountered
 *	it checks the kind of the array at the top of this stack and if it is
 *	K_ARRAY or K_DARRAY, requires an INT valued expr; if it is K_HASHED,
 *	it requires a STRING valued expr.
 *
 *	IDENT requires special handling because the action associated
 *	with the IDENT depends on whether it is appearing a declaration or
 *	a statement.  In a statement, the IDENT must have been previously
 *	declared and the action is to introduce its attr's into the expr;
 *	in a decl, it must not have been previously declared and when 1st
 *	encountered will have only its type known;  The class and kind 
 *	can not be determined until the next symbol - '(' for func's, '['
 *	for arrays and ',' or ';' for vars has been seen.  The solution here
 *	is to provide a special case entry for IDENT that check's that it 
 *	has been declared.
 *
 *	ATTRIBUTE also requires special handling.  The actions here are to
 *	look up the attribute and if it is known, return its type, kind and
 *	class and to print an errormsg if the attribute is unknown.  This
 *	is done by passing the ATTRIBUTE node to the procedure checkattr()
 *	which does the lookup and set the n_type, n_class and n_kind fields
 *	of the ATTRIBUTE node.  The procedure checkattr() is generated from
 *	attribute.tab by the auxilliary proc mk_checkattr.c
 *
 */


#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdbool.h>

#ifdef WIN32
# define WIN32_LEAN_AND_MEAN
# include "Windows.h"
#else
# include <libgen.h> // for basename
# include <unistd.h>
#endif

#include "parse2.h"

#define	FALSE	0
#define	TRUE	1
#define	UNDEF	(-1)

#define	MAXFIELDS	50
	static	char	*fields[ MAXFIELDS ];
	int	n_fields;

	/* symbols (#defines) from nab.h	*/


#define	TVTAB_SIZE	30
static	DEF_T	tvtab[ TVTAB_SIZE ];
static	int	n_tvtab;
static	int	v_error;
static	int	l_utype;

#define	CVTAB_SIZE	10
static	DEF_T	cvtab[ CVTAB_SIZE ];
static	int	n_cvtab;

#define	KVTAB_SIZE	10
static	DEF_T	kvtab[ KVTAB_SIZE ];
static	int	n_kvtab;

#define	AVTAB_SIZE	5
static	DEF_T	avtab[ AVTAB_SIZE ];
static	int	n_avtab;

#define	SYMTAB_SIZE	300
static	DEF_T	symtab[ SYMTAB_SIZE ];
static	int  n_symtab;

	/* rule files reserved words	*/

static	RESWORD	rwords[] = {
	{ "akind", TOK_AKIND },
	{ "class", TOK_CLASS },
	{ "code", TOK_CODE },
	{ "if", TOK_IF },
	{ "kind", TOK_KIND },
	{ "left", TOK_LEFT },
	{ "operator", TOK_OPERATOR },
	{ "output", TOK_OUTPUT },
	{ "prec", TOK_PREC },
	{ "print", TOK_PRINT },
	{ "right", TOK_RIGHT },
	{ "sym", TOK_SYM },
	{ "type", TOK_TYPE },
	{ "use", TOK_USE }
	
};
static	int	n_rwords = sizeof( rwords ) / sizeof( RESWORD );

	/* attribute information:	*/
#define	ATTAB_SIZE	50
static	ATREC_T	attab[ ATTAB_SIZE ];
static	int	n_attab;

static	char	*rfname; // name of current file being read
static	int	rfeof = FALSE;
static	char	rfline[ 256 ] = "";
static	char	*rfp = rfline;

static	int	readahead = FALSE;
static	int	ra_tok;
static	int	ra_tokval;
static	char	ra_tokstr[ 256 ];

static	int	tok;
static	int	tokval;
static	char	tokstr[ 256 ];

#define	KIND_SIZE	10
static	int	a_kind[ KIND_SIZE ];
static	int	l_kind[ KIND_SIZE ];
static	int	r_kind[ KIND_SIZE ];
static	int	o_ktok[ KIND_SIZE ];
static	int	o_kval[ KIND_SIZE ];
static	int	o_kqtok[ KIND_SIZE ];
static	int	o_kqval[ KIND_SIZE ];

#define	TTAB_SIZE	5
static	TTAB_T	ttabs[ TTAB_SIZE ];

static	RULE_T	**rules;

	/* HTML versions of the rules:	*/
#define	TWID1	100
#define	TWID2	400
#define	SELSIZE	6

#define	FRFILE	"frames.html"
#define	CTLFILE	"control.html"
#define	CFSIZE	30
static	char	r1fname[ 256 ] = "";

// variables which are set in main() but need to be global
char* nab_h_path;
char* progname;
int next_type_index;

int	main( int argc, char *argv[] )
{
	FILE	*fp;
	int	err = 0, err_count = 0;

	progname = argv[ 0 ];

	const char* usage_string = "usage: %s --amber-data-dir <path to amber data dir> --y-tab-h <path to nab2c's y.tab.h> --attributes-tab <path to atrfile.tab>\
 [--html-output-dir <dir to write HTML files to>] --nab-h-path <path to nab.h> -o <output file path> <rule-files...> \n";

	// we need at least 11 arguments, plus the argv[0] program name
	if( argc < 12 ){
		fprintf( stderr, usage_string, progname );
		exit( 1 );
	}

	// variables that will be set by the arguments
	char* amber_data_dir_path = NULL;
	char* y_tab_h_path = NULL;
	char* attributes_tab_path = NULL;
	nab_h_path = NULL;
	char* output_file_path = NULL;

	bool generate_html = false;
	char* html_output_path = NULL;

	int first_rules_file_in_argv_index;

	for(int arg_index = 1; arg_index < argc; ++arg_index)
	{
        fprintf( stderr, "%d  %s\n", arg_index, argv[arg_index] );
		if(strcmp(argv[arg_index], "--amber-data-dir") == 0)
		{
			amber_data_dir_path = argv[++arg_index];
		}
		else if(strcmp(argv[arg_index], "--y-tab-h") == 0)
		{
			y_tab_h_path = argv[++arg_index];
            fprintf( stderr, "%d  %s\n", arg_index, argv[arg_index] );
		}
		else if(strcmp(argv[arg_index], "--attributes-tab") == 0)
		{
			attributes_tab_path = argv[++arg_index];
		}
		else if(strcmp(argv[arg_index], "--nab-h-path") == 0)
		{
			nab_h_path = argv[++arg_index];
		}
		else if(strcmp(argv[arg_index], "-o") == 0)
		{
			output_file_path = argv[++arg_index];
		}
		else if(strcmp(argv[arg_index], "--html-output-dir") == 0)
		{
			html_output_path = argv[++arg_index];
			generate_html = true;
		}
		else
		{
			// since we have seen something that isn't an option, then we must be at the start of the source file list
			first_rules_file_in_argv_index = arg_index;
			break;
		}
	}

	// check for missing options
	if(y_tab_h_path == NULL || attributes_tab_path == NULL || nab_h_path == NULL || output_file_path == NULL)
	{
		fprintf( stderr, "error: missing arguments\n" );
		fprintf( stderr, "y_tab_h_path: %s\n", y_tab_h_path );
		fprintf( stderr, "attributes_tab_path: %s\n", attributes_tab_path );
		fprintf( stderr, "nab_h_path: %s\n", nab_h_path );
		fprintf( stderr, "output_file_path: %s\n", nab_h_path );
		fprintf( stderr, usage_string, progname );
		exit( 1 );
	}


	if( ( fp = fopen( y_tab_h_path, "r" ) ) == NULL )
	{
		fprintf( stderr, "%s: can't read symbol file %s\n",
			progname, y_tab_h_path );
		exit( 1 );
	}

	if( ( n_symtab = getsyms( fp, symtab ) ) == 0 )
	{
		fprintf( stderr, "%s: no symbols found in symbol header %s.\n", progname, y_tab_h_path );
		exit( 1 );
	}

	fclose( fp );

	fprintf( stderr, "%s: %3d symbols.\n", progname, n_symtab );

	rules = ( RULE_T ** )malloc( n_symtab * sizeof( RULE_T * ) );
	if( rules == NULL )
	{
		fprintf( stderr, "%s: can't allocate rules.\n", progname );
		exit( 1 );
	}
	
	if( ( fp = fopen( nab_h_path, "r" ) ) == NULL )
	{
		fprintf( stderr, "%s: can't read type/class/kind file %s\n",
			progname, nab_h_path );
		exit( 1 );
	}

	if( gettcka( fp, &l_utype, &v_error, &n_tvtab, tvtab, &n_cvtab, cvtab, &n_kvtab, kvtab, &n_avtab, avtab ) )
	{
		exit( 1 );
	}

	fclose( fp );
	fprintf( stderr, "%s: %3d types, lastuser = %d, v_error = %d.\n", progname, n_tvtab, l_utype, v_error );
	fprintf( stderr, "%s: %3d classes.\n", progname, n_cvtab ); 
	fprintf( stderr, "%s: %3d kinds.\n", progname, n_kvtab ); 
	fprintf( stderr, "%s: %3d accesses.\n", progname, n_avtab ); 

	/* experimental stuff: do attribute checkin here,	*/
	/* create the file checkattr() from the file ATRFILE	*/

	if( ( fp = fopen( attributes_tab_path, "r" ) ) == NULL )
	{
		fprintf( stderr, "%s: can't read attibute file %s.\n", progname, attributes_tab_path );
		exit( 1 );
	}
	rfname = attributes_tab_path;
	rfeof = FALSE;
	err = parse_attrs( fp, &n_attab, attab );
	fclose( fp );

	if( n_attab == 0 )
	{
		fprintf( stderr, "%s: No attributes found in %s.\n", progname, attributes_tab_path );
		exit( 1 );
	}
	else
	{
		fprintf( stderr, "%s: %3d attributes.\n", progname, n_attab ); 
	}

	int num_rules = 0;
	static	RULE_T	c_rule; // we re-use one rule struct to hold each file in turn
	memset(&c_rule, 0, sizeof(RULE_T));

	next_type_index = 0;

	for(int rule_index = 0; rule_index < (argc - first_rules_file_in_argv_index) ; ++rule_index )
	{
		err = 0;

		int arg_index = rule_index + first_rules_file_in_argv_index;

		if( ( fp = fopen( argv[ arg_index ], "r" ) ) == NULL )
		{
			fprintf( stderr, "%s: can't read rule-file %s.\n", progname, argv[ arg_index ] );
			err_count++;
			continue;
		}
		else
		{
			rfname = argv[ arg_index ];
			rfeof = FALSE;
		}


		clear_rule( &c_rule );
		err = parse_rule( fp, &c_rule );
		fclose( fp );
		if( err )
		{
			err_count++;
			continue;
		}
		if((err = check_rule( &c_rule )))
		{
			err_count++;
			continue;
		}

		rules[ rule_index ] = save_rule( &c_rule );
		if( rules[ rule_index ] == NULL )
		{
			exit( 1 );	/* end of the world!	*/
		}

		if( ( fp = fopen( argv[ arg_index ], "wb" ) ) == NULL )
		{
			fprintf( stderr, "%s: can't rewrite rule-file %s.\n", progname, argv[ arg_index ] );
			err_count++;
			continue;
		}
		fprint_rule( fp, &c_rule );
		fclose( fp );

		if(generate_html)
		{
			if((err = hprint_rule( argv[ arg_index ], &c_rule, html_output_path )))
			{
				err_count++;
			}
		}

		++num_rules;
	}

	fprintf( stderr, "%s: %3d rules.\n", progname, num_rules );

	if(generate_html)
	{
		if( mk_frameset( html_output_path ) )
		{
			err_count++;
		}

		if( mk_controlfile( html_output_path ) )
		{
			err_count++;
		}
	}

	if(err_count == 0)
	{
		mk_checkexpr(output_file_path, num_rules);
	}

	fprintf( stderr, "%s: %3d used symbols.\n", progname, countudefs( n_symtab, symtab ) );
	fprintf( stderr, "%s: %3d used types.\n", progname, countudefs( n_tvtab, tvtab ) );
	fprintf( stderr, "%s: %3d used classes.\n", progname, countudefs( n_cvtab, cvtab ) );
	fprintf( stderr, "%s: %3d used kinds.\n", progname, countudefs( n_kvtab, kvtab ) );
	fprintf( stderr, "%s: %3d errors.\n", progname, err_count );

	exit( err_count > 0 );
}

static	int	getsyms( FILE *fp, DEF_T symtab[] )
{
	char	line[ 256 ];
	DEF_T	*sp;

	for( sp = symtab; fgets( line, sizeof( line ), fp ); sp++ ){
		sscanf( line, "     %s = %d,", sp->d_name, &sp->d_val );
		sp->d_used = FALSE;
	}
	return( (int)(sp - symtab) );
}

static	int	gettcka( FILE* fp, int *l_utype, int *v_error, int *n_tvtab, DEF_T tvtab[],
	int *n_cvtab, DEF_T cvtab[], int *n_kvtab, DEF_T kvtab[], int *n_avtab, DEF_T avtab[] )
{
	char	line[ 256 ];
	DEF_T	*tvp, *tvp1, *cvp, *kvp, *avp;
	int	found, err;
	int	f, i;

	err = FALSE;
	*l_utype = *n_tvtab = *n_cvtab = *n_kvtab = *n_avtab = 0;
	*v_error = UNDEF;
	tvp = tvtab;
	cvp = cvtab;
	kvp = kvtab;
	avp = avtab;
	while( fgets( line, sizeof( line ), fp ) ){
		n_fields = split( line, fields, " \t\n" );
		if( n_fields < 3 ){ 
			for( f = 0; f < n_fields; f++ )
				free( fields[ f ] );
			continue;
		}
		if( strcmp( fields[ 0 ], "#define" ) ){
			for( f = 0; f < n_fields; f++ )
				free( fields[ f ] );
			continue;
		}
		if( !strncmp( fields[ 1 ], "T_", 2 ) ){
			if( !strcmp( fields[ 1 ], "T_LASTUSER" ) ){
				tvp1 = tvtab;
				found = 0;
				for( i = 0; i < *n_tvtab; i++, tvp1++ ){
					if(!strcmp( tvp1->d_name, fields[2] )){
						found = 1;
						break;
					}
				}
				if( !found ){
					fprintf( stderr, 
			"%s: gettcka: T_LASTUSER value %s not defined.\n",
						progname, fields[2] );
					return( TRUE );
				}else
					*l_utype = tvp1->d_val;
			}else{
				strcpy( tvp->d_name, fields[ 1 ] );
				tvp->d_val = atoi( fields[ 2 ] );
				tvp->d_used = 0;
				tvp++;
				( *n_tvtab )++;
			}
		}else if( !strncmp( fields[ 1 ], "C_", 2 ) ){
			strcpy( cvp->d_name, fields[ 1 ] );
			cvp->d_val = atoi( fields[ 2 ] );
			cvp->d_used = 0;
			cvp++;
			( *n_cvtab )++;
		}else if( !strncmp( fields[ 1 ], "K_", 2 ) ){
			strcpy( kvp->d_name, fields[ 1 ] );
			kvp->d_val = atoi( fields[ 2 ] );
			kvp->d_used = 0;
			kvp++;
			( *n_kvtab )++;
		}else if( !strncmp( fields[ 1 ], "A_", 2 ) ){
			strcpy( avp->d_name, fields[ 1 ] );
			avp->d_val = atoi( fields[ 2 ] );
			avp->d_used = 0;
			avp++;
			( *n_avtab )++;
		}
		for( f = 0; f < n_fields; f++ )
			free( fields[ f ] );
	}
	if( *n_tvtab == 0 ){
		fprintf( stderr,
			"%s: gettcka: no types (T_* defines) found.\n",
			progname );
		err = TRUE;
	}
	for( tvp = tvtab, i = 0; i < *n_tvtab; i++, tvp++ ){
		if( tvp->d_val != i ){
			fprintf( stderr,
		"%s: gettcka: type %s has value %d should have value %d.\n",
				progname, tvp->d_name, tvp->d_val, i );
			err = TRUE;
		}
		if( !strcmp( tvp->d_name, "T_ERROR" ) )
			*v_error = tvp->d_val;
	}
	if( *v_error == UNDEF ){
		fprintf( stderr,
			"%s: gettcka: no #define for T_ERROR.\n", progname );
		err = TRUE;
	}
	if( *n_cvtab == 0 ){
		fprintf( stderr, "gettcka: no classes (C_* defines) found.\n" );
		err = TRUE;
	}
	for( cvp = cvtab, i = 0; i < *n_cvtab; i++, cvp++ ){
		if( cvp->d_val != i ){
			fprintf( stderr,
		"%s: gettcka: class %s has value %d should have value %d.\n",
				progname, cvp->d_name, cvp->d_val, i );
			err = TRUE;
		}
	}
	if( *n_kvtab == 0 ){
		fprintf( stderr, "gettcka: no kinds (K_* defines) found.\n" );
		err = TRUE;
	}
	for( kvp = kvtab, i = 0; i < *n_kvtab; i++, kvp++ ){
		if( kvp->d_val != i ){
			fprintf( stderr,
		"%s: gettcka: kind %s has value %d should have value %d.\n",
				progname, kvp->d_name, kvp->d_val, i );
			err = TRUE;
		}
	}
	if( *n_avtab == 0 ){
		fprintf( stderr, "gettcka: no accesses (A_* defines) found.\n" );
		err = TRUE;
	}
	for( avp = avtab, i = 0; i < *n_avtab; i++, avp++ ){
		if( avp->d_val != i ){
			fprintf( stderr,
		"%s: gettcka: access %s has value %d should have value %d.\n",
				progname, avp->d_name, avp->d_val, i );
			err = TRUE;
		}
	}
	return( err );
}

static	DEF_T	*finddef( char name[], int n_tab, DEF_T tab[] )
{
	int	i;
	DEF_T	*dp;

	for( dp = tab, i = 0; i < n_tab; i++, dp++ ){
		if( !strcmp( name, dp->d_name ) )
			return( dp );
	}
	return( NULL );
}

static	DEF_T	*finddef_withval( int val, int n_tab, DEF_T tab[] )
{
	int	i;
	DEF_T	*dp;

	for( dp = tab, i = 0; i < n_tab; i++, dp++ ){
		if( val == dp->d_val )
			return( dp );
	}
	return( NULL );
}

static	int	countudefs( int n_tab, DEF_T tab[] )
{
	int	i, cnt;
	DEF_T	*dp;

	for( cnt = 0, dp = tab, i = 0; i < n_tab; i++, dp++ )
		cnt = dp->d_used ? cnt + 1 : cnt;
	return( cnt );
}

static	int	defcmp( void const *v1, void const *v2 )
{
	DEF_T const	* const *d1 = v1;
	DEF_T const	* const *d2 = v2;
	return( strcmp( ( *d1 )->d_name, ( *d2 )->d_name ) );
}

static	int	parse_attrs( FILE *fp, int *n_attab, ATREC_T attab[] )
{
	int	err;
	int	na;
	ATREC_T	*ap;
	char	*sp;
	DEF_T	*dp;
	int	val, t_idx;

	err = FALSE;
	*n_attab = 0;
	for( na = 0; gettoken( fp ); ){
		if( tok == TOK_IDENT ){
			ap = &attab[ na ]; 
			na++;
			sp = ( char * )malloc( strlen( tokstr ) + 1 );
			strcpy( sp, tokstr );
			ap->a_name = sp;
			ap->a_type = 0;
			ap->a_class = 0;
			ap->a_kind = 0;
			ap->a_in = 0;		/* set of types	*/
			ap->a_access = 0;
			gettoken( fp );
			if( tok == TOK_IDENT ){	
				dp = finddef( tokstr, n_tvtab, tvtab );
				if( dp == NULL ){
					fprintf( stderr,
						"%s: undefined type %s.\n", 
						rfname, tokstr );
					err = TRUE;
					goto skip;
				}else
					ap->a_type = dp->d_val;
			}else{
				fprintf( stderr,
					"%s: attribute %s -- bad type value.\n",
					rfname, ap->a_name );
				err = TRUE;
				goto skip;
			}
			gettoken( fp );
			if( tok == TOK_IDENT ){	
				dp = finddef( tokstr, n_cvtab, cvtab );
				if( dp == NULL ){
					fprintf( stderr,
						"%s: undefined class %s.\n", 
						rfname, tokstr );
					err = TRUE;
					goto skip;
				}else
					ap->a_class = dp->d_val;
			}else{
				fprintf( stderr,
				"%s: attribute %s -- bad class value.\n",
					rfname, ap->a_name );
				err = TRUE;
				goto skip;
			}
			gettoken( fp );
			if( tok == TOK_IDENT ){	
				dp = finddef( tokstr, n_kvtab, kvtab );
				if( dp == NULL ){
					fprintf( stderr,
						"%s: undefined kind %s.\n", 
						rfname, tokstr );
					err = TRUE;
					goto skip;
				}else
					ap->a_kind = dp->d_val;
			}else{
				fprintf( stderr,
					"%s: attribute %s -- bad kind value.\n",
					rfname, ap->a_name );
				err = TRUE;
				goto skip;
			}
			err |= getsetvalue( fp, "attribute in vals",
				n_tvtab, tvtab, &val, &t_idx );
			ap->a_in = val;
			if( tok == TOK_IDENT ){
				dp = finddef( tokstr, n_avtab, avtab );
				if( dp == NULL ){
					fprintf( stderr,
						"%s: undefined access %s.\n", 
						rfname, tokstr );
					err = TRUE;
					goto skip;
				}else
					ap->a_access = dp->d_val;
			}else{
				fprintf( stderr,
					"%s: attribute %s -- bad access value.\n",
					rfname, ap->a_name );
				err = TRUE;
				goto skip;
			}
		}else if( tok != TOK_COMMENT && tok != TOK_NL ){
			fprintf( stderr, "%s: syntax error '%s'.\n",
				rfname, tokstr );
		}
	skip : ;
		skiptonl( fp );
	}
	*n_attab = na;
	return( err );
}

static	int	parse_rule( FILE *fp, RULE_T *rp )
{
	int	err;
	
	err = FALSE;
	while( gettoken( fp ) ){
		if( tok == TOK_OPERATOR ){
			skiptonl( fp );
			err |= operator( fp, rp );
		}else if( tok == TOK_CLASS ){
			skiptonl( fp );
			err |= class( fp, rp );
		}else if( tok == TOK_KIND ){
			skiptonl( fp );
			err |= kind( fp, rp );
		}else if( tok == TOK_TYPE ){
			err |= type( fp, rp );
		}else if( tok == TOK_CODE ){
			skiptonl( fp );
		}else if( tok == TOK_ERROR ){
			fprintf( stderr,
				"%s: syntax error: tokstr = '%s'\n",
				rfname, tokstr );
			skiptonl( fp );
			err = TRUE;
		}
	}
	return( err );
}

static	int	operator( FILE *fp, RULE_T *rp )
{
	int	err;
	char	*sp;
	DEF_T	*dp;
	int	g_sym, g_prec, g_print;

	err = FALSE;
	g_sym = g_prec = g_print = FALSE;
	while( gettoken( fp ) ){
		if( tok == TOK_SYM ){
			if( g_sym ){
				fprintf( stderr,
				"%s: operator part: mulitple sym entries.\n",
					rfname );
				err = TRUE;
				skiptonl( fp );
				continue;
			}else
				g_sym = TRUE;
			gettoken( fp );
			if( tok == TOK_EQUAL ){
				gettoken( fp );
				if( tok == TOK_IDENT ){
					if((dp=finddef(tokstr,n_symtab,symtab))){
						rp->r_sym = dp->d_val;
						dp->d_used = TRUE;
					}else{
						fprintf( stderr,
				"%s: operator part: undefined sym %s.\n",
							rfname, tokstr );
						return( TRUE );
					}
				}else
					err = TRUE;
			}else
				err = TRUE;
			skiptonl( fp );
		}else if( tok == TOK_PREC ){
			if( g_prec ){
				fprintf( stderr,
				"%s: operator part: multiple prec entries.\n",
					rfname );
				err = TRUE;
				skiptonl( fp );
				continue;
			}else
				g_prec = TRUE;
			gettoken( fp );
			if( tok == TOK_EQUAL ){
				gettoken( fp );
				if( tok == TOK_INT )
					rp->r_prec = tokval;
				else
					err = TRUE;
			}else
				err = TRUE;
			skiptonl( fp );
		}else if( tok == TOK_PRINT ){
			if( g_print ){
				fprintf( stderr,
				"%s: operator part: mulitple print entries.\n",
					rfname );
				err = TRUE;
				skiptonl( fp );
				continue;
			}else
				g_print = TRUE;
			gettoken( fp );
			if( tok == TOK_EQUAL ){
				gettoken( fp );
				if( tok == TOK_STRING ){
					sp = ( char * )malloc(
						strlen( tokstr ) + 1 );
					strcpy( sp, tokstr );
					rp->r_print = sp;
				}else
					err = TRUE;
			}else
				err = TRUE;
			skiptonl( fp );
		}else if( tok != TOK_NL && tok != TOK_COMMENT ){
			ungettoken();
			break;
		}
	}
	if( !g_sym ){
		fprintf( stderr,
			"%s: operator part: no sym entry.\n", rfname );
		err = TRUE;
	}
	if( !g_prec ){
		fprintf( stderr,
			"%s: operator: no prec entry.\n", rfname );
		err = TRUE;
	}
	if( !g_prec ){
		fprintf( stderr,
			"%s: operator: no print entry.\n", rfname );
		err = TRUE;
	}
	return( err );
}

static	int	class( FILE *fp, RULE_T *rp )
{
	int	err;
	int	val, t_idx;
	int	otok, oval, qtok, qval;
	int	g_left, g_right, g_output;
	
	err = FALSE;
	g_left = g_right = g_output = FALSE;
	while( gettoken( fp ) ){
		if( tok == TOK_LEFT ){
			if( g_left ){
				fprintf( stderr,
				"%s: class part: multiple left entries.\n",
					rfname );
					err = TRUE;
				skiptonl( fp );
				continue;
			}else
				g_left = TRUE;
			gettoken( fp );
			if( tok == TOK_EQUAL ){
				err |= getsetvalue( fp, "class part: left",
					n_cvtab, cvtab, &val, &t_idx );
				if( !err )
					rp->r_lclass = val;
			}else
				err = TRUE;
			skiptonl( fp );
		}else if( tok == TOK_RIGHT ){
			if( g_right ){
				fprintf( stderr,
			"%s: class part: multiple right class entries.\n",
					rfname );
					err = TRUE;
				skiptonl( fp );
				continue;
			}else
				g_right = TRUE;
			gettoken( fp );
			if( tok == TOK_EQUAL ){
				err |= getsetvalue( fp, "class part: right",
					n_cvtab, cvtab, &val, &t_idx );
				if( !err )
					rp->r_rclass = val;
			}else
				err = TRUE;
			skiptonl( fp );
		}else if( tok == TOK_OUTPUT ){
			if( g_output ){
				fprintf( stderr,
			"%s: class part: multiple output class entries.\n",
					rfname );
					err = TRUE;
				skiptonl( fp );
				continue;
			}else
				g_output = TRUE;
			gettoken( fp );
			if( tok == TOK_EQUAL ){
				err |= getovalue( fp, "class part: output",
					n_cvtab, cvtab,
					&otok, &oval, &qtok, &qval );
				if( !err ){
					rp->r_octok = otok;
					rp->r_ocval = oval;
				}
			}else
				err = TRUE;
			skiptonl( fp );
		}else if( tok != TOK_NL && tok != TOK_COMMENT ){
			ungettoken();
			break;
		}
	}
	if( !g_left && !g_right ){
		fprintf( stderr,
			"%s: class part: no input class entries.\n",
			rfname );
		err = TRUE;
	}
	if( !g_output ){
		fprintf( stderr,
			"%s: class part: no output class entry.\n",
			rfname );
		err = TRUE;
	}
	return( err );
}

static	int	kind( FILE *fp, RULE_T *rp )
{
	int	err, mkind;
	int	val, t_idx;
	int	otok, oval;
	int	g_akind, g_left, g_right, g_output;
	int	qtok, qval;
	int	i, j;
	
	err = FALSE;
	mkind = FALSE;
	g_akind = g_left = g_right = g_output = FALSE;
	while( gettoken( fp ) ){
		if( tok == TOK_AKIND ){
			g_akind = TRUE;
			gettoken( fp );
			if( tok == TOK_EQUAL ){
				err |= getsetvalue( fp, "kind part: array",
					n_kvtab, kvtab, &val, &t_idx );
				if( !err ){
					rp->r_akind[rp->r_nakind] = val;
					rp->r_nakind++;
				}
			}else
				err = TRUE;
			skiptonl( fp );
		}else if( tok == TOK_LEFT ){
			g_left = TRUE;
			gettoken( fp );
			if( tok == TOK_EQUAL ){
				err |= getsetvalue( fp, "kind part: left",
					n_kvtab, kvtab, &val, &t_idx );
				if( !err ){
					rp->r_lkind[rp->r_nlkind] = val;
					rp->r_nlkind++;
				}
			}else
				err = TRUE;
			skiptonl( fp );
		}else if( tok == TOK_RIGHT ){
			g_right = TRUE;
			gettoken( fp );
			if( tok == TOK_EQUAL ){
				err |= getsetvalue( fp, "kind part: right",
					n_kvtab, kvtab, &val, &t_idx );
				if( !err ){
					rp->r_rkind[rp->r_nrkind] = val;
					rp->r_nrkind++;
				}
			}else
				err = TRUE;
			skiptonl( fp );
		}else if( tok == TOK_OUTPUT ){
			g_output = TRUE;
			gettoken( fp );
			if( tok == TOK_EQUAL ){
				err |= getovalue( fp, "kind part: output",
					n_kvtab, kvtab, &otok, &oval,
					&qtok, &qval );
				if( !err ){
					rp->r_oktok[rp->r_nokind] = otok;
					rp->r_okval[rp->r_nokind] = oval;
					rp->r_okqtok[rp->r_nokind] = qtok;
					rp->r_okqval[rp->r_nokind] = qval;
					rp->r_nokind++;
				}
			}else
				err = TRUE;
			skiptonl( fp );
		}else if( tok != TOK_NL && tok != TOK_COMMENT ){
			ungettoken();
			break;
		}
	}
	if( g_akind ){
		if( !g_left ){
			fprintf( stderr, "%s: kind part: no left entry.\n",
				rfname );
			err = TRUE;
		}
		if( !g_right ){
			fprintf( stderr, "%s: kind part: no right entry.\n",
				rfname );
			err = TRUE;
		}
	}else if( !g_left && !g_right ){
		fprintf( stderr, "%s: kind part: no left/right entries.\n",
			rfname );
		err = TRUE;
	}
	if( !g_output ){
		fprintf( stderr,
			"%s: kind part: no output entry.\n",
			rfname );
		err = TRUE;
	}
	for( i = 0; i < rp->r_nakind - 1; i++ ){
		for( j = i + 1; j < rp->r_nakind; j++ ){
			if( rp->r_akind[i] & rp->r_akind[j] ){
				err = TRUE;
				fprintf( stderr, 
		"%s: kind part: akind entries %d and %d have common members.\n",
					rfname, i+1, j+1 );
			}
		}
	}
	for( i = 0; i < rp->r_nlkind - 1; i++ ){
		for( j = i + 1; j < rp->r_nlkind; j++ ){
			if( rp->r_lkind[i] & rp->r_lkind[j] ){
				err = TRUE;
				fprintf( stderr, 
		"%s: kind part: left entries  %d and %d have common members.\n",
					rfname, i+1, j+1 );
			}
		}
	}
	for( i = 0; i < rp->r_nrkind - 1; i++ ){
		for( j = i + 1; j < rp->r_nrkind; j++ ){
			if( rp->r_rkind[i] & rp->r_rkind[j] ){
				err = TRUE;
				fprintf( stderr, 
		"%s: kind part: right entries %d and %d have common members.\n",
					rfname, i+1, j+1 );
			}
		}
	}
	if( rp->r_nakind > 1 )
		mkind = TRUE;
	if( rp->r_nlkind > 1 ){
		if( mkind ){
			fprintf( stderr,
	"%s: kind part: only 1 of akind, left, right can have >1 entry.\n",
				rfname );
			err = TRUE;
		}
		mkind = TRUE;
	}
	if( rp->r_nrkind > 1 ){
		if( mkind ){
			fprintf( stderr,
	"%s: kind part: only 1 of akind, left, right can have >1 entry.\n",
				rfname );
			err = TRUE;
		}
		mkind = TRUE;
	}
	return( err );
}

static	int	getsetvalue( FILE *fp, char ctxt[], int n_tab, DEF_T tab[], int *val, int *t_idx )
{
	int	cnt;
	DEF_T	*dp;
	
	*val = 0;
	*t_idx = 1;
	gettoken( fp );
	if( tok == TOK_IDENT ){
		if((dp = finddef( tokstr, n_tab, tab ))){
			*val = 1 << dp->d_val;
			dp->d_used = TRUE;
		}else{
			fprintf( stderr, "%s: %s entry %s undefined.\n",
				rfname, ctxt, tokstr );
			return( TRUE );
		}
		gettoken( fp );
	}else if( tok == TOK_LBRACE ){
		for( cnt = 0; ; ){
			gettoken( fp );
			if( tok != TOK_IDENT ){
				fprintf( stderr,
				"%s: %s set -- syntax error.\n",
					rfname, ctxt );
				return( TRUE );
			}
			if( (dp = finddef( tokstr, n_tab, tab )) ){
				*val |= 1 << dp->d_val;
				dp->d_used = TRUE;
				cnt++;
			}else{
				fprintf( stderr, "%s: %s entry %s undefined.\n",
					rfname, ctxt, tokstr );
				return( TRUE );
			}
			gettoken( fp );
			if( tok == TOK_COMMA )
				continue;
			else if( tok == TOK_RBRACE ){
				gettoken( fp );
				break;
			}else{
				fprintf( stderr,
				"%s: %s set -- syntax error.\n",
					rfname, ctxt );
				return( TRUE );
			}
		}
		if( cnt == 0 ){
			fprintf( stderr,
				"%s: %s set is empty.\n",
				rfname, ctxt );
			return( TRUE );
		}
	}
	if( tok == TOK_USE ){
		gettoken( fp );
		if( tok == TOK_TYPE ){
			gettoken( fp );
			if( tok == TOK_EQUAL ){
				gettoken( fp );
				if( tok == TOK_INT )
					*t_idx = tokval;
				else{
					fprintf( stderr,
					"%s: %s use clause -- syntax error.\n",
						rfname, ctxt );
					return( TRUE );
				}
			}else{
				fprintf( stderr,
					"%s: %s use clause -- syntax error.\n",
					rfname, ctxt );
				return( TRUE );
			}
		}else{
			fprintf( stderr, "%s: %s use clause -- syntax error.\n",
				rfname, ctxt );
			return( TRUE );
		}
/*
	}else if( tok != TOK_NL && tok != TOK_EOF ){
		fprintf( stderr, "%s: %s set -- syntax error.\n",
			rfname, ctxt );
		return( TRUE );
*/
	}
	return( FALSE );
}

static	int	getovalue( FILE *fp, char ctxt[], int n_tab, DEF_T tab[], int *otok, int *oval, int *qtok, int *qval )
{
	int	err;
	DEF_T	*dp;
	int	t_idx;

	err = FALSE;
	gettoken( fp );
	if( tok == TOK_IDENT ){
		*otok = tok;
		if( (dp = finddef( tokstr, n_tab, tab )) ){
			dp->d_used = TRUE;
			*oval = dp->d_val;
		}else{
			fprintf( stderr, "%s: %s entry %s undefined.\n",
				rfname, ctxt, tokstr );
			return( TRUE );
		}
		gettoken( fp );
	}else if( tok == TOK_LEFT || tok == TOK_RIGHT ){
		*otok = tok;
		*oval = 0;
		gettoken( fp );
	}
	*qtok = TOK_NL;
	*qval = 0;
	if( tok == TOK_IF ){
		gettoken( fp );
		if( tok == TOK_LEFT || tok == TOK_RIGHT ){
			*qtok = tok;
			gettoken( fp );
			if( tok == TOK_EQUAL ){
				err |= getsetvalue( fp, ctxt, n_tab, tab,
					qval, &t_idx );
			}else
				err = TRUE;
		}else
			err = TRUE;
		skiptonl( fp );
	}
	if( tok != TOK_NL && tok != TOK_EOF ){
		fprintf( stderr, "%s: %s set -- syntax error.\n",
			rfname, ctxt );
		err = TRUE;
	}
	return( err );
}

static	int	type( FILE *fp, RULE_T *rp )
{
	int	err;
	int	t, tnum;
	int	n_rows;
	TROW_T	*trows[ TVTAB_SIZE ];
	char	*ip;
	TROW_T	*trp, **trp1;

	err = FALSE;
	tnum = 1;
	gettoken( fp );
	if( tok == TOK_EQUAL ){
		gettoken( fp );
		if( tok == TOK_INT ){
			tnum = tokval;
		}
	}else if( tok != TOK_NL ){
		fprintf( stderr, "%s: type: syntax error in type header.\n",
			rfname );
		err = TRUE;
	}
	skiptonl( fp );
	for( n_rows = 0; gettoken( fp ); ){
		if( tok == TOK_IDENT ){
			trp = ( TROW_T * )malloc( sizeof( TROW_T ) );
			ip = ( char * )malloc( strlen( tokstr ) + 1 );
			strcpy( ip, tokstr );
			trp->t_name = ip;
			trp->t_nval = UNDEF;
			trp->t_nrvals = 0;
			trp->t_rvals = NULL;
			gettoken( fp );
			if( tok == TOK_EQUAL ){
				gettoken( fp );
				err |= getvecvalue( fp, trp );
				trows[ n_rows ] = trp;
				n_rows++;
			}
			skiptonl( fp );
		}else if( tok != TOK_NL && tok != TOK_COMMENT ){
			ungettoken();
			break;
		}
	}

	trp1 = ( TROW_T ** )malloc( n_rows * sizeof( TROW_T * ) );
	for( t = 0; t < n_rows; t++ )
		trp1[ t ] = trows[ t ];
	ttabs[ tnum - 1 ].t_index = next_type_index;
	++next_type_index;
	ttabs[ tnum - 1 ].t_nrows = n_rows;
	ttabs[ tnum - 1 ].t_rows = trp1;
	rp->r_ntypes++;

	return( err );
}

static	int	getvecvalue( FILE *fp, TROW_T *trp )
{
	int	err;
	int	i, n_col;
	int	trow[ TVTAB_SIZE ];
	int	*rp;

	err = FALSE;
	if( tok != TOK_LBRACK )
		return( TRUE );
	gettoken( fp );
	for( n_col = 0; tok != TOK_RBRACK; ){
		if( tok == TOK_INT ){
			trow[ n_col ] = tokval;
			n_col++;
			gettoken( fp );
			if( tok == TOK_COMMA )
				gettoken( fp );
		}
	}
	if( n_col == 0 ){
		fprintf( stderr, "%s: getvecvalue: empty vector.\n", rfname );
		err = TRUE;
	}else{
		trp->t_nrvals = n_col;
		rp = ( int * )malloc( n_col * sizeof( int ) );
		for( i = 0; i < n_col; i++ )
			rp[ i ] = trow[ i ];
		trp->t_rvals = rp;
	}
	return( err );
}

static	int	gettoken( FILE *fp )
{
	char	*tsp;
	int	i;

	if( readahead ){
		tok = ra_tok;
		tokval = ra_tokval;
		strcpy( tokstr, ra_tokstr );
		readahead = FALSE;
		return( TRUE );
	}

	tok = TOK_EOF;
	tokval = 0;
	*tokstr = '\0';
	if( rfeof )
		return( FALSE );
	else if( *rfp == '\0' ){
		if( !fgets( rfline, sizeof( rfline ), fp ) ){
			rfeof = TRUE;
			return( FALSE );
		}else
			rfp = rfline;
	}
	while( *rfp == ' ' || *rfp == '\t' )
		rfp++;
	tsp = tokstr;
	if( isalpha( *rfp ) ){
		tok = TOK_IDENT;
		*tsp++ = *rfp++; 
		for( ; isalnum( *rfp ) || *rfp == '_' ; rfp++ )
			*tsp++ = *rfp;
		*tsp = '\0';
		for( i = 0; i < n_rwords; i++ ){
			if( !strcmp( tokstr, rwords[i].r_name ) ){
				tok = rwords[i].r_token;
				break;
			}
		}
	}else if( *rfp == '"' ){
		tok = TOK_STRING;
		rfp++;
		while( *rfp != '\0' && *rfp != '"' && *rfp != '\n' ){
			if( *rfp == '\\' ){
				if( rfp[1] == '\n' || rfp[1] == '\0' ){
					fprintf( stderr,
					"%s: gettoken: Unterminated string.",
						rfname );
					break;
				}else
					rfp++;
			}
			*tsp++ = *rfp++;
		}
		*tsp = '\0';
		if( *rfp == '"' )
			rfp++;
		tok = TOK_STRING;
	}else if( isdigit( *rfp ) ){
		tok = TOK_INT;
		*tsp++ = *rfp++; 
		for( ; isdigit( *rfp ); rfp++ )
			*tsp++ = *rfp;
		*tsp = '\0';
		tokval = atoi( tokstr );
	}else if( *rfp == '#' ){
		tok = TOK_COMMENT;
		*tsp++ = *rfp++; 
		for( ; *rfp && *rfp != '\n'; rfp++ )
			*tsp++ = *rfp;
		*tsp = '\0';
	}else if( *rfp == '=' ){
		tok = TOK_EQUAL;
		*tsp++ = *rfp++; 
		*tsp = '\0';
	}else if( *rfp == '{' ){
		tok = TOK_LBRACE;
		*tsp++ = *rfp++; 
		*tsp = '\0';
	}else if( *rfp == '}' ){
		tok = TOK_RBRACE;
		*tsp++ = *rfp++; 
		*tsp = '\0';
	}else if( *rfp == '[' ){
		tok = TOK_LBRACK;
		*tsp++ = *rfp++; 
		*tsp = '\0';
	}else if( *rfp == ']' ){
		tok = TOK_RBRACK;
		*tsp++ = *rfp++; 
		*tsp = '\0';
	}else if( *rfp == ',' ){
		tok = TOK_COMMA;
		*tsp++ = *rfp++; 
		*tsp = '\0';
	}else if( *rfp == '\r' ){
		*tsp++ = *rfp++;
		*tsp = '\0';
		if(*rfp == '\n'){
			tok = TOK_NL;
		}
		else{
			tok = TOK_ERROR;
		}
		*tsp++ = *rfp++;
	}else if( *rfp == '\n' ){
		tok = TOK_NL;
		*tsp++ = *rfp++; 
		*tsp = '\0';
	}else{
		tok = TOK_ERROR;
		*tsp++ = *rfp++; 
		*tsp = '\0';
	}
	return( TRUE );
}

static	void	ungettoken()
{

	ra_tok = tok;
	ra_tokval = tokval;
	strcpy( ra_tokstr, tokstr );
	readahead = TRUE;
}

static	void	skiptonl( FILE *fp )
{

	while( tok != TOK_NL && tok != TOK_EOF )
		gettoken( fp );
}

static	void	clear_rule( RULE_T *rp )
{
	int	i;

	rp->r_sym = UNDEF;
	rp->r_prec = UNDEF;
	rp->r_print = NULL;
	rp->r_lclass = 0;
	rp->r_rclass = 0;
	rp->r_octok = 0;
	rp->r_ocval = 0;
	rp->r_nakind = 0;
	rp->r_akind = a_kind;
	rp->r_nlkind = 0;
	rp->r_lkind = l_kind;
	rp->r_nrkind = 0;
	rp->r_rkind = r_kind;
	rp->r_nokind = 0;
	rp->r_oktok = o_ktok;
	rp->r_okval = o_kval;
	rp->r_okqtok = o_kqtok;
	rp->r_okqval = o_kqval;
	for( i = 0; i < KIND_SIZE; i++ ){
		rp->r_akind[i] = 0;
		rp->r_lkind[i] = 0;
		rp->r_rkind[i] = 0;
		rp->r_oktok[i] = 0;
		rp->r_okval[i] = 0;
		rp->r_okqtok[i] = 0;
		rp->r_okqval[i] = 0;
	}
	rp->r_ntypes = 0;
	rp->r_types = ttabs;
}

static	RULE_T	*save_rule( RULE_T *rp )
{
	RULE_T	*nrp;
	char	*sp;
	int	i, *ip;

	nrp = ( RULE_T * )malloc( sizeof( RULE_T ) );
	if( nrp == NULL ){
		fprintf( stderr, "%s: can't allocate new rule.\n", progname );
		return( NULL );
	}
	nrp->r_sym = rp->r_sym;
	nrp->r_prec = rp->r_prec;
	sp = ( char * )malloc( strlen( rp->r_print ) + 1 );
	if( sp == NULL ){
		fprintf( stderr, "%s: can't allocaate r_print for new rule.\n",
			progname );
		return( NULL );
	}
	strcpy( sp, rp->r_print );
	nrp->r_print = sp;

	nrp->r_lclass = rp->r_lclass;
	nrp->r_rclass = rp->r_rclass;
	nrp->r_octok = rp->r_octok;
	nrp->r_ocval = rp->r_ocval;

	nrp->r_nakind = rp->r_nakind;
	if( nrp->r_nakind > 0 ){
		ip = ( int * )malloc( nrp->r_nakind * sizeof( int ) );
		if( ip == NULL ){
			fprintf( stderr,
				"%s: can't allocate r_akind for new rule.\n", 
				progname );
			return( NULL );
		}
		nrp->r_akind = ip;
		for( i = 0; i < nrp->r_nakind; i++ )
			nrp->r_akind[ i ] = rp->r_akind[ i ];
		
	}else
		nrp->r_akind = NULL;

	nrp->r_nlkind = rp->r_nlkind;
	if( nrp->r_nlkind > 0 ){
		ip = ( int * )malloc( nrp->r_nlkind * sizeof( int ) );
		if( ip == NULL ){
			fprintf( stderr,
				"%s: can't allocate r_lkind for new rule.\n", 
				progname );
			return( NULL );
		}
		nrp->r_lkind = ip;
		for( i = 0; i < nrp->r_nlkind; i++ )
			nrp->r_lkind[ i ] = rp->r_lkind[ i ];
		
	}else
		nrp->r_lkind = NULL;

	nrp->r_nrkind = rp->r_nrkind;
	if( nrp->r_nrkind > 0 ){
		ip = ( int * )malloc( nrp->r_nrkind * sizeof( int ) );
		if( ip == NULL ){
			fprintf( stderr,
				"%s: can't allocate r_rkind for new rule.\n", 
				progname );
			return( NULL );
		}
		nrp->r_rkind = ip;
		for( i = 0; i < nrp->r_nrkind; i++ )
			nrp->r_rkind[ i ] = rp->r_rkind[ i ];
		
	}else
		nrp->r_rkind = NULL;

	nrp->r_nokind = rp->r_nokind;
	if( nrp->r_nokind > 0 ){
		ip = ( int * )malloc( nrp->r_nokind * sizeof( int ) );
		if( ip == NULL ){
			fprintf( stderr,
				"%s: can't allocate r_oktok for new rule.\n", 
				progname );
			return( NULL );
		}
		nrp->r_oktok = ip;
		for( i = 0; i < nrp->r_nokind; i++ )
			nrp->r_oktok[ i ] = rp->r_oktok[ i ];
		
		ip = ( int * )malloc( nrp->r_nokind * sizeof( int ) );
		if( ip == NULL ){
			fprintf( stderr,
				"%s: can't allocate r_okval for new rule.\n", 
				progname );
			return( NULL );
		}
		nrp->r_okval = ip;
		for( i = 0; i < nrp->r_nokind; i++ )
			nrp->r_okval[ i ] = rp->r_okval[ i ];

		ip = ( int * )malloc( nrp->r_nokind * sizeof( int ) );
		if( ip == NULL ){
			fprintf( stderr,
				"%s: can't allocate r_okqtok for new rule.\n", 
				progname );
			return( NULL );
		}
		nrp->r_okqtok = ip;
		for( i = 0; i < nrp->r_nokind; i++ )
			nrp->r_okqtok[ i ] = rp->r_okqtok[ i ];
		
		ip = ( int * )malloc( nrp->r_nokind * sizeof( int ) );
		if( ip == NULL ){
			fprintf( stderr,
				"%s: can't allocate r_okqval for new rule.\n", 
				progname );
			return( NULL );
		}
		nrp->r_okqval = ip;
		for( i = 0; i < nrp->r_nokind; i++ )
			nrp->r_okqval[ i ] = rp->r_okqval[ i ];
	}

	nrp->r_ntypes = rp->r_ntypes;
	nrp->r_types = ( TTAB_T * )malloc( nrp->r_ntypes * sizeof( TTAB_T ) );
	if( nrp->r_types == NULL ){
		fprintf( stderr, "%s: can't allocate r_types for new rule.\n",
			progname );
		return( NULL );
	}
	for( i = 0; i < nrp->r_ntypes; i++ ){
		nrp->r_types[i].t_index = rp->r_types[i].t_index;
		nrp->r_types[i].t_nrows = rp->r_types[i].t_nrows;
		nrp->r_types[i].t_rows = rp->r_types[i].t_rows;
	}

	return( nrp );
}

static	int	check_rule( RULE_T *rp )
{
	int	err;
	int	t, nkind;
	TTAB_T	*ttp;

	err = FALSE;
	if( ( nkind = rp->r_nakind ) > 1 && nkind != rp->r_ntypes ){
		fprintf( stderr, "%d akind values requires %d types.\n",
			nkind, nkind );
		return( TRUE );
	}
	if( ( nkind = rp->r_nlkind ) > 1 && nkind != rp->r_ntypes ){
		fprintf( stderr, "%d left kind values requires %d types.\n",
			nkind, nkind );
		return( TRUE );
	}
	if( ( nkind = rp->r_nrkind ) > 1 && nkind != rp->r_ntypes ){
		fprintf( stderr,
			"%d right kind values requires %d types.\n",
			nkind, nkind );
		return( TRUE );
	}

	for( t = 0; t < rp->r_ntypes; t++ ){
		ttp = &rp->r_types[ t ];
		err |= check_type( ttp );
	}

	return( err );
}

static	int	check_type( TTAB_T *ttp )
{
	int	err;
	int	t, n_rows;
	int	i, j;
	int	t1, ot, nt;
	TROW_T	*trp, *ntrp, *otrp;
	DEF_T	*dp;
	int	ins, del;
	int	lval;
	int	dv2rv[ TVTAB_SIZE ];
	int	rv2dv[ TVTAB_SIZE ];
	TROW_T	**newrows;
	char	*tnp;
	int	*nrp;

	err = FALSE;

	/* Square up the ttab array by adding/removing		*/
	/* elements:						*/

	n_rows = ttp->t_nrows;
	for( t = 0; t < n_rows; t++ ){
		trp = ttp->t_rows[ t ];
		trp->t_nval = t;
		if( trp->t_nrvals < n_rows ){
			fprintf( stderr,
				"%s: short type row: add %d T_ERROR's.\n",
				rfname, n_rows - trp->t_nrvals );
			nrp = ( int * )malloc( n_rows * sizeof( int ) );
			for( i = 0; i < trp->t_nrvals; i++ )
				nrp[i] = trp->t_rvals[i];
			for( i = trp->t_nrvals; i < n_rows; i++ )
				nrp[i] = v_error;
			free( trp->t_rvals );
			trp->t_rvals = nrp;
		}else if( trp->t_nrvals > n_rows ){
			fprintf( stderr,
				"%s: long type row: truncate last %d values.\n",
				rfname, trp->t_nrvals - n_rows );
			nrp = ( int * )malloc( n_rows * sizeof( int ) );
			for( i = 0; i < n_rows; i++ )
				nrp[i] = trp->t_rvals[i];
			free( trp->t_rvals );
			trp->t_rvals = nrp;
		}
	}

	/* construct the mapping from values in nab.h to the	*/
	/* values from the rows of the current type		*/
	/* types of interest run from 0 .. l_utype (the last	*/
	/* type)						*/

	for( i = 0; i <= l_utype; i++ )
		dv2rv[ i ] = UNDEF;
	ins = FALSE;
	for( t = 0; t < n_tvtab; t++ ){
		dp = &tvtab[ t ];
		if((trp = findtrow( dp->d_name, n_rows, ttp->t_rows ))){
			dv2rv[ dp->d_val ] = trp->t_nval;
			dp->d_used = TRUE;
		}else
			ins = TRUE;
	}

	/* ensure that after any insert/delete ops, that the 	*/
	/* symbols defined in the type table have the same 	*/
	/* order as the symbols in TCK_FNAME			*/

	for( lval = UNDEF, i = 0; i <= l_utype; i++ ){
		if( dv2rv[ i ] == UNDEF )
			continue;
		else{
			if( lval == UNDEF )
				lval = dv2rv[ i ];
			else if( dv2rv[ i ]  < lval ){
				fprintf( stderr,
		"%s: type table row order differs from that of %s.\n",
					rfname, nab_h_path );
				return( TRUE );
			}else
				lval = dv2rv[ i ];
		}
	}

	/* construct the mapping from values in the type table	*/
	/* back to the values in nab.h				*/
	for( i = 0; i < n_rows; i++ )
		rv2dv[ i ] = UNDEF;
	del = FALSE;
	for( i = 0; i < n_rows; i++ ){
		trp = ttp->t_rows[ i ];
		if( (dp = finddef( trp->t_name, n_tvtab, tvtab)) )
			rv2dv[ i ] = dp->d_val;
		else
			del = TRUE;
	}

	/* no inserts or deletes, table is ok as far as	*/
	/* this program is concerned.  Adding/removing	*/
	/* or changing the permitted operations and the	*/
	/* result types is beyond the scope of this 	*/
	/* program					*/

	if( !ins && !del )
		return( err );

	newrows = ( TROW_T ** )malloc( (l_utype+1) * sizeof( TROW_T * ) );
	for( i = 0; i <= l_utype; i++ ){
		ntrp = ( TROW_T * )malloc( sizeof( TROW_T ) );
		tnp = ( char * )malloc( strlen( tvtab[i].d_name ) + 1 );
		strcpy( tnp, tvtab[i].d_name );
		ntrp->t_name = tnp;
		ntrp->t_nval = i;
		ntrp->t_nrvals = l_utype + 1;
		nrp = ( int * )malloc( (l_utype+1) * sizeof( int ) );
		ntrp->t_rvals = nrp;
		for( j = 0; j <= l_utype; j++ )
			ntrp->t_rvals[j] = v_error;	
		newrows[ i ] = ntrp;
	}

	/* fix the table:					*/
	for( i = 0; i <= l_utype; i++ ){
		ntrp = newrows[ i ];
		/* UNDEF means new type, add row.  new row was	*/
		/* to all T_UNDEF's above, so its done		*/
		if( ( t1 = dv2rv[ i ] ) == UNDEF ){
			fprintf( stderr, "%s: add row for type = %s.\n",
				rfname, tvtab[ i ].d_name );
			continue;
		}
		otrp = ttp->t_rows[ t1 ];
		for( j = 0; j <= l_utype; j++ ){
			if( ( t1 = dv2rv[ j ] ) == UNDEF )
				continue;
			ot = otrp->t_rvals[ t1 ];
			nt = rv2dv[ ot ];
			if( nt == UNDEF ){
				fprintf( stderr,
					"%s: row uses undefind value %d\n",
					rfname, nt );
				return( TRUE );
			}else
				ntrp->t_rvals[ j ] = nt;
		}
	}

	ttp->t_nrows = l_utype + 1;
	ttp->t_rows = newrows;

	return( err );
}
static	TROW_T	*findtrow( char name[], int n_trows, TROW_T *trows[] )
{
	int	i;
	TROW_T	*trp;

	for( i = 0; i < n_trows; i++ ){
		trp = trows[ i ];
		if( !strcmp( name, trp->t_name ) )
			return( trp );
	}
	return( NULL );
}

static	void	fprint_rule( FILE *fp, RULE_T *rp )
{
	int	i, j, k;
	DEF_T	*dp;

	fprintf( fp, "operator\n" );
	dp = finddef_withval( rp->r_sym, n_symtab, symtab );
	fprintf( fp, "\tsym    = %s\n", dp->d_name );
	fprintf( fp, "\tprec   = %d\n", rp->r_prec );
	fprintf( fp, "\tprint  = \"%s\"\n", rp->r_print );
	fprintf( fp, "\n" );
	fprintf( fp, "class\n" );
	if( rp->r_lclass != 0 ){
		fprintf( fp, "\tleft   = " );
		fprint_set( fp, rp->r_lclass, n_cvtab, cvtab );
		fprintf( fp, "\n" );
	}
	if( rp->r_rclass != 0 ){
		fprintf( fp, "\tright  = " );
		fprint_set( fp, rp->r_rclass, n_cvtab, cvtab );
		fprintf( fp, "\n" );
	}
/*
	fprintf( fp, "\toutput = %s\n", rp->r_oclass );
*/
	fprintf( fp, "\toutput = " );
	if( rp->r_octok == TOK_LEFT )
		fprintf( fp, "left\n" );
	else if( rp->r_octok == TOK_RIGHT )
		fprintf( fp, "right" );
	else{
		dp = finddef_withval( rp->r_ocval, n_cvtab, cvtab );
		fprintf( fp, "%s", dp->d_name );
	}
	fprintf( fp, "\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "kind\n" );
	if( rp->r_nakind == 1 ){
		fprintf( fp, "\takind  = " );
		fprint_set( fp, *rp->r_akind, n_kvtab, kvtab );
		fprintf( fp, "\n" );
	}else{
		for( i = 0; i < rp->r_nakind; i++ ){
			fprintf( fp, "\takind  = " );
			fprint_set( fp, rp->r_akind[i], n_kvtab, kvtab );
			fprintf( fp, " use type = %d\n", i + 1 );
		}
	}
	if( rp->r_nlkind == 1 ){
		fprintf( fp, "\tleft   = " );
		fprint_set( fp, *rp->r_lkind, n_kvtab, kvtab );
		fprintf( fp, "\n" );
	}else{
		for( i = 0; i < rp->r_nlkind; i++ ){
			fprintf( fp, "\tleft   = " );
			fprint_set( fp, rp->r_lkind[i], n_kvtab, kvtab );
			fprintf( fp, " use type = %d\n", i + 1 );
		}
	}
	if( rp->r_nrkind == 1 ){
		fprintf( fp, "\tright  = " );
		fprint_set( fp, *rp->r_rkind, n_kvtab, kvtab );
		fprintf( fp, "\n" );
	}else{
		for( i = 0; i < rp->r_nrkind; i++ ){
			fprintf( fp, "\tright  = " );
			fprint_set( fp, rp->r_rkind[i], n_kvtab, kvtab );
			fprintf( fp, " use type = %d\n", i + 1 );
		}
	}
/*
	fprintf( fp, "\toutput = %s\n", rp->r_okind );
	fprintf( fp, "\n" );
*/
	for( i = 0; i < rp->r_nokind; i++ ){
		fprintf( fp, "\toutput = " );
		if( rp->r_oktok[i] == TOK_IDENT ){
			dp = finddef_withval( rp->r_okval[i], n_kvtab, kvtab );
			fprintf( fp, "%s", dp->d_name );
		}else if( rp->r_oktok[i] == TOK_LEFT )
			fprintf( fp, "left" );
		else if( rp->r_oktok[i] == TOK_RIGHT )
			fprintf( fp, "right" );
		if( rp->r_okqtok[i] == TOK_LEFT ){
			fprintf( fp, " if left = " );
			fprint_set( fp, rp->r_okqval[i], n_kvtab, kvtab );
		}else if( rp->r_okqtok[i] == TOK_LEFT ){
			fprintf( fp, " if right = " );
			fprint_set( fp, rp->r_okqval[i], n_kvtab, kvtab );
		}
		fprintf( fp, "\n" );
	}
	fprintf( fp, "\n" );

	for( i = 0; i < rp->r_ntypes; i++ ){
		fprintf( fp, "type" );
		if( rp->r_ntypes > 1 )
			fprintf( fp, " = %d", i + 1 );
		fprintf( fp, "\n" );
		for( j = 0; j < rp->r_types[i].t_nrows; j++ ){
			fprintf( fp, "\t%-12s = [",
				rp->r_types[i].t_rows[j]->t_name );
			for(k = 0; k < rp->r_types[i].t_rows[j]->t_nrvals; k++){
				fprintf( fp, " %2d", 
					rp->r_types[i].t_rows[j]->t_rvals[k] );
				if( k < rp->r_types[i].t_rows[j]->t_nrvals-1 )
					fprintf( fp, "," );
			}
			fprintf( fp, " ]\n" );
		}
		fprintf( fp, "\n" );
	}

}

static	void	fprint_set( FILE *fp, int set, int n_tab, DEF_T tab[] )
{
	int	i, c, card;

	for( card = 0, i = 0; i < n_tab; i++ ){
		if( set & ( 1 << i ) )
			card++;
	}
	if( card == 1 ){
		for( i = 0; i < n_tab; i++ ){
			if( set & ( 1 << i ) ){
				fprintf( fp, "%s", tab[i].d_name );
				break;
			}
		}
	}else{
		fprintf( fp, "{" );
		for( c = 0, i = 0; i < n_tab; i++ ){
			if( set & ( 1 << i ) ){
				fprintf( fp, " %s", tab[i].d_name );
				c++;
				if( c < card )
					fprintf( fp, "," );
			}
		}
		fprintf( fp, " }" );
	}
}

static	void	fprint_ovalue( FILE *fp, char *oval, int n_tab, DEF_T tab[] )
{

	fputs( oval, fp );
		
}

static	int	mk_frameset( char* html_output_folder )
{
	FILE	*fp;
	int	err;
	char	hstr[ 256 ];
	
	err = FALSE;

	// figure out file path
	// trying to avoid fixed-length buffers, at least here -JS
	char* frameset_html_file_path = malloc(strlen(html_output_folder) + strlen(FRFILE) + 2); // +1 for the slash, and +1 for the null terminator

	sprintf(frameset_html_file_path, "%s/%s", html_output_folder, FRFILE);

	if( ( fp = fopen( frameset_html_file_path, "w" ) ) == NULL ){
		fprintf( stderr, "%s: can't write frame file %s.\n",
			progname, frameset_html_file_path );
		return( TRUE );
	}

	free(frameset_html_file_path);

	fprintf( fp, "<HTML>\n" );
	str2hstr( "NAB Operator Rules", hstr );
	fprintf( fp, "<HEAD><TITLE>%s</TITLE></HEAD>\n", hstr );
	fprintf( fp, "</HTML>\n" );
	fprintf( fp, "<FRAMESET rows=%d%%,*>\n", CFSIZE );
	fprintf( fp, "<FRAME src=\"%s\" name=\"ControlFrame\">\n", CTLFILE );
	fprintf( fp, "<FRAME src=\"%s\" name=\"RuleFrame\">", r1fname );
	fprintf( fp, "</FRAMESET>\n" );
	fprintf( fp, "</HTML>\n" );
	fclose( fp );
	return( err );
}

static	int	mk_controlfile( char* html_output_folder )
{
	int	err;
	FILE	*fp;
	char	hstr[ 256 ];
	int	i, s;
	DEF_T	**s_symtab;

	err = FALSE;
	s_symtab = ( DEF_T ** )malloc( n_symtab * sizeof( DEF_T * ) );
	if( s_symtab == NULL ){
		fprintf( stderr, "%s: can't allocate s_symtab.\n", progname );
		return( TRUE );
	}

	// figure out file path
	// trying to avoid fixed-length buffers, at least here -JS
	char* control_html_file_path = malloc(strlen(html_output_folder) + strlen(CTLFILE) + 2); // +1 for the slash, and +1 for the null terminator

	sprintf(control_html_file_path, "%s/%s", html_output_folder, CTLFILE);

	if( ( fp = fopen( control_html_file_path, "w" ) ) == NULL){
		fprintf( stderr, "%s: can't write control file %s.\n",
			progname, control_html_file_path );
		return( TRUE );
	}

	free(control_html_file_path);

	str2hstr( "ControlFile", hstr );
	fprintf( fp, "<HTML>\n" );
	fprintf( fp, "<HEAD>\n" );
	fprintf( fp, "<TITLE>%s</TITLE>\n", hstr );
	fprintf( fp, "<SCRIPT language=\"JavaScript\" src=\"update.js\">\n" );
	fprintf( fp, "</SCRIPT>\n" );
	fprintf( fp, "</HEAD>\n" );

	fprintf( fp, "<BODY>\n" );
	fprintf( fp, "<TABLE>\n" );

	fprintf( fp, "<TR>\n" );

	fprintf( fp, "<TD>\n" );
	fprintf( fp, "<FORM name=\"f_dop\">\n" );
	fprintf( fp, "<TABLE border>\n" );
	fprintf( fp, "<CAPTION>Defined Operators</CAPTION>\n" );
	fprintf( fp, "<TR>\n" );
	fprintf( fp, "<TD colspan=2>\n" );
	for( s = 0, i = 0; i < n_symtab; i++ ){ 
		if( symtab[ i ].d_used ){
			s_symtab[ s ] = &symtab[ i ];
			s++;
		}
	}
	qsort( s_symtab, s, sizeof( DEF_T * ), defcmp );
	fprintf( fp, "Current<INPUT type=\"text\" name=\"d_current\"" );
	fprintf( fp, " value=\"%s\"", s_symtab[0]->d_name );
	fprintf( fp, " onChange=\"update('d','f');\">\n" );
	fprintf( fp, "</TD>\n" );
	fprintf( fp, "</TR>\n" );
	fprintf( fp, "<TR>\n" );
	fprintf( fp, "<TD colspan=2>\n" );
	fprintf( fp, "<SELECT name=\"s_dop\" size=%d", SELSIZE );
	fprintf( fp, " onChange=\"update('d','s');\">\n" );
	for( i = 0; i < s; i++ )
		fprintf( fp, "<OPTION value=\"%s\">%s\n", 
			s_symtab[i]->d_name,
			s_symtab[i]->d_name );
	fprintf( fp, "</SELECT>\n" );
	fprintf( fp, "</TD>\n" );
	fprintf( fp, "</TR>\n" );
	fprintf( fp, "<TR>\n" );
	fprintf( fp, "<TD>\n" );
	fprintf( fp, "<iNPUT type=\"button\" value=\" &lt;- Prev\"" );
	fprintf( fp, " onClick=\"update('d','p');\">\n" );
	fprintf( fp, "</TD>\n" );
	fprintf( fp, "<TD>\n" );
	fprintf( fp, "<INPUT type=\"button\" value=\"Next -&gt;\"" );
	fprintf( fp, " onClick=\"update('d','n');\">\n" );
	fprintf( fp, "</TD>\n" );
	fprintf( fp, "</TR>\n" );
	fprintf( fp, "</TABLE>\n" );
	fprintf( fp, "</FORM>\n" );
	fprintf( fp, "</TD>\n" );

	fprintf( fp, "<TD>\n" );
	fprintf( fp, "<FORM name=\"f_uop\">\n" );
	fprintf( fp, "<TABLE border>\n" );
	fprintf( fp, "<CAPTION>Undefined Operators</CAPTION>\n" );
	fprintf( fp, "<TR>\n" );
	fprintf( fp, "<TD colspan=2>\n" );
	for( s = 0, i = 0; i < n_symtab; i++ ){ 
		if( !symtab[ i ].d_used ){
			s_symtab[ s ] = &symtab[ i ];
			s++;
		}
	}
	qsort( s_symtab, s, sizeof( DEF_T * ), defcmp );
	fprintf( fp, "Current<INPUT type=\"text\" name=\"u_current\"" );
	fprintf( fp, " value=\"%s\"", s_symtab[0]->d_name );
	fprintf( fp, " onChange=\"update('u','f');\">\n" );
	fprintf( fp, "</TD>\n" );
	fprintf( fp, "</TR>\n" );
	fprintf( fp, "<TR>\n" );
	fprintf( fp, "<TD colspan=2>\n" );
	fprintf( fp, "<SELECT name=\"s_uop\" size=%d", SELSIZE );
	fprintf( fp, " onChange=\"update('u','s');\">\n" );
	for( i = 0; i < s; i++ )
		fprintf( fp, "<OPTION value=\"%s\">%s\n", 
			s_symtab[i]->d_name,
			s_symtab[i]->d_name );
	fprintf( fp, "</SELECT>\n" );
	fprintf( fp, "</TD>\n" );
	fprintf( fp, "</TR>\n" );
	fprintf( fp, "<TR>\n" );
	fprintf( fp, "<TD>\n" );
	fprintf( fp, "<INPUT type=\"button\" value=\" &lt;- Prev\"" );
	fprintf( fp, " onClick=\"update('u','p');\">\n" );
	fprintf( fp, "</TD>\n" );
	fprintf( fp, "<TD>\n" );
	fprintf( fp, "<INPUT type=\"button\" value=\"Next -&gt;\"" );
	fprintf( fp, " onClick=\"update('u','n');\">\n" );
	fprintf( fp, "</TD>\n" );
	fprintf( fp, "</TR>\n" );
	fprintf( fp, "</TABLE>\n" );
	fprintf( fp, "</FORM>\n" );
	fprintf( fp, "</TD>\n" );

	fprintf( fp, "</TR>\n" );

	fprintf( fp, "</TABLE>\n" );
	fprintf( fp, "</BODY>\n" );
	fprintf( fp, "</HTML>\n" );

	fclose( fp );
	return( err );
}

static	int	hprint_rule(char* rulePath, RULE_T * rp, char* outputFolder)
{
	int	err;
	char	htmlFilename[ 256 ];
	char	hstr[ 256 ];
	FILE	*hfp;
	DEF_T	*dp, *dp1;
	int	t, i, j, tv;

	err = FALSE;

	/* Figure out the filename of the HTML file */

#ifdef WIN32

#define RULE_FILENAME_BUFSZ 64
	 char ruleFilename[RULE_FILENAME_BUFSZ];

	_splitpath_s(rulePath, NULL, 0, NULL, 0, ruleFilename, RULE_FILENAME_BUFSZ, NULL, 0);
#else
	char* ruleFilename = basename(rulePath);

	 //attempt to reterminate ruleFilename at the .
	char* rnp = ruleFilename; //JS: if this line is put inside the for loop on the next line, the Intel compiler derps out and can't compile it
	for(; *rnp; rnp++ ){
		if( *rnp == '.' )
		{
			*rnp = '\0';
		}
	}
#endif



	sprintf(htmlFilename, "%s/%s.html", outputFolder, ruleFilename);

	if( ( hfp = fopen( htmlFilename, "w" ) ) == NULL ){
		fprintf( stderr, "%s: can't write html file %s.\n",
			progname, htmlFilename );
		return( TRUE );
	}
	if( *r1fname == '\0' )
		strcpy( r1fname, htmlFilename );
	str2hstr( rulePath, hstr );
	fprintf( hfp, "<HTML><HEAD><TITLE>%s</TITLE></HEAD>\n", hstr ); 
	fprintf( hfp, "<BODY>\n" );

	fprintf( hfp, "<B>operator</B><BR>\n" );
	fprintf( hfp,
		"<TABLE border bgcolor=white cellpadding=2 cellspacing=2>\n" );
	fprintf( hfp, "<TR>\n" );
	fprintf( hfp, "<TD width=%d>sym</TD>\n", TWID1 );
	dp = finddef_withval( rp->r_sym, n_symtab, symtab );
	str2hstr( dp->d_name, hstr );
	fprintf( hfp, "<TD width=%d>%s</TD>\n", TWID2, hstr );
	fprintf( hfp, "</TR>\n" );
	fprintf( hfp, "<TR>\n" );
	fprintf( hfp, "<TD>prec</TD>\n" );
	fprintf( hfp, "<TD>%d</TD>\n", rp->r_prec );
	fprintf( hfp, "</TR>\n" );
	fprintf( hfp, "<TR>\n" );
	fprintf( hfp, "<TD>print</TD>\n" );
	str2hstr( rp->r_print, hstr );
	fprintf( hfp, "<TD>%s</TD>\n", hstr );
	fprintf( hfp, "</TR>\n" );
	fprintf( hfp, "</TABLE>\n" );
	fprintf( hfp, "<BR>\n" );

	fprintf( hfp, "<B>class</B><BR>\n" );
	fprintf( hfp,
		"<TABLE border bgcolor=white cellpadding=2 cellspacing=2>\n" );
	if( rp->r_lclass != 0 ){
		fprintf( hfp, "<TR>\n" );
		fprintf( hfp, "<TD width=%d>left</TD>\n", TWID1 );
		fprintf( hfp, "<TD width=%d>\n", TWID2 );
		fprint_set( hfp, rp->r_lclass, n_cvtab, cvtab );
		fprintf( hfp, "</TD>\n" );
		fprintf( hfp, "</TR>\n" );
	}
	if( rp->r_rclass != 0 ){
		fprintf( hfp, "<TR>\n" );
		fprintf( hfp, "<TD width=%d>right</TD>\n", TWID1 );
		fprintf( hfp, "<TD width=%d>\n", TWID2 );
		fprint_set( hfp, rp->r_rclass, n_cvtab, cvtab );
		fprintf( hfp, "</TD>\n" );
		fprintf( hfp, "</TR>\n" );
	}
	fprintf( hfp, "<TR>\n" );
	fprintf( hfp, "<TD width=%d>output</TD>\n", TWID1 );
	fprintf( hfp, "<TD width=%d>\n", TWID2 );
/*
	fprint_ovalue( hfp, rp->r_oclass, n_cvtab, cvtab );
*/
	if( rp->r_octok == TOK_LEFT )
		fprintf( hfp, "left" );
	else if( rp->r_octok == TOK_RIGHT )
		fprintf( hfp, "right" );
	else{
		dp = finddef_withval( rp->r_ocval, n_cvtab, cvtab );
		fprintf( hfp, "%s", dp->d_name );
	}
	fprintf( hfp, "</TD>\n" );
	fprintf( hfp, "</TR>\n" );
	fprintf( hfp, "</TABLE>\n" );
	fprintf( hfp, "<BR>\n" );

	fprintf( hfp, "<B>kind</B><BR>\n" );
	fprintf( hfp,
		"<TABLE border bgcolor=white cellpadding=2 cellspacing=2>\n" );
	if( rp->r_nakind == 1 ){
		fprintf( hfp, "<TR>\n" );
		fprintf( hfp, "<TD width=%d>akind</TD>\n", TWID1 );
		fprintf( hfp, "<TD width=%d>\n", TWID2 );
		fprint_set( hfp, *rp->r_akind, n_kvtab, kvtab );
		fprintf( hfp, "</TD>\n" );
		fprintf( hfp, "</TR>\n" );
	}else{
		for( i = 0; i < rp->r_nakind; i++ ){
			fprintf( hfp, "<TR>\n" );
			fprintf( hfp, "<TD width=%d>akind</TD>\n", TWID1 );
			fprintf( hfp, "<TD width=%d>\n", TWID2 );
			fprint_set( hfp, rp->r_akind[i], n_kvtab, kvtab );
			fprintf( hfp, " use type = %d", i + 1 );
			fprintf( hfp, "</TD>\n" );
			fprintf( hfp, "</TR>\n" );
		}
	}
	if( rp->r_nlkind == 1 ){
		fprintf( hfp, "<TR>\n" );
		fprintf( hfp, "<TD width=%d>left</TD>\n", TWID1 );
		fprintf( hfp, "<TD width=%d>\n", TWID2 );
		fprint_set( hfp, *rp->r_lkind, n_kvtab, kvtab );
		fprintf( hfp, "</TD>\n" );
		fprintf( hfp, "</TR>\n" );
	}else{
		for( i = 0; i < rp->r_nlkind; i++ ){
			fprintf( hfp, "<TR>\n" );
			fprintf( hfp, "<TD width=%d>left</TD>\n", TWID1 );
			fprintf( hfp, "<TD width=%d>\n", TWID2 );
			fprint_set( hfp, rp->r_lkind[i], n_kvtab, kvtab );
			fprintf( hfp, " use type = %d", i + 1 );
			fprintf( hfp, "</TD>\n" );
			fprintf( hfp, "</TR>\n" );
		}
	}
	if( rp->r_nrkind == 1 ){
		fprintf( hfp, "<TR>\n" );
		fprintf( hfp, "<TD width=%d>right</TD>\n", TWID1 );
		fprintf( hfp, "<TD width=%d>\n", TWID2 );
		fprint_set( hfp, *rp->r_rkind, n_kvtab, kvtab );
		fprintf( hfp, "</TD>\n" );
		fprintf( hfp, "</TR>\n" );
	}else{
		for( i = 0; i < rp->r_nrkind; i++ ){
			fprintf( hfp, "<TR>\n" );
			fprintf( hfp, "<TD width=%d>right</TD>\n", TWID1 );
			fprintf( hfp, "<TD width=%d>\n", TWID2 );
			fprint_set( hfp, rp->r_rkind[i], n_kvtab, kvtab );
			fprintf( hfp, " use type = %d", i + 1 );
			fprintf( hfp, "</TD>\n" );
			fprintf( hfp, "</TR>\n" );
		}
	}
	fprintf( hfp, "<TR>\n" );
/*
	fprintf( hfp, "<TD width=%d>output</TD>\n", TWID1 );
	fprintf( hfp, "<TD width=%d>\n", TWID2 );
	fprint_ovalue( hfp, rp->r_okind, n_kvtab, kvtab );
	fprintf( hfp, "</TD>\n" );
	fprintf( hfp, "</TR>\n" );
*/
	for( i = 0; i < rp->r_nokind; i++){
		fprintf( hfp, "<TR>\n" );
		fprintf( hfp, "<TD width=%d>output</TD>\n", TWID1 );
		fprintf( hfp, "<TD width=%d>\n", TWID2 );
		if( rp->r_oktok[i] == TOK_LEFT )
			fprintf( hfp, "left" );
		else if( rp->r_oktok[i] == TOK_RIGHT )
			fprintf( hfp, "right" );
		else{
			dp = finddef_withval( rp->r_okval[i], n_kvtab,kvtab );
			fprintf( hfp, "%s", dp->d_name );
			if( rp->r_okqtok[i] == TOK_LEFT ){
				fprintf( hfp, " if left = " );
				fprint_set( hfp, rp->r_okqval[i],
					n_kvtab, kvtab );
			}else if( rp->r_okqtok[i] == TOK_RIGHT ){
				fprintf( hfp, " if left = " );
				fprint_set( hfp, rp->r_okqval[i],
					n_kvtab, kvtab );
			}
		}
		fprintf( hfp, "</TD>\n" );
		fprintf( hfp, "</TR>\n" );
	}
	fprintf( hfp, "</TD>\n" );
	fprintf( hfp, "</TR>\n" );

	fprintf( hfp, "</TABLE>\n" );
	fprintf( hfp, "<BR>\n" );

	for( t = 0; t < rp->r_ntypes; t++ ){
		fprintf( hfp, "<B>type" );
		if( rp->r_ntypes > 1 )
			fprintf( hfp, " = %d", t + 1 );
		fprintf( hfp, ", t_index = %d</B><BR>\n",
			rp->r_types[t].t_index );
		fprintf( hfp,
		"<TABLE border bgcolor=white cellpadding=2 cellspacing=2>\n" );
		fprintf( hfp, "<TR>\n" );
		fprintf( hfp, "<TD width=%d>left/right</TD>\n", TWID1 );
		fprintf( hfp, "<TD></TD>\n" );
		for( i = 0; i <= l_utype; i++ ){
			dp = finddef_withval( i, n_tvtab, tvtab );
			fprintf( hfp, "<TD>%s</TD>\n", dp->d_name );
		}
		fprintf( hfp, "</TR>\n" );
		fprintf( hfp, "<TR>\n" );
		for( i = 0; i <= l_utype + 1; i++ )
			fprintf( hfp, "<TD></TD>\n" );
		fprintf( hfp, "</TR>\n" );
		for( i = 0; i <= l_utype; i++ ){
			fprintf( hfp, "<TR>\n" );
			dp = finddef_withval( i, n_tvtab, tvtab );
			fprintf( hfp, "<TD width=%d>%s</TD>\n",
				TWID1, dp->d_name );
			fprintf( hfp, "<TD></TD>\n" );
			for( j = 0; j <= l_utype; j++ ){
				tv = rp->r_types[t].t_rows[i]->t_rvals[j];
				if( tv == v_error )
					fprintf( hfp, "<TD><BR></TD>\n" );
				else{
					dp1 = finddef_withval(tv,n_tvtab,tvtab);
					fprintf( hfp, "<TD>%s</TD>\n",
					dp1->d_name );
				}
			}
			fprintf( hfp, "</TR>\n" );
		}
		fprintf( hfp, "</TABLE>\n" );
		fprintf( hfp, "<BR>\n" );
	}

	fprintf( hfp, "</BODY>\n" );
	fprintf( hfp, "</HTML>\n" );
	fclose( hfp );

	return( err );
}

static	char	*str2hstr( char str[], char hstr[] )
{
	char	*sp, *hp;

	for( sp = str, hp = hstr; *sp; sp++ ){
		if( *sp == '<' ){
			strcpy( hp, "&lt;" );
			hp += 4;
		}else if( *sp == '>' ){
			strcpy( hp, "&gt;" );
			hp += 4;
		}else if( *sp == '&' ){
			strcpy( hp, "&amp;" );
			hp += 5;
		}else
			*hp++ = *sp;
	}
	*hp = '\0';
	return( hstr );
}

static	int	mk_checkexpr( char* output_file_path, int n_rules )
{
	int	err;
	FILE	*cfp;
	int	r;

	err = FALSE;
	if( ( cfp = fopen( output_file_path, "w" ) ) == NULL ){
		fprintf( stderr, "%s: can't write code file %s.\n", progname, output_file_path );
		return( TRUE );
	}
	fprintf( cfp, "#include <stdio.h>\n" );
	fprintf( cfp, "#include <string.h>\n" );
	fprintf( cfp, "\n" );
	fprintf( cfp, "#include \"nab.h\"\n" );
	fprintf( cfp, "#include \"nabgrm.tab.h\"\n" );
	fprintf( cfp, "#include \"errormsg.h\"\n" );
	fprintf( cfp, "#include \"symbol.h\"\n" );
	fprintf( cfp, "\n" );
	fprintf( cfp, "#define BIT(e) (1<<(e))\n" );
	fprintf( cfp, "#define INSET(e,s) (BIT(e)&(s))\n" );
	fprintf( cfp, "\n" );
	fprintf( cfp, "extern int cg_emsg_lineno;\n" );
	fprintf( cfp, "extern SYMREC_T *astk[];\n" );
	fprintf( cfp, "extern int astkp;\n" );
	fprintf( cfp, "static int akind;\n" );	
	fprintf( cfp, "SYMREC_T *findsym();\n" );
	fprintf( cfp, "static void checkattr();\n" );
	fprintf( cfp, "void checkid();\n" );
	fprintf( cfp, "\n" );

	fprintf( cfp, "static char typetab[][%d][%d] = {\n",
		l_utype+1, l_utype+1 );
	for( r = 0; r < n_rules; r++ )
		mk_typetabent( cfp, r == 0, rules[ r ] );
	fprintf( cfp, "};\n" );

	fprintf( cfp, "\n" );
	fprintf( cfp, "void checkexpr( expr )\n" );
	fprintf( cfp, "NODE_T *expr;\n" );
	fprintf( cfp, "{\n" );
	fprintf( cfp, "  NODE_T *npl, *npr, *npd;\n" );
	fprintf( cfp, "  char *fname;\n" );
	fprintf( cfp, "  int class, nd;\n" );
	fprintf( cfp, "  SYMREC_T *s_array;\n" );
	fprintf( cfp, "  int l_class, r_class, o_class;\n" );
	fprintf( cfp, "  int l_kind, r_kind, o_kind;\n" );
	fprintf( cfp, "  int t_idx, l_type, r_type, o_type;\n" );
	fprintf( cfp, "\n" );
	fprintf( cfp, "  if( expr ){\n" );
	fprintf( cfp, "    npl = expr->n_left;\n" );
	fprintf( cfp, "    npr = expr->n_right;\n" );
	fprintf( cfp, "    if( expr->n_sym==SYM_LBRACK ){\n" );
	fprintf( cfp,
		"      s_array = findsym( npl->n_val.v_value.v_cval );\n" );
	fprintf( cfp, "      akind = s_array ? s_array->s_kind : K_UNDEF;\n" );
	fprintf( cfp, "      astk[ astkp ] = s_array;\n" );
	fprintf( cfp, "      astkp++;\n" );
	fprintf( cfp, "    }\n" );
	fprintf( cfp, "    checkexpr( npl );\n" );
	fprintf( cfp, "    checkexpr( npr );\n" );
	fprintf( cfp, "\n" );
	fprintf( cfp, "    cg_emsg_lineno = expr->n_lineno;\n" );
	fprintf( cfp, "\n" );
	fprintf( cfp, "    l_class = npl ? npl->n_class : C_UNDEF;\n" ); 
	fprintf( cfp, "    r_class = npr ? npr->n_class : C_UNDEF;\n" ); 
	fprintf( cfp, "    l_kind = npl ? npl->n_kind : K_UNDEF;\n" ); 
	fprintf( cfp, "    r_kind = npr ? npr->n_kind : K_UNDEF;\n" ); 
	fprintf( cfp, "    l_type = npl ? npl->n_type : T_UNDEF;\n" ); 
	fprintf( cfp, "    r_type = npr ? npr->n_type : T_UNDEF;\n" ); 
	fprintf( cfp, "\n" );
	fprintf( cfp, "    switch( expr->n_sym ){\n" );
	fprintf( cfp, "\n" );

	/* NOTE: code for SYM_IDENT is not from a rule-file!	*/	
	fprintf( cfp, "    case SYM_IDENT :\n" );
	fprintf( cfp, "      checkid( expr );\n" );
	fprintf( cfp, "      break;\n" );
	fprintf( cfp, "\n" );

	/* NOTE: code for SYM_PERIOD is not from a rule-file!	*/	
	fprintf( cfp, "    case SYM_PERIOD :\n" );
	fprintf( cfp, "      checkattr( expr, npl, npr );\n" );
	fprintf( cfp, "      break;\n" );
	fprintf( cfp, "\n" );

	for( r = 0; r < n_rules; r++ )
	{
		mk_rulecase( cfp, rules[ r ] );
	}
	
	fprintf( cfp, "\n" );
	fprintf( cfp, "    }\n" );
	fprintf( cfp, "\n" );
	fprintf( cfp, "    if( expr->n_sym==SYM_LBRACK ){\n" );
	fprintf( cfp, "      astkp--;\n" );
	fprintf( cfp, "      if( astkp > 0 )\n" );
	fprintf( cfp, "        akind = astk[astkp-1]->s_kind;\n" );
	fprintf( cfp, "    }\n" );
	fprintf( cfp, "\n" );
	fprintf( cfp, "  }\n" );
	fprintf( cfp, "}\n" );
	mk_checkattr( cfp );
	fclose( cfp );
	return( err );
}

void
mk_typetabent( FILE *fp, int first, RULE_T *rp )
{
	int	t, i, j;
	TTAB_T	*ttp;
	TROW_T	*trp;
	DEF_T	*dp;

	for( t = 0; t < rp->r_ntypes; t++ ){
		if( first ){
			fprintf( fp, "  {\n" );
			first = FALSE;
		}else
			fprintf( fp, "  ,{\n" );
		ttp = &rp->r_types[ t ];
		for( i = 0; i < ttp->t_nrows; i++ ){
			trp = ttp->t_rows[i];
			if( i == 0 )
				fprintf( fp, "    { " );
			else
				fprintf( fp, "   ,{ " );
			for( j = 0; j < trp->t_nrvals; j++ ){
				dp = finddef_withval( trp->t_rvals[j],
					n_tvtab, tvtab );
				fprintf( fp, "%s", dp->d_name );
				if( j < trp->t_nrvals-1 )
					fprintf( fp, "," );
			}
			fprintf( fp, " }\n" );
		}
		fprintf( fp, "  }\n" );
	}
}

void mk_rulecase( FILE *fp, RULE_T *rp )
{
	DEF_T	*dp;
	int	i;

	dp = finddef_withval( rp->r_sym, n_symtab, symtab );
	fprintf( fp, "    case %s :\n", dp->d_name );
	fprintf( fp, "      expr->n_class = C_UNDEF;\n" );
	fprintf( fp, "      expr->n_kind = K_UNDEF;\n" );
	fprintf( fp, "      expr->n_type = T_UNDEF;\n" );
	fprintf( fp, "      if( l_type == T_ERROR || r_type == T_ERROR ){\n" );
	fprintf( fp, "        expr->n_type = T_ERROR;\n" );
	fprintf( fp, "        break;\n" );
	fprintf( fp, "      }\n" );
	if( rp->r_lclass != 0 ){
		fprintf( fp, "      if( !INSET( l_class, 0%o ) ){\n",
			rp->r_lclass );
		fprintf( fp, "        expr->n_type = T_ERROR;\n" );
		fprintf( fp, "        errormsg( FALSE," );
		fprintf( fp,
			"\"operator '%s': left operand has wrong class.\\n\"",
			rp->r_print );
		fprintf( fp, ");\n" );
		fprintf( fp, "        break;\n" );
		fprintf( fp, "      }\n" );
	}
	if( rp->r_rclass != 0 ){
		fprintf( fp, "      if( !INSET( r_class, 0%o ) ){\n",
			rp->r_rclass );
		fprintf( fp, "        expr->n_type = T_ERROR;\n" );
		fprintf( fp, "        errormsg( FALSE," );
		fprintf( fp,
			"\"operator '%s': right operand has wrong class.\\n\"",
			rp->r_print );
		fprintf( fp, ");\n" );
		fprintf( fp, "        break;\n" );
		fprintf( fp, "      }\n" );
	}
	fprintf( fp, "      o_class = " );
/*
	if( !strcmp( rp->r_oclass, "left" ) )
		fprintf( fp, "l_class;\n" );
	else if( !strcmp( rp->r_oclass, "right" ) )
		fprintf( fp, "r_class;\n" );
	else 
		fprintf( fp, "%s;\n", rp->r_oclass );
*/
	if( rp->r_octok == TOK_LEFT )
		fprintf( fp, "l_class;\n" );
	else if( rp->r_octok == TOK_RIGHT )
		fprintf( fp, "r_class;\n" );
	else{
		dp = finddef_withval( rp->r_ocval, n_cvtab, cvtab );
		fprintf( fp, "%s;\n", dp->d_name );
	}
	fprintf( fp, "      t_idx = %d;\n", rp->r_types[0].t_index );
	if( rp->r_nakind > 0 ){
		for( i = 0; i < rp->r_nakind; i++ ){
			if( i == 0 )
				fprintf( fp, "      " );
			else
				fprintf( fp, "      else " );
			fprintf( fp, "if( INSET( akind, 0%o ) )\n",
				rp->r_akind[i] );
			fprintf( fp, "        t_idx += %d;\n", i );
		}
		fprintf( fp, "      else{\n" );
		fprintf( fp, "        expr->n_type = T_ERROR;\n" );
		fprintf( fp, "        errormsg( FALSE," );
		fprintf( fp,
			"\"operator '%s': array operand has wrong kind.\\n\"",
			rp->r_print );
		fprintf( fp, ");\n" );
		fprintf( fp, "        break;\n" );
		fprintf( fp, "      }\n" );
	}
	if( rp->r_nlkind > 0 ){
		for( i = 0; i < rp->r_nlkind; i++ ){
			if( i == 0 )
				fprintf( fp, "      " );
			else
				fprintf( fp, "      else " );
			fprintf( fp, "if( INSET( l_kind, 0%o ) )\n",
				rp->r_lkind[i] );
			fprintf( fp, "        t_idx += %d;\n", i );
		}
		fprintf( fp, "      else{\n" );
		fprintf( fp, "        expr->n_type = T_ERROR;\n" );
		fprintf( fp, "        errormsg( FALSE," );
		fprintf( fp,
			"\"operator '%s': left operand has wrong kind.\\n\"",
			rp->r_print );
		fprintf( fp, ");\n" );
		fprintf( fp, "        break;\n" );
		fprintf( fp, "      }\n" );
	}
	if( rp->r_nrkind > 0 ){
		for( i = 0; i < rp->r_nrkind; i++ ){
			if( i == 0 )
				fprintf( fp, "      " );
			else
				fprintf( fp, "      else " );
			fprintf( fp, "if( INSET( r_kind, 0%o ) )\n",
				rp->r_rkind[i] );
			fprintf( fp, "        t_idx += %d;\n", i );
		}
		fprintf( fp, "      else{\n" );
		fprintf( fp, "        expr->n_type = T_ERROR;\n" );
		fprintf( fp, "        errormsg( FALSE," );
		fprintf( fp,
			"\"operator '%s': right operand has wrong kind.\\n\"",
			rp->r_print );
		fprintf( fp, ");\n" );
		fprintf( fp, "        break;\n" );
		fprintf( fp, "      }\n" );
	}
/*
	fprintf( fp, "      o_kind = " );
	if( !strcmp( rp->r_okind, "left" ) )
		fprintf( fp, "l_kind;\n" );
	else if( !strcmp( rp->r_okind, "right" ) )
		fprintf( fp, "r_kind;\n" );
	else
		fprintf( fp, "%s;\n", rp->r_okind );
*/

	for( i = 0; i < rp->r_nokind; i++ ){
		if( rp->r_oktok[i] == TOK_LEFT )
			fprintf( fp, "      o_kind = l_kind;\n" );
		else if( rp->r_oktok[i] == TOK_RIGHT )
			fprintf( fp, "      o_kind = r_kind;\n" );
		else if( rp->r_okqtok[i] == TOK_LEFT ){
			fprintf( fp, "      if( INSET( l_kind, 0%o ) )\n", rp->r_okqval[i] );
			dp = finddef_withval( rp->r_okval[i], n_kvtab, kvtab );
			fprintf( fp, "        o_kind = %s;\n", dp->d_name );
		}else if( rp->r_okqtok[i] == TOK_RIGHT ){
			fprintf( fp, "      if( INSET( r_kind, 0%o ) )\n", rp->r_okqval[i] );
			dp = finddef_withval( rp->r_okval[i], n_kvtab, kvtab );
			fprintf( fp, "        o_kind = %s;\n", dp->d_name );
		}else{
			dp = finddef_withval( rp->r_okval[i], n_kvtab, kvtab );
			fprintf( fp, "      o_kind = %s;\n", dp->d_name );
		}
	}

	fprintf( fp, "      o_type = typetab[t_idx][l_type][r_type];\n" );
	fprintf( fp, "      if( o_type != T_ERROR ){\n" );
	fprintf( fp, "        expr->n_class = o_class;\n" );
	fprintf( fp, "        expr->n_kind = o_kind;\n" );
	fprintf( fp, "        expr->n_type = o_type;\n" );
	fprintf( fp, "      }else{\n" );
	fprintf( fp, "        expr->n_type = T_ERROR;\n" );
	fprintf( fp, "        errormsg( FALSE," );
	fprintf( fp, "\"operator '%s': operands have wrong types.\\n\"",
		rp->r_print );
	fprintf( fp, ");\n" );
	fprintf( fp, "      }\n" );
	fprintf( fp, "      break;\n" );
}

static	void	mk_checkattr( FILE *fp )
{
	int	i;
	ATREC_T	*ap;
	DEF_T	*dp;

	fprintf( fp, "\f\n" );
	fprintf( fp, "typedef struct atrec_t {\n" );
	fprintf( fp, "  char *a_name;\n" );
	fprintf( fp, "  int a_type;\n" );
	fprintf( fp, "  int a_class;\n" );
	fprintf( fp, "  int a_kind;\n" );
	fprintf( fp, "  int a_in;\n" );
	fprintf( fp, "  int a_access;\n" );
	fprintf( fp, "} ATREC_T;\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "static ATREC_T attab[] = {\n" );
	for( ap = attab, i = 0; i < n_attab; i++, ap++ ){
		if( i == 0 )
			fprintf( fp, "  {" ); 
		else
			fprintf( fp, " ,{" ); 
		fprintf( fp, " \"%s\"", ap->a_name );
		dp = finddef_withval( ap->a_type, n_tvtab, tvtab );
		fprintf( fp, ", %s", dp->d_name );
		dp = finddef_withval( ap->a_class, n_cvtab, cvtab );
		fprintf( fp, ", %s", dp->d_name );
		dp = finddef_withval( ap->a_kind, n_kvtab, kvtab );
		fprintf( fp, ", %s", dp->d_name );
		fprintf( fp, ", 0%o", ap->a_in );
		dp = finddef_withval( ap->a_access, n_avtab, avtab );
		fprintf( fp, ", %s", dp->d_name );
		fprintf( fp, " }\n" );
	}
	fprintf( fp, "};\n" );
	fprintf( fp,
		"static int n_attab = sizeof( attab ) / sizeof( ATREC_T );\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "int CG_get_attr_access( expr, npl, npr )\n" );
	fprintf( fp, "NODE_T *expr;\n" );
	fprintf( fp, "NODE_T *npl;\n" );
	fprintf( fp, "NODE_T *npr;\n" );
	fprintf( fp, "{\n" );
	fprintf( fp, "  int i;\n" );
	fprintf( fp, "  ATREC_T *ap;\n" );
	fprintf( fp, "  VALUE_T *vp;\n" );
	fprintf( fp, "  char e_msg[ 256 ];\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "  if( npr->n_sym == SYM_IDENT )\n" );
	fprintf( fp, "    return( A_STRUCT );\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "  vp = &npr->n_val;\n" );
	fprintf( fp, "  ap = attab;\n" );
	fprintf( fp,
		"  for( i = 0; i < n_attab; i++, ap++ ){\n" );
	fprintf( fp, "    if( !strcmp( ap->a_name, vp->v_value.v_cval ) )\n" );
	fprintf( fp, "      return(ap->a_access);\n" );
	fprintf( fp, "  }\n" );
	fprintf( fp, "  sprintf( e_msg, \"unknown attribute '%%s'\\n\"," );
	fprintf( fp, " vp->v_value.v_cval );\n" );
	fprintf( fp, "  expr->n_type = T_ERROR;\n" );
	fprintf( fp, "  errormsg( FALSE, e_msg );\n" );
	dp = finddef( "A_UNDEF", n_avtab, avtab );
	fprintf( fp, "  return(%d);\n", dp->d_val );
	fprintf( fp, "}\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "static void checkattr( expr, npl, npr )\n" );
	fprintf( fp, "NODE_T *expr;\n" );
	fprintf( fp, "NODE_T *npl;\n" );
	fprintf( fp, "NODE_T *npr;\n" );
	fprintf( fp, "{\n" );
	fprintf( fp, "  int i, found;\n" );
	fprintf( fp, "  ATREC_T *ap;\n" );
	fprintf( fp, "  VALUE_T *vp;\n" );
	fprintf( fp, "  NODE_T *npv;\n" );
	fprintf( fp, "  SYMREC_T *str, *fld;\n" );
	fprintf( fp, "  char e_msg[ 256 ];\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "  if( npl->n_type == T_USER ){\n" );
	fprintf( fp, "    for( npv = npl; npv; npv = npv->n_left ){\n" );
	fprintf( fp, "      if( npv->n_sym == SYM_IDENT )\n" );
	fprintf( fp, "        break;\n" );
	fprintf( fp, "    }\n" );
	fprintf( fp, "    str = findsym( npv->n_val.v_value.v_cval );\n" );
	fprintf( fp, "    if( openusyms( str ) ){\n" );
	fprintf( fp, "      fld = findsym( npr->n_val.v_value.v_cval );\n" );
	fprintf( fp, "      if( fld == NULL || fld->s_scope != S_USER ){\n" );
	fprintf( fp,
		"        sprintf(e_msg,\"fld '%%s' not in struct '%%s'\\n\"," );
	fprintf( fp,
		" npr->n_val.v_value.v_cval, npv->n_val.v_value.v_cval );\n" );
	fprintf( fp, "        errormsg(FALSE,e_msg);\n" );
	fprintf( fp, "        return;\n" );
	fprintf( fp, "      }\n" );
	fprintf( fp, "      expr->n_type = fld->s_type;\n" );
	fprintf( fp, "      expr->n_class = fld->s_class;\n" );
	fprintf( fp, "      expr->n_kind = fld->s_kind;\n" );
	fprintf( fp, "      npr->n_sym = SYM_IDENT;\n" );
	fprintf( fp, "      closeusyms( NULL );\n" );
	fprintf( fp, "    }else{\n" );
	fprintf( fp, "    }\n" );
	fprintf( fp, "  }else{\n" );
	fprintf( fp, "    vp = &npr->n_val;\n" );
	fprintf( fp, "    ap = attab;\n" );
	fprintf( fp,
		"    for( found = FALSE, i = 0; i < n_attab; i++, ap++ ){\n" );
	fprintf( fp, "      if( !strcmp( ap->a_name, vp->v_value.v_cval ) ){\n" );
	fprintf( fp, "        found = TRUE;\n" );
	fprintf( fp, "        break;\n" );
	fprintf( fp, "      }\n" );
	fprintf( fp, "    }\n" );
	fprintf( fp, "    if( !found ){\n" );
	fprintf( fp, "      sprintf( e_msg, \"unknown attribute '%%s'\\n\"," );
	fprintf( fp, " vp->v_value.v_cval );\n" );
	fprintf( fp, "      expr->n_type = T_ERROR;\n" );
	fprintf( fp, "      errormsg( FALSE, e_msg );\n" );
	fprintf( fp, "      return;\n" );
	fprintf( fp, "    }\n" );
	fprintf( fp, "    if( npl->n_class != C_VAR ){\n" );
	fprintf( fp, "      expr->n_type = T_ERROR;\n" );
	fprintf( fp,
"      errormsg( FALSE,\"operator '.': left operand must be a variable.\\n\" );\n" );
	fprintf( fp, "      return;\n" );
	fprintf( fp, "    }\n" );
	fprintf( fp, "    if( npl->n_kind != K_SCALAR &&\n" );
	fprintf( fp, "        npl->n_kind != K_DARRAYEL &&\n" );
	fprintf( fp, "        npl->n_kind != K_HASHEL ){\n" );
	fprintf( fp, "      expr->n_type = T_ERROR;\n" );
	fprintf( fp,
"      errormsg( FALSE,\"operator '.': left operand must be a scalar.\\n\" );\n" );
	fprintf( fp, "      return;\n" );
	fprintf( fp, "    }\n" );
	fprintf( fp, "    if( !INSET( npl->n_type, ap->a_in ) ){\n" );
	fprintf( fp, "      expr->n_type = T_ERROR;\n" );
	fprintf( fp,
"      errormsg_s( FALSE,\"operator '.': left operand has no attribute '%%s'.\\n\"," );
	fprintf( fp, " ap->a_name );\n" );
	fprintf( fp, "      return;\n" );
	fprintf( fp, "    }\n" );
	fprintf( fp, "    expr->n_type = ap->a_type;\n" );
	fprintf( fp, "    expr->n_class = ap->a_class;\n" );
	fprintf( fp, "    expr->n_kind = ap->a_kind;\n" );
	fprintf( fp, "  }\n" );
	fprintf( fp, "}\n" );
}

static	int	split( char str[], char *fields[], char *fsep )
{
	int	nf, white;
	size_t flen;
	char	*sp, *fp, *efp, *nfp;

	if( !str )
		return( 0 );

	/* Use fsep of white space is special				*/ 
	/* although a \n is a white space character, a fsep of a	*/
	/* single \n is not considered white space, allowing multi line	*/
	/* strings to be broken into lines				*/

	if( !strcmp( fsep, "\n" ) )
		white = 0;
	else if( strspn( fsep, " \t\n" ) == strlen( fsep ) ){
		white = 1;
	}else
		white = 0;

	if( white ){
		for( nf = 0, sp = str; ; ){
			fp = sp + strspn( sp, fsep );
			if( !*fp )
				return( nf );
			if((efp = strpbrk( fp, fsep ))){
				if( !( flen = efp - fp ) )
					return( nf );
				nfp = (char *)malloc( (flen + 1) * sizeof(char) );
				strncpy( nfp, fp, flen );
				nfp[ flen ] = '\0';
				fields[ nf ] = nfp;
				sp = efp;
				nf++;
			}else{
				flen = strlen( fp );
				nfp = (char *)malloc( (flen + 1) * sizeof(char) );
				strcpy( nfp, fp );
				fields[ nf ] = nfp;
				nf++;
				return( nf );
			}
		}
	}else{
		for( nf = 0, sp = str; ; ){
			if((fp = strchr( sp, *fsep ))){
				flen = fp - sp;
				nfp = (char *)malloc( (flen + 1) * sizeof(char) );
				strncpy( nfp, sp, flen );
				nfp[ flen ] = '\0';
				fields[ nf ] = nfp;
				nf++;
				sp = fp + 1;
			}else{
				flen = strlen( sp );
				nfp = (char *)malloc( (flen + 1) * sizeof(char) );
				strcpy( nfp, sp );
				fields[ nf ] = nfp;
				nf++;
				return( nf );
			}
		}
	}
}
