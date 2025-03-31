#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

#include "nab.h"
#include "nabgrm.tab.h"
#include "cgen.h"
#include "errormsg.h"
#include "symbol.h"

#ifdef WIN32
 #include <windows.h>
#else
 #include <unistd.h>
#endif

	/* literals	*/
VALUE_T	v_0;
VALUE_T	v_1;
VALUE_T	v_NULL;
VALUE_T	v_hfirst;
VALUE_T	v_hnext;
VALUE_T	v_mfirst;
VALUE_T	v_mnext;
VALUE_T	v_rfirst;
VALUE_T	v_rnext;
VALUE_T	v_afirst;
VALUE_T	v_anext;
static	VALUE_T	v_ainit;
static	VALUE_T v_sinit;
static	VALUE_T	v_hinit;

#define		ASTKSIZE	20
SYMREC_T	*astk[ ASTKSIZE ];
int		astkp;

extern	char	c_funcname[];

extern	int	cg_dval;
extern	char	cg_nfname[];
extern	int	cg_lineno;
extern	int	cg_noassert;
extern	int	cg_nodebug;

extern  int  errors();

static	char	e_msg[ 256 ];

static	char	codefname[ 256 ];
static	char	ctempfname[ 256 ];
FILE	*cg_cfp = NULL;
static	FILE	*cfp1 = NULL;
static	FILE	*cfp2 = NULL;
static	int	aoption;

static	int	needspace = FALSE;
static	int	needmainstart = TRUE;
static	int	needmainend = FALSE;
static	int	infunc = FALSE;
static	int	idtype;
static	int	idscope = S_GLOBAL;
static	SYMREC_T	*curfunc;
static	NODE_T	*g_hlist = NULL;	/* global hashed arrays for init */
static	NODE_T	*l_hlist = NULL;	/* local hashed arrays for init  */

static	NODE_T	*l_alist = NULL;
static	NODE_T	*l_slist = NULL;
static	NODE_T	*plist = NULL;

	/* temporary expression stuff */
static	int	maxtemp[ N_TYPES ];
static	int	curtemp[ N_TYPES ];
static	NODE_T	*l_dabinit[ 200 ];	/* initializers for parm arrays */
static	NODE_T	*g_dabinit[ 200 ];	/* initializers for parm arrays */
					/* l_dabinit for locals, 	*/
					/* g_dabinit for globals. 	*/
					/* ht = maxtemp[ T_DABOUND ]    */

void	checkexpr( NODE_T * );

SYMREC_T	*CG_genstructdecl( NODE_T *, int, int );
static	void	CG_dumpvars( NODE_T * );
static	char	*CG_genlit( char *, VALUE_T * );
static	char	*CG_genid( char *, VALUE_T *, int );
void	CG_genpids( NODE_T * );

static	int	getdims( NODE_T *, int, int, int *, NODE_T *[] );
static	NODE_T	*mk_arraysize( int, NODE_T *[], int );
static	char	*CG_genexprcode( char *, NODE_T * );
static	void	freeexpr( NODE_T * );
static	char	*CG_gentype( char *, NODE_T * );
static	int	needparen( NODE_T *, NODE_T * );
static	int	op_prec( int );
static	void	fixassert( NODE_T * );
static	void	CG_write_cstr( FILE *, char [] );
static	void	CG_debug( char [], NODE_T *, int );

FILE	*tmpfile();


int	CG_init( char cfname[], int aopt )
{
	aoption = aopt;
	v_0.v_type = T_INT;
	v_0.v_value.v_ival = 0;
	v_1.v_type = T_INT;
	v_1.v_value.v_ival = 1;
	v_NULL.v_type = T_STRING;
	v_NULL.v_value.v_cval = "NULL";
	v_ainit.v_type = T_STRING;
	v_ainit.v_value.v_cval = "NAB_AINIT";
	v_sinit.v_type = T_STRING;
	v_sinit.v_value.v_cval = "NAB_SINIT";
	v_hinit.v_type = T_STRING;
	v_hinit.v_value.v_cval = "NAB_hinit";
	v_hfirst.v_type = T_STRING;
	v_hfirst.v_value.v_cval = "NAB_hfirst";
	v_hnext.v_type = T_STRING;
	v_hnext.v_value.v_cval = "NAB_hnext";
	v_mfirst.v_type = T_STRING;
	v_mfirst.v_value.v_cval = "NAB_mfirst";
	v_mnext.v_type = T_STRING;
	v_mnext.v_value.v_cval = "NAB_mnext";
	v_rfirst.v_type = T_STRING;
	v_rfirst.v_value.v_cval = "NAB_rfirst";
	v_rnext.v_type = T_STRING;
	v_rnext.v_value.v_cval = "NAB_rnext";
	v_afirst.v_type = T_STRING;
	v_afirst.v_value.v_cval = "NAB_afirst";
	v_anext.v_type = T_STRING;
	v_anext.v_value.v_cval = "NAB_anext";
	if( ( cfp1 = fopen( cfname, "w" ) ) == NULL ){
		fprintf( stderr, "can't open c-code file %s\n", cfname );
		return( 1 );
	}

#ifdef WIN32
	char tempFilePath[MAX_PATH];

	GetTempPath(MAX_PATH, ctempfname);
	GetTempFileName(ctempfname, "nab", 0, (LPSTR)(&tempFilePath));
	
    if( ( cfp2 = fopen( tempFilePath, "w+" ) ) == NULL ){
		fprintf( stderr, "can't create c-temp file %s, error: %s\n", tempFilePath, strerror(errno));
		fclose( cfp1 );
		return( 1 );
	}

#else
	strcpy( ctempfname, "/tmp/NABctemp_XXXXXX" );
	int cfd2 = mkstemp( ctempfname );

	if( cfd2 < 0 ){
		fprintf( stderr, "can't create c-temp file %s, error: %s\n", ctempfname, strerror(errno));
		return( 1 );
	}

	if( ( cfp2 = fdopen( cfd2, "w+" ) ) == NULL ){
		fprintf( stderr, "can't open c-temp file %s, error: %s\n", ctempfname, strerror(errno));
		fclose( cfp1 );
		return( 1 );
	}
#endif

	cg_cfp = cfp1;
	strcpy( codefname, cfname );
	fprintf( cg_cfp, "#include <stdio.h>\n" );
	fprintf( cg_cfp, "#include <string.h>\n" );
	fprintf( cg_cfp, "#include <stdlib.h>\n" );
	fprintf( cg_cfp, "#include <math.h>\n" );
	fprintf( cg_cfp, "#include <assert.h>\n" );
	fprintf( cg_cfp, "#include \"nabcode.h\"\n" );
	fprintf( cg_cfp, "extern char NAB_rsbuf[];\n" );
	fprintf( cg_cfp, "static int mytaskid, numtasks;\n" );
	putc( '\n', cg_cfp );
	return( 0 );
}

void	CG_exit( int status )
{
	int	c;

	rewind( cfp2 );
	while( ( c = getc( cfp2 ) ) != EOF )
		putc( c, cfp1 );
/*
	if( status )
		unlink( codefname );
*/
	fclose( cfp1 );
	fclose( cfp2 );
	unlink( ctempfname );
	exit( status );
}

void	CG_genend( void )
{

	if( needmainend ){
		CG_genestmts( FALSE );
#if defined(SCALAPACK) || defined(MPI)
		/* add in an mpifinalize() call:  */
		fprintf( cg_cfp, "\tmpifinalize();\n" );
#endif
		fprintf( cg_cfp, "\n" );
		fprintf( cg_cfp, "\texit( 0 );\n" );
		fprintf( cg_cfp, "}\n" );
	}
	CG_exit( errors() );
}

void	CG_genedefs( int local )
{
	int	t, maxt, type;
	VALUE_T	v_type;
	NODE_T	*npl, *npt, *npd;

	type = local ? T_LDABOUND : T_GDABOUND;
	maxt = maxtemp[ type ];
	curtemp[ type ] = 0;
	maxtemp[ type ] = 0;
	for( t = 0; t < maxt; t++ ){
		if( local )
			npl = node( SYM_LIST, NULL, l_dabinit[ t ], NULL );
		else
			npl = node( SYM_LIST, NULL, g_dabinit[ t ], NULL );
/*
		npt = node( SYM_TYPE, T_INT, NULL, NULL );
*/
		v_type.v_type = T_INT;
		v_type.v_value.v_ival = T_INT;
		npt = node( SYM_TYPE, &v_type, NULL, NULL );
		npd = node( SYM_DECL, 0, npt, npl );
		CG_genvardecl( npd, FALSE, FALSE, FALSE );
	}
	maxtemp[ type ] = 0;
	curtemp[ type ] = 0;
	cg_cfp = cfp2;
}

void	CG_genestmts( int local )
{
	int	c;
	int	t, maxt, type;
	VALUE_T	val, v_type;
	NODE_T	*npi, *npl, *npt, *npd;

	cg_cfp = cfp1;
	val.v_type = T_STRING;
	for( type = 0; type < N_TYPES; type++ ){
		if( local && type == T_GDABOUND )
			continue;
		curtemp[ type ] = 0;
		maxt = maxtemp[ type ];
		for( t = 0; t < maxt; t++ ){
			val.v_value.v_cval = CG_gentemp( type );
			npi = node( SYM_IDENT, &val, NULL, NULL );
			npl = node( SYM_LIST, 0, npi, NULL );
/*
			npt = node( SYM_TYPE, type, NULL, NULL );
*/
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = type;
			npt = node( SYM_TYPE, &v_type, NULL, NULL );
			npd = node( SYM_DECL, 0, npt, npl );
			CG_genvardecl( npd, FALSE, FALSE, FALSE );
		}
		maxtemp[ type ] = 0;
	}

	rewind( cfp2 );
	while( ( c = getc( cfp2 ) ) != EOF )
		putc( c, cfp1 );
	rewind( cfp2 );

	// ftruncate() does not exist in Windows's posix compatibility API
#ifdef WIN32
	HANDLE cfp2Handle = (HANDLE) _get_osfhandle(_fileno(cfp2));

	if(SetFilePointer(cfp2Handle, 0, NULL, FILE_BEGIN) == INVALID_SET_FILE_POINTER)
	{
		fprintf( stderr, "can't set pointer of temp file, error: %lu\n", GetLastError());
	}

	if(!SetEndOfFile(cfp2Handle))
	{
		fprintf( stderr, "can't truncate temp file, error: %lu\n", GetLastError());
	}

#else
	ftruncate( fileno( cfp2 ), (off_t)0 );
#endif
}

void	CG_gennl( void )
{

	cg_lineno++;
	putc( '\n', cg_cfp );
	fflush( cg_cfp );
	needspace = FALSE;
}

void	CG_genmain( void )
{
	static	NODE_T	*npai = NULL, *npsi = NULL;
	static	NODE_T	*npa, *nps;
	NODE_T	*np, *np1;
	SYMREC_T *sp;
	int	t;

	if( !infunc ){
		if( needmainstart ){
			needmainstart = FALSE;
			needmainend = TRUE;
			fprintf( cg_cfp, "int main( argc, argv )\n" );
			fprintf( cg_cfp, "\tint\targc;\n" );
			fprintf( cg_cfp, "\tchar\t*argv[];\n" );
			fprintf( cg_cfp, "{\n" );
			fprintf( cg_cfp, "\tnabout = stdout; /*default*/\n\n" );
#if defined(SCALAPACK) || defined(MPI)
			/* add in an mpiinit() call:  */
			fprintf( cg_cfp, "\tif ( mpiinit(&argc, argv, &mytaskid, &numtasks) != 0 ) {\n" );
			fprintf( cg_cfp, "\t\tprintf(\"Error in mpiinit!\");\n" );
			fprintf( cg_cfp, "\t\tfflush(stdout);\n\t\texit(1);\n\t}\n" );
#else
			fprintf( cg_cfp, "\tmytaskid=0; numtasks=1;\n" );
#endif
			CG_genedefs( FALSE );
		}
	}else {
		if( !npai ){
			/*	Build static (reused) parse tree for
			 *	NAB_AINIT( ar );
			 *	for each local string array
			 */
			npa = node( SYM_IDENT, &v_NULL, NULL, NULL );
			np = node( SYM_LIST, 0, npa, NULL );
			np1 = node( SYM_IDENT, &v_ainit, NULL, NULL );
			npai = node( SYM_CALL, 0, np1, np );
		}
		if( !npsi ){
			/*	Build static (reused) parse tree for
			 *	NAB_SINIT( st );
			 *	for each local struct w/strings
			 */
			nps = node( SYM_IDENT, &v_NULL, NULL, NULL );
			np = node( SYM_LIST, 0, nps, NULL );
			np1 = node( SYM_IDENT, &v_sinit, NULL, NULL );
			npsi = node( SYM_CALL, 0, np1, np );
		}
		/*
		 *	generate code to initialize local string arrays
		 *	to NULL:
		 */

		for( np = l_alist; np; np = np->n_right ){
			np1 = np->n_left->n_right;
			sp = findsym( np1->n_val.v_value.v_cval );
			npa->n_val.v_value.v_cval = sp->s_name;
			CG_genexprcode( NULL, npai );
			CG_genop( NULL, SYM_SEMICOLON );
			CG_gennl();
		}
		l_alist = NULL;

		for( np = l_slist; np; np = np->n_right ){
			np1 = np->n_left;
			sp = findsym( np1->n_val.v_value.v_cval );
			nps->n_val.v_value.v_cval = sp->s_name;
			CG_genexprcode( NULL, npsi );
			CG_genop( NULL, SYM_SEMICOLON );
			CG_gennl();
		}
		l_slist = NULL;
	}

	if( infunc )
		l_hlist = NULL;
	else
		g_hlist = NULL;

	for( t = 0; t < N_TYPES; t++ ){
		if( t != T_CURHASH )
			curtemp[ t ] = 0;
	}
}

void	CG_genvardecl( NODE_T *expr, int isparm, int instruct, int lastparm )
{
	NODE_T	*npt, *npl, *npl0, *npv, *np1;
	NODE_T	*nphd, *npind, *npnull;
	NODE_T	*nplb;
	int	kind, nd, ndim, isdyn;
	SYMREC_T	*sym, *stag = NULL;
	NODE_T	*dims[ 20 ];
	VALUE_T	v_type;

	if( cg_dval & CGD_VARDECL )
		CG_debug( "VAR", expr, 0 );
	npt = expr->n_left;
	idtype = npt->n_type = npt->n_val.v_value.v_ival;
	npt->n_class = C_VAR;
	if( idtype == T_USER ){
		stag = CG_genstructdecl( npt, isparm, lastparm );
		idtype = stag != NULL ? stag->s_type : T_UNDEF;
	}
	nphd = NULL;
	npl0 = expr;
	for( npl = expr->n_right; npl; npl = npl->n_right ){
		kind = K_SCALAR;
		nplb = NULL;
		npv = npl->n_left;
		isdyn = FALSE;
		if( npv->n_sym == SYM_ASSIGN )
			npv = npv->n_left;
		if( npv->n_sym == SYM_LBRACK ){
			nplb = npv;
			np1 = npv->n_right;
			if( np1->n_sym == SYM_LIST ){
				ndim = getdims( np1,
					idscope, isparm, &isdyn, dims );
				nplb->n_right = mk_arraysize(ndim,dims,isparm);
				kind = isdyn ? K_DARRAY : K_ARRAY;
			}else if( np1->n_sym == SYM_HASHED ){
				kind = K_HASHED;
			}
			npv = npv->n_left;
		}
		if( (sym = findsym( npv->n_val.v_value.v_cval )) ){
			if( sym->s_scope == S_SYSTEM ){
				errormsg_s( 0, E_REDEF_BUILTIN_S, sym->s_name );
				continue;
			}else if( sym->s_scope == idscope ){
				errormsg_s( 0, E_REDEF_SYM_S, sym->s_name );
				continue;
			}
		}
		if( !( sym = entersym( idscope, npv->n_val.v_value.v_cval,
			idtype, C_VAR, kind, CC_UNDEF ) ) ){
			sprintf( e_msg, "new symbol %s\n",
				npv->n_val.v_value.v_cval );
			errormsg_s( 1, E_NOMEM_FOR_S, e_msg );
			continue;
		}
		if( idtype == T_USER && stag != NULL ){
			sym->s_user = stag->s_user;
			sym->s_uname = stag->s_name;
			sym->s_init = stag->s_init;
		}
		sym->s_isdyn = isdyn;
		sym->s_isparm = isparm;
		if( kind == K_ARRAY || kind == K_DARRAY ){
			sym->s_parts = ( NODE_T ** )
				malloc( ndim * sizeof( NODE_T * ) );
			for( nd = 0; nd < ndim; nd++ )
				sym->s_parts[ nd ] = dims[ nd ];
			sym->s_pcount = ndim;
		}
		npv->n_type = idtype;
		npv->n_class = C_VAR;
		npv->n_kind = kind;

		/* NOTE: all hashed variables are the same thing	*/  
		/* at the C level: HASH_T *ha, where ha is the hashed	*/
		/* array name. This next if stmt removes the subtree	*/
		/* corresponding to the ha and places it on a hash decl	*/
		/* code for which is generated after the rest (if any)	*/
		/* of the decls in this tree are created.		*/

		if( npv->n_kind == K_HASHED ){
			if( infunc )
				l_hlist = node( SYM_LIST, 0, npv, l_hlist );
			else
				g_hlist = node( SYM_LIST, 0, npv, g_hlist );
			npind = node( SYM_INDIRECT, 0, NULL, npv );
			if( isparm ){
				npind = node( SYM_INDIRECT, 0, NULL, npind );
			}else if( infunc ){
				npnull = node( SYM_IDENT, &v_NULL, NULL, NULL );
				npind = node( SYM_ASSIGN, NULL, npind, npnull );
			}
			nphd = node( SYM_LIST, 0, npind, nphd );
			npl0->n_right = npl->n_right;
			continue;
		}else if( npv->n_kind == K_ARRAY && npv->n_type == T_STRING ){
			if( infunc )
				l_alist = node( SYM_LIST, 0, npv, l_alist );
			npl0 = npl;
		}else
			npl0 = npl;

		if( npv->n_type == T_USER && stag->s_init && infunc )
			l_slist = node( SYM_LIST, 0, npv, l_slist );

		if( npv->n_type == T_STRING ||
			(npv->n_type >= T_FILE &&
			 npv->n_type != T_CURHASH &&
			 npv->n_type != T_USER) )
		{
			np1 = node( SYM_IDENT, &npv->n_val, NULL, NULL );
			np1->n_type = idtype;
			np1->n_class = C_VAR;
			np1->n_kind = kind;
			npv->n_sym = SYM_INDIRECT;
			npv->n_left = NULL;
			npv->n_right = np1;
		}
		if( npv->n_type == T_STRING && npv->n_kind == K_SCALAR &&
			!isparm && !instruct )
		{
			npnull = node( SYM_IDENT, &v_NULL, NULL, NULL );
			np1 = node( SYM_INDIRECT, NULL, NULL, npv->n_right );
			npv->n_sym = SYM_ASSIGN;
			npv->n_left = np1;
			npv->n_right = npnull;
		}
		if( isdyn ){
			nplb->n_sym = SYM_INDIRECT;
			nplb->n_right = nplb->n_left;
			nplb->n_left = NULL;
		}else if( kind != K_ARRAY && isparm ){
			if( npv->n_type != T_POINT && npv->n_type != T_MATRIX ){
				if( npv->n_sym == SYM_INDIRECT ){
					np1 =node(SYM_INDIRECT,0,NULL,
						npv->n_right);
					npv->n_right = np1;
				}else{
					np1 = node(SYM_IDENT,&npv->n_val,
						NULL,NULL);
					np1->n_type = idtype;
					np1->n_class = C_VAR;
					np1->n_kind = kind;
					npv->n_sym = SYM_INDIRECT;
					npv->n_left = NULL;
					npv->n_right = np1;
				}
			}
		}
		if( nplb ){
			nplb->n_type = idtype;
			nplb->n_class = C_VAR;
			nplb->n_kind = kind;
		}
	}
	if( !instruct ){
		if( expr->n_right ){
			if( idscope == S_GLOBAL )
				fprintf( cg_cfp, "static " );
			CG_genexprcode( NULL, expr );
			if( !isparm ){
				CG_genop( NULL, SYM_SEMICOLON );
				CG_gennl();
			}else if( !lastparm )
				CG_genop( NULL, SYM_COMMA );
		}else if( expr->n_left->n_type == T_USER ){
			CG_genexprcode( NULL, expr );
			CG_genop( NULL, SYM_SEMICOLON );
			CG_gennl();
		}
	}
	if( nphd ){
		v_type.v_type = T_INT;
		v_type.v_value.v_ival = T_HASH;
		npt = node( SYM_TYPE, &v_type, NULL, NULL );
		npt->n_type = T_HASH;
		np1 = node( SYM_DECL, 0, npt, nphd );
		CG_genexprcode( NULL, np1 );
		if( !isparm ){
			CG_genop( NULL, SYM_SEMICOLON );
			CG_gennl();
		}else if( !lastparm )
			CG_genop( NULL, SYM_COMMA );
	}
	freeexpr( expr );
}

SYMREC_T	*CG_genstructdecl( NODE_T *expr, int isparm, int lastparm )
{
	NODE_T	*npstag, *npflds, *npl, *npd;
	SYMREC_T	*stag = NULL;
	int		l_idscope;

/*
fprintf( stderr, "CG_genstructdecl: struct decl\n" );
dumpexpr( stderr, expr, 0 );
*/

	npstag = expr->n_left;		/* SYM_STRUCT					*/
	npflds = npstag->n_right;	/* SYM_LIST						*/
	npstag = npstag->n_left;	/* SYM_ID, the struct tag		*/

/*
dumpexpr( stderr, npstag, 0 );
dumpexpr( stderr, npflds, 0 );
*/

	stag = findsym( npstag->n_val.v_value.v_cval );
	if( !stag ){
		if( npflds == NULL ){
			errormsg_s( 0, E_MISSING_FIELDS_S, npstag->n_val.v_value.v_cval );
		}else{	/* this is good, create the stag and do decl the fields	*/
			if( !( stag = entersym( idscope, npstag->n_val.v_value.v_cval,
				idtype, C_STRTAG, K_SCALAR, CC_UNDEF ) ) ){
				sprintf( e_msg, "new symbol %s\n",
					npstag->n_val.v_value.v_cval );
				errormsg_s( 1, E_NOMEM_FOR_S, e_msg );
			}else{
				l_idscope = idscope;
				idscope = S_USER;
				openusyms( stag );
				for( npl = npflds; npl; npl = npl->n_right ){
					npl->n_val.v_type = T_INT; 
					npl->n_val.v_value.v_ival = SYM_SEMICOLON;
					npd = npl->n_left;
/*
fprintf( stderr, "FLD DECL:\n" );
dumpexpr( stderr, npd, 0 );
*/
					CG_genvardecl( npd, isparm, TRUE, lastparm );
				}
				closeusyms( stag );
				idscope = l_idscope;

/*
dsym( stderr, "", stag );
*/

			}
		}
	}else if( stag->s_class != C_STRTAG ){
		errormsg_s( 0, E_STRUCT_TAG_EXPECTED_S, npstag->n_val.v_value.v_cval );
	}else if( npflds != NULL ){
		errormsg_s( 0, E_NOFIELDS_ALLOWED_S, npstag->n_val.v_value.v_cval );
	}else{	/* this is good	*/
	}

	return( stag );
}

void	CG_genassert( NODE_T *expr )
{
	NODE_T	*npr;
	static	char	cbuf[ 10000 ];
	char	*csp;

	if( cg_noassert )
		return;
	if( expr ){
		for( npr = expr->n_right; npr->n_sym == SYM_LPAREN; )
			npr = npr->n_right;
		checkexpr( npr );
		csp = CG_genexprcode( cbuf, npr );
		fixassert( npr );
		fixexpr( npr, FALSE, CC_UNDEF );
		fprintf( cg_cfp, "if( !( " );
		CG_genexprcode( NULL, npr );
		fprintf( cg_cfp, " ) ){\n" );
		fprintf( cg_cfp, "  fprintf( stderr, " );
		fprintf( cg_cfp, "\"%s:%d assert '", cg_nfname, expr->n_lineno );
		CG_write_cstr( cg_cfp, cbuf );
		fprintf( cg_cfp, "' failed.\\n\"" );
		fprintf( cg_cfp, " );\n" );
		CG_dumpvars( npr );
		fprintf( cg_cfp, "  exit( 1 );\n" );
		fprintf( cg_cfp, "}\n" );
	}
}

void	CG_gendebug( NODE_T *expr )
{
	NODE_T	*npl, *npe;
	int	lno;
	static	char	cbuf[ 10000 ];
	char	*csp;
	char	fmt[ 8 ];

	if( cg_nodebug)
		return;
	if( expr ){
		for( lno = UNDEF, npl = expr; npl; npl = npl->n_right ){
			npe = npl->n_left;
			if( npe == NULL )
				continue;
			if( lno == UNDEF )
				lno = npe->n_lineno;
			checkexpr( npe );
			if( npe->n_type == T_INT )
				strcpy( fmt, "%d" );
			else if( npe->n_type == T_SIZE_T )
				strcpy( fmt, "%ld" );
			else if( npe->n_type == T_FLOAT )
				strcpy( fmt, "%lf" );
			else if( npe->n_type == T_STRING )
				strcpy( fmt, "%s" );
			else
				strcpy( fmt, "%p" );

			csp = CG_genexprcode( cbuf, npe );
			fixassert( npe );
			fixexpr( npe, FALSE, CC_UNDEF );
			fprintf( cg_cfp, "fprintf( stderr, " );
			fprintf( cg_cfp, "\"%s:%d '", cg_nfname, lno );
			CG_write_cstr( cg_cfp, cbuf );
			fprintf( cg_cfp, "' = '%s'\\n\", ", fmt );
			CG_genexprcode( NULL, npe );
			fprintf( cg_cfp, " );\n" );
		}
	}
}

void	CG_genexpr( NODE_T *expr )
{

	if( expr ){
		astkp = 0;
		if( cg_dval & CGD_CHKEXPR )
			CG_debug( "before CHK", expr, 0 );
		checkexpr( expr );
		astkp = 0;
		if( cg_dval & CGD_FIXEXPR )
			CG_debug( "before FIX", expr, 0 );
		fixexpr( expr, FALSE, CC_UNDEF );
		astkp = 0;
		if( cg_dval & CGD_EXPRCODE )
			CG_debug( "before GEN", expr, 0 );
		CG_genexprcode( NULL, expr );
		freeexpr( expr );
	}
}

static	char	*CG_genlit( char *csp, VALUE_T *vp )
{

	switch( vp->v_type ){
	case T_INT :
		if( csp ){
			sprintf( csp, "%d", vp->v_value.v_ival );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, "%d", vp->v_value.v_ival );
		break;
	case T_FLOAT :
		if( csp ){
			sprintf( csp, "%E", vp->v_value.v_fval );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, "%E", vp->v_value.v_fval );
		break;
	case T_STRING :
		if( csp ){
			sprintf( csp, "\"%s\"", vp->v_value.v_cval );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, "\"%s\"", vp->v_value.v_cval );
		break;
	}
	return( csp );
}

static	char	*CG_genid( char *csp, VALUE_T *vp, int lookup )
{

	if( needspace ){
		if( csp ){
			*csp++ = ' ';
			*csp = '\0';
		}else
			putc( ' ', cg_cfp );
	}
	if( csp ){
		sprintf( csp, "%s", vp->v_value.v_cval );
		csp += strlen( csp );
	}else
		fprintf( cg_cfp, "%s", vp->v_value.v_cval );
	needspace = TRUE;
	return( csp );
}

void	checkid( NODE_T *np )
{
	SYMREC_T	*sp;

	if( (sp = findsym( np->n_val.v_value.v_cval )) ){
		np->n_type = sp->s_type;
		np->n_class = sp->s_class;
		np->n_kind = sp->s_kind;
	}else
		errormsg_s( TRUE, E_NO_DECL_S, np->n_val.v_value.v_cval );
}

char	*CG_genop( char *csp, int op )
{

	switch( op ){
	case SYM_ADDRESS :
		if( csp ){
			sprintf( csp, " &" );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " &" );
		break;
	case SYM_AND :
		if( csp ){
			sprintf( csp, " && " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " && " );
		break;
	case SYM_ASSIGN :
		if( csp ){
			sprintf( csp, " = " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " = " );
		break;
	case SYM_ATSIGN :
		if( csp ){
			sprintf( csp, " @ " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " @ " );
		break;
	case SYM_COMMA :
		if( csp ){
			sprintf( csp, ", " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, ", " );
		break;
	case SYM_DONT_MATCH :
		if( csp ){
			sprintf( csp, " !~ " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " !~ " );
		break;
	case SYM_EQUAL :
		if( csp ){
			sprintf( csp, " == " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " == " );
		break;
	case SYM_GREATER_EQUAL :
		if( csp ){
			sprintf( csp, " >= " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " >= " );
		break;
	case SYM_GREATER :
		if( csp ){
			sprintf( csp, " > " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " > " );
		break;
	case SYM_IN :
		if( csp ){
			sprintf( csp, " in " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " in " );
		break;
	case SYM_INDEX :	/* FIX ME!!! */
		break;
	case SYM_INDIRECT :
		if( csp ){
			sprintf( csp, " *" );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " *" );
		break;
	case SYM_LBRACE :
		if( csp ){
			sprintf( csp, "{" );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, "{" );
		break;
	case SYM_LBRACK :
		if( csp ){
			sprintf( csp, "[" );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, "[" );
		break;
	case SYM_LESS_EQUAL :
		if( csp ){
			sprintf( csp, " <= " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " <= " );
		break;
	case SYM_LESS :
		if( csp ){
			sprintf( csp, " < " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " < " );
		break;
	case SYM_LPAREN :
		if( csp ){
			sprintf( csp, "( " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, "( " );
		break;
	case SYM_MATCH :
		if( csp ){
			sprintf( csp, " =~ " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " =~ " );
		break;
	case SYM_MINUS :
		if( csp ){
			sprintf( csp, " - " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " - " );
		break;
	case SYM_MINUS_ASSIGN :
		if( csp ){
			sprintf( csp, " -= " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " -= " );
		break;
	case SYM_MINUS_MINUS :
		if( csp ){
			sprintf( csp, " -- " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " -- " );
		break;
	case SYM_MODULUS :
		if( csp ){
			sprintf( csp, " %% " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " %% " );
		break;
	case SYM_MODULUS_ASSIGN :
		if( csp ){
			sprintf( csp, " %%= " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " %%= " );
		break;
	case SYM_NEGATE :
		if( csp ){
			sprintf( csp, " - " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " - " );
		break;
	case SYM_NOT :
		if( csp ){
			sprintf( csp, " !" );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " !" );
		break;
	case SYM_NOT_EQUAL :
		if( csp ){
			sprintf( csp, " != " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " != " );
		break;
	case SYM_OR :
		if( csp ){
			sprintf( csp, " || " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " || " );
		break;
	case SYM_PARM :	/* no code required	*/
		break;
	case SYM_PERIOD :
		if( csp ){
			sprintf( csp, " . " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " . " );
		break;
	case SYM_PLUS :
		if( csp ){
			sprintf( csp, " + " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " + " );
		break;
	case SYM_PLUS_ASSIGN :
		if( csp ){
			sprintf( csp, " += " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " += " );
		break;
	case SYM_PLUS_PLUS :
		if( csp ){
			sprintf( csp, " ++ " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " ++ " );
		break;
	case SYM_POINTS_TO :
		if( csp ){
			sprintf( csp, "->" );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, "->" );
		break;
	case SYM_RBRACE :
		if( csp ){
			sprintf( csp, "}" );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, "}" );
		break;
	case SYM_RBRACK :
		if( csp ){
			sprintf( csp, "]" );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, "]" );
		break;
	case SYM_RPAREN :
		if( csp ){
			sprintf( csp, " )" );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " )" );
		break;
	case SYM_SEMICOLON :
		if( csp ){
			sprintf( csp, ";" );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, ";" );
		break;
	case SYM_STAR :
		if( csp ){
			sprintf( csp, " * " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " * " );
		break;
	case SYM_STAR_ASSIGN :
		if( csp ){
			sprintf( csp, " *= " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " *= " );
		break;
	case SYM_SLASH :
		if( csp ){
			sprintf( csp, " / " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " / " );
		break;
	case SYM_SLASH_ASSIGN :
		if( csp ){
			sprintf( csp, " /= " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " /= " );
		break;
	case SYM_UPARROW :
		if( csp ){
			sprintf( csp, " ^ " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " ^ " );
		break;
	case SYM_UPARROW_ASSIGN :
		if( csp ){
			sprintf( csp, " ^= " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " ^= " );
		break;
	case SYM_STRUCT :
		if( csp ){
			sprintf( csp, " struct " );
			csp += strlen( csp );
		}else
			fprintf( cg_cfp, " struct " );
		break;
	default :
		fprintf( stderr, "CG_genop: unexpected operator %d\n", op );
		exit( 1 );
		break;
	}
	needspace = FALSE;
	return( csp );
}

void	CG_genrword( int sym )
{

	if( needspace )
		putc( ' ', cg_cfp );
		
	switch( sym ){
	case SYM_STRUCT :
		fprintf( cg_cfp, "struct" );
		break;
	case SYM_BREAK :
		fprintf( cg_cfp, "break" );
		break;
	case SYM_CONTINUE :
		fprintf( cg_cfp, "continue" );
		break;
	case SYM_IF :
		fprintf( cg_cfp, "if" );
		break;
	case SYM_ELSE :
		fprintf( cg_cfp, "else" );
		break;
	case SYM_WHILE :
		fprintf( cg_cfp, "while" );
		break;
	case SYM_FOR :
		fprintf( cg_cfp, "for" );
		break;
	case SYM_RETURN :
		fprintf( cg_cfp, "return" );
		break;
	default :
		fprintf( stderr,
			"CG_genrword: symbol %d is not a reserved word\n",
			sym );
		exit( 1 );
		break;
	}
	needspace = TRUE;
}

void	CG_genfhdr( NODE_T *ftype, NODE_T *fname )
{
	char	*fnp;
	int	type, kind;
	SYMREC_T	*sp;
	NODE_T	*np;

	/*
	 *	This code will need to worked on some more to work 
 	 * 	correctly with the apparant double declaration of a
	 *	function that must be declared and defined.
	 */
	type = ftype->n_type = ftype->n_val.v_value.v_ival;
	kind = K_SCALAR;
	CG_genexprcode( NULL, ftype );
	if( ( sp = findsym( fname->n_val.v_value.v_cval ) ) ){
		/* work here required */
	}else if( !( sp = entersym( S_GLOBAL, fname->n_val.v_value.v_cval,
		type, C_FUNC, kind, CC_UNDEF ) ) ){
		/* ERROR */
	}
	curfunc = sp;
	fname->n_type = type;
	fname->n_class = C_FUNC;
	fname->n_kind = kind;
	if( type == T_INT || type == T_SIZE_T || type == T_FLOAT )
		CG_genexprcode( NULL, fname );
	else{
		np = node( SYM_INDIRECT, 0, NULL, fname );
		CG_genexprcode( NULL, np );
	}
	fnp = fname->n_val.v_value.v_cval;
	strcpy( c_funcname, fnp );
	idscope = S_LOCAL;
}

void	CG_genplist( NODE_T *np )
{

	plist = copynode( np );
	CG_genop( NULL, SYM_LPAREN );
/*	CG_genpids( np );	*/
	CG_genpdecl( plist );
	CG_genop( NULL, SYM_RPAREN );
}

void	CG_genpdecls( void )
{

/*
	CG_genpdecl( plist );
*/
}

void	CG_genfstart( void )
{

	curfunc->s_cconv = CC_NAB;
	infunc = TRUE;
}

void	CG_genfend( NODE_T *expr )
{
	char	*ip;

	if( expr ){
		ip = expr->n_val.v_value.v_cval;
		if( !strcmp( ip, "nab" ) )
			curfunc->s_cconv = CC_NAB;
		else if( !strcmp( ip, "c" ) )
			curfunc->s_cconv = CC_CC;
		else if( !strcmp( ip, "fortran" ) )
			curfunc->s_cconv = CC_FORTRAN;
		else{
			fprintf( stderr,
				"CG_genfend: unknown function class %s\n", ip );
		}
	}else if( curfunc->s_cconv == CC_UNDEF )
		curfunc->s_cconv = CC_NAB;
	idscope = S_GLOBAL;
	freelsyms();
	plist = NULL;
	infunc = FALSE;
	curfunc = NULL;
	strcpy( c_funcname, "main" );
}

void	CG_genpids( NODE_T *np )
{

	if( np ){
		CG_genpids( np->n_left );
		if( np->n_sym == SYM_IDENT )
			CG_genid( NULL, &np->n_val, FALSE );
		else if( np->n_sym == SYM_LIST && np->n_right )
			CG_genop( NULL, SYM_COMMA );
		if( np->n_sym != SYM_LBRACK )
			CG_genpids( np->n_right );
	}
}

void	CG_genpdecl( NODE_T *np )
{
	NODE_T	*npp;

	/*
	 *	np points to the parameter list, or at list what remains
	 * 	of.  The tree holding the plist behaves like a list: np points
	 *	to a node of type SYM_LIST; the right subtree is the next parm
	 *	(if it exists); the left subtree is SYM_TYPE decl.  Thus a loop
	 *	can be used to tranverse the plist. 
	 *
	 */
	for( npp = np; npp; npp = npp->n_right ){
		CG_genvardecl( npp->n_left, TRUE, FALSE, npp->n_right == NULL );
	}
}

char	*CG_gentemp( int type )
{
	char	*tnp;
	int	tn;

	tn = ++curtemp[ type ];
	if( tn > maxtemp[ type ] )
		maxtemp[ type ] = tn;
	tnp = ( char * )malloc( 16 * sizeof(char) );

	switch( type ){
	case T_INT :
		sprintf( tnp, "__it%04d__", tn );
		break;
	case T_SIZE_T :
		sprintf( tnp, "__szt%04d__", tn );
		break;
	case T_FLOAT :
		sprintf( tnp, "__ft%04d__", tn );
		break;
	case T_STRING :
		sprintf( tnp, "__st%04d__", tn );
		break;
	case T_FILE :
		sprintf( tnp, "__fpt%04d__", tn );
		break;
	case T_CURHASH :
		sprintf( tnp, "__cht%04d__", tn );
		break;
	case T_POINT :
		sprintf( tnp, "__pt%04d__", tn );
		break;
	case T_MATRIX :	
		sprintf( tnp, "__mat%04d__", tn );
		break;
	case T_LDABOUND :
		sprintf( tnp, "__ldab%04d__", tn );
		break;
	case T_GDABOUND :
		sprintf( tnp, "__gdab%04d__", tn );
		break;
	case T_UNDEF :
	case T_ATOM :
	case T_RESIDUE :
	case T_MOLECULE :
	case T_BOUNDS :
	case T_HASH :
	default :
		fprintf( stderr,
		"CG_gentemp: internal error: can't create temp of type %d\n",
			type );
		exit( 1 );
		break;
	}
	return( tnp );
}

static	int	getdims( NODE_T * expr, int scope, int isparm, int *isdyn,
	NODE_T *dims[] )
{
	int	nd;
	NODE_T	*npl, *npd, *npda, *np1;
	VALUE_T	dab;

	*isdyn = FALSE;
	for( nd = 0, npl = expr; npl; npl = npl->n_right ){
		npd = npl->n_left;
		if( npd->n_sym == SYM_INT_LIT )
			dims[ nd ] = node( SYM_INT_LIT, &npd->n_val, 0, 0 );
		else if( npd->n_sym == SYM_IDENT ){
			np1 = node( SYM_IDENT, &npd->n_val, 0, 0 );
			np1 = node( SYM_INDIRECT, 0, 0, np1 );
			dab.v_value.v_cval = CG_gentemp( T_LDABOUND );
			npda = node( SYM_IDENT, &dab, 0, 0 );
			dims[ nd ] = npda;
			np1 = node( SYM_ASSIGN, 0, npda, np1 );
			l_dabinit[ maxtemp[ T_LDABOUND ] - 1 ] = np1;
		}else if( npd->n_sym == SYM_DYNAMIC ){
			*isdyn = TRUE;
			if( scope == S_LOCAL )
				dab.v_value.v_cval = CG_gentemp( T_LDABOUND );
			else
				dab.v_value.v_cval = CG_gentemp( T_GDABOUND );
			dims[ nd ] = npda = node( SYM_IDENT, &dab, 0, 0 );
			npda->n_type = T_INT;
			npda->n_class = C_VAR;
			npda->n_kind = K_SCALAR;
			if( scope == S_LOCAL )
				l_dabinit[ maxtemp[ T_LDABOUND ] - 1 ] = npda;
			else
				g_dabinit[ maxtemp[ T_GDABOUND ] - 1 ] = npda;
			if( isparm ){
				/* ERROR: parm arrays can not be dynamic */
			}
		}else{
			/* ERROR: array bounds must be int/id/dynamic */
		}
		nd++;
	}
	return( nd );
}

static	NODE_T	*mk_arraysize( int ndim, NODE_T *dims[], int isparm )
{
	int	nd;
	NODE_T	*npas;

	if( isparm )
		return( NULL );
	for( nd = 0; nd < ndim; nd++ ){
		if( dims[ nd ]->n_sym == SYM_IDENT ){
			return( NULL );
		}
	}
	if( ndim == 1 )
		return( dims[ 0 ] );
	npas = dims[ 0 ];
	for( nd = 1; nd < ndim; nd++ ){
		npas = node( SYM_STAR, 0, npas, dims[ nd ] );
	}
	return( npas );
}

static	char	*CG_genexprcode( char *csp, NODE_T *expr )
{
	int	n_rparen;

	if( expr ){
		n_rparen = 0;
		if( expr->n_sym == SYM_CALL ){
			csp = CG_genid( csp, &expr->n_left->n_val, FALSE );
			csp = CG_genop( csp, SYM_LPAREN );
		}else if( expr->n_sym == SYM_LPAREN )
			csp = CG_genop( csp, SYM_LPAREN );
		else if( expr->n_sym == SYM_LBRACK ){
			csp = CG_genexprcode( csp, expr->n_left );
			csp = CG_genop( csp, SYM_LBRACK );
		}else if( expr->n_sym == SYM_STRUCT ){
			csp = CG_genop( csp, SYM_STRUCT );
			csp = CG_genexprcode( csp, expr->n_left );
			if( expr->n_right != NULL )
				csp = CG_genop( csp, SYM_LBRACE );
		}else{
			if( (n_rparen = needparen( expr, expr->n_left )) )
				csp = CG_genop( csp, SYM_LPAREN );
			csp = CG_genexprcode( csp, expr->n_left );
			if( n_rparen )
				csp = CG_genop( csp, SYM_RPAREN );
		}
		switch( expr->n_sym ){

		/* symbols for the statement reserved words:	*/
		case SYM_BREAK :
		case SYM_CONTINUE :
		case SYM_ELSE :
		case SYM_FOR :
		case SYM_FOREACH :
		case SYM_IF :
		case SYM_RETURN :
		case SYM_WHILE :
			break;

		/* special nab stmts that are both stmts &	*/
		/* exprs but have no direct C xlation:		*/
		case SYM_ALLOCATE :
		case SYM_ASSERT :
		case SYM_DEALLOCATE :
		case SYM_DEBUG :
		case SYM_DELETE :
			break;

		/* function call: is handled both pre & postfix	*/
		/* but not infix.				*/
		case SYM_CALL :
			break;

		/* this symbol generates no code, unlike LIST	*/
		/* which may generate a comma			*/
		case SYM_STMTLIST :
			break;

		/* symbols involved in declarations:		*/
		case SYM_DECL :
			break;

		/* symbols that represent the nab types:	*/ 
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
			break;

		/* special array type symbols:			*/
		case SYM_DYNAMIC :
		case SYM_HASHED :
			break;

		case SYM_TYPE :
			csp = CG_gentype( csp, expr );
			break;

		case SYM_IDENT :
		case SYM_ATTRIBUTE :
			csp = CG_genid( csp, &expr->n_val, FALSE );
			break;

		case SYM_INT_LIT :
		case SYM_FLOAT_LIT :
		case SYM_STRING_LIT :
			csp = CG_genlit( csp, &expr->n_val );
			break;

		case SYM_ADDRESS :
		case SYM_AND :
		case SYM_ASSIGN :
		case SYM_COMMA :
		case SYM_DONT_MATCH :
		case SYM_EQUAL :
		case SYM_GREATER :
		case SYM_GREATER_EQUAL :
		case SYM_IN :
		case SYM_INDIRECT :
		case SYM_LESS :
		case SYM_LESS_EQUAL :
		case SYM_LBRACE :
		case SYM_MATCH :
		case SYM_MINUS :
		case SYM_MINUS_ASSIGN :
		case SYM_MINUS_MINUS :
		case SYM_MODULUS :
		case SYM_MODULUS_ASSIGN :
		case SYM_NEGATE :
		case SYM_NOT :
		case SYM_NOT_EQUAL :
		case SYM_OR :
		case SYM_PARM :
		case SYM_PERIOD :
		case SYM_PLUS :
		case SYM_PLUS_ASSIGN :
		case SYM_PLUS_PLUS :
		case SYM_POINTS_TO :
		case SYM_RBRACE :
		case SYM_SEMICOLON :
		case SYM_SLASH :
		case SYM_SLASH_ASSIGN :
		case SYM_STAR :
		case SYM_STAR_ASSIGN :
		case SYM_UPARROW :
		case SYM_UPARROW_ASSIGN :
			csp = CG_genop( csp, expr->n_sym );
			break;

		case SYM_INDEX :
			if( expr->n_right )
				csp = CG_genop( csp, SYM_COMMA );
			break;
		case SYM_LIST :
			if( expr->n_val.v_type == T_INT )
				csp = CG_genop( csp, expr->n_val.v_value.v_ival );
			else if( expr->n_right )
				csp = CG_genop( csp, SYM_COMMA );
			break;

		/* bracketing symbols are generated pre/post	*/
		case SYM_LPAREN :
		case SYM_LBRACK :
			break;

		/* struct is generated pre/post */
		case SYM_STRUCT :
			break;

		/* top level operator in if(), while() and for( ;; )	*/
		/* generates no code but provides a `hook' to attach	*/
		/* the rewrite rule for tests of point & matrix exprs	*/
		case SYM_TEST :
			break;

		case SYM_ERROR :
			break;
		default :
			fprintf( stderr,
				"CG_genexprcode: unexpected symbol: %d\n",
				expr->n_sym );
			exit( 1 );
			break;
		}
		if( expr->n_sym == SYM_CALL ){
			csp = CG_genexprcode( csp, expr->n_right );
			csp = CG_genop( csp, SYM_RPAREN );
		}else if( expr->n_sym == SYM_LPAREN ){
			csp = CG_genexprcode( csp, expr->n_right );
			csp = CG_genop( csp, SYM_RPAREN );
		}else if( expr->n_sym == SYM_LBRACK ){
			csp = CG_genexprcode( csp, expr->n_right );
			csp = CG_genop( csp, SYM_RBRACK );
		}else if( expr->n_sym == SYM_STRUCT ){
			if( expr->n_right ){
				csp = CG_genexprcode( csp, expr->n_right );
				csp = CG_genop( csp, SYM_RBRACE );
			}
		}else{
			if( (n_rparen = needparen( expr, expr->n_right )) )
				csp = CG_genop( csp, SYM_LPAREN );
			csp = CG_genexprcode( csp, expr->n_right );
			if( n_rparen )
				csp = CG_genop( csp, SYM_RPAREN );
		}
	}
	return( csp );
}

static	char	*CG_gentype( char *csp, NODE_T *npt )
{
	char	*tname;

	tname = "";
	switch( npt->n_type ){
	case T_INT :
		tname = "INT_T";
		break;
	case T_SIZE_T :
		tname = "SIZE_T";
		break;
	case T_FLOAT :
		tname = "REAL_T";
		break;
	case T_STRING :
		tname = "STRING_T";
		break;
	case T_POINT :
		tname = "POINT_T";
		break;
	case T_MATRIX :
		tname = "MATRIX_T";
		break;
	case T_FILE :
		tname = "FILE_T";
		break;
	case T_ATOM :
		tname = "ATOM_T";
		break;
	case T_RESIDUE :
		tname = "RESIDUE_T";
		break;
	case T_MOLECULE :
		tname = "MOLECULE_T";
		break;
	case T_BOUNDS :
		tname = "BOUNDS_T";
		break;
	case T_HASH :
		tname = "HASH_T";
		break;
	case T_CURHASH :
		tname = "CURHASH_T";
		break;
	case T_LDABOUND :
		tname = "INT_T";
		break;
	case T_GDABOUND :
		tname = "INT_T";
		break;
	}
	if( csp ){
		strcpy( csp, tname );
		csp += strlen( csp );
	}else
		fputs( tname, cg_cfp );
	needspace = TRUE;
	return( csp );
}

static	int	needparen( NODE_T *expr, NODE_T *s_expr )
{
	int	p, sp;

	if( !expr || !s_expr )
		return( FALSE );
	else if( (p=op_prec( expr->n_sym )) > (sp=op_prec( s_expr->n_sym )) )
		return( TRUE );
	else
		return( FALSE );
}

static	int	op_prec( int sym )
{

	switch( sym ){
	case SYM_ASSIGN :
	case SYM_PLUS_ASSIGN :
	case SYM_MINUS_ASSIGN :
	case SYM_STAR_ASSIGN :
	case SYM_SLASH_ASSIGN :
	case SYM_MODULUS_ASSIGN :
		return( 0 );
	case SYM_COMMA :
		return( 0 );

	case SYM_OR :
		return( 1 );

	case SYM_AND :
		return( 2 );

	case SYM_EQUAL :
	case SYM_NOT_EQUAL :
		return( 3 );

	case SYM_LESS :
	case SYM_LESS_EQUAL :
	case SYM_GREATER_EQUAL :
	case SYM_GREATER :
		return( 4 );

	case SYM_PLUS :
	case SYM_MINUS :
		return( 5 );

	case SYM_STAR :
	case SYM_SLASH :
	case SYM_MODULUS :
	case SYM_ATSIGN :
		return( 6 );

	case SYM_UPARROW :
		return( 7 );

	case SYM_PLUS_PLUS :
	case SYM_MINUS_MINUS :
	case SYM_NEGATE :
	case SYM_NOT :

/* proposed:*/
	case SYM_ADDRESS :
	case SYM_INDIRECT :
/* proposed:*/
		return( 8 );

/* orig
	case SYM_PERIOD :
*/
	case SYM_LPAREN :
	case SYM_RPAREN :
	case SYM_LBRACK :
	case SYM_RBRACK :
	case SYM_LBRACE :
	case SYM_RBRACE :
/* orig
	case SYM_ADDRESS :
	case SYM_INDIRECT :
*/
		return( 9 );

/* proposed:*/
	case SYM_PERIOD :
/* proposed:*/
	case SYM_POINTS_TO :
		return( 10 );

	case SYM_INT_LIT :
	case SYM_FLOAT_LIT :
	case SYM_STRING_LIT :
	case SYM_IDENT :
	case SYM_ATTRIBUTE :
		return( 11 );

	default :
		return( 0 );
	}
}

static	void	freeexpr( NODE_T *expr )
{

	return;
/*
	if( expr ){
		freeexpr( expr->n_left );
		freeexpr( expr->n_right );
		if( expr->n_val.v_type == T_STRING )
			free( expr->n_val.v_value.v_cval );	
		free( expr );
	}
*/
}

static	void	fixassert( NODE_T *expr )
{
	static	char	cbuf[ 10000 ];
	char	*csp;

	if( expr ){
		if( expr->n_sym == SYM_PERIOD ){
			CG_genexprcode( cbuf, expr );
			csp = ( char * )malloc( strlen( cbuf ) + 1 );
			strcpy( csp, cbuf );
			expr->n_val.v_type = T_STRING;
			expr->n_val.v_value.v_cval = csp;
		}else if( expr->n_sym == SYM_LBRACK ){
			CG_genexprcode( cbuf, expr );
			csp = ( char * )malloc( strlen( cbuf ) + 1 );
			strcpy( csp, cbuf );
			expr->n_val.v_type = T_STRING;
			expr->n_val.v_value.v_cval = csp;
			fixassert( expr->n_right );
		}else{
			fixassert( expr->n_left );
			fixassert( expr->n_right );
		}
	}
}

static	void	CG_write_cstr( FILE *fp, char str[] )
{
	char	*sp;

	for( sp = str; *sp; sp++ ){
		if( *sp == '\\' )
			putc( '\\', cg_cfp );
		else if( *sp == '"' )
			putc( '\\', cg_cfp );
		putc( *sp, cg_cfp );
	}
}

static	void	CG_dumpvars( NODE_T *expr )
{
	char	fmt[ 8 ];

	if( expr ){
		if( expr->n_type == T_INT )
			strcpy( fmt, "%d" );
		else if( expr->n_type == T_SIZE_T )
			strcpy( fmt, "%ld" );
		else if( expr->n_type == T_FLOAT )
			strcpy( fmt, "%lf" );
		else if( expr->n_type == T_STRING )
			strcpy( fmt, "'%s'" );
		else
			*fmt = '\0';
		if( expr->n_sym == SYM_LBRACK ){
			if( *fmt ){
				fprintf( cg_cfp, "  fprintf( stderr, " );
				fprintf( cg_cfp, "\"%s:%d        '",
					cg_nfname, expr->n_lineno );
				CG_write_cstr( cg_cfp,
					expr->n_val.v_value.v_cval );
				fprintf( cg_cfp, "' = %s\\n\", ", fmt );
				CG_genexprcode( NULL, expr );
				fprintf( cg_cfp, " );\n" );
			}
			CG_dumpvars( expr->n_right );
		}else if( expr->n_sym == SYM_INDIRECT ){
			if( *fmt ){
				fprintf( cg_cfp, "  fprintf( stderr, " );
				fprintf( cg_cfp, "\"%s:%d        '",
					cg_nfname, expr->n_lineno );
				CG_write_cstr( cg_cfp,
					expr->n_val.v_value.v_cval );
				fprintf( cg_cfp, "' = %s\\n\", ", fmt );
				CG_genexprcode( NULL, expr );
				fprintf( cg_cfp, " );\n" );
			}
		}else if( expr->n_sym == SYM_CALL ){
			if( expr->n_class == C_VAR ){
				if( *fmt ){
					fprintf( cg_cfp, "  fprintf( stderr, " );
					fprintf( cg_cfp, "\"%s:%d        '",
						cg_nfname, expr->n_lineno );
					CG_write_cstr( cg_cfp,
						expr->n_val.v_value.v_cval );
					fprintf( cg_cfp, "' = %s\\n\", ", fmt );
					CG_genexprcode( NULL, expr );
					fprintf( cg_cfp, " );\n" );
				}
				CG_dumpvars( expr->n_right );
			}
		}else if( expr->n_sym == SYM_IDENT ){
			if( expr->n_kind == K_SCALAR ){
				if( *fmt ){
					fprintf( cg_cfp, "  fprintf( stderr, " );
					fprintf( cg_cfp, "\"%s:%d        '",
						cg_nfname, expr->n_lineno );
					CG_write_cstr( cg_cfp,
						expr->n_val.v_value.v_cval );
					fprintf( cg_cfp, "' = %s\\n\", ", fmt );
					CG_genexprcode( NULL, expr );
					fprintf( cg_cfp, " );\n" );
				}
			}
		}else{
			CG_dumpvars( expr->n_left );
			CG_dumpvars( expr->n_right );
		}
	}
}

static	void	CG_debug( char msg[], NODE_T *expr, int offset )
{
	static	FILE	*cg_dfile = NULL;

	if( !cg_dfile ){
		if( !( cg_dfile = fopen( "debug.file", "w" ) ) )
			return;
	}
	fprintf( cg_dfile, "%s:%d %s\n", cg_nfname, cg_lineno, msg );
	fflush( cg_dfile );
	dumpexpr( cg_dfile, expr, offset );
	fflush( cg_dfile );
}
