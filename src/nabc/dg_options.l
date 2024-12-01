%{

#ifdef flex

int dgo_dgoinput( char*, int );

#undef YY_INPUT
#define YY_INPUT(b, r, ms) (r = dgo_yyinput(b, ms))

#else

#undef input
#undef unput
#define input() (*gopts++)
#define unput(c) (*--gopts = c)

#endif

#undef ECHO
#define ECHO yyout=stdout;fprintf(yyout, "\tdg_options:  %s\n",yytext)

%}
%%

ddm[=\ ][^\ \n\t,]+			{ ECHO; ddm = 1; 
								strncpy( ddmn, &yytext[4], yyleng - 4 );
								ddmn[yyleng-4] = '\0'; }
rdm[=\ ][^\ \n\t,]+			{ ECHO; rdm = 1; 
								strncpy( rdmn, &yytext[4], yyleng - 4 );
								rdmn[yyleng-4] = '\0'; }
metrize[=\ ][^\ \n\t]+			{ ECHO; mxp = 1; 
								strncpy( mxpn, &yytext[8], yyleng - 8 );
								mxpn[yyleng-8] = '\0'; }
dmm[=\ ][^\ \n\t,]+			{ ECHO; dmm = 1; 
								strncpy( dmmn, &yytext[4], yyleng - 4 );
								dmmn[yyleng-4] = '\0'; }
rmm[=\ ][^\ \n\t,]+			{ ECHO; rmm = 1; 
								strncpy( rmmn, &yytext[4], yyleng - 4 );
								rmmn[yyleng-4] = '\0'; }
majorize[=\ ][0-9]+			{ ECHO; majorize_opt = 1; niter = atoi( &yytext[9] ); }
gdist[=\ ][0-9]+			{ ECHO; gdist = atoi( &yytext[6] ); }
seed[=\ ]-?[0-9]+			{ ECHO; seed = atoi( &yytext[5] ); }
randpair[=\ ][0-9eE\.\+\-]+	{ ECHO; randpair = atof( &yytext[9] ); }
ntpr[=\ ][0-9]+				{ ECHO; ntpr = atoi( &yytext[5] ); }
eamax[=\ ][0-9]+			{ ECHO; eamax = atoi( &yytext[6] ); }
kchi[=\ ][0-9eE\.\+\-]+		{ ECHO; kchi  = atof( &yytext[5] ); }
k4d[=\ ][0-9eE\.\+\-]+		{ ECHO; k4d  = atof( &yytext[4] ); }
lbpen[=\ ][0-9eE\.\+\-]+	{ ECHO; lbpen  = atof( &yytext[6] ); }
pencut[=\ ][0-9eE\.\+\-]+	{ ECHO; pencut  = atof( &yytext[7] ); }
slearn[=\ ][0-9eE\.\+\-]+	{ ECHO; slearn  = atof( &yytext[7] ); }
sqviol[=\ ][01]				{ ECHO; sqviol = atof( &yytext[7] ); }
pembed[=\ ][01]				{ ECHO; pembed = atoi( &yytext[7] ); }
shuffle[=\ ][01]			{ ECHO; shuffle = atoi( &yytext[8] ); }
rbox[=\ ][0-9eE\.\+\-]+		{ ECHO; rbox  = atof( &yytext[5] ); }
riter[=\ ][0-9]+			{ ECHO; riter = atoi( &yytext[6] ); }

[\ \n\t,]+						;

.						{	
		fprintf( stderr, "unknown option: %s\n", --gopts );
		exit( 1 );
	}

%%

int yywrap(){
	return( 1 );
}

#ifdef flex


int  dgo_yyinput( buf, max_size )

char *buf;
int  max_size;
{
	int n;

	if( max_size > &gopts[dgoinputlim] - dgoinputptr ) 
		n =  &gopts[dgoinputlim] - dgoinputptr;
	else
		n = max_size;

	if( n>0 ){
		memcpy( buf, dgoinputptr, n );
		dgoinputptr += n;
	}
	return n;
}

#endif
