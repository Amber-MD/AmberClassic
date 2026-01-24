%option prefix="dgo"
%option noyywrap

%{

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
    char c;
    fprintf( stderr, "Error: dg_options: unknown keyword or invalid value: %s", yytext );
    for(c=input();c != ' ' && c != '=' && c != EOF;c=input())
      fprintf( stderr, "%c",c);
    fprintf( stderr,"\n");
    exit( 1 );
}

%%

int dg_options(BOUNDS_T * bp, char *opts)
{
    dgo_scan_string( opts );
    dgolex();
    yy_delete_buffer(YY_CURRENT_BUFFER);

    natoms = bp->b_natoms;
    bmat = bp->b_bmat;
    nchi = bp->b_nchiral;
    chi = bp->b_chiral;

    return (0);
}
