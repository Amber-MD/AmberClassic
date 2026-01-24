%option prefix="mmo"
%option noyywrap

%{
extern int mytaskid;

#undef ECHO
#define ECHO yyout=nabout;if(mytaskid==0) fprintf(yyout, "\tmm_options:  %s\n",yytext)

char* setVarStr ( char *str);
int setVarInt ( char *str);
double setVarDbl ( char *str);
void setArrayStr ( char *str, char *a, int n, int len);
void setArrayInt ( char *str, int *a, int n);
void setArrayDbl ( char *str, double *a, int n);

%}
Float [+-]?((([0-9]+\.?[0-9]*)|([0-9]*\.?[0-9]+))([eE][-+]?[0-9]+)?)
Int [+-]?[0-9]+
PosInt [+]?[1-9][0-9]*
NegInt [-]?[1-9][0-9]*
NonNegInt [+]?[0-9]+
Str [^,\ =]+
Assign (\ *=\ *|\ +)
ListSep (\ *,\ *)
%%

cut{Assign}{Float}                 { ECHO; cut = setVarDbl(yytext); }
cutnp{Assign}{Float}               { ECHO; cutnp = setVarDbl(yytext); }
scnb{Assign}{Float}                { ECHO; scnb = setVarDbl(yytext); }
scee{Assign}{Float}                { ECHO; scee = setVarDbl(yytext); }
ntpr{Assign}[1-9][0-9]*            { ECHO; ntpr = setVarInt(yytext); 
                                       if( ntpr==0) ntpr = 9999999; }
hcp{Assign}[0124]                  { ECHO; hcp = setVarInt(yytext); }
hcp_h1{Assign}{Float}              { ECHO; hcp_h1 = setVarDbl(yytext); }
hcp_h2{Assign}{Float}              { ECHO; hcp_h2 = setVarDbl(yytext); }
hcp_h3{Assign}{Float}              { ECHO; hcp_h3 = setVarDbl(yytext); }
dhcp{Assign}{Float}                { ECHO; dhcp = setVarDbl(yytext); }

nscm{Assign}{NonNegInt}            { ECHO; nscm = setVarInt(yytext); }
nsnb{Assign}{NonNegInt}            { ECHO; nsnb = setVarInt(yytext); }
nsnp{Assign}{NonNegInt}            { ECHO; nsnp = setVarInt(yytext); }
ntwx{Assign}{NonNegInt}            { ECHO; ntwx = setVarInt(yytext); }
nchk{Assign}{NonNegInt}            { ECHO; nchk = setVarInt(yytext); }
nchk2{Assign}{NonNegInt}           { ECHO; nchk2 = setVarInt(yytext); }

dt{Assign}{Float}                  { ECHO; dt = setVarDbl(yytext); }
t{Assign}{Float}                   { ECHO; t = setVarDbl(yytext); }
tautp{Assign}{Float}               { ECHO; tautp = setVarDbl(yytext); }
temp0{Assign}{Float}               { ECHO; temp0 = setVarDbl(yytext); }
tempi{Assign}{Float}               { ECHO; tempi = setVarDbl(yytext); }
vlimit{Assign}{Float}              { ECHO; vlimit = setVarDbl(yytext); }
gamma_ln{Assign}{Float}            { ECHO; gamma_ln = setVarDbl(yytext); }
genmass{Assign}{Float}             { ECHO; genmass = setVarDbl(yytext); }
wcons{Assign}{Float}               { ECHO; wcons = setVarDbl(yytext); }
ntpr_md{Assign}[1-9][0-9]*         { ECHO; ntpr_md = setVarInt(yytext); 
                                        if( ntpr_md==0) ntpr_md = 9999999; }
zerov{Assign}{NonNegInt}           { ECHO; zerov = setVarInt(yytext); }
k4d{Assign}{Float}                 { ECHO; k4d = setVarDbl(yytext); }

ips{Assign}[0123]                  { ECHO; ips = setVarInt(yytext); }

gb{Assign}[0123578]                { ECHO; gb = setVarInt(yytext); }
gb_debug{Assign}[01]               { ECHO; gb_debug = setVarInt(yytext); }
nr_debug{Assign}[01]               { ECHO; nr_debug = setVarInt(yytext); }
gb2_debug{Assign}[01]              { ECHO; gb2_debug = setVarInt(yytext); }
gbsa_debug{Assign}[01]             { ECHO; gbsa_debug = setVarInt(yytext); }
e_debug{Assign}[0123]              { ECHO; e_debug = setVarInt(yytext); }
epsext{Assign}{Float}              { ECHO; epsext = setVarDbl(yytext); }
kappa{Assign}{Float}               { ECHO; kappa = setVarDbl(yytext); }
rgbmax{Assign}{Float}              { ECHO; rgbmax = setVarDbl(yytext); }

surften{Assign}{Float}             { ECHO; surften = setVarDbl(yytext); }
dradius{Assign}{Float}             { ECHO; dradius = setVarDbl(yytext); }
deltar{Assign}{Float}              { ECHO; deltar = setVarDbl(yytext); }
rwater{Assign}{Float}              { ECHO; rwater = setVarDbl(yytext); }
alphanp{Assign}{Float}             { ECHO; alphanp = setVarDbl(yytext); }
kappanp{Assign}{Float}             { ECHO; kappanp = setVarDbl(yytext); }
min_volume{Assign}{Float}          { ECHO; min_volume = setVarDbl(yytext); }
max_set_size{Assign}{NonNegInt}    { ECHO; max_set_size = setVarInt(yytext); }
MPI_min_tasks{Assign}{NonNegInt}   { ECHO; MPI_min_tasks = setVarInt(yytext); }
cull_np_lists{Assign}[01]          { ECHO; cull_np_lists = setVarInt(yytext); }
use_lower_tri{Assign}[01]          { ECHO; use_lower_tri = setVarInt(yytext); }
dynamic_loops{Assign}[01]          { ECHO; dynamic_loops = setVarInt(yytext); }
blocksize{Assign}{NonNegInt}       { ECHO; blocksize = setVarInt(yytext); }
gbsa{Assign}[01]                   { ECHO; gbsa = setVarInt(yytext); }
dim{Assign}[34]                    { ECHO; dim = setVarInt(yytext); }

diel{Assign}Yammp               { ECHO; dield = -5; }
diel{Assign}R94                 { ECHO; dield = -4; }
diel{Assign}RL94                { ECHO; dield = -3; }
diel{Assign}RL                  { ECHO; dield = -2; }
diel{Assign}R                   { ECHO; dield =  0; }
diel{Assign}C                   { ECHO; dield =  1; }
dielc{Assign}{Float}            { ECHO; dielc = setVarDbl(yytext); }

rattle{Assign}[012]             { ECHO; irattle = setVarInt(yytext); }

[\ \n\t,]+                                                ;

.                               {        
  if(mytaskid==0){
    char c;
    fprintf( stderr, "Error: mm_option: unknown keyword or invalid value: %s", yytext );
    for(c=input();c != ' ' && c != '=' && c != EOF;c=input())
      fprintf( stderr, "%c",c);
    fprintf( stderr,"\n");
  }
  
  exit( 1 );
}
%%

/**
 * Extract a single value from an option statement
 * Takes a string of the form
 * <keyword>=<val1>
 * with an unlimited amount of white space.
 *IN:
 *   str : string 
 *OUT:
 *   returns the string after the keyword
 */
char* setVarStr ( char *str){
  char *copy,*token;
  copy = strdup(str);
  /* the first token is the variable name */
  token=strtok(copy,"= ");
  /* the remaining characters should be the string */
  token = strtok(NULL," ,=");
  return token;
}

/**
 * Extract a single value from an option statement
 * Takes a string of the form
 * <keyword>=<val1>
 * with an unlimited amount of white space.
 *IN:
 *   str : string 
 *OUT:
 *   returns the double after the keyword
 */
double setVarDbl ( char *str){
  double v[1];
  setArrayDbl(str,v,1);
  return v[0];
}

/**
 * Extract a single value from an option statement
 * Takes a string of the form
 * <keyword>=<val1>
 * with an unlimited amount of white space.
 *IN:
 *   str : string 
 *OUT:
 *   returns the integer after the keyword
 */
int setVarInt ( char *str){
  int v[1];
  setArrayInt(str,v,1);
  return v[0];
}

/**
 * Extract mulitple values from an option statement
 * Takes a string of the form
 * <keyword>=<val1>,<val2>,...
 * with an unlimited amount of white space.  Up to n Values are stored
 * in array a.  A 1D character array is used and treated as a 2D array
 * that is [len][n] and must be preallocated.  Each string in the
 * array is set to empty on entry.
 *IN:
 *   str : string 
 *   a   : array of char (len,n)
 *   len : length of 
 *   n   : number of elements of a
 *SIDE EFFECTS:
 *   Values from string are stored in a
 */
void setArrayStr ( char *str, char *a, int len, int n){
  int i;
  char *copy, *token;
  copy = malloc(strlen(str)+1);
  if(copy == NULL){
    fprintf(stderr,"Error: mm_option: failed to allocate memory parsing:\n");
    fprintf(stderr,"'%s'\n",str);
    exit(1);
  }
  strcpy(copy,str);
  /* set all strings in the array to empty */
  for(i=0; i<n; i++){
     *(a+len*i)='\0';
   }

  /* the first token is the variable name */
  token = strtok(copy,"= ");
  /* the remaining tokens should all be floating point numbers */
  for(i = 0; i< n && (token = strtok(NULL," ,=")) != NULL; i++){
        fflush(stderr);
    strncpy((a+len*i), token, len);
  }
  free(copy);  
}

/**
 * Extract mulitple values from an option statement
 * Takes a string of the form
 * <keyword>=<val1>,<val2>,...
 * with an unlimited amount of white space.  Up to n Values are stored in array a.
 *IN:
 *   str : string 
 *   a   : array of doubles
 *   n   : length of a
 *SIDE EFFECTS:
 *   Values from string are stored in a
 */
void setArrayDbl ( char *str, double *a, int n){
  int i;
  char *copy, *token;

  copy = malloc(strlen(str)+1);
  if(copy == NULL){
    fprintf(stderr,"Error: mm_option: failed to allocate memory parsing:\n");
    fprintf(stderr,"'%s'\n",str);
    exit(1);
  }
  strcpy(copy,str);
  /* the first token is the variable name */
  token = strtok(copy,"= ");
  /* the remaining tokens should all be floating point numbers */
  for(i = 0; i< n && (token = strtok(NULL," ,=")) != NULL; i++){
    sscanf(token,"%lg",&(a[i]));
  }
  free(copy);  
}

/**
 * Extract mulitple values from an option statement
 * Takes a string of the form
 * <keyword>=<val1>,<val2>,...
 * with an unlimited amount of white space.  Up to n Values are stored in array a.
 *IN:
 *   str : string 
 *   a   : array of doubles
 *   n   : length of a
 *SIDE EFFECTS:
 *   Values from string are stored in a
 */
void setArrayInt ( char *str, int *a, int n){
  int i;
  char *copy, *token;
  copy = malloc(strlen(str)+1);
  if(copy == NULL){
    fprintf(stderr,"Error: mm_option: failed to allocate memory parsing:\n");
    fprintf(stderr,"'%s'\n",str);
    exit(1);
  }
  strcpy(copy,str);
  /* the first token is the variable name */
  token = strtok(copy,"= ");
  /* the remaining tokens should all be integer numbers */
  for(i = 0; i< n && (token = strtok(NULL," ,=")) != NULL; i++){
    sscanf(token,"%i",&(a[i]));
  }
  free(copy);
}

/***********************************************************************
                            MM_OPTIONS()
************************************************************************/

/* Set the options for mme, md, etc. */

int mm_options(char *opts)
{

   mmo_scan_string( opts );
   mmolex();
   yy_delete_buffer(YY_CURRENT_BUFFER);
   return (0);
}
