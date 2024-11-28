/*  Written by Garry Gippert, The Scripps Research Institute              */

/*   #include <sys/times.h>  */
#include <stdio.h>
#include <string.h>
#include <math.h>
#define _STANDALONE /* SGI R10K bug ... */
#include <time.h>
#undef  _STANDALONE
#include <assert.h>

/* shorthand definitions and "macros" here */
#define ABS( A )	( ( A ) > 0 ? ( A ) : - ( A ) )
#define MIN( A, B )	( ( A ) < ( B ) ? ( A ) : ( B ) )
#define MAX( A, B )	( ( A ) > ( B ) ? ( A ) : ( B ) )
#define SIGN( A, B )	( ( B ) <  0  ? -fabs( A ) : fabs( A ) )
#define SQR( A )	( ( A ) * ( A ) )

/* architecture and similar dependencies here */
#ifdef CRAY
#define FSIZE 8
#else
#define FSIZE 4
#endif

/* local ugly stuff */
#define getline        	(fgets(line,sizeof line,fp)!=NULL)
#define lineis(fword)	(strncmp(line,(fword),strlen((fword)))==0)
#define	blankline	((slen=strlen( line )) <= 1)

/* molecular, constraints and real-space limits here */
#define MAXATOMS	50000
#define	MAXDEGEN	400

#define BIGDIST		9999.9999
#define	TSTHRESHOLD	1.0e-1
#define	SMALLDIST	1.0e-4
#define	VERYSMALLDIST	1.0e-8
#define PI		 3.1415926535897932384626433832795028841972

/* dimensionality here */
#define MAXDIMEN	3
#define X		0
#define Y		1
#define Z		2

/* bookkeeping, other non-molecular and trivial definitions here */
#define FALSE		0
#define TRUE		1
#define	ONE		1
#define TWO		2
#define	THREE		3
#define	FOUR		4
#define	MAXWORDS	200
#define WORDSIZE	12
#define	LINESIZE	1024
#define	COMMENTSIZE	1024

#define identical(fw,sw)	(strncmp(fw,sw,strlen(sw))==0 && strncmp(sw,fw,strlen(fw))==0)
#define ew(fw,sw)    (strcmp((fw),(sw)) == strcmp((sw),(fw)))

