// wc_complement() - create W/C string to seq.

#include "nabc.h"

char*	wc_complement( char* seq, char* rlt )
{
	char	acbase, base;
    char*   wcseq=strdup("");
    char    wcbase[2];
	int	i, len;

	if( !strcmp(rlt,"dna") )
		acbase = 't';
	else if( !strcmp(rlt,"rna") )
		acbase = 'u';
	else{
		fprintf( stderr,
		"wc_complement: rlt (%s) is not dna/rna, has no W/C complement\n",
			rlt );
		exit( 1 );
	}
	
	len = strlen( seq );
    wcbase[1] = '\0';
	for( i=0; i<len; i++ ){
		base = seq[i];
		if( base == 'a' || base == 'A' )      wcbase[0] = acbase;
		else if( base == 'c' || base == 'C' ) wcbase[0] = 'g';
		else if( base == 'g' || base == 'G' ) wcbase[0] = 'c';
		else if( base == 't' || base == 'T' ) wcbase[0] = 'a';
		else if( base == 'u' || base == 'U' ) wcbase[0] = 'a';
		else{
			fprintf( stderr, "wc_complement: unknown base %c\n",
				base );
			exit( 1 );
		}
		// wcseq = wcseq + wcbase; ignore memory leaks here
        asprintf(&wcseq, "%s%s", wcseq, wcbase );
	}
	return( wcseq );
}

#if 0
//  simple driver to test
int main( int argc, char* argv[] )
{
    printf( "wc_complement returns %s\n", 
       wc_complement( argv[1], argv[2] ) );
}
#endif
