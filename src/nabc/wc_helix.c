// wc_helix() - create Watson/Crick duplex

#include "nabc.h"

//   see user's manual for explanation of the code.

MOLECULE_T *wc_helix(
	char *seq,  char *snatype, char *aseq, char *anatype,
	double xoff, double incl, double twist, double rise,
	char *opts )
{
	MOLECULE_T  *m1, *m2, *m3;
	MATRIX_T xomat, inmat, mat;
	char *arname, *srname;
	RESIDUE_T *sres, *ares;
	int	has_s, has_a;
	int i, slen;
	double	ttwist, trise;

	has_s = 1; has_a = 1;

	if( seq == NULL && aseq == NULL ){
		fprintf( stderr, "wc_helix: no sequence\n" );
		return( NULL );
	}else if( seq == NULL ){
		seq = wc_complement( aseq, snatype );
		has_s = 0;
	}else if( aseq == NULL ){
		aseq = wc_complement( seq, anatype );
		has_a = 0;
	}

	slen = strlen( seq );

    i = 1;
	setreslibkind( snatype );
	if( !strcmp(snatype, "dna"))  
        asprintf( &srname, "D%c", toupper(seq[i-1]) );
    else
        asprintf( &srname, "%c", toupper(seq[i-1]) );
	if( strstr( opts, "s5" ) )
        asprintf( &srname, "%s%c", srname, '5' );
	else if( strstr( opts, "s3" ) && slen == 1 )
        asprintf( &srname, "%s%c", srname, '3' );
	sres = getresidue( srname );

	setreslibkind( anatype );
	if( !strcmp(anatype, "dna"))  
        asprintf( &arname, "D%c", toupper(aseq[i-1]) );
    else
        asprintf( &arname, "%c", toupper(aseq[i-1]) );
	if( strstr( opts, "a3" ) )
        asprintf( &arname, "%s%c", arname, '3' );
	else if( strstr( opts, "a5" ) && slen == 1 )
        asprintf( &arname, "%s%c", arname, '5' );
	ares = getresidue( arname );

	m1 = wc_basepair( sres, ares );
	freeresidue( sres );
	freeresidue( ares );

	copy_mat( newtransform(xoff, 0., 0., 0., 0., 0. ), xomat );
	transformmol( xomat, m1, NULL );
	copy_mat( newtransform( 0., 0., 0., incl, 0., 0.), inmat );
	transformmol( inmat, m1, NULL );

	trise = rise; ttwist = twist;
	for( i=2; i<=slen-1; i++ ){

        setreslibkind( snatype );
        if( !strcmp(snatype, "dna"))  
            asprintf( &srname, "D%c", toupper(seq[i-1]) );
        else
            asprintf( &srname, "%c", toupper(seq[i-1]) );
        sres = getresidue( srname );

        setreslibkind( anatype );
        if( !strcmp(anatype, "dna"))  
            asprintf( &arname, "D%c", toupper(aseq[i-1]) );
        else
            asprintf( &arname, "%c", toupper(aseq[i-1]) );
        ares = getresidue( arname );

		m2 = wc_basepair( sres, ares );
		freeresidue( sres );
		freeresidue( ares );
		transformmol( xomat, m2, NULL );
		transformmol( inmat, m2, NULL );
		copy_mat( newtransform( 0., 0., trise, 0., 0., ttwist ), mat );
		transformmol( mat, m2, NULL );
		mergestr( m1, "sense", "last", m2, "sense", "first" );
		connectres( m1, "sense", i-1, "O3'", i, "P" );
		mergestr( m1, "anti", "first", m2, "anti", "last" );
		connectres( m1, "anti", 1, "O3'", 2, "P" );
		trise = trise + rise;
		ttwist = ttwist + twist;
		freemolecule( m2 );
	}

	i = slen;       // add in final residue pair
	if( i > 1 ){

        setreslibkind( snatype );
        if( !strcmp(snatype, "dna"))  
            asprintf( &srname, "D%c", toupper(seq[i-1]) );
        else
            asprintf( &srname, "%c", toupper(seq[i-1]) );
	    if( strstr( opts, "s3" ) )
            asprintf( &srname, "%s%c", srname, '3' );
        sres = getresidue( srname );

        setreslibkind( anatype );
        if( !strcmp(anatype, "dna"))  
            asprintf( &arname, "D%c", toupper(aseq[i-1]) );
        else
            asprintf( &arname, "%c", toupper(aseq[i-1]) );
	    if( strstr( opts, "a5" ) )
            asprintf( &arname, "%s%c", arname, '5' );
        ares = getresidue( arname );

		m2 = wc_basepair( sres, ares );
		freeresidue( sres );
		freeresidue( ares );
		transformmol( xomat, m2, NULL );
		transformmol( inmat, m2, NULL );
		copy_mat( newtransform( 0., 0., trise, 0., 0., ttwist ), mat );
		transformmol( mat, m2, NULL );
		mergestr( m1, "sense", "last", m2, "sense", "first" );
		connectres( m1, "sense", i-1, "O3'", i, "P" );
		mergestr( m1, "anti", "first", m2, "anti", "last" );
		connectres( m1, "anti", 1, "O3'", 2, "P" );
		trise = trise + rise;
		ttwist = ttwist + twist;
		freemolecule( m2 );
	}

	m3 = newmolecule();
	addstrand( m3, "sense" );
	addstrand( m3, "anti" );
	if( has_s )
		mergestr( m3, "sense", "last", m1, "sense", "first" );
	if( has_a )
		mergestr( m3, "anti", "last", m1, "anti", "first" );

	freemolecule( m1 );

	return( m3 );
};

