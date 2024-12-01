// wc_basepair() - create Watson/Crick base pair

#include "nabc.h"

#define AT_SEP 8.29
#define CG_SEP 8.27
MOLECULE_T *wc_basepair( RESIDUE_T *sres, RESIDUE_T *ares )
{
	MOLECULE_T *m, *m_sense, *m_anti;
	double	sep;
	char *arname, *srname;
	char *xtail, *xhead;
	char *ytail, *yhead;
	MATRIX_T	mat;

	m = newmolecule();
	m_sense = newmolecule();
	m_anti = newmolecule();
	addstrand( m, "sense" );
	addstrand( m, "anti" );
	addstrand( m_sense, "sense" );
	addstrand( m_anti, "anti" );

    srname = strdup( getresname(  sres ) );
    arname = strdup( getresname(  ares ) );
    ytail = strdup( "sense::C1'" );
    yhead = strdup( "anti::C1'" );
    if( strchr(srname,'A') ){
		sep = AT_SEP;
		xtail = strdup( "sense::C5" );
		xhead = strdup( "sense::N3" );
		addresidue( m_sense, "sense", sres );
		setframe( 2, m_sense, "::C4", "::C5", "::N3", "::C4", "::N1" );
    }else if( strchr(srname,'C') ){
		sep = CG_SEP;
		xtail = strdup( "sense::C6" );
		xhead = strdup( "sense::N1" );
		addresidue( m_sense, "sense", sres );
		setframe( 2, m_sense, "::C6", "::C5", "::N1", "::C6", "::N3" );
    }else if( strchr(srname,'G') ){
		sep = CG_SEP;
		xtail = strdup( "sense::C5" );
		xhead = strdup( "sense::N3" );
		addresidue( m_sense, "sense", sres );
		setframe( 2, m_sense, "::C4", "::C5", "::N3", "::C4", "::N1" );
    }else if( strchr(srname,'T') ){
		sep = AT_SEP;
		xtail = strdup( "sense::C6" );
		xhead = strdup( "sense::N1" );
		addresidue( m_sense, "sense", sres );
		setframe( 2, m_sense, "::C6", "::C5", "::N1", "::C6", "::N3" );
    }else if( strchr(srname,'U') ){
		sep = AT_SEP;
		xtail = strdup( "sense::C6" );
		xhead = strdup( "sense::N1" );
		addresidue( m_sense, "sense", sres );
		setframe( 2, m_sense, "::C6", "::C5", "::N1", "::C6", "::N3" );
	}else{
		fprintf( stderr,"wc_basepair : unknown sres %s\n",srname );
 		exit( 1 );
	}
    if( strchr(arname,'A') ){
		addresidue( m_anti, "anti", ares );
		setframe( 2, m_anti, "::C4", "::C5", "::N3", "::C4", "::N1" );
    }else if( strchr(arname,'C') ){
		addresidue( m_anti, "anti", ares );
		setframe( 2, m_anti, "::C6", "::C5", "::N1", "::C6", "::N3" );
    }else if( strchr(arname,'G') ){
		addresidue( m_anti, "anti", ares );
		setframe( 2, m_anti, "::C4", "::C5", "::N3", "::C4", "::N1" );
    }else if( strchr(arname,'T') ){
		addresidue( m_anti, "anti", ares );
		setframe( 2, m_anti, "::C6", "::C5", "::N1", "::C6", "::N3" );
    }else if( strchr(arname,'U') ){
		addresidue( m_anti, "anti", ares );
		setframe( 2, m_anti, "::C6", "::C5", "::N1", "::C6", "::N3" );
	}else{
		fprintf( stderr,"wc_basepair : unknown ares %s\n",arname );
 		exit( 1 );
	}
	
	alignframe( m_sense, NULL );
	alignframe( m_anti, NULL );
	NAB_matcpy( mat, newtransform( 0., 0., 0., 180., 0., 0. ) );
	transformmol( mat, m_anti, NULL );
	NAB_matcpy( mat, newtransform( 0., sep, 0., 0., 0., 0. ) );
	transformmol( mat, m_anti, NULL );
	mergestr( m, "sense", "last", m_sense, "sense", "first" );
	mergestr( m, "anti", "last", m_anti, "anti", "first" );

	freemolecule( m_sense );
	freemolecule( m_anti );

	setframe( 2, m, "::C1'", xtail, xhead, ytail, yhead );
	alignframe( m, NULL );

	return( m );
}
