#include <stdio.h>
#include <string.h>
#include <avs/avs.h>
#include <avs/port.h>
#include <avs/field.h>

#include <math.h>

#include "nabcode.h"

#define	MAT_ALLOC(n)	( ( MATRIX_T * )malloc( ( n )*sizeof( MATRIX_T ) ) )
#define	BPN		14
char	*MAT_getsyminfo();

#define	D_CONST		0
#define	D_UNIFORM	1
#define	D_GAUSS		2

#define	S_CONST		"Constant"
#define	S_CONST_P1	"value"
#define	S_CONST_P2	"ignore"
#define	S_UNIFORM	"Uniform"
#define	S_UNIFORM_P1	"lower"
#define	S_UNIFORM_P2	"upper"
#define	S_GAUSS		"Gaussian"
#define	S_GAUSS_P1	"mu"
#define	S_GAUSS_P2	"sigma"
 
static	s_xtdist = D_GAUSS;
static	s_xrdist = D_GAUSS;
static	s_ytdist = D_GAUSS;
static	s_yrdist = D_GAUSS;
static	s_ztdist = D_GAUSS;
static	s_zrdist = D_GAUSS;

static	void	get_dist_names();
static	int	check_dist_parms();

static	float	mk_rand_val();

static	int	r_init = 1;

/*long	random();*/
#define	MAXRAND	2147483647
#define	FRAND()	(1.*random()/MAXRAND)

static	float	g_rdist();

/* *****************************************/
/*  Module Description                     */
/* *****************************************/
int mg_rand_desc()
{

	int in_port, out_port, param, iresult;
	extern int mg_rand_compute();

	AVSset_module_name( "mg_rand", MODULE_DATA );

	/* Input Port Specifications              */
	in_port = AVScreate_input_port( "InMatrices", 
		"field 1D 1-space 1-vector uniform byte", OPTIONAL );

	/* Output Port Specifications              */
	out_port = AVScreate_output_port( "OutMatrices", 
		"field 1D 1-space 1-vector uniform byte" );

	/* Parameter Specifications                */
	param = AVSadd_parameter( "count", "integer", 1, INT_UNBOUND,
		INT_UNBOUND );
	AVSconnect_widget( param, "typein_integer" );

	param = AVSadd_parameter( "xtrn", "string", "X Translation", "", ":" );
	AVSconnect_widget( param, "text" );

	param = AVSadd_parameter( "xtdist", "tristate", 2, 0, 2 );
	AVSconnect_widget( param, "tristate" );
	AVSadd_parameter_prop( param, "title", "string", S_GAUSS );

	param = AVSadd_float_parameter( "xt1", 0.00000, FLOAT_UNBOUND,
		FLOAT_UNBOUND );
	AVSconnect_widget( param, "typein_real" );
	AVSadd_parameter_prop( param, "title", "string", S_GAUSS_P1 );

	param = AVSadd_float_parameter( "xt2", 1.00000, FLOAT_UNBOUND,
		FLOAT_UNBOUND );
	AVSconnect_widget( param, "typein_real" );
	AVSadd_parameter_prop( param, "title", "string", S_GAUSS_P2 );

	param = AVSadd_parameter( "xrot", "string", "X Rotation", "", ":" );
	AVSconnect_widget( param, "text" );
	AVSadd_parameter_prop( param, "layout", "string_block",
	"panel $Module -w panel -p \"Top Level Stack\" \n\
		manipulator $Module:xrot -w text -xy 128,32" );

	param = AVSadd_parameter( "xrdist", "tristate", 2, 0, 2 );
	AVSconnect_widget( param, "tristate" );
	AVSadd_parameter_prop( param, "title", "string", S_GAUSS );
	AVSadd_parameter_prop( param, "layout", "string_block",
	"panel $Module -w panel -p \"Top Level Stack\" \n\
		manipulator $Module:xrdist -w tristate -xy 128,54" );

	param = AVSadd_float_parameter( "xr1", 0.00000, FLOAT_UNBOUND,
		FLOAT_UNBOUND );
	AVSconnect_widget( param, "typein_real" );
	AVSadd_parameter_prop( param, "title", "string", S_GAUSS_P1 );
	AVSadd_parameter_prop( param, "layout", "string_block",
	"panel $Module -w panel -p \"Top Level Stack\" \n\
		manipulator $Module:xr1 -w typein_real -xy 128,76" );

	param = AVSadd_float_parameter( "xr2", 1.00000, FLOAT_UNBOUND,
		FLOAT_UNBOUND );
	AVSconnect_widget( param, "typein_real" );
	AVSadd_parameter_prop( param, "title", "string", S_GAUSS_P2 );
	AVSadd_parameter_prop( param, "layout", "string_block",
	"panel $Module -w panel -p \"Top Level Stack\" \n\
		manipulator $Module:xr2 -w typein_real -xy 128,98" );

	param = AVSadd_parameter( "ytrn", "string", "Y Translation", "", ":" );
	AVSconnect_widget( param, "text" );
	AVSadd_parameter_prop( param, "layout", "string_block",
	"panel $Module -w panel -p \"Top Level Stack\" \n\
		manipulator $Module:ytrn -w text -xy 10,120" );

	param = AVSadd_parameter( "ytdist", "tristate", 2, 0, 2 );
	AVSconnect_widget( param, "tristate" );
	AVSadd_parameter_prop( param, "title", "string", S_GAUSS );
	AVSadd_parameter_prop( param, "layout", "string_block",
	"panel $Module -w panel -p \"Top Level Stack\" \n\
		manipulator $Module:ytdist -w tristate -xy 10,142" );

	param = AVSadd_float_parameter( "yt1", 0.00000, FLOAT_UNBOUND,
		FLOAT_UNBOUND );
	AVSconnect_widget( param, "typein_real" );
	AVSadd_parameter_prop( param, "title", "string", S_GAUSS_P1 );
	AVSadd_parameter_prop( param, "layout", "string_block",
	"panel $Module -w panel -p \"Top Level Stack\" \n\
		manipulator $Module:yt1 -w typein_real -xy 10,164" );

	param = AVSadd_float_parameter( "yt2", 1.00000, FLOAT_UNBOUND,
		FLOAT_UNBOUND );
	AVSconnect_widget( param, "typein_real" );
	AVSadd_parameter_prop( param, "title", "string", S_GAUSS_P2 );
	AVSadd_parameter_prop( param, "layout", "string_block",
	"panel $Module -w panel -p \"Top Level Stack\" \n\
		manipulator $Module:yt2 -w typein_real -xy 10,186" );

	param = AVSadd_parameter( "yrot", "string", "Y Rotation", "", ":" );
	AVSconnect_widget( param, "text" );
	AVSadd_parameter_prop( param, "layout", "string_block",
	"panel $Module -w panel -p \"Top Level Stack\" \n\
		manipulator $Module:yrot -w text -xy 128,120" );

	param = AVSadd_parameter( "yrdist", "tristate", 2, 0, 2 );
	AVSconnect_widget( param, "tristate" );
	AVSadd_parameter_prop( param, "title", "string", S_GAUSS );
	AVSadd_parameter_prop( param, "layout", "string_block",
	"panel $Module -w panel -p \"Top Level Stack\" \n\
		manipulator $Module:yrdist -w tristate -xy 128,142" );

	param = AVSadd_float_parameter( "yr1", 0.00000, FLOAT_UNBOUND,
		FLOAT_UNBOUND );
	AVSconnect_widget( param, "typein_real" );
	AVSadd_parameter_prop( param, "title", "string", S_GAUSS_P1 );
	AVSadd_parameter_prop( param, "layout", "string_block",
	"panel $Module -w panel -p \"Top Level Stack\" \n\
		manipulator $Module:yr1 -w typein_real -xy 128,164" );

	param = AVSadd_float_parameter( "yr2", 1.00000, FLOAT_UNBOUND,
		FLOAT_UNBOUND );
	AVSconnect_widget( param, "typein_real" );
	AVSadd_parameter_prop( param, "title", "string", S_GAUSS_P2 );
	AVSadd_parameter_prop( param, "layout", "string_block",
	"panel $Module -w panel -p \"Top Level Stack\" \n\
		manipulator $Module:yr2 -w typein_real -xy 128,186" );

	param = AVSadd_parameter( "ztrn", "string", "Z Translation", "", ":" );
	AVSconnect_widget( param, "text" );
	AVSadd_parameter_prop( param, "layout", "string_block",
	"panel $Module -w panel -p \"Top Level Stack\" \n\
		manipulator $Module:ztrn -w text -xy 10,208" );

	param = AVSadd_parameter( "ztdist", "tristate", 2, 0, 2 );
	AVSconnect_widget( param, "tristate" );
	AVSadd_parameter_prop( param, "title", "string", S_GAUSS );
	AVSadd_parameter_prop( param, "layout", "string_block",
	"panel $Module -w panel -p \"Top Level Stack\" \n\
		manipulator $Module:ztdist -w tristate -xy 10,230" );

	param = AVSadd_float_parameter( "zt1", 0.00000, FLOAT_UNBOUND,
		FLOAT_UNBOUND );
	AVSconnect_widget( param, "typein_real" );
	AVSadd_parameter_prop( param, "title", "string", S_GAUSS_P1 );
	AVSadd_parameter_prop( param, "layout", "string_block",
	"panel $Module -w panel -p \"Top Level Stack\" \n\
		manipulator $Module:zt1 -w typein_real -xy 10,252" );

	param = AVSadd_float_parameter( "zt2", 1.00000, FLOAT_UNBOUND,
		FLOAT_UNBOUND );
	AVSconnect_widget( param, "typein_real" );
	AVSadd_parameter_prop( param, "title", "string", S_GAUSS_P2 );
	AVSadd_parameter_prop( param, "layout", "string_block",
	"panel $Module -w panel -p \"Top Level Stack\" \n\
		manipulator $Module:zt2 -w typein_real -xy 10,274" );

	param = AVSadd_parameter( "zrot", "string", "Z Rotation", "", ":" );
	AVSconnect_widget( param, "text" );
	AVSadd_parameter_prop( param, "layout", "string_block",
	"panel $Module -w panel -p \"Top Level Stack\" \n\
		manipulator $Module:zrot -w text -xy 128,208" );

	param = AVSadd_parameter( "zrdist", "tristate", 2, 0, 2 );
	AVSconnect_widget( param, "tristate" );
	AVSadd_parameter_prop( param, "title", "string", S_GAUSS );
	AVSadd_parameter_prop( param, "layout", "string_block",
	"panel $Module -w panel -p \"Top Level Stack\" \n\
		manipulator $Module:zrdist -w tristate -xy 128,230" );

	param = AVSadd_float_parameter( "zr1", 0.00000, FLOAT_UNBOUND,
		FLOAT_UNBOUND );
	AVSconnect_widget( param, "typein_real" );
	AVSadd_parameter_prop( param, "title", "string", S_GAUSS_P1 );
	AVSadd_parameter_prop( param, "layout", "string_block",
	"panel $Module -w panel -p \"Top Level Stack\" \n\
		manipulator $Module:zr1 -w typein_real -xy 128,252" );

	param = AVSadd_float_parameter( "zr2", 1.00000, FLOAT_UNBOUND,
		FLOAT_UNBOUND );
	AVSconnect_widget( param, "typein_real" );
	AVSadd_parameter_prop( param, "title", "string", S_GAUSS_P2 );
	AVSadd_parameter_prop( param, "layout", "string_block",
	"panel $Module -w panel -p \"Top Level Stack\" \n\
		manipulator $Module:zr2 -w typein_real -xy 128,274" );

	AVSset_compute_proc( mg_rand_compute );
	return( 1 );
}
 
/* *****************************************/
/* Module Compute Routine                  */
/* *****************************************/
int mg_rand_compute( InMatrices, OutMatrices, count,
	xtrn, xtdist, xt1, xt2,
	xrot, xrdist, xr1, xr2,
	ytrn, ytdist, yt1, yt2,
	yrot, yrdist, yr1, yr2,
	ztrn, ztdist, zt1, zt2,
	zrot, zrdist, zr1, zr2 )
AVSfield_char *InMatrices;
AVSfield_char **OutMatrices;
int count;
char *xtrn;
int xtdist;
float *xt1;
float *xt2;
char *xrot;
int xrdist;
float *xr1;
float *xr2;
char *ytrn;
int ytdist;
float *yt1;
float *yt2;
char *yrot;
int yrdist;
float *yr1;
float *yr2;
char *ztrn;
int ztdist;
float *zt1;
float *zt2;
char *zrot;
int zrdist;
float *zr1;
float *zr2;
{
	int	create;
	int	i, j, k;
	POINT_T pts[ 4 ], axes[ 4 ];
	REAL_T	angs[ 3 ];
	int	m, ms, mi, mo;
	int	nmats, simats, nimats, nomats;
	MATRIX_T	*mats, *imats, *omats;
	MATRIX_T	mat;
	int	dims0[ 1 ];
	char	*sip;
	char	m_si[ 1000 ];
	int	s_sip, s_o_sip;
	char	*dp;
	char	*s_dist, *s_dp1, *s_dp2;

	if( *OutMatrices )
		AVSfield_free( *OutMatrices );

	if( r_init ){
		srandom(0);
		r_init = 0;
	}

	create = 1;
	angs[ 0 ] = angs[ 1 ] = angs[ 2 ] = 0.0;

	nimats = nomats = 0;
	imats = omats = NULL;

	for( i = 0; i < 4; i++ ){
		for( j = 0; j < 3; j++ ){
			pts[ i ][ j ] = 0.0;
			axes[ i ][ j ] = 0.0;
		}
	}

	axes[ 1 ][ 0 ] = 1.0;
	axes[ 2 ][ 1 ] = 1.0;
	axes[ 3 ][ 2 ] = 1.0;

	if( count <= 0 ){
		AVSerror( "Count must > 0." );
		return( 0 );
	}

	if( xtdist != s_xtdist ){
		s_xtdist = xtdist;
		get_dist_names( xtdist, &s_dist, &s_dp1, &s_dp2 );
		AVSmodify_parameter_prop( "xtdist", "title", "string", s_dist );
		AVSmodify_parameter_prop( "xt1", "title", "string", s_dp1 );
		AVSmodify_parameter_prop( "xt2", "title", "string", s_dp2 );
	}
	if( !check_dist_parms("X Translation", xtdist, xt1, "xt1", xt2, "xt2") )
		return( 0 );
	if( xrdist != s_xrdist ){
		s_xrdist = xrdist;
		get_dist_names( xrdist, &s_dist, &s_dp1, &s_dp2 );
		AVSmodify_parameter_prop( "xrdist", "title", "string", s_dist );
		AVSmodify_parameter_prop( "xr1", "title", "string", s_dp1 );
		AVSmodify_parameter_prop( "xr2", "title", "string", s_dp2 );
	}
	if( !check_dist_parms( "X Rotation", xrdist, xr1, "xr1", xr2, "xr2" ) )
		return( 0 );
	if( ytdist != s_ytdist ){
		s_ytdist = ytdist;
		get_dist_names( ytdist, &s_dist, &s_dp1, &s_dp2 );
		AVSmodify_parameter_prop( "ytdist", "title", "string", s_dist );
		AVSmodify_parameter_prop( "yt1", "title", "string", s_dp1 );
		AVSmodify_parameter_prop( "yt2", "title", "string", s_dp2 );
	}
	if( !check_dist_parms("Y Translation", ytdist, yt1, "yt1", yt2, "yt2") )
		return( 0 );
	if( yrdist != s_yrdist ){
		s_yrdist = yrdist;
		get_dist_names( yrdist, &s_dist, &s_dp1, &s_dp2 );
		AVSmodify_parameter_prop( "yrdist", "title", "string", s_dist );
		AVSmodify_parameter_prop( "yr1", "title", "string", s_dp1 );
		AVSmodify_parameter_prop( "yr2", "title", "string", s_dp2 );
	}
	if( !check_dist_parms( "Y Rotation", yrdist, yr1, "yr1", yr2, "yr2" ) )
		return( 0 );
	if( ztdist != s_ztdist ){
		s_ztdist = ztdist;
		get_dist_names( ztdist, &s_dist, &s_dp1, &s_dp2 );
		AVSmodify_parameter_prop( "ztdist", "title", "string", s_dist );
		AVSmodify_parameter_prop( "zt1", "title", "string", s_dp1 );
		AVSmodify_parameter_prop( "zt2", "title", "string", s_dp2 );
	}
	if( !check_dist_parms("Z Translation", ztdist, zt1, "zt1", zt2, "zt2") )
		return( 0 );
	if( zrdist != s_zrdist ){
		s_zrdist = zrdist;
		get_dist_names( zrdist, &s_dist, &s_dp1, &s_dp2 );
		AVSmodify_parameter_prop( "zr1", "title", "string", s_dp1 );
		AVSmodify_parameter_prop( "zr2", "title", "string", s_dp2 );
		AVSmodify_parameter_prop( "zrdist", "title", "string", s_dist );
	}
	if( !check_dist_parms( "Z Rotation", zrdist, zr1, "zr1", zr2, "zr2" ) )
		return( 0 );

	nmats = count;
	mats = MAT_ALLOC( nmats );
	if( mats == NULL ){
		AVSerror( "Allocation of mats failed." );
		return( 0 );
	}

	for( i = 0; i < nmats; i++ ){
		angs[0] = mk_rand_val( xrdist, xr1, xr2 );
		angs[1] = mk_rand_val( yrdist, yr1, yr2 );
		angs[2] = mk_rand_val( zrdist, zr1, zr2 );
		MAT_orient( axes, angs, mats[i] );
		mats[i][3][0] = mk_rand_val( xtdist, xt1, xt2 );
		mats[i][3][1] = mk_rand_val( ytdist, yt1, yt2 );
		mats[i][3][2] = mk_rand_val( ztdist, zt1, zt2 );
	}

	if( !InMatrices ){
		nimats = 0;
		sip = NULL;
		s_sip = 0;
	}else{
		simats = MAT_count( InMatrices->data );
		imats = MAT_ALLOC( simats );
		if( imats == NULL ){
			AVSerror( "Allocation of imats failed." );
			free( mats );
			return( 0 );
		}
		if( ( nimats = MAT_sscan(InMatrices->data,simats,imats))==0 ){
			AVSerror( "InMatrices is not a valid matrix field." );
			free( imats );
			free( mats );
			return( 0 );
		}
		sip = MAT_getsyminfo( InMatrices->data );
		s_sip = strlen( sip );
		create = 0;
	}

	if( create || nimats == 0 ){
		nomats = nmats;
		omats = ( MATRIX_T * )mats;
	}else{
		nomats = nmats;
		omats = MAT_ALLOC( nomats );
		if( omats == NULL ){
			AVSerror( "Allocation of omats failed." );
			free( imats );
			free( mats );
			return( 0 );
		}
		mo = 0;
		for( m = 0; m < nmats; m++, mo++ ){
			NAB_matcpy( omats[ mo ],
				MAT_concat( imats[ mi ], mats[ m ] ) );
		}
	}

	dims0[ 0 ] = 16 * BPN * nomats;
	*OutMatrices = ( AVSfield_char * )AVSdata_alloc(
		"field 1D 1-space 1-vector uniform byte", dims0 );
	if( *OutMatrices == NULL ){
		AVSerror( "Allocation of output field failed." );
		if( !create && nimats != 0 )
			free( omats );
		if( imats )
			free( imats );
		free( mats );
		return( 0 );
	}

	dp = ( *OutMatrices )->data;

	MAT_sprint( dp, nomats, omats );

	if( !create && nimats != 0 )
		free( omats );
	if( imats )
		free( imats );
	free( mats );

	return(1);
}
 
/* ***********************************************************************/
/* Initialization for modules contained in this file.                    */
/* ***********************************************************************/
static int ((*mod_list[])()) = {
	mg_rand_desc
};
#define NMODS (sizeof(mod_list) / sizeof(char *))

AVSinit_modules()
{
	AVSinit_from_module_list(mod_list, NMODS);
}

static	void	get_dist_names( dist_type, d_name, p1_name, p2_name )
int	dist_type;
char	**d_name;
char	**p1_name;
char	**p2_name;
{

	switch( dist_type ){
	case D_CONST :
		*d_name = S_CONST;
		*p1_name = S_CONST_P1;
		*p2_name = S_CONST_P2;
		break;
	case D_UNIFORM :
		*d_name = S_UNIFORM;
		*p1_name = S_UNIFORM_P1;
		*p2_name = S_UNIFORM_P2;
		break;
	case D_GAUSS :
		*d_name = S_GAUSS;
		*p1_name = S_GAUSS_P1;
		*p2_name = S_GAUSS_P2;
		break;
	}
}

static	int	check_dist_parms( pname, dist_type, dp1, pn1, dp2, pn2 )
char	*pname;
int	dist_type;
float	*dp1;
char	*pn1;
float	*dp2;
char	*pn2;
{
	char	e_msg[ 256 ];

	switch( dist_type ){
	case D_CONST :
		break;
	case D_UNIFORM :
		if( *dp1 > *dp2 ){
			sprintf( e_msg, "Uniform %s: %s must be <= %s.",
				pname, pn1, pn2 );
			AVSerror( e_msg );
			return( 0 );
		}
		break;
	case D_GAUSS :
		if( *dp2 < 0.0 ){
			sprintf( e_msg, "Gaussian %s: %s must be >= 0.",
				pname, pn2 );
			AVSerror( e_msg );
			return( 0 );
		}
		break;
	}
	return( 1 );
}

static	float	mk_rand_val( dist_type, dp1, dp2 )
int	dist_type;
float	*dp1;
float	*dp2;
{

	switch( dist_type ){
	case D_CONST :
		return( *dp1 );
		break;
	case D_UNIFORM :
		if( *dp1 == *dp2 )
			return( *dp1 );
		else
			return( ( *dp2 - *dp1 ) * random() / MAXRAND + *dp1 );
			return( ( *dp2 - *dp1 ) * FRAND() + *dp1 );
		break;
	case D_GAUSS :
		if( *dp2 == 0.0 )
			return( *dp1 );
		else
			return( g_rdist( dp1, dp2 ) );
		break;
	}
}

static	float	g_rdist( mu, sigma ) 
float	*mu;
float	*sigma;
{
	static	int	iset = 0;
	static	float	gset;
	float	fac, rsq, v1, v2, r;

	if( !iset ){
		do{
			v1 = 2.0 * FRAND() - 1.0;
			v2 = 2.0 * FRAND() - 1.0;
			rsq = v1 * v1 + v2 * v2;
		}while( rsq >= 1.0 || rsq == 0.0 );
		fac = sqrt( -2.0 * log( rsq ) / rsq ); 
		gset = v1 * fac;
		iset = 1;
		r = v2 * fac;
	}else{
		iset = 0;
		r = gset;
	}
	return( r * *sigma + *mu );
}
