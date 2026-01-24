#include <stdio.h>
#include <string.h>
#include <avs/avs.h>
#include <avs/port.h>
#include <avs/field.h>

#include "nabcode.h"

#define	MAT_ALLOC(n)	( ( MATRIX_T * )malloc( ( n )*sizeof( MATRIX_T ) ) )
#define	WHITE	" \t\n"
#define	BPN	14
typedef	char	NAME_T[ 32 ];
static	void	mk_setmatcmd();
static	void	mk_hidecmd();

#define	OV_ZERO		0
#define	SV_ZERO		"0"
#define	OV_ONE		1
#define	SV_ONE		"1"
#define	OV_LAST		2
#define	SV_LAST		"Last Value"

static	int	xform_obj();
 
/* *****************************************/
/*  Module Description                     */
/* *****************************************/
int ma_desc()
{

	int in_port, out_port, param, iresult;
	extern int ma_compute();

	AVSset_module_name("ma", MODULE_RENDER);

	/* Input Port Specifications               */
	in_port = AVScreate_input_port("Matrices", 
		"field 1D 1-space 1-vector uniform byte", REQUIRED);
	in_port = AVScreate_input_port("ObjNames", 
		"field 1D 1-space uniform byte", REQUIRED);
	in_port = AVScreate_input_port("Occupancy", 
		"field 1D 1-space 1-vector uniform integer", OPTIONAL);

	param = AVSadd_parameter( "hide", "boolean", 0, 0, 1 );
	AVSadd_parameter_prop( param, "title", "string",
		"Hide Untransformed Objects" );
	AVSadd_parameter_prop( param, "width", "integer", 4 );

	param = AVSadd_parameter( "ig_occ", "boolean", 0, 0, 1 );
	AVSadd_parameter_prop( param, "title", "string",
		"Ignore Occupancy Field" );
	AVSadd_parameter_prop( param, "width", "integer", 4 );

	param = AVSadd_parameter( "t_ovals", "string",
		"Set missing Occupancy values to", "", "" );
	AVSadd_parameter_prop( param, "width", "integer", 4 );
	AVSconnect_widget( param, "text" );

	param = AVSadd_parameter( "m_ovals", "choice",
		"False", "False:True:Last Value", ":" );
	AVSconnect_widget( param, "radio_buttons" );
	AVSadd_parameter_prop( param, "width", "integer", 4 );

	AVSset_compute_proc(ma_compute);
	return(1);
}
 
/* *****************************************/
/* Module Compute Routine                  */
/* *****************************************/
int ma_compute( Matrices, ObjNames, Occupancy, hide, ig_occ, t_ovals, m_ovals )
AVSfield_char *Matrices;
AVSfield_char *ObjNames;
AVSfield_int  *Occupancy;
int	hide;
int	ig_occ;	
char	*t_ovals;
char	*m_ovals;
{
	int	m, nmats;
	MATRIX_T	*mats;
	int	n, nw, nnames;
	char	*np, *enp;
	NAME_T	*names;
	char	cmdbuf[ 1024], *outbuf, *errbuf;

	nmats = MAT_count( Matrices->data );
	if( nmats == 0 )
		return( 0 );

	mats = MAT_ALLOC( nmats );
	if( mats == NULL ){
		AVSerror( "Allocation of mats failed." );
		return( 0 );
	}

	for( nnames = 0, np = ObjNames->data; *np; nnames++ ){
		np += strspn( np, WHITE );
		enp = strpbrk( np, WHITE );
		np = enp + 1;
	}
	if( nnames == 0 ){
		free( mats );
		return( 0 );
	}
	names = ( NAME_T * )malloc( nnames * sizeof( NAME_T ) ); 
	if( names == NULL ){
		AVSerror( "Allocation of names failed." );
		free( mats );
		return( 0 );
	}
	for( n = 0, np = ObjNames->data; *np; n++ ){
		np += strspn( np, WHITE );
		enp = strpbrk( np, WHITE );
		strncpy( names[ n ], np, enp - np );
		names[ n ][ enp - np ] = '\0';
		np = enp + 1;
	}

	MAT_sscan( Matrices->data, nmats, mats );

	if( nmats <= nnames ){
		for( m = 0; m < nmats; m++ ){
			if( xform_obj( m, ig_occ, m_ovals, Occupancy ) ){
				mk_setmatcmd( names[ m ], mats[ m ], cmdbuf );
				AVScommand( "kernel", cmdbuf,
					&outbuf, &errbuf );
			}else if( hide ){
				sprintf( cmdbuf,
					"geom_set_visibility -object %s 0",
					names[ m ] );
				AVScommand( "kernel", cmdbuf,
					&outbuf, &errbuf );
			}else{
				sprintf( cmdbuf,
					"geom_set_visibility -object %s 1",
					names[ m ] );
				AVScommand( "kernel", cmdbuf,
					&outbuf, &errbuf );
			}
		}
		if( hide ){
			for( m = nmats; m < nnames; m++ ){
				sprintf( cmdbuf,
					"geom_set_visibility -object %s 0",
					names[ m ] );
				AVScommand( "kernel", cmdbuf,
					&outbuf, &errbuf );
			}
		}else{
			for( m = nmats; m < nnames; m++ ){
				sprintf( cmdbuf,
					"geom_set_visibility -object %s 1",
					names[ m ] );
				AVScommand( "kernel", cmdbuf,
					&outbuf, &errbuf );
			}
		}
	}else{
		for( n = 0; n < nnames; n++ ){
			if( xform_obj( n, ig_occ, m_ovals, Occupancy ) ){
				mk_setmatcmd( names[ n ], mats[ n ], cmdbuf );
				AVScommand( "kernel", cmdbuf,
					&outbuf, &errbuf );
			}else if( hide ){
				sprintf( cmdbuf,
					"geom_set_visibility -object %s 0",
					names[ n ] );
				AVScommand( "kernel", cmdbuf,
					&outbuf, &errbuf );
			}else{
				sprintf( cmdbuf,
					"geom_set_visibility -object %s 1",
					names[ n ] );
				AVScommand( "kernel", cmdbuf,
					&outbuf, &errbuf );
			}
		}
	}
	AVScommand( "kernel", "geom_refresh", &outbuf, &errbuf );

	free( mats );
	free( names );
	return(1);
}

static	void	mk_setmatcmd( name, mat, cmd )
char		name[];
MATRIX_T	mat;
char		cmd[];
{

	sprintf( cmd, "geom_set_matrix -object %s \
		-mat %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
		name,
		mat[ 0 ][ 0 ], mat[ 0 ][ 1 ], mat[ 0 ][ 2 ], mat[ 0 ][ 3 ],
		mat[ 1 ][ 0 ], mat[ 1 ][ 1 ], mat[ 1 ][ 2 ], mat[ 1 ][ 3 ],
		mat[ 2 ][ 0 ], mat[ 2 ][ 1 ], mat[ 2 ][ 2 ], mat[ 2 ][ 3 ],
		mat[ 3 ][ 0 ], mat[ 3 ][ 1 ], mat[ 3 ][ 2 ], mat[ 3 ][ 3 ] );
}

static	int	xform_obj( n, ig_occ, m_oval, occ )
int	n;
int	ig_occ;
char	*m_oval;
AVSfield_int	*occ;
{
	int	n_occ;
	int	val;

	if( occ == NULL || ig_occ )
		return( 1 );
	
	n_occ = occ->dimensions[ 0 ];
	if( n < n_occ ){
		return( occ->data[ n ] );
	}else if( !strcmp( m_oval, "0" ) )
		return( 0 );
	else if( !strcmp( m_oval, "1" ) )
		return( 1 );
	else
		return( occ->data[ n_occ - 1 ] );
}

/* ***********************************************************************/
/* Initialization for modules contained in this file.                    */
/* ***********************************************************************/
static int ((*mod_list[])()) = {
	ma_desc
};
#define NMODS (sizeof(mod_list) / sizeof(char *))

AVSinit_modules()
{
	AVSinit_from_module_list(mod_list, NMODS);
}
