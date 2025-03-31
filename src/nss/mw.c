#include <stdio.h>
#include <avs/avs.h>
#include <avs/port.h>
#include <avs/field.h>

#include "nabcode.h"

#define	MAT_ALLOC(n)	( ( MATRIX_T * )malloc( ( n )*sizeof( MATRIX_T ) ) )
#define	BPN	14
 
char	*MAT_getsyminfo();

/* *****************************************/
/*  Module Description                     */
/* *****************************************/
int mw_desc()
{

	int in_port, out_port, param, iresult;
	extern int mw_compute();

	AVSset_module_name( "mw", MODULE_RENDER );

	/* Input Port Specifications               */
	in_port = AVScreate_input_port( "InMatrices", 
		"field 1D 1-space 1-vector uniform byte", REQUIRED );

	/* Parameter Specifications                */
	param = AVSadd_parameter( "MatFile", "string", "", "", ":" );
	AVSconnect_widget( param, "browser" );

	AVSset_compute_proc( mw_compute );
	return( 1 );
}
 
/* *****************************************/
/* Module Compute Routine                  */
/* *****************************************/
int mw_compute( InMatrices, MatFile )
	AVSfield_char *InMatrices;
	char *MatFile;
{
	FILE	*fp;
	char	e_msg[ 256 ];
	int	nmats, smats;
	MATRIX_T *mats;
	char	*sip;

	if( !*MatFile )
		return( 0 );

	if( ( fp = fopen( MatFile, "w" ) ) == NULL ){
		sprintf( e_msg, "Can't open MatFile %s.", MatFile );
		AVSerror( e_msg );
		return( 0 );
	}

	smats = MAT_count( InMatrices->data );
	mats = MAT_ALLOC( smats );
	if( mats == NULL ){
		AVSerror( "Allocation of mats failed." );
		return( 0 );
	}

	sip = MAT_getsyminfo( InMatrices->data );
	fputs( sip, fp );
	if( ( nmats = MAT_sscan( InMatrices->data, smats, mats ) ) == 0 ){
		AVSerror( "InMatrices is not a valid matrix stream." );
		return( 0 );
	}
	MAT_fprint( fp, nmats, mats );
	fclose( fp );

	free( mats );

	return( 1 );
}
 
/* ***********************************************************************/
/* Initialization for modules contained in this file.                    */
/* ***********************************************************************/
static int ((*mod_list[])()) = {
	mw_desc
};
#define NMODS (sizeof(mod_list) / sizeof(char *))

AVSinit_modules()
{
	AVSinit_from_module_list(mod_list, NMODS);
}
