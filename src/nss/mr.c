#include <stdio.h>
#include <avs/avs.h>
#include <avs/port.h>
#include <avs/field.h>

#include "nabcode.h"

#define	MAT_ALLOC(n)	( ( MATRIX_T * )malloc( ( n )*sizeof( MATRIX_T ) ) )
#define	BPN	14

static	int	pid = -1;
static	char	mname[ 256 ] = "";
 
/* *****************************************/
/*  Module Description                     */
/* *****************************************/
int mr_desc()
{

	int in_port, out_port, param, iresult;
	extern int mr_compute();

	AVSset_module_name( "mr", MODULE_DATA );

	/* Output Port Specifications              */
	out_port = AVScreate_output_port( "OutMatrices", 
		"field 1D 1-space 1-vector uniform byte" );

	/* Parameter Specifications                */
	param = AVSadd_parameter( "MatFile", "string", "", "", ":" );
	AVSconnect_widget( param, "browser" );

	param = AVSadd_parameter( "Name", "string", "", "", "" );
	AVSadd_parameter_prop( param, "width", "integer", 4 );
	AVSconnect_widget( param, "typein" );

	AVSset_compute_proc( mr_compute );
	return( 1 );
}
 
/* *****************************************/
/* Module Compute Routine                  */
/* *****************************************/
int mr_compute( OutMatrices, MatFile, Name)
	AVSfield_char **OutMatrices;
	char *MatFile;
	char *Name;
{
	FILE	*fp;
	char	e_msg[ 256 ];
	int	smats, nmats;
	MATRIX_T	*mats;
	int	dims0[ 1 ];
	int	sisize;
	char	sibuf[ 1000 ];

	if( pid == -1 )
		pid = getpid();

	if( !*Name ){
		if( !*mname )
			sprintf( mname, "m%d", pid );
	}else
		strcpy( mname, Name );
	AVSmodify_parameter( "Name", AVS_VALUE, mname, "", "" );

	if( *OutMatrices )
		AVSfield_free( *OutMatrices );

	if( !*MatFile )
		return( 0 );
	if( ( fp = fopen( MatFile, "r" ) ) == NULL ){
		sprintf( e_msg, "Can't read MatFile %s.", MatFile );
		AVSerror( e_msg );
		return( 0 );
	}
	smats = 1000;
	mats = MAT_ALLOC( smats );
	if( mats == NULL ){
		AVSerror( "Allocation of mats failed." );
		return( 0 );
	}

	if( ( nmats = MAT_fscan( fp, smats, mats ) ) == 0 ){
		AVSerror( "MatFile is not a valid matrix file." );
		fclose( fp );
		return( 0 );
	}

	sprintf( sibuf, "#S{ read %d %s\n#S+   file     %s\n#S}\n",
		pid, mname, MatFile );
	sisize = strlen( sibuf );

	dims0[ 0 ] = 16 * BPN * nmats + sisize;
	*OutMatrices = ( AVSfield_char * )AVSdata_alloc(
		"field 1D 1-space 1-vector uniform byte", dims0 );
	if( *OutMatrices == NULL ){
		AVSerror( "Allocation of output field failed." );
		return( 0 );
	}

	strcpy( ( *OutMatrices )->data, sibuf );

	MAT_sprint( &( *OutMatrices )->data[ sisize ], nmats, mats );

	free( mats );

	return( 1 );
}
 
/* ***********************************************************************/
/* Initialization for modules contained in this file.                    */
/* ***********************************************************************/
static int ((*mod_list[])()) = {
	mr_desc
};
#define NMODS (sizeof(mod_list) / sizeof(char *))

AVSinit_modules()
{
	AVSinit_from_module_list(mod_list, NMODS);
}
