#include <stdio.h>
#include <string.h>
#include <avs/avs.h>
#include <avs/port.h>
#include <avs/field.h>

#define	NBUFSIZE	10240
static	char	nbuf[ NBUFSIZE ];
static	int	getobjnames();

/* *****************************************/
/*  Module Description                     */
/* *****************************************/
int ron_desc()
{

	int in_port, out_port, param, iresult;
	extern int ron_compute();

	AVSset_module_name("ron", MODULE_DATA);

	/* Output Port Specifications              */
	out_port = AVScreate_output_port("ObjNames", 
		"field 1D 1-space 1-vector uniform byte");

	/* Parameter Specifications                */
	param = AVSadd_parameter("ObjNameFile", "string", "", "", ":");
	AVSconnect_widget(param, "browser");

	AVSset_compute_proc(ron_compute);
	return(1);
}
 
/* *****************************************/
/* Module Compute Routine                  */
/* *****************************************/
int ron_compute( ObjNames, ObjNameFile)
	AVSfield_char **ObjNames;
	char *ObjNameFile;
{
	FILE	*fp;
	char	e_msg[ 256 ];
	int	dims0[ 1 ];

	if( *ObjNames )
		AVSfield_free( *ObjNames );

	if( !*ObjNameFile )
		return( 0 );
	if( ( fp = fopen( ObjNameFile, "r" ) ) == NULL ){
		sprintf( e_msg, "Can't read ObjNameFile %s.", ObjNameFile );
		AVSerror( e_msg );
		return( 0 );
	}
	dims0[ 0 ] = getobjnames( fp, nbuf ) + 1;	/* add \0 */
	fclose( fp );

	if( dims0[ 0 ] == 0 ){
		AVSerror( "ObjNamesFile is invalid." );
		return( 0 );
	}
	
	*ObjNames = ( AVSfield_char * )AVSdata_alloc(
		"field 1D 1-space 1-vector uniform byte", dims0 );
	if( *ObjNames == NULL ){
		AVSerror( "Allocation of output field failed." );
		return( 0 );
	}
	strcpy( ( *ObjNames )->data, nbuf );

	return(1);
}

static	int	getobjnames( fp, nbuf )
FILE	*fp;
char	nbuf[];
{
	char	*nbp, *lp, line[ 256 ];

	for( nbp = nbuf; fgets( line, sizeof( line ), fp ); ){
		if( *line == '#' )
			continue;
		lp = line + strspn( line, " \t" );
		strcpy( nbp, lp );
		nbp += strlen( lp );
	}
	*nbp = '\0';
	return( nbp - nbuf );
}
 
/* ***********************************************************************/
/* Initialization for modules contained in this file.                    */
/* ***********************************************************************/
static int ((*mod_list[])()) = {
	ron_desc
};
#define NMODS (sizeof(mod_list) / sizeof(char *))

AVSinit_modules()
{
	AVSinit_from_module_list(mod_list, NMODS);
}
