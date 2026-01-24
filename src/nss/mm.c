#include <stdio.h>
#include <string.h>
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
int mm_desc()
{

	int in_port, out_port, param, iresult;
	extern int mm_compute();

	AVSset_module_name("mm", MODULE_FILTER);

	/* Input Port Specifications               */
	in_port = AVScreate_input_port( "InMatrices1", 
		"field 1D 1-space 1-vector uniform byte", REQUIRED );
	in_port = AVScreate_input_port( "InMatrices2", 
		"field 1D 1-space 1-vector uniform byte", REQUIRED );
	in_port = AVScreate_input_port( "InMatrices3", 
		"field 1D 1-space 1-vector uniform byte", OPTIONAL );
	in_port = AVScreate_input_port( "InMatrices4", 
		"field 1D 1-space 1-vector uniform byte", OPTIONAL );

	/* Output Port Specifications              */
	out_port = AVScreate_output_port( "OutMatrices", 
		"field 1D 1-space 1-vector uniform byte" );

	AVSset_compute_proc( mm_compute );
	return( 1 );
}
 
/* *****************************************/
/* Module Compute Routine                  */
/* *****************************************/
int mm_compute( InMatrices1, InMatrices2, InMatrices3, InMatrices4, 
	OutMatrices
)
	AVSfield_char *InMatrices1;
	AVSfield_char *InMatrices2;
	AVSfield_char *InMatrices3;
	AVSfield_char *InMatrices4;
	AVSfield_char **OutMatrices;
{
	int	nimats1, nimats2, nimats3, nimats4, nomats;
	MATRIX_T	*imats1, *imats2, *imats3, *imats4;
	int	dims0[ 1 ];
	char	*sip1, *sip2, *sip3, *sip4, *dp;
	char	m_si[ 100 ];
	int	s_sip1, s_sip2, s_sip3, s_sip4, s_o_sip;

	if( *OutMatrices )
		AVSfield_free( *OutMatrices );

	nimats1 = MAT_count( InMatrices1->data );
	sip1 = MAT_getsyminfo( InMatrices1->data );
	s_sip1 = strlen( sip1 );
	nimats2 = MAT_count( InMatrices2->data );
	sip2 = MAT_getsyminfo( InMatrices2->data );
	s_sip2 = strlen( sip2 );
	if( InMatrices3 ){
		nimats3 = MAT_count( InMatrices3->data );
		sip3 = MAT_getsyminfo( InMatrices3->data );
		s_sip3 = strlen( sip3 );
	}else{
		nimats3 = 0;
		sip3 = NULL;
		s_sip3 = 0;
	}
	if( InMatrices4 ){
		nimats4 = MAT_count( InMatrices4->data );
		sip4 = MAT_getsyminfo( InMatrices4->data );
		s_sip4 = strlen( sip4 );
	}else{
		nimats4 = 0;
		sip4 = NULL;
		s_sip4 = 0;
	}
	nomats = nimats1 + nimats2 + nimats3 + nimats4;

	imats1 = MAT_ALLOC( nimats1 );
	if( imats1 == NULL ){
		AVSerror( "Allocation of imats1 failed." );
		return( 0 );
	}
	imats2 = MAT_ALLOC( nimats2 );
	if( imats2 == NULL ){
		AVSerror( "Allocation of imats2 failed." );
		free( imats2 );
		return( 0 );
	}
	if( nimats3 != 0 ){
		imats3 = MAT_ALLOC( nimats3 );
		if( imats3 == NULL ){
			AVSerror( "Allocation of imats3 failed." );
			free( imats2 );
			free( imats1 );
			return( 0 );
		}
	}else
		imats3 = NULL;
	if( nimats4 != 0 ){
		imats4 = MAT_ALLOC( nimats4 );
		if( imats4 == NULL ){
			AVSerror( "Allocation of imats4 failed." );
			if( imats3 != NULL )
				free( imats3 );
			free( imats2 );
			free( imats1 );
			return( 0 );
		}
	}else
		imats4 = NULL;

	MAT_sscan( InMatrices1->data, nimats1, imats1 );
	MAT_sscan( InMatrices2->data, nimats2, imats2 );
	if( imats3 )
		MAT_sscan( InMatrices3->data, nimats3, imats3 );
	if( imats4 )
		MAT_sscan( InMatrices4->data, nimats4, imats4 );

	sprintf( m_si, "#S{ merge %d\n", getpid() );
	s_o_sip = strlen( m_si ) + s_sip1 + s_sip2 + s_sip3 + s_sip4 +
		strlen( "#S}\n" );

	dims0[ 0  ] = ( 16 * BPN * nomats ) + s_o_sip;
	*OutMatrices = (AVSfield_char * )AVSdata_alloc(
		"field 1D 1-space 1-vector uniform byte", dims0 );
	if( *OutMatrices == NULL ){
		AVSerror( "Allocation of output field failed." );
		free( imats1 );
		free( imats2 );
		if( imats3 != NULL )
			free( imats3 );
		if( imats4 != NULL )
			free( imats4 );
		return( 0 );
	}

	dp = ( *OutMatrices )->data;

	strcpy( dp, m_si );
	dp += strlen( m_si );
	strcpy( dp, sip1 );
	dp += s_sip1;
	strcpy( dp, sip2 );
	dp += s_sip2;
	if( sip3 ){
		strcpy( dp, sip3 );
		dp += s_sip3;
	}
	if( sip4 ){
		strcpy( dp, sip4 );
		dp += s_sip4;
	}
	strcpy( dp, "#S}\n" );
	dp += strlen( "#S}\n" );

	MAT_sprint( dp, nimats1, imats1 );
	dp += ( 16 * BPN * nimats1 );
	dp[ -1 ] = '\n';
	MAT_sprint( dp, nimats2, imats2 );
	dp += ( 16 * BPN * nimats2 );
	if( nimats3 > 0 ){
		dp[ -1 ] = '\n';
		MAT_sprint( dp, nimats3, imats3 );
		dp += ( 16 * BPN * nimats3 );
	}
	if( nimats4 > 0 ){
		dp[ -1 ] = '\n';
		MAT_sprint( dp, nimats4, imats4 );
	}

	free( imats1 );
	free( imats2 );
	if( imats3 != NULL )
		free( imats3 );
	if( imats4 != NULL )
		free( imats4 );

	return(1);
}
 
/* ***********************************************************************/
/* Initialization for modules contained in this file.                    */
/* ***********************************************************************/
static int ((*mod_list[])()) = {
	mm_desc
};
#define NMODS (sizeof(mod_list) / sizeof(char *))

AVSinit_modules()
{
	AVSinit_from_module_list(mod_list, NMODS);
}
