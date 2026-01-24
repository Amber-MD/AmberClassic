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
int mx_desc()
{

	int in_port, out_port, param, iresult;
	extern int mx_compute();

	AVSset_module_name( "mx", MODULE_FILTER );

	/* Input Port Specifications               */
	in_port = AVScreate_input_port( "RightMatrices", 
		"field 1D 1-space 1-vector uniform byte", REQUIRED );
	in_port = AVScreate_input_port( "LeftMatrices", 
		"field 1D 1-space 1-vector uniform byte", REQUIRED );

	/* Output Port Specifications              */
	out_port = AVScreate_output_port( "OutMatrices", 
		"field 1D 1-space 1-vector uniform byte" );

	AVSset_compute_proc( mx_compute );
	return( 1 );
}
 
/* *****************************************/
/* Module Compute Routine                  */
/* *****************************************/
int mx_compute( RightMatrices, LeftMatrices, OutMatrices )
	AVSfield_char *RightMatrices;
	AVSfield_char *LeftMatrices;
	AVSfield_char **OutMatrices;
{
	int	m, nrmats, nlmats, nomats;
	MATRIX_T	*rmats, *lmats, *omats;
	int	dims0[ 1 ];
	char	*sipr, *sipl;
	char	m_si[ 100 ];
	int	s_sipr, s_sipl, s_o_sip;
	char	*dp;

	if( *OutMatrices )
		AVSfield_free( *OutMatrices );

	nrmats = MAT_count( RightMatrices->data );
	sipr = MAT_getsyminfo( RightMatrices->data );
	s_sipr = strlen( sipr );
	nlmats = MAT_count( LeftMatrices->data );
	sipl = MAT_getsyminfo( LeftMatrices->data );
	s_sipl = strlen( sipl );
	nomats = ( nrmats > nlmats ) ? nrmats : nlmats;

	rmats = MAT_ALLOC( nrmats );
	if( rmats == NULL ){
		AVSerror( "Allocation of rmats failed." );
		return( 0 );
	}
	lmats = MAT_ALLOC( nlmats );
	if( lmats == NULL ){
		AVSerror( "Allocation of lmats failed." );
		free( rmats );
		return( 0 );
	}
	omats = MAT_ALLOC( nomats );
	if( omats == NULL ){
		AVSerror( "Allocation of omats failed." );
		free( rmats );
		free( lmats );
		return( 0 );
	}

	MAT_sscan( RightMatrices->data, nrmats, rmats );
	MAT_sscan( LeftMatrices->data, nlmats, lmats );

	sprintf( m_si, "#S{ multiply %d\n", getpid() );
	s_o_sip = strlen( m_si ) + strlen( "#S}\n" ) + s_sipr + s_sipl;

	if( nlmats == nrmats ){
		for( m = 0; m < nlmats; m++ ){
			NAB_matcpy( omats[ m ],
				MAT_concat( lmats[ m ], rmats[ m ] ) );
		}
	}else if( nlmats < nrmats ){
		for( m = 0; m < nlmats; m++ ){
			NAB_matcpy( omats[ m ],
				MAT_concat( lmats[ m ], rmats[ m ] ) );
		}
		for( m = nlmats; m < nrmats; m++ ){
			NAB_matcpy( omats[ m ],
				MAT_concat( lmats[ nlmats - 1 ], rmats[ m ] ) );
		}
	}else{
		for( m = 0; m < nrmats; m++ ){
			NAB_matcpy( omats[ m ],
				MAT_concat( lmats[ m ], rmats[ m ] ) );
		}
		for( m = nrmats; m < nlmats; m++ ){
			NAB_matcpy( omats[ m ],
				MAT_concat( lmats[ m ], rmats[ nrmats - 1 ] ) );
		}
	}

	dims0[ 0  ] = ( 16 * BPN * nomats ) + s_o_sip;
	*OutMatrices = (AVSfield_char * )AVSdata_alloc(
		"field 1D 1-space 1-vector uniform byte", dims0 );
	if( *OutMatrices == NULL ){
		AVSerror( "Allocation of output field failed." );
		free( omats );
		free( lmats );
		free( rmats );
		return( 0 );
	}
	dp = ( *OutMatrices )->data;
	strcpy( dp, m_si );
	dp += strlen( m_si );
	strcpy( dp, sipl );
	dp += s_sipl;
	strcpy( dp, sipr );
	dp += s_sipr;

	MAT_sprint( dp, nomats, omats );
	dp += ( 16 * BPN * nomats );
	strcpy( dp, "#S}\n" );

	free( omats );
	free( lmats );
	free( rmats );
	
	return( 1 );
}
 
/* ***********************************************************************/
/* Initialization for modules contained in this file.                    */
/* ***********************************************************************/
static int ((*mod_list[])()) = {
	mx_desc
};
#define NMODS (sizeof(mod_list) / sizeof(char *))

AVSinit_modules()
{
	AVSinit_from_module_list(mod_list, NMODS);
}
