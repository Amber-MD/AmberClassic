#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "nab.h"

#define	HE_SIZE_0	251

/* call by NAB_href() as hashes are init on first ref	*/
void	NAB_hinit( HASH_T **htab, int type )
{
	int		i;
	H_ENTRY_T	**hep;

	if( ( *htab = ( HASH_T * )malloc( sizeof( HASH_T ) ) ) == NULL ){
		fprintf( stderr, "NAB_hinit: can't alloc new hashed array\n" );
		exit( 1 );
	}
	if( ( hep = ( H_ENTRY_T ** )malloc( HE_SIZE_0*sizeof( H_ENTRY_T * ) ) )
		== NULL ){
		fprintf( stderr, "NAB_hinit: can't alloc new he tab\n" );
		exit( 1 );
	}
	for( i = 0; i < HE_SIZE_0; i++ )
		hep[ i ] = NULL;

	( *htab )->h_type = type;
	( *htab )->h_llen = 0;
	( *htab )->h_n_entries = 0;
	( *htab )->h_he_size = HE_SIZE_0;
	( *htab )->h_entries = hep;
}

VALUE_T	*NAB_href( HASH_T **htab, char key[], int type, int size )
{
	int	klen;
	int	hash;
	char	*kp;
	int	cv, ll;
	H_ENTRY_T	*hep, *hepl, *hepk;

	if( *htab == NULL )
		NAB_hinit( htab, type );

	klen = strlen( key );
	hash = ( key[ 0 ] + key[ klen - 1 ] + 7 * klen ) % (*htab)->h_he_size;

	hepk = NULL;	/* key if found	*/
	hepl = NULL;	/* prev entry for maintaining ordered list	*/
	ll = 0;		/* longest list	*/
	for( hep = (*htab)->h_entries[ hash ]; hep; hep = hep->he_next ){
		if( ( cv = strcmp( hep->he_key, key ) ) > 0 )
			break;
		else if( cv == 0 ){
			hepk = hep;
			break;
		}
		hepl = hep;
		ll++;
	}

	if( !hepk ){
		if( ( hepk = ( H_ENTRY_T * )malloc( sizeof( H_ENTRY_T ) ) )
			== NULL ){
			fprintf( stderr, "NAB_href: can't add %s to hash\n",
				key );
			exit( 1 );
		}
		if(( kp = (char *)malloc( (klen+1) * sizeof(char))) == NULL ){
			fprintf( stderr, "NAB_href: can't save %s in hash\n",
				key );
			exit( 1 );
		}
		strcpy( kp, key );

		if( ++ll > (*htab)->h_llen )
			(*htab)->h_llen = ll;
		(*htab)->h_n_entries++;
		hepk->he_next = hep;
		if( !hepl )
			(*htab)->h_entries[ hash ] = hepk;
		else
			hepl->he_next = hepk;
		hepk->he_key = kp;
		hepk->he_val.v_type = (*htab)->h_type;
		switch( (*htab)->h_type ){
		case T_INT :
			hepk->he_val.v_value.v_ival = 0;
			break;
		case T_SIZE_T :
			hepk->he_val.v_value.v_szval = 0;
			break;
		case T_FLOAT :
			hepk->he_val.v_value.v_fval = 0.0;
			break;
		case T_STRING :
			hepk->he_val.v_value.v_cval = NULL; /*strdup("");*/
			break;
		case T_POINT :
			hepk->he_val.v_value.v_ptval[ 0 ] = 0.0;
			hepk->he_val.v_value.v_ptval[ 1 ] = 0.0;
			hepk->he_val.v_value.v_ptval[ 2 ] = 0.0;
			break;
		case T_MATRIX :
			hepk->he_val.v_value.v_matval = NULL;
			break;
		case T_FILE :
			hepk->he_val.v_value.v_fpval = NULL;
			break;
		case T_ATOM :
			hepk->he_val.v_value.v_atomval = NULL;
			break;
		case T_RESIDUE :
			hepk->he_val.v_value.v_resval = NULL;
			break;
		case T_MOLECULE :
			hepk->he_val.v_value.v_molval = NULL;
			break;
		case T_BOUNDS :
			hepk->he_val.v_value.v_bval = NULL;
			break;
		case T_USER :
			hepk->he_val.v_value.v_uval =
				( char * )calloc( (size_t)size, sizeof(char));
			if( hepk->he_val.v_value.v_uval == NULL ){
				fprintf( stderr, "NAB_href: can't calloc space for %s in hash\n",
					key );
				exit( 1 );
			}
		}
	}

	return( &hepk->he_val );
}

int	NAB_hfirst( HASH_T *htab, CURHASH_T *c_hash )
{
	c_hash->c_index = 0;
	c_hash->c_entry = NULL;
	return( 0 );
}

char	*NAB_hnext( HASH_T *htab, CURHASH_T *c_hash )
{
	int	i;
	H_ENTRY_T	*hep;

	if( htab == NULL )
		return( NULL );

	if( !c_hash->c_entry ){
		for( i = 0; i < htab->h_he_size; i++ ){
			if( (hep = htab->h_entries[ i ]) ){
				c_hash->c_index = i;
				c_hash->c_entry = hep;
				return( hep->he_key );
			}
		}
		return( NULL );
	}else if( (hep = c_hash->c_entry->he_next) ){
		c_hash->c_entry = hep;
		return( hep->he_key );
	}else{
		for( i = c_hash->c_index + 1; i < htab->h_he_size; i++ ){
			if( (hep = htab->h_entries[ i ]) ){
				c_hash->c_index = i;
				c_hash->c_entry = hep;
				return( hep->he_key );
			}
		}
		return( NULL );
	}
}

int	NAB_hin( HASH_T *htab, char key[] )
{
	int	klen;
	int	cv, hash;
	H_ENTRY_T	*hep;

	if( htab == NULL )
		return( 0 );

	klen = strlen( key );
	hash = ( key[ 0 ] + key[ klen - 1 ] + 7 * klen ) % htab->h_he_size;

	for( hep = htab->h_entries[ hash ]; hep; hep = hep->he_next ){
		if( ( cv = strcmp( hep->he_key, key ) ) > 0 )
			return( 0 );	/* not in tab	*/
		else if( cv == 0 ){
			return( 1 );	/* found it!	*/
		}
	}
	return( 0 );	/* not in tab	*/
}

void	NAB_hdelete( HASH_T *htab, char key[] )
{
	int	klen;
	int	cv, hash;
	H_ENTRY_T	*hep, *hep1;

	if( htab == NULL )
		return;

	/* clear the array */
	if( !key ){ 
		for( hash = 0; hash < htab->h_he_size; hash++ ){
			for( hep = htab->h_entries[ hash ]; hep; hep = hep1 ){
				hep1 = hep->he_next;
				free( hep );
			}
			htab->h_entries[ hash ] = NULL;
		}
		htab->h_n_entries = 0;
		return;
	}

	/* delete a key */
	klen = strlen( key );
	hash = ( key[ 0 ] + key[ klen - 1 ] + 7 * klen ) % htab->h_he_size;
	hep1 = NULL;
	for( hep = htab->h_entries[ hash ]; hep; hep = hep->he_next ){
		if( ( cv = strcmp( hep->he_key, key ) ) > 0 )
			return;		/* not in table, return */
		else if( cv == 0 ){
			if( !hep1 )
				htab->h_entries[ hash ] = hep->he_next;
			else
				hep1->he_next = hep->he_next;
			free( hep );
			return;		/* found it!	*/
		}
		hep1 = hep;
	}
}

void	NAB_hdump( FILE *fp, char *msg, HASH_T *htab, int pzero )
{
	int	i, ec;
	H_ENTRY_T	*hep;
	
	if( msg != NULL && *msg != '\0' )
		fprintf( fp, "HD: %s\n", msg );
	if( htab == NULL ){
		fprintf( fp, "htab == NULL\n" );
		return;
	}
	fprintf( fp, "HD: type\t%d\n", htab->h_type );
	fprintf( fp, "HD: llen\t%d\n", htab->h_llen );
	fprintf( fp, "HD: #ent\t%d\n", htab->h_n_entries );
	fprintf( fp, "HD: size\t%d\n", htab->h_he_size );
	for( i = 0; i < htab->h_he_size; i++ ){
		for( ec = 0, hep = htab->h_entries[ i ]; hep;
			hep = hep->he_next )
		{
			ec++;
		}
		if( ec > 0 || pzero ){
			fprintf( fp, "HD: htab[%3d] %3d entries\n", i, ec );
			for( ec = 0, hep = htab->h_entries[ i ]; hep;
				hep = hep->he_next )
			{
				fprintf( fp, "HD:\t%p %p\t'%s'", hep,
					hep->he_next, hep->he_key );
				switch( hep->he_val.v_type ){
				case T_UNDEF :
					fprintf( fp, "\tU\t%d",
						hep->he_val.v_value.v_ival );
					break;
				case T_INT :
					fprintf( fp, "\tI\t%d",
						hep->he_val.v_value.v_ival );
					break;
				case T_SIZE_T :
					fprintf( fp, "\tS\t%ld",
						hep->he_val.v_value.v_szval );
					break;
				case T_FLOAT :
					fprintf( fp, "\tF\t%g",
						hep->he_val.v_value.v_fval );
					break;
				case T_STRING :
					fprintf( fp, "\ts\t'%s'",
						hep->he_val.v_value.v_cval ?
						hep->he_val.v_value.v_cval :
						"(NULL)" );
					break;
				case T_POINT :
					fprintf( fp, "\tP" );
					break;
				case T_MATRIX :
					fprintf( fp, "\tM" );
					break;
				case T_FILE :
					fprintf( fp, "\tf" );
					break;
				case T_ATOM :
					fprintf( fp, "\ta" );
					break;
				case T_RESIDUE :
					fprintf( fp, "\tr" );
					break;
				case T_MOLECULE :
					fprintf( fp, "\tm" );
					break;
				case T_BOUNDS :
					fprintf( fp, "\tb" );
					break;
				default :
					fprintf( fp, "\t?" );
					break;
				}
				fprintf( fp, "\n" );
			}
		}
	}
}
