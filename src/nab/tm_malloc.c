#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "tm_malloc.h"

#define TM_MEM_SIZE     20000
static  TM_MEM_T        tm_mem[ TM_MEM_SIZE ];
static  int     tm_n_mem;

char    *tm_malloc( size_t size, char id[] )
{
        char    *ptr;

        ptr = ( char * )malloc( size * sizeof( char ) );
#ifdef  MEMDEBUG
        if( tm_n_mem < TM_MEM_SIZE ){
                tm_mem[ tm_n_mem ].t_freed = 0;
                tm_mem[ tm_n_mem ].t_id = id;
                tm_mem[ tm_n_mem ].t_addr = ptr;
                tm_n_mem++;
        }
#endif
        return( ptr );
}

void    tm_free( void *ptr )
{
#ifdef  MEMDEBUG
        int     i;
#endif

        free( ptr );
#ifdef  MEMDEBUG
        for( i = 0;i < tm_n_mem; i++ ){
                if( tm_mem[ i ].t_addr == ptr ){
                        tm_mem[ i ].t_freed = 1;
                }
        }
#endif
}

void    tm_report( void )
{
        int     i;
        static  int     cnt = 0;

        for( i = 0; i < tm_n_mem; i++ ){
                fprintf( stderr, "tm_mem[%4d]:", i );
                fprintf( stderr, " addr = %8p, freed = %d, id = %s\n",
                        tm_mem[i].t_addr, tm_mem[i].t_freed,
                        tm_mem[i].t_id );
        }
        cnt++;
/*
        if( cnt > 1 )
                exit( 1 );
*/
}

