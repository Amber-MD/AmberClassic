#ifndef	TM_MALLOC_H
#define	TM_MALLOC_H

typedef struct  tm_mem_t        {
        int     t_freed;
        char    *t_id;
        void    *t_addr;
} TM_MEM_T;

char    *tm_malloc( size_t, char [] );
void    tm_free( void * );
void    tm_report( void );

#endif
