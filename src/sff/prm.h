#ifndef  PRM_H
#define  PRM_H

/*  prmtop routine public interfaces: */
INT_T    free_prm( PARMSTRUCT_T* prm );
PARMSTRUCT_T*  rdparm( STRING_T* name );
INT_T    wrparm( PARMSTRUCT_T* prm, STRING_T* name );

/*  mme routine public interfaces: */
INT_T    mme_init_sff( PARMSTRUCT_T*, INT_T*, INT_T*, REAL_T*, STRING_T* );

#endif
