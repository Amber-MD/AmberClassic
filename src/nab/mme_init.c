/*
 * mme_init(): just a wrapper for the mme_init_sff() function in libsff
 */

#include <stdio.h>
#include <stdlib.h>

#include "nab.h"
#include "memutil.h"
#include "molutil.h"
#include "../sff/prm.h"

INT_T get_mytaskid(void);

/* these pointers will be used by various sff routines: */
static int *frozen = NULL;
static int *constrained = NULL;

/* output from molecular mechanics routines goes to nabout  */
/* nabout global should be defined only once */
FILE *nabout;

/***********************************************************************
                            MME_INIT()
************************************************************************/

/* Initialize many variables and data structures for mme, mme2, md, etc. */

int mme_init(MOLECULE_T * m, char *aexp, char *aexp2,
             REAL_T * x0i, char *ncfilename )
{

   int ier, nfrozen, nconstrained;
   PARMSTRUCT_T *prm;

   prm = m->m_prm;

   /* allocate the frozen and constrained vectors. */

   free_ivector(frozen, 0, prm->Natom);
   free_ivector(constrained, 0, prm->Natom);

   frozen = ivector(0, prm->Natom);
   constrained = ivector(0, prm->Natom);

   nconstrained = set_cons_mask(m, aexp2, constrained);
   if(get_mytaskid()==0){ 
     if (nconstrained) {
       if (aexp2 == NULL) {
         fprintf(nabout, "constrained all %d atoms\n", nconstrained);
       } else {
         fprintf(nabout, "constrained %d atoms using expression %s\n",
                 nconstrained, aexp2);
       }
     }
   }
   nfrozen = set_belly_mask(m, aexp, frozen);
   if(get_mytaskid()==0){ 
     if (nfrozen) {
       fprintf(nabout,
               "freezing %d atoms using expression %s for moving atoms\n",
               nfrozen, aexp);
     }
   }

   ier = mme_init_sff( prm, frozen, constrained, x0i, ncfilename );

   return (ier);
}

