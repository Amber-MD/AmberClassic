#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "nab.h"
#include "molutil.h"
#include "../sff/prm.h"

/*
 * These routines are NAB wrappers for the corresponding SFF routines.
 */

/***************************************************************************
			     FREEPARM()
****************************************************************************/

/*
 * freeparm() - free a PARMSTRUCT_T
 */

int	freeparm( MOLECULE_T *mol )
{
	if( !mol )
		return( 0 );

	if( mol->m_prm )
		free_prm( mol->m_prm );

	return( 0 );
}

/***************************************************************************
			     READPARM()
****************************************************************************/

/*
 * readparm() - read a prmtop file and instantiate a PARMSTRUCT_T
 */

int readparm(MOLECULE_T * mol, char *name)
{

   int ai;
   ATOM_T *a;

   mol->m_prm = rdparm( name );

   /*  fill in the charge and radii arrays in the molecule with the values
      we just got from the prmtop file:                                  */

   for (ai = 0, a = NULL; (a = NAB_mnext(mol, a)); ai++) {
      a->a_charge = mol->m_prm->Charges[ai] / 18.2223;
      a->a_radius = mol->m_prm->Rborn[ai];
   }

   return (0);
}
#if 0

/***************************************************************************
			     WRITEPARM()
****************************************************************************/

/*
 * writeparm() - write a PARMSTRUCT_T but only if task zero
 */

int writeparm(MOLECULE_T * mol, char *name)
{
   int ier;

   ier = wrparm( mol->m_prm, name );

   return (ier);
}
#endif
