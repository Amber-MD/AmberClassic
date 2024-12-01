#ifndef  SELECT_ATOMS_H
#define  SELECT_ATOMS_H

int	setpoint( MOLECULE_T *mol, char aexpr[], POINT_T point );
int	select_atoms( MOLECULE_T *mol, char aex[] );
int	atom_in_aexpr( ATOM_T *ap, char aex[] );
void	set_attr_if( MOLECULE_T *mol, int attr, int i_attr );
void	clear_attr( MOLECULE_T *mol, int attr );

#endif
