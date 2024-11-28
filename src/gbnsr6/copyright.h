!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! GBNSR6
!
! An analysis program for solvent-mediated electrostatic and non-electrostatic
! interactions in biomolecules. Can be used for efficient estimates of 
! solvation free energies. 
!
! Please acknowledge your use of GBNSR6 by citing:
!
! S. Izadi, B. Aguilar, and A. V. Onufriev,  J. Chem. Theory and Comput. 11(9), 4450–4459 (2015)
! B. Aguilar and A.V. Onufriev. J. Chem. Theory and Comput., 8, 2404-2411 (2012)
! B. Aguilar, R. Shadrach, and A. V. Onufriev, J. Chem. Theory and Comput. 6, 3613(2010)
! A.  Mukhopadhyay, B. Aguilar, Igor S. Tolokh, and A. V. Onufriev, J. Chem. Theory and Comput., 10 1788–1794 (2014)
!
! Major Developers:
!
! Saeed Izadi and Boris Aguilar, with some helps from Abhishek Mukhopadhyay
! Alexey V. Onufriev, coordinator of overall development
!
! Departments of Computer Science and Physics,
! Virginia Tech, Blacksburg, Virginia
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Overview:
!
! GBNSR6 is an implementation of the Generalized Born (GB) model in which
! the efective Born radii are computed numerically, via the so-called R6 integration 
! over molecular surface of the solute. For most structures, solvation energies computed 
! via the GB formula based on the numerical R6 radii are virtually as accuarte as the 
! energies based on the gold standard perfect effective radii, which can in 
! principle be obtained from numerical solution of the PB equation.
! As a result, the numerical R6 formulation is generally more accurate than
! the fast analytical approaches available in AMBER.
! In contrast to most GB practical models, NSR6 model is parameter-free in
! the same sense as the numerical PB framework is.
! Thus, accuracy of NSR6 relative to the PB standard is virtually unaffacted
! by the choice of input atomic radii.
! Three versions of the GB equation are available in this implementation:
! (1) the canonical (Still 1990) GB
! (2) the canonical GB with the ALPB correction
! (3) the charge hydration asymmetric generalized Born model.
! The first two models are described in more detail in the GB section of
! the main manual. For more information about the methods please refer to the 
! original references (above).
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! GBNSR6 is free software; you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! GBNSR6 is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! GBNSR6; if not, write to the Free Software Foundation, Inc., 59 Temple Place,
! Suite 330, Boston, MA 02111-1307, USA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
