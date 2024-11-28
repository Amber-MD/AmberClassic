!<compile=optimized>

! The 3D-RISM-KH software found here is copyright (c) 2010-2012 by Andriy Kovalenko,
! Tyler Luchko and David A. Case.
! 
! This program is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.
! 
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
! for more details.
! 
! You should have received a copy of the GNU General Public License in the
! ../../LICENSE file.  If not, see <http://www.gnu.org/licenses/>.
! 
! Users of the 3D-RISM capability found here are requested to acknowledge
! use of the software in reports and publications.  Such acknowledgement
! should include the following citations:
! 
! 1) A. Kovalenko and F. Hirata. J. Chem. Phys., 110:10095-10112  (1999);
! ibid. 112:10391-10417 (2000).
! 
! 2) A. Kovalenko,  in:  Molecular  Theory  of  Solvation,  edited  by
! F. Hirata  (Kluwer Academic Publishers, Dordrecht, 2003), pp.169-275.
! 
! 3) T. Luchko, S. Gusarov, D.R. Roe, C. Simmerling, D.A. Case, J. Tuszynski,
! and  A. Kovalenko, J. Chem. Theory Comput., 6:607-624 (2010).

#include "../include/dprec.fh"

!> Handles quaternions.  Currently quaterions are just 1D arrays w/ 4
!! indicies.  This may be expanded to define a quaterion type with
!! overloaded operators.
module quaternion
  implicit none

contains

  !> New quaternion from Euler angle and rotation axis. The rotation axis is
  !! automatically normalized. If dir = (0, 0, 0), then no
  !! normalization is done and the imaginary parts are left as zero.
  !! This does not change the angular part.
  !! @param[in,out] angle Angle, in radians, to rotate by.
  !! @param[in,out] dir Vector to rotate about.
  !! @return quaternion (4-vector) corresponding to the angle and direction input
  function quaternionFromEulerAxis(angle, direction) result(quaternion)
    implicit none
    _REAL_, intent(in) :: angle, direction(3)
    _REAL_ :: quaternion(4), magnitude

    magnitude  = sqrt(sum(direction**2))
    quaternion(1) = cos(angle / 2d0)
    if (magnitude == 0d0) then
       quaternion(2:4) = 0d0
    else
       quaternion(2:4) = direction / magnitude * sin(angle / 2d0)
    end if
  end function quaternionFromEulerAxis


  !> New quaternion from Euler angles. Right-handed coordinate system
  !! is assumed.
  !! @param[in] a The yaw angle in radians.
  !! @param[in] b The pitch angle in radians.
  !! @param[in] c The roll angle in radians.
  !! @return quaternion (4-vector) corresponding to Euler angles
  function quaternionFromEulerAngles(a, b, c) result(quaternion)
    implicit none
    _REAL_, intent(in) :: a, b, c
    _REAL_ :: quaternion(4)
    
    ! quaternion(1) = cos((a - c) / 2) * sin(b / 2)
    ! quaternion(2) = sin((a - c) / 2) * sin(b / 2)
    ! quaternion(3) = sin((a + c) / 2) * cos(b / 2)
    ! quaternion(4) = cos((a + c) / 2) * cos(b / 2)
    quaternion(1) = cos(a/2) * cos(b/2) * cos(c/2) + sin(a/2) * sin(b/2) * sin(c/2)
    quaternion(2) = sin(a/2) * cos(b/2) * cos(c/2) - cos(a/2) * sin(b/2) * sin(c/2)
    quaternion(3) = cos(a/2) * sin(b/2) * cos(c/2) + sin(a/2) * cos(b/2) * sin(c/2)
    quaternion(3) = cos(a/2) * cos(b/2) * sin(c/2) - sin(a/2) * sin(b/2) * cos(c/2)
  end function quaternionFromEulerAngles

  
  !> Converts a quaternion to an angle and rotation axis. The rotation
  !! axis is automatically normalized.
  !! @param[in] quat A quaternion.
  !! @return the angle and direction corresponding to the input quaternion
  !!     (4-vector)
  function quat_2_angle_axis(quat)
    implicit none
    _REAL_, intent(in) :: quat(4)
    _REAL_ :: quat_2_angle_axis(4), sin_ang
    sin_ang = sqrt(sum(quat(2:4)**2))
    quat_2_angle_axis(2:4) = quat(2:4) / sin_ang
    quat_2_angle_axis(1) = 2d0 * atan2(sin_ang, quat(1))
  end function quat_2_angle_axis


  !> Rotates a vector by the given quaterion.
  !! @param[in,out] vec 3D vector to be rotated.
  !! @param[in] quat Quaternion that defines the rotation.
  subroutine rotate_quat(vec, quat)
    implicit none
    _REAL_,intent(inout) :: vec(3)
    _REAL_,intent(in) :: quat(4)
    _REAL_ :: rot(3),t(10),qvec(4),qtemp(4),q(4)
    t(2) =   quat(1)*quat(2)
    t(3) =   quat(1)*quat(3)
    t(4) =   quat(1)*quat(4)
    t(5) =  -quat(2)*quat(2)
    t(6) =   quat(2)*quat(3)
    t(7) =   quat(2)*quat(4)
    t(8) =  -quat(3)*quat(3)
    t(9) =   quat(3)*quat(4)
    t(10)=  -quat(4)*quat(4)
    rot(1) =2*((t(8)+t(10))*vec(1)+(t(6)- t(4))*vec(2)+(t(3)+t(7))*vec(3))+vec(1)
    rot(2) =2*((t(4)+ t(6))*vec(1)+(t(5)+t(10))*vec(2)+(t(9)-t(2))*vec(3))+vec(2)
    rot(3) =2*((t(7)- t(3))*vec(1)+(t(2)+ t(9))*vec(2)+(t(5)+t(8))*vec(3))+vec(3)
    qvec(1) = 0
    qvec(2:4) = vec
    call quat_mult(quat,qvec,qtemp)
    q(1) = quat(1)
    q(2:4) = -quat(2:4)
    call quat_mult(qtemp,q,qvec)
    vec = qvec(2:4)
  end subroutine rotate_quat


  !> Multiplies quaternion a and b to produce c.
  !! @param[in] a Quaternion a.
  !! @param[in] b Quaternion b.
  !! @param[out] c Product quaternion.
  subroutine quat_mult(a, b, c)
    implicit none
    _REAL_, intent(in) :: a(4), b(4)
    _REAL_, intent(out) :: c(4)
    _REAL_ :: d(4), t(3)
    d(1) = a(1) * b(1) - dot_product(a(2:4), b(2:4))
    d(2) = a(1) * b(2) + a(2) * b(1) +a(3) * b(4) - a(4) * b(3)
    d(3) = a(1) * b(3) -a(2) * b(4) + a(3) * b(1) + a(4) * b(2)
    d(4) = a(1) * b(4) + a(2) * b(3) - a(3) * b(2) + a(4) * b(1)
    c = d
  end subroutine quat_mult

  
end module quaternion
