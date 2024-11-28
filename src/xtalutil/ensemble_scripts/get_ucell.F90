#define _REAL_ double precision
!---------------------------------------------------------------------

!     --- GET_UCELL ---

!     ...this routine produces the direct and reciprocal lattice
!     vectors from the unit cell edge lengths and angles which are
!     passed to it.  It is assumed that the 1st vector (length a)
!     lies along the cartesian x-axis the 2nd vector (length b) is
!     in the x-y plane with positive y, and that the direct lattice
!     vectors are a non-degenerate right handed system.  Thus the 3rd
!     vector has positive z component.  Alpha is the angle (in degrees)
!     between 2nd and 3rd vectors, beta is the angle (in degrees)
!     between 1st and 3rd vectors, and gamma is the angle (in degrees)
!     between 1st and 2nd vectors.  The lengths of the direct lattice
!     vectors are given by dirlng(1), dirlng(2) and dirlng(3),
!     whereas the reciprocal lengths of the reciprocal vectors are
!     reclng(1),reclng(2)  and reclng(3).


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine get_ucell here]
program get_ucell
   implicit none
   _REAL_ ::  a,b,c,alpha,beta,gamma
   _REAL_ ::  sphere,volume
   _REAL_ ::  ucell(3,3),recip(3,3),dirlng(3),reclng(3)

   _REAL_ u23(3),u31(3),u12(3)
   _REAL_ result,distance,onevolume
   integer i,j
   _REAL_, parameter :: PI      = 3.1415926535897932384626433832795d0
   _REAL_, parameter :: DEG_TO_RAD = PI / 180.0d0

   ! input unit cell parameters:
   read(5,*) a,b,c,alpha,beta,gamma
   write(6,*) 'input unit cell parmameters:'
   write(6,'(3f14.7)') a,b,c
   write(6,'(3f14.7)') alpha,beta,gamma

   ucell(1,1) = a
   ucell(2,1) = 0.d0
   ucell(3,1) = 0.d0
   ucell(1,2) = b*cos(DEG_TO_RAD*gamma)
   ucell(2,2) = b*sin(DEG_TO_RAD*gamma)
   ucell(3,2) = 0.d0
   ucell(1,3) = c*cos(DEG_TO_RAD*beta)
   ucell(2,3) = &
         (b*c*cos(DEG_TO_RAD*alpha)-ucell(1,3)*ucell(1,2))/ucell(2,2)
   ucell(3,3) = sqrt( c*c - ucell(1,3)*ucell(1,3) - &
         ucell(2,3)*ucell(2,3) )
   dirlng(1) = a
   dirlng(2) = b
   dirlng(3) = c

   !  now get reciprocal vectors

   call cross(ucell(1,2),ucell(1,3),u23)
   call cross(ucell(1,3),ucell(1,1),u31)
   call cross(ucell(1,1),ucell(1,2),u12)
   call dot(ucell(1,1),u23,volume)
   onevolume=1.0d0/volume
   do j = 1,3
      recip(j,1) = u23(j)*onevolume
      recip(j,2) = u31(j)*onevolume
      recip(j,3) = u12(j)*onevolume
   end do

   reclng(1) = 1.d0/sqrt( recip(1,1)*recip(1,1) + &
         recip(2,1)*recip(2,1) + &
         recip(3,1)*recip(3,1) )
   reclng(2) = 1.d0/sqrt( recip(1,2)*recip(1,2) + &
         recip(2,2)*recip(2,2) + &
         recip(3,2)*recip(3,2) )
   reclng(3) = 1.d0/sqrt( recip(1,3)*recip(1,3) + &
         recip(2,3)*recip(2,3) + &
         recip(3,3)*recip(3,3) )

   ! interfacial distances given by dot of direct,recip
   ! sphere is radius of largest sphere inscribed in unit cell
   ! the minimum image cutoff must be less than or equal to this

   sphere = a+b+c
   do i = 1,3
      call dot(recip(1,i),ucell(1,i),result)
      distance = result*reclng(i)
      if ( distance < sphere )sphere = distance
   end do
   sphere = 0.5d0*sphere

   write(6, '(a,f9.3)') &
    '| Largest sphere to fit in unit cell has radius = ', sphere
   write(6, '(a)')  '# recip matrix, for conversion to fractional coords:'
   write(6, '(a,f14.7,a,f14.7,a,f14.7,a)') &
            '$r11 =',recip(1,1),'; $r12 =',recip(1,2),'; $r13 =',recip(1,3),';'
   write(6, '(a,f14.7,a,f14.7,a,f14.7,a)') &
            '$r21 =',recip(2,1),'; $r22 =',recip(2,2),'; $r23 =',recip(2,3),';'
   write(6, '(a,f14.7,a,f14.7,a,f14.7,a)') &
            '$r31 =',recip(3,1),'; $r32 =',recip(3,2),'; $r33 =',recip(3,3),';'
   write(6, '(a)')  '# ucell matrix, for conversion from fractional coords:'
   write(6, '(a,f14.7,a,f14.7,a,f14.7,a)') &
            '$u11 =',ucell(1,1),'; $u12 =',ucell(1,2),'; $u13 =',ucell(1,3),';'
   write(6, '(a,f14.7,a,f14.7,a,f14.7,a)') &
            '$u21 =',ucell(2,1),'; $u22 =',ucell(2,2),'; $u23 =',ucell(2,3),';'
   write(6, '(a,f14.7,a,f14.7,a,f14.7,a)') &
            '$u31 =',ucell(3,1),'; $u32 =',ucell(3,2),'; $u33 =',ucell(3,3),';'

end program get_ucell 

!     --- DOT ---

subroutine dot(v1,v2,result)
   _REAL_ v1(3),v2(3),result
   result = v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
   return
end subroutine dot 

!     --- CROSS ---

subroutine cross(v1,v2,v12)
   
   !    v12 is cross product of v1 and v2
   
   _REAL_ v1(3),v2(3),v12(3)
   v12(1) = v1(2)*v2(3)-v1(3)*v2(2)
   v12(2) = v1(3)*v2(1)-v1(1)*v2(3)
   v12(3) = v1(1)*v2(2)-v1(2)*v2(1)
   return
end subroutine cross 

