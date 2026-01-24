!<compile=optimized>

#include "../include/dprec.fh"

!> Utility functions and subroutines used throughout the RISM code
!! base.
module rism_util
  use rism_report_c
  ! use, intrinsic :: iso_c_binding

  interface r2c_pointer
     module procedure r2c_pointer_1d
     module procedure r2c_pointer_2d
  end interface r2c_pointer

  public r2c_pointer
  
contains

  !> use to remove any overall force component
  subroutine corr_drift(ff,mass,numAtoms &
#ifdef MPI
         ,rank,size,comm &
#endif /*MPI*/
         )
#ifdef MPI
    use mpi
#endif
      implicit none
#if MPI
      integer, intent(in) :: rank,size,comm
#endif /*MPI*/
      integer, intent(in) :: numAtoms
      _REAL_,intent(inout) :: ff(3,numAtoms)
      _REAL_,intent(in) :: mass(numAtoms)
      
      integer :: id,iatu
      _REAL_ :: totmass,totfrc(3)
      integer :: err
#ifdef RISM_DEBUG
      write(6,*) "CORR_DRIFT"
#endif /*RISM_DEBUG*/
      do id = 1, 3
         totfrc(id) = sum(ff(id,1:numAtoms))
#ifdef RISM_DEBUG
         write (6,*) "ID=",id
         write (6,*) ff(id,1:numAtoms)
#endif /*RISM_DEBUG*/
      end do
#if defined(MPI)
      call MPI_ALLREDUCE(MPI_IN_PLACE,totfrc,3,MPI_DOUBLE_PRECISION,&
           MPI_SUM,comm,err)
      if (err /=0) call rism_report_error&
           ("RISM3D CORR_DRIFT: could not reduce TOTFRC")
#endif /*defined(MPI)*/

#ifdef RISM_DEBUG
      write(6,*)"checking drift..."
      write(6,*)"x",totfrc(1)
      write(6,*)"y",totfrc(2)
      write(6,*)"z",totfrc(3)
      write(6,*)"corrected drift..."
#endif /*RISM_DEBUG*/

      
      ! The force is distributed across all so the correction must be
      ! too.  That is, the correction is divided by the number of
      ! processors.
      
      totmass = sum(mass)
      do iatu=1,numAtoms
#if defined(MPI)
         ff(:,iatu) = ff(:,iatu) - totfrc*mass(iatu)/totmass/size
#else
         ff(:,iatu) = ff(:,iatu) - totfrc*mass(iatu)/totmass
#endif /*defined(MPI)*/
      end do
#ifdef RISM_DEBUG
      do id = 1, 3
         totfrc(id) = sum(ff(id,1:numAtoms))
      end do      
#  if defined(MPI)
      call MPI_ALLREDUCE(MPI_IN_PLACE,totfrc,3,MPI_DOUBLE_PRECISION,&
           MPI_SUM,comm,err)
      if (err /=0) call rism_report_error&
           ("RISM3D CORR_DRIFT: could not reduce TOTFRC")
#  endif /*defined(MPI)*/
      write(6,*)"x",totfrc(1)
      write(6,*)"y",totfrc(2)
      write(6,*)"z",totfrc(3)
      write(6,*)'outputting ff to file:  ala.ff'
      !#endif
#endif /*RISM_DEBUG*/
    end subroutine corr_drift

!> Calculates the vector cross product of A and B and places it in C.
!! It is assumed that A, B and C are of length three.
!! @param[in] a Three element array.
!! @param[in] b Three element array.
!! @param[in] c Three element array, on return contains the cross product of A and B.
!! @return Vector cross product of a and b.
function cross(a, b) result(c)
  implicit none
  _REAL_, intent(in) :: a(3), b(3)
  _REAL_ :: c(3)
  c(1) =  a(2) * b(3) - a(3) * b(2)
  c(2) = -a(1) * b(3) + a(3) * b(1)
  c(3) =  a(1) * b(2) - a(2) * b(1)  
end function cross

!> Calculates the magnitude of a 3D vector.
!! @param[in] a Three element array.
!! @return Magnitude of a.
function magnitude(a)
  implicit none
  _REAL_, intent(in) :: a(3)
  _REAL_ :: magnitude
  magnitude = sqrt(abs(dot_product(a, a)))
end function magnitude

!! Tests if number is prime. Note: this is a slow algorithm.  It checks if 2 or
!! any odd number (except 1) less than the square root of number is a factor.
!! IN:
!!    number : number to test
!! OUT:
!!     .true. if prime, otherwise .false.
function isPrime(number)
  implicit none
  logical :: isprime
  integer :: i,number,j
  isprime = .false.
  i = abs(number)
  if (mod(i,2)==0) return 
  j=3
  do while (j**2 <= i)
     if (mod(i,j)==0) return 
     j = j+2
  end do
  isprime = .true.
end function isPrime


!! Checks if the the number is factorizable by the list of numbers given. 
!! Works best for prime numbers.  It is also best to pre-sort the list from 
!! highest to lowest.
!! IN:
!!    number : number to test
!!    factors: potential factors
!! OUT:
!!     .true. if factorizable, otherwise .false.
function isFactorable(number, factor)
  implicit none
  integer, intent(in) :: number, factor(:)
  logical :: isFactorable
  integer :: i, num, numold
  isFactorable = .false.
  num = number
  do i = 1, ubound(factor, 1)
     numold = num + 1
     do while (numold > num)
        numold = num
        if (mod(num, factor(i)) == 0) then
           num = num / factor(i)
        end if
     end do
  end do
  if (num == 1) isFactorable=.true.
end function isFactorable

!> Find the largest prime factor of a number.
!! @param[in] number Number to test.
!! @return Largest prime factor of number.
function largestPrimeFactor(number)
  implicit none
  integer, intent(in) :: number
  integer :: i, largestPrimeFactor
  largestPrimeFactor = 1
  do i = 2, number
     if (isPrime(i) .and. mod(number, i) == 0) then
        largestPrimeFactor = i
     end if
  end do
end function largestPrimeFactor

!! Fills PTR with indices of VAL to give the values of VAL in accending order.
!! I.e.  VAL(PTR) will be in accending order.  The reproduces the functionality
!! of INDEXX from Numerical Recipes.
!! IN:
!!    val :: array of values
!!    ptr :: will be an array of VAL indices.  I.e. it is a 'pointer' to the VAL
!!           elements
!!    n   :: length of the arrays

subroutine indexArray(val, ptr, n)
  implicit none
  integer, intent(in) :: n
  _REAL_, intent(in) :: val(n)
  integer, intent(out) :: ptr(n)
  integer :: i, temp

  !Following Numerical Recipes, we are using heapsort to sort the array in place.
  !For the actual heapsort algorithm we have follow Wikipedia heapsort article
  !from 2010/02/13.

  !initialize ptr
  ptr =  (/(i, i=1,n)/)

  !build the heap
  do i =  (n+1)/2, 1, -1
     call siftDown(val,ptr,i,n,n)
  end do

  !Traverse the heap to produce the sorted order
  do i = n, 1, -1
     !move the largest element (found in the first index) to the last index
     temp = ptr(i)
     ptr(i) = ptr(1)
     ptr(1) = temp
     !bring up the next largest element to the first index
     call siftDown(val(1:i-1),ptr(1:i-1),1,i-1,n)
  end do
end subroutine indexArray


!! Used in the indexArray heapsort.  Brings the largest value between I and M
!! to index I while ensuring that values at 2*I and 2*I+1 are less than the value
!! at I.
!! IN:
!!    val :: array of values
!!    ptr :: will be an array of VAL indices.  I.e. it is a 'pointer' to the VAL
!!           elements
!!    i   :: the start element to use for the arrays
!!    m   :: the end element to use for the arrays
!!    n   :: length of the arrays

subroutine siftDown(val, ptr, i, n, m)
  implicit none
  integer, intent(in) :: i,n,m
  _REAL_, intent(in) :: val(m)
  integer, intent(inout) :: ptr(m)
  integer:: root,child
  integer :: temp
  root=i
  do while (root*2 <= n)
     child = root*2
     if (child+1 <= n) then
        if (val(ptr(child+1))> val(ptr(child))) child = child+1
     end if
     if (child <= n) then
        if (val(ptr(child))> val(ptr(root))) then
           temp = ptr(child) 
           ptr(child) = ptr(root)
           ptr(root) = temp
           root = child
        else 
           exit
        end if
     end if
  end do
end subroutine siftDown


!! 'Progressive' polynomial interpolation using Neville's algorithm. Adds
!! one point at a time, starting from the lowest indicies, to the
!! interpolation until either the data is exhausted or the relative
!! difference begins to increase.  Can be used as a drop in
!! replacement for polynomialInterpolation.
!! IN:
!!    xa :: array of x values
!!    ya :: array of y values
!!    n  :: number of array elements, must be >2
!!    x  :: argument to the polynomial
!!    y  :: (out) polynomial value at x
!!    dy :: (out) error estimate in y
subroutine polynomialInterpolation_progressive(xa, ya,n, x, y, dy)
  implicit none
  integer, intent(in) :: n
  _REAL_, intent(in) :: xa(n),ya(n),x
  _REAL_, intent(out) :: y, dy
  _REAL_ :: P(n),rel_diff,rel_diff0, y0
  integer :: i,m
  rel_diff=huge(1d0)
  rel_diff0=huge(1d0)
  y0=huge(1d0)
  P=ya
  do m = 1,n-2
     !Neville's algorithm
     do i = 1, n-m
        p(i) = (x-xa(i+m))*p(i) - (x-xa(i))*p(i+1)
        p(i) = p(i)/(xa(i)-xa(i+m))
     end do
     y=((x-xa(n))*p(1) - (x-xa(1))*p(2))/(xa(1)-xa(n))

     !compute relative difference between this and the previous iteration
     if (y0 /= huge(1d0)) then
        rel_diff = abs((y0-y)/y)
     end if
     !estimated error
     dy = (p(1)+p(2)-2d0*y)/2d0
     !if the relative difference is not decreasing, we're done
     if (rel_diff0 < rel_diff) then
        y=y0
        return
     end if
     !update old values
     y0=y
     rel_diff0=rel_diff
  end do
  !if we get here, we are just returning the y from using all points
end subroutine polynomialInterpolation_progressive


!! Polynomial interpolation using Neville's algorithm.
!! Serves as a drop in replacement for Numerical Recipes POLINT subroutine; 
!! however, the error estimate is done differently but is of the same order.
!! IN:
!!    xa :: array of x values
!!    ya :: array of y values
!!    n  :: number of array elements
!!    x  :: argument to the polynomial
!!    y  :: polynomial value at x
!!    error :: Error estimate in y.
subroutine polynomialInterpolation(xa, ya, n, x, y, error)
  implicit none
  integer, intent(in) :: n
  _REAL_, intent(in) :: xa(n), ya(n), x
  _REAL_, intent(out) :: y, error
  _REAL_ :: P(n)
  integer :: i, m
  P = ya
  do m = 1, n - 2
     do i = 1, n - m
        p(i) = (x - xa(i + m)) * p(i) - (x - xa(i)) * p(i + 1)
        p(i) = p(i) / (xa(i) - xa(i + m))
     end do
  end do
  y = ((x - xa(n)) * p(1) - (x - xa(1)) * p(2)) / (xa(1) - xa(n))
  error = (p(1) + p(2) - 2d0 * y) / 2d0
end subroutine polynomialInterpolation

!! Does a MPI sum

function checksum(a,n,comm)
#ifdef MPI
    use mpi
#endif
  implicit none
  integer, intent(in) :: n,comm
  _REAL_, intent(in) :: a(n)
  _REAL_ :: checksum
  integer :: err
  checksum = sum(a)
#ifdef MPI
    call mpi_allreduce(MPI_IN_PLACE,checksum,1,MPI_DOUBLE_PRECISION,mpi_sum,comm,err)
    if (err /=0) call rism_report_error&
         ("RISM3D CHECKSUM: could not reduce CHECKSUM")
#endif /*MPI*/
  
end function checksum


!> Least common multiple of two integers.
!! @param[in] a Integer.
!! @param[in] b Integer.
!! @return Least common multiple of a and b.
function lcm(a, b)
  implicit none
  integer, intent(in) :: a, b
  integer :: lcm
  lcm = a * b / gcd(a,b)
end function lcm


!> Greatest common divisor of two integers.
!! @param[in] s Integer.
!! @param[in] t Integer.
!! @return Greatest common divisor of s and t.
recursive function gcd(s,t) result(res)
  implicit none
  integer, intent(in) :: s, t
  integer :: a, b, c, res
  a = s
  b = t
  if (b > a) then
     c = a
     a = b
     b = c
  end if
  if (b == 0) then
     res = a
     return
  end if
  a = mod(a, b)
  res = gcd(b, a)
end function gcd


!! Convert string to upper case
!! IN:
!!    cc : string will be converted to upper case in place
subroutine  caseup (cc)
  implicit none
  character(len=*)  cc
  integer :: i,bi
  do i=1,len(cc)
     bi = iachar(cc(i:i))
     if (bi >= iachar('a') .AND. bi <= iachar('z')) &
          cc(i:i) = achar(iachar(cc(i:i)) - 32 )
  end do
  return
end subroutine caseup


!! Convert string to lower case
!! IN:
!!    cc : string will be converted to lower case in place
subroutine  caselow (cc)
  implicit none
  character(len=*)  cc
  integer :: i,bi
  do i=1,len(cc)
     bi = iachar(cc(i:i))
     if (bi >= iachar('A') .AND. bi <= iachar('Z')) &
          cc(i:i) = achar(iachar(cc(i:i)) + 32 )
  end do
  return
end subroutine caselow


!> Find a free Fortran unit for temporary use.
function freeUnit(start) result(unit)
  implicit none
  integer, optional, intent(in) :: start !< Minimum unit number that is acceptible.
  integer :: unit
  logical :: opened
  opened = .true.
  ! Skip commonly used unit numbers.
  unit = 10
  ! User defined minimum unit number.
  if (present(start)) unit = start
  ! Search for free unit.
  unit = unit - 1
  do while (opened)
     unit = unit + 1
     inquire(unit, opened = opened)
  end do
end function freeUnit


!> Intel CPUs typically use extended (80-bit) precision when possible.
!! However, different compilers may or may not use the extended
!! precision value when printing.  To ensure consistency, we test each
!! value and values less than TINY() are printed as zero.  These very
!! small values are actually meaningless and break our testing
!! procedures.
elemental function rmExPrec(ep) result(rp)
  implicit none
  _REAL_, intent(in) :: ep
  _REAL_ :: rp
  if (abs(ep) < tiny(ep)) then
     rp = 0d0
  else
     rp = ep
  end if
end function rmExPrec

!> Approximate rounding a Fortran floating-point number to a certain
!! number of decimals.
!! This is useful for getting rid of insignificant digits in code
!! sensitive to the exact value of the float. Unfortunately Fortran
!! does not have a built-in utility for this.
function round(a, decimal)
  implicit none
  _REAL_, intent(in) :: a !< Fortran floating-point number.
  integer, intent(in) :: decimal !< Decimal place to round to.
  _REAL_ :: round

  !TODO: Is there a better, perhaps more robust way to achieve this?
  round = anint(a * 10**decimal) / 10**decimal
end function round


!> Takes a _REAL_ pointer and returns a complex pointer of the same
!! size with the leading dimension cut in half. This is the memory
!! layout required for FFTW. This uses Fortan 2003 features.
function r2c_pointer_1d(rdata) result (cdata)
  use, intrinsic :: iso_c_binding
  implicit none
  _REAL_, pointer, intent(in) :: rdata(:) !< real data (ndata).
  complex(kind(1d0)), pointer :: cdata(:) !< pointer to complex data
                                          !! (ndata/2).
  type(c_ptr) :: cptr = c_null_ptr
  nullify(cdata)
  cptr = c_loc(rdata(1))
  call c_f_pointer(cptr, cdata, [ubound(rdata, 1) / 2])
end function r2c_pointer_1d


!> Takes a _REAL_ pointer and returns a complex pointer of the same
!! size with the leading dimension cut in half. This is the memory
!! layout required for FFTW. This uses Fortan 2003 features.
function r2c_pointer_2d(rdata) result (cdata)
  use, intrinsic :: iso_c_binding
  implicit none
  _REAL_, pointer, intent(in) :: rdata(:,:) !< real data (ndata, narray).
  complex(kind(1d0)), pointer :: cdata(:,:) !< pointer to complex data
                                            !! (ndata/2, narray).
  type(c_ptr) :: cptr = c_null_ptr
  nullify(cdata)
  cptr = c_loc(rdata(1, 1))
  call c_f_pointer(cptr, cdata, [ubound(rdata, 1)/2, ubound(rdata, 2)])
end function r2c_pointer_2d

!> Heaviside function with offset
!! @param x value to apply Heaviside to.  May be an array
!! @param offset offset from 0. for the transition point
elemental _REAL_ function heaviside(x,offset)
  implicit none
  _REAL_, intent(in)::x,offset
  if (x==offset) then
     heaviside=0.5d0
  elseif (x.gt.offset) then
     heaviside = 1.d0
  else
     heaviside = 0d0
  end if
end function heaviside
  
!! Computes zeroth order spherical Bessel function.  Adapted from the GNU 
!! Scientific library.
!! IN:
!!    x :: calculate the spherical Bessel function a x
!!    o_err :: (optional) absolute error in the result.  For abs(x)>0.5d0 this is
!!             an approximation of the GSL version since we don't know the error
!!             for the intrinsic sin(x)
!! OUT:
!!     the value of the zeroth spherical Bessel function at x

function spherical_bessel_j0(x,o_err) result(val)
  implicit none
  _REAL_, intent(in) :: x
  _REAL_, optional, intent(out) :: o_err
  _REAL_ :: val
  _REAL_ :: ax
  _REAL_, parameter :: c1 = -1d0/6d0, &
       c2 =  1d0/120d0,&
       c3 = -1d0/5040d0,&
       c4 =  1d0/362880d0,&
       c5 = -1d0/39916800d0,&
       c6 =  1d0/6227020800d0
  ax = abs(x)
  if (ax < 0.5d0) then
     ax = ax*ax
     val = 1d0 + ax*(c1 + ax*(c2 + ax*(c3 + ax*(c4 + ax*(c5 + ax*c6)))))
     if (present(o_err)) o_err = 2d0 * epsilon(1d0) * abs(val)
  else
     val  = sin(x)/x;
     if (present(o_err)) then
        !we don't know what the error is for intrinsic sin(x)!
        o_err  = abs(epsilon(1d0)/x)  + 2.0 * epsilon(1d0) * abs(val)
     end if
  end if
end function spherical_bessel_j0


!! Computes first order spherical Bessel function.  Adapted from the GNU 
!! Scientific library.
!! IN:
!!    x :: calculate the spherical Bessel function a x
!!    o_err :: (optional) absolute error in the result.  For abs(x)>0.25d0 this is
!!             an approximation of the GSL version since we don't know the error
!!             for the intrinsics sin(x) and cos(x)
!! OUT:
!!     the value of the first spherical Bessel function at x

function spherical_bessel_j1(x,o_err) result(val)
  implicit none
  _REAL_, intent(in) :: x
  _REAL_, optional, intent(out) :: o_err
  _REAL_ :: val
  _REAL_ :: ax
  _REAL_, parameter :: c1 = -1d0/10d0,&
       c2 =  1d0/280d0,&
       c3 = -1d0/15120d0,&
       c4 =  1d0/1330560d0,&
       c5 = -1d0/172972800d0

  ax = abs(x)
  if (x==0) then
     val=0d0
     if (present(o_err)) o_err=0d0
  else if (ax < 3.1d0*tiny(1d0)) then
     !GSL  gives an underflow error here.  We'll just approximate it as zero
     val=0d0
     if (present(o_err)) o_err=0d0
  else if (ax < 0.25d0) then
     ax = ax*ax
     val = x/3d0*(1d0 + ax*(c1 + ax*(c2 + ax*(c3 + ax*(c4 + ax*c5)))))
     if (present(o_err)) o_err = 2d0 * epsilon(1d0) * abs(val)
  else
     val  = (sin(x)/x - cos(x))/x;
     if (present(o_err)) then
        !we don't know what the error is for intrinsic sin(x)!
        o_err  = (abs(epsilon(1d0)/x)+epsilon(1d0))/abs(x)  + 2.0 * epsilon(1d0) * abs(val)
        o_err = o_err+ 2.0 * epsilon(1d0)*(abs(sin(x)/(x*x)) + abs(cos(x)/x))
     end if
  end if
end function spherical_bessel_j1

! Calculate the principal axes of the atom distribution.  
! This follows, in part, mofi() in nmode/thermo.f.
!IN:
!   ratu   :: the x,y,z position of each solute atom. (3,natom)
!   mass   :: mass of each atom
!   pa     :: the three prinicpal axes
subroutine principalAxes(ratu,mass,pa)
  implicit none
  _REAL_,intent(in) :: ratu(:,:),mass(:)
  _REAL_, intent(out) :: pa(3,3)
  integer :: id,ier

  !t        : moment of inertia tensor in 1-d, in upper triangular, 
  !           column-major format (xx, xy, yy, xz, yz, zz)
  !eigenval : This will be the moment of interia in each direction
  !work     : temp space for the algorithm
  _REAL_ :: t(6),eigenval(3),eigenvec(3,3),work(3*3)
#ifdef RISM_DEBUG
  write(6,*)"CALC PA"
  call flush(6)
#endif /*RISM_DEBUG*/
  
  call dspev('V','U',3,momentOfInertia(ratu,mass),eigenval,pa,3,work,ier)

end subroutine principalAxes

! Calculate the moment of inertia in upper triangular form.
! Uses column-major format (xx, xy, yy, xz, yz, zz).
! This follows, in part, mofi() in nmode/thermo.f.
!IN:
!   ratu   :: the x,y,z position of each solute atom. (3,numAtoms)
!   mass   :: mass of each atom
function momentOfInertia(ratu,mass)
  implicit none
  _REAL_,intent(in) :: ratu(:,:),mass(:)

  !> Moment of inertia tensor in 1-d, in upper triangular form.
  _REAL_ :: momentOfInertia(6)
#ifdef RISM_DEBUG
  write(6,*)"CALC momentOfInertia"
  call flush(6)
#endif /*RISM_DEBUG*/
  momentOfInertia(1) = sum(mass * (ratu(2, :)**2 + ratu(3, :)**2))
  momentOfInertia(3) = sum(mass * (ratu(1, :)**2 + ratu(3, :)**2))
  momentOfInertia(6) = sum(mass * (ratu(1, :)**2 + ratu(2, :)**2))
  momentOfInertia(2) = -sum(mass * (ratu(1, :) * ratu(2, :)))
  momentOfInertia(4) = -sum(mass * (ratu(1, :) * ratu(3, :)))
  momentOfInertia(5) = -sum(mass * (ratu(2, :) * ratu(3, :)))
  
end function momentOfInertia

!> Assumes that the principal axes have been aligned to coincide with
!! the x-, y- and z- axes for both structures. The structures are then
!! compared to each other and the first rotated to best fit the second
!! while maintaining the orientation of the PA.
!! @param[in,out] ratu The x,y,z position of each solute atom.
!! @param[in,out] ratu2 The x,y,z position of each solute atom.
!! @param[in] numAtoms The number of solute atoms.
subroutine alignorient(ratu, ratu2, numAtoms, backquat)
  use constants, only : PI
  use quaternion, only : rotate_quat, quat_mult
  use rism_report_c
  implicit none
  integer,intent(in) :: numAtoms
  _REAL_, intent(inout) :: ratu(3, numAtoms)
  _REAL_, intent(in) :: ratu2(3, numAtoms)
  _REAL_, intent(out) :: backquat(4)
  _REAL_ :: rotaxis(3), quat(4), tempquat(4)
  integer :: axis(3)
  integer :: id,iatu

  ! Here we calculate the product of each coordinate of each atom in
  ! the two structures.  If both have had their PA aligned to the
  ! coordinate system there are zero to three 180 degree rotations
  ! that will align the two molecules.  The rotation axes are those
  ! with over all positive coordinate products (some atom movement may
  ! distort the molecule).

#ifdef RISM_DEBUG
  write(6,*) "ALIGNORIENT"
#endif /* RISM_DEBUG*/
  do id = 1, 3
     rotaxis(id) = sum(ratu(id, :) * ratu2(id, :))
     if (rotaxis(id) < 0) then
        axis(id) = 0
     else
        axis(id) = 1
     end if
  end do
#ifdef RISM_DEBUG
  write(6,*) "rotate axis", axis, rotaxis
#endif /*RISM_DEBUG*/
  if (sum(axis) == 0) then
     call rism_report_error("ALIGNORIENT failed")
  else if (sum(axis) == 3) then
     ! Already aligned.
     backquat = (/ 1d0, 0d0, 0d0, 0d0 /)
     return
  end if
  backquat = 0d0
  do id = 1, 3
     if (axis(id) == 1) then
        rotaxis = 0d0
        rotaxis(id) = 1d0

!! $        quat(2:4) = rotaxis * sin(PI / 2d0)
!! $        quat(1) = cos(PI / 2d0)
        quat = 0d0
        quat(id + 1) = 1d0

        do iatu = 1, numAtoms
           call rotate_quat(ratu(1:3, iatu), quat)
        end do
#ifdef RISM_DEBUG
        write(6,*) "QUAT", quat, backquat
#endif /*RISM_DEBUG*/
        if (sum(backquat) == 0d0) then
           backquat(1) = quat(1)
           backquat(2:4) = -1d0 * quat(2:4)
!           backquat(2:4) = quat(2:4)
        else
           quat(2:4) = -1d0 * quat(2:4)
           call quat_mult(quat, backquat, tempquat)
           backquat = tempquat
        end if
#ifdef RISM_DEBUG
        write(6,*) "QUAT", quat, backquat
#endif /*RISM_DEBUG*/
     end if
  end do
#ifdef RISM_DEBUG
  write(6,*) ratu(:,1)
#endif /*RISM_DEBUG*/
end subroutine alignorient


!Rotates the system such that the previously calculated pricipal axes
!coincide with the x-,y- and z-axes.
!IN:
!   ratu   :: the x,y,z position of each solute atom.  (3,natom) This is modified.
!   pa     :: the three prinicpal axes
subroutine orientToPrincipalAxes(ratu, pa, backquat)
  use constants, only : PI
  use quaternion, only : quaternionFromEulerAxis, rotate_quat, quat_mult
  implicit none
  _REAL_,intent(inout) :: ratu(:,:)
  _REAL_, intent(inout) :: pa(3,3)
  _REAL_, intent(out) :: backquat(4)
  _REAL_ :: angle, quat(4),dir(3),xaxis(3)=(/1d0,0d0,0d0/),yaxis(3)=(/0d0,1d0,0d0/),&
       checkv(3), tempquat(4)
  integer :: numAtoms

  integer :: iatu,ipa
#ifdef RISM_DEBUG
  write(6,*) "ORIENT_PA"
#endif /*RISM_DEBUG*/
  numAtoms=ubound(ratu,2)

  ! Rotate first principal axis the x-axis.

  ! Get the angle.  This is always positive.
  angle = acos(min(1d0,max(-1d0,dot_product(pa(1:3,1),xaxis))))
  if (angle < PI - 1d-6 .and. angle > -PI+1d-6) then
     dir = cross(pa(1:3,1), xaxis)
  else
     dir = yaxis
  end if
  
  ! Get the cross product between pa(:,1) and dir.  This will determine
  ! if we should rotate using +ive or -ive angle.
  checkv = cross(dir,pa(1:3,1))
  if (dot_product(checkv,xaxis) < 0) then
     angle = -angle
  end if
  quat = quaternionFromEulerAxis(angle,dir)
  do iatu = 1 ,numAtoms
     call rotate_quat(ratu(1:3,iatu),quat)
  end do
  do ipa = 1,3
     call rotate_quat(pa(1:3,ipa),quat)
  end do
  backquat = quat

  ! Next the second PA (this also places the third one as well)
  ! max/min protects against round-off errors in normalized vectors.
  angle = acos(min(1d0,max(-1d0,dot_product(pa(1:3,2),yaxis))))
  dir = xaxis
  checkv = cross(dir, pa(1:3,2))
  
  if (dot_product(checkv,yaxis) < 0) then
     angle = -angle
  end if

  quat = quaternionFromEulerAxis(angle,dir)
  do iatu = 1 ,numAtoms
     call rotate_quat(ratu(1:3,iatu),quat)
  end do
  do ipa = 1,3
     call rotate_quat(pa(1:3,ipa),quat)
  end do

  call quat_mult(quat,backquat,tempquat)
  backquat = tempquat
  backquat(2:4) = -1d0*backquat(2:4)

#ifdef RISM_DEBUG
  write(6,*) ratu(:,1)
#endif /*RISM_DEBUG*/

end subroutine orientToPrincipalAxes

end module rism_util
