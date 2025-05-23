! <compile=optimized>
#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ compute penalty function for residual dipolar couplings
subroutine align1( natom, x, f, amass )

   use file_io_dat
   
   implicit none
   integer   natom
   _REAL_    x(*),f(*),amass(*)
#  include "nmr.h"
#  include "extra.h"
   character(len=1) uplo,jobz
   _REAL_    com(3)
   _REAL_    work(27), vect_al(3,3), almat(6), vect_ref(3,3)
   integer   ipear, iof, i, idip, j, iset, k
   _REAL_    dx, dy, dz, dr
   _REAL_    csx, csy, csz
   _REAL_    csx2, csy2, csz2
   _REAL_    djcst, dcalc
   _REAL_    ealigni
   _REAL_    temp1, temp2
   _REAL_    fx, fy, fz
   _REAL_    dtarget
   integer   ier
   _REAL_    tmass
   _REAL_    rootsum
   _REAL_    axlen
   _REAL_    news(3,3), Qnum, Qden
   _REAL_    root_tmp(3), vect_tmp(3,3)
   
   data uplo,jobz / 'U','V' /
   
   ipear = 0
   if (master .and. iprint /= 0) then
      rewind (57)
      write(57,9032)
      write(57,9031)
      write(57,9032)
   end if
   9031 format( &
         ' First atom        Last atom    curr. value target', &
         ' deviation  penalty  distance')
   9032 format(' ',78('-'))
   
   iof = 0
   do i=1,num_datasets
      s11(i)  = x(3*natom + iof + 1)
      s12(i)  = x(3*natom + iof + 2)
      s22(i)  = x(3*natom + iof + 3)
      s13(i)  = x(3*natom + iof + 4)
      s23(i)  = x(3*natom + iof + 5)
      s33(i)  = -s11(i) - s22(i)
      iof = iof + 5
   end do

   if( itarget .gt. 0 ) then

      !  ---first, diagonalize the input alignment tensor(s):
   
      do iset=1,num_datasets
         almat(1) = s11(iset)
         almat(2) = s12(iset)
         almat(3) = s22(iset)
         almat(4) = s13(iset)
         almat(5) = s23(iset)
         almat(6) = s33(iset)
         call D_OR_S()spev(jobz,uplo,3,almat,root(1,iset),vect_al,3,work,ier)

         if( itarget .eq. 1 ) then

            ! order the roots so by absolute value, so that |root(3)| >
            !    |root(2)| > |root(1)|:
            root_tmp = root(:,iset)
            vect_tmp = vect_al
            if( abs(root(1,iset)) > abs(root(3,iset)) ) then
               root(3,iset) = root_tmp(1)
               root(2,iset) = root_tmp(3)
               root(1,iset) = root_tmp(2)
               vect_al(:,3) = vect_tmp(:,1)
               vect_al(:,2) = vect_tmp(:,3)
               vect_al(:,1) = vect_tmp(:,2)
            else
               root(3,iset) = root_tmp(3)
               root(2,iset) = root_tmp(1)
               root(1,iset) = root_tmp(2)
               vect_al(:,3) = vect_tmp(:,3)
               vect_al(:,2) = vect_tmp(:,1)
               vect_al(:,1) = vect_tmp(:,2)
            end if

            ! construct the new eigenvalues from da_target and r_target:

            root(1,iset) = 0.5d0*da_target(iset) * (3.d0*r_target(iset) - 2.d0)
            root(2,iset) = -0.5d0*da_target(iset) * (3.d0*r_target(iset) + 2.d0)
            root(3,iset) =  2.d0 * da_target(iset)

         else if( itarget .eq. 2 ) then

            ! dac note: don't order the roots here to force an order of
            !   the eigenvalues: that can lead to a "jump" that reverse
            !   the sign of the tensor, and can mess up minimization or
            !   MD.  Below, we do order the roots for printing out of the
            !   final tensor.

            ! force the alignment tensor to have particular eigenvectors:
            ! uses vectors from dataset 1 for the other two:

            if( iset .eq. 1 ) then
                vect_ref = vect_al
            else
                vect_al = vect_ref
            end if

         end if

         ! reconstruct the alignment tensor:
         do i=1,3
            do j=1,3
               news(i,j) = 0.d0
               do k=1,3
                  news(i,j) = news(i,j) + vect_al(i,k)*root(k,iset)*vect_al(j,k)
               end do
            end do
         end do

         ! put news back into s:
         s11(iset) = news(1,1)
         s12(iset) = news(1,2)
         s22(iset) = news(2,2)
         s13(iset) = news(1,3)
         s23(iset) = news(2,3)
         s33(iset) = news(3,3)
      end do

      ! also, put back in the end of the coordinate array:
      iof = 0
      do i=1,num_datasets
         x(3*natom + iof + 1) = s11(i)
         x(3*natom + iof + 2) = s12(i)
         x(3*natom + iof + 3) = s22(i)
         x(3*natom + iof + 4) = s13(i)
         x(3*natom + iof + 5) = s23(i)
         iof = iof + 5
      end do

   end if  ! itarget .gt. 0
   
   !    Loop over obvserved dipolar couplings:
   
   Qnum = 0.0
   Qden = 0.0
   do idip=1,ndip
      
      i = id(idip)
      j = jd(idip)
      iset = dataset(idip)
      iof = 5*(iset-1)
      !========================================================================
      
      !    DJcst == -10**(-5)*gammai*gammaj*h/(2*pi*pi*r**3), in Hz
      
      !          = -0.077004 * gi * gj / r**3,  with r in Ang.
      
      !             Here, "gamma" is the nuclear gyromagneto ratio, and "g"
      !             is the nuclear "g" value (dimensionless).  The connection
      !             is gamma = g * betaN / hbar, where "betaN" is the "nuclear
      !             Bohr magneton" and hbar is Planck's constant/2*pi.
      
      !             Some values:  g( 1H) = 5.5856
      !                           g(13C) = 1.4048
      !                           g(15N) =-0.5663   note: negative!
      !                           g(31P) = 2.2632
      
      !              Examples:
      
      !                    DJcst=  0.2295      for N-H with rNH=1.02 A
      !                    DJcst=  0.2165      for N-H with rNH=1.04 A
      !                    DJcst= -0.4797      for C-H with rCH=1.08 A
      !                    DJcst= -0.4666      for C-H with rCH=1.09 A
      
      !   In earlier versions, DJcst was directly entered, assuming that
      !   the distance was constant.  Now we enter the product of the
      !   nuclear g-factors as gigj(idip) and the distance in Angstroms as
      !   dij(idip).  This allows us to also include gradients with respect
      !   to the interatomic distance.  If the value entered for dij is zero,
      !   then the actual distance is used, and derivatives with respect to
      !   this distance are included in the refinement.
      
      !   The anisotropy of the alignment tensor is defined by
      !   S11,  S22, S33, S12, S23, S13, in the  coordinate
      !   system X={1,0,0}, Y={0,1,0}, Z={0,0,1}
      
      !   In order to have the order of magnitude of the S values be
      !   commensurate with coordinates, we will mutliply them by 10**5.
      
      !   dcalc: calculated resudual D-D splitting for i-j spin pair.
      !   dobsu: upper bound for the observed resudual D-D splitting for
      !          i-j spin pair; if calc. value is larger than this, a
      !          penalty will be imposed
      !   dobsl: lower bound for the observed resudual D-D splitting for
      !          i-j spin pair; if calc. value is smaller than this, a
      !          penalty will be imposed
      !   (No penalty if dobsl < dobs < dobsu).
      !   ealign: Resudual D-D splitting restrant energy term is
      !        EJ=dwt*(DJ-OJ)**2
      
      !========================================================================
      
      dx = x(3*i-2)-x(3*j-2)
      dy = x(3*i-1)-x(3*j-1)
      dz = x(3*i  )-x(3*j  )
      dr = sqrt(dx**2+dy**2+dz**2)
      
      csx = dx/dr
      csy = dy/dr
      csz = dz/dr
      csx2 = csx*csx
      csy2 = csy*csy
      csz2 = csz*csz
      
      if( dij(idip) == 0.d0 ) then
         djcst = -0.077004d0 * gigj(idip) / (dr*dr*dr)
      else
         djcst = -0.077004d0 * gigj(idip) / dij(idip)**3
      end if

      dcalc = djcst*( s11(iset)*csx2 + &
            s22(iset)*csy2 + &
            s33(iset)*csz2 + &
            2.d0*s12(iset)*csx*csy + &
            2.d0*s13(iset)*csx*csz + &
            2.d0*s23(iset)*csy*csz  )

      if( dcalc > dobsu(idip) ) then
         ealigni = dwt(idip)*(dcalc - dobsu(idip))**2
      else if( dcalc < dobsl(idip) ) then
         ealigni = dwt(idip)*(dcalc - dobsl(idip))**2
      else
         ealigni = 0.d0
      end if
      ealign = ealign + ealigni
      
      !---   Derivatives with respect to Sij, assuming dij is fixed:
      
      if( dcalc > dobsu(idip) ) then
         temp1=2.0*dwt(idip)*(dcalc-dobsu(idip))*djcst
      else if( dcalc < dobsl(idip) ) then
         temp1=2.0*dwt(idip)*(dcalc-dobsl(idip))*djcst
      else
         temp1 = 0.d0
      end if
      
      if( ifreezes .eq. 0 ) then
         f(3*natom+iof+1) = f(3*natom+iof+1) - temp1*(csx2 - csz2)
         f(3*natom+iof+3) = f(3*natom+iof+3) - temp1*(csy2 - csz2)
      
         f(3*natom+iof+2) = f(3*natom+iof+2) - 2.d0*temp1*csx*csy
         f(3*natom+iof+4) = f(3*natom+iof+4) - 2.d0*temp1*csx*csz
         f(3*natom+iof+5) = f(3*natom+iof+5) - 2.d0*temp1*csy*csz
      end if
      
      !---   Derivatives with respect to x,y,z
      
      temp2 = -2.d0*temp1/dr
      fx = temp2*(  s11(iset)*csx*(csx2 - 1.d0) &
            + s22(iset)*csy2*csx &
            + s33(iset)*csz2*csx &
            + s12(iset)*csy*(2.d0*csx2 - 1.d0) &
            + s13(iset)*csz*(2.d0*csx2 - 1.d0) &
            + s23(iset)*2.d0*csx*csy*csz        )
      
      fy = temp2*(  s11(iset)*csx2*csy &
            + s22(iset)*csy*(csy2 - 1.d0) &
            + s33(iset)*csz2*csy &
            + s12(iset)*csx*(2.d0*csy2 - 1.d0) &
            + s13(iset)*2.d0*csx*csy*csz &
            + s23(iset)*csz*(2.d0*csy2 - 1.d0)  )
      
      fz = temp2*(  s11(iset)*csx2*csz &
            + s22(iset)*csy2*csz &
            + s33(iset)*csz*(csz2 - 1.d0) &
            + s12(iset)*2.d0*csx*csy*csz &
            + s13(iset)*csx*(2.d0*csz2 - 1.d0) &
            + s23(iset)*csy*(2.d0*csz2 - 1.d0)  )
      
      !---   Now, if desired, we need to add in the derivatives with respect
      !      to the (1/r**3) part of the dipolar coupling formula:
      
      if( dij(idip) == 0.d0 ) then
         if( dcalc > dobsu(idip) ) then
            temp1=2.0*dwt(idip)*(dcalc-dobsu(idip))*dcalc
         else if( dcalc < dobsl(idip) ) then
            temp1=2.0*dwt(idip)*(dcalc-dobsl(idip))*dcalc
         else
            temp1 = 0.d0
         end if
         temp2 = -3.d0*temp1/(dr*dr)
         fx = fx + temp2*dx
         fy = fy + temp2*dy
         fz = fz + temp2*dz
      end if
      
      f(3*i-2) = f(3*i-2) - fx
      f(3*i-1) = f(3*i-1) - fy
      f(3*i  ) = f(3*i  ) - fz
      f(3*j-2) = f(3*j-2) + fx
      f(3*j-1) = f(3*j-1) + fy
      f(3*j  ) = f(3*j  ) + fz
      
      if ( master .and. iprint /= 0 ) then
         if ( ealigni > dcut ) then
            dtarget = 0.5d0*(dobsu(idip)+dobsl(idip))
            if( dij(idip) == 0.d0 ) then
               write(57,9073) resat(i)(1:13),resat(j)(1:13), &
                     dcalc,dtarget,dcalc-dtarget,ealigni,dr
            else
               write(57,9073) resat(i)(1:13),resat(j)(1:13), &
                     dcalc,dtarget,dcalc-dtarget,ealigni,dij(idip)
            end if
         end if
      end if
      9073 format(' ',a13,' -- ',a13,':',5f9.3)

      ! stats for Qfactor:
      Qnum = Qnum + (dcalc - dtarget)**2
      Qden = Qden + dtarget**2
      
   end do  !  idip=1,ndip
   
   !     ---printing of results:
   
   if (master .and. iprint /= 0) then
      write(57,44) ealign, sqrt(Qnum/Qden)
      44 format(20x,'Total align    constraint:',f8.2, ' Q = ', f10.5)
      
      !       ---diagonalize the alignment tensor:
      
      do iset=1,num_datasets
         write(57,'(a,i3)') 'alignment tensor', iset
         write(57, &
        '(a,i2,a,f8.3,a,i2,a,f8.3,a,i2,a,f8.3,a,/,a,i2,a,f8.3,a,i2,a,f8.3,a)') & 
               's11(',iset,')=',s11(iset), &
             ', s12(',iset,')=',s12(iset), &
             ', s13(',iset,')=',s13(iset), ', ', &
             '  s22(',iset,')=',s22(iset), &
             ', s23(',iset,')=',s23(iset),', '
         almat(1) = s11(iset)
         almat(2) = s12(iset)
         almat(3) = s22(iset)
         almat(4) = s13(iset)
         almat(5) = s23(iset)
         almat(6) = s33(iset)
         call D_OR_S()spev(jobz,uplo,3,almat,root(1,iset),vect_al,3,work,ier)

         ! order the roots so by absolute value, so that |root(3)| >
         !    |root(2)| > |root(1)|:
         root_tmp = root(:,iset)
         vect_tmp = vect_al
         if( abs(root(1,iset)) > abs(root(3,iset)) ) then
            root(3,iset) = root_tmp(1)
            root(2,iset) = root_tmp(3)
            root(1,iset) = root_tmp(2)
            vect_al(:,3) = vect_tmp(:,1)
            vect_al(:,2) = vect_tmp(:,3)
            vect_al(:,1) = vect_tmp(:,2)
         else
            root(3,iset) = root_tmp(3)
            root(2,iset) = root_tmp(1)
            root(1,iset) = root_tmp(2)
            vect_al(:,3) = vect_tmp(:,3)
            vect_al(:,2) = vect_tmp(:,1)
            vect_al(:,1) = vect_tmp(:,2)
         end if

         write(57,*) 'Diagonalization:'
         do i=1,3
            write(57,'(f15.5,5x,3f12.5)') root(i,iset), (vect_al(j,i),j=1,3)
         end do
         write(57,'(a,f10.4,a,f10.4)') ' Da = ', root(3,iset)/2.d0, &
            ' x 10^-5;  R =', &
            2.d0*(root(1,iset)-root(2,iset))/(3.d0*root(3,iset))
         write(57,*)
      end do

   end if  ! (master .and. iprint /= 0)
   
   return
end subroutine align1 

!----------------------------------------------------------------------


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine alignread here]
subroutine alignread(natom,x)

   use file_io_dat

   implicit none
   integer   natom
   _REAL_    x(*)
#  include "nmr.h"
#  include "../include/md.h"
   character(len=80) line
   logical freezemol,freezes
   integer   iin, i, iof
   namelist /align/ ndip, id, jd, dobsu, dobsl, dataset, &
         num_datasets, s11,s22,s12,s13,s23,dcut,gigj,dij,dwt, &
         freezemol,freezes,da_target,r_target,itarget
   
   ! If restraint input has been redirected, open the appropriate file
   
   call amopen(37,redir(8)(1:iredir(8)),'O','F','R')
   iin = 37
   write(6,10) redir(8)(1:iredir(8))
   10 format(' Alignment info will be read from file: ',a)
   
   !  --- read and echo title from alignment file:
   
   write(6,*) 'Here are comments from the alignment input file:'
   do i=1,20
      read(iin,'(a)') line
      write(6,*) line(1:78)
   end do
   rewind (iin)
   write(6,*)
   
   ! read the namelist align, first setting up defaults:
   
   dcut = 0.1
   num_datasets = 1
   freezemol = .false.
   freezes = .false.
   itarget = 0
   do i=1,maxdip
      dwt(i) = 1.0d0
      dataset(i) = 1
      dij(i) = 0.d0
      dobsu(i) = 0.d0
      dobsl(i) = 0.d0
   end do
   read(iin,nml=align,err=30)
   if( ndip > maxdip) then
      write(6,*) 'ndip is too big: ',ndip,maxdip
      call mexit(6,1)
   end if
   if( num_datasets > maxdipsets) then
      write(6,*) 'num_datasets is too big: ',num_datasets,maxdipsets
      call mexit(6,1)
   end if
   ifreeze = 0
   if( freezemol ) ifreeze = 1
   ifreezes = 0
   if( freezes ) ifreezes = 1

   if( itarget .eq. 1 ) then
      do i=1,num_datasets
         if( (r_target(i) .lt. 0.0) .or. (r_target(i) .gt. 0.667) ) then
            write(6,*) 'Error: r_target must be between 0 and 2/3: ', i, &
                r_target(i)
            call mexit(6,1)
         end if
      end do
   end if
   
   iof = 0
   do i=1,num_datasets
      x(3*natom + iof + 1) = s11(i)
      x(3*natom + iof + 2) = s12(i)
      x(3*natom + iof + 3) = s22(i)
      x(3*natom + iof + 4) = s13(i)
      x(3*natom + iof + 5) = s23(i)
      iof = iof + 5
   end do
   iscale = 5*num_datasets
   
   return
   30 write(6,*) 'namelist reports error reading &align'
   call mexit(6,1)
end subroutine alignread 
