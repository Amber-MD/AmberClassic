! <compile=optimized>
#include "../include/assert.fh"
#include "../include/dprec.fh"

module sgld

! MODULE: sgld
! ================== SELF-GUIDED MOLECULAR/LANGEVIN DYNAMICS ====================
! Xiongwu Wu, 2024

#ifdef MPI
   use mpi
#endif

implicit none

!     Head file for the self-guided Langevin Dynamics simulation  
!
!      variables for SGLD simulation
!


! ... integers:
!
!
!    SGMD/SGLD applying range
!      ISGSTA      Begining atom index applying SGLD
!      ISGEND      Ending atom index applying SGLD
!
!
      integer,save:: ISGSTA,ISGEND,NSGATOM
!
! ... floats:
!
!
!    SGMD/SGLD VARIABLES
!     SGFT    !  Guiding factor 
!     TSGAVG  !  Local average time, ps
!     TSGAVP  !  Convergence time, ps
!
!
!     SGAVG0  !  Local average remains
!     SGAVG1  !  Local average factor, SGAVG1=1-SGAVG0
!     SGAVP0  !  Convergence average remains
!     SGAVP1  !  Convergency average factor, SGAVP1=1-SGAVP0
!     GAMMAS  !  friction coefficient
!     TEMPLF  !  Low frequency motion temperature 
!     TEMPHG  !  High frequency temperature
!

      _REAL_  SGFT,SGFF,SGFG,TEMPSG,TSGAVG,TSGSET,TSGAVP,GAMMAS,SGGAMMA,SGSIZE, &
          SGAVG0,SGAVG1,SGAVP0,SGAVP1,SGFTI,SGFFI,SGFGI,EPOTLF,EPOTHF,EPOTLLF, &
          TEMPLF,TEMPHF, &
          SGWT, sgrndf, sgmsum, com0sg(3),com1sg(3),com2sg(3), &
          MYSCALSG,SGSCALE,FSGLDG,PSGLDG,tsgfac

!  common block for parallel broadcast
      common/sgldr/SGFT,SGFF,SGFG,TEMPSG,TSGAVG,TSGSET,TSGAVP,GAMMAS,SGGAMMA,SGSIZE, &
      SGAVG0,SGAVG1,SGAVP0,SGAVP1,SGFTI,SGFFI,SGFGI,FSGLDG,PSGLDG, &
      EPOTLF,EPOTHF,EPOTLLF,TEMPLF,TEMPHF, &
      SGWT, sgrndf, sgmsum, com0sg,com1sg,com2sg, &
      MYSCALSG,SGSCALE,tsgfac
!  Number of broadcasting variables
integer, parameter :: nsgld_real=39
       
!
! ... flags:
!
!
!     isgld   ! input control; positive values activate SGLD
!     TSGLDGLE ! Perform SGLD-GLE to mantain canonical ensemble distribution
!     trxsgld   !  replica exchange sgld (isgld>0 & rem>0)
!     tsgbond   !  local average over bonded atoms
!     tsgmap    !  local average over atoms within sgsize distance
!
!
      integer, save :: isgld,sgtype
      LOGICAL, save, private :: TSGLD,TSGLDGLE,TSGBOND,TSGMAP,tsggamma
      LOGICAL, save :: TRXSGLD
      character(8), save:: sglabel
      character(256), save:: sgmask
      
!*******************************************************************************
!
! ...allocatable arrays:
!
!     avgx1     ! local averages of position
!     avgx2     ! local averages of local averages of position
!     avgp     ! local averages of momentum
!     avgr     ! local average of random forces
!     sgfps   !  averages of force-momentum product
!     sgpps   ! averages of momentum-momentum product
!*******************************************************************************

      double precision, dimension(:,:), allocatable, save :: avgx0,avgx1,avgx2,avgr
      double precision, dimension(:), allocatable, private,save :: sgfps,sgpps,sgmass

      type listdata_rec
      integer             :: offset
      integer             :: cnt
      end type listdata_rec
  

!*******************************************************************************
!
! spatial average variables:
!
!     atm_sg_maskdata     ! SG local structure of each atom (count, offset)
!     atm_sg_mask         ! atom id in SG local structures
!     NGRIDX,NGRIDY,NGRIDZ     ! SG map dimensions in x, y, z directions
!     NGRIDXY,NGRIDYZ,NGRIDZX     ! SG map dimensions in x-y, y-z, z-x plan
!     NGRIDXYZ     ! SG map dimensions in x-y-z
!     GXMIN,GXMAX,GYMIN,GYMAX,GZMIN,GZMAX     ! SG map boundary
!     RHOM   !  mass map
!     RHOVX,RHOVY,RHOVZ   ! velocity vx, vy, vz maps
!     RHOAX,RHOAY,RHOAZ   ! acceleration ax, ay, az maps
!*******************************************************************************
      type(listdata_rec), allocatable, save :: atm_sg_maskdata(:)
      integer, allocatable, save            :: atm_sg_mask(:),sgatoms(:)

      integer,private,save:: NGRIDX,NGRIDY,NGRIDZ,NGRIDXY,NGRIDYZ,NGRIDZX,NGRIDXYZ
      double precision,private,save::  GXMIN,GXMAX,GYMIN,GYMAX,GZMIN,GZMAX
      double precision, dimension(:), allocatable, private,save :: RHOM,RHOVX,RHOVY,RHOVZ,RHOAX,RHOAY,RHOAZ
    
       
   
contains


    SUBROUTINE PSGLD(atm_cnt,ix,ih,numex,natex,AMASS,crd,vel, rem)
!-----------------------------------------------------------------------
!     This routine performs initiation for the Self-Guided        
!       Langevin Dynamcs (SGLD) simulaiton                 
!
      use md_scheme, only: gamma_ln
      use findmask,only:atommask
      implicit none
#include "../include/md.h"
#include "../include/memory.h"
#ifdef MPI
#  include "parallel.h"
#endif
#  include "extra_pts.h"
      INTEGER atm_cnt,ix(*),numex(*),natex(*)
      character(len=4) :: ih(*)
      _REAL_ AMASS(*),crd(3,*),vel(3,*)
      integer rem
      INTEGER I,I3,M,ierror,j,nsgsubi,idx_nbex,jatm
      _REAL_ AMASSI,XI3,VI3,ekin,ekinsg,GAMM,FACT1,FACT2
      logical is_langevin  ! Is this a Langevin dynamics simulation
!
#ifdef MPI
#ifdef VSCODE
!     mpi debug
      i=1
      do while (i==0)
        call sleep(5)
      enddo
#endif /* VSCODE */
#endif /* MPI */
!
      is_langevin = gamma_ln > 0.0d0
      TSGLD = (isgld >0)
      tsggamma = (isgld == 1)
      tsgbond = (sgtype > 1).and.(sgtype < 4)
      tsgmap = (sgtype > 3)
      tsgldgle = (abs(sgfg) > 1.0d-6).and.is_langevin
      trxsgld = rem > 0 .and. tsgld
!  Check for invalid sgld setting
      IF(ISGSTA < 1)ISGSTA=1
      IF(ISGEND > NATOM .OR. ISGEND < 1)ISGEND=NATOM
      IF(TSGAVG.LT.DT)TSGAVG=DT
      IF(TSGAVP.LT.DT)TSGAVP=10.0D0*TSGAVG
      SGAVG1=DT/TSGAVG
      SGAVG0=1.0D0-SGAVG1
      SGAVP1=DT/TSGAVP
      SGAVP0=1.0D0-SGAVP1
      gammas=gamma_ln/20.455d0
      tsgfac=1.0d0/20.455d0/tsgavg
      sglabel="SGMD: "
      if(is_langevin)sglabel="SGLD: "
      tsgset=temp0
      if(tsgset<1.0d-6)tsgset=300.0d0
      if(isgld==3)then
        if(sgft<=1.0d0)then
          if(ABS(sgft)>ABS(sgff))then
            sgfti=sgft
            FACT1=9.0d0-SQRT(81.0d0-12.0d0*SGFTI*SGFTI*SGFTI)
            FACT2=(ABS(FACT1)*1.5d0)**(1.0/3.0d0)
            PSGLDG=SIGN(1.0d0,FACT1)*(FACT2/3.0d0+SGFTI/FACT2)-1.0d0
            sgffi=psgldg
          else
            sgffi=sgff
            sgfti=(1.0d0+SGFFI)*(1.0d0+SGFFI)-1.0d0/(1.0d0+SGFFI)
          endif
        else
          sgffi=sgff
          sgfti=(1.0d0+SGFFI)*(1.0d0+SGFFI)-1.0d0/(1.0d0+SGFFI)
        endif
      else 
        sgfti=sgft 
        sgffi=sgff
      endif
      if(ABS(SGFTI)>1.0d-8.and. sgfti<=1.0d0)then
        FACT1=9.0d0-SQRT(81.0d0-12.0d0*SGFTI*SGFTI*SGFTI)
        FACT2=(ABS(FACT1)*1.5d0)**(1.0/3.0d0)
        PSGLDG=SIGN(1.0d0,FACT1)*(FACT2/3.0d0+SGFTI/FACT2)-1.0d0
      else if(sgfti>-1.0d0)then
        PSGLDG=SQRT(1.0d0+SGFTI)-1.0d0
      else
        PSGLDG=-1.0d0
      endif
      if(tsgldgle)then
        sgfgi=sgfg
        if(sgfgi<=1.0d0)then
          FSGLDG=SQRT(1.0d0-SGFGI)-1.0d0
        else 
          fsgldg=-1.0d0
        endif
      else
        sgfgi=0.0d0
        FSGLDG=0.0d0
      endif
      IF(TEMPSG>1.0d-6)THEN
        ! when tempsg is set
        FACT1=1.0d0-TSGSET/TEMPSG
        IF(SGFTI*SGFTI>1.0d-8)THEN
          SGFFI=PSGLDG-FACT1
        ELSE 
          PSGLDG=SGFFI+FACT1
          SGFTI=(1.0d0+PSGLDG)*(1.0d0+PSGLDG)-1.0d0/(1.0d0+PSGLDG)
        ENDIF
      ELSE
        TEMPSG=TSGSET/(1.0d0-PSGLDG+SGFFI)
      ENDIF
#ifdef MPI
      if(numtasks>1)then
        call mpi_bcast(sgmask,256,MPI_CHARACTER,0,commsander,ierror)
      endif
#endif
      if(allocated(sgatoms))deallocate(sgatoms)
      allocate(sgatoms(NATOM),stat=ierror)
      REQUIRE( ierror == 0 )
      call atommask( natom, nres, 0, ih(m04), ih(m06), &
                ix(i02), ih(m02), crd, sgmask, sgatoms )
      !     allocate working arrays
      allocate( avgx0(3,natom),avgx1(3,natom),avgx2(3,natom),  &
      avgr(3,natom ),sgfps(natom ),sgpps(natom ),&
                stat=ierror)
      REQUIRE( ierror == 0 )
     ! build bidireectional  exclusion lists
      if(tsgbond)then 
        allocate(atm_sg_maskdata(natom), &
                 atm_sg_mask(nnb*2), &
                 sgmass(natom),stat = ierror)
        REQUIRE( ierror == 0 )

        call make_sgavg_mask_list(natom,nnb,ix(i08), ix(i10),numnb14,ix(inb_14))
        REQUIRE( ierror == 0 )
      endif
      if(tsgmap)then
        gxmin=1.0d8
        gxmax=-1.0d8
        gymin=1.0d8
        gymax=-1.0d8
        gzmin=1.0d8
        gzmax=-1.0d8
      endif
 !  !    Initialize arrays
      gamm=sqrt(dt/tsgavg)
      ekin=0.0d0
      ekinsg=0.0d0
      sgmsum=0.0d0
      nsgatom=0
      DO I=1,NATOM
        AMASSI = AMASS(I)
        !IF((I>=ISGSTA).AND.(I<=ISGEND))THEN
        !ENDIF
        DO M=1,3
          xi3=crd(m,i)
          vi3=vel(m,i)
          ekin=ekin+amassi*vi3*vi3
          if(tsgbond)then 
            nsgsubi=atm_sg_maskdata(i)%cnt
            xi3=amassi*xi3
            vi3=amassi*vi3
            sgmass(i)=amassi
            do j=1,nsgsubi
              idx_nbex=atm_sg_maskdata(i)%offset + j 
              jatm=atm_sg_mask(idx_nbex)
              xi3=xi3+amass(jatm)*crd(m,jatm)
              vi3=vi3+amass(jatm)*vel(m,jatm)
              sgmass(i)=sgmass(i)+amass(jatm)
            enddo
            xi3=xi3/sgmass(i)
            vi3=vi3/sgmass(i)
          endif
          avgx0(m,i)=xi3
          avgx1(m,i)=xi3-vi3*sqrt(sgavg1)/tsgfac
          avgx2(m,i)=avgx1(m,i)-vi3*sgavg1/tsgfac
          avgr(m,i)=0.0d0
          ekinsg=ekinsg+amassi*vi3*vi3
        END DO
        sgpps(i)=amassi*0.001987*tsgset*sgavg1
        sgfps(i)=-0.01*sgpps(i)/20.455d0
        if(i>=isgsta.and.i<=isgend.and.sgatoms(i)>0)then
          nsgatom=nsgatom+1
          sgmsum=sgmsum+amassi
        else
          sgatoms(i)=0
        endif
      END DO
      epotlf=2.0d10
      epotllf=2.0d10
      templf=tsgset*sgavg1
      temphf=tsgset-templf
      sgwt=0.0d0
      sggamma=0.01d0/20.455d0
      com0sg=0.0d0
      com1sg=0.0d0
      com2sg=0.0d0
#ifdef MPI
      if(mytaskid.eq.0)THEN
#endif
      write(6,910)isgsta,isgend,nsgatom
      write(6,915)tsgavg,tsgavp
      if(is_langevin)then
        if(tsgldgle)then
          write(6,928)
          write(6,927)sgfgi,fsgldg
        else
          write(6,940)
        endif
        write(6,930)gamma_ln
      else
          write(6,941)
      endif
      write(6,925)sgfti,psgldg
      write(6,926)sgffi
      write(6,932)tempsg
      if(tsgbond)then
        if(sgtype==2)then
          write(6,942)
        else
          write(6,943)
        endif
      endif      
      if(tsgmap)then
        write(6,948)sgsize
      endif      
      write(6,935)
#ifdef MPI
    ENDIF
#endif
910   format("  _________________ SGMD/SGLD parameters _________________"/  &
      "  Parameters for self-guided Molecular/Langevin dynamics (SGMD/SGLD) simulation"//  &
          "  Guiding range from ",i5,"  to ",i8, " with ",i8," guiding atoms")
915   format("  Local averaging time: tsgavg: ",f10.4," ps,  tsgavp: ",f10.4," ps")
925   format("  sgfti: ",f8.4," psgldg: ",f8.4)
926   format("  sgffi: ",f8.4)
927   format("  momentum factor sgfgi= ",f8.4," random force factor fsgldg=",f8.4)
928   format("  SGLD-GLE method is used to mantain a canonical distribution. ")
932   format("  Guided sampling effective temperature (TEMPSG): ",f8.2)
940   format("  SGLDg  method is used to enhance conformational search. ")
941   format("  SGMDg  method is used to enhance conformational search. ")
942   format("  SGTYPE=2, Guiding forces are averaged over 1-2,1-3 bonded structures" )
943   format("  SGTYPE=3, Guiding forces are averaged over 1-2,1-3,1-4 bonded structures" )
948   format("  SGTYPE=4, Guiding forces are averaged over a cutoff: ",F8.4 )
930   format("  Collision frequency:",f8.2," /ps" )
935   format("  Output properties:"    /  &
             "  SGMD/SGLD:  SGGAMMA TEMPLF  TEMPHF  EPOTLF EPOTHF EPOTLLF SGWT" /  &
             "         SGMD/SGLD weighting factor =exp(SGWT)"/  &
              " _______________________________________________________"/)
      RETURN
      END SUBROUTINE PSGLD


      subroutine make_sgavg_mask_list(atm_cnt, nnb, numex, natex,numnb14,nb_14_list)

     
        implicit none
      
      ! Formal arguments:
      
        integer               :: atm_cnt,nnb,numnb14
        ! Excluded atom count for each atom.
        integer               :: numex(atm_cnt)
        ! Excluded atom concatenated list:
        integer               :: natex(nnb),nb_14_list(3,numnb14)
      
      ! Local variables:
      
        integer               :: atm_i, atm_j
        integer               :: lst_idx, sublst_idx, num_sublst
        integer               :: mask_idx
        integer               :: offset
        integer               :: total_excl
        ! nb14 variables
        integer               :: j,idx14,cnt14,list14(20)
        logical               :: addj
      
      ! Double the mask to deal with our list generator
      
      ! Pass 1: get pointers, check size
      
        lst_idx = 0
      
        atm_sg_maskdata(:)%cnt = 0       ! array assignment
      
        do atm_i = 1, atm_cnt - 1         ! last atom never has any...
          if(sgtype==2)then 
            cnt14=0
            do idx14=1,numnb14
              if(nb_14_list(1,idx14)==atm_i)then 
                cnt14=cnt14+1
                list14(cnt14)=nb_14_list(2,idx14)
              endif
            end do
          endif
          num_sublst = numex(atm_i)
          do sublst_idx = 1, num_sublst
            atm_j = natex(lst_idx + sublst_idx)
            if (atm_j .gt. 0 ) then
              addj=.true.
              if(sgtype==2)then 
                do j=1,cnt14
                  if(list14(j)==atm_j)addj=.false.
                enddo
              endif 
              if(addj)then
                atm_sg_maskdata(atm_i)%cnt = atm_sg_maskdata(atm_i)%cnt + 1
                atm_sg_maskdata(atm_j)%cnt = atm_sg_maskdata(atm_j)%cnt + 1
              endif 
            end if
          end do
          lst_idx = lst_idx + num_sublst
        end do
      
        total_excl = 0
      
        do atm_i = 1, atm_cnt
          total_excl = total_excl + atm_sg_maskdata(atm_i)%cnt
        end do
        !write(6,*)"total_excl, nnb: ",total_excl, nnb
        if (total_excl .gt. nnb*2) then
          write(6, '(a,a)') "SGBOND: ", &
               'The total number of sg substructure exceeds that stipulated by the'
          write(6, '(a,a)') "SGBOND: ", &
               'prmtop.  This is likely due to a very high density of added extra points.'
          write(6, '(a,a)') "SGBOND: ", &
               'Scale back the model detail, or contact the developers for a workaround.'
          call mexit(6, 1)
        end if
      
        offset = 0
      
        do atm_i = 1, atm_cnt
          atm_sg_maskdata(atm_i)%offset = offset
          offset = offset + atm_sg_maskdata(atm_i)%cnt
        end do
      
      ! Pass 2: fill mask array
      
        lst_idx = 0
      
        atm_sg_maskdata(:)%cnt = 0       ! array assignment
        
          do atm_i = 1, atm_cnt - 1
            if(sgtype==2)then 
              cnt14=0
              do idx14=1,numnb14
                if(nb_14_list(1,idx14)==atm_i)then 
                  cnt14=cnt14+1
                  list14(cnt14)=nb_14_list(2,idx14)
                endif
              end do
            endif
            num_sublst = numex(atm_i)
            do sublst_idx = 1, num_sublst
              atm_j = natex(lst_idx + sublst_idx)
              if (atm_j .gt. 0 ) then
                addj=.true.
                if(sgtype==2)then 
                  do j=1,cnt14
                    if(list14(j)==atm_j)addj=.false.
                  enddo
                endif 
                if(addj)then
                  atm_sg_maskdata(atm_j)%cnt = atm_sg_maskdata(atm_j)%cnt + 1
                  mask_idx = atm_sg_maskdata(atm_j)%offset + &
                           atm_sg_maskdata(atm_j)%cnt
                  atm_sg_mask(mask_idx) = atm_i
      
                  atm_sg_maskdata(atm_i)%cnt = atm_sg_maskdata(atm_i)%cnt + 1
                  mask_idx = atm_sg_maskdata(atm_i)%offset + &
                           atm_sg_maskdata(atm_i)%cnt
                  atm_sg_mask(mask_idx) = atm_j
                endif
              end if
            end do
            lst_idx = lst_idx + num_sublst
          end do
        return
      
      end subroutine make_sgavg_mask_list
      
      
    subroutine sg_fix_degree_count(sgsta_rndfp, sgend_rndfp, ndfmin, rndf)
!-----------------------------------------------------------------------
!     Correct the total number of degrees of freedom for a translatable COM,
!       and compute the number of degrees of freedom in the SGLD part.
!       The latter is mostly done by the caller to avoid passing the long
!       argument list needed by routine degcnt which also differs between
!       sander and pmemd.
!
      implicit none
      _REAL_, intent(in)    :: sgsta_rndfp, sgend_rndfp
      integer, intent(in)   :: ndfmin
      _REAL_, intent(inout) :: rndf

      sgrndf = sgend_rndfp - sgsta_rndfp
      return
      end subroutine sg_fix_degree_count



      SUBROUTINE SGLDW(NATOM,ISTART,IEND, &
             DTX,TEMP0,ENER,AMASS,WINV,crd,Frc,Vel)
!-----------------------------------------------------------------------
!     This routine perform SGLD integration        
!
      use state
      use random, only: GAUSS
      implicit none
#ifdef MPI
      integer ierr
# include "parallel.h"
      _REAL_ temp1(20)
# ifndef USE_MPI_IN_PLACE
      _REAL_ :: temp2(20)
# endif
#endif
      INTEGER NATOM,ISTART,IEND
      _REAL_ DTX,TEMP0
      type(state_rec) :: ener
      _REAL_ AMASS(*),WINV(*),crd(3,*),Frc(3,*),Vel(3,*)
!
      INTEGER I,M,JSTA,JEND,j,nsgsubi,idx_nbex,jatm
      _REAL_ BOLTZ,AMASSI,TEMPI
      _REAL_ FACT,WFAC,GAM,RSD,FLN
      _REAL_ EKIN,EKINSG,SGBETA
      _REAL_ dcom1(3),dcom2(3),com0(3),com1(3),com2(3)
      _REAL_ sumgam,sumfp,sumpp,sumgv,sumpv
      _REAL_ sggammai,avgpi3,pi3t,avgdfi3,avgri3,fsgpi,fsgfi,fsgi3,frici
      _REAL_ xi3,x1i3,x2i3,vi3t,vi3,fi3,x0i3,dx0i3,psgi3
      _REAL_ avgpi(3),avgfi(3)
      PARAMETER (BOLTZ = 1.987192d-3)
!
    !
        if(isgld.eq.2)then
          dcom1=0.0d0
          dcom2=0.0d0
        else
          dcom1=com0sg-com1sg
          dcom2=com0sg-com2sg
        endif
    ! build spatial average maps
        if(tsgmap)then
          call mapbuild(ISTART,IEND,dtx,amass,crd,vel,dcom1,dcom2)
        endif
            !
        gam=gammas*dtx
        JSTA=ISTART
        JEND=IEND
        IF(JSTA < ISGSTA)JSTA=ISGSTA
        IF(JEND > ISGEND)JEND=ISGEND
        sumgam=0.0d0
        EKIN=0.0D0
        EKINSG=0.0D0
        com0=0.0d0
        com1=0.0d0
        com2=0.0d0
        DO  I = 1,NATOM 
          AMASSI = AMASS(I)
          WFAC =  2.0D0*DTX*WINV(I)
          RSD = SQRT(2.D0*GAMMAS*BOLTZ*TEMP0*AMASSI/DTX)
          IF(I>=JSTA.AND.I<=JEND.AND.SGATOMS(I)>0)THEN
            ! sggamma
          sggammai=-sgfps(i)/sgpps(i)
          sumgam=sumgam+sggammai
          if(tsggamma)sggammai=sggamma
          if(tsgmap)then
            ! spatial average using map
            call mapvalues(amassi,crd(:,i),avgpi,avgfi)
          endif
          sumfp=0.0d0
          sumpp=0.0d0
          sumgv=0.0d0
          sumpv=0.0d0
          DO  M = 1,3
!   Keep random number series the same as that in a single cpu simulation
            CALL GAUSS( 0.D0, RSD, FLN )
            if(tsgmap)then
              ! spatial average using map
                avgpi3=avgpi(m)
                avgdfi3=avgfi(m)
            else
              ! avg(x)
              xi3=crd(m,i)
              if(tsgbond)then 
                nsgsubi=atm_sg_maskdata(i)%cnt
                xi3=amassi*xi3
                do j=1,nsgsubi
                  idx_nbex=atm_sg_maskdata(i)%offset + j 
                  jatm=atm_sg_mask(idx_nbex)
                  xi3=xi3+amass(jatm)*crd(m,jatm)
                enddo
                xi3=xi3/sgmass(i)
              endif
                ! adjust previous averages
                x0i3=xi3-vel(m,i)*dtx
                dx0i3=x0i3-avgx0(m,i)
                x1i3=avgx1(m,i)+dx0i3
                x2i3=avgx2(m,i)+dx0i3
                psgi3=tsgfac*amassi*(x0i3-x1i3)
                avgx0(m,i)=xi3
                ! Obtain local averages
                x1i3=sgavg0*(x1i3+dcom1(m))+sgavg1*xi3
                avgx1(m,i)=x1i3
                ! avgavg(x)
                x2i3=sgavg0*(x2i3+dcom2(m))+sgavg1*x1i3
                avgx2(m,i)=x2i3
              ! avg(p)
              avgpi3=tsgfac*amassi*(xi3-x1i3)
              pi3t=(avgpi3-sgavg0*psgi3)/sgavg1
              ! avg(f-avg(f))
              avgdfi3=tsgfac*(pi3t-2.0d0*avgpi3+tsgfac*amassi*(x1i3-x2i3))
              com0(m)=com0(m)+amassi*xi3
              com1(m)=com1(m)+amassi*x1i3
              com2(m)=com2(m)+amassi*x2i3
            endif /* tsgmap */
            ! sum(avg(f-avg(f))avg(p))
              sumfp=sumfp+avgdfi3*avgpi3
              ! sum(avg(p)avg(p))
              sumpp=sumpp+avgpi3*avgpi3
              ! average random forces
              avgri3=sgavg0*avgr(m,i)+sgavg1*fln
              avgr(m,i)=avgri3
              ! guiding forces
              fsgpi=(sgfti*sggammai)*avgpi3
              fsgfi=sgffi*avgdfi3
              fsgi3=sgfgi*gammas*avgpi3+(fsgldg-sgffi)*avgri3+fsgpi+fsgfi
              fi3=frc(m,i)+fln+fsgi3
              frc(m,i)=fi3
              ! estimate velocity at t+dt/2
              ! Using volocities at t avoid SHAKE complication
              vi3t=vel(m,i)
              ! sum(g*v)
              sumgv=sumgv+fsgi3*vi3t
              ! sum(p*v)
              sumpv=sumpv+amassi*vi3t*vi3t
              ekin=ekin+amassi*vi3t*vi3t
              ekinsg=ekinsg+avgpi3*avgpi3/amassi
          end do
            ! <(avg(f-avg(f))avg(v))>
            sgfps(i)=sgavp0*sgfps(i)+sgavp1*sumfp
            ! <(avg(p)avg(v))>
            sgpps(i)=sgavp0*sgpps(i)+sgavp1*sumpp
            ! energy conservation friction constant
            if(sumpv<1.0d-8)then
              sgbeta=0.0d0
            else
              sgbeta=sumgv/(2.0d0*sumpv/(2.0d0+gam)-sumgv*dtx)
            endif
            !sgbeta=0.0d0
            fact=dtx*(gammas+sgbeta)
            do  m = 1,3
              fi3=frc(m,i)
              vi3t=((2.0d0-fact)*vel(m,i)+fi3*wfac)/(2.0d0+fact)
              vel(m,i)=vi3t
            end do
          ELSE 
               ! without guiding forces
            do  m = 1,3
              !   generate random number 
              call gauss( 0.d0, rsd, fln )
              FI3=FRC(m,I)+fln
              Frc(m,I)=FI3
              vi3t=((2.0d0-gam)*vel(m,i)+fi3*wfac)/(2.0d0+gam)
              vel(m,i)=vi3t
          END DO
        ENDIF
      END DO
#ifdef MPI
        IF(SANDERSIZE > 1)THEN
!  Combining all node results
!
          TEMP1(1)=sumgam
          TEMP1(2)=EKIN
          TEMP1(3)=EKINSG
          temp1(4:6)=com0
          temp1(7:9)=com1
          temp1(10:12)=com2

# ifdef USE_MPI_IN_PLACE
          call mpi_allreduce(MPI_IN_PLACE,temp1,12,&
             MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ierr)
          sumgam=TEMP1(1)
          EKIN=TEMP1(2)
          EKINSG=TEMP1(3)
          com0=temp1(4:6)
          com1=temp1(7:9)
          com2=temp1(10:12)

#else
          CALL MPI_ALLREDUCE(TEMP1,TEMP2,12, &
          MPI_DOUBLE_PRECISION,MPI_SUM,COMMSANDER,IERR)
          sumgam=TEMP2(1)
          EKIN=TEMP2(2)
          EKINSG=TEMP2(3)
          com0=temp2(4:6)
          com1=temp2(7:9)
          com2=temp2(10:12)

# endif
        ENDIF
#endif
    ! Estimate low frequency temperatures
        !TEMPI=EKIN/sgrndf/BOLTZ
        TEMPI=tsgset
        !TEMPLF=SGAVP0*TEMPLF+SGAVP1*EKINSG/sgrndf/BOLTZ
        TEMPLF=SGAVP0*TEMPLF+SGAVP1*EKINSG*tsgset/EKIN
        TEMPHF=Tempi-TEMPLF
        sggamma=sumgam/nsgatom
        sgscale=20.455d0*sggamma
        com0sg=com0/sgmsum
        com1sg=com1/sgmsum
        com2sg=com2/sgmsum
        CALL SGENERGY(ENER)
        RETURN
        END SUBROUTINE SGLDW

        SUBROUTINE SGMDW(NATOM,ISTART,IEND, &
             DTX,ENER,AMASS,WINV,crd,Frc,Vel)
!-----------------------------------------------------------------------
!     This routine calculate guiding force using SGLD method 
!     for MD simulation        
!
      use state
      use random, only: GAUSS
      implicit none
#ifdef MPI
      integer ierr
# include "parallel.h"
      _REAL_ temp1(20)
# ifndef USE_MPI_IN_PLACE
      _REAL_ :: temp2(20)
# endif
#endif
      INTEGER NATOM,ISTART,IEND,NTP
      _REAL_ DTX
      type(state_rec) :: ener
      _REAL_ AMASS(*),WINV(*),crd(3,*),Frc(3,*),Vel(3,*)
!
      INTEGER JSTA,JEND,I,M,j,nsgsubi,idx_nbex,jatm
      _REAL_ BOLTZ,AMASSI,TEMPI
      _REAL_ FACT,WFAC
      _REAL_ EKIN,EKINSG,SGBETA
      _REAL_ dcom1(3),dcom2(3),com0(3),com1(3),com2(3)
      _REAL_ sumgam,sumfp,sumpp,sumgv,sumpv
      _REAL_ sggammai,avgpi3,pi3t,avgdfi3,avgri3,fsgpi,fsgfi,fsgi3,frici
      _REAL_ avgpi(3),avgfi(3)
      _REAL_ xi3,x1i3,x2i3,vi3t,vi3,fi3,x0i3,dx0i3,psgi3
      PARAMETER (BOLTZ = 1.987192d-3)
!
    !
        if(isgld.eq.2)then
          dcom1=0.0d0
          dcom2=0.0d0
        else
          dcom1=com0sg-com1sg
          dcom2=com0sg-com2sg
        endif
    ! build spatial average maps
        if(tsgmap)then
            call mapbuild(ISTART,IEND,dtx,amass,crd,vel,dcom1,dcom2)
        endif
        !
        JSTA=ISTART
        JEND=IEND
        IF(JSTA < ISGSTA)JSTA=ISGSTA
        IF(JEND > ISGEND)JEND=ISGEND
        sumgam=0.0d0
        EKIN=0.0D0
        EKINSG=0.0D0
        com0=0.0d0
        com1=0.0d0
        com2=0.0d0
        DO  I = jsta,jend
          if(sgatoms(i).eq.0)cycle
          AMASSI = AMASS(I)
          WFAC =  DTX*0.5D0*WINV(I)
          sggammai=-sgfps(i)/sgpps(i)
          sumgam=sumgam+sggammai
          if(tsggamma)sggammai=sggamma
          if(tsgmap)then
            ! spatial average using map
            call mapvalues(amassi,crd(:,i),avgpi,avgfi)
          endif

          sumfp=0.0d0
          sumpp=0.0d0
          sumgv=0.0d0
          sumpv=0.0d0
          DO  M = 1,3
            if(tsgmap)then
              ! spatial average using map
                avgpi3=avgpi(m)
                avgdfi3=avgfi(m)
            else
            ! avg(x)
            xi3=crd(m,i)
            if(tsgbond)then 
              nsgsubi=atm_sg_maskdata(i)%cnt
              xi3=amassi*xi3
              do j=1,nsgsubi
                idx_nbex=atm_sg_maskdata(i)%offset + j 
                jatm=atm_sg_mask(idx_nbex)
                xi3=xi3+amass(jatm)*crd(m,jatm)
              enddo
              xi3=xi3/sgmass(i)
            endif
                ! adjust previous averages
                x0i3=xi3-vel(m,i)*dtx
                dx0i3=x0i3-avgx0(m,i)
                x1i3=avgx1(m,i)+dx0i3
                x2i3=avgx2(m,i)+dx0i3
                psgi3=tsgfac*amassi*(x0i3-x1i3)
                avgx0(m,i)=xi3
                ! Obtain local averages
                x1i3=sgavg0*(x1i3+dcom1(m))+sgavg1*xi3
                avgx1(m,i)=x1i3
                ! avgavg(x)
                x2i3=sgavg0*(x2i3+dcom2(m))+sgavg1*x1i3
                avgx2(m,i)=x2i3
              ! avg(p)
              avgpi3=tsgfac*amassi*(xi3-x1i3)
              pi3t=(avgpi3-sgavg0*psgi3)/sgavg1
              ! avg(f-avg(f))
              avgdfi3=tsgfac*(pi3t-2.0d0*avgpi3+tsgfac*amassi*(x1i3-x2i3))
                com0(m)=com0(m)+amassi*xi3
                com1(m)=com1(m)+amassi*x1i3
                com2(m)=com2(m)+amassi*x2i3
            endif /* tsgmap */
            ! sum(avg(f-avg(f))avg(p))
            sumfp=sumfp+avgdfi3*avgpi3
            ! sum(avg(p)avg(p))
            sumpp=sumpp+avgpi3*avgpi3
            ! guiding forces
            fsgpi=sgfti*sggammai*avgpi3
            fsgfi=sgffi*avgdfi3
            fsgi3=fsgpi+fsgfi
            fi3=frc(m,i)+fsgi3
            frc(m,i)=fi3
            ! estimate velocity at t+dt/2
            !vi3t=vel(m,i)+fi3*wfac
            ! Using volocities at t avoid SHAKE complication
            vi3t=vel(m,i)
            ! sum(g*v)
            sumgv=sumgv+fsgi3*vi3t
            ! sum(p*v)
            sumpv=sumpv+amassi*vi3t*vi3t
            ekin=ekin+amassi*vi3t*vi3t
            ekinsg=ekinsg+avgpi3*avgpi3/amassi
          enddo
            ! <(avg(f-avg(f))avg(v))>
          sgfps(i)=sgavp0*sgfps(i)+sgavp1*sumfp
          ! <(avg(p)avg(v))>
          sgpps(i)=sgavp0*sgpps(i)+sgavp1*sumpp
          ! energy conservation friction constant
          if(sumpv<1.0d-8)then
            sgbeta=0.0d0
          else
            sgbeta=2.0d0*sumgv/(2.0d0*sumpv-sumgv*dtx)
          endif
          fact=sgbeta/(1.0d0+0.5d0*sgbeta*dtx)
          do  m = 1,3
            fi3=frc(m,i)
            vi3t = vel(m,i) + fi3*wfac
            frici=fact*amassi*vi3t
            frc(m,i)=fi3-frici
          end do
        END DO
#ifdef MPI
        IF(SANDERSIZE > 1)THEN
!  Combining all node results
!
          TEMP1(1)=sumgam
          TEMP1(2)=EKIN
          TEMP1(3)=EKINSG
          temp1(4:6)=com0
          temp1(7:9)=com1
          temp1(10:12)=com2
# ifdef USE_MPI_IN_PLACE
          call mpi_allreduce(MPI_IN_PLACE,temp1,12, &
            MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ierr)
          sumgam=TEMP1(1)
          EKIN=TEMP1(2)
          EKINSG=TEMP1(3)
          com0=temp1(4:6)
          com1=temp1(7:9)
          com2=temp1(10:12)
#else
          CALL MPI_ALLREDUCE(TEMP1,TEMP2,12, &
          MPI_DOUBLE_PRECISION,MPI_SUM,COMMSANDER,IERR)
          sumgam=TEMP2(1)
          EKIN=TEMP2(2)
          EKINSG=TEMP2(3)
          com0=temp2(4:6)
          com1=temp2(7:9)
          com2=temp2(10:12)
# endif
        ENDIF
#endif
    ! Estimate low frequency temperatures
        !TEMPI=EKIN/sgrndf/BOLTZ
        TEMPI=tsgset
        !TEMPLF=SGAVP0*TEMPLF+SGAVP1*EKINSG/sgrndf/BOLTZ
        TEMPLF=SGAVP0*TEMPLF+SGAVP1*EKINSG*tsgset/EKIN
        TEMPHF=TEMPI-TEMPLF
        sggamma=sumgam/nsgatom
        sgscale=20.455d0*sggamma
        com0sg=com0/sgmsum
        com1sg=com1/sgmsum
        com2sg=com2/sgmsum
        ! update accumulators
        CALL SGENERGY(ENER)
        RETURN
        END SUBROUTINE SGMDW
     
        SUBROUTINE SGENERGY(ENER)
!-----------------------------------------------------------------------
!     This routine set the ener fields of SGLD variables
!
      use state
      implicit none
      type(state_rec) :: ener
      _REAL_ BOLTZ
      PARAMETER (BOLTZ = 1.987192d-3)
      _REAL_ EPOTI
    ! Weighting accumulators
        EPOTI=ENER%POT%TOT
        IF(EPOTLF>1.0D10)THEN
          EPOTLF=EPOTI
          EPOTLLF=EPOTI
        ELSE
          EPOTLF=SGAVG0*EPOTLF+SGAVG1*EPOTI
          EPOTLLF=SGAVG0*EPOTLLF+SGAVG1*EPOTLF
        ENDIF
        EPOTHF=EPOTI-EPOTLF
        sgwt=(psgldg-sgffi)*(epotlf-epotllf)/(boltz*tsgset)
    ! Update ENER structure
        ENER%SGLD%SGSCALE=SGSCALE
        ENER%SGLD%TEMPLF=TEMPLF
        ENER%SGLD%TEMPHF=TEMPHF
        ENER%SGLD%EPOTLF=EPOTLF
        ENER%SGLD%EPOTHF=EPOTHF
        ENER%SGLD%EPOTLLF=EPOTLLF
        ENER%SGLD%SGWT=SGWT
       RETURN
       END SUBROUTINE SGENERGY

#ifdef MPI

!*********************************************************************
!               SUBROUTINE SGLD_EXCHG
!*********************************************************************
!  exchange all data in the SGLDR common block
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine sgld_exchg(irep)

   implicit none
#  include "parallel.h"

   integer, intent(in) :: irep
   integer  ierror, istatus(mpi_status_size)
   call mpi_sendrecv_replace(sgft,nsgld_real, mpi_double_precision, &
                   irep, 511, irep, 511, commmaster, istatus, ierror)
    !call mpi_barrier(commmaster, ierror)
   return
   end subroutine sgld_exchg

!*********************************************************************
!               SUBROUTINE REMD_SCALE_VELO
!*********************************************************************
! Scale velocities based on new temps after exchange
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine rxsgld_scale(stagid,nr,myscaling,amass,crd, vel)

   implicit none
#  include "parallel.h"

   integer, intent(in) :: stagid,nr
   _REAL_, intent(in) :: myscaling
   _REAL_,  intent(in) :: crd(3,nr)
   _REAL_,  intent(inout) :: vel(3,nr)
   _REAL_, dimension (*), intent(in)    :: amass
   integer ierror
   integer i,j,jsta,jend
   _REAL_ amassi,xi,x0i,x1i,x2i
   _REAL_ temp1(10),temp2(10)
!--------------------
         !if (sanderrank==0) then
         !   write (6,'(a,i4,2x,f8.3,a,2f8.3)') &
         !      "RXSGLD: stagid, scalsg  ",stagid,myscalsg,&
         !      " to match sgft,tempsg,: ",sgft,tempsg
         !endif
      call mpi_bcast(stagid,1,mpi_integer,0,commsander,ierror)

         ! All processes scale velocities.
         ! DAN ROE: This could potentially be divided up as in runmd
         !  since when there are mutiple threads per group each thread 
         !  only ever knows about its own subset of velocities anyway.
#ifdef VERBOSE_REMD
         if (sanderrank==0) then
            write (6,'(a,f8.3,a,f8.3)') &
               "| RXSGLD: scaling guiding properties by ",myscalsg,&
               " to match a new guiding factors sgfti, sgffi ",sgfti,sgffi
         endif
#endif
! ---=== Broadcast RXSGLD guiding effect ===---
      IF(SANDERSIZE > 1)call mpi_bcast(sgft,nsgld_real,mpi_double_precision,&
                                              0,commsander,ierror)
      if (myscaling > 0.0d0) then
        vel(:,:)=myscaling*vel(:,:)
      endif                         
      if (myscalsg > 0.0d0) then
        templf=myscalsg*myscalsg*templf
        avgr(:,:)=myscalsg*avgr(:,:)
        sgfps(:)=myscalsg*sgfps(:)
        sgpps(:)=myscalsg*myscalsg*sgpps(:)
        jsta = iparpt(mytaskid) + 1
        jend = iparpt(mytaskid+1)
        IF(JSTA < ISGSTA)JSTA=ISGSTA
        IF(JEND > ISGEND)JEND=ISGEND
        do i = jsta,jend
            amassi=amass(i)
            do j=1,3
              xi=crd(j,i)
              x0i=avgx0(j,i)
              x1i=avgx1(j,i)
              x2i=avgx2(j,i)
              avgx1(j,i)=xi+myscalsg*(x1i-x0i)
              avgx2(j,i)=xi+myscalsg*(x2i-x0i)
           enddo
         enddo
      endif
   return

end subroutine rxsgld_scale

!*********************************************************************
!               FUNCTION TEMPSGLOOKUP
!*********************************************************************
! lookup temp in templist and return its index
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
integer function tempsglookup(numreps,tempsg,tempsgs)

   implicit none
#  include "parallel.h"

   integer numreps
   _REAL_, intent(in) :: tempsg
   _REAL_, dimension(numreps), intent(in) :: tempsgs

   integer i
   
   tempsglookup=0
   do i = 1, numreps
      if(abs(tempsg-tempsgs(i)) < 1.0d-6 ) then
         if(tempsglookup>0)then
            write (6,*) "================================"
            write (6,*) "Two replicas are the same: ",tempsglookup,i
            write (6,*) "================================"
            call mexit(6,1)
         endif
         tempsglookup = i
      end if
   end do
   return
end function tempsglookup

!*********************************************************************
!               FUNCTION STAGIDLOOKUP
!*********************************************************************
! lookup stagid in stagidlist and return its neighboring index
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
integer function stagidlookup(numreps,id,idtable)

   implicit none
#  include "parallel.h"

   integer numreps
   INTEGER, dimension(numreps), intent(in) :: idtable
   INTEGER, intent(in) :: id

   integer i
   
   stagidlookup=-1
   do i = 1, numreps
      if(id==idtable(i)) stagidlookup=i
   end do
   return
end function stagidlookup

!*********************************************************************
!               SUBROUTINE SORTTEMPSG
!*********************************************************************
! sort temp ascendingly
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine sorttempsg(numreps,temps,tempsgs,sgfts,sgffs,psgldgs)

   implicit none

#  include "parallel.h"

   integer numreps
   _REAL_, dimension(numreps), intent(inout) :: temps,tempsgs,sgfts,sgffs,psgldgs
   _REAL_, dimension(numreps) :: tmp
   INTEGER, dimension(numreps) :: tmpid

   _REAL_ tempt
   integer i, j, ii

   do i = 1, numreps
     tmp(i)=tempsgs(i)
     tmpid(i)=i
   enddo
   do i = 1, numreps
      do j = i + 1, numreps
         if(tmp(j) < tmp(i)) then
            tempt = tmp(i)
            tmp(i) = tmp(j)
            tmp(j) = tempt
            ii=tmpid(i)
            tmpid(i)=tmpid(j)
            tmpid(j)=ii
         end if
      end do
   end do
   do i = 1, numreps
     ii=tmpid(i)
     tmp(i)=tempsgs(ii)
   enddo
   do i = 1, numreps
     ii=tmpid(i)
     tempsgs(i)=tmp(i)
     tmp(i)=temps(ii)
   enddo
   do i = 1, numreps
     ii=tmpid(i)
     temps(i)=tmp(i)
     tmp(i)=sgfts(ii)
   enddo
   do i = 1, numreps
     ii=tmpid(i)
     sgfts(i)=tmp(i)
     tmp(i)=psgldgs(ii)
   enddo
   do i = 1, numreps
     ii=tmpid(i)
     psgldgs(i)=tmp(i)
     tmp(i)=sgffs(ii)
   enddo
   do i = 1, numreps
     sgffs(i)=tmp(i)
   enddo
   return
end subroutine sorttempsg

#endif /* MPI */


  SUBROUTINE MAPBUILD(ISTART,IEND,DELTA,AMASS,CRD,VEL,DCOM1,DCOM2)
  !-----------------------------------------------------------------------
  !     This routine build velocity and acceleration maps
  !
  !-----------------------------------------------------------------------
  
  implicit none
#ifdef MPI
     integer ierr
# include "parallel.h"
#endif
 INTEGER ISTART,IEND
  double precision CRD(3,*),VEL(3,*),AMASS(*),DCOM1(3),DCOM2(3)
  double precision DELTA,XI,YI,ZI,AMASSI
  INTEGER X0,Y0,Z0,X1,Y1,Z1
  INTEGER I000,I100,I010,I001,I110,I101,I011,I111
  double precision GX,GY,GZ,A0,B0,C0,A1,B1,C1
  double precision ABC000,ABC100,ABC010,ABC001,ABC110,ABC101,ABC011,ABC111
  INTEGER IA,I
    double precision X1I,Y1I,Z1I,X2I,Y2I,Z2I,XT,YT,ZT
    double precision GX0,GY0,GZ0,GXI,GYI,GZI,PXT,PYT,PZT,FXI,FYI,FZI
    double precision FACT1,FACT2,RHOMI
  !
  ! check mapsize
  CALL MAPINIT(ISTART,IEND,CRD)
  ! interpolate structure to protein map
  !WXW Calculate constraint factor to eliminate energy input from guiding force
        do i = ISTART,IEND
          if(i>=isgsta.and.i<=isgend)then

          AMASSI=AMASS(I)
          XI=CRD(1,I)
          YI=CRD(2,I)
          ZI=CRD(3,I)
            XT=CRD(1,I)-VEL(1,I)*DELTA-AVGX0(1,I)
            YT=CRD(2,I)-VEL(2,I)*DELTA-AVGX0(2,I)
            ZT=CRD(3,I)-VEL(3,I)*DELTA-AVGX0(3,I)
            X1I=AVGX1(1,I)+XT
            Y1I=AVGX1(2,I)+YT
            Z1I=AVGX1(3,I)+ZT
            X2I=AVGX2(1,I)+XT
            Y2I=AVGX2(2,I)+YT
            Z2I=AVGX2(3,I)+ZT
            FACT1=TSGFAC*AMASSI
            GX0=FACT1*(XI-X1I)
            GY0=FACT1*(YI-Y1I)
            GZ0=FACT1*(ZI-Z1I)
            X1I=SGAVG0*(X1I+DCOM1(1))+SGAVG1*XI
            Y1I=SGAVG0*(Y1I+DCOM1(2))+SGAVG1*YI
            Z1I=SGAVG0*(Z1I+DCOM1(3))+SGAVG1*ZI
            X2I=SGAVG0*(X2I+DCOM2(1))+SGAVG1*X1I
            Y2I=SGAVG0*(Y2I+DCOM2(2))+SGAVG1*Y1I
            Z2I=SGAVG0*(Z2I+DCOM2(3))+SGAVG1*Z1I
            AVGX0(1,I)=XI
            AVGX0(2,I)=YI
            AVGX0(3,I)=ZI
            AVGX1(1,I)=X1I
            AVGX1(2,I)=Y1I
            AVGX1(3,I)=Z1I
            AVGX2(1,I)=X2I
            AVGX2(2,I)=Y2I
            AVGX2(3,I)=Z2I
            ! avg(p)
            FACT1=TSGFAC*AMASSI
            GXI=FACT1*(XI-X1I)
            GYI=FACT1*(YI-Y1I)
            GZI=FACT1*(ZI-Z1I)
            PXT=(GXI-SGAVG0*GX0)/SGAVG1
            PYT=(GYI-SGAVG0*GY0)/SGAVG1
            PZT=(GZI-SGAVG0*GZ0)/SGAVG1
            ! average force deviation 
            FXI=tsgfac*(pxt-2.0d0*gxi+fact1*(x1i-x2i))
            FYI=tsgfac*(pyt-2.0d0*gyi+fact1*(y1i-y2i))
            FZI=tsgfac*(pzt-2.0d0*gzi+fact1*(z1i-z2i))
            ! grid distribution
            GX=(XI-GXMIN)/SGSIZE+1.0d0
            GY=(YI-GYMIN)/SGSIZE+1.0d0
            GZ=(ZI-GZMIN)/SGSIZE+1.0d0
            X0=INT(GX)
            Y0=INT(GY)
            Z0=INT(GZ)
            IF(X0<1.OR.Y0<1.OR.Z0<1)THEN
              STOP 'Position outside of lower boundary'
            ENDIF
            X1=X0+1
            Y1=Y0+1
            Z1=Z0+1
            IF(X1>NGRIDX.OR.Y1>NGRIDY.OR.Z1>NGRIDZ)THEN
              !write(mdout,*)'I= ',i,xi,yi,zi,x1,y1,z1
              STOP 'Position outside of higher boundary'
            ENDIF
            A0=X1-GX
            B0=Y1-GY
            C0=Z1-GZ
            A1=1.0d0-A0
            B1=1.0d0-B0
            C1=1.0d0-C0
            ABC000=A0*B0*C0
            ABC100=A1*B0*C0
            ABC010=A0*B1*C0
            ABC001=A0*B0*C1
            ABC110=A1*B1*C0
            ABC011=A0*B1*C1
            ABC101=A1*B0*C1
            ABC111=A1*B1*C1
            I000=X0+NGRIDX*(Y0-1+NGRIDY*(Z0-1))
            I100=I000+1
            I010=I000+NGRIDX
            I001=I000+NGRIDXY
            I110=I010+1
            I101=I001+1
            I011=I001+NGRIDX
            I111=I011+1
            RHOM(I000)=RHOM(I000)+ABC111*AMASSI
            RHOM(I100)=RHOM(I100)+ABC011*AMASSI
            RHOM(I010)=RHOM(I010)+ABC101*AMASSI
            RHOM(I001)=RHOM(I001)+ABC110*AMASSI
            RHOM(I110)=RHOM(I110)+ABC001*AMASSI
            RHOM(I101)=RHOM(I101)+ABC010*AMASSI
            RHOM(I011)=RHOM(I011)+ABC100*AMASSI
            RHOM(I111)=RHOM(I111)+ABC000*AMASSI
            RHOVX(I000)=RHOVX(I000)+ABC111*GXI
            RHOVY(I000)=RHOVY(I000)+ABC111*GYI
            RHOVZ(I000)=RHOVZ(I000)+ABC111*GZI
            RHOVX(I100)=RHOVX(I100)+ABC011*GXI
            RHOVY(I100)=RHOVY(I100)+ABC011*GYI
            RHOVZ(I100)=RHOVZ(I100)+ABC011*GZI
            RHOVX(I010)=RHOVX(I010)+ABC101*GXI
            RHOVY(I010)=RHOVY(I010)+ABC101*GYI
            RHOVZ(I010)=RHOVZ(I010)+ABC101*GZI
            RHOVX(I001)=RHOVX(I001)+ABC110*GXI
            RHOVY(I001)=RHOVY(I001)+ABC110*GYI
            RHOVZ(I001)=RHOVZ(I001)+ABC110*GZI
            RHOVX(I110)=RHOVX(I110)+ABC001*GXI
            RHOVY(I110)=RHOVY(I110)+ABC001*GYI
            RHOVZ(I110)=RHOVZ(I110)+ABC001*GZI
            RHOVX(I101)=RHOVX(I101)+ABC010*GXI
            RHOVY(I101)=RHOVY(I101)+ABC010*GYI
            RHOVZ(I101)=RHOVZ(I101)+ABC010*GZI
            RHOVX(I011)=RHOVX(I011)+ABC100*GXI
            RHOVY(I011)=RHOVY(I011)+ABC100*GYI
            RHOVZ(I011)=RHOVZ(I011)+ABC100*GZI
            RHOVX(I111)=RHOVX(I111)+ABC000*GXI
            RHOVY(I111)=RHOVY(I111)+ABC000*GYI
            RHOVZ(I111)=RHOVZ(I111)+ABC000*GZI

            RHOAX(I000)=RHOAX(I000)+ABC111*FXI
            RHOAY(I000)=RHOAY(I000)+ABC111*FYI
            RHOAZ(I000)=RHOAZ(I000)+ABC111*FZI
            RHOAX(I100)=RHOAX(I100)+ABC011*FXI
            RHOAY(I100)=RHOAY(I100)+ABC011*FYI
            RHOAZ(I100)=RHOAZ(I100)+ABC011*FZI
            RHOAX(I010)=RHOAX(I010)+ABC101*FXI
            RHOAY(I010)=RHOAY(I010)+ABC101*FYI
            RHOAZ(I010)=RHOAZ(I010)+ABC101*FZI
            RHOAX(I001)=RHOAX(I001)+ABC110*FXI
            RHOAY(I001)=RHOAY(I001)+ABC110*FYI
            RHOAZ(I001)=RHOAZ(I001)+ABC110*FZI
            RHOAX(I110)=RHOAX(I110)+ABC001*FXI
            RHOAY(I110)=RHOAY(I110)+ABC001*FYI
            RHOAZ(I110)=RHOAZ(I110)+ABC001*FZI
            RHOAX(I101)=RHOAX(I101)+ABC010*FXI
            RHOAY(I101)=RHOAY(I101)+ABC010*FYI
            RHOAZ(I101)=RHOAZ(I101)+ABC010*FZI
            RHOAX(I011)=RHOAX(I011)+ABC100*FXI
            RHOAY(I011)=RHOAY(I011)+ABC100*FYI
            RHOAZ(I011)=RHOAZ(I011)+ABC100*FZI
            RHOAX(I111)=RHOAX(I111)+ABC000*FXI
            RHOAY(I111)=RHOAY(I111)+ABC000*FYI
            RHOAZ(I111)=RHOAZ(I111)+ABC000*FZI
    ENDIF  /* ISGSTA */
  enddo 
    !
#ifdef MPI
        if(numtasks > 1)then
!  combining all node results
!
          call mpi_allreduce(MPI_IN_PLACE,rhom,ngridxyz, &
          mpi_double_precision,mpi_sum,commsander,ierr)
          call mpi_allreduce(MPI_IN_PLACE,rhovx,ngridxyz, &
          mpi_double_precision,mpi_sum,commsander,ierr)
          call mpi_allreduce(MPI_IN_PLACE,rhovy,ngridxyz, &
          mpi_double_precision,mpi_sum,commsander,ierr)
          call mpi_allreduce(MPI_IN_PLACE,rhovz,ngridxyz, &
          mpi_double_precision,mpi_sum,commsander,ierr)
          call mpi_allreduce(MPI_IN_PLACE,rhoax,ngridxyz, &
          mpi_double_precision,mpi_sum,commsander,ierr)
          call mpi_allreduce(MPI_IN_PLACE,rhoay,ngridxyz, &
          mpi_double_precision,mpi_sum,commsander,ierr)
          call mpi_allreduce(MPI_IN_PLACE,rhoaz,ngridxyz, &
          mpi_double_precision,mpi_sum,commsander,ierr)
        endif
#endif
  ! remove net velocity and acceleration
    AMASSI=SUM(RHOM)
    GXI=SUM(RHOVX)/AMASSI
    GYI=SUM(RHOVY)/AMASSI
    GZI=SUM(RHOVZ)/AMASSI
    FXI=SUM(RHOAX)/AMASSI
    FYI=SUM(RHOAY)/AMASSI
    FZI=SUM(RHOAZ)/AMASSI
    DO I=1,NGRIDXYZ
      RHOMI=RHOM(I)
      IF(RHOMI>1.0d-6)THEN
        RHOVX(I)=RHOVX(I)/RHOMI
        RHOVY(I)=RHOVY(I)/RHOMI
        RHOVZ(I)=RHOVZ(I)/RHOMI
        RHOAX(I)=RHOAX(I)/RHOMI
        RHOAY(I)=RHOAY(I)/RHOMI
        RHOAZ(I)=RHOAZ(I)/RHOMI
      ENDIF
    ENDDO
    RHOVX=RHOVX-GXI
    RHOVY=RHOVY-GYI
    RHOVZ=RHOVZ-GZI
    RHOAX=RHOAX-FXI
    RHOAY=RHOAY-FYI
    RHOAZ=RHOAZ-FZI
    RETURN
  END SUBROUTINE MAPBUILD

  SUBROUTINE mapvalues(amassi,crdi,avgpi,avgfi)

  !-----------------------------------------------------------------------
  !     This routine to find map values at input position
  !     throuth linear interpolation
  !-----------------------------------------------------------------------

  implicit none
  double precision AMASSI,TMASS
  double precision crdi(3),avgpi(3),avgfi(3)
  INTEGER X0,Y0,Z0,X1,Y1,Z1
  INTEGER I000,I100,I010,I001,I110,I101,I011,I111
  double precision XI,YI,ZI,VXI,VYI,VZI,AXI,AYI,AZI
  double precision GX,GY,GZ,A0,B0,C0,A1,B1,C1
  double precision ABC000,ABC100,ABC010,ABC001,ABC110,ABC101,ABC011,ABC111
     GX=(crdi(1)-GXMIN)/SGSIZE+1.0d0
     GY=(crdi(2)-GYMIN)/SGSIZE+1.0d0
     GZ=(crdi(3)-GZMIN)/SGSIZE+1.0d0
     X0=INT(GX)
     Y0=INT(GY)
     Z0=INT(GZ)
     IF(X0<1.OR.Y0<1.OR.Z0<1)THEN
       STOP 'Position outside of lower boundary'
     ENDIF
     X1=X0+1
     Y1=Y0+1
     Z1=Z0+1
     IF(X1>NGRIDX.OR.Y1>NGRIDY.OR.Z1>NGRIDZ)THEN
       STOP 'Position outside of higher boundary'
     ENDIF
     A0=X1-GX
     B0=Y1-GY
     C0=Z1-GZ
     A1=1.0d0-A0
     B1=1.0d0-B0
     C1=1.0d0-C0
     ABC000=A0*B0*C0
     ABC100=A1*B0*C0
     ABC010=A0*B1*C0
     ABC001=A0*B0*C1
     ABC110=A1*B1*C0
     ABC011=A0*B1*C1
     ABC101=A1*B0*C1
     ABC111=A1*B1*C1
     I000=X0+NGRIDX*(Y0-1+NGRIDY*(Z0-1))
     I100=I000+1
     I010=I000+NGRIDX
     I001=I000+NGRIDXY
     I110=I010+1
     I101=I001+1
     I011=I001+NGRIDX
     I111=I011+1
     TMASS=ABC111*RHOM(I000)+ABC110*RHOM(I001)+ABC101*RHOM(I010)+ABC011*RHOM(I100)+     &
     ABC100*RHOM(I011)+ABC010*RHOM(I101)+ABC001*RHOM(I110)+ABC000*RHOM(I111)
     VXI=ABC111*RHOVX(I000)+ABC110*RHOVX(I001)+ABC101*RHOVX(I010)+ABC011*RHOVX(I100)+     &
     ABC100*RHOVX(I011)+ABC010*RHOVX(I101)+ABC001*RHOVX(I110)+ABC000*RHOVX(I111)
     VYI=ABC111*RHOVY(I000)+ABC110*RHOVY(I001)+ABC101*RHOVY(I010)+ABC011*RHOVY(I100)+     &
     ABC100*RHOVY(I011)+ABC010*RHOVY(I101)+ABC001*RHOVY(I110)+ABC000*RHOVY(I111)
     VZI=ABC111*RHOVZ(I000)+ABC110*RHOVZ(I001)+ABC101*RHOVZ(I010)+ABC011*RHOVZ(I100)+     &
     ABC100*RHOVZ(I011)+ABC010*RHOVZ(I101)+ABC001*RHOVZ(I110)+ABC000*RHOVZ(I111)
     AXI=ABC111*RHOAX(I000)+ABC110*RHOAX(I001)+ABC101*RHOAX(I010)+ABC011*RHOAX(I100)+     &
     ABC100*RHOAX(I011)+ABC010*RHOAX(I101)+ABC001*RHOAX(I110)+ABC000*RHOAX(I111)
     AYI=ABC111*RHOAY(I000)+ABC110*RHOAY(I001)+ABC101*RHOAY(I010)+ABC011*RHOAY(I100)+     &
     ABC100*RHOAY(I011)+ABC010*RHOAY(I101)+ABC001*RHOAY(I110)+ABC000*RHOAY(I111)
     AZI=ABC111*RHOAZ(I000)+ABC110*RHOAZ(I001)+ABC101*RHOAZ(I010)+ABC011*RHOAZ(I100)+     &
     ABC100*RHOAZ(I011)+ABC010*RHOAZ(I101)+ABC001*RHOAZ(I110)+ABC000*RHOAZ(I111)
     avgpi(1)=AMASSI*VXI
     avgpi(2)=AMASSI*VYI
     avgpi(3)=AMASSI*VZI
     avgfi(1)=AMASSI*AXI
     avgfi(2)=AMASSI*AYI
     avgfi(3)=AMASSI*AZI
     RETURN
  END SUBROUTINE MAPVALUES
  

  SUBROUTINE MAPINIT(ISTART,IEND,CRD)
  !-----------------------------------------------------------------------
  !     This routine build and initialize maps
  !
  !-----------------------------------------------------------------------

  implicit none
  integer ISTART,IEND
  double precision crd(3,*)
#ifdef MPI
     integer ierr
# include "parallel.h"
#endif

  integer I,IA,alloc_err,dealloc_err
  double precision XI,YI,ZI
  double precision XMAX,XMIN,YMAX,YMIN,ZMAX,ZMIN,TEMP1(6),TEMP2(6)
  logical rmap

  XMIN=1.0D8
  XMAX=-1.0D8
  YMIN=1.0D8
  YMAX=-1.0D8
  ZMIN=1.0D8
  ZMAX=-1.0D8
  do i = ISTART,IEND

       XI=CRD(1,I)
       YI=CRD(2,I)
       ZI=CRD(3,I)
       IF(XI>XMAX)XMAX=XI
       IF(XI<XMIN)XMIN=XI
       IF(YI>YMAX)YMAX=YI
       IF(YI<YMIN)YMIN=YI
       IF(ZI>ZMAX)ZMAX=ZI
       IF(ZI<ZMIN)ZMIN=ZI
  ENDDO
#ifdef MPI
        if(numtasks > 1)then
!  combining all node results
!
          temp1(1)=XMAX
          temp1(2)=-XMIN
          temp1(3)=YMAX
          temp1(4)=-YMIN
          temp1(5)=ZMAX
          temp1(6)=-ZMIN

          call mpi_allreduce(temp1,temp2,6, &
          mpi_double_precision,mpi_max,commsander,ierr)

          xmax=temp2(1)
          xmin=-temp2(2)
          ymax=temp2(3)
          ymin=-temp2(4)
          zmax=temp2(5)
          zmin=-temp2(6)
        endif
#endif
  RMAP=(XMAX>GXMAX.OR.XMIN<GXMIN.OR.YMAX>GYMAX.OR.YMIN<GYMIN.OR.ZMAX>GZMAX.OR.ZMIN<GZMIN)
  IF(RMAP)THEN

    GXMAX=SGSIZE*(AINT(XMAX/SGSIZE)+1.0d0)
    GXMIN=SGSIZE*(AINT(XMIN/SGSIZE)-1.0d0)
    GYMAX=SGSIZE*(AINT(YMAX/SGSIZE)+1.0d0)
    GYMIN=SGSIZE*(AINT(YMIN/SGSIZE)-1.0d0)
    GZMAX=SGSIZE*(AINT(ZMAX/SGSIZE)+1.0d0)
    GZMIN=SGSIZE*(AINT(ZMIN/SGSIZE)-1.0d0)
    NGRIDX=INT((GXMAX-GXMIN)/SGSIZE)+1
    NGRIDY=INT((GYMAX-GYMIN)/SGSIZE)+1
    NGRIDZ=INT((GZMAX-GZMIN)/SGSIZE)+1
    NGRIDXY=NGRIDX*NGRIDY
    NGRIDYZ=NGRIDY*NGRIDZ
    NGRIDZX=NGRIDZ*NGRIDX
    NGRIDXYZ=NGRIDX*NGRIDY*NGRIDZ
    ! deallocate maps
    if(allocated(rhom))deallocate(rhom)  ! deallocate sgmaps
    if(allocated(rhovx))deallocate(rhovx)  ! deallocate sgmaps
    if(allocated(rhovy))deallocate(rhovy)  ! deallocate sgmaps
    if(allocated(rhovz))deallocate(rhovz)  ! deallocate sgmaps
    if(allocated(rhoax))deallocate(rhoax)  ! deallocate sgmaps
    if(allocated(rhoay))deallocate(rhoay)  ! deallocate sgmaps
    if(allocated(rhoaz))deallocate(rhoaz)  ! deallocate sgmaps
    !if(allocated(rhoaz))deallocate(rhoaz,stat=dealloc_err)  ! deallocate sgmaps
    !if(dealloc_err /= 0 ) then
    !  stop "unable to deallocate SG maps "
    !endif
    allocate( rhom(ngridxyz),rhovx(ngridxyz),rhovy(ngridxyz),rhovz(ngridxyz),rhoax(ngridxyz),rhoay(ngridxyz),rhoaz(ngridxyz),stat=alloc_err)  ! allocate sgmaps
    !if (alloc_err .ne. 0) then
    !  stop "unable to allocate SG maps "
    !endif
  ENDIF
    RHOM=0.0d0
    RHOVX=0.0d0
    RHOVY=0.0d0
    RHOVZ=0.0d0
    RHOAX=0.0d0
    RHOAY=0.0d0
    RHOAZ=0.0d0
    RETURN
  END SUBROUTINE MAPINIT

end module sgld

