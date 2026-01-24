! <compile=optimized>
module emap

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  EMAP module  for flexible fitting macromolecules
!  into experimental maps
!         By Xiongwu Wu:  wuxw@nhlbi.nih.gov
!            updated 10/1/2011
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use, intrinsic :: iso_fortran_env
  implicit none

#include "../include/assert.fh"
#include "../include/dprec.fh"

public   

    TYPE EMAPOBJECT           ! map object type
      integer :: mapid              ! map object identity
      character(len=80) mapfile
      integer :: type               ! map object type
      integer :: MODE               ! map object mode
      character(LEN=80) :: name     ! map object name
      character(LEN=516) :: info    ! information about the map object
      _REAL_,pointer :: xatom(:)=> null(),yatom(:)=> null(),zatom(:)=> null()
      integer :: lx,ly,lz     !  map size in grid numbers (voxels)
      integer :: mx,my,mz     !  starting indexes (voxels)
      integer :: nx,ny,nz     !  map size for a unit cell (voxels)
      _REAL_    :: avg,std,max,min   ! statistics
      _REAL_    :: alpha,beta,gamma   ! lattice angles (degrees)
      _REAL_    :: cx,cy,cz  !  map center in the map coordinate (vexols)
      _REAL_    :: dx,dy,dz  !  grid sizes (angstroms/vexol)
      _REAL_    :: ux,uy,uz           !  unit cell sizes (angstroms)
      _REAL_    :: ox,oy,oz           !  absolute coordinates of map center (angstroms)
      _REAL_    :: sx,sy,sz           !  parameters for reduced coordinates (angstroms)
      integer(int32),pointer :: bdata(:)=> null()     !   integer data
      integer(int64),pointer :: idata(:)=> null()     !   integer data
      real(real32),pointer:: rdata(:)=> null()        !   real data
      complex(real64), pointer :: Cdata(:)=> null()  !   double precision data
      type (emapobject), pointer :: next   !  pointer to next emap object
      type (emapobject), pointer :: ref    !  pointer to a reference emap object
    END TYPE EMAPOBJECT

    TYPE EMAPRIGID        ! rigid domain type
      integer :: rigid      ! rigid domain identity
      integer :: mapid      ! rigid domain identity
      _REAL_ fcons
      _REAL_ :: RESOLUTION        !  resolution of the map object
      character(len=256) atmask
      character(LEN=80) :: name     ! rigid domain name
      character(LEN=80) :: mapfit     ! output map file name
      character(LEN=80) :: molfit     ! output molecular file name
      integer :: natc             !  number of atoms represented by the map
      integer, pointer, dimension(:) :: idxatom => null()  !  pointer to atom index array
      logical :: movable,fitting,tomap        ! flag indicating a fixed rigid domain
      _REAL_ :: mass,winv,tcm         ! mass, invert
      _REAL_ :: energy         ! energy
      INTEGER,dimension(6) :: grids   ! grid numbers in x,y,z,phi,psi,theta
      integer :: nminim             !  number of minimization steps
      _REAL_,dimension(3) :: prinv   ! inverse of principle inertia
      _REAL_,dimension(3) :: center         ! Position center coordinates 
      _REAL_,dimension(3) :: eular ! Eular angles:phi,theta,psi  
      _REAL_,dimension(4) :: quad   ! quarterions  
      _REAL_,dimension(3) :: v         ! velocities 
      _REAL_,dimension(3) :: omiga        ! angular velocities 
      _REAL_,dimension(3) :: angm        ! angular momentum 
      _REAL_,dimension(3) :: TRAN,TRANREF ! Translation vectors
      _REAL_,dimension(3,3) :: ROT,ROTREF ! rotational matrix
      _REAL_,dimension(3) :: force ! force on center of mass
      _REAL_,dimension(3) :: torq ! torque
      type (emapobject), pointer :: emap ! pointer to the associated emap object
      type (emaprigid), pointer :: next ! pointer to next rigid domain
      type (emaprigid), pointer :: ref  ! pointer to a reference rigid domain
      integer :: nc      ! number of constraint points
      real, pointer :: xc,yc,zc  ! coordinate array of constraint points
    END TYPE EMAPRIGID


    LOGICAL TEMAP
    INTEGER NEMAP,NRIGID,BORDER
    _REAL_ GAMMAMAP,EKMAP,MAPDT,XGAMMA,XAVG0,XAVG1,EMSIGMA,SCEMAP
    type(EMAPRIGID), save :: rigmin
    type(EMAPOBJECT), allocatable,dimension(:),save :: EMAPS
    type(EMAPRIGID), allocatable,dimension(:),save  :: EMRIGS

contains



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE EMAP_OPTIONS(UNIT)
!_________________________________________________________________
!  read in map constraint information
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
      INTEGER UNIT
      character(len=80) mapfile
      character(len=80) mapfit,molfit
      character(len=256) atmask
      _REAL_ fcons,resolution
      INTEGER move,ifit,minim,grids(6)
      namelist /emap/mapfile,mapfit,molfit,atmask,fcons,move,ifit,grids,minim,resolution
      INTEGER I,J,IFIND,mapid,rigid,alloc_err
      mapid=0
      nrigid=0
      nemap=0
!  search for namelist EMAP
      rewind UNIT
      DO
        call nmlsrc('emap',UNIT,ifind)
        if (ifind == 0) exit
        read(UNIT,nml=emap,err=999)
        nrigid=nrigid+1
      ENDDO
      if(allocated(EMRIGS))deallocate(EMRIGS)
      allocate(EMRIGS(NRIGID),stat=alloc_err)
      if(alloc_err /= 0 ) write(6,*)"unable to allocate EMRIGS"
      rewind UNIT
      DO I=1,NRIGID
        mapfile=''
        atmask=':*'
        mapfit=''
        molfit=''
        move=0
        ifit=0
        minim=100
        grids=1
        fcons=0.05d0
        resolution=0.0d0
        rigid=i
        read(UNIT,nml=emap,err=999)
        if(len_trim(mapfile)==0)then
          nemap=nemap+1
          mapid=nemap
        else
          DO J=1,I-1
            if(emrigs(j)%name==mapfile)then
              mapid=emrigs(j)%mapid
              exit
            endif
          ENDDO
          if(j==i)then
            nemap=nemap+1
            mapid=nemap
          endif
        endif
        call rigid_init(emrigs(i))
        emrigs(i)%name=mapfile
        emrigs(i)%mapfit=mapfit
        emrigs(i)%molfit=molfit
        emrigs(i)%rigid=rigid
        emrigs(i)%atmask=atmask
        emrigs(i)%mapid=mapid
        emrigs(i)%fcons=fcons
        emrigs(i)%resolution=resolution
        emrigs(i)%movable=move>0
        emrigs(i)%fitting=ifit>0
        emrigs(i)%tomap=ifit>1
        emrigs(i)%grids=grids
        emrigs(i)%nminim=minim
      ENDDO
      if(allocated(EMAPS))deallocate(EMAPS)
      allocate(EMAPS(NEMAP),stat=alloc_err)
      if(alloc_err /= 0 ) write(6,*)"unable to allocate EMAPS"
      J=0
      DO I=1,NRIGID
        if(emrigs(i)%mapid>j)then
          j=emrigs(i)%mapid
          emaps(j)%mapfile=emrigs(i)%name
        endif
      ENDDO
999   return
      end subroutine emap_options


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE PEMAP(dtmol,tempmol,x,ix,ih)
!_________________________________________________________________
!  Preprocessing EMAP module
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
      use constants, only: KB
      use findmask,only:atommask
      use md_scheme, only: gamma_ln
#ifdef MPI
      use mpi
#endif
      implicit none
#ifdef MPI
#  include "parallel.h"
#endif
#  include "../include/memory.h"
#  include "../include/md.h"
      _REAL_ dtmol,tempmol
      _REAL_ x(*)
      integer ix(*)
      character(len=4) :: ih(*)

      INTEGER I,J,nm,iat,i3,j3,idemp
      INTEGER alloc_err,UNITOUT
      _REAL_ RESO
      INTEGER, allocatable, dimension(:)::mask
      _REAL_, allocatable, dimension(:)::TMP_CRD,TMP_AMASS
#ifdef MPI
      INTEGER ierr
#endif
!
      BORDER=6
      UNITOUT=6
      !  map moving time step
      mapdt=20.455d0*dtmol
      ! Map moving kinetic energy
      EKMAP=1.5d0*KB*300.0d0
      IF(TEMPMOL>0.0d0)EKMAP=1.5d0*KB*TEMPMOL
      ! fraction of map and molecule frictions
      xgamma=gamma_ln/(gamma_ln+gammamap)
      XAVG1=0.01d0
      XAVG0=1.0d0-XAVG1
#ifdef MPI
      if ( sandersize > 1 ) then
        if ( sanderrank > 0 ) then
          UNITOUT=-1
          if(allocated(EMRIGS))deallocate(EMRIGS)
          allocate(EMRIGS(NRIGID),stat=alloc_err)
          if(alloc_err /= 0 ) write(6,*)"unable to allocate EMRIGS"
          if(allocated(EMAPS))deallocate(EMAPS)
          allocate(EMAPS(NEMAP),stat=alloc_err)
          if(alloc_err /= 0 ) write(6,*)"unable to allocate EMAPS"
        end if
        DO I=1,NRIGID
          if ( sanderrank > 0 )call rigid_init(emrigs(i))
          call mpi_bcast(emrigs(i)%mapid,1,MPI_INTEGER,0,commsander,ierr)
          call mpi_bcast(emrigs(i)%rigid,1,MPI_INTEGER,0,commsander,ierr)
          call mpi_bcast(emrigs(i)%name,80,MPI_CHARACTER,0,commsander,ierr)
          call mpi_bcast(emrigs(i)%mapfit,80,MPI_CHARACTER,0,commsander,ierr)
          call mpi_bcast(emrigs(i)%atmask,256,MPI_CHARACTER,0,commsander,ierr)
          call mpi_bcast(emrigs(i)%fcons,1,MPI_DOUBLE_PRECISION,0,commsander,ierr)
          call mpi_bcast(emrigs(i)%resolution,1,MPI_DOUBLE_PRECISION,0,commsander,ierr)
          call mpi_bcast(emrigs(i)%movable,1,MPI_LOGICAL,0,commsander,ierr)
          call mpi_bcast(emrigs(i)%fitting,1,MPI_LOGICAL,0,commsander,ierr)
          call mpi_bcast(emrigs(i)%nminim,1,MPI_INTEGER,0,commsander,ierr)
          call mpi_bcast(emrigs(i)%grids,6,MPI_INTEGER,0,commsander,ierr)
        ENDDO
        DO I=1,NEMAP
          call mpi_bcast(emaps(i)%mapfile,80,MPI_CHARACTER,0,commsander,ierr)
        ENDDO
      end if
#endif

      DO I=1,NEMAP
        CALL RDEMAP(EMAPS(I)%MAPFILE,I,UNITOUT)
      ENDDO
      if(allocated(mask))deallocate(mask)
      allocate(mask(NATOM),stat=alloc_err)
      DO I=1,NRIGID
        call atommask( natom, nres, 0, ih(m04), ih(m06), &
                ix(i02), ih(m02), x(lcrd), emrigs(i)%atmask, mask )
        nm=0
        DO J=1,NATOM
          if(mask(J)>0)nm=nm+1
        ENDDO
        emrigs(i)%natc=nm
        if(associated(EMRIGS(I)%IDXATOM))deallocate(EMRIGS(I)%IDXATOM)
        allocate(EMRIGS(I)%IDXATOM(NM),stat=alloc_err)
        if(alloc_err /= 0 ) write(6,*)"unable to allocate IDXATOM"
        nm=0
        DO J=1,NATOM
          if(mask(J)>0)then
            nm=nm+1
            emrigs(i)%idxatom(nm)=j
          endif
        ENDDO
        call rigid_prop(x(lmass),x(lcrd),emrigs(i))
        if(len_trim(emrigs(i)%name)==0)then
          ! create map from the coordinates of masked atoms
          if(allocated(TMP_CRD))deallocate(TMP_CRD)
          allocate(TMP_CRD(3*NM),stat=alloc_err)
          if(alloc_err /= 0 ) write(6,*)"unable to allocate TMP_CRD"
          if(allocated(TMP_AMASS))deallocate(TMP_AMASS)
          allocate(TMP_AMASS(NM),stat=alloc_err)
          if(alloc_err /= 0 ) write(6,*)"unable to allocate TMP_AMASS"
          DO IAT=1,NM
            J=emrigs(i)%idxatom(iat)
            tmp_amass(iat)=x(lmass+j-1)
            i3=3*iat-2
            j3=lcrd+3*j-3
            tmp_crd(i3)=x(j3)
            tmp_crd(i3+1)=x(j3+1)
            tmp_crd(i3+2)=x(j3+2)
          enddo
          IDEMP=emrigs(i)%mapid
          RESO=emrigs(i)%resolution
          if(RESO<0.1D0)RESO=2.0D0
          if(unitout>0)then
            write(unitout,'("Map ",I4," is created from ",i6," constrained atoms with resolution: ",F4.1)') &
              IDEMP,nm, RESO
          endif
          CALL COR2MAP(IDEMP,RESO,nm,tmp_amass,tmp_crd,unitout)      
        endif
        if(unitout>0)then
           write(unitout,'("Rigid ",I4," has ",i6," constrained atoms with mask: ",A)') &
              I,nm, emrigs(i)%atmask
           if(emrigs(i)%fitting)then
             write(unitout,'("    Rigid-fitting is requested for the constrained atoms with a grid:")')
             write(unitout,'("      x: ",I2," y: ",I2," z: ",I2," theta: ",I2," phi: ",I2," psi: ",I2)')(emrigs(i)%grids(j),j=1,6)
             if(emrigs(i)%tomap) &
             write(unitout,'("    Atoms for Rigid ",I4," will be moved to fit the map")')I
           else
             write(unitout,'("    Rigid-fitting is not requested for the constrained atoms.")')
           endif
           if(emrigs(i)%movable)then
             write(unitout,'("    The constraint map is movable to follow the constrained atoms.")')
           else
             write(unitout,'("    The constraint map is fixed in space.")')
           endif
        endif
        ! perform map to atom fit if required
        if(emrigs(i)%fitting)then
           j=emrigs(i)%mapid
           call emapfit(emrigs(i),emaps(j),UNITOUT)
           if(unitout>0)then
              write(unitout,'("Rigid ",I4," has been fitted to: ",A," with energy:",E14.6)') &
              I, trim(emaps(j)%mapfile),emrigs(i)%energy
           endif
           if(emrigs(i)%tomap)then
             ! transform the atoms for this rigid dormain
             CALL TRANSCRD(emrigs(i)%NATC,emrigs(i)%IDXATOM,EMAPS(J)%CX,&
            EMAPS(J)%CY,EMAPS(J)%CZ,emrigs(i)%TRAN,emrigs(i)%ROT)
             call rigid_init(emrigs(i))
           endif
        endif
      ENDDO
      if(allocated(MASK))deallocate(mask)
      if(allocated(TMP_AMASS))deallocate(tmp_amass)
      if(allocated(TMP_CRD))deallocate(tmp_crd)
      return
      end subroutine pemap

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE QEMAP()
!_________________________________________________________________
!  Clean and quit EMAP module
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
#ifdef MPI
      use mpi
#endif
      implicit none
#ifdef MPI
#  include "parallel.h"
#endif
#  include "../include/memory.h"

      INTEGER I,dealloc_err,UNITOUT
!
      UNITOUT=6
#ifdef MPI
      if ( sanderrank > 0 ) then
        UNITOUT=-1
      end if
#endif

      DO I=1,NRIGID
        if(UNITOUT>0)call RIGMAP(emrigs(i),UNITOUT)
        if(associated(EMRIGS(I)%IDXATOM))deallocate(EMRIGS(I)%IDXATOM,stat=dealloc_err)
        if(dealloc_err /= 0 ) write(6,*)"unable to deallocate IDXATOM for RIGID ",I
      ENDDO
      if(allocated(EMRIGS))deallocate(EMRIGS)
      if(allocated(EMAPS))deallocate(EMAPS)
      return
      end subroutine qemap

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE EMAPFIT(RIGOBJ,MAPOBJ,OUTU)
!_________________________________________________________________
!  fit map to a rigid with grid-threading minimization
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
   use constants, only: PI,TWOPI
      type(EMAPRIGID) :: rigobj
      type(EMAPOBJECT) :: mapobj
      INTEGER OUTU

#  include "../include/memory.h"
#  include "../include/md.h"
      INTEGER I
      INTEGER ktx,kty,ktz,krx,kry,krz
      INTEGER it,jt,kt,ir,jr,kr,iter
      _REAL_ dtx,gx,gy,gz,dphi,dpsi,dtht
      _REAL_ xl,yl,zl,pos0(6),pos(6),posi(6),gtol
      _REAL_ p(7,6),engs(7)
!
      dtx = dt*20.455d+00
      gtol=1.0d-20
!
      ktx=rigobj%grids(1)
      kty=rigobj%grids(2)
      ktz=rigobj%grids(3)
      krx=rigobj%grids(4)
      kry=rigobj%grids(5)
      krz=rigobj%grids(6)
      xl=mapobj%dx*mapobj%mx
      yl=mapobj%dy*mapobj%my
      zl=mapobj%dz*mapobj%mz
      gx=mapobj%dx*mapobj%lx/(ktx)
      gy=mapobj%dy*mapobj%ly/(kty)
      gz=mapobj%dz*mapobj%lz/(ktz)
      dphi=TWOPI/krx
      dpsi=TWOPI/kry
      dtht=PI/krz
      rigmin=rigobj
      pos0(1)=rigobj%tran(1)+anint((rigobj%center(1)-rigobj%tran(1)-xl)/gx-0.5d0)*gx
      pos0(2)=rigobj%tran(2)+anint((rigobj%center(2)-rigobj%tran(2)-yl)/gy-0.5d0)*gy
      pos0(3)=rigobj%tran(3)+anint((rigobj%center(3)-rigobj%tran(3)-zl)/gz-0.5d0)*gz
      pos0(4:6)=rigobj%eular
      do it=1,ktx
        pos(1)=pos0(1)-(it-1)*gx
        do jt=1,kty
          pos(2)=pos0(2)-(jt-1)*gy
          do kt=1,ktz
            pos(3)=pos0(3)-(kt-1)*gz
            do ir=1,krx
              pos(4)=pos0(4)+(ir-1)*dphi
              do jr=1,kry
                pos(5)=pos0(5)+(jr-1)*dpsi
                do kr=1,krz
                  pos(6)=pos0(6)+(kr-1)*dtht
                  !call dfpmin(pos,6,gtol,iter,enrig,rigeng)
                  do i=1,7
                    p(i,1:6)=pos
                  enddo                  
                  p(1,1)=p(1,1)+gx/2
                  p(2,2)=p(2,2)+gy/2
                  p(3,3)=p(3,3)+gz/2
                  p(4,4)=p(4,4)+dphi/2
                  p(5,5)=p(5,5)+dpsi/2
                  p(6,6)=p(6,6)+dtht/2
                  do i=1,7
                    posi=p(i,1:6)
                    engs(i)=rigeng(posi)
                  enddo
                  if(outu>0)write(outu,'("initial:",6f8.2,e14.6)')(p(7,i),i=1,6),engs(7)
                  call simplex(p,engs,7,6,6,gtol,rigeng,iter)  
                  if(outu>0)write(outu,'("result: ",6f8.2,e14.6)')(p(1,i),i=1,6),engs(1)
                  if(rigmin%energy<rigobj%energy)then
                     rigobj=rigmin
                 endif
                enddo  !kr
              enddo  !jr
            enddo  !ir
          enddo  !kt
        enddo  !jt
      enddo  !it
      RETURN
      END SUBROUTINE EMAPFIT


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE rigid_init(rigobj)
!_________________________________________________________________
!  initialize rigid domain
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
      type(EMAPRIGID) :: rigobj

      INTEGER I,J

      rigobj%energy=1.0D9
      DO I=1,3
        DO J=1,3
          rigobj%ROT(I,J)=0.0D0
        ENDDO
        rigobj%ROT(I,I)=1.0D0
        rigobj%TRAN(I)=0.0D0
        rigobj%center(I)=0.0D0
        rigobj%FORCE(I)=0.0D0
        rigobj%TORQ(I)=0.0D0
        rigobj%QUAD(I+1)=0.0D0
      ENDDO
      rigobj%QUAD(1)=1.0D0
      return
      end subroutine rigid_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE rigid_prop(amass,crd,rigobj)
!_________________________________________________________________
!  initialize rigid domain
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
      _REAL_ amass(*),CRD(*)
      type(EMAPRIGID) :: rigobj

      INTEGER I,N,IA
      _REAL_ AMASSI,X,Y,Z,XC,YC,ZC,TCM,TMASS

      N=rigobj%natc

      XC=0.0D0
      YC=0.0D0
      ZC=0.0D0
      DO I=1,N
        IA=rigobj%idxatom(i)
        AMASSI=AMASS(IA)
        X=crd(3*IA-2)
        Y=crd(3*IA-1)
        Z=crd(3*IA)
        XC=XC+X
        YC=YC+Y
        ZC=ZC+Z
      ENDDO
      XC=XC/N
      YC=YC/N
      ZC=ZC/N
      TCM=0.0D0
      TMASS=0.0D0
      DO I=1,N
        IA=rigobj%idxatom(i)
        AMASSI=AMASS(IA)
        X=crd(3*IA-2)-XC
        Y=crd(3*IA-1)-YC
        Z=crd(3*IA)-ZC
        TCM=TCM+AMASSI*(X*X+Y*Y+Z*Z)
        TMASS=TMASS+AMASSI
      ENDDO
      rigobj%mass=TMASS
      rigobj%tcm=TCM
      rigobj%center(1)=xc
      rigobj%center(2)=yc
      rigobj%center(3)=zc
      return
      end subroutine rigid_prop

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE EMAP_move()
!_________________________________________________________________
!  Perform motion step for all rigid domains
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
      _REAL_ enemap

      INTEGER I,J !,OUTU
      enemap=0.0d0
      DO I=1,NRIGID
        enemap=enemap+emrigs(i)%energy
        if(emrigs(i)%movable)call rigid_move(emrigs(i))
        DO J=1,3
          emrigs(i)%force(j)=0.0d0
          emrigs(i)%torq(j)=0.0d0
        ENDDO
        emrigs(i)%energy=0.0d0
      ENDDO
      !IF(OUTU>0)WRITE(OUTU,'("EMAP = ",F14.6)')ENEMAP
      return
      end subroutine EMAP_move

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE rigid_move(rigobj)
!_________________________________________________________________
!  moving rigid domain based on position Langevin motion
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!
!     -dEmap/dr=f=gamma*m*v    v=f/(gamma*m)
!     -dEmap/da=T=gamma*I*w    w=T/(gamma*I)
!
!     v0=sqrt(3kT/m)   w0=sqrt(3kT/I)
!
!  Scale down based on the thermo motion
!     dr=dt*v/(1+sqrt(v/v0))
!     da=dt*w/(1+sqrt(w/w0))
!
!_______________________________________________________________________
!
      use constants,only: PI,TWOPI
      
      type(EMAPRIGID) :: rigobj

      INTEGER I
      _REAL_  ff,tt,scalf,scalt,ek
      _REAL_  DRF,DGT
      _REAL_ a(3,3),b(3,3)

      ff=dot_product(rigobj%force,rigobj%force)
      tt=dot_product(rigobj%torq,rigobj%torq)
      scalf=gammamap*gammamap*rigobj%mass/20.455d0/20.455d0
      scalt=gammamap*gammamap*rigobj%tcm/20.455d0/20.455d0
      ! Map moving kinetic energy
      EK=2.0d0*EKMAP*rigobj%natc
      DRF=mapdt/sqrt(rigobj%mass*(scalf+ff/ek))
      DGT=mapdt*sqrt(tt/(rigobj%tcm*(scalt+tt/ek)))

      ! rotation matrix due to the torque
      call rotaxisphi(a,rigobj%torq,DGT)
      
      ! update total rotation matrix
      b= matmul(a,rigobj%rot)
      
      ! calculate quaterion
      rigobj%quad=rot2quad(b)
       ! Update rotation matrix
      rigobj%rot=quad2rot(rigobj%quad)

       ! Cacluate Eular angles
      call rot2eular(rigobj%rot,rigobj%eular(1),rigobj%eular(2),rigobj%eular(3))
      
      DO I=1,3
        rigobj%TRAN(I)=rigobj%TRAN(I)+rigobj%FORCE(I)*DRF
      ENDDO
      return
      end subroutine rigid_move

      pure FUNCTION rot2quad(u)
!_______________________________________________________________________
!     This routine convert rotation matrix to quaterion
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_______________________________________________________________________
    !     quaterion: q(4)  
    !     matrix: u
    !
    !     q4=sqrt(1+u11+u22+u33)/2
    !     q1=(u32-u23)/4/q4
    !     q2=(u13-u31)/4/q4
    !     q3=(u21-u12)/4/q4
    !
    !
      implicit none
    !
      _REAL_, intent(in) :: U(3,3)
      _REAL_ rot2quad(4),q1,q2,q3,q4,qsum
    !
      q4=2.0d0*sqrt(1.0d0+U(1,1)+U(2,2)+U(3,3))
      q1=(U(3,2)-U(2,3))/q4
      q2=(U(1,3)-U(3,1))/q4
      q3=(U(2,1)-U(1,2))/q4
      q4=q4/4.0d0
      qsum=sqrt(q1*q1+q2*q2+q3*q3+q4*q4)
      rot2quad(1)=q1/qsum
      rot2quad(2)=q2/qsum
      rot2quad(3)=q3/qsum
      rot2quad(4)=q4/qsum
      END FUNCTION ROT2QUAD

      pure FUNCTION quad2rot(q)
!_______________________________________________________________________
!     This routine convert quaterion to rotation matrix 
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_______________________________________________________________________
    !     quaterion: q(4)  
    !     matrix: u
    !    
    !       { 1-2(q2^2+q3^2)      2(q1q2-q3q4)      2(q1q3+q2q4) }
    !     U={ 2(q1q2+q3)q4        1-2(q1^2+q3^2)    2(q2q3-q1q4)    }
    !       { 2(q1q3-q2q4)        2(q1q4+q2q3)      1-2(q1^2+q2^2)   }
    !
      implicit none
    !
      _REAL_ , intent(in) :: q(4)
      _REAL_ quad2rot(3,3)
    !
      quad2rot(1,1)=2.0d0*(q(1)*q(1)+q(4)*q(4))-1.0d0
      quad2rot(1,2)=2.0d0*(q(1)*q(2)-q(3)*q(4))
      quad2rot(1,3)=2.0d0*(q(1)*q(3)+q(2)*q(4))
      quad2rot(2,1)=2.0d0*(q(1)*q(2)+q(3)*q(4))
      quad2rot(2,2)=2.0d0*(q(2)*q(2)+q(4)*q(4))-1.0d0
      quad2rot(2,3)=2.0d0*(q(2)*q(3)-q(1)*q(4))
      quad2rot(3,1)=2.0d0*(q(1)*q(3)-q(2)*q(4))
      quad2rot(3,2)=2.0d0*(q(2)*q(3)+q(1)*q(4))
      quad2rot(3,3)=2.0d0*(q(3)*q(3)+q(4)*q(4))-1.0d0
      END FUNCTION QUAD2ROT

      SUBROUTINE eularrot(U,PHI,PSI,THT)
!_______________________________________________________________________
!     This routine convert Eular angles (z,x,z)
!     to rotational matrix
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!     emap energy and derivatives to x,y,z, and phi,psi,theta
!_______________________________________________________________________
    !     axis: ax(3)  
    !     angle: phir in rad
    !     matrix: u
    !     lOK: validate
    !
    !       { cosFcosY-sinFsinYcosS   sinFcosY+cosFsinYcosS   sinSsinY  }
    !     u={ -cosFsinY-sinFcosYcosS  -sinFsinY+cosFcosYcosS   sinScosY  }
    !       { sinFsinS                 -cosFsinS                cosS     }
    !
    !
    !
      implicit none
    !
      _REAL_ U(3,3),PHI,PSI,THT
    !
      _REAL_ sinF,cosF,sinY,cosY,sinS,cosS
    !
      sinF=sin(phi)
      cosF=cos(phi)
      sinY=sin(psi)
      cosY=cos(psi)
      sinS=sin(tht)
      cosS=cos(tht)
      u(1,1)=cosF*cosY-sinF*sinY*cosS 
      u(1,2)=sinF*cosY+cosF*sinY*cosS
      u(1,3)=sinS*sinY
      u(2,1)=-cosF*sinY-sinF*cosY*cosS
      u(2,2)=-sinF*sinY+cosF*cosY*cosS
      u(2,3)=sinS*cosY 
      u(3,1)=sinF*sinS  
      u(3,2)=-cosF*sinS
      u(3,3)=cosS
      RETURN
      END SUBROUTINE EULARROT


      SUBROUTINE rot2eular(U,PHI,PSI,THT)
!_______________________________________________________________________
!     This routine convert Eular angles
!     to rotational matrix
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!     emap energy and derivatives to x,y,z, and phi,psi,theta
!_______________________________________________________________________
    !     axis: ax(3)  
    !     angle: phir in rad
    !     matrix: u
    !     lOK: validate
    !
    !       { cosFcosY-sinFsinYcosS   sinFcosY+cosFsinYcosS   sinSsinY  }
    !     u={ -cosFsinY-sinFcosYcosS  -sinFsinY+cosFcosYcosS   sinScosY  }
    !       { sinFsinS                 -cosFsinS                cosS     }
    !
    !
    !
      implicit none
    !
      _REAL_ U(3,3),PHI,PSI,THT
    !
      phi=datan2(U(3,1),-U(3,2))
      tht=dacos(U(3,3))
      psi=datan2(U(1,3),U(2,3))
      RETURN
      END SUBROUTINE ROT2EULAR
      
      SUBROUTINE eulardrot(A,PHI,PSI,THT,daf,day,das)
!_______________________________________________________________________
!     This routine convert Eular angles
!     to rotational matrix and it derivative to Eular angles
!
!       ( cosf cosy-sinf siny coss,  sinf cosy+cosf siny coss,   siny sins  )
!     A=(-cosf siny-sinf cosy coss, -sinf siny+cosf cosy coss,   cosy sins  )
!       (   sinf sins             ,    -cosf sins            ,    coss      )
!
!       (-sinf cosy-cosf siny coss,  cosf cosy-sinf siny coss,   0  ) (-a12,a11,0)
! dA/df=( sinf siny-cosf cosy coss, -cosf siny-sinf cosy coss,   0  )=(-a22,a21,0)
!       (   cosf sins             ,    sinf sins            ,    0  ) (-a32,a31,0)
!
!       (-cosf siny-sinf cosy coss,-sinf siny+cosf cosy coss, cosy sins) ( a21, a22, a23)
! dA/dy=(-cosf cosy+sinf siny coss,-sinf cosy-cosf siny coss,-siny sins)=(-a11,-a12,-a31)
!       (   0                     ,    0                    ,   0      )=(  0 ,  0 ,  0 )
!
!       ( sinf siny sins,  -cosf siny sins,   siny coss  )    (a13a31,a32a13,a13a33   )
! dA/ds=( sinf cosy sins, -cosf cosy sins,    cosy coss  )=   (a31a23,a32a23,a23a33   )/b33
!       (   sinf coss   ,    -cosf coss   ,   -sins      )    (a31a33,a32a33,a33*a33-1)
!    here b33=sqrt(1-a33*a33)
!_________________________________________________________________
    !
    !
      implicit none
    !
      _REAL_ PHI,PSI,THT
    !
      _REAL_ a(3,3),daf(3,3),day(3,3),das(3,3)
      _REAL_ sinf,cosf,siny,cosy,sins,coss
    !
      sinf=sin(phi)
      siny=sin(psi)
      sins=sin(tht)
      cosf=cos(phi)
      cosy=cos(psi)
      coss=cos(tht)
      a(1,1)=cosF*cosY-sinF*sinY*cosS 
      a(1,2)=sinF*cosY+cosF*sinY*cosS
      a(1,3)=sinS*sinY
      a(2,1)=-cosF*sinY-sinF*cosY*cosS
      a(2,2)=-sinF*sinY+cosF*cosY*cosS
      a(2,3)=sinS*cosY 
      a(3,1)=sinF*sinS  
      a(3,2)=-cosF*sinS
      a(3,3)=cosS
      daf=0.0d0
      daf(1,1)=-a(1,2)
      daf(1,2)=a(1,1)
      daf(2,1)=-a(2,2)
      daf(2,2)=a(2,1)
      daf(3,1)=-a(3,2)
      daf(3,2)=a(3,1)
      day=0.0d0
      day(1,1)=a(2,1)
      day(1,2)=a(2,2)
      day(1,3)=a(2,3)
      day(2,1)=-a(1,1)
      day(2,2)=-a(1,2)
      day(2,3)=-a(1,3)
      das(1,1)=sinf*siny*sins
      das(1,2)=-cosf*siny*sins
      das(1,3)=siny*coss 
      das(2,1)=sinf*cosy*sins 
      das(2,2)=-cosf*cosy*sins 
      das(2,3)=cosy*coss 
      das(3,1)=sinf*coss 
      das(3,2)=-cosf*coss 
      das(3,3)= -sins 
      RETURN
      END SUBROUTINE EULARDROT


      SUBROUTINE rotaxisphi(U,RN,PHIR)
    !-------------  --------------------------------------------------------
    !     This routine convert rotation axis and angle
    !     to rotational matrix
    !     axis: ax(3)  
    !     angle: phir in rad
    !     matrix: u
    !     lOK: validate
    !
    !       { cos+ux2(1-cos)   uxuy(1-cos)-uzsin   uxuz(1-cos)+uysin   }
    !     u={uyux(1-cos)+uzsin   cos+uy2(1-cos)   uyuz(1-cos)-uxsin    }
    !       {uzux(1-cos)-uysin   uzuy(1-cos)+uxsin   cos+uz2(1-cos)    }
    !
    !     Xiongwu Wu
    !
    !
      implicit none
    !
      _REAL_ U(3,3),RN(3),PHIR
    !
      _REAL_ rnorm,rx,ry,rz,cosp,cosp1,sinp
    !
      cosp=cos(phir)
      cosp1=1.0D0-cosp
      sinp=sin(phir)
      rnorm=sqrt(rn(1)*rn(1)+rn(2)*rn(2)+rn(3)*rn(3))
      if(rnorm<1.0d-8)then
        u(1,1)=1.0d0
        u(1,2)=0.0d0
        u(1,3)=0.0d0
        u(2,1)=0.0d0
        u(2,2)=1.0d0
        u(2,3)=0.0d0
        u(3,1)=0.0d0
        u(3,2)=0.0d0
        u(3,3)=1.0d0
        return
      endif
      rx=rn(1)/rnorm
      ry=rn(2)/rnorm
      rz=rn(3)/rnorm
      u(1,1)=cosp+rx*rx*cosp1
      u(1,2)=rx*ry*cosp1-rz*sinp
      u(1,3)=rx*rz*cosp1+ry*sinp
      u(2,1)=rx*ry*cosp1+rz*sinp
      u(2,2)=cosp+ry*ry*cosp1
      u(2,3)=ry*rz*cosp1-rx*sinp
      u(3,1)=rx*rz*cosp1-ry*sinp
      u(3,2)=ry*rz*cosp1+rx*sinp
      u(3,3)=cosp+rz*rz*cosp1
      RETURN
      END SUBROUTINE rotaxisphi



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE RIGMAP(RIGOBJ,OUTU)
!_________________________________________________________________
!  write out map and molecule of a rigid 
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
      type(EMAPRIGID) :: rigobj
      type(EMAPOBJECT) :: mapobj
      INTEGER OUTU

      INTEGER IDEMP,UNIT
      CHARACTER*80 FNAME
!
      IF(rigobj%mapfit/='')then
        FNAME=rigobj%mapfit
        IDEMP=rigobj%mapid
        !  Transform map to rigid position
        IF(INDEX(FNAME,'.ccp4')>0.or.INDEX(FNAME,'.CCP4')>0 .or. &
           INDEX(FNAME,'.map')>0.or.INDEX(FNAME,'.MAP')>0 .or. &
           INDEX(FNAME,'.mrc')>0.or.INDEX(FNAME,'.MRC')>0)THEN
           CALL MAPCAST(MAPOBJ,EMAPS(IDEMP),RIGOBJ%TRAN,RIGOBJ%ROT,0.0D0)
           CALL WRTCCP4(FNAME,MAPOBJ,OUTU)      
           deallocate(MAPOBJ%RDATA)
        ELSE
          write(outu,*)' Map format is not supported. Skip output:',fname
        ENDIF
      ENDIF
      IF(rigobj%molfit/='')then
        FNAME=rigobj%molfit
        IDEMP=rigobj%mapid
        IF(INDEX(FNAME,'.pdb')>0.or.INDEX(FNAME,'.PDB')>0 )THEN
          write(outu,'(/" Fit structure is written to: ",A/)')fname
          UNIT=99
          OPEN(UNIT,FILE=FNAME,STATUS='UNKNOWN',FORM='FORMATTED')
          CALL WRTFITPDB(UNIT,RIGOBJ%NATC,RIGOBJ%IDXATOM,EMAPS(IDEMP)%CX,&
                         EMAPS(IDEMP)%CY,EMAPS(IDEMP)%CZ,RIGOBJ%TRAN,RIGOBJ%ROT)
          CLOSE(UNIT)
        ELSE
          write(outu,*)' structure format is not supported. Skip output:',fname
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE RIGMAP

      SUBROUTINE MAPCAST(MAPN,MAP,T,U,RESO)
!_________________________________________________________________
!  cast map according transform
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
      use ew_bspline,only:fill_bspline_0
!
      type(EMAPOBJECT) :: mapn,map
      _REAL_  T(3),U(3,3),RESO
!
      INTEGER*4 LX,LY,LZ,MODE,MNX,MNY,MNZ
      _REAL_ CX,CY,CZ,DX,DY,DZ,ALPHA,BETA,GAMMA
      INTEGER*4 MAPC,MAPR,MAPS
      REAL*4 AMAX,AMIN,AMEAN,ARMS
      INTEGER*4 ISPG,NSYMBT,LSKFLG,NNOTE
      CHARACTER*4 MAPLABEL,MACHST
      CHARACTER*80 NOTES(10)
      INTEGER*4 LXN,LYN,LZN,MXN,MYN,MZN,NDATA
      _REAL_ CXN,CYN,CZN,DXN,DYN,DZN
!
      INTEGER I,J,K,IN,JN,KN,M1,M2,M3,alloc_err
      INTEGER IPT1,IPT2,IPT3,ITH1,ITH2,ITH3,MORDER
      _REAL_ X,Y,Z,XN,YN,ZN,XM,YM,ZM
      _REAL_ FR1,FR2,FR3,W,RHOI
      REAL*8 BS1(8),BS2(8),BS3(8)
      REAL*8 VAL0,VAL0A
!
      LX=MAP%LX
      LY=MAP%LY
      LZ=MAP%LZ
      MODE=MAP%MODE
      MNX=MAP%MX
      MNY=MAP%MY
      MNZ=MAP%MZ
      DX=MAP%DX
      DY=MAP%DY
      DZ=MAP%DZ
      CX=MAP%CX
      CY=MAP%CY
      CZ=MAP%CZ
      XN=LX*DX
      YN=LY*DY
      ZN=LZ*DZ
      XM=ABS(U(1,1)*XN+U(1,2)*YN+U(1,3)*ZN)
      YM=ABS(U(2,1)*XN+U(2,2)*YN+U(2,3)*ZN)
      ZM=ABS(U(3,1)*XN+U(3,2)*YN+U(3,3)*ZN)
      X=ABS(U(1,1)*XN-U(1,2)*YN+U(1,3)*ZN)
      Y=ABS(U(2,1)*XN-U(2,2)*YN+U(2,3)*ZN)
      Z=ABS(U(3,1)*XN-U(3,2)*YN+U(3,3)*ZN)
      IF(XM<X)XM=X
      IF(YM<Y)YM=Y
      IF(ZM<Z)ZM=Z
      X=ABS(U(1,1)*XN+U(1,2)*YN-U(1,3)*ZN)
      Y=ABS(U(2,1)*XN+U(2,2)*YN-U(2,3)*ZN)
      Z=ABS(U(3,1)*XN+U(3,2)*YN-U(3,3)*ZN)
      IF(XM<X)XM=X
      IF(YM<Y)YM=Y
      IF(ZM<Z)ZM=Z
      X=ABS(U(1,1)*XN-U(1,2)*YN-U(1,3)*ZN)
      Y=ABS(U(2,1)*XN-U(2,2)*YN-U(2,3)*ZN)
      Z=ABS(U(3,1)*XN-U(3,2)*YN-U(3,3)*ZN)
      IF(XM<X)XM=X
      IF(YM<Y)YM=Y
      IF(ZM<Z)ZM=Z
      IF(RESO>0.0D0)THEN
        MORDER=BORDER
        DXN=2.0D0*RESO/MORDER
        DYN=2.0D0*RESO/MORDER
        DZN=2.0D0*RESO/MORDER
        DXN=DX
        DYN=DY
        DZN=DZ
      ELSE
        MORDER=2
        DXN=DX
        DYN=DY
        DZN=DZ
      ENDIF
      CXN=CX+T(1)
      CYN=CY+T(2)
      CZN=CZ+T(3)
      MXN=NINT((CXN-XM/2.0)/DXN)
      MYN=NINT((CYN-YM/2.0)/DYN)
      MZN=NINT((CZN-ZM/2.0)/DZN)
      LXN=NINT((XM)/DXN)
      LYN=NINT((YM)/DYN)
      LZN=NINT((ZM)/DZN)
      MAPN=MAP
      NDATA=LXN*LYN*LZN
      allocate(MAPN%RDATA(NDATA),stat=alloc_err)
      if(alloc_err /= 0 ) write(6,*)"unable to allocate RDATA"
      MAPN%LX=LXN
      MAPN%LY=LYN
      MAPN%LZ=LZN
      MAPN%MX=MXN
      MAPN%MY=MYN
      MAPN%MZ=MZN
      MAPN%NX=NINT(MAPN%NX*DX/DXN)
      MAPN%NY=NINT(MAPN%NY*DY/DYN)
      MAPN%NZ=NINT(MAPN%NZ*DZ/DZN)
      MAPN%CX=CXN
      MAPN%CY=CYN
      MAPN%CZ=CZN
      MAPN%DX=DXN
      MAPN%DY=DYN
      MAPN%DZ=DZN
      MAPN%MODE=2
      ALPHA=MAP%ALPHA
      BETA=MAP%BETA
      GAMMA=MAP%GAMMA
      MAPC=1
      MAPR=2
      MAPS=3
      AMIN=REAL(MAP%MAX)
      AMAX=REAL(MAP%MIN)
      AMEAN=0.0D0
      ARMS=0.0D0
      ISPG=0
      NSYMBT=0
      LSKFLG=0
      MAPLABEL='MAP '
      MACHST='ALL '
      NNOTE=3
      NOTES(1)=" This map is created with the emap module "
      NOTES(2)=" Report questions to Dr. Xiongwu Wu  "
      NOTES(3)="             Email: wuxw@nhlbi.nih.gov "
      DO KN=1,LZN
        ZN=(MZN+KN-1)*DZN-CZN
        DO JN=1,LYN
          YN=(MYN+JN-1)*DYN-CYN
          DO IN=1,LXN
            XN=(MXN+IN-1)*DXN-CXN
            X=U(1,1)*XN+U(2,1)*YN+U(3,1)*ZN+CX
            Y=U(1,2)*XN+U(2,2)*YN+U(3,2)*ZN+CY
            Z=U(1,3)*XN+U(2,3)*YN+U(3,3)*ZN+CZ
! calculate B-spline parameters
            FR1=(X)/DX
            FR2=(Y)/DY
            FR3=(Z)/DZ
            M1=INT(ANINT(FR1-0.5D0))
            W=FR1-M1
            call fill_bspline_0(w,MORDER,BS1)
            M2=INT(ANINT(FR2-0.5D0))
            W=FR2-M2
            call fill_bspline_0(w,MORDER,BS2)
            M3=INT(ANINT(FR3-0.5D0))
            W=FR3-M3
            call fill_bspline_0(w,MORDER,BS3)
! Calculate B-spline interception
            RHOI=0.0D0
            K = M3-MNZ+1 - MORDER/2
            DO 100 ITH3 = 1,MORDER
               K=K+1
               IF(K<1)THEN
                 IPT1=0
               ELSE IF(K>LZ)THEN
                 IPT1=LX*LY*(LZ-1)
               ELSE
                 IPT1=(K-1)*LY*LX
               ENDIF
               VAL0A =  BS3(ITH3)
!
               J = M2-MNY+1 - MORDER/2
!
               DO 200 ITH2 = 1,MORDER
                 J=J+1
                 IF(J<1)THEN
                   IPT2=IPT1
                 ELSE IF(J>LY)THEN
                   IPT2=IPT1+LX*(LY-1)
                 ELSE
                   IPT2=IPT1+(J-1)*LX
                 ENDIF
!
                 VAL0= VAL0A * BS2(ITH2)
!
                 I = M1-MNX+1 - MORDER/2
!
                 DO 300 ITH1 = 1,MORDER
                   I=I+1
                   IF(I<1)THEN
                     IPT3=IPT2+1
                   ELSE IF(I>LX)THEN
                     IPT3=IPT2+LX
                   ELSE
                     IPT3=IPT2+I
                   ENDIF
                   RHOI=RHOI+ VAL0 * BS1(ITH1)*MAP%RDATA(IPT3)
!
300              CONTINUE
200            CONTINUE
100         CONTINUE
!
            MAPN%RDATA(IN+LXN*(JN-1+LYN*(KN-1)))=REAL(RHOI)
            IF(RHOI<AMIN)AMIN=REAL(RHOI)
            IF(RHOI>AMAX)AMAX=REAL(RHOI)
            AMEAN=AMEAN+REAL(RHOI)
            ARMS=ARMS+REAL(RHOI*RHOI)
          ENDDO
        ENDDO
      ENDDO
      AMEAN=AMEAN/NDATA
      ARMS=ARMS/NDATA
      ARMS=ARMS-AMEAN*AMEAN
      IF(ARMS>0.0D0)ARMS=SQRT(ARMS)
      MAPN%AVG=AMEAN
      MAPN%STD=ARMS
      MAPN%MIN=AMIN
      MAPN%MAX=AMAX
      RETURN
      END SUBROUTINE MAPCAST


      SUBROUTINE TRANSCRD(NATRIG,IDXATM,CX,CY,CZ,T,U)
!_________________________________________________________________
!  transform atom coordinates to their fitting position
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
      use memory_module, only: memory_init,  &
            residue_pointer,coordinate
!
#  include "../include/memory.h"
      INTEGER NATRIG,IDXATM(*)
      _REAL_  CX,CY,CZ,T(3),U(3,3)
!
      INTEGER IRES,IATOM,I,IA
      _REAL_ XN,YN,ZN
      LOGICAL ATSKIP
!
      call memory_init()
      do ires = 1,nres
         do iatom = residue_pointer(ires), residue_pointer(ires+1)-1
           ATSKIP=.TRUE.
           DO I=1,NATRIG
             IA=IDXATM(I)
             IF(IATOM==IA)THEN
               ATSKIP=.FALSE.
               EXIT
             ENDIF
           ENDDO
           IF(ATSKIP)CYCLE
           XN=coordinate(1,ia)-T(1)-CX
           YN=coordinate(2,ia)-T(2)-CY
           ZN=coordinate(3,ia)-T(3)-CZ
           coordinate(1,ia)=U(1,1)*XN+U(2,1)*YN+U(3,1)*ZN+CX
           coordinate(2,ia)=U(1,2)*XN+U(2,2)*YN+U(3,2)*ZN+CY
           coordinate(3,ia)=U(1,3)*XN+U(2,3)*YN+U(3,3)*ZN+CZ
         end do
      end do
      RETURN
      END SUBROUTINE TRANSCRD


      SUBROUTINE WRTFITPDB(UNIT,NATRIG,IDXATM,CX,CY,CZ,T,U)
!_________________________________________________________________
!  transform atom coordinates to their fitting position
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
      use memory_module, only: memory_init,  &
            residue_pointer,residue_label,atom_name,coordinate
!
#  include "../include/memory.h"
      INTEGER UNIT,NATRIG,IDXATM(*)
      _REAL_  CX,CY,CZ,T(3),U(3,3)
      _REAL_ CRDN(3)
!
      INTEGER IRES,IATOM,I,IA
      _REAL_ XN,YN,ZN
      LOGICAL ATSKIP
      character(len=4) :: name
      character(len=4) :: resName
!
      call memory_init()
      write(unit,'(A)') 'REMARK  Written by Amber 14, SANDER, EMAP'
      do ires = 1,nres
         do iatom = residue_pointer(ires), residue_pointer(ires+1)-1
           ATSKIP=.TRUE.
           DO I=1,NATRIG
             IA=IDXATM(I)
             IF(IATOM==IA)THEN
               ATSKIP=.FALSE.
               EXIT
             ENDIF
           ENDDO
           IF(ATSKIP)CYCLE
           XN=coordinate(1,ia)-T(1)-CX
           YN=coordinate(2,ia)-T(2)-CY
           ZN=coordinate(3,ia)-T(3)-CZ
           CRDN(1)=U(1,1)*XN+U(2,1)*YN+U(3,1)*ZN+CX
           CRDN(2)=U(1,2)*XN+U(2,2)*YN+U(3,2)*ZN+CY
           CRDN(3)=U(1,3)*XN+U(2,3)*YN+U(3,3)*ZN+CZ
            name = atom_name(iatom)
            resName=residue_label(ires)
            resName=adjustr(resName)
            write(unit,'(A6,I5,1X,A4,A1,A4,A1,I4,A1,3X,3F8.3)')&
                  'ATOM  ', &
                  iatom,name,' ', &
                  resName,' ', &
                  ires,' ', &
                  CRDN(1:3)
         end do
      end do
      write(unit,'(A)') 'END'
      RETURN
      END SUBROUTINE WRTFITPDB


      SUBROUTINE WRTCCP4(FNAME,MAP,OUTU)
!_________________________________________________________________
!  write out CCP4 map information and its density array
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
      CHARACTER*80 FNAME
      INTEGER UNIT,OUTU
      type(EMAPOBJECT) :: map
!
      INTEGER*4 LX,LY,LZ,MODE,MNX,MNY,MNZ,NX,NY,NZ
      REAL*4 DX,DY,DZ,XL,YL,ZL,ALPHA,BETA,GAMMA
      REAL*4 X0,Y0,Z0
      INTEGER*4 MAPC,MAPR,MAPS
      REAL*4 AMAX,AMIN,AMEAN,ARMS
      INTEGER*4 ISPG,NSYMBT,LSKFLG,NNOTE
      REAL*4 SKWMAT(9),SKWTRN(3),EXTRA(15),EXTRA25(25)
      CHARACTER*4 MAPLABEL,MACHST
      CHARACTER*80 NOTES(10)
      INTEGER*1 BDATAI
      INTEGER*2 IDATAI
      REAL*4 RDATAI
      COMPLEX(int32) CDATAI
!
      INTEGER I,NDATA
      logical ismrc

!
      ismrc=(INDEX(FNAME,'.mrc')>0.or.INDEX(FNAME,'.MRC')>0)
!
      LX=MAP%LX
      LY=MAP%LY
      LZ=MAP%LZ
      NDATA=LX*LY*LZ
      MODE=MAP%MODE
      MNX=MAP%MX
      MNY=MAP%MY
      MNZ=MAP%MZ
      NX=MAP%NX
      NY=MAP%NY
      NZ=MAP%NZ
      DX=REAL(MAP%DX)
      DY=REAL(MAP%DY)
      DZ=REAL(MAP%DZ)
      XL=NX*DX
      YL=NY*DY
      ZL=NZ*DZ
      ALPHA=REAL(MAP%ALPHA)
      BETA=REAL(MAP%BETA)
      GAMMA=REAL(MAP%GAMMA)
      MAPC=1
      MAPR=2
      MAPS=3
      AMIN=REAL(MAP%MIN)
      AMAX=REAL(MAP%MAX)
      AMEAN=REAL(MAP%AVG)
      ARMS=REAL(MAP%STD)
      ISPG=0
      NSYMBT=0
      if(ismrc)then
        EXTRA25=0.0D0
        X0=MNX*DX
        Y0=MNY*DY
        Z0=MNZ*DZ
        EXTRA25=0.0D0
      else
        LSKFLG=0
        SKWMAT=0.0d0
        SKWMAT(1)=1.0D0
        SKWMAT(5)=1.0D0
        SKWMAT(9)=1.0D0
        SKWTRN=0.0d0
        EXTRA=0.0D0
      endif
      MAPLABEL='EMAP'
      MACHST='ALL '
      NNOTE=3
      NOTES=""
      NOTES(1)=" This map is created with the emap module "
      NOTES(2)=" Report questions to Dr. Xiongwu Wu  "
      NOTES(3)="             Email: wuxw@nhlbi.nih.gov "
      UNIT=99
      OPEN(UNIT,FILE=FNAME,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='UNKNOWN')
      WRITE(UNIT)LX,LY,LZ               !1,2,3
      WRITE(UNIT)MODE                   ! 4
      WRITE(UNIT)MNX,MNY,MNZ            ! 5,6,7
      WRITE(UNIT)NX,NY,NZ               ! 8,9,10
      WRITE(UNIT)XL,YL,ZL               ! 11,12,13
      WRITE(UNIT)ALPHA,BETA,GAMMA       ! 14,15,16
      WRITE(UNIT)MAPC,MAPR,MAPS         ! 17,18,19
      WRITE(UNIT)AMIN,AMAX,AMEAN        ! 20,21,22
      WRITE(UNIT)ISPG,NSYMBT            ! 23,24
      if(ismrc)then
        WRITE(UNIT)EXTRA25              ! 25-49
        WRITE(UNIT)X0,Y0,Z0             ! 50-52
      else
        WRITE(UNIT)LSKFLG               ! 25
        WRITE(UNIT)SKWMAT               ! 26-34
        WRITE(UNIT)SKWTRN               ! 35-37
        WRITE(UNIT)EXTRA                ! 38-52
      endif
      WRITE(UNIT)MAPLABEL               ! 53
      WRITE(UNIT)MACHST                 ! 54
      WRITE(UNIT)ARMS                   ! 55
      WRITE(UNIT)NNOTE                  ! 56
      WRITE(UNIT)NOTES                  ! 57-256
! write data
      IF(MODE==0)THEN
        DO I=1,NDATA
          BDATAI=MAP%BDATA(I)
          WRITE(UNIT)BDATAI
        ENDDO
      ELSE IF(MODE==1)THEN
        DO I=1,NDATA
          IDATAI=MAP%IDATA(I)
          WRITE(UNIT)IDATAI
        ENDDO
      ELSE IF(MODE==2)THEN
        DO I=1,NDATA
          RDATAI=MAP%RDATA(I)
          WRITE(UNIT)RDATAI
        ENDDO
      ELSE IF(MODE==3)THEN
        DO I=1,NDATA
          CDATAI=CMPLX(MAP%CDATA(I))
          WRITE(UNIT)CDATAI
        ENDDO
      ELSE IF(MODE==4)THEN
        DO I=1,NDATA
          CDATAI=CMPLX(MAP%CDATA(I))
          WRITE(UNIT)CDATAI
        ENDDO
      ELSE IF(MODE==5)THEN
        DO I=1,NDATA
          BDATAI=MAP%BDATA(I)
          WRITE(UNIT)BDATAI
        ENDDO
      ENDIF
      CLOSE(UNIT)
! print map information
      IF(OUTU > 0)THEN
      WRITE(OUTU,1000)TRIM(MAP%NAME), TRIM(FNAME)  !map name to be written
      WRITE(OUTU,1010)LX,LY,LZ               !1,2,3
      WRITE(OUTU,1020)MODE                   ! 4
      WRITE(OUTU,1030)MNX,MNY,MNZ            ! 5,6,7
      WRITE(OUTU,1040)NX,NY,NZ               ! 8,9,10
      WRITE(OUTU,1050)XL,YL,ZL               ! 11,12,13
      WRITE(OUTU,1060)ALPHA,BETA,GAMMA       ! 14,15,16
      WRITE(OUTU,1070)MAPC,MAPR,MAPS         ! 17,18,19
      WRITE(OUTU,1080)AMIN,AMAX,AMEAN,ARMS       ! 20,21,22,55
      WRITE(OUTU,1090)ISPG,NSYMBT     ! 23,24
      if(ismrc)then
        WRITE(OUTU,1102)EXTRA25                  ! 38-52
        WRITE(OUTU,1105)x0,y0,z0                  ! 38-52
      else
        WRITE(OUTU,1100)LSKFLG,NNOTE     ! 25,56
        WRITE(OUTU,1110)SKWMAT                 ! 26-34
        WRITE(OUTU,1120)SKWTRN                 ! 35-37
        WRITE(OUTU,1130)EXTRA                  ! 38-52
      endif
      WRITE(OUTU,1140)MAPLABEL               ! 53
      WRITE(OUTU,1150)MACHST                 ! 54
      WRITE(OUTU,1160)(I,NOTES(I),I=1,NNOTE)   ! 57-256
      WRITE(OUTU,1210)NDATA                  ! DATA NUMBER
      !WRITE(OUTU,1220)(MAP%RDATA(I),I=1,NDATA)           ! DATA 
      ENDIF
1000  FORMAT(" map object: ",A," is written to: ",A)
1010  FORMAT(" LX, LY, LZ              = ",3I8)
1020  FORMAT(" MODE                    = ",I8)
1030  FORMAT(" MX, MY, MZ              = ",3I8)
1040  FORMAT(" NX, NY, NZ              = ",3I8)
1050  FORMAT(" XL, YL, ZL              = ",3F8.2)
1060  FORMAT(" ALPHA,BETA,GAMMA        = ",3F8.2)
1070  FORMAT(" MAPC, MAPR, MAPS        = ",3I8)
1080  FORMAT(" MIN,MAX,MEAN,STD        = ",4E12.4)
1090  FORMAT(" ISPG,NSYMBT= ",2I8)
1102  FORMAT(" EXTRA                   = ",5F8.2)
1105  FORMAT(" X0,Y0,Z0                = ",3F8.2)
1100  FORMAT(" LSKFLG,NNOTE= ",2I8)
1110  FORMAT(" SKWMAT                  = ",3F8.2)
1120  FORMAT(" SKWTRN                  = ",3F8.2)
1130  FORMAT(" EXTRA                   = ",5F8.2)
1140  FORMAT(" MAPLABEL                = ",A)
1150  FORMAT(" MACHST                  = ",A)
1160  FORMAT(" NOTES ",I2,": ",A80)
1210  FORMAT(" DATA POINT NUMBER       = ",I8)
!1220  FORMAT(8F10.4)
      RETURN
      END SUBROUTINE WRTCCP4


      SUBROUTINE RDEMAP(FNAME,IDEMP,OUTU)
!-----------------------------------------------------------------------
!  read in map information and its density array
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!                        
      INTEGER IDEMP,OUTU
      CHARACTER*(*)FNAME
!
      INTEGER I
      _REAL_ RESO
!
      EMAPS(IDEMP)%name=FNAME
      IF(LEN_TRIM(FNAME)==0)THEN
        ! map will be created from input coordinates based on rigid mask
        RETURN      
      ELSE IF(INDEX(FNAME,'.ccp4')>0.or.INDEX(FNAME,'.CCP4')>0 .or.  &
         INDEX(FNAME,'.map')>0.or.INDEX(FNAME,'.MAP')>0 .or. &
         INDEX(FNAME,'.mrc')>0.or.INDEX(FNAME,'.MRC')>0)THEN
        CALL READCCP4(FNAME,IDEMP,OUTU)      
      ELSE IF(INDEX(FNAME,'.pdb')>0.or.INDEX(FNAME,'.PDB')>0)THEN
        DO I=1,NRIGID
          IF(EMRIGS(I)%MAPID==IDEMP)RESO=EMRIGS(I)%resolution
        ENDDO
        IF(RESO<1.0D-1)RESO=2.0D0
        CALL PDB2MAP(FNAME,IDEMP,RESO,OUTU)      
      ELSE
        write(outu,*)' Constraint map format is not supported:',fname
        stop
      ENDIF
      RETURN
      END SUBROUTINE RDEMAP

      SUBROUTINE READCCP4(FNAME,IDEMP,OUTU)
!_________________________________________________________________
!  read in CCP4 map information and its density array
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
     use constants
!
     CHARACTER*(*)FNAME
     INTEGER UNIT,IDEMP,OUTU
!
      INTEGER*4 LX,LY,LZ,MODE,MNX,MNY,MNZ,NX,NY,NZ
      REAL*4 XL,YL,ZL,ALPHA,BETA,GAMMA
      REAL*4 X0,Y0,Z0
      INTEGER*4 MAPC,MAPR,MAPS
      REAL*4 AMAX,AMIN,AMEAN,ARMS
      INTEGER*4 ISPG,NSYMBT,LSKFLG,NNOTE
      REAL*4 SKWMAT(9),SKWTRN(3),EXTRA(15),EXTRA25(25)
      CHARACTER*4 MAPLABEL,MACHST
      CHARACTER*80 NOTES(10)
      INTEGER*1 BDATAI
      INTEGER*2 IDATAI
      REAL*4 RDATAI
      COMPLEX(real32) CDATAI
!
      INTEGER*4 LXC,LYC,LZC
      INTEGER*4 IX,IY,IZ,JX,JY,JZ,JX1,JY1,JZ1
      REAL*8 DXU,DYU,DZU,DXC,DYC,DZC,DDX,DDY,DDZ,DDX1,DDY1,DDZ1
      REAL*8 COSA,COSB,COSG,SING,VUC
      REAL*8 MAT11,MAT12,MAT13,MAT22,MAT23,MAT33
      REAL*8 BG,XI,YI,ZI
      INTEGER I,IXYZ,NDATA,MDATA,alloc_err
      real(real32),pointer:: rdata(:)        !   real data
      logical ismrc

!
      ismrc=(INDEX(FNAME,'.mrc')>0.or.INDEX(FNAME,'.MRC')>0)
      X0=0.0d0
      y0=0.0d0
      z0=0.0d0
      IX=0
      IY=0
      IZ=0

      UNIT=99
      OPEN(UNIT,FILE=FNAME,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='OLD')
      !
      READ(UNIT)LX,LY,LZ               !1,2,3
      READ(UNIT)MODE                   ! 4
      READ(UNIT)MNX,MNY,MNZ            ! 5,6,7
      READ(UNIT)NX,NY,NZ               ! 8,9,10
      READ(UNIT)XL,YL,ZL               ! 11,12,13
      READ(UNIT)ALPHA,BETA,GAMMA       ! 14,15,16
      READ(UNIT)MAPC,MAPR,MAPS         ! 17,18,19
      READ(UNIT)AMIN,AMAX,AMEAN        ! 20,21,22
      READ(UNIT)ISPG,NSYMBT            ! 23,24
      if(ismrc)then
        READ(UNIT)EXTRA25              ! 25-49
        READ(UNIT)x0,y0,z0             ! 50-52
     else
        READ(UNIT)LSKFLG                 ! 25
        READ(UNIT)SKWMAT                 ! 26-34
        READ(UNIT)SKWTRN                 ! 35-37
        READ(UNIT)EXTRA                  ! 38-52
      endif
      READ(UNIT)MAPLABEL               ! 53
      READ(UNIT)MACHST                 ! 54
      READ(UNIT)ARMS                   ! 55
      READ(UNIT)NNOTE                  ! 56
      READ(UNIT)(NOTES(I),I=1,10)      ! 57-256
! write data
      NDATA=LX*LY*LZ
      allocate(RDATA(NDATA))
!     if(alloc_err /= 0 ) write(outu,*)"unable to allocate RDATA"
      IF(MODE==0)THEN
        DO I=1,NDATA
            READ(UNIT)BDATAI
            RDATA(I)=BDATAI
        ENDDO
      ELSE IF(MODE==1)THEN
        DO I=1,NDATA
          READ(UNIT)IDATAI
          RDATA(I)=REAL(IDATAI,real32)
        ENDDO
      ELSE IF(MODE==2)THEN
        DO I=1,NDATA
          READ(UNIT)RDATAI
          RDATA(I)=RDATAI
        ENDDO
      ELSE IF(MODE==3)THEN
        DO I=1,NDATA
          READ(UNIT)CDATAI
          RDATA(I)=REAL(CDATAI)
        ENDDO
      ELSE IF(MODE==4)THEN
        DO I=1,NDATA
          READ(UNIT)CDATAI
          RDATA(I)=REAL(CDATAI)
        ENDDO
      ELSE IF(MODE==5)THEN
        DO I=1,NDATA
          READ(UNIT)BDATAI
          RDATA(I)=BDATAI
        ENDDO
      ENDIF
      CLOSE(UNIT)
! convert to orthogonal system
!   v=sqrt(1-cos(alpha)^2-cos(beta)^2-cos(gamma)^2+2cos(alpha)cos(beta)cos(gamma))
! from cartisian to fraction
! {xu}={xm-ym*cos(gamma)/sin(gamma)+zm*(cos(alpha)cos(gamma)-cos(beta))/sin(gamma)/v)/a}
! {yu}={ym/sin(gamma)+zm*(cos(beta)cos(gamma)-cos(a))/sin(gamma)/v)/b}
! {zu}={zm*sin(gamma)/v/c}
! from fraction to cartisian
! {xc}={a*xu+b*cos(gamma)*yu+c*cos(beta)*zu}
! {yc}={b*sin(gamma)*yu+c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma)*zu}
! {zc}={zu*c*v/sin(gamma)}
!
!
      cosa=cos(alpha*DEG_TO_RAD)
      cosb=cos(beta*DEG_TO_RAD)
      cosg=cos(gamma*DEG_TO_RAD)
      sing=sin(gamma*DEG_TO_RAD)
      vuc=SQRT(one-cosa*cosa-cosb*cosb-cosg*cosg+two*cosa*cosb*cosg)
      mat11=one
      mat12=cosg
      mat13=cosb
      mat22=sing
      mat23=(cosa-cosb*cosg)/sing
      mat33=vuc/sing
      dxu=xl/nx
      dyu=yl/ny
      dzu=zl/nz
      dxc=dxu*mat11+dyu*mat12+dzu*mat13
      dyc=dyu*mat22+dzu*mat23
      dzc=dzu*mat33
      
      LXC=NINT((LX*dxu*mat11+LY*dyu*mat12+LZ*dzu*mat13)/dxc)
      LYC=NINT((LY*dyu*mat22+LZ*dzu*mat23)/dyc)
      LZC=NINT((LZ*dzu*mat33)/dzc)
      EMAPS(IDEMP)%LX=LXC
      EMAPS(IDEMP)%LY=LYC
      EMAPS(IDEMP)%LZ=LZC
      EMAPS(IDEMP)%MX=NINT((MNX*dxu*mat11+MNY*dyu*mat12+MNZ*dzu*mat13)/dxc)
      EMAPS(IDEMP)%MY=NINT((MNY*dyu*mat22+MNZ*dzu*mat23)/dyc)
      EMAPS(IDEMP)%MZ=NINT((MNZ*dzu*mat33)/dzc)
      EMAPS(IDEMP)%NX=NINT((NX*dxu*mat11+NY*dyu*mat12+NZ*dzu*mat13)/dxc)
      EMAPS(IDEMP)%NY=NINT((NY*dyu*mat22+NZ*dzu*mat23)/dyc)
      EMAPS(IDEMP)%NZ=NINT((NZ*dzu*mat33)/dzc)
      EMAPS(IDEMP)%MODE=2
      EMAPS(IDEMP)%DX=dxc
      EMAPS(IDEMP)%DY=dyc
      EMAPS(IDEMP)%DZ=dzc
      EMAPS(IDEMP)%ALPHA=90.0d0
      EMAPS(IDEMP)%BETA=90.0d0
      EMAPS(IDEMP)%GAMMA=90.0d0
!      MAPC=1
!      MAPR=2
!C      MAPS=3
!      ISPG=0
!      NSYMBT=0
!      LSKFLG=0
!      MAPLABEL='MAP '
!      MACHST='ALL '
!      NNOTE=3
!      NOTES(1)="          "
!     &    //" This map is created with the emap module of charmm "
!      NOTES(2)="          "
!     &    //" Report questions to Dr. Xiongwu Wu  "
!      NOTES(3)="          "
!     &    //"             Email: wuxw@nhlbi.nih.gov "
      EMAPS(IDEMP)%CX=(EMAPS(IDEMP)%MX+EMAPS(IDEMP)%LX/2.0D0)*EMAPS(IDEMP)%DX
      EMAPS(IDEMP)%CY=(EMAPS(IDEMP)%MY+EMAPS(IDEMP)%LY/2.0D0)*EMAPS(IDEMP)%DY
      EMAPS(IDEMP)%CZ=(EMAPS(IDEMP)%MZ+EMAPS(IDEMP)%LZ/2.0D0)*EMAPS(IDEMP)%DZ
      BG=SUM(RDATA)/NDATA
      MDATA=LXC*LYC*LZC
      if(associated(EMAPS(IDEMP)%RDATA))deallocate(EMAPS(IDEMP)%RDATA)
      allocate(EMAPS(IDEMP)%RDATA(MDATA),stat=alloc_err)
      if (alloc_err /= 0 ) then
        write(outu,'( "Error:  Unable to allocate EMAP rdata !" )')
        write(outu,'(" Work map size and background: ",3I6,F10.4)')LXC,LYC,LZC,BG
      endif
      EMAPS(IDEMP)%RDATA(:)=REAL(BG)
      DO I=1,NDATA
        IF(MAPC==1.and.MAPR==2)THEN
          IX=MOD(I-1,LX)
          IY=MOD(I-1,LX*LY)/LX
          IZ=(I-1)/(LX*LY)
        ELSE IF(MAPC==2.and.MAPR==3)THEN
          IZ=MOD(I-1,LZ)
          IX=MOD(I-1,LX*LZ)/LZ
          IY=(I-1)/(LX*LZ)
        ELSE IF(MAPC==1.and.MAPR==3)THEN
          IX=MOD(I-1,LX)
          IZ=MOD(I-1,LX*LZ)/LX
          IY=(I-1)/(LX*LZ)
        ELSE IF(MAPC==2.and.MAPR==1)THEN
          IY=MOD(I-1,LY)
          IX=MOD(I-1,LX*LY)/LY
          IZ=(I-1)/(LX*LY)
        ELSE IF(MAPC==3.and.MAPR==1)THEN
          IY=MOD(I-1,LY)
          IZ=MOD(I-1,LY*LZ)/LY
          IX=(I-1)/(LY*LZ)
        ELSE IF(MAPC==3.and.MAPR==2)THEN
          IZ=MOD(I-1,LZ)
          IY=MOD(I-1,LY*LZ)/LZ
          IX=(I-1)/(LY*LZ)
        ENDIF
        XI=IX*DXU*MAT11+IY*DYU*MAT12+IZ*DZU*MAT13
        YI=IY*DYU*MAT22+IZ*DZU*MAT23
        ZI=IZ*DZU*MAT33
        JX=MOD(FLOOR(XI/DXC),LXC)+1
        JY=MOD(FLOOR(YI/DYC),LYC)+1
        JZ=MOD(FLOOR(ZI/DZC),LZC)+1
        JX1=MOD(JX,LXC)+1
        JY1=MOD(JY,LYC)+1
        JZ1=MOD(JZ,LZC)+1
        DDX=FLOOR(XI/DXC)+ONE-XI/DXC
        DDY=FLOOR(YI/DYC)+ONE-YI/DYC
        DDZ=FLOOR(ZI/DZC)+ONE-ZI/DZC
        DDX1=ONE-DDX
        DDY1=ONE-DDY
        DDZ1=ONE-DDZ
        RDATAI=RDATA(I)-REAL(BG)
        IXYZ=JX+LXC*((JY-1+(JZ-1)*LYC))
        EMAPS(IDEMP)%RDATA(IXYZ)=EMAPS(IDEMP)%RDATA(IXYZ)+REAL(DDX*DDY*DDZ)*RDATAI
        IXYZ=JX+LXC*((JY1-1+(JZ-1)*LYC))
        EMAPS(IDEMP)%RDATA(IXYZ)=EMAPS(IDEMP)%RDATA(IXYZ)+REAL(DDX*DDY1*DDZ)*RDATAI
        IXYZ=JX+LXC*((JY-1+(JZ1-1)*LYC))
        EMAPS(IDEMP)%RDATA(IXYZ)=EMAPS(IDEMP)%RDATA(IXYZ)+REAL(DDX*DDY*DDZ1)*RDATAI
        IXYZ=JX+LXC*((JY1-1+(JZ1-1)*LYC))
        EMAPS(IDEMP)%RDATA(IXYZ)=EMAPS(IDEMP)%RDATA(IXYZ)+REAL(DDX*DDY1*DDZ1)*RDATAI
        IXYZ=JX1+LXC*((JY-1+(JZ-1)*LYC))
        EMAPS(IDEMP)%RDATA(IXYZ)=EMAPS(IDEMP)%RDATA(IXYZ)+REAL(DDX1*DDY*DDZ)*RDATAI
        IXYZ=JX1+LXC*((JY1-1+(JZ-1)*LYC))
        EMAPS(IDEMP)%RDATA(IXYZ)=EMAPS(IDEMP)%RDATA(IXYZ)+REAL(DDX1*DDY1*DDZ)*RDATAI
        IXYZ=JX1+LXC*((JY-1+(JZ1-1)*LYC))
        EMAPS(IDEMP)%RDATA(IXYZ)=EMAPS(IDEMP)%RDATA(IXYZ)+REAL(DDX1*DDY*DDZ1)*RDATAI
        IXYZ=JX1+LXC*((JY1-1+(JZ1-1)*LYC))
        EMAPS(IDEMP)%RDATA(IXYZ)=EMAPS(IDEMP)%RDATA(IXYZ)+REAL(DDX1*DDY1*DDZ1)*RDATAI
        !if(outu>0)write(outu,*)'jc=',ix,iy,iz,jx,jy,jz,xi,yi,zi,rdatai
      ENDDO
!     statistics
      CALL MAPSTAT(IDEMP)
! print map information
      IF(OUTU > 0 )THEN
      WRITE(OUTU,1000)IDEMP,TRIM(EMAPS(IDEMP)%NAME), UNIT  !map name to be written
      WRITE(OUTU,1010)LX,LY,LZ               !1,2,3
      WRITE(OUTU,1020)MODE                   ! 4
      WRITE(OUTU,1030)MNX,MNY,MNZ            ! 5,6,7
      WRITE(OUTU,1040)NX,NY,NZ               ! 8,9,10
      WRITE(OUTU,1050)XL,YL,ZL               ! 11,12,13
      WRITE(OUTU,1060)ALPHA,BETA,GAMMA       ! 14,15,16
      WRITE(OUTU,1070)MAPC,MAPR,MAPS         ! 17,18,19
      WRITE(OUTU,1080)AMIN,AMAX,AMEAN,ARMS       ! 20,21,22,55
      WRITE(OUTU,1090)ISPG,NSYMBT     ! 23,24
      if(ismrc)then
        WRITE(OUTU,1102)EXTRA25                  ! 38-52
        WRITE(OUTU,1105)x0,y0,z0                  ! 38-52
      else
        WRITE(OUTU,1100)LSKFLG,NNOTE     ! 25,56
        WRITE(OUTU,1110)SKWMAT                 ! 26-34
        WRITE(OUTU,1120)SKWTRN                 ! 35-37
        WRITE(OUTU,1130)EXTRA                  ! 38-52
      endif
      CALL ASCIIMASK(MAPLABEL)
      WRITE(OUTU,1140)MAPLABEL               ! 53
      CALL ASCIIMASK(MACHST)
      WRITE(OUTU,1150)MACHST                 ! 54
      DO I=1,NNOTE
        CALL ASCIIMASK(NOTES(I))
        IF(LEN(trim(NOTES(I)))>0)WRITE(OUTU,1160)I,NOTES(I)   ! 57-256
      ENDDO
      WRITE(OUTU,1210)NDATA                  ! DATA NUMBER
      !WRITE(OUTU,1220)(RDATA(I),EMAPS(IDEMP)%RDATA(I),I=1,NDATA)           ! DATA 
      WRITE(OUTU,1230)IDEMP                    ! DATA NUMBER
      ENDIF
      if(associated(RDATA))deallocate(RDATA)
1000  FORMAT(" ------------------EMAP IMAGE ",I4," INPUT ---------------" / &
            "map file: ",A," is read from unit: ",I4)
1010  FORMAT(" LX, LY, LZ              = ",3I8)
1020  FORMAT(" MODE                    = ",I8)
1030  FORMAT(" MX, MY, MZ              = ",3I8)
1040  FORMAT(" NX, NY, NZ              = ",3I8)
1050  FORMAT(" XL, YL, ZL              = ",3F8.2)
1060  FORMAT(" ALPHA,BETA,GAMMA        = ",3F8.2)
1070  FORMAT(" MAPC, MAPR, MAPS        = ",3I8)
1080  FORMAT(" MIN,MAX,MEAN,STD        = ",4E12.4)
1090  FORMAT(" ISPG,NSYMBT= ",2I8)
1102  FORMAT(" EXTRA                   = ",5F8.2)
1105  FORMAT(" X0,Y0,Z0                = ",3F8.2)
1100  FORMAT(" LSKFLG,NNOTE= ",2I8)
1110  FORMAT(" SKWMAT                  = ",3F8.2)
1120  FORMAT(" SKWTRN                  = ",3F8.2)
1130  FORMAT(" EXTRA                   = ",5F8.2)
1140  FORMAT(" MAPLABEL                = ",A)
1150  FORMAT(" MACHST                  = ",A)
1160  FORMAT(" NOTES ",I2,": ",A80)
1210  FORMAT(" DATA POINT NUMBER       = ",I8)
!1220  FORMAT(8f10.4)
1230  FORMAT(" ----------------------- END OF EMAP IMAGE ",I4, &
             "  -------------------------- ")
      RETURN
      END SUBROUTINE READCCP4
     
      SUBROUTINE ASCIIMASK(STRING)
      CHARACTER*(*) STRING
      INTEGER I,J
      I=1
      DO WHILE(I.LT.LEN(STRING))
        J=ICHAR(STRING(I:I))
        IF(J<=0.or.J>127)THEN
          STRING(I:)=STRING(I+1:)
        ELSE
          I=I+1
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE ASCIIMASK

      SUBROUTINE PDB2MAP(FNAME,IDEMP,RESO,OUTU)
!_________________________________________________________________
!  read in PDB coordinates to create a map object
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
      use ew_bspline,only:fill_bspline_0
!
      CHARACTER*(*)FNAME
      INTEGER UNIT,IDEMP,OUTU
      _REAL_ RESO
!
      CHARACTER*80 pdblin
      CHARACTER*4 CHAINI,RESIN,ATOMIN,RID,ATOMNM
      _REAL_ XIN,YIN,ZIN,AMASSI
      INTEGER NAT,IAT,ISEQ
!
      INTEGER alloc_err
      _REAL_, allocatable, dimension(:)::tmp_crd,tmp_amass
!
      UNIT=99
      OPEN(UNIT,FILE=FNAME,ACCESS='SEQUENTIAL',STATUS='OLD')
      NAT=0
99    READ(UNIT,'(A)',END=120,ERR=120) PDBLIN
      PDBLIN=TRIM(PDBLIN)
      IF(LEN_TRIM(PDBLIN)<44)GOTO 99
      IF (PDBLIN(1:6).NE.'ATOM  '.AND.PDBLIN(1:6).NE.'HETATM') GOTO 99
      NAT=NAT+1
      GOTO 99
120   CONTINUE
      if(allocated(tmp_crd))deallocate(tmp_crd)
      allocate(tmp_crd(3*NAT),stat=alloc_err)
      if(alloc_err /= 0 ) write(6,*)"unable to allocate CRD !"
      if(allocated(tmp_amass))deallocate(tmp_amass)
      allocate(tmp_amass(NAT),stat=alloc_err)
      if(alloc_err /= 0 ) write(6,*)"unable to allocate AMASS !"
      CLOSE(UNIT)
      OPEN(UNIT,FILE=FNAME,ACCESS='SEQUENTIAL',STATUS='OLD')
      IAT=0
299   READ(UNIT,'(A)',END=320,ERR=320) PDBLIN
      PDBLIN=TRIM(PDBLIN)
      IF(LEN(PDBLIN)<44)GOTO 299
      IF (PDBLIN(1:6).NE.'ATOM  '.AND.PDBLIN(1:6).NE.'HETATM') GOTO 299
      READ(PDBLIN, &
          '(6X,I5,1X,A4,1X,A3,1X,A1,A5,3X,3F8.3)') &
               ISEQ,ATOMIN,RESIN,CHAINI,RID,XIN,YIN,ZIN     
!  determine atom mass
      ATOMNM=ADJUSTL(ATOMIN)
      IF(INDEX(ATOMNM,'CAL')==1)THEN
        AMASSI=40.0D0
      ELSE IF(INDEX(ATOMNM,'CL')==1)THEN
        AMASSI=35.5D0
      ELSE IF(INDEX(ATOMNM,'MG')==1)THEN
        AMASSI=24.3D0
      ELSE IF(INDEX(ATOMNM,'S')==1)THEN
        AMASSI=32.06D0
      ELSE IF(INDEX(ATOMNM,'O')==1)THEN
        AMASSI=16.0D0
      ELSE IF(INDEX(ATOMNM,'N')==1)THEN
        AMASSI=14.01D0
      ELSE IF(INDEX(ATOMNM,'P')==1)THEN
        AMASSI=31.0D0
      ELSE IF(INDEX(ATOMNM,'F')==1)THEN
        AMASSI=19.0D0
      ELSE IF(INDEX(ATOMNM,'C')==1)THEN
        AMASSI=12.01D0
      ELSE IF(INDEX(ATOMNM,'H')==1)THEN
        AMASSI=1.008D0
      ELSE 
        AMASSI=1.008D0
      ENDIF
      IAT=IAT+1
      tmp_amass(IAT)=AMASSI
      tmp_crd(3*IAT-2)=XIN
      tmp_crd(3*IAT-1)=YIN
      tmp_crd(3*IAT)=ZIN
      goto 299
320   CONTINUE
!     assign map object
      CALL COR2MAP(IDEMP,RESO,NAT,tmp_amass,tmp_crd,OUTU)
      if(allocated(tmp_crd))deallocate(tmp_crd)
      if(allocated(tmp_amass))deallocate(tmp_amass)
      CLOSE(UNIT)
      RETURN
      END SUBROUTINE PDB2MAP


      SUBROUTINE COR2MAP(IDEMP,RESO,NATC,AMASS,CRD,OUTU)
!_________________________________________________________________
!  Blur input coordinates to a map object
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
      use ew_bspline,only:fill_bspline_0
!
      INTEGER IDEMP,NATC,OUTU
      _REAL_ AMASS(*)
      _REAL_ CRD(*)
      _REAL_ RESO
!
      real(real32),pointer:: rho(:)=> null()        !   real data
!
      INTEGER*4 LX,LY,LZ,MODE,MNX,MNY,MNZ
      REAL*4 XC,YC,ZC
      CHARACTER*80 NOTES(10)
      REAL*4 XIN,YIN,ZIN,AMASSI,DG
      INTEGER IAT
!
      INTEGER I,J,K,M1,M2,M3,NDATA
      INTEGER IPT1,IPT2,IPT3,ITH1,ITH2,ITH3
      _REAL_ XMIN,YMIN,ZMIN,XMAX,YMAX,ZMAX
      _REAL_ FR1,FR2,FR3,W
      REAL*8 BS1(8),BS2(8),BS3(8)
      REAL*8 VAL0,VAL0A
      INTEGER alloc_err

!
      MODE=2
      XC=0.0D0
      YC=0.0D0
      ZC=0.0D0
      XMIN=1.0D9
      YMIN=1.0D9
      ZMIN=1.0D9
      XMAX=-1.0D9
      YMAX=-1.0D9
      ZMAX=-1.0D9
      DO IAT=1,NATC
        XIN=REAL(CRD(3*IAT-2))
        YIN=REAL(CRD(3*IAT-1))
        ZIN=REAL(CRD(3*IAT))
        XC=XC+XIN
        YC=YC+YIN
        ZC=ZC+ZIN
        IF(XMIN>XIN)XMIN=XIN
        IF(YMIN>YIN)YMIN=YIN
        IF(ZMIN>ZIN)ZMIN=ZIN
        IF(XMAX<XIN)XMAX=XIN
        IF(YMAX<YIN)YMAX=YIN
        IF(ZMAX<ZIN)ZMAX=ZIN
      ENDDO
      XC=XC/NATC
      YC=YC/NATC
      ZC=ZC/NATC
      DG=REAL(2.0D0*RESO/DBLE(BORDER))
      MNX=nint((XMIN-RESO)/DG-0.5d0)
      MNY=nint((YMIN-RESO)/DG-0.5d0)
      MNZ=nint((ZMIN-RESO)/DG-0.5d0)
      LX=nint((XMAX-XMIN+2.0D0*RESO)/DG+0.5d0)+1
      LY=nint((YMAX-YMIN+2.0D0*RESO)/DG+0.5d0)+1
      LZ=nint((ZMAX-ZMIN+2.0D0*RESO)/DG+0.5d0)+1
      NDATA=LX*LY*LZ
      if(associated(EMAPS(IDEMP)%RDATA))deallocate(EMAPS(IDEMP)%RDATA)
      allocate(EMAPS(IDEMP)%RDATA(NDATA),stat=alloc_err)
      if(alloc_err /= 0 ) write(6,*)"unable to allocate RDATA"
      RHO=>EMAPS(IDEMP)%RDATA
      RHO=0.0D0
!   generate map with b-spline
      DO IAT=1,NATC
!  determine atom mass
        AMASSI=REAL(AMASS(IAT))
        XIN=REAL(CRD(3*IAT-2))
        YIN=REAL(CRD(3*IAT-1))
        ZIN=REAL(CRD(3*IAT))
! calculate B-spline parameters
        FR1=XIN/DG
        FR2=YIN/DG
        FR3=ZIN/DG
        M1=INT(ANINT(FR1-0.5D0))
        W=FR1-M1
        call fill_bspline_0(w,BORDER,BS1)
        M2=INT(ANINT(FR2-0.5D0))
        W=FR2-M2
        call fill_bspline_0(w,BORDER,BS2)
        M3=INT(ANINT(FR3-0.5D0))
        W=FR3-M3
        call fill_bspline_0(w,BORDER,BS3)
! Calculate B-spline interception
         K = M3-MNZ+1 - BORDER/2
         DO 410 ITH3 = 1,BORDER
           K=K+1
           IF(K<1)THEN
             IPT1=0
             ELSE IF(K>LZ)THEN
               IPT1=LX*LY*(LZ-1)
             ELSE
               IPT1=(K-1)*LY*LX
             ENDIF
             VAL0A =  AMASSI*BS3(ITH3)
!
             J = M2-MNY+1 - BORDER/2
!
             DO 420 ITH2 = 1,BORDER
               J=J+1
               IF(J<1)THEN
                 IPT2=IPT1
               ELSE IF(J>LY)THEN
                 IPT2=IPT1+LX*(LY-1)
               ELSE
                 IPT2=IPT1+(J-1)*LX
               ENDIF
!
               VAL0= VAL0A * BS2(ITH2)
!
               I = M1-MNX+1 - BORDER/2
!
               DO 430 ITH1 = 1,BORDER
                 I=I+1
                 IF(I<1)THEN
                   IPT3=IPT2+1
                 ELSE IF(I>LX)THEN
                   IPT3=IPT2+LX
                 ELSE
                   IPT3=IPT2+I
                 ENDIF
                 RHO(IPT3)=RHO(IPT3)+ REAL(VAL0 * BS1(ITH1))
!
430          CONTINUE
420        CONTINUE
410      CONTINUE
!
      ENDDO
!     assign map object
      
      EMAPS(IDEMP)%LX=LX
      EMAPS(IDEMP)%LY=LY
      EMAPS(IDEMP)%LZ=LZ
      EMAPS(IDEMP)%MODE=MODE
      EMAPS(IDEMP)%MX=MNX
      EMAPS(IDEMP)%MY=MNY
      EMAPS(IDEMP)%MZ=MNZ
      EMAPS(IDEMP)%NX=LX
      EMAPS(IDEMP)%NY=LY
      EMAPS(IDEMP)%NZ=LZ
      EMAPS(IDEMP)%DX=DG
      EMAPS(IDEMP)%DY=DG
      EMAPS(IDEMP)%DZ=DG
      EMAPS(IDEMP)%ALPHA=90.0D0
      EMAPS(IDEMP)%BETA=90.0D0
      EMAPS(IDEMP)%GAMMA=90.0D0
!      MAPC=1
!      MAPR=2
!C      MAPS=3
!      ISPG=0
!      NSYMBT=0
!      LSKFLG=0
!      MAPLABEL='MAP '
!      MACHST='ALL '
!      NNOTE=3
      NOTES(1)="          " &
         //" This map is created with the emap module of charmm "
      NOTES(2)="          " &
         //" Report questions to Dr. Xiongwu Wu  "
      NOTES(3)="          " &
        //"             Email: wuxw@nhlbi.nih.gov "
      EMAPS(IDEMP)%CX=XC
      EMAPS(IDEMP)%CY=YC
      EMAPS(IDEMP)%CZ=ZC
!     statistics
      CALL MAPSTAT(IDEMP)
! print map information
      IF(OUTU > 0 )THEN
      WRITE(OUTU,1000)IDEMP  !map id
      WRITE(OUTU,1010)LX,LY,LZ               !1,2,3
      WRITE(OUTU,1020)MODE                   ! 4
      WRITE(OUTU,1030)MNX,MNY,MNZ            ! 5,6,7
      WRITE(OUTU,1050)LX*DG,LY*DG,LZ*DG            ! 5,6,7
      WRITE(OUTU,1080)EMAPS(IDEMP)%MIN,EMAPS(IDEMP)%MAX,EMAPS(IDEMP)%AVG,EMAPS(IDEMP)%STD       ! 20,21,22,55
!      WRITE(OUTU,1090)ISPG,NSYMBT,LSKFLG,NNOTE     ! 23,24,25,56
!      WRITE(OUTU,1110)SKWMAT                 ! 26-34
!      WRITE(OUTU,1120)SKWTRN                 ! 35-37
!      WRITE(OUTU,1130)EXTRA                  ! 38-52
!      WRITE(OUTU,1140)MAPLABEL               ! 53
!      WRITE(OUTU,1150)MACHST                 ! 54
!      WRITE(OUTU,1160)(I,NOTES(I),I=1,NNOTE)   ! 57-256
      WRITE(OUTU,1210)NDATA                  ! DATA NUMBER
      WRITE(OUTU,1220)IDEMP                    ! DATA NUMBER
      ENDIF
1000  FORMAT(" ------------------EMAP ID ",I4," CREATED ---------------" )
1010  FORMAT(" LX, LY, LZ              = ",3I8)
1020  FORMAT(" MODE                    = ",I4)
1030  FORMAT(" MX, MY, MZ              = ",3I8)
!1040  FORMAT(" NX, NY, NZ              = ",3I8)
1050  FORMAT(" XL, YL, ZL              = ",3F8.2)
!1060  FORMAT(" ALPHA,BETA,GAMMA        = ",3F8.2)
!1070  FORMAT(" MAPC, MAPR, MAPS        = ",3I8)
1080  FORMAT(" MIN,MAX,MEAN,STD        = ",4E12.4)
!1090  FORMAT(" ISPG,NSYMBT,LSKFLG,NNOTE= ",4I8)
!1110  FORMAT(" SKWMAT                  = ",3F8.2)
!1120  FORMAT(" SKWTRN                  = ",3F8.2)
!1130  FORMAT(" EXTRA                   = ",5F8.2)
!1140  FORMAT(" MAPLABEL                = ",A)
!1150  FORMAT(" MACHST                  = ",A)
!1160  FORMAT(" NOTES ",I2,":"/A)
1210  FORMAT(" DATA POINT NUMBER       = ",I8)
1220  FORMAT(" ----------------------- END OF EMAP IMAGE ",I4, &
             "  -------------------------- ")
      RETURN
      END SUBROUTINE COR2MAP


      SUBROUTINE MAPSTAT(IDEMP)
!_________________________________________________________________
!  calculate statistics of a map
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
!
      INTEGER IDEMP
      real(real32),pointer:: rho(:)=> null()        !   real data
!
      INTEGER NDATA,I
      real(real32) RHOI
      REAL*8 AMAX,AMIN,AMEAN,ARMS
!
      NDATA=EMAPS(IDEMP)%LX*EMAPS(IDEMP)%LY*EMAPS(IDEMP)%LZ
      RHO=>EMAPS(IDEMP)%RDATA
!     map statistics
      AMIN=1.0D9
      AMAX=-1.0D9
      AMEAN=0.0D0
      ARMS=0.0D0
      DO I=1,NDATA
        RHOI=RHO(I)
        IF(AMIN>RHOI)AMIN=RHOI
        IF(AMAX<RHOI)AMAX=RHOI
        AMEAN=AMEAN+RHOI
        ARMS=ARMS+RHOI*RHOI
      ENDDO
      AMEAN=AMEAN/NDATA
      ARMS=ARMS/NDATA-AMEAN*AMEAN
      IF(ARMS>0.0D0)ARMS=SQRT(ARMS)
      EMAPS(IDEMP)%MIN=AMIN
      EMAPS(IDEMP)%MAX=AMAX
      EMAPS(IDEMP)%AVG=AMEAN
      EMAPS(IDEMP)%STD=ARMS
      RETURN
      END SUBROUTINE MAPSTAT

!+ Read maps for constraint
subroutine emapforce(natom,enemap,amass,x,f)
!________________________________________________________________
!  Calculate map constraint forces
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
      _REAL_ ENEMAP,AMASS(*),X(*),F(*)
      INTEGER NATOM,I
      _REAL_ ENRIG,SCRIG,TMASS
      ENEMAP=0.0D0
      SCEMAP=0.0D0
      TMASS=0.0D0
!  Loop over all active rigid domains
      DO I=1,NRIGID
        CALL RIGIDFORCE(NATOM,ENRIG,SCRIG,AMASS,X,F,EMRIGS(I))
        ENEMAP=ENEMAP+ENRIG
        SCEMAP=SCEMAP+SCRIG
        TMASS=TMASS+EMRIGS(I)%MASS
      ENDDO
      SCEMAP=SCEMAP/TMASS
      return
      end subroutine emapforce

      SUBROUTINE RIGIDENG(NATOM,ENRIG,SCRIG,AMASS,X,rigobj)
!_________________________________________________________________
!   This routine transform reference atoms to map frame and
!     calculate B-splie parameters and the energy and force
!     then transform the force back to simulation frame
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
!
      use ew_bspline,only:fill_bspline_1
#ifdef MPI
      use mpi
#endif
      implicit none
!
#ifdef MPI
#  include "parallel.h"
#endif
      INTEGER NATOM
      REAL*8 ENRIG,SCRIG,AMASS(*),X(*)
      type(EMAPRIGID) :: rigobj
      type(EMAPOBJECT) :: mapobj
      real(real32),pointer:: rho(:)=> null()        !   real data
      integer,pointer:: idx(:)=> null()        !   real data
      INTEGER NATC,MNX,MNY,MNZ,LX,LY,LZ
      REAL*8 TR(3),U(3,3)
      REAL*8 RHO0,CX,CY,CZ,GX,GY,GZ
      REAL*8 BS1(8),BS2(8),BS3(8)
      REAL*8 DBS1(8),DBS2(8),DBS3(8)
      REAL*8 XI,YI,ZI,XJ,YJ,ZJ
      REAL*8 AMASSI,RHOI,ENRIGM
      REAL*8 CFACT,CFACTX,CFACTY,CFACTZ,FR1,FR2,FR3,W
      REAL*8 VAL0,VAL0A
      INTEGER ID,IAT,IATX,IATY,IATZ,I,J,K,M1,M2,M3,NSTA,NEND
      INTEGER ITH1,ITH2,ITH3,IPT1,IPT2,IPT3
#ifdef MPI
      INTEGER IERR
      _REAL_ TEMP1(10),TEMP2(10)

      nsta = iparpt(mytaskid) + 1
      nend = iparpt(mytaskid+1)
#else
      nsta = 1
      nend = natom
#endif
      mapobj=emaps(rigobj%mapid)
      MNX=MAPOBJ%MX
      MNY=MAPOBJ%MY
      MNZ=MAPOBJ%MZ
      LX=MAPOBJ%LX
      LY=MAPOBJ%LY
      LZ=MAPOBJ%LZ
      GX=MAPOBJ%DX
      GY=MAPOBJ%DY
      GZ=MAPOBJ%DZ
      CX=MAPOBJ%CX
      CY=MAPOBJ%CY
      CZ=MAPOBJ%CZ
      RHO=>MAPOBJ%RDATA
      RHO0=mapobj%AVG
      CFACT=rigobj%fcons/mapobj%std
      CFACTX=CFACT/GX
      CFACTY=CFACT/GY
      CFACTZ=CFACT/GZ
      NATC=RIGOBJ%NATC
      U=RIGOBJ%ROT
      TR=RIGOBJ%TRAN
      IDX=>RIGOBJ%IDXATOM
      SCRIG=0.0D0
!  Calculate selected atoms
      DO ID=1,NATC
        IAT=IDX(ID)
        IF(IAT<NSTA.OR.IAT>NEND)CYCLE
        IATX=3*IAT-2
        IATY=IATX+1
        IATZ=IATY+1
        AMASSI=AMASS(IAT)
        XI=X(IATX)-CX-TR(1)
        YI=X(IATY)-CY-TR(2)
        ZI=X(IATZ)-CZ-TR(3)
! transform to map frame
        XJ=XI*U(1,1)+YI*U(2,1)+ZI*U(3,1)+CX
        YJ=XI*U(1,2)+YI*U(2,2)+ZI*U(3,2)+CY
        ZJ=XI*U(1,3)+YI*U(2,3)+ZI*U(3,3)+CZ
! calculate B-spline parameters
        FR1=XJ/GX
        FR2=YJ/GY
        FR3=ZJ/GZ
        M1=INT(ANINT(FR1-0.5D0))
        W=FR1-M1
        call fill_bspline_1(w,BORDER,BS1,dBS1)
        M2=INT(ANINT(FR2-0.5D0))
        W=FR2-M2
        call fill_bspline_1(w,BORDER,BS2,dBS2)
        M3=INT(ANINT(FR3-0.5D0))
        W=FR3-M3
        call fill_bspline_1(w,BORDER,BS3,dBS3)
! Calculate B-spline interception
        K = M3-MNZ+1 - BORDER/2
        DO 100 ITH3 = 1,BORDER
           K=K+1
           IF(K<=0.OR.K>LZ)GOTO 100
           IPT1=(K-1)*LY*LX
            VAL0A = AMASSI  * BS3(ITH3)
!
            J = M2-MNY+1 - BORDER/2
!
            DO 200 ITH2 = 1,BORDER
               J=J+1
               IF(J<=0.OR.J>LY)GOTO 200
               IPT2=IPT1+(J-1)*LX
!
               VAL0= VAL0A * BS2(ITH2)
!
               I = M1-MNX+1 - BORDER/2
!
               DO 300 ITH1 = 1,BORDER
                  I=I+1
                  IF(I<=0.OR.I>LX)GOTO 300
                  IPT3=IPT2+I
                  RHOI=RHO(IPT3)-RHO0
                  SCRIG = SCRIG + VAL0 * RHOI * BS1(ITH1)
!
300          CONTINUE
200        CONTINUE
100      CONTINUE
!
      ENDDO
      ENRIG=-CFACT*SCRIG
#ifdef MPI
        IF(SANDERSIZE > 1)THEN
!  Combining all node results
!
          TEMP1(1)=SCRIG
          call mpi_allreduce(MPI_IN_PLACE,temp1,1,MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ierr)
          ENRIGM=TEMP1(1)
        ENDIF
#endif
      ENRIGM=-CFACT*SCRIG
      rigobj%energy=enrigm
      RETURN
      END SUBROUTINE RIGIDENG

      SUBROUTINE RIGIDFORCE(NATOM,ENRIG,SCRIG,AMASS,X,F,rigobj)
!_________________________________________________________________
!   This routine transform reference atoms to map frame and
!     calculate B-splie parameters and the energy and force
!     then transform the force back to simulation frame
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
!
      use ew_bspline,only:fill_bspline_1
#ifdef MPI
      use mpi
#endif
      implicit none
!
#ifdef MPI
#  include "parallel.h"
#endif
      INTEGER NATOM
      REAL*8 ENRIG,SCRIG,AMASS(*),X(*),F(*)
      type(EMAPRIGID) :: rigobj
      type(EMAPOBJECT) :: mapobj
      real(real32),pointer:: rho(:)=> null()        !   real data
      integer,pointer:: idx(:)=> null()        !   real data
      INTEGER NATC,MNX,MNY,MNZ,LX,LY,LZ
      REAL*8 TR(3),U(3,3)
      REAL*8 RHO0,CX,CY,CZ,GX,GY,GZ
      REAL*8 BS1(8),BS2(8),BS3(8)
      REAL*8 DBS1(8),DBS2(8),DBS3(8)
      REAL*8 XI,YI,ZI,XJ,YJ,ZJ,XC,YC,ZC
      REAL*8 AMASSI,RHOI,ENRIGM,F1,F2,F3,FX,FY,FZ
      REAL*8 CFACT,CFACTX,CFACTY,CFACTZ,FR1,FR2,FR3,T2,W
      REAL*8 rfx,rfy,rfz,rtx,rty,rtz
      REAL*8 VAL0,VAL1,VAL2,VAL3,VAL0A,VAL1A,VAL2A,VAL3A
      INTEGER ID,IAT,IATX,IATY,IATZ,I,J,K,M1,M2,M3,NSTA,NEND
      INTEGER ITH1,ITH2,ITH3,IPT1,IPT2,IPT3
#ifdef MPI
      INTEGER IERR
      _REAL_ TEMP1(10),TEMP2(10)
#endif
!
#ifdef MPI
      nsta = iparpt(mytaskid) + 1
      nend = iparpt(mytaskid+1)
#else
      nsta = 1
      nend = natom
#endif
      mapobj=emaps(rigobj%mapid)
      MNX=MAPOBJ%MX
      MNY=MAPOBJ%MY
      MNZ=MAPOBJ%MZ
      LX=MAPOBJ%LX
      LY=MAPOBJ%LY
      LZ=MAPOBJ%LZ
      GX=MAPOBJ%DX
      GY=MAPOBJ%DY
      GZ=MAPOBJ%DZ
      CX=MAPOBJ%CX
      CY=MAPOBJ%CY
      CZ=MAPOBJ%CZ
      RHO=>MAPOBJ%RDATA
      RHO0=mapobj%AVG
      CFACT=rigobj%fcons/mapobj%std
      CFACTX=CFACT/GX
      CFACTY=CFACT/GY
      CFACTZ=CFACT/GZ
      NATC=RIGOBJ%NATC
      U=RIGOBJ%ROT
      TR=RIGOBJ%TRAN
      IDX=>RIGOBJ%IDXATOM
      SCRIG=0.0D0
      rfx=0.0D0
      rfy=0.0D0
      rfz=0.0D0
      rtx=0.0D0
      rty=0.0D0
      rtz=0.0D0
      xc=0.0d0
      yc=0.0d0
      zc=0.0d0
!  Calculate selected atoms
      DO ID=1,NATC
        IAT=IDX(ID)
        IATX=3*IAT-2
        IATY=IATX+1
        IATZ=IATY+1
        XI=X(IATX)
        YI=X(IATY)
        ZI=X(IATZ)
        xc=xc+xi
        yc=yc+yi
        zc=zc+zi
        IF(IAT<NSTA.OR.IAT>NEND)CYCLE
        AMASSI=AMASS(IAT)
        XI=XI-CX-TR(1)
        YI=YI-CY-TR(2)
        ZI=ZI-CZ-TR(3)
! transform to map frame
        XJ=XI*U(1,1)+YI*U(2,1)+ZI*U(3,1)+CX
        YJ=XI*U(1,2)+YI*U(2,2)+ZI*U(3,2)+CY
        ZJ=XI*U(1,3)+YI*U(2,3)+ZI*U(3,3)+CZ
! calculate B-spline parameters
        FR1=XJ/GX
        FR2=YJ/GY
        FR3=ZJ/GZ
        M1=INT(ANINT(FR1-0.5D0))
        W=FR1-M1
        call fill_bspline_1(w,BORDER,BS1,dBS1)
        M2=INT(ANINT(FR2-0.5D0))
        W=FR2-M2
        call fill_bspline_1(w,BORDER,BS2,dBS2)
        M3=INT(ANINT(FR3-0.5D0))
        W=FR3-M3
        call fill_bspline_1(w,BORDER,BS3,dBS3)
! Calculate B-spline interception
        F1 = 0.0D0 
        F2 = 0.0D0 
        F3 = 0.0D0 
        K = M3-MNZ+1 - BORDER/2
        DO 100 ITH3 = 1,BORDER
           K=K+1
           IF(K<=0.OR.K>LZ)GOTO 100
           IPT1=(K-1)*LY*LX
            VAL0A = AMASSI  * BS3(ITH3)
            VAL1A = AMASSI * CFACTX * BS3(ITH3)
            VAL2A = AMASSI * CFACTY * BS3(ITH3)
            VAL3A = AMASSI * CFACTZ * DBS3(ITH3)
!
            J = M2-MNY+1 - BORDER/2
!
            DO 200 ITH2 = 1,BORDER
               J=J+1
               IF(J<=0.OR.J>LY)GOTO 200
               IPT2=IPT1+(J-1)*LX
!
               VAL0= VAL0A * BS2(ITH2)
               VAL1= VAL1A * BS2(ITH2)
               VAL2= VAL2A * DBS2(ITH2)
               VAL3= VAL3A * BS2(ITH2)

!
               I = M1-MNX+1 - BORDER/2
!
               DO 300 ITH1 = 1,BORDER
                  I=I+1
                  IF(I<=0.OR.I>LX)GOTO 300
                  IPT3=IPT2+I
                  RHOI=RHO(IPT3)-RHO0
                  SCRIG = SCRIG + VAL0 * RHOI * BS1(ITH1)
!                                   force is negative of grad
                  F1 = F1 - VAL1 * RHOI * DBS1(ITH1)
                  F2 = F2 - VAL2 * RHOI * BS1(ITH1)
                  F3 = F3 - VAL3 * RHOI * BS1(ITH1)
!
300          CONTINUE
200        CONTINUE
100      CONTINUE
!
         FX = (U(1,1)*F1+U(1,2)*F2+U(1,3)*F3)
         FY = (U(2,1)*F1+U(2,2)*F2+U(2,3)*F3)
         FZ = (U(3,1)*F1+U(3,2)*F2+U(3,3)*F3)
         F(IATX) = F(IATX)-FX
         F(IATY) = F(IATY)-FY
         F(IATZ) = F(IATZ)-FZ
        ! Calculate force and torque on the rigid
         rfx=rfx+fx
         rfy=rfy+fy
         rfz=rfz+fz
         rtx=rtx+yi*fz-zi*fy
         rty=rty+zi*fx-xi*fz
         rtz=rtz+xi*fy-yi*fx
      ENDDO
      ENRIG=-CFACT*SCRIG
#ifdef MPI
        IF(SANDERSIZE > 1)THEN
!  Combining all node results
!
          TEMP1(1)=SCRIG
          TEMP1(2)=RFX
          TEMP1(3)=RFY
          TEMP1(4)=RFZ
          TEMP1(5)=RTX
          TEMP1(6)=RTY
          TEMP1(7)=RTZ
          call mpi_allreduce(MPI_IN_PLACE,temp1,7,MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ierr)
          SCRIG=TEMP1(1)
          RFX=TEMP1(2)
          RFY=TEMP1(3)
          RFZ=TEMP1(4)
          RTX=TEMP1(5)
          RTY=TEMP1(6)
          RTZ=TEMP1(7)
        ENDIF
#endif
      ENRIGM=-CFACT*SCRIG
      rigobj%energy=enrigm
      rigobj%force(1)=rfx
      rigobj%force(2)=rfy
      rigobj%force(3)=rfz
      rigobj%torq(1)=rtx
      rigobj%torq(2)=rty
      rigobj%torq(3)=rtz
      ! If the map object is movable, remove partial net force/torque to constraint atoms
      if(rigobj%movable)then
        ! Center of the constrained atoms
        xc=xc/natc
        yc=yc/natc
        zc=zc/natc
        ! Center difference between map and atoms
        XI=XC-CX-TR(1)
        YI=YC-CY-TR(2)
        ZI=ZC-CZ-TR(3)
        !  Torque on molecule corrected for center difference
        rtx=(rtx-yi*rfz+zi*rfy)
        rty=(rty-zi*rfx+xi*rfz)
        rtz=(rtz-xi*rfy+yi*rfx)
        t2=rtx*rtx+rty*rty+rtz*rtz
        ! Roataional inertial around the torque
        w=0.0d0
        if(t2>1.0d-6)then
          DO ID=1,NATC
            IAT=IDX(ID)
            IATX=3*IAT-2
            IATY=IATX+1
            IATZ=IATY+1
            AMASSI=AMASS(IAT)
            XI=X(IATX)-XC
            YI=X(IATY)-YC
            ZI=X(IATZ)-ZC
            ! distance to axis by out product
            xj=yi*rtz-zi*rty
            yj=zi*rtx-xi*rtz
            zj=xi*rty-yi*rtx
            w=w+amassi*(xj*xj+yj*yj+zj*zj)
          enddo
        endif
        if(w>1.0d-6)then
          ! correct for torque in product w=w/t2
          ! partial angular accelaration to be removal
          rtx=-xgamma*t2*rtx/w
          rty=-xgamma*t2*rty/w
          rtz=-xgamma*t2*rtz/w
        else
          rtx=0.0d0
          rty=0.0d0
          rtz=0.0d0
        endif
        ! partial acceleration to be removal
        rfx=-xgamma*rfx/rigobj%mass
        rfy=-xgamma*rfy/rigobj%mass
        rfz=-xgamma*rfz/rigobj%mass
        DO ID=1,NATC
          IAT=IDX(ID)
          IF(IAT<NSTA.OR.IAT>NEND)CYCLE
          IATX=3*IAT-2
          IATY=IATX+1
          IATZ=IATY+1
          AMASSI=AMASS(IAT)
          XI=X(IATX)-XC
          YI=X(IATY)-YC
          ZI=X(IATZ)-ZC
          fx=rty*zi-rtz*yi
          fy=rtz*xi-rtx*zi
          fz=rtx*yi-rty*xi
          F(IATX) = F(IATX)-AMASSI*(RFX+fx)
          F(IATY) = F(IATY)-AMASSI*(RFY+fy)
          F(IATZ) = F(IATZ)-AMASSI*(RFZ+fz)
        enddo
      endif
      RETURN
      END SUBROUTINE RIGIDFORCE

      SUBROUTINE EMAPGUID(IDEMP,AMASS,NAT,IDX)
!_________________________________________________________________
!     Initialize constraint parameters
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
!
!
!
      INTEGER IDEMP,NAT,IDX(*)
      _REAL_ AMASS(*)
      INTEGER IAT,I,NDATA
      _REAL_ TMASS,DELT
!
!
      TMASS=0.0D0
      DELT=0.d0
!  Loop over all active rigid domains
      DO I=1,NAT
        IAT=IDX(I)
        TMASS=TMASS+AMASS(IAT)
      ENDDO
!  Statistics of map object
      !CALL EMAPSTAT(IDEMP)
      NDATA=EMAPS(IDEMP)%LX*EMAPS(IDEMP)%LY*EMAPS(IDEMP)%LZ
      !DELT=(RREMAP(IDEMP)-AREMAP(IDEMP)*AREMAP(IDEMP)/NDATA)/NDATA
      IF(DELT>0.0D0)DELT=SQRT(DELT)
      !EMAPS(IDEMP)%FCONS=EMGUID/DELT/TMASS
      RETURN
      END SUBROUTINE EMAPGUID



      FUNCTION rigeng(pos)
!_______________________________________________________________________
!     emap energy for rigid fitting
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
!
      use memory_module, only: xx =>x
#  include "../include/memory.h"
      _REAL_ rigeng,rigengi,RIGSCI
      _REAL_ pos(6)
!
      integer i
      
       do i=1,3
         rigmin%tran(i)=pos(i)
         rigmin%eular(i)=pos(i+3)
       enddo
       call eularrot(rigmin%rot,pos(4),pos(5),pos(6))
       CALL RIGIDENG(NATOM,RIGENGI,RIGSCI,xx(lmass),xx(lcrd),RIGMIN)
       RIGENG=RIGMIN%ENERGY
      return
      end function rigeng

      FUNCTION rigengd(pos,dpos)
!_______________________________________________________________________
!     emap energy and derivatives to x,y,z, and phi,psi,theta
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!
!     dEmap/dr=-f    r=(x,y,z)
!     dEmap/da=(dEmap/drc)(drc/da)=-T(d(A*r0)/da)=-T(dA/da)(A-)rc
!     a=(phi,psi,theta)=(f,y,s), rc=(1,1,1), T=(Tx,Ty, Tz)
!_________________________________________________________________
!
!

      use memory_module, only: xx =>x
#  include "../include/memory.h"
      _REAL_ rigengd,rigengi,rigsci
      _REAL_ pos(6),dpos(6)
!
      integer i
      _REAL_ a(3,3),daf(3,3),day(3,3),das(3,3),b(3),t(3),c(3)
      
       do i=1,3
         rigmin%tran(i)=pos(i)
         rigmin%eular(i)=pos(i+3)
       enddo
       call eulardrot(rigmin%rot,pos(4),pos(5),pos(6),daf,day,das)
       CALL RIGIDFORCE(NATOM,RIGENGI,RIGSCI,xx(lmass),xx(lcrd),xx(lforce),RIGMIN)
       RIGENGD=RIGMIN%ENERGY
       do i=1,3
         dpos(i)=-rigmin%force(i)
       enddo
       a=rigmin%rot
       t=rigmin%torq
      b(1)=a(1,1)+a(2,1)+a(3,1)
      b(2)=a(1,2)+a(2,2)+a(3,2)
      b(3)=a(1,3)+a(2,3)+a(3,3)
      c=matmul(daf,b)
      dpos(4)=-dot_product(t,c)
      c=matmul(day,b)
      dpos(5)=-dot_product(t,c)
      c=matmul(das,b)
      dpos(6)=-dot_product(t,c)
      return
      end function rigengd


      subroutine amotry(simtry,p,y,psum,mp,np,ndim,funk,ihi,fac)
!_____________________________________________
!     Trypoint converted from F77 numerical recipes: amotry.for
!______________________________________________
      INTEGER ihi,mp,ndim,np,NMAX  
      _REAL_ simtry,fac,p(mp,np),psum(np),y(mp),funk
      PARAMETER (NMAX=20)  
      EXTERNAL funk  
!    USES funk  
      INTEGER j  
      _REAL_ fac1,fac2,ytry,ptry(NMAX)  
      fac1=(1.-fac)/ndim  
      fac2=fac1-fac  
      do 11 j=1,ndim  
        ptry(j)=psum(j)*fac1-p(ihi,j)*fac2  
11    continue  
      ytry=funk(ptry)  
      if (ytry.lt.y(ihi)) then  
        y(ihi)=ytry  
        do 12 j=1,ndim  
          psum(j)=psum(j)-p(ihi,j)+ptry(j)  
          p(ihi,j)=ptry(j)  
12      continue  
      endif  
      simtry=ytry  
      return  
      END subroutine amotry
 

      SUBROUTINE simplex(p,y,mp,np,ndim,ftol,funk,iter)  
!_____________________________________________
!     Downhill simplex method converted from F77 numerical recipes: amoeba.for
!______________________________________________
      INTEGER iter,mp,ndim,np,NMAX,ITMAX  
      _REAL_ ftol,p(mp,np),y(mp),funk
      PARAMETER (NMAX=20,ITMAX=10000)  
      EXTERNAL funk
!    USES simtry,funk  
      INTEGER i,ihi,ilo,inhi,j,m,n  
      _REAL_ rtol,sum,swap,ysave,ytry,psum(NMAX)
      iter=0  
1     do 12 n=1,ndim  
        sum=0.  
        do 11 m=1,ndim+1  
          sum=sum+p(m,n)  
11      continue  
        psum(n)=sum  
12    continue  
2     ilo=1  
      if (y(1).gt.y(2)) then  
        ihi=1  
        inhi=2  
      else  
        ihi=2  
        inhi=1  
      endif  
      do 13 i=1,ndim+1  
        if(y(i).le.y(ilo)) ilo=i  
        if(y(i).gt.y(ihi)) then  
          inhi=ihi  
          ihi=i  
        else if(y(i).gt.y(inhi)) then  
          if(i.ne.ihi) inhi=i  
        endif  
13    continue  
      rtol=2.*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))  
      if (rtol.lt.ftol) then  
        swap=y(1)  
        y(1)=y(ilo)  
        y(ilo)=swap  
        do 14 n=1,ndim  
          swap=p(1,n)  
          p(1,n)=p(ilo,n)  
          p(ilo,n)=swap  
14      continue  
        return  
      endif  
      if (iter.ge.ITMAX) then
         !WRITE(6,*) 'ITMAX exceeded in simplex'  
         return
      endif
      iter=iter+2  
      !ytry = simtry(p,y,psum,mp,np,ndim,funk,ihi,-1.0d0)  
      call amotry(ytry,p,y,psum,mp,np,ndim,funk,ihi,-1.0d0)  
      if (ytry.le.y(ilo)) then  
        !ytry = simtry(p,y,psum,mp,np,ndim,funk,ihi,2.0d0)  
        call amotry(ytry,p,y,psum,mp,np,ndim,funk,ihi,2.0d0)  
      else if (ytry.ge.y(inhi)) then  
        ysave=y(ihi)  
        !ytry = simtry(p,y,psum,mp,np,ndim,funk,ihi,0.5d0)  
        call amotry(ytry,p,y,psum,mp,np,ndim,funk,ihi,0.5d0)  
        if (ytry.ge.ysave) then  
          do 16 i=1,ndim+1  
            if(i.ne.ilo)then  
              do 15 j=1,ndim  
                psum(j)=0.5*(p(i,j)+p(ilo,j))  
                p(i,j)=psum(j)  
15            continue  
              y(i)=funk(psum)  
            endif  
16        continue  
          iter=iter+ndim  
          goto 1  
        endif  
      else  
        iter=iter-1  
      endif  
      goto 2  
      RETURN
      END SUBROUTINE simplex

  

end module emap
