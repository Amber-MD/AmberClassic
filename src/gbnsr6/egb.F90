! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
#include "pb_def.h"
#include "timer.h"

module genborn

   implicit none
#  include "pb_constants.h"
   _REAL_, parameter :: alpb_alpha = 0.571412d0 !Alpha prefactor for ALPB
   integer gb_alpb 
   _REAL_  gb_B ! B parameter 
   integer gb_dgij ! print pair wise dgij if  

   _REAL_, allocatable :: reff(:)
   _REAL_, allocatable :: onereff(:)
  
   logical, allocatable :: skipv(:)
 
   _REAL_ molecule_mass
   _REAL_ x_cm
   _REAL_ y_cm
   _REAL_ z_cm
   _REAL_, allocatable :: Xcm(:), Ycm(:), Zcm(:) 
   _REAL_ gb_arad

   ! Computation of chunks used for AR6
   logical gb_chunk
   integer gb_depth
   integer copy_natom 
   _REAL_, allocatable  ::  chk_packing(:)  ! this array stores the chunk factors
   _REAL_, allocatable  ::  chk_crd(:,:)
   _REAL_, allocatable  ::  copyacrd(:,:)
   _REAL_, allocatable  ::  copyradi(:)
   _REAL_, allocatable  ::  chk_radi(:)
   integer, allocatable ::  chk_nex(:)
   integer, allocatable ::  chk_iex(:,:) 
   integer, allocatable ::  chk_natex(:) 
   integer,  allocatable :: chk_nshrt(:)
  
   _REAL_ , allocatable :: ar6_neck(:) ! neck parameters for ar6
   _REAL_ , allocatable :: ar6_facts(:)   ! facts parameters for ar6
   _REAL_ gb_Rs, gb_tau, gb_ROH ! Rs is the shift to the VDW surf. and
                                !subtracted from probe, tau is the scaling 
                                !for neighbors participation in CHA,
                                !ROH is propensity of explicit model CHA 
  integer  gb_cha ! CHAGB flag = 1, for canonical GB =0 
contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ generalized Born nonbonded routine place holder
subroutine egb( natom,nres,ntypes,npdec,ipres,iac,ico,natex,cn1,cn2,cg,x,f,enb,eel,epol, nbonh, ntbon, ibh, jbh ,isymbl )


   !use solvent_accessibility, only : dprob, radi, radip, radip2, radip3, nzratm, &
   !                                 sa_init, sa_driver, sa_free, sa_free_mb

   use solvent_accessibility, only : dprob, iprob, radi, radip, radip2, radip3,nzratm, &
         narcdot,maxarc,marc,m2narc,fstarc,arcatm,arccrd,savarc,dotarc, &
         sa_init, sa_driver, sa_free, sa_free_mb, cha_rad, radiopt

   use variable_module 
   use decomp, only : irespw, jgroup
   use pbtimer_module
 

#  include "pb_md.h"
#  include "md.h"
#  include "box.h"

   integer, parameter :: mytaskid = 0
   integer, parameter :: numtasks = 1

#  include "extra.h"
#  include "files.h"
 
   ! passed variables
   character (len=4) :: isymbl(*)
 
   integer natom, nres, ntypes, npdec, ipres(*), iac(*), ico(*), natex(*)
   _REAL_ cn1(*), cn2(*), cg(natom), x(3,natom), f(3,natom)
   _REAL_ enb, eel, epol    
   integer nbonh, ntbon, ibh(ntbon), jbh(ntbon)
   ! Local variables
   integer cnt, i, j, k
   integer iatm, jatm,  proatm, atmfirst, atmlast
   integer atmind(natom)
   _REAL_ acg(natom)
   _REAL_ pbcutcap, pbxcap, pbycap, pbzcap
   _REAL_ eelrffd, eelrfmp
   _REAL_ pbfrc(3,natom)
   _REAL_ ionene
   _REAL_ rh

   ! used to store chunk radii
   _REAL_ chk_rinvi 
   _REAL_ xi,yi,zi, dij, ri, rj, r2, ri2, ri4, rj2
   _REAL_ xij,yij,zij,v0, v1,v2,v3,v4, uij
   _REAL_ Ivdw, Iexact, ri3
   character(100):: line
   integer io

   ! Local multiblock variables
   integer ipermute     !
   integer iblock       !
!   integer boolsorted
   integer mynblock     !
   integer lvl2begin    !
   integer orphanblkn   !
   integer guess_int    !
!#endif
   integer ierr         !
   integer taskpiece    !
   integer myblkstart   !
   integer myblkend     !
   integer tasknatom    !
   integer ihavedone    !
   ! readin info verification
   integer tmp_nsatm    !
   integer mingrdblk    !
   integer myldim       !
   logical indexmatched !

   !This is not an efficient way to use the memory, need to estimate
   !the size of the array.
   integer, allocatable :: tmpindex(:) !

   _REAL_ myh           !
   _REAL_ guess_float   !
   _REAL_ blk_eelrffd   ! temporary storage for eelrffd
   _REAL_ blk_eel       ! temporary storage for eel
   _REAL_ blk_enb       ! temporary storage for enb

    ! end of multi-block
   logical localpbgrid


   ! Allocate, this needs to be done outside outside (gb_init)
   call gb_init(natom) 


   ! Variables initialization
   enb = ZERO; eel = ZERO; epol = ZERO
   eelrffd = ZERO; eelrfmp = ZERO
   pbfrc = ZERO; ionene = ZERO
   atmind = 0 

  ! Variables initialization, multi-block
   myblkstart = 0
   myblkend   = 0
   ierr = 0
   acg = ZERO !making sure clien(s) have clean acg

   firstleveldone = .false.

   ! assuming ifcap == 0;
   pbcutcap = ZERO; pbxcap = ZERO; pbycap = ZERO; pbzcap = ZERO


   ! split atoms into internal/external and ipdate nblis
   call pbtimer_start(PBTIME_PBLIST)
   outflag = 0;

 
   ! save coordinates in acrd  and multiply chages by 18.
   call gb_atmconv( natom,x,cg,acg )
   !call pb_atmconv(mpopt,ifcap,natom,ibgwat,ienwat,ibgion,ienion,atmind,ipres,x,cg,acg)

   ! I do not know if we need this part, seems that generate a lis of
   ! non bonded atoms for md and for docking???
   if ( ntnba == 1 .and. max(cutnb,cutsa,cutfd) > ZERO ) call pb_atmlist(pbverbose,pbprint,&
      maxnba,natom,ntypes,iac,ico,natex,nshrt,nex,iex,iar1pb,iprshrt,&
      cutnb,cutsa,cutfd,cn1,cn2,cn1pb,cn2pb,cn3pb,cg,acrd(1,1))
   if ( ntnbr == 1 ) ntnbr = 0
   if ( ntnba == 1 ) ntnba = 0
   call pbtimer_stop(PBTIME_PBLIST)

    
   ! set up the grid, ligand should false, no multiblock and ifcap==0
   if ( mpopt /=2 .and. pbgrid ) then
         call pb_setgrd(ipb,pbverbose,pbprint,pbinit,pbgrid,ifcap,1,natom,pbxcap,pbycap,pbzcap,pbcutcap)
   end if


   if ( srsas  ) then  ! ifcapp == 0
         call sa_init(pbverbose,pbprint,natom,natom,ifcap,dprob,radi,radip,radip2,outflag)
   !write(600+mytaskid,*)ipb,inp,natom,natom,dosas,ndosas,npbstep,nsaslag,ligand,multiblock,outflag;call mexit(0,0)
   !write(6,*)ipb,inp,natom,natom,dosas,ndosas,npbstep,nsaslag,ligand,multiblock,outflag;call mexit(0,0)
         call sa_driver(pbverbose,pbprint,ipb,inp,natom,natom,dosas,ndosas,npbstep,nsaslag,&
                        (ligand .or. multiblock), outflag,&
                        acrd(1,1),iar1pb(1,0),iprshrt,nex,iex,.false.)
   end if


   
   call gb_numerical_r6(pbverbose,pbprint,pbgrid,ifcap,ipb,imin,natom,natom,npdec,idecomp,irespw, &
                       ipres,jgroup,ibgwat,ibgion,pbfrc,eelrffd,ionene,npbstep,npbgrid,nstlim)

   ! Add the offset paramber B
if (gb_cha==1) then
      ! Define B
call gb_elsize( natom, acrd, radi, gb_arad  )
      if (gb_arad < 10) then
        gb_B = 0.0
      else
!        gb_B = 0.04
       gb_B = 0.0015* gb_arad + 0.01
      endif
endif
   do i=1,natom
      onereff(i) = onereff(i)  + gb_B
      reff(i) = ONE/onereff(i)  
   enddo

   ! print if necessary 
   if (rbornstat .eq. 1 ) then
     write(6,'(a)') "Print:    Atom number      Inverse Effective Born Radii"
      do i=1,natom
         write(6,'(a,i8,a,f14.8)') "rinv    " , i,"        ",  onereff(i)
      end do
      write(6, '(a)'        )
   end if

   if (gb_cha == 0) then
      ! computing electrostatic size
      if(gb_alpb == 1) call gb_elsize( natom, acrd, radi, gb_arad  )
     ! ALPB or Canonical GB 
      call gb_equation( gb_alpb, gb_arad, epol, eel, natex,nshrt  ) 
   
   elseif (gb_cha == 1) then
      ! computing electrostatic size
      if(gb_alpb == 1) call gb_elsize( natom, acrd, radi, gb_arad  )
      !THis is AM et. al. CHA-GB model
      call chagb_equation(gb_arad, epol, eel, natex,nshrt  )
   endif


   ! Compute Solvent Accesible Surface Area  (SASA)
   !level =1 ;
   !bcopt = savbcopt(level)
   !xm = savxm(level); ym = savym(level); zm = savzm(level)
   !xmym = xm*ym; xmymzm = xmym*zm
   !h = savh(level)
   !gox = savgox(level); goy = savgoy(level); goz = savgoz(level)
  
   !rh = ONE/h  ! rh is local variable
   !do iatm = atmfirst, atmlast
   !   gcrd(1,iatm) = (acrd(1,iatm) - gox)*rh
   !   gcrd(2,iatm) = (acrd(2,iatm) - goy)*rh
   !   gcrd(3,iatm) = (acrd(3,iatm) - goz)*rh
   !end do
   
   ! sasopt = 1, for solvent accesible option   
   !call pb_exmol_ses( .false., 0, ipb,savbcopt,saopt, 1,natom,&
   !     smoothopt,dprob,epsin,epsout,&
   !     h,gox,goy,goz,xm,ym,zm,xmymzm,&
   !     level,nfocus,&
   !     narcdot,maxarc,nbnd,nbndx,nbndy,nbndz,&
   !     outflag,gcrd,acrd,radi,radip3,&
   !     marc,m2narc,fstarc,arcatm,dotarc,arccrd,savarc,&
   !     atmsas,insas,lvlset,zv,epsx,epsy,epsz,&
   !     iepsav,iepsavx,iepsavy,iepsavz,fedgex,fedgey,fedgez)

   !call  calc_sa1(acrd,xm,ym,zm,xmymzm,nbndx,nbndy,nbndz,iepsavx,iepsavy,&
   !                 iepsavz,gox,goy,goz,h,fedgex,fedgey,fedgez,1)


   
   call pbtimer_start(PBTIME_PBSETUP)
   if ( srsas .and. (ifcap == 0 .or. ifcap == 5) ) then
      call sa_free( dosas,ndosas,.false. )
   end if
   call pbtimer_stop(PBTIME_PBSETUP)

   atmlast = natom

   !write(6,*) "atom name ", newtop, mdout 
   if ( newtop  .ne. ' ') then 

      ! change natex and nshrt, for chunk computation
      call gb_chunk_setup(gb_depth,natom,ntbon,ibh,jbh,nex,iex,natex,nshrt)  
      ! setting up necks and facts parameter per atoms 
      call ar6_necks_facts(natom,nbonh,ntbon,ibh,jbh,radi,ar6_neck,ar6_facts)

      copyacrd(1,:) =   acrd(1,:)    
      copyacrd(2,:) =   acrd(2,:)    
      copyacrd(3,:) =   acrd(3,:)    
      copyradi = radi
      copy_natom = natom 

      do iatm=1,copy_natom
      !iatm = 1
      natom = nex(iatm) + 1 
      
      ! step 1 : get coordinates
      ! loop over exclusion pairs:
      acrd(1,1)= x(1,iatm); acrd(2,1)=x(2,iatm); acrd(3,1)=x(3,iatm)
      radi(1) = copyradi(iatm )
      do j = 1, nex(iatm) 
         jatm = iex(j,iatm)
         acrd(1,j+1) = x(1,jatm ) 
         acrd(2,j+1) = x(2,jatm )
         acrd(3,j+1) = x(3,jatm )
         radi(j+1) = copyradi(jatm )
      end do
      
      ! step 2: setup the nonbonded lists
      cnt = 0;
      chk_nshrt(0) = 0 
      do i = 1, natom
          do j = i+1 , natom
             cnt = cnt + 1
             chk_natex(cnt) = j
          end do 
          chk_nshrt(i) = cnt 
      end do
      chk_nshrt(i-1) = cnt+1;
      chk_natex(cnt+1) = 0;

      call pb_atmlist(pbverbose,pbprint,&
      maxnba,natom,ntypes,iac,ico,chk_natex, chk_nshrt,chk_nex,chk_iex,iar1pb,iprshrt,&
      cutnb,cutsa,cutfd,cn1,cn2,cn1pb,cn2pb,cn3pb,cg,acrd(1,1))  

      ! step 3:  setup the grid of the small chunk    
      pbinit = .false.
      call pb_setgrd(ipb,pbverbose,pbprint,pbinit,pbgrid,ifcap,1,natom,pbxcap,pbycap,pbzcap,pbcutcap)
   
      ! step 4, create arcs and points of sasa
      call sa_init(.false.,.false.,natom,natom,ifcap,dprob,radi,radip,radip2,outflag)

      call sa_driver(.false.,.false.,ipb,inp,natom,natom,dosas,ndosas,npbstep,nsaslag,&
                        .false., outflag,&
                        acrd(1,1),iar1pb(1,0),iprshrt,chk_nex,chk_iex,.false.)

      ! step 5: find the boundary points and boundar edges
      call gb_chunk_radius( pbverbose,pbprint,pbgrid,ifcap,ipb,imin,natom,natom,npdec,idecomp,irespw, &
                     ipres,jgroup,ibgwat,ibgion,pbfrc,eelrffd,ionene,npbstep,npbgrid,nstlim,chk_rinvi ) 
      call sa_free( dosas,ndosas,.false. )

     
      ! find Ivdw : R6 integral over the atomic spheres. 
      xi = acrd(1,1) ! acrd is defined by the module variable_module
      yi = acrd(2,1)
      zi = acrd(3,1)
      ri = radi(1)
      ri2 =  ri*ri
      ri4 = ri2 * ri2
      Ivdw = 0.0d0
      do j = 2, natom

         xij = xi - acrd(1,j)
         yij = yi - acrd(2,j)
         zij = zi - acrd(3,j)
         rj = radi(j)
         rj2 = rj*rj
        
         r2 = xij*xij + yij*yij + zij*zij         
         dij=sqrt(r2)
        
         !Here it goes Grycuk equations for chunks
         if( dij > ri + rj ) then
            uij =  rj / ( r2 - rj2)
            Ivdw   = Ivdw  + uij * uij * uij

         else if ( dij > abs(ri - rj) ) then
            v0 = 16.0d0 * dij
            v1 = dij + 3.0d0 * rj
            v2 = dij + rj
            v3 = 3.0d0*(rj2 -2.0d0*ri2 -r2) + 8.0d0*dij*ri
            v4 = ri4
            Ivdw = Ivdw  +  (1.0d0/v0)*(v1/(v2*v2*v2)+v3/v4)

         else if (rj > ri ) then
            uij =  rj / ( r2 - rj2)
            Ivdw = Ivdw + one/(ri2*ri)  +  uij*uij*uij
         end if

      enddo
      ri3 = ri2 * ri
      Iexact  =   1.0d0/ri3 - chk_rinvi
      chk_packing(iatm) =  Iexact/Ivdw
!     write(70,*) "chunk", iatm ,   Iexact/Ivdw   
      enddo

      open(8,file=parm) !open topology file 
      open(70,file=newtop) ! open new topology file

      ! read all lines
      do
        read(8,'(A)',iostat=io) line
        if (io<0) exit
        write(70,'(A)')  trim(line)
      enddo
      close(8)

     
      write(70,'(A16)') '%FLAG CHUNKS_AR6'
      write(70,'(A15)') '%FORMAT(5E16.8)'
      write(70,100)   chk_packing
      
      write(70,'(A15)') '%FLAG NECKS_AR6'
      write(70,'(A15)') '%FORMAT(5E16.8)'
      write(70,100)   ar6_neck

      write(70,'(A16)') '%FLAG VOLUME_AR6'
      write(70,'(A15)') '%FORMAT(5E16.8)'
      write(70,100)   ar6_facts

      
      100 FORMAT(5ES16.8)

      close(70)

   end if

   call gb_free()


contains


! output arrays : arcd (atomic cordinates) 
!                 acrg (atomic charge cg/18.22223  charge units?  )
!                 acg  (charges in amber formats ) 
subroutine gb_atmconv( natom,x,cg,acg )

   ! Passed variables
   integer natom 
   _REAL_ x(3,natom), cg(natom), acg(natom)

   !local variables
   integer iatm

   ! acrd, acrg, and mapout are defined in poisson boltzmann module
   do iatm = 1, natom
       acrd(1,iatm)=x(1,iatm); acrd(2,iatm)=x(2,iatm); acrd(3,iatm)=x(3,iatm)
       acrg(iatm) = cg(iatm)/18.2223d0; acg(iatm) = cg(iatm)
       mapout(iatm) = iatm
   end do
end subroutine gb_atmconv



subroutine pb_atmconv( mpopt,ifcap,natom,ibgwat,ienwat,ibgion,ienion,atmind,ipres,x,cg,acg )

   ! output arrays : arcd (atomic cordinates) 
   !                 acrg (atomic charge cg/18.22223  charge units?  )
   !                 acg  (charges in amber formats ) 
   !                 mapout  = iatm ?

   ! Passed variables

   integer mpopt, ifcap, natom, ibgwat, ienwat, ibgion, ienion
   integer atmind(natom), ipres(*)
   _REAL_ x(3,natom), cg(natom), acg(natom)

   ! Local variables

   integer i, j, ifirst, ilast, iatm, ires, num
 
      do iatm = 1, natom
         acrd(1,iatm) = x(1,iatm); acrd(2,iatm) = x(2,iatm); acrd(3,iatm) = x(3,iatm)
         acrg(iatm) = cg(iatm)/18.2223d0; acg(iatm) = cg(iatm)
         mapout(iatm) = iatm
      end do

end subroutine pb_atmconv

! Computation of epol  using ALPB (alpb==1) or original still equation 
subroutine gb_equation( alpb, Arad, epol, eel  ,natex,nshrt   )

    !use gb_constants, only: alpb_alpha 

    ! Passed variables
    integer alpb
    _REAL_ Arad, epol

    integer nshrt(0:natom)
    integer natex(*)
    _REAL_ eel

    ! Local variables
    integer i, j, iatm, k, jjv 

    _REAL_ extdiel, intdiel, extdieli,intdieli
    _REAL_ onekappa, kappa
    _REAL_ xi, yi, zi, qi, ri 
    _REAL_ xij, yij, zij, r2
    _REAL_ alpb_beta, one_Arad_beta
    _REAL_ four_ri 
    _REAL_ expmkf, dl, qi2h, qid2h, temp1, fgb
    _REAL_ qiqj
    _REAL_ self_e, e

    ! epsin, epsout are multiplyied by eps0
    intdiel = epsin/eps0 
    extdiel = epsout/eps0
    kappa = pbkappa 
 
    if(alpb == 1) then
!     Sigalov Onufriev ALPB (epsilon-dependent GB):
      alpb_beta = alpb_alpha*(intdiel/extdiel)
      extdieli = ONE/(extdiel*(ONE + alpb_beta))
      intdieli = ONE/(intdiel*(ONE + alpb_beta))
      one_Arad_beta = alpb_beta/Arad
      if (kappa/=ZERO) onekappa = ONE/kappa
   else
   !  Standard Still's GB - alpb=0
      extdieli = ONE/extdiel
      intdieli = ONE/intdiel
      one_Arad_beta = ZERO
   end if 
   
   epol = ZERO 
   do i = 1, natom
      xi = acrd(1,i) ! acrd is defined by the module variable_module
      yi = acrd(2,i)
      zi = acrd(3,i)
      qi =  acg(i)   ! charge in amber format (poisson_boltzman)
      ri =  reff(i)  ! effective radii ( calculated by gb_numerical_r6 )
      four_ri = FOURTH*onereff(i)

      expmkf = exp( -kappa * reff(i) )*extdieli
      dl = intdieli - expmkf
      qi2h = HALF*qi*qi
      qid2h = qi2h * dl

      temp1 = (onereff(i) + one_Arad_beta)
      self_e = qid2h*temp1
      
      epol = epol - self_e

      ! print self terms
      if ( gb_dgij > 0) then
          write(6,'(a,2i8,f16.10)') "DGij ", i, i , -self_e
      endif

      
      !check the exclusion list for eel and vdw:
      do k=i+1,natom
         skipv(k) = .false.
      end do
      ! mark excluded atoms
      do jjv  = nshrt(i-1) + 1, nshrt(i)
         if ( natex(jjv) == 0 ) cycle
         skipv(natex(jjv))= .true.
         !write(6,*) natex(jjv) 
      end do


      do j=i+1, natom

         xij = xi - acrd(1,j)
         yij = yi - acrd(2,j)
         zij = zi - acrd(3,j)

         r2 = xij*xij + yij*yij + zij*zij
         qiqj = qi * acg(j) 

         temp1 = exp(-r2*four_ri*onereff(j) )
         fgb = sqrt( r2 + reff(i)*reff(j)*temp1 )

         if( kappa == ZERO ) then
            expmkf = extdieli
         else
            expmkf = exp(-kappa * fgb) *extdieli
         end if
         
         dl = intdieli - expmkf
         temp1 = -dl*( ONE/fgb + one_Arad_beta)
         e = qiqj*temp1
         
         epol = epol + e 

         ! Gas coulomb interaction
         if ( .not.  skipv(j)  ) then
             ! fixed intdiele_inv --> intdiel by Dan Parkin from Waseda Uni.
             temp1 =  intdiel * sqrt(r2) 
             eel = eel +  qiqj/temp1

             ! This is for printing 
             if ( gb_dgij .eq. 2 ) then
                 e = e + qiqj/temp1
             endif

         end if

         ! print the pair wise terms of elstrostatic 
         if ( gb_dgij > 0) then
            write(6,'(a,2i8,f12.4)') "DGij ", i, j , e
         endif
 
      end do 
   end do

end subroutine gb_equation


! Computation of epol  using CHAGB (alpb==1) or original still equation 
subroutine chagb_equation (Arad, epol, eel  ,natex,nshrt   )


    ! Passed variables
    _REAL_ Arad, epol
    integer nshrt(0:natom)
    integer natex(*)
    _REAL_ eel

    ! Local variables
    integer i, j, iatm, k, jjv

    _REAL_ extdiel, intdiel, extdieli,intdieli
    _REAL_ onekappa, kappa
    _REAL_ xi, yi, zi, qi, ri 
    _REAL_ xij, yij, zij, r2
    _REAL_ alpb_beta, one_Arad_beta
    _REAL_ four_ri 
    _REAL_ expmkf, dl, qi2h, qid2h, temp1, fgb
    _REAL_ qiqj
    _REAL_ self_e, e
    _REAL_ qeff(natom)
    _REAL_ mu_in_i,mu_in_j
    _REAL_ B


    ! Zero out
    mu_in_i = 0.0d0
    mu_in_j = 0.0d0
    qeff = 0.0d0

    ! epsin, epsout are multiplyied by eps0
    intdiel = epsin/eps0 
    extdiel = epsout/eps0
    kappa = pbkappa



   if (gb_alpb == 1) then
      ! Define B
 !     if (Arad < 140) then    
 !       B = 0.0 
 !     else
 !       B = 0.00000323* Arad +0.023
 !     endif
 !     reff = reff + B
!     Sigalov Onufriev ALPB (epsilon-dependent GB):
      alpb_beta = alpb_alpha*(intdiel/extdiel)
      extdieli = ONE/(extdiel*(ONE + alpb_beta))
      intdieli = ONE/(intdiel*(ONE + alpb_beta))
      one_Arad_beta = alpb_beta/Arad
      if (kappa/=ZERO) onekappa = ONE/kappa
   else
   !  Standard Still's GB - alpb=0
      extdieli = ONE/extdiel
      intdieli = ONE/intdiel
      one_Arad_beta = ZERO
   end if 

! Compute Q_eff
   do i = 1, natom
   qeff(i) = 0.0d0
    do j = 1, natom
    xij = acrd(1,i)-acrd(1,j)
    yij = acrd(2,i)-acrd(2,j)
    zij = acrd(3,i)-acrd(3,j)
    r2 = xij*xij + yij*yij +zij*zij
    qeff(i) = qeff(i)+acg(j)*exp(-gb_tau*r2/(reff(i)*reff(j)))
    enddo
   enddo
   epol = ZERO 
   do i = 1, natom
      xi = acrd(1,i) ! acrd is defined by the module variable_module
      yi = acrd(2,i)
      zi = acrd(3,i)
      qi =  acg(i)   ! charge in amber format (poisson_boltzman)
      ri =  reff(i)  ! effective radii ( calculated by gb_numerical_r6 )
      four_ri = FOURTH*onereff(i)

      expmkf = exp( -kappa * reff(i) )*extdieli
      dl = intdieli - expmkf
      qi2h = HALF*qi*qi
      qid2h = qi2h * dl
       call mu_inv(qeff(i), reff(i), mu_in_i)
      temp1 = (onereff(i)/mu_in_i + one_Arad_beta)
      self_e = qid2h*temp1
      
      epol = epol - self_e

      ! print self terms
      if ( gb_dgij > 0) then
          write(6,'(a,2i7,f12.4)') "DGij ", i, i , -self_e
      endif


      !check the exclusion list for eel and vdw:
      do k=i+1,natom
         skipv(k) = .false.
      end do
      ! mark excluded atoms
      do jjv  = nshrt(i-1) + 1, nshrt(i)
         if ( natex(jjv) == 0 ) cycle
         skipv(natex(jjv))= .true.
         !write(6,*) natex(jjv) 
      end do


      do j=i+1, natom

         xij = xi - acrd(1,j)
         yij = yi - acrd(2,j)
         zij = zi - acrd(3,j)

         r2 = xij*xij + yij*yij + zij*zij
         qiqj = qi * acg(j) 

         temp1 = exp(-r2*four_ri*onereff(j) )
         call mu_inv(qeff(j), reff(j), mu_in_j)
         fgb = sqrt( r2 + reff(i)*reff(j)*temp1*mu_in_i* mu_in_j)

         if( kappa == ZERO ) then
            expmkf = extdieli
         else
            expmkf = exp(-kappa * fgb) *extdieli
         end if
         
         dl = intdieli - expmkf
         temp1 = -dl*( ONE/fgb + one_Arad_beta)
         e = qiqj*temp1
         
         epol = epol + e 


         ! Gas coulomb interaction
         if ( .not.  skipv(j)  ) then
             ! fixed intdiele_inv --> intdiel by Dan Parkin from Waseda Uni.
             temp1 =  intdiel * sqrt(r2)
             eel = eel +  qiqj/temp1
             
             ! This is for printing 
             if ( gb_dgij .eq. 2 ) then
                 e = e + qiqj/temp1
             endif

         end if
         
         ! print the pair wise terms of elstrostatic 
         if ( gb_dgij > 0) then
            write(6,'(a,2i7,f12.4)') "DGij ", i, j , e
         endif

      end do 
   end do

end subroutine chagb_equation

! SUBROTINE FOR MU
subroutine mu_inv (q, R, mu_in)
!   use genborn, only: gb_tau, gb_ROH
  use solvent_accessibility, only: dprob 
     _REAL_ mu_in, q, R
    if(q .lt. 0) then
    mu_in = 1.0d0 - gb_ROH/(R+dprob)
    else if (q .gt. 0) then
    mu_in = 1.0d0 + gb_ROH/(R+dprob)
    else
    mu_in =1.0d0
    end if
end subroutine mu_inv

    
end subroutine egb 



subroutine gb_chunk_radius( pbverbose,pbprint,pbgrid,ifcap,ipb,imin,natom,atmlast,npdec,idecomp,irespw, &
                     ipres,jgroup,ibgwat,ibgion,pbfrc,eelrf,ionene,npbstep,npbgrid,nstlim, rinvi )
   use variable_module
   use solvent_accessibility, only : dprob,iprob,radi,radip3,nzratm, &
       narcdot,maxarc,marc,m2narc,fstarc,arcatm,arccrd,savarc,dotarc     

    ! passed variables

   logical pbverbose, pbprint, pbgrid
   integer ifcap, ipb, imin, natom, atmlast, npdec, idecomp, ibgwat, ibgion
   integer npbstep, npbgrid, nstlim
   integer irespw(*), ipres(*), jgroup(*)
   _REAL_ ionene, eelrf, pbfrc(3,natom)!, fnet(3)
   _REAL_ rinvi

   !local variables
   integer atmfirst
   integer iatm, lastatm, mpdec, i
   integer cnt, j, k ! phiform = 2
   _REAL_ rh
   !integer atmid
   !_REAL_ rinv
 
   
   atmfirst = 1
   ! atmfirst is supposedly to be passed
   m = 1
   level = 1

   ! retrieving saved grid data into working variables
   bcopt = savbcopt(level)
   xm = savxm(level); ym = savym(level); zm = savzm(level)
   xmym = xm*ym; xmymzm = xmym*zm
   h = savh(level)
   gox = savgox(level); goy = savgoy(level); goz = savgoz(level)

   rh = ONE/h
   do iatm = atmfirst, atmlast
      gcrd(1,iatm) = (acrd(1,iatm) - gox)*rh
      gcrd(2,iatm) = (acrd(2,iatm) - goy)*rh
      gcrd(3,iatm) = (acrd(3,iatm) - goz)*rh
   end do

   call pb_exmol_ses( pbverbose,ifcap,ipb,savbcopt,saopt,sasopt,natom,&
        smoothopt,dprob,epsin,epsout,&
        h,gox,goy,goz,xm,ym,zm,xmymzm,&
        level,nfocus,&
        narcdot,maxarc,nbnd,nbndx,nbndy,nbndz,&
        outflag,gcrd,acrd,radi,radip3,&
        marc,m2narc,fstarc,arcatm,dotarc,arccrd,savarc,&
        atmsas,insas,lvlset,zv,epsx,epsy,epsz,&
        iepsav,iepsavx,iepsavy,iepsavz,fedgex,fedgey,fedgez)

   call calc_chunk_inv(acrd,xm,ym,zm,xmymzm,nbndx,nbndy,nbndz,iepsavx,iepsavy,&
                    iepsavz,gox,goy,goz,h,fedgex,fedgey,fedgez,sasopt,natom, atmfirst , rinvi )

end subroutine gb_chunk_radius


subroutine gb_numerical_r6( pbverbose,pbprint,pbgrid,ifcap,ipb,imin,natom,atmlast,npdec,idecomp,irespw, &
                     ipres,jgroup,ibgwat,ibgion,pbfrc,eelrf,ionene,npbstep,npbgrid,nstlim )

   use variable_module
   use solvent_accessibility, only : dprob,iprob,radi,radip3,nzratm, &
       narcdot,maxarc,marc,m2narc,fstarc,arcatm,arccrd,savarc,dotarc
!#ifndef SANDER
   use pbtimer_module
!#endif


#  include "flocntrl.h"
#include "parallel.h"

   ! passed variables

   logical pbverbose, pbprint, pbgrid
   integer ifcap, ipb, imin, natom, atmlast, npdec, idecomp, ibgwat, ibgion
   integer npbstep, npbgrid, nstlim
   integer irespw(*), ipres(*), jgroup(*)
   _REAL_ ionene, eelrf, pbfrc(3,natom)!, fnet(3)

   ! local variables
   integer atmfirst
   integer iatm, lastatm, mpdec, i
   integer cnt, j, k ! phiform = 2
   _REAL_ eelself, eelcoul, rh, fcrd(3,atmlast)
   _REAL_ aa, bb, cc, aa1, bb1, cc1
   _REAL_ bb1cc1, bb_cc1, bb1cc, bb_cc
   _REAL_ acrgtmp, factor
   _REAL_ grdreac, drcreac

   !XP: test variables
   _REAL_  px0,py0,pz0,pu,pex,pey,pez,pcirreg(1:137440,1:6)
   integer ipp,iip
   character(len=20):: filename = "pcrd"
   logical alive
   integer :: status = 0

   character(len=12) phifilename
   character(len=23) phidataname
   integer           phifilenum

   if (phiform == 0 .OR. phiform == 1) then
      phifilename="pbsa_phi.phi"
   else
      phifilename="pbsa_phi.dx "
   end if
   phidataname="Electrostatic Potential"
   phifilenum=64
    
   !write(6,*) "FALGGGGGGGGG" ,  do_pbfd 
   if ( do_pbfd == 0 ) return

   mpdec = 1
   if ( idecomp > 2 ) mpdec = npdec

   !-- PB pairwise decomp <<1,mpdec>>

   atmfirst = 1
   ! atmfirst is supposedly to be passed
   m = 1 
   eelself = ZERO; eelcoul = ZERO 
   if ( .not. firstleveldone ) xsoffset = 1
   
   level = 1 

   call pbtimer_start(PBTIME_PBBUILDSYS)

   ! retrieving saved grid data into working variables

   bcopt = savbcopt(level)
   xm = savxm(level); ym = savym(level); zm = savzm(level)
   xmym = xm*ym; xmymzm = xmym*zm
   h = savh(level)
   gox = savgox(level); goy = savgoy(level); goz = savgoz(level)

   rh = ONE/h
   do iatm = atmfirst, atmlast
      gcrd(1,iatm) = (acrd(1,iatm) - gox)*rh
      gcrd(2,iatm) = (acrd(2,iatm) - goy)*rh
      gcrd(3,iatm) = (acrd(3,iatm) - goz)*rh
   end do

   do iatm = atmfirst, atmlast
      icrd(1,iatm) = floor(gcrd(1,iatm))
      icrd(2,iatm) = floor(gcrd(2,iatm))
      icrd(3,iatm) = floor(gcrd(3,iatm))
      fcrd(1,iatm) = dble(icrd(1,iatm))
      fcrd(2,iatm) = dble(icrd(2,iatm))
      fcrd(3,iatm) = dble(icrd(3,iatm))

      if ( icrd(1,iatm) > xm-1 .or. icrd(2,iatm) > ym-1 .or. icrd(3,iatm) > zm-1 ) then
         write(6,'(3f12.4)') acrd(1:3,iatm)
         write(6,'(a,3i6)') "pb_fdfrc(): Atom out of focusing box",icrd(1:3,iatm)
         call mexit(6,1)
      end if
   end do

   gcrg=0d0!shouldn't be necessary

   do iatm = atmfirst, atmlast
      !if ( level > 1 .and. ligand     .and. realflag(iatm) == 0 ) cycle
      !if ( level > 1 .and. multiblock .and. realflag(iatm) == 0 ) cycle
      aa = gcrd(1,iatm) - fcrd(1,iatm)
      bb = gcrd(2,iatm) - fcrd(2,iatm)
      cc = gcrd(3,iatm) - fcrd(3,iatm)
      bb1 = ONE - bb; cc1 = ONE - cc
      !-- PB decomp
      if (idecomp < 3 ) then
         acrgtmp = acrg(iatm)
      else if (iatm >= ipres(irespw(m)) .and. iatm < ipres(irespw(m)+1)) then
         acrgtmp = acrg(iatm)
      else
          acrgtmp = ZERO
      end if
      aa  = acrgtmp*aa; aa1 = acrgtmp - aa
      bb1cc1 = bb1*cc1; bb_cc1 = bb *cc1
      bb1cc  = bb1*cc ; bb_cc  = bb *cc
      if ( (ifcap == 2 .or. ifcap == 5) .and. outflag(iatm) == 1 ) then
         gcrg(1,iatm) = ZERO; gcrg(2,iatm) = ZERO
         gcrg(3,iatm) = ZERO; gcrg(4,iatm) = ZERO
         gcrg(5,iatm) = ZERO; gcrg(6,iatm) = ZERO
         gcrg(7,iatm) = ZERO; gcrg(8,iatm) = ZERO
      else
         gcrg(1,iatm) = aa1*bb1cc1; gcrg(2,iatm) = aa *bb1cc1
         gcrg(3,iatm) = aa1*bb_cc1; gcrg(4,iatm) = aa *bb_cc1
         gcrg(5,iatm) = aa1*bb1cc ; gcrg(6,iatm) = aa *bb1cc
         gcrg(7,iatm) = aa1*bb_cc ; gcrg(8,iatm) = aa *bb_cc
      end if
   end do
  
   call pb_exmol_ses( pbverbose,ifcap,ipb,savbcopt,saopt,sasopt,natom,&
        smoothopt,dprob,epsin,epsout,&
        h,gox,goy,goz,xm,ym,zm,xmymzm,&
        level,nfocus,&
        narcdot,maxarc,nbnd,nbndx,nbndy,nbndz,&
        outflag,gcrd,acrd,radi,radip3,&
        marc,m2narc,fstarc,arcatm,dotarc,arccrd,savarc,&
        atmsas,insas,lvlset,zv,epsx,epsy,epsz,&
        iepsav,iepsavx,iepsavy,iepsavz,fedgex,fedgey,fedgez)

!print *, "before NSR6"
!call mexit(6,1)

   call calc_NSR6(acrd,xm,ym,zm,xmymzm,nbndx,nbndy,nbndz,iepsavx,iepsavy,&
                iepsavz,gox,goy,goz,h,fedgex,fedgey,fedgez,sasopt,natom)  
 
   call pbtimer_stop(PBTIME_PBBUILDSYS) 


end subroutine gb_numerical_r6

subroutine gb_elsize( natom, acrd, radi, A_det  )

   ! passed variables
   integer natom
   _REAL_ acrd(3,*), radi(*)
   _REAL_ A_det

   ! local variables
   integer i
   integer ierr
   _REAL_ rad3
   _REAL_ I11, I12, I13, I22, I23, I33
   _REAL_ atom_mass, atom_MI
   _REAL_ x2, y2, z2 
   _REAL_ det_I
   _REAL_ sqrt_25_mass
   _REAL_ sixthroot

   ! allocation of arrays 
   allocate(Xcm(natom),stat=ierr); if (ierr /= 0) call mexit(6,1) 
   allocate(Ycm(natom),stat=ierr); if (ierr /= 0) call mexit(6,1) 
   allocate(Zcm(natom),stat=ierr); if (ierr /= 0) call mexit(6,1) 
   
   molecule_mass = 0.0d0
   x_cm = 0.0d0 
   y_cm = 0.0d0
   z_cm = 0.0d0

   do i = 1, natom
      
      rad3 = radi(i)*radi(i)*radi(i) 

      x_cm = x_cm + rad3*acrd(1,i)
      y_cm = y_cm + rad3*acrd(2,i)
      z_cm = z_cm + rad3*acrd(3,i)

      molecule_mass = molecule_mass + rad3
      
   end do

   !test uf molecule_mass is grater than zero
   x_cm  = x_cm /  molecule_mass 
   y_cm  = y_cm /  molecule_mass 
   z_cm  = z_cm /  molecule_mass

   do i = 1, natom
      Xcm(i) = acrd(1,i) - x_cm
      Ycm(i) = acrd(2,i) - y_cm
      Zcm(i) = acrd(3,i) - z_cm

   end do 
  

   I11 = 0.0d0
   I12 = 0.0d0
   I13 = 0.0d0
   I22 = 0.0d0
   I23 = 0.0d0
   I33 = 0.0d0
   do i = 1, natom
      atom_MI  = radi(i)*radi(i)
      atom_mass = atom_MI*radi(i)
      atom_MI = atom_MI / 5.0d0

      x2 = Xcm(i)*Xcm(i) + atom_MI 
      y2 = Ycm(i)*Ycm(i) + atom_MI
      z2 = Zcm(i)*Zcm(i) + atom_MI
  
      I11 = I11 + atom_mass*(y2 + z2)
      I22 = I22 + atom_mass*(z2 + x2)
      I33 = I33 + atom_mass*(x2 + y2)

      I12 = I12 - atom_mass*Xcm(i)*Ycm(i)
      I13 = I13 - atom_mass*Xcm(i)*Zcm(i)
      I23 = I23 - atom_mass*Ycm(i)*Zcm(i)


   end do

   det_I = I11*I22*I33 + 2*I12*I23*I13 - I11*I23*I23 &
           - I22*I13*I13 - I33*I12*I12
  
    A_det = 0;
    if (det_I  >  0.0d0 ) then
        sqrt_25_mass = sqrt(2.5d0/molecule_mass);
        sixthroot = det_I**(1.0d0/6.0d0)
        A_det = sqrt_25_mass*sixthroot;
 
    else    
        write(6,'(a,a,f12.4)') "problem with the determinant of the ", &
                    "tensor of intertia: Det I= ", det_I 

        call mexit(6,1) 
    end if


   deallocate(Xcm, stat = ierr )
   deallocate(Ycm, stat = ierr )
   deallocate(Zcm, stat = ierr )

end subroutine gb_elsize

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ computation of the chunks atoms
!+ inputs: 
!+    depth : chunk depth
!+    natom : number of atoms of the molecule
!+    ntbon :  number of covalent bonds
!+    ibh, jbh: arrays with indexes of the atoms of the bonds 
!+ outputs 
!+    nex, idex : arrays to store the chunks 
!+    natex     : array with "excluded atoms" similar to amber format
!+    nshrt     : array with indexes to read natex   
subroutine gb_chunk_setup( depth, natom, ntbon, ibh, jbh, nex,iex, natex, nshrt)

   implicit none

   !passed variables
   integer natom , ntbon, ibh(ntbon) , jbh(ntbon)  
   integer nex(natom) , iex(64,natom)  ! these will be updated
   integer depth ! this determina the chunk size
                 ! depth is the maximun distance allowed betwen
                 ! the atoms of the chunk and atom i
   integer natex(*), nshrt(0:natom)  ! will be updated

   !local variables
   logical chunks(natom)
   logical overlap(natom)
   integer iatm, ii, jj, i,j,k, jjv, cnt
   integer counter

   ! this will be chunk = 0
   nex(1:natom) = 1
   do iatm =1,natom  ; iex(1,iatm) = iatm ; enddo
  
   do  counter = 1,depth
   do iatm =1, natom ! loop over atoms

      chunks = .false.
      overlap = .false.
      
      do jjv = 1, nex(iatm) !loop over the chunk i
         chunks( iex(jjv, iatm) ) = .true.
      end do 

      do k = 1, ntbon  ! loop over covalent bonds
         !get the atomic indices of the bond k
         ii = (ibh(k) + 3 )/3
         jj = (jbh(k) + 3 )/3
         
         ! chek if ii or jj is bonded to the chunk
         if ( chunks(ii) .AND. (.NOT. chunks(jj)) ) then
             overlap(jj) = .true.
         elseif ( chunks(jj) .AND. (.NOT. chunks(ii)) ) then
             overlap(ii) = .true.
         end if
      end do

      ! loop over atoms to increment the chunk size
      nex(iatm) = 0 ;
      do  j = 1, natom
         if ( chunks(j) ) then ! original chunk atoms
            nex(iatm) = nex(iatm) + 1
            iex( nex(iatm) , iatm) = j
         else if( overlap(j) ) then ! new atoms 
            nex(iatm) = nex(iatm) + 1
            iex( nex(iatm) , iatm) = j
         end if
      end do  
       
   end do 
   end do

   ! change natex and nshrt accordint to chunk size
   ! the chunk size is determined by depth
   cnt = 0;
   nshrt(0) = 0
   do iatm = 1, natom
      i = cnt
      do j = 1 , nex(iatm) ! loop over atoms of the chunk_i
         if ( iex(j,iatm) > iatm ) then 
             cnt = cnt + 1
             natex(cnt) = iex(j,iatm)
         end if
      end do
      ! update nshrt (sizes) 
      if ( i .eq. cnt) then !if no atom found for chunk of i
            cnt = cnt + 1
            natex(cnt) = 0
            nshrt(iatm) = cnt
      else
            nshrt(iatm) = cnt
      end if
   end do

end subroutine gb_chunk_setup 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ setting up necks and factos for AR6 
!+ inputs: 
!+    natom : number of atoms of the molecule
!+    ntbon :  number of covalent bonds with hydrogens
!+    ibh, jbh: arrays with indexes of the atoms of the bonds 
!+    radi : the radii of atomos (just bondi and mbondi2 are expected
!+ outputs 
!+    ar6_necks 
!+    ar6_facts 
subroutine ar6_necks_facts(natom,nbonh, ntbon,  ibh,jbh,  radi,ar6_neck,ar6_facts)

   implicit none
#  include "ar6_constants.h"

   !passed variables
   integer natom , nbonh, ntbon,  ibh(ntbon) , jbh(ntbon) 

   _REAL_ ar6_neck(natom) , ar6_facts(natom)   ! will be updated
   _REAL_ radi(natom) ! atomic radii    

   !local variables
   integer iatm, k, ii, jj, cnt
   integer nha(natom)
   integer nbo(natom)
   _REAL_ temp
  
   nha = ZERO
   nbo = ZERO
   ! determine the number of bonded hydrogen atoms 
   do k = 1, nbonh  ! loop over covalent bonds
      !get the atomic indices of the bond k
      ii = (ibh(k) + 3 )/3
      jj = (jbh(k) + 3 )/3

      ! determine hydrogen ( with smaller radius)
      if ( radi(ii) > radi(jj) ) then
         nha(ii) = nha(ii) + 1
      else
         nha(jj) = nha(jj) + 1   
      end if

   end do

   ! determine the number of covalent bonds
   do k = 1, ntbon  ! loop over covalent bonds
      !get the atomic indices of the bond k
      ii = (ibh(k) + 3 )/3
      jj = (jbh(k) + 3 )/3
      
      !update total number of bonds
      nbo(ii) = nbo(ii) + 1
      nbo(jj) = nbo(jj) + 1
   end do


   ! put tha values to of necks and facts to each atom
   ! aton type are dinstinguised by its radii, needs improvement
   do iatm = 1 , natom
      if (radi(iatm) == 1.7d0) then  ! Carbon
         if (nbo(iatm) == 4) then
            if (nha(iatm) == 3) then
               ar6_neck(iatm)  = AR6_NC_3H
               ar6_facts(iatm) = AR6_SC_3H
            else
               ar6_neck(iatm)  = AR6_NC
               ar6_facts(iatm) = AR6_SC
            end if
         else 
            ar6_neck(iatm)  = AR6_NC_3B
            ar6_facts(iatm) = AR6_SC_3B
         end if

      else if (radi(iatm) == 1.55d0) then !nitrogen
         ar6_neck(iatm)  = AR6_NN
         ar6_facts(iatm) = AR6_SN 
      else if (radi(iatm) == 1.5d0) then !oxygen
         ar6_neck(iatm)  = AR6_NO
         ar6_facts(iatm) = AR6_SO
      else if (radi(iatm) == 1.3d0) then !H -N
         ar6_neck(iatm)  = AR6_NHH
         ar6_facts(iatm) = AR6_SHH
      else if (radi(iatm) == 1.2d0) then ! hydrogens
         ar6_neck(iatm)  = AR6_NH
         ar6_facts(iatm) = AR6_SH
      else if (radi(iatm) == 1.8d0) then ! sulfur
         ar6_neck(iatm)  = AR6_NS
         ar6_facts(iatm) = AR6_SS
      else if (radi(iatm) == 1.85d0) then ! phosphorous
         ar6_neck(iatm)  = AR6_NP
         ar6_facts(iatm) = AR6_SP
      else
         ar6_neck(iatm)  = 0.4056d0
         ar6_facts(iatm) = 0.0d0
      end if
   end do 

end subroutine ar6_necks_facts



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ GB cleanup routine
subroutine gb_free

   !use genborn


   implicit none

   integer  ierr

   deallocate(reff , stat = ierr )
   deallocate(onereff , stat = ierr )

   deallocate(ar6_neck, stat = ierr )
   deallocate(ar6_facts, stat = ierr )
   deallocate(skipv, stat = ierr )
 
 
   if (gb_chunk) then
      deallocate(chk_crd  , stat = ierr )
      deallocate(chk_radi  , stat = ierr )
      deallocate(chk_nex  , stat = ierr )
      deallocate(chk_iex , stat = ierr )
      deallocate(chk_packing , stat = ierr )

      deallocate( copyacrd , stat = ierr ) 
      deallocate( copyradi, stat = ierr )
      deallocate( chk_natex, stat = ierr ) 
      deallocate( chk_nshrt, stat = ierr )
   end if


end subroutine gb_free


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine gb_init(natom)

   ! Module variables

   !use genborn
   implicit none

   integer natom

   integer ierr 

   gb_chunk = .true. 

   allocate(reff(natom),    stat=ierr); if (ierr /= 0) call mexit(6,1)
   allocate(onereff(natom),    stat=ierr); if (ierr /= 0) call mexit(6,1)

   allocate(ar6_neck(natom), stat=ierr); if (ierr /= 0) call mexit(6,1) 
   allocate(ar6_facts(natom), stat=ierr); if (ierr /= 0) call mexit(6,1)    
   allocate(skipv(natom), stat=ierr); if (ierr /= 0) call mexit(6,1)    


   ! chunk computation maximun number fo atoms in chunk = 64
   if (gb_chunk ) then

     allocate(chk_crd(3,64),  stat=ierr); if (ierr /= 0) call mexit(6,1)
     allocate(chk_radi(64),  stat=ierr); if (ierr /= 0) call mexit(6,1)
     allocate( chk_nex(64),    stat=ierr); if (ierr /= 0) call mexit(6,1)
     allocate( chk_iex(64,64), stat=ierr); if (ierr /= 0) call mexit(6,1)
     allocate( chk_natex(2048), stat=ierr); if (ierr /= 0) call mexit(6,1)
     allocate( chk_nshrt(0:64), stat=ierr); if (ierr /= 0) call mexit(6,1) 

     allocate( copyacrd(3,natom) , stat=ierr); if (ierr /= 0) call mexit(6,1)
     allocate( copyradi(natom), stat=ierr); if (ierr /= 0) call mexit(6,1) 
     allocate( chk_packing(natom), stat=ierr); if (ierr /= 0) call mexit(6,1) 


     chk_nex = ZERO
     chk_iex =  ZERO

   end if 
 
   ! default parameters of genborn module
   !gb_alpb = 1
   !gb_arad = 30

end subroutine gb_init





end module genborn

