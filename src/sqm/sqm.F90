#include "copyright.h"
#include "../include/dprec.fh"
#include "../include/assert.fh"
!
!        SQM  stand-alone quantum program
!

program sqm

   use qmmm_module, only : qmmm_nml, qmmm_struct, qmmm_mpi, qm2_struct, &
                           qmmm_vsolv, qm2_params, deallocate_qmmm
   use sqm_qmmm_read_and_alloc, only : read_qmmm_nm_and_alloc
   use qm2_dftb_module, only : ks_struct
   use constants, only : KCAL_TO_EV, EV_TO_KCAL
   use qm2_pm6_hof_module, only : cct, nsp2, print, strlen
   use file_io_dat, only : MAX_FN_LEN

   use UtilitiesModule, only : print

   implicit none

   _REAL_ x(3000), f(3000), escf
   character(len=8) atnam(1000)
   _REAL_ born_radii(1000), one_born_radii(1000)
   _REAL_ intdiel, extdiel, Arad
   integer natom, ier, atnum(1000), xmin_iter
   character(len=80) arg ! temp for each of the command line arguments
   integer iarg !         index of the current argument
   integer last_arg_index !   index of the last argument
   integer ntpr
   character owrite
   character(len=MAX_FN_LEN) mdin, mdout 
   ! external charge
   _REAL_ excharge(4000)
   integer chgatnum(1000)
   character(len=8) chgnam(1000)
   integer ncharge

   integer :: igb, maxcyc
   _REAL_  :: grms_tol
   _REAL_  :: total_energy
   logical :: master=.true.


   character(len=strlen) :: string

!   interface
!      subroutine qm2_calc_dipole(coord,mass,ipres,lbres,nres)
!        integer, optional, intent(in) :: nres, ipres(*)
!        _REAL_, optional, intent(inout) :: mass(*)
!        _REAL_, intent(inout) :: coord(*)
!        character(len=4), optional, intent(in) :: lbres(*)
!      end subroutine qm2_calc_dipole
!   end interface


   ! ==== Initialise first_call flags for QMMM ====
   qmmm_struct%qm_mm_first_call = .true.
   qmmm_struct%fock_first_call = .true.
   qmmm_struct%fock2_2atm_first_call = .true.
   ! qmmm_struct%qm2_deriv_qm_analyt_first_call = .true.
   qmmm_struct%qm2_allocate_e_repul_first_call = .true.
   ! qmmm_struct%qm2_rotate_qmqm_first_call = .true.
   qmmm_struct%qm2_calc_rij_eqns_first_call = .true.
   qmmm_struct%qm2_scf_first_call = .true.
   qmmm_struct%zero_link_charges_first_call = .true.
   qmmm_struct%adj_mm_link_pair_crd_first_call = .true.

   !     --- default file names ---
   
   mdin   = 'mdin'
   mdout  = 'mdout'
   iarg = 0
   owrite = 'N'  ! output status: New
   last_arg_index = command_argument_count()
   do while (iarg < last_arg_index)

      iarg = iarg + 1
      call getarg(iarg,arg)

      if (arg == '-i') then
         iarg = iarg + 1
         call getarg(iarg,mdin)
      else if (arg == '-o') then
         iarg = iarg + 1
         call getarg(iarg,mdout)
      else if (arg == '-O') then
         owrite = 'R'   ! output status: Replace
      else if (arg == '-h') then
         write(6,'(a)') 'sqm [-O] -i <input> -o <output>'
         write(6,'(a)') '   -O           Overwrite output file if it exists'
         write(6,'(a)') '   -i <input>   Input file'
         write(6,'(a)') '   -o <output>  Output file'
         write(6,'(a)') '   -h           Show this message'
         call mexit(6, 0)
      else if (arg == ' ') then
         continue
      else
         write(0,'(/,5x,a,a)') 'Error unknown flag: ',arg
         call mexit(6, 1)
      end if 
   end do  !  while (iarg < last_arg_index)

   igb = 0
   call amopen(5,mdin,'O','F','R')
   call amopen(6,mdout,owrite,'F','W')

   write(6,*) '           --------------------------------------------------------'
   write(6,*) '                            AMBER SQM VERSION 19'
   write(6,*) ''
   write(6,*) '                                    By'
   write(6,*) '             Ross C. Walker, Michael F. Crowley, Scott Brozell,'
   write(6,*) '                        Tim Giese, Andreas W. Goetz,'
   write(6,*) '                       Tai-Sung Lee and David A. Case'
   write(6,*) ''              
   write(6,*) '           --------------------------------------------------------'
   write(6,*) ''                  

   call getsqmx( natom, x, atnam, atnum, ncharge, excharge, chgnam, chgatnum )
   call read_qmmm_nm_and_alloc(natom,igb,atnum,maxcyc,grms_tol,ntpr, &
                               ncharge,excharge,chgatnum )
   call qm_assign_atom_types

   ! Set default QMMM MPI parameters - for single cpu operation.
   ! These will get overwritten by qmmm_mpi_setup if MPI is on.
   ! qmmm_mpi%master = master
   qmmm_mpi%commqmmm_master = master
   qmmm_mpi%numthreads = 1
   qmmm_mpi%mytaskid = 0
   qmmm_mpi%natom_start = 1
   qmmm_mpi%natom_end = natom
   qmmm_mpi%nquant_nlink_start = 1
   qmmm_mpi%nquant_nlink_end = qmmm_struct%nquant_nlink
   call allocate_qmgb(qmmm_struct%nquant_nlink)

   allocate( qmmm_struct%dxyzqm(3, qmmm_struct%nquant_nlink), stat = ier )
   REQUIRE(ier == 0)

   allocate ( qm2_struct%scf_mchg(qmmm_struct%nquant_nlink), stat = ier )
   REQUIRE(ier == 0)

   if (maxcyc < 1) then
      ! ------------------------
      ! Single point calculation
      ! ------------------------
      call sqm_energy(natom, x, escf, born_radii, one_born_radii, &
                 intdiel, extdiel, Arad, qm2_struct%scf_mchg )
   else
      ! ---------------------
      ! Geometry optimization
      ! ---------------------
      call xmin(natom, x, f, escf, xmin_iter, maxcyc, born_radii, &
           one_born_radii, intdiel, extdiel, Arad, qm2_struct%scf_mchg, grms_tol, ntpr)
   end if

   ! --------------------
   ! print MO eigenvalues
   ! --------------------
   if ( (qmmm_nml%print_eigenvalues > 0) .and. qmmm_mpi%commqmmm_master) then
      write (6,*) ''
      if (qmmm_nml%qmtheory%DFTB) then
         call print(' Final MO eigenvalues (au)', ks_struct%ev(1:ks_struct%ind(qmmm_struct%nquant_nlink+1)))
      else
         call print(' Final MO eigenvalues (eV)', qm2_struct%eigen_values)
      end if
   end if

   ! ----------------
   ! print SCF energy
   ! ----------------
   string = 'CC triple bond correction (unpublished)'
   call print(cct, string)
   string = 'Nitrogen pyramidalization correction'
   call print(nsp2, string)
   ! at present sqm does only pure QM, need to update this for QM/MM
   total_energy = qmmm_struct%elec_eng +  qmmm_struct%enuclr_qmqm + qmmm_struct%enuclr_qmmm
   write(6,'(/,a,f20.8,a,f18.8,a)') ' Heat of formation   =', &
        escf, ' kcal/mol  (', escf*KCAL_TO_EV, ' eV)'
   write(6,'(/a,f20.8,a,f18.8,a)')   ' Total SCF energy    =', &
        total_energy*EV_TO_KCAL, ' kcal/mol  (', total_energy, ' eV)'
   write(6,'(a,f20.8,a,f18.8,a)')   ' Electronic energy   =', &
        qmmm_struct%elec_eng*EV_TO_KCAL, ' kcal/mol  (', qmmm_struct%elec_eng, ' eV)'
   write(6,'(a,f20.8,a,f18.8,a)')   ' Core-core repulsion =', &
        (qmmm_struct%enuclr_qmqm+qmmm_struct%enuclr_qmmm)*EV_TO_KCAL, ' kcal/mol  (', &
        (qmmm_struct%enuclr_qmqm+qmmm_struct%enuclr_qmmm), ' eV)'
    if (qmmm_nml%qmtheory%DISPERSION .or. qmmm_nml%qmtheory%DISPERSION_HYDROGENPLUS) then
       write(6,'(/a,f20.8,a,f18.8,a)')  ' Dispersion energy   =', &
            qmmm_struct%dCorrection, ' kcal/mol  (', qmmm_struct%dCorrection*KCAL_TO_EV, ' eV)'
    end if
    if (qmmm_nml%qmtheory%DISPERSION_HYDROGENPLUS) then
       write(6,'(a,f20.8,a,f18.8,a)')   ' H-bond energy       =', &
            qmmm_struct%hCorrection, ' kcal/mol  (', qmmm_struct%hCorrection*KCAL_TO_EV, ' eV)'
    end if

   write(6,*) ''
   call qm2_print_charges(1,qmmm_nml%dftb_chg,qmmm_struct%nquant_nlink, &
                            qm2_struct%scf_mchg,qmmm_struct%iqm_atomic_numbers)

   write(6,*) ''
   
   call qm2_calc_dipole(x)
   
   write(6,*) ''
   
   write(6,*) 'Final Structure'
   call qm_print_coords(0,.true.)
   if ( qmmm_nml%printbondorders ) then
      write(6,*) ''
      write(6,*) 'Bond Orders'
      call qm2_print_bondorders()
   end if

   if (qmmm_nml%verbosity > 3) then
      ! Calculate and print also forces in final step
      call sqm_forces(natom, f)
   end if

   write(6,*)
   
   write(6,*) '          --------- Calculation Completed ----------'
   write(6,*)

   call deallocate_qmmm(qmmm_nml, qmmm_struct, qmmm_vsolv, qm2_params)

   call mexit(6,0)

end program sqm

subroutine sqm_energy(natom,coords,escf, &
                 born_radii,one_born_radii, &
                 intdiel, extdiel, Arad, scf_mchg ) 
!
!     Argument list variables:
!
!     coords(natom*3)                 - Cartesian coordinates for all atoms.
!                                       Amber array
!     natom                           - Total number of REAL atoms.
!     qmmm_struct%nquant              - Number of REAL quantum atoms as specified in mdin.
!     iqmatoms(qmmm_struct%nquant)
!                                     - Atom numbers for quantum atoms link atoms given values of -1
!     qmmm_struct%iqm_atomic_numbers(qmmm_struct%nquant) - Atomic numbers for qm atoms.
!     qmmm_struct%nlink               - Number of link atoms.
!     escf                            - Heat of formation from QM.
!     qmmm_struct%qm_coords(3,qmmm_struct%nquant+qmmm_struct%nlink)  
!                                     - Cartesian coordinates of quantum atoms.
!                                       (Extracted from coords by qm_extract_coords)

!    Locally defined arrays:
!    born_radii(1->natom)      - Effective GB radii - only used when doing qm with gb (and qm_gb==2)
!                                Calculated via an initial call to egb.
!    one_born_radii(1->natom)  - 1.0d0/born_radii(i)
!    scf_mchg                  - nquant long, gets filled with the mulliken charges during scf.

   use qmmm_module, only : qmmm_nml,qmmm_struct, qm_gb, qmmm_mpi
   use constants, only : EV_TO_KCAL, KCAL_TO_EV, zero, one, alpb_alpha
  
   implicit none

#include "../include/assert.fh"

#ifdef MPI
   include 'mpif.h'
#endif

! Passed in
   integer, intent(in) :: natom
   _REAL_ , intent(inout)  :: coords(natom*3) !Amber array - adjusted for link atoms
   _REAL_ , intent(out) :: escf
   _REAL_ , intent(in) :: born_radii(natom), one_born_radii(natom)
   _REAL_ , intent(in) :: intdiel, extdiel, Arad
   _REAL_ , intent(inout) :: scf_mchg(qmmm_struct%nquant_nlink)

!Locals
   _REAL_ :: alpb_beta

   integer :: ier=0
   integer i, i3

  ! interface
  !   subroutine qm2_calc_dipole(coord,mass,ipres,lbres,nres)
  !   integer, optional, intent(in) :: nres, ipres(*)
  !   _REAL_, optional, intent(inout) :: mass(*)
  !   _REAL_, intent(inout) :: coord(*)
  !   character(len=4), optional, intent(in) :: lbres(*)
  !   end subroutine qm2_calc_dipole
  ! end interface

!=============================================================================
!                   START OF QMMM SETUP: allocate list memory
!=============================================================================

!  If this is the first call to the routine, do some initial allocation
!  that has not been done elsewhere.
   if (qmmm_struct%qm_mm_first_call) then

     allocate ( qmmm_struct%qm_coords(3,qmmm_struct%nquant_nlink), stat=ier )
                !Stores the REAL and link atom qm coordinates
     REQUIRE(ier == 0)

     !Allocation for QM_GB (qmgb==2)
     if (qmmm_nml%qmgb == 2) then
       !Calculate dielectric factor
       if (qm_gb%alpb_on) then
         alpb_beta=alpb_alpha*(intdiel/extdiel)
         qm_gb%intdieli = one/(intdiel*(one + alpb_beta))
         qm_gb%extdieli = one/(extdiel*(one + alpb_beta))
         qm_gb%one_Arad_beta = alpb_beta/Arad
       else
         qm_gb%intdieli = 1.0d0/intdiel
         qm_gb%extdieli = 1.0d0/extdiel
       end if
       qm_gb%mmcut2 = 999.d0
     end if
   end if ! ---- first call endif ----------

   ! call qm_extract_coords(coords)
   i3 = 0
   !do i=1,natom
   do i=1,qmmm_struct%nquant_nlink
      qmmm_struct%qm_coords(1,i) = coords(i3+1)
      qmmm_struct%qm_coords(2,i) = coords(i3+2)
      qmmm_struct%qm_coords(3,i) = coords(i3+3)
      i3 = i3 + 3
   end do

!=============================================================================
!                   START OF REST OF QMMM SETUP
!=============================================================================
   if(qmmm_struct%qm_mm_first_call) then 
       if (qmmm_mpi%commqmmm_master) then
          write(6,'(/80(1H-)/''  QM CALCULATION INFO'',/80(1H-))')
       end if

       call qm2_load_params_and_allocate(.false.) !Load the parameters
             !Also does a lot of memory allocation and pre-calculates all
             !the STO-6G orbital expansions.

       if (qmmm_mpi%commqmmm_master) then
          ! call qm_print_dyn_mem(natom,qmmm_struct%qm_mm_pairs)
          call qm_print_coords(0,.true.)
          !Finally print the result header that was skipped in sander.
          write(6,'(/80(1H-)/''  RESULTS'',/80(1H-)/)')
       end if
   end if !if (qmmm_struct%qm_mm_first_call)

!======================END OF QMMM SETUP ======================================

   !Calculate RIJ and many related equations here. Necessary memory allocation
   !is done inside the routine.
!Parallel
   call qm2_calc_rij_and_eqns(qmmm_struct%qm_coords, qmmm_struct%nquant_nlink, &
          qmmm_struct%qm_xcrd, natom, qmmm_struct%qm_mm_pairs)
                                !and store them in memory to save time later.

   !============================
   ! Calculate SCF Energy
   !============================
   call qm2_energy(escf, scf_mchg, natom, born_radii, one_born_radii)

   !=============================
   !   Print Mulliken Charges
   !=============================

   if (qmmm_nml%printcharges .and. qmmm_mpi%commqmmm_master) then
     call qm2_print_charges(1,qmmm_nml%dftb_chg,qmmm_struct%nquant_nlink, &
                            scf_mchg,qmmm_struct%iqm_atomic_numbers)
   end if
   !=============================
   !   Print Dipole Charges
   !=============================

   select case (qmmm_nml%printdipole)
      case (1)
         call qm2_calc_dipole(coords)
      case (2)
         write (6,'("QMMM: Not MM part; please check your selection")')
      case default
   end select

   ! Print some extra informatiom about energy contributions
   ! (This is really only required in sander since we print the energies anyways
   !  but kept here for historical reasons)

   if (qmmm_mpi%commqmmm_master) then
      call qm2_print_energy(qmmm_nml%verbosity, qmmm_nml%qmtheory, escf, qmmm_struct)
   end if

   qmmm_struct%qm_mm_first_call = .false.

end subroutine sqm_energy

subroutine sqm_forces(natom, forces)

   !=============================
   ! Calculation of QM Forces
   !=============================

   use qmmm_module, only : qmmm_nml, qmmm_struct, qmmm_mpi
   use constants, only : zero
  
  implicit none

  integer, intent(in) :: natom
  _REAL_ , intent(out) :: forces(natom*3)

  integer :: i, j, m

  qmmm_struct%dxyzqm=zero
   if (qmmm_nml%qmtheory%DFTB) then
     call qm2_dftb_get_qm_forces(qmmm_struct%dxyzqm)
   else
     !standard semi-empirical
     call qm2_get_qm_forces(qmmm_struct%dxyzqm)
   end if

   !NOW PUT THE CALCULATED gradient (not force!) INTO THE SANDER FORCE ARRAY
   do i=1,qmmm_struct%nquant
     m = qmmm_struct%iqmatoms(i)
     m = (m-1)*3
     forces(m+1) = qmmm_struct%dxyzqm(1,i)
     forces(m+2) = qmmm_struct%dxyzqm(2,i)
     forces(m+3) = qmmm_struct%dxyzqm(3,i)
   enddo

   if (qmmm_mpi%commqmmm_master .AND. qmmm_nml%verbosity > 3) then
      
      !If verbosity level is greater than 3 we also print the force array on the QM atoms
      write (6,'("QMMM:")')
      write (6,'("QMMM: Forces on QM atoms from SCF calculation")')
      write (6,'("QMMM: Atm ",i6,": ",3f20.14)') (j,qmmm_struct%dxyzqm(1,j), qmmm_struct%dxyzqm(2,j), &
           qmmm_struct%dxyzqm(3,j), j=1,qmmm_struct%nquant_nlink)
      if (qmmm_nml%verbosity > 4) then
         !Also print info in KJ/mol
         write (6,'("QMMM:")')
         write (6,'("QMMM: Forces on QM atoms from SCF calculation (KJ/mol)")')
         write (6,'("QMMM: Atm ",i6,": ",3f20.14)') (j,qmmm_struct%dxyzqm(1,j)*4.184d0, &
              qmmm_struct%dxyzqm(2,j)*4.184d0, qmmm_struct%dxyzqm(3,j)*4.184d0, &
              j=1,qmmm_struct%nquant_nlink)
      end if
   end if

end subroutine sqm_forces

!======================END OF QM_MM ======================================

!-------------------------------------------------
!     --- FLOAT_LEGAL_RANGE ---
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Check the range of a float; abort on illegal values.
subroutine float_legal_range(string,param,lo,hi)
   implicit none
   _REAL_ param,lo,hi
   character(len=*)string

   if ( param < lo .or. param > hi )then
      write(6,59)
      write(6,60)string,param
      write(6,61)
      write(6,62)lo,hi
      write(6,63)
      call mexit(6,1)
   end if
   59 format(/,1x,'Ewald PARAMETER RANGE CHECKING: ')
   60 format(1x,'parameter ',a,' has value ',e12.5)
   61 format(1x,'This is outside the legal range')
   62 format(1x,'Lower limit: ',e12.5,' Upper limit: ',e12.5)
   63 format(1x,'Check ew_legal.h')
   return
end subroutine float_legal_range 

!-------------------------------------------------
!     --- INT_LEGAL_RANGE ---
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Check the range of an integer; abort on illegal values.
subroutine int_legal_range(string,param,lo,hi)
   implicit none
   integer param,lo,hi
   character(len=*)string

   if ( param < lo .or. param > hi )then
      write(6,59)
      write(6,60)string,param
      write(6,61)
      write(6,62)lo,hi
      write(6,63)
      call mexit(6,1)
   end if
   59 format(/,1x,'PARAMETER RANGE CHECKING: ')
   60 format(1x,'parameter ',a,' has value ',i8)
   61 format(1x,'This is outside the legal range')
   62 format(1x,'Lower limit: ',i8,' Upper limit: ',i8)
   63 format(1x,'The limits may be adjustable; search in the .h files ')
   return
end subroutine int_legal_range 

!-------------------------------------------------
!     --- SANDER_BOMB ---
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Print an error message and quit
subroutine sander_bomb(routine,string1,string2)
   implicit none
   character(len=*) routine,string1,string2

   write(6, '(1x,2a)') &
         'SANDER BOMB in subroutine ', routine
   write(6, '(1x,a)') string1
   write(6, '(1x,a)') string2
   call mexit(6,1)
end subroutine sander_bomb
!-------------------------------------------------

subroutine getsqmx(natom,x,atnam,atnum,ncharge,excharge,chgnam,chgatnum)
   
   !     --- reads initial coords,

   implicit none
   _REAL_ x(*)
   integer i,i3,lun
   integer natom,atnum(*)
   character(len=8) atnam(*)
   character(len=80) line
   ! test-local
   _REAL_ excharge(*)
   integer chgatnum(*)
   character(len=8) chgnam(*)
   integer ncharge
   integer ia, ic, ihead, iend
   logical mdin_qm_atom
   logical mdin_external_charge

   lun = 5 
   mdin_qm_atom = .false.
   mdin_external_charge = .false.
   ncharge = 0

   ! check header names
   ihead=0
   iend=0
   do i=1,9999
      read(lun,'(a)',end=10) line
      if (line(1:1) == "#") then
         if (line(1:10) == "#EXCHARGES") then
            mdin_external_charge = .true.
            ihead = ihead + 1
            ! write(0,*) 'Header "#EXCHARGES" found'
         else if (line(1:4) == "#END") then
            iend = iend + 1
            ! write(0,*) 'Header "#END" found'
         else
            write(0,*) 'Unrecognized header name ', line(1:10)
            call mexit(6,1)
         end if
      end if
   end do

   10 if (iend < ihead) then
      write(0,*) 'Missing "#END" termination sign, exit program'
      call mexit(6,1)
   end if

   rewind(lun)

   !  skip over the &qmmm namelist at the beginning:
   do i=1,20
      read(5,'(a)') line
      if( line(1:2) == " /" ) go to 11
   end do
   write(0,*) 'Error in finding end of qmmm namelist'
   call mexit(6,1)
   
   ! reading QM atoms
   11 i3=0
   ia=0
   do i=1,999
      read(lun,'(a)',end=12) line
      if (line(1:80) /= "") then
         if (line(1:80) /= "#EXCHARGES") then
            ia = ia + 1
            read(line,*,err=15) atnum(ia),atnam(ia),x(i3+1),x(i3+2),x(i3+3)
            i3 = i3 + 3
         else
            go to 12
         end if
      end if
   end do

   12 natom = ia
   !write(0,*) 'finish reading QM atoms, natom =', natom

   ! reading external charges
   if (mdin_external_charge) then
      i3=0
      ic=0
      do i=1,999
         read(lun,'(a)',end=14) line
         if (line(1:80) /= "") then
            if (line(1:80) /= "#END") then
               ic = ic + 1
               read(line,*,err=16) chgatnum(ic), chgnam(ic),excharge(i3+1:i3+4) 
               i3 = i3 + 4
            else
               go to 13
            end if
         end if
      end do
   13 ncharge = ic

      write(6,'(/80(1H-)/''  EXTERNAL CHARGES FOUND IN INPUT'',/80(1H-))')
      write(6,'(2x,"QMMM: External Charge Info")')
      write(6,'(2x,"QMMM:",1x,"ATOMIC",3x,"NAME",8x,"X",9x,"Y",9X,"Z",8X,"CHARGE")')

      i3=0
      do i=1,ncharge
         write(6,'(2x,"QMMM:",3x,i2,6x,a6,4f10.4)') chgatnum(i), chgnam(i), excharge(i3+1:i3+4)
         i3=i3+4
      end do
   end if

   return

   14 write(0,*) 'The termination sign "#END" is missing'
   call mexit(6,1)

   15 write(0,*) 'Error in reading QM atoms'
   call mexit(6,1)

   16 write(0,*) 'Error in reading external charges'
   call mexit(6,1)

end subroutine getsqmx 

!  following stub routine to avoid changing qmmm_module.f for sqm:
#ifdef MPI
subroutine qmmm_vsolv_mpi_setup
   return
end subroutine qmmm_vsolv_mpi_setup
#endif
