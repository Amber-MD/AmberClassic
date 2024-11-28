#include "../include/dprec.fh"
#include "../include/assert.fh"
#include "nfe-config.h"

!------------------------------------------------------------------------------
! sander: the Molecular Dynamics/NMR Refinement/Modeling Module of AMBER
!         (Assisted Model Building with Energy Refinement).  This main driver
!         routine will manage all operations of the program.
!------------------------------------------------------------------------------
subroutine sander()

  use state
  use runmd_module, only : runmd
#if !defined(DISABLE_NFE)
  use nfe_sander_hooks, only : &
      nfe_on_sander_init => on_sander_init, &
      nfe_on_sander_exit => on_sander_exit
  use nfe_sander_proxy, only : infe
#endif /* DISABLE_NFE */

  use lmod_driver
  use constants, only : INV_AMBER_ELECTROSTATIC

#ifdef LES
  use genbornles
#  ifdef MPI
  use les_data, only : cnum
#  endif /* MPI */
#else
  use genborn
#endif /* LES */
  use les_data, only : temp0les
  use fastwt
  use relax_mat
  use nmr, only: nmrrad, impnum
  use ew_recip, only: deallocate_m1m2m3,first_pme
  use parms
  use molecule, only : mol_info, &
                       allocate_molecule, deallocate_molecule
  use nblist, only: cutoffnb, skinnb, nblist_allocate, nblist_deallocate, &
                    nblist_allreal, nblist_allint, num_calls_nblist, &
                    first_list_flag
  use stack

  use sander_rism_interface, only: rism_setparam, rism_init, rism_finalize

  use xray_cpu_module, only: xray_init=>init, xray_read_parm, &
           xray_read_mdin, xray_write_options, xray_write_pdb, xray_write_fmtz
  use xray_globals_module, only: xray_active,pdb_read_coordinates,pdb_outfile, &
         fmtz_outfile

#ifdef MPI /* SOFT CORE */
  use softcore, only: setup_sc, cleanup_sc, ifsc, extra_atoms, sc_sync_x, &
                      summarize_ti_changes, sc_check_perturbed_molecules, &
                      ti_check_neutral, tishake
  use mbar, only: setup_mbar, cleanup_mbar, ifmbar
#endif

  ! For Linear Interaction Energy calculations
  use linear_response, only: ilrt, setup_linear_response, &
                              cleanup_linear_response

#if defined(MPI)
  ! Replica Exchange Molecular Dynamics
  use remd, only : rem, mdloop, remd1d_setup, remd_exchange, reservoir_remd_exchange, &
    remd_cleanup, hremd_exchange, &
    multid_remd_setup, multid_remd_exchange, setup_pv_correction, rremd, stagid
#ifdef VERBOSE_REMD
  use remd, only : repnum
#endif
  use bintraj, only: setup_remd_indices
#else
#  define rem 0
#endif /* MPI */

  use trajenemod, only: trajene

  ! CHARMM support
  use charmm_mod, only : charmm_active, charmm_deallocate_arrays
  use ff11_mod, only : cmap_active, deallocate_cmap_arrays
  use memory_module

  ! Self-Guided molecular/Langevin Dynamics (SGLD)
  use sgld, only: isgld, psgld

  use emap,only: temap,pemap,qemap

  use file_io_dat
  use barostats, only: mcbar_setup
  use random, only: amrset

  use music_module, only: read_music_nml, print_music_settings

  use commandline_module, only: cpein_specified

#ifdef MPI
   use mpi
#endif
  implicit none

  logical belly, erstop
  integer ier,ncalls,xmin_iter,ntbond
  logical ok
  logical newstyle
#  include "nmr.h"
#  include "box.h"
#  include "../include/md.h"
#  include "extra.h"
#  include "tgtmd.h"
#  include "multitmd.h"

#  include "parallel.h"

  ! AMBER/MPI
#ifdef MPI
#  ifdef MPI_DOUBLE_PRECISION
#    undef MPI_DOUBLE_PRECISION
#  endif
!  REMD: loop is the current exchange. runmd is called numexchg times.
  integer loop

  _REAL_ emtmd
#endif /* MPI */
  ! End of declarations and inclues for AMBER/MPI

#  include "ew_pme_recip.h"
#  include "ew_frc.h"
#  include "ew_erfc_spline.h"
#ifdef MPI
#  include "ew_parallel.h"
#endif
#  include "ew_mpole.h"
#  include "ew_cntrl.h"
#  include "def_time.h"
   character(len=30) omp_num_threads

  type(state_rec) ::  ene
  integer native,nr3,nr

  ! nmrcal vars
  _REAL_ enmr(6), devdis(4), devang(4), devtor(4), &
         devplpt(4), devpln(4), devgendis(4)
  _REAL_ ag(1), bg(1), cg(1)

  integer numphi, nhb

  ! trajene var
  _REAL_ carrms

  character(len=8) initial_date, setup_end_date, final_date
  character(len=10) initial_time, setup_end_time, final_time

  ! Needed for final time printout
  integer nstlim_total
  integer, dimension(:), allocatable :: dummy

  _REAL_ time0, time1

  integer i

  integer, dimension(:), allocatable :: ipairs
  logical qsetup
#ifdef MPI
  integer :: n_force_calls
  logical :: do_list_update=.false.
#endif
  _REAL_ :: box_center(3)

  ! Here begin the executable statements.  Start by initializing the cpu
  ! timer. Needed for machines where returned cpu times are relative.
  call date_and_time( initial_date, initial_time )
  call wallclock( time0 )
  call init_timers()

  ! Initialize the printing of ongoing time and performance summaries.
  call print_ongoing_time_summary(0,0,0.0d0,7)

  ! Initialize the number of copies -- always assume 1 until we know otherwise
  ncopy = 1

  ! Flag to tell list builder to print size of list on first call
  first_list_flag = .true.

  ! Flag to tell recip space routines to allocate on first call
  first_pme = .true.

#ifdef MPI
  ! Parallel initialization (setup is done in multisander.F90).
  ! Make PE 0 the master.
  master = (mytaskid == 0)
  master_master = (masterrank == 0)
  if (master .and. numtasks > MPI_MAX_PROCESSORS) then
    write(0, '(a,i4,a,i4)') &
         'Error: the number of processors must not be greater than ', &
         MPI_MAX_PROCESSORS, ', but is ', numtasks
    call mexit(6,1)
  end if
#else /* not MPI follows */
  ! In the single-threaded version, the one process is master
  master = .true.
#endif /* MPI */
  temp0les = 0.d0
  erstop = .false.
  qsetup = .true.

  ! Generic packing scheme
  nwdvar = 1
  native = 32
  numpk = nwdvar
  nbit = native/numpk

  ! Only the master node (only node when single-process)
  ! performs the initial setup and reading/writing
  call timer_start(TIME_TOTAL)

    masterwork: if (master) then

        ! First, initial reads to determine memory sizes
        call mdread1()
        call amopen(8,parm,'O','F','R')
        call rdparm1(8)
        if (mtmd /= 'mtmd' .or. itgtmd == 2) then
          call mtmdlx(natom)
        end if

        ! Now, we can allocate memory
        call locmem()
        write(6,'(/,a,5x,a)') '|','Memory Use     Allocated'
        write(6,'(a,5x,a,i14)') '|', 'Real      ', lastr
        write(6,'(a,5x,a,i14)') '|', 'Hollerith ', lasth
        write(6,'(a,5x,a,i14)') '|', 'Integer   ', lasti
        write(6,'(a,5x,a,i14)') '|', 'Max Pairs ', lastpr

        ! Dynamic memory allocation: allocate space for
        ! module molecule in the master node
        mol_info%natom = natom
        mol_info%nres  = nres
        call allocate_molecule()

        ! Allocate all global arrays
        allocate( x(lastr), ix(lasti), ipairs(lastpr), ih(lasth), stat=ier)
        REQUIRE(ier == 0)
        ix(1:lasti) = 0

        ! This sets up pointer arrays in MEMORY_MODULE to match array-offsets
        ! into the shared X, IX, and IH arrays. Eventually, LOCMEM code should
        ! be merged with MEMORY_MODULE to allocate individual allocatable
        ! arrays, but that will also require updating the MPI code to handle
        ! individual arrays.
        call memory_init()

        ! Allocate the parm arrays
        call allocate_parms()

        if (igb .ne. 0 .and. igb .ne. 10 .and. ipb == 0)  then
          call allocate_gb( natom, ncopy )
        end if

        write(6,'(a,5x,a,i14)'  ) '|', 'nblistReal', nblist_allreal
        write(6,'(a,5x,a,i14)'  ) '|', 'nblist Int', nblist_allint
        write(6,'(a,5x,a,i14,a)') '|', '  Total   ', &
              (8*(lastr+lastrst+nblist_allreal)  &
              + 4*(lasth+lasti+lastpr+lastist+nblist_allint))/1024, &
              ' kbytes'

        ! Finish reading the prmtop file and other user input:
        call rdparm2(x, ix, ih, 8)

      if( xray_active ) then
        call xray_read_mdin(mdin_lun=5)
        call amopen(8,parm,'O','F','R')
        call xray_read_parm(8,6)
        close(8)
        call xray_write_options()
      end if

      call mdread2(x, ix, ih)
      call read_music_nml()
      call print_music_settings()
#ifdef MPI
      if (ifmbar .ne. 0) then
        call setup_mbar(nstlim)
      end if
#endif

!$    call set_omp_num_threads()
!$    call set_omp_num_threads_rism()
      call rism_setparam(mdin, commsander, natom, ntypes, x(L15:L15+natom-1), &
                         x(LMASS:LMASS+natom-1), cn1, cn2, &
                         ix(i04:i04+ntypes**2-1), ix(i06:i06+natom-1))

      ! Evaluate constants frommderead settings
      nr = nrp
      nr3 = 3*nr
      belly = (ibelly > 0)

      ! Seed the random number generator (master node only)
#ifdef MPI
        if (rem == 0) then
          call amrset(ig)
        else

          ! Set the random seed different for different replicas,
          ! but keep same seed for CPUs in the same replica since
          ! we want data from diff numbers of cpus to match.
          !
          ! The variable nodeid is declared through md.h
          ! and is equal to repnum-1
          call amrset(ig + (17 * nodeid))
        end if
#else
        call amrset(ig)
#endif /* MPI */
        if (ntp > 0.and.iabs(ntb) /= 2) then
          write(6,*) 'Input of NTP/NTB inconsistent'
          call mexit(6, 1)
        end if

      ! Read coordinates and velocities.  This begins the qmstep = 1 block
        call timer_start(TIME_RDCRD)
#ifdef LES
#  ifdef MPI
        call getcor(nr, x(lcrd), x(lvel), x(lforce), ntx, box, irest, t, &
                    temp0les, .TRUE., solvph, solve)
#  else
        call getcor(nr, x(lcrd), x(lvel), x(lforce), ntx, box, irest, t, &
                    .TRUE.)
#  endif
#else
#  ifdef MPI
        call getcor(nr, x(lcrd), x(lvel), x(lforce), ntx, box, irest, &
                      t, temp0, .TRUE., solvph, solve)
#  else
        call getcor(nr, x(lcrd), x(lvel), x(lforce), ntx, box, irest, t, &
                      .TRUE.)
#  endif
#endif /* LES */

        ! Set the initial velocities
        if (ntx <= 3) then
          call setvel(nr, x(lvel), x(lwinv), tempi, iscale, scalm)

          ! Random numbers may have been "used up" in setting the intial
          ! velocities; re-set the generator so that all nodes are back in
          ! sync.  This again happens only on the master node.
#ifdef MPI
          if (rem == 0) then
            call amrset(ig)
          end if
#else
          call amrset(ig)
#endif
        else if (irest>0) then
          ! still need to set initial velocities for extra variables:
          if( iscale>0 ) x(lvel+3*nr:lvel+3*nr+iscale-1) = 0.d0
        end if
        if (belly) then
          call bellyf(natom, ix(ibellygp), x(lvel))
        end if
        call timer_stop(TIME_RDCRD)

        ! If we are reading NMR restraints/weight changes, read them now:
        if (nmropt >= 1) then
          call nmrcal(x(lcrd), x(lforce), ih(m04), ih(m02), ix(i02), &
                      x(lwinv), enmr, devdis, devang, devtor, devplpt, &
                      devpln, devgendis, temp0, tautp, cut, x(lnmr01), &
                      ix(inmr02), x(l95), 5, 6, rk, tk, pk, cn1, cn2, ag, &
                      bg, cg, numbnd, numang, numphi, nimprp, nhb, natom, &
                      natom, ntypes, nres, rad, wel, radhb, welhb, rwell, &
                      tgtrmsd, temp0les, -1, 'READ')

          ! Updated 9/2007 by Matthew Seetin to enable plane-point and
          ! plane-plane restraints.  Determine how many of the torsional
          ! parameters are impropers
          call impnum(ix(i46), ix(i56), ix(i48), ix(i58), nphih, nphia, &
                      0, nptra, nimprp)
        end if

        ! Set up info related to weight changes for the non-bonds:
        call nmrrad(rad, wel, cn1, cn2, ntypes, 0, 0.0d0)
        call decnvh(asol, bsol, nphb, radhb, welhb)
        if (iredir(4) > 0) then
          call noeread(x, ix, ih)
        end if
        if (iredir(8) > 0) then
          call alignread(natom, x(lcrd))
        end if
        if (iredir(9) > 0) then
          call csaread
        end if

      ! need to delay xray_init call until here, so that xray_read_pdb
      ! can over-write the coordinates that came from getcor() call above:
      if( xray_active) call xray_init()

      ! Call the fastwat subroutine to tag those bonds which are part
      ! of 3-point water molecules. Constraints will be performed for
      ! these waters using a fast analytic routine.
      call timer_start(TIME_FASTWT)
      call fastwat(ih(m04), nres, ix(i02), ih(m02), nbonh, nbona, ix(iibh), &
                   ix(ijbh), ibelly, ix(ibellygp), iwtnm, iowtnm, ihwtnm, &
                   jfastw, ix(iifstwt), ix(iifstwr), ibgwat, ienwat, ibgion, &
                   ienion, iorwat, 6, natom)
      call timer_stop(TIME_FASTWT)
      call getwds(ih(m04), nres, ix(i02), ih(m02), nbonh, nbona, 0, ix(iibh), &
                  ix(ijbh), iwtnm, iowtnm, ihwtnm, jfastw, ix(iicbh), req, &
                  x(lwinv), rbtarg, ibelly, ix(ibellygp), 6, imin)

      ! Open the data dumping files and position them
      ! depending on the type of run:
      call open_dump_files
      call flush(6)

    end if masterwork
    ! End of master process setup

    ! rism initialization
    call rism_init(commsander)

#ifdef MPI
    call mpi_barrier(commsander,ier)

    ! AMBER/MPI
    !
    ! NOTE: in the current AMBER/MPI implementation, two means of
    ! running in parallel within sander are supported. The value
    ! of mpi_orig determines which approach is used.
    ! This is turned on when minimization (imin .ne. 0) is requested,
    ! and is otherwise off.
    !
    ! When running the mpi_orig case, a variable notdone is now
    ! set by the master and determines when to exit the force()
    ! loop.  When the master has finished calling force, the
    ! master changes notdone to 0 and broadcasts the data one more
    ! time to signal end of the loop.  force() is modified so that
    ! in the mpi_orig case, an initial broadcast is done to receive
    ! the value from the master to decide whether to do the work or
    ! simply exit.
    
    ! Set up initial data and send all needed data to other nodes, {{{
    ! now that the master has it
    !
    ! First, broadcast parameters in memory.h, so that all processors
    ! will know how much memory to allocate:
    call mpi_bcast(natom, BC_MEMORY, mpi_integer, 0, commsander, ier)

    ! Set up integer stack initial size
    call mpi_bcast(lastist, 1, mpi_integer, 0, commsander, ier)
    call mpi_bcast(lastrst, 1, mpi_integer, 0, commsander, ier)

    call stack_setup()
    call mpi_bcast(plumed, 1, MPI_INTEGER, 0, commsander, ier)
    call mpi_bcast(plumedfile, MAX_FN_LEN, MPI_CHARACTER, 0, commsander, ier)

    ! GMS: Broadcast parameters from module 'molecule'
    call mpi_bcast(mol_info%natom, 1, mpi_integer, 0, commsander, ier)
    call mpi_bcast(mol_info%nres, 1, mpi_integer, 0, commsander, ier)
    call mpi_barrier(commsander, ier)
    ! }}}

    ! Allocate memory on the non-master nodes {{{
    if (.not. master) then

      allocate(x(1:lastr), stat=ier)
      REQUIRE(ier == 0)

      allocate(ix(1:lasti), stat=ier)
      REQUIRE(ier == 0)

      allocate(ipairs(1:lastpr), stat=ier)
      REQUIRE(ier == 0)

      allocate(ih(1:lasth), stat=ier)
      REQUIRE(ier == 0)

      ! AWG Set up pointer arrays also on non-master nodes
      call memory_init()

      ! Allocate space for molecule module arrays in the other nodes
      call allocate_molecule()

    end if
    ! End memory allocation on non-master nodes  }}}

    ! Broadcast arrays from module 'molecule'
    call mpi_bcast(mol_info%natom_res, mol_info%nres, MPI_INTEGER, &
                   0, commsander, ier)
    call mpi_bcast(mol_info%atom_to_resid_map, mol_info%natom, MPI_INTEGER, &
                   0, commsander, ier)
    call mpi_bcast(mol_info%atom_mass, mol_info%natom, MPI_DOUBLE_PRECISION, &
                   0, commsander, ier)
    call startup_groups(ier)
    call startup(x, ix, ih)

    call mpi_bcast(xray_active , 1, MPI_LOGICAL, 0, commworld, ier)
    call mpi_bcast (mdin, MAX_FN_LEN, MPI_CHARACTER, 0, commworld, ier)
    call mpi_bcast (parm, MAX_FN_LEN, MPI_CHARACTER, 0, commworld, ier)
    call mpi_bcast (inpcrd, MAX_FN_LEN, MPI_CHARACTER, 0, commworld, ier)

#  if defined(LES)
    call mpi_bcast (ncopy, 1, MPI_INTEGER, 0, commworld, ier)
    call mpi_bcast (cnum(1:natom), natom, MPI_INTEGER, 0, commworld, ier)
    call mpi_bcast (evbin, MAX_FN_LEN, MPI_CHARACTER, 0, commworld, ier)
#  endif /* LES */

    if (ifsc .ne. 0) then

      ! Multi-CPU minimization does not work with soft core
      if (imin > 0 .and. numtasks > 1) then
        call sander_bomb('imin > 0 and numtasks > 1', &
                     'TI minimizations cannot be performed with > 2 CPUs','')
      end if
      call setup_sc(natom, nres, ih(m04), ih(m06), &
                    ix(i02), ih(m02), x(lcrd), ntypes, clambda, nstlim)
      if (ntp > 0 .and. master) then

        ! Check which molecules are perturbed in NPT runs
        call sc_check_perturbed_molecules(nspm, ix(i70))
      end if

      ! Make sure all common atoms have the same v (that of V0) in TI runs
      if (ifsc .ne. 2) then
        if (master) call sc_sync_x(x(lvel), nr3)
        if (numtasks > 1) then
          call mpi_bcast(nr3, 1, MPI_INTEGER, 0, commsander, ier)
          call mpi_bcast(x(lvel), nr3, MPI_DOUBLE_PRECISION, &
                         0, commsander, ier)
        end if
      end if
      if (tishake .ne. 0) then
        call setnoshake_sc(ix, ntc, num_noshake, master)
      end if
    else
      extra_atoms=0
    end if
    if (.not. master .and. igb == 0 .and. ipb == 0) then
      call nblist_allocate(natom,ntypes,num_direct,numtasks)
    end if

    ! Allocate memory for GB on the non-master nodes:
    if (.not. master) then
      if (igb /= 0 .and. igb /= 10 .and. ipb == 0) then
          call allocate_gb( natom, ncopy )
      end if
    end if
    ! End non-msater process GB memory allocation

    nr = nrp
    nr3 = 3 * nr
    belly = (ibelly > 0)

    ! All nodes are calling this. amrset(ig) has already been called
    ! by the master node (twice if initial velocities are set, ntx <= 3).
    ! Replica Exchange MD now requires need to call amrset on all child
    ! threads, masters have called it above before initial coord read.
      if (rem == 0) then
        call amrset(ig+1)
      else
        if (.not. master) then
          call amrset(ig + 17*nodeid)
        end if
      endif
      if (nmropt >= 1) then
        call nmrcal(x(lcrd), x(lforce), ih(m04), ih(m02), ix(i02), x(lwinv), &
                    enmr, devdis, devang, devtor, devplpt, devpln, devgendis, &
                    temp0, tautp, cut, x(lnmr01), ix(inmr02), x(l95), 5, 6, &
                    rk, tk, pk, cn1, cn2, ag, bg, cg, numbnd, numang, numphi, &
                    nimprp, nhb, natom, natom, ntypes, nres, rad, wel, radhb, &
                    welhb, rwell, tgtrmsd, temp0les, -1, 'MPI ')
      end if
    call mpi_bcast(lmtmd01, 1, mpi_integer, 0, commsander, ier)
    call mpi_bcast(imtmd02, 1, mpi_integer, 0, commsander, ier)
    if (itgtmd == 2) then
      call mtmdcall(emtmd, x(lmtmd01), ix(imtmd02), x(lcrd), x(lforce), &
                    ih(m04), ih(m02), ix(i02), ih(m06), x(lmass), natom, &
                    nres, 'MPI ')
    end if

    ! Check that the system is neutral and print warning message
    ! if not.  Adjust charges for roundoff error.

    !  DAC note: to do mixed GB/RISM calcs, both sides need to call this
    ! if (igb == 0 .and. ipb == 0 .and. iyammp == 0) then
      if (icfe == 0) then
        call check_neutral(x(l15),natom)
      else
        call ti_check_neutral(x(l15),natom)
      end if
    !end if

    ! Use old parallelism for energy minimization
    if (imin .ne. 0) then
      mpi_orig = .true.
      notdone = 1
    else
      mpi_orig = .false.
    end if

    if (mpi_orig .and. .not. master) then
      ! All nodes only do the force calculations.  Minimization  {{{
      ! so only master gets past the loop below
      if (igb == 7 .or. igb == 8) then

        ! x(l97) is rborn(), the Generalized Born radii
        call igb7_init(natom, x(l97))

        ! TODO: add igb == 8 here
      end if
      n_force_calls = 0
      do while(notdone == 1)
        n_force_calls = n_force_calls+1
        call force(x, ix, ih, ipairs, x(lcrd), x(lforce), ene,ene%vir, &
                   x(l96), x(l97), x(l98), x(l99), qsetup, &
                   do_list_update, n_force_calls)
      end do

      ! Deallocate and return
      goto 999
      ! }}}
    end if

    ! Report the parallel configuration {{{
    if (master) then
      write(6, '(a,i4,a,/)') '|  Running AMBER/MPI version on ', &
                               numtasks, ' nodes'

      ! BTREE is selected by default if cpu is a power of two.
      ! The number of processes is required to be a power of two for Btree
      ! Print a warning about inefficiency with cpus not being a power of 2.
      if (numtasks > 1 .and. logtwo(numtasks) <= 0) then
          write(6, '(a,i4,a,/)') &
                '|  WARNING: The number of processors is not a power of 2'
          write(6, '(a,i4,a,/)') &
                '|           this may be inefficient on some systems.'
      end if
    end if
    ! End reporting of parallel configuration }}}

    if (master .and. numgroup > 1) then
      write(6, '(a,i4,a,i4,a,i4,a)') '|  MULTISANDER: ', numgroup, &
           ' groups. ',  numtasks, ' processors out of ', worldsize, &
           ' total.'
    end if

    if (master) then
      call flush(6)
    end if

    ! End of AMBER/MPI work in this block

#else

    ! For debugging, the charges must be copied at the start so that
    ! they can't change later.  Check that the system is neutral and
    ! print warning message adjust charges for roundoff error.
    call check_neutral(x(l15),natom)
    call amrset(ig+1)
    call stack_setup()

#endif /* MPI */

    ! Initialize LIE module if used
    if ( ilrt /= 0 ) then
      call setup_linear_response(natom, nres, ih(m04), ih(m06), ix(i02), &
                                 ih(m02), x(lcrd), x(l15), ntypes, ix(i04), &
                                 ix(i06), cn1, cn2, master)
    end if
    call date_and_time(setup_end_date, setup_end_time)

    ! Initialize the printing of ongoing time and performance summaries.
    ! We do this quite late here to avoid including all the setup time.
    if (master) then
      call print_ongoing_time_summary(0, 0, 0.0d0, 7)
    end if

    ! Now do the dynamics or minimization.

    if (igb == 7 .or. igb == 8 ) then
      call igb7_init(natom, x(l97)) !x(l97) is rborn()
    end if

    ! Use the debugf namelist to activate
    call debug_frc(x, ix, ih, ipairs, x(lcrd), x(lforce), cn1, cn2, qsetup)

    ! Prepare for SGLD simulation
    if (isgld > 0) then
      call psgld(natom,ix(i08), ix(i10), x(lmass),x(lcrd),x(lvel), rem)
    end if

    ! Prepare for EMAP constraints
    if (temap) then
      call pemap(dt,temp0,x,ix,ih)
    end if

    if (master) then
      write(6,'(/80(''-'')/,''   4.  RESULTS'',/80(''-'')/)')
    end if

    ! Set up the MC barostat if requested
    if (ntp > 0 .and. barostat == 2) then
      call mcbar_setup(ig)
    end if

    ! Input flag imin determines the type of calculation: MD, minimization, ...
    select case (imin)
      case (0)

        ! Dynamics:  {{{
        call timer_start(TIME_RUNMD)

#ifdef MPI
        ! Replica Exchange Molecular Dynamics.  If this is not a REMD run,
        ! runmd is called only once.  If this is a REMD run, runmd is
        ! called 0 to numexchg times, where the 0th runmd is just for getting
        ! initial PEs (no dynamics).
        if (rem == 0) then

          ! Not a REMD run. runmd will be called once.
          loop = 0
        else

          ! This is a REMD run. runmd will be called numexchg times.
          loop = numexchg
          if (rem < 0) then

            ! Multi-D REMD run
            call multid_remd_setup(numexchg, numwatkeep, temp0, mxvar, &
                                   natom, ig, solvph, solve, irest)
          else

            ! 1D REMD. Set up temptable, open remlog, etc.
            call remd1d_setup(numexchg, hybridgb, numwatkeep, &
                              temp0, mxvar, natom, ig, solvph, solve,stagid)
          end if
          call setup_pv_correction(6, ntp, sanderrank)

          ! Now set up REMD indices for traj/restart writes. Only do this on
          ! master since only master handles writes.
          if (master) then
            call setup_remd_indices
          end if
        end if ! Replica run setup

        ! Loop over REMD exchanges
        do mdloop = 0, loop

          ! REMD exchange handling.  mdloop == 0 is just used
          ! to get initial energies for the first exchange.
          if (rem < 0 .and. mdloop > 0) then
            call multid_remd_exchange(x, ix, ih, ipairs, qsetup, &
                                      do_list_update, temp0, solvph, solve, &
                                      reservoir_exchange_step)
          else if ((rem == 1 .or. rem == 2) .and. mdloop > 0) then
             if(rremd>0) then
                call reservoir_remd_exchange(1, rem, x(lcrd), x(lvel), x(lmass), &
                               nr3, natom, nr, temp0, reservoir_exchange_step)
             else
                call remd_exchange(1, rem, x(lcrd), x(lvel), x(lmass), &
                               nr3, natom, nr, temp0)
             end if
          else if (rem == 3 .and. mdloop > 0) then
            call hremd_exchange(1, x, ix, ih, ipairs, qsetup, do_list_update)

            ! Force was called inside hremd_exchange, so call nmrdcp
            ! to decrement the NMR counter, since this should not count
            ! as a real step. This is OK, since the counter got
            ! incremented at the _very_ end of nmrcal, so we haven't already
            ! printed an unwanted value (JMS 2/12)
            if (nmropt .ne. 0) then
              call nmrdcp
            end if
          end if
          ! End of block for replica exchange handling

#ifdef VERBOSE_REMD
          if (rem > 0 .and. mdloop .eq. 0 .and. master) then
            write (6,'(a,i4)') "| REMD: Getting initial energy for replica ", &
                               repnum
          end if
#endif /* VERBOSE_REMD */
#endif  /* MPI */

#if !defined(DISABLE_NFE)
#ifdef MPI
          call mpi_bcast(infe, 1, mpi_integer, 0, commsander, ier)
#endif
          if ( infe == 1) &
            call nfe_on_sander_init(ih, x(lmass), x(lcrd), rem)
#endif /* DISABLE_NFE */

          ntbond = nbonh  + nbona + nbper
          call runmd(x, ix, ih, ipairs, coord3, massinv, mass, &
                       force3, vel3, vel3_old, rest_coord3, &
                       conp, skip, atoms_per_molecule, erstop, qsetup)

#if !defined(DISABLE_NFE)
          if (infe == 1) call nfe_on_sander_exit()
#endif /* DISABLE_NFE */
          ! End of Replica Exchange MD block

#ifdef MPI
        end do
        ! End of loop over REMD exchanges

        ! Cleanup REMD files.
        if (rem .ne. 0) then
          call remd_cleanup()
        end if
#endif
        call timer_stop(TIME_RUNMD)
        if (master) then
          call flush(6)
        end if

        ! The erstop error condition stems from subroutine shake;
        ! furthermore, it seems that erstop can never be true since shake
        ! can never return with its third last argument, niter, equal to 0.
        if (erstop) then
          if (master) then
            write(6, *) 'FATAL ERROR'
          end if
          call mexit(6,1)
        end if
      ! }}}
      case (1)

        ! Minimization:  {{{
        !    input flag ntmin determines the method of minimization
        select case (ntmin)

          case (LMOD_NTMIN_XMIN)
            write(6,'(a)') '  LMOD XMIN Minimization.'
            write(6,'(a)') ''
            write(6,'(a)') '  Note: Owing to the behaviour of the XMIN &
                            &algorithm,'
            write(6,'(a)') '        coordinates in the trajectory and &
                            &intermediate'
            write(6,'(a)') '        restart files will not match up with &
                            &energies'
            write(6,'(a)') '        in the mdout and mdinfo files. The final &
                            &energy'
            write(6,'(a)') '        and final coordinates do match.'
            write(6,'(a)') ''
            xmin_iter = 0
            call run_xmin(x, ix, ih, ipairs, x(lcrd), x(lforce), &
                          ene, qsetup, xmin_iter, ntpr, iscale)
            if (master) then

              ! Write the restart file
              call minrit(0,nrp,ntxo,x(lcrd))
            end if
          case (LMOD_NTMIN_LMOD)
            write(6,'(a)') '  LMOD LMOD Minimization.'
            write(6,'(a)') ''
            write(6,'(a)') '  Note: Owing to the behaviour of the XMIN &
                            &algorithm,'
            write(6,'(a)') '        coordinates in the trajectory and &
                            &intermediate'
            write(6,'(a)') '        restart files will not match up with &
                            &energies'
            write(6,'(a)') '        in the mdout and mdinfo files. The final &
                            &energy'
            write(6,'(a)') '        and final coordinates do match.'
            write(6,'(a)') ''
            call run_lmod(x, ix, ih, ipairs, x(lcrd), x(lforce), ene, qsetup)
            if (master) then

              ! Write the restart file
              call minrit(0, nrp, ntxo, x(lcrd))
            end if
          case default

            ! invalid ntmin: input validation occurs in mdread.f
            ASSERT(.false.)
        end select
      ! }}}

      case (5)

        ! Modified for reading trajectories (trajene option) {{{
        write (6,*) "POST-PROCESSING OF TRAJECTORY ENERGIES"

        ! Read trajectories and calculate energies for each frame
        call trajene(x, ix, ih, ipairs, ene, ok, qsetup)
        if (.not. ok) then
          write (6,*) 'error in trajene()'
          call mexit(6, 1)
        end if
      ! }}}

      case default

        ! Invalid imin: input validation should be transferred to mdread.f
        write(6,'(/2x,a,i3,a)') 'Error: Invalid IMIN (', imin, ').'
        ASSERT(.false.)
    end select

#ifdef MPI /* SOFT CORE */
    if (master) then
      if (icfe .ne. 0 .and. ifsc == 1) then
        call summarize_ti_changes(natom, resat)
      end if
    end if
#endif

    ! Finish up EMAP
    if (temap) call qemap()

  ! Calc time spent running vs setup
  call timer_stop(TIME_TOTAL)
  call wallclock( time1 )
  call date_and_time( final_date, final_time )
#ifdef MPI
  call profile_time(time1 - time0, num_calls_nblist, profile_mpi)
#else
  call profile_time(time1 - time0, num_calls_nblist)
#endif

#ifdef MPI
  ! Set and broadcast notdone in mpi_orig case to inform
  ! other nodes that we are finished calling force().
  if (mpi_orig) then
    notdone = 0
    call mpi_bcast(notdone, 1, mpi_integer, 0, commsander, ier)
  end if
#endif

  if (master) then
    call close_dump_files
  end if

    ! Write out the final timings, taking Replica Exchange MD into account
#ifdef MPI
    if (rem .ne. 0) then
      nstlim_total = nstlim * numexchg
    else
      nstlim_total = nstlim
    end if
#else
    nstlim_total = nstlim
#endif
    if( master ) then
       if (imin == 0) then
         call print_ongoing_time_summary(nstlim_total, nstlim_total, dt, 6)
       end if
       write(6,'(12(a))') '|           Job began  at ', initial_time(1:2), &
             ':', initial_time(3:4), ':', initial_time(5:10), '  on ',&
             initial_date(5:6), '/', initial_date(7:8), '/', initial_date(1:4)
       write(6,'(12(a))') '|           Setup done at ', setup_end_time(1:2),  &
             ':', setup_end_time(3:4), ':', setup_end_time(5:10), '  on ', &
             setup_end_date(5:6), '/', setup_end_date(7:8), '/', &
             setup_end_date(1:4)
       write(6,'(12(a))') '|           Run   done at ', final_time(1:2),  &
             ':', final_time(3:4), ':', final_time(5:10), '  on ', &
             final_date(5:6), '/', final_date(7:8), '/', final_date(1:4)
       call nwallclock( ncalls )
       write(6, '(''|'',5x,''wallclock() was called'',I8,'' times'')') ncalls
       call flush(6)
    end if

#ifdef MPI
   ! --- dynamic memory deallocation:
   999 continue
#endif

#ifdef MPI /* SOFT CORE */
  if (ifsc .ne. 0) then
    call cleanup_sc()
  end if
  if (ifmbar .ne. 0) then
    call cleanup_mbar()
  end if
#endif

  ! Finalize Linear Interaction Energy module if initiated above
  if (ilrt .ne. 0) then
    call cleanup_linear_response(master)
  end if

   call rism_finalize()

   if( xray_active .and. master ) then
      if (pdb_outfile /= '') then
         call xray_write_pdb(pdb_outfile)
      end if
      if (fmtz_outfile /= '') then
         call xray_write_fmtz(fmtz_outfile)
      endif
   end if

  call nblist_deallocate()
  call deallocate_stacks()
  if (igb /= 0 .and. igb /= 10 .and. ipb == 0)  then
    call deallocate_gb()
  end if
  deallocate(ih, stat=ier)
  REQUIRE(ier == 0)
  deallocate(ipairs, stat=ier)
  REQUIRE(ier == 0)
  deallocate(ix, stat=ier)
  REQUIRE(ier == 0)
  deallocate(x, stat=ier)
  REQUIRE(ier == 0)
  if (ntb > 0 .and. ifbox == 1) then
    call deallocate_m1m2m3()
  end if
  call deallocate_molecule()

  if (master .and. mdout .ne. 'stdout') then
    close(6)
  end if

  return

end subroutine sander
