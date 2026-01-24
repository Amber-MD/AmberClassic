! <compile=optimized>
#include "../include/assert.fh"
#include "../include/dprec.fh"

module remd

! MODULE: REMD
! ================== REPLICA EXCHANGE MOLECULAR DYNAMICS ====================
! Daniel R. Roe, 2007
! Jason M. Swails, 2013 (multi-D/pH REMD)
! Vinicius Wilian D. Cruzeiro, 2017 (E REMD/added pH and E on multi-D)
! Based on original implementation in multisander.F90 by:
!    Guanglei Cui
!    Carlos Simmerling
!
! ---=== Subroutine list ===---
!   Called from outside remd.F90:
!     remd1d_setup (from sander.F90)
!     multid_remd_setup (from sander.F90)
!     remd_cleanup (from sander.F90)
!     remd_exchange (from sander.F90)
!     hremd_exchange (from sander.F90)
!     reservoir_remd_exchange (from sander.F90)
!     hybrid_remd_ene (from runmd.F90)
!     multid_remd_exchange (from sander.F90)
!     setup_pv_correction (from sander.F90)
!   Internal:
!     subrem
!     subrem_reservoir
!     remd_scale_velo
!     set_partners
!     sorttable
!     load_reservoir_structure
!     stripwat
!     calc_rremd_cluster
!     load_reservoir_files
!     templookup *This is a function.
!     pv_correction
!     remd_bcast_cell

#ifdef MPI

use file_io_dat, only : MAX_FN_LEN, REMLOG_UNIT, REMTYPE_UNIT, &
                        REMSTRIPCOORD_UNIT, REMIN_UNIT, RESERVOIR_UNIT, &
                        REMD_DIM_UNIT, &
                        remlog, remtype, remstripcoord, saveenefile, &
                        clusterinfofile, reservoirname, remd_dimension_file, &
                        inpcrd, numgroup, owrite, target_ph, target_e
use commandline_module, only : cpein_specified
use random
use mpi

implicit none

! DAN ROE: None of this should be available to sander if no MPI - REMD only
!          works in parallel.

#include "extra.h"

! ... MPI Communicators and variables:
!  remd_comm         - Communicator that links all masters in a REMD group
!  group_master_comm - Communicator of all masters of remd_comm
!  group_master_rank - Rank in group_master_comm
!  group_master_size - Size of the group_master_comm
!  num_rem_grp       - Number of REMD groups
!  remd_rank         - Rank in remd_comm
!  remd_size         - Size of remd_comm
!  remd_master       - Am I the master of remd_comm?
!  master_master     - Am I the master of the group_master_comm?
integer, save     :: remd_comm
integer, save     :: group_master_comm
integer, save     :: group_master_rank
integer, save     :: group_master_size
integer, save     :: num_rem_grp
integer, save     :: remd_rank
integer, save     :: remd_size
logical, save     :: remd_master

! ... Constants:
!   maxreservoir - largest number of structures allowed in a reservoir.
!                  Currently this is only necessary since reservoir restart
!                  files are assumed to have a 6 digit extension (see
!                  load_reservoir_structure()).
integer, parameter :: maxreservoir=999999

! ... logical variables:
!   hybridwritetraj - if true a trajectory of stripped coords will be written
!                     during hybrid remd.
!   jumpright       - if true replica that controls the exchange will attempt
!                     to exchange with the replica at next highest temperature.
!                     (one value for each dimension)
!   initmodwt       - If nmropt>0 modwt will reset temp0 to its initial value
!                     after each exchange. initmodwt will force a re-read of
!                     temp0 after each exchange.
!   reserv_velo     - If the reservoir restart files have velocities
!   even_replica    - If this replica has an even index. Helps decide whether or
!                     not this replica will exchange
!   even_exchanges  - Whether even-numbered replicas exchange here or not.
!   exchange_ucell  - If true, exchange unit cells during Hamiltonian exchange.
!   use_pv          - If true, use the pressure/volume correction to exch. calc.
logical, allocatable, save :: jumpright(:)
logical,              save :: hybridwritetraj
logical,              save :: initmodwt
logical,              save :: reserv_velo
logical, allocatable, save :: even_replica(:)
logical, allocatable, save :: even_exchanges(:)
logical,              save :: exchange_ucell = .false.
logical,              save :: use_pv = .false.

! ... integer variables:
! mdloop      - Current exchange index. mdloop will be always 0 unless rem>0.
! rem         - The type of replica exchange:
!              -1, multi-D REMD
!               0, regular MPI/MD
!               1, conventional REM
!               2, partial REM
!               3, h-remd
!               4, pH-REMD
!               5, E-REMD
!               6, Reservoir REMD
! next_rem_method - same as "rem", but is the "rem" for the next dimension we
!                   are exchanging in.
! rremd       - The type of reservoir replica exchange:
!               0, no reservoir
!               1, Boltzmann weighted reservoir
!               2, 1/N weighted reservoir
!               3, reservoir with weights defined by dihedral clustering
! repnum          - Replica number.
! numreps         - Total # of replicas, set equal to numgroup.
! remd_dimension  - The total number of exchange dimensions we have.
! remd_repidx     - Overall replica index for this replica.
! remd_crdidx     - Overall coordinate index for this replica.
! exchsuccess     - The number of successful exchanges between each neighbor
!                   pair for printing acceptance %.
! group_num       - Which group we are a part of in each dimension.
! replica_indexes - A collection of which ranks this replica is in each of its
!                   REMD dimensions.
! partners    - repnums of replica to the left and replica to the right.
! index_list  - Ordered list of repnum for each replica in the active dimension.
!               This is updated in subroutine set_partners(), which should be
!               called at the beginning of each exchange subroutine.
! remd_types  - What type of exchange attempt we make in each dimension:
!               TEMPERATURE : 1
!               HAMILTONIAN : 3
!               PH          : 4
!               REDOX       : 5
!               RESERVOIR   : 6
integer, save :: mdloop
integer, save :: rem
integer, save :: next_rem_method
integer, save :: rremd
integer, save :: repnum
integer, save :: numreps
integer, save :: remd_dimension
integer, save :: remd_repidx
integer, save :: remd_crdidx
integer, allocatable, save  :: exchsuccess(:,:)
integer, allocatable, save  :: group_num(:)
integer, allocatable, save  :: replica_indexes(:)
integer, dimension(2), save :: partners
integer, allocatable, save  :: index_list(:)
integer, allocatable, save  :: remd_types(:)
integer, allocatable, save  :: num_right_exchg(:,:,:)
integer, allocatable, save  :: num_left_exchg(:,:,:)

! Data structure for holding properties and bookkeeping
! Optimized for MPI communications
type remd_data
   sequence
   _REAL_   :: mytemp         ! My current instantaneous temperature
   _REAL_   :: myEptot        ! My current potential energy
   _REAL_   :: mytargettemp   ! My current target temperature (temp0)
   _REAL_   :: myscaling      ! My scaling factor
   _REAL_   :: newtargettemp  ! My new target temperature after exchange attempt
end type remd_data

integer, parameter :: SIZE_REMD_DATA = 5
type (remd_data), parameter :: NULL_REMD_DATA = &
         remd_data(0.d0, 0.d0, 0.d0, -1.d0, -99.d0)

type (remd_data), save      :: my_remd_data

! Data for pressure/volume correction.
! TODO incorporate into remd_data?
_REAL_, save :: my_pressure = 0.d0
_REAL_, save :: my_volume = 0.d0

! Data structure for outputting stuff to the REM.log file

type :: remlog_data
   sequence
   _REAL_   :: scaling        ! T-REMD
   _REAL_   :: real_temp      ! T-REMD
   _REAL_   :: new_temp0      ! T-REMD
   _REAL_   :: struct_num     ! T-REMD
   _REAL_   :: pot_ene_tot    ! T-REMD & H-REMD
   _REAL_   :: temp0          ! T-REMD & H-REMD
   _REAL_   :: success_ratio  ! T-REMD & H-REMD & pH-REMD & E-REMD
   _REAL_   :: num_rep        ! T-REMD & H-REMD & pH-REMD & E-REMD
   _REAL_   :: group_num      ! T-REMD & H-REMD & pH-REMD & E-REMD
   _REAL_   :: repnum         !          H-REMD
   _REAL_   :: neighbor_rep   !          H-REMD
   _REAL_   :: nei_pot_ene    !          H-REMD
   _REAL_   :: left_fe        !          H-REMD
   _REAL_   :: right_fe       !          H-REMD
   _REAL_   :: left_exchg     !          H-REMD
   _REAL_   :: right_exchg    !          H-REMD
   _REAL_   :: success        !          H-REMD & pH-REMD & E-REMD
   _REAL_   :: my_ph          !                   pH-REMD
   _REAL_   :: nei_ph         !                   pH-REMD
   _REAL_   :: nprot          !                   pH-REMD
   _REAL_   :: my_e           !                             E-REMD
   _REAL_   :: nei_e          !                             E-REMD
   _REAL_   :: nelec          !                             E-REMD
end type remlog_data

integer, parameter :: SIZE_REMLOG_DATA = 23

type(remlog_data), allocatable, save :: multid_print_data(:)
type(remlog_data), allocatable, save :: multid_print_data_buf(:,:)
type(remlog_data), parameter         :: NULL_REMLOG_DATA = &
    remlog_data(0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0, &
                0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0, &
                0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0, &
                0.d0,0.d0)

! ... floats:
! remd_ekmh     - Store ekmh variable between runmd calls (exchanges). ekmh is
!                 the KE from old velocities.
! statetable()  - Store sorted list of state values. Only used for RXSGLD
!                 outside of setup routines.
!                 TODO Remove dependence of RXSGLD on this array and treat
!                 indexing the same way (so it will work with multi-D REMD)
_REAL_,              save :: remd_ekmh
_REAL_, allocatable, save :: statetable(:)
_REAL_, allocatable, save :: total_left_fe(:,:,:)
_REAL_, allocatable, save :: total_right_fe(:,:,:)

! ... Reservoir REMD (RREMD):
! saveene()      - Energy of each structure in the reservoir.
! restemp0       - reservoir temperature
! rremd_idx      - # of random structure from reservoir
! reservoirsize  - total # of structures in reservoir
! reservoir_ncid - NCID of netcdf reservoir; -1 when not netcdf
_REAL_, allocatable, save :: saveene(:)
_REAL_,              save :: restemp0
integer,             save :: rremd_idx, reservoirsize
integer,             save :: reservoir_ncid = -1
! Required for netcdf reservoir
integer,             save :: coordVID, velocityVID

! ... RREMD Dihedral Clustering:
! clusternum()   - cluster that reservoir structure belongs to
! clusterid(,)   - 2d array of bin values for each dihedral, 1 set per cluster
! clustersize()  - # structures from reservoir in this cluster (note 0 is unknown
!                  cluster, currently set to 0)
! dihclustat(,)  - 2d, 4 atoms that define each of the dihedrals used for clustering
!     TODO: dihclustnbin is only used to hold values for output - could be axed.
! dihclustnbin() - #bins used for each of the dihedrals (each is dihclustnbin/360
!                  in width)
! dihclustmin()  - For each dihedral, torsion value that should get bin 0
! dihcluststep() - Step size for each dihedral (calculated during read of
!                  clusterinfo as 360 / Nbins).
! currdihid()    - dih bin values for the current MD structure (to assign it to a
!                  cluster)
! incluster      - cluster for the current MD structure
! nclust         - # clusters
! nclustdih      - # dihedrals used for clustering
integer, allocatable, save :: clusternum(:)
integer, allocatable, save :: clusterid(:,:)
integer, allocatable, save :: clustersize(:)
integer, allocatable, save :: currdihid(:)
integer, allocatable, save :: dihclustnbin(:)
integer, allocatable, save :: dihclustat(:,:)
_REAL_,  allocatable, save :: dihclustmin(:)
_REAL_,  allocatable, save :: dihcluststep(:)
integer,              save :: incluster
integer,              save :: nclust
integer,              save :: nclustdih

! Hybrid REMD Temp. Storage
_REAL_, dimension(:), allocatable, private :: hybrid_coord
_REAL_, dimension(:), allocatable, private :: hybrid_force
_REAL_, dimension(:), allocatable, private :: hybrid_refc

! H-remd temporary coordinates, force
_REAL_, dimension(:), allocatable, private :: xtemp, ftemp

! RXSGLD variables
! stagid       - current stage of this replica
! myscalsg     - scaling factor for guiding properties after exchange. (in sgld
!                module)
! sgfttable()  - Store sorted list of replica self-guiding momentum factors.
! sgfftable()   - Store sorted list of replica self-guiding force factors.
! psgldgtable()   - Store sorted list of replica self-guiding effective force factors.
integer, save :: stagid
_REAL_, allocatable, save :: sgfttable(:),sgfftable(:),psgldgtable(:)
! REMD RNG
type(rand_gen_state), save :: remd_gen

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! The following subroutines are used for REMD


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Determine if the pressure-volume correction should be used.
subroutine setup_pv_correction(outu, ntpin, rank)
  ! USE STATEMENTS
  ! ARGUMENTS
  implicit none
  integer, intent(in) :: outu, ntpin, rank

  if (ntpin > 0) then
    exchange_ucell = .true. ! Only matters for hremd_exchange
    use_pv = .true.
    if (rank .eq. 0) &
      write(outu,'(a)') '| REMD: Pressure/volume correction to exchange calc active for TREMD/HREMD.'
  else
    exchange_ucell = .false.
    use_pv = .false.
  endif
end subroutine setup_pv_correction

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Set up REMD run - open log, set up temperature table, etc
subroutine remd1d_setup(numexchg, hybridgb, numwatkeep, temp0, &
                        mxvar, natom, ig, solvph, solve, stagid)

   use sgld, only : isgld, sgfti, sgffi, psgldg, sorttempsg, tempsglookup

implicit none
#  include "parallel.h"

! Passed variables
   integer, intent(in) :: numexchg
   integer, intent(in) :: hybridgb
   integer, intent(in) :: numwatkeep
   integer, intent(in) :: mxvar
   integer, intent(in) :: natom
   integer, intent(in) :: ig
   integer, intent(in out) :: stagid

   _REAL_, intent(in) :: temp0
   _REAL_, intent(in) :: solvph
   _REAL_, intent(in) :: solve

! Local variables
   integer i, j, ierror, add_fac

!--------------------

   ! ---=== INITIALIZE VARIABLES ===---

   ! Allocate jumpright (only 1 dimension)
   allocate(jumpright(1), stat=ierror)
   REQUIRE(ierror == 0)
   jumpright(1) = .true.
   remd_master = .false.

   my_remd_data = NULL_REMD_DATA
   initmodwt = .true.
   rremd_idx = -1
   remd_ekmh = 0.0d0

   ! Nullify communicators
   remd_comm = mpi_comm_null
   group_master_comm = mpi_comm_null

   ! next_rem_method is always rem
   next_rem_method = rem

   ! This is a 1-D REMD setup. Set up variables for this
   remd_dimension = 1
   if (master) then
      call mpi_comm_dup(commmaster, remd_comm, ierror)
      call mpi_comm_rank(remd_comm, remd_rank, ierror)
      call mpi_comm_size(remd_comm, remd_size, ierror)
      remd_master = remd_rank == 0
      group_master_rank = 0
      group_master_size = 1
      num_rem_grp = 1
      master_master = remd_master
   end if

   ! Initialize the RNG (nodeid == repnum - 1).
   add_fac = mod(repnum-1, numreps) / 2 * 2 ! For compatibility with pmemd
   call amrset_gen(remd_gen, ig + add_fac)

   ! Allocate memory for statetable and exchsuccess
   allocate(exchsuccess(1, numreps), &
            replica_indexes(1),      &
            remd_types(1),            &
            group_num(1),            &
            index_list(numreps),     &
            stat=ierror)
   REQUIRE(ierror==0)
   group_num(1) = 1

   ! Setting remd_types
   remd_types(1) = rem

   ! allocate memory for statetable

   if (rem /= 3) then
      allocate(statetable(numreps), stat=ierror)
      REQUIRE(ierror == 0)
   else
      allocate(total_left_fe(1, 1, numreps), &
               total_right_fe(1, 1, numreps), &
               num_right_exchg(1, 1, numreps), &
               num_left_exchg(1, 1, numreps), &
               stat=ierror)
      REQUIRE(ierror == 0)
      total_left_fe(:,:,:)   = 0.d0
      total_right_fe(:,:,:)  = 0.d0
      num_right_exchg(:,:,:) = 0
      num_left_exchg(:,:,:)  = 0
   end if

   ! allocate memory for self-guiding temperature table: tempsgtable
   if (isgld > 0) then
      allocate(sgfttable(numreps), sgfftable(numreps), psgldgtable(numreps), stat=ierror)
      REQUIRE(ierror==0)
   end if

   ! initialize the exchange success array to zero (for printing acceptance
   ! rate exchsuccess(i) give the number of times that temperature i
   ! has exch with i+1
   exchsuccess(:,:) = 0

   ! ---=== OPEN REMD FILES ===---
   !  Open remtype_unit to record initial REMD information
   if (master_master .and. rem > 0) then
      call amopen(remtype_unit, remtype, 'U', 'F', 'W')
   end if

   ! R-REMD: load in the file containing energies for the reservoir.
   ! Only sander masters perform the read since they do the exchanges.
   ! Coordinates/velocities for reservoir are loaded as needed during the
   ! RREMD run.
   if (rremd > 0) call load_reservoir_files()


   ! ---=== HYBRID REMD SETUP ===---
   hybridwritetraj=.false.
   if (numwatkeep >= 0) then
      ! 1a- Allocate memory for temp. coord/force storage
      ! Note: This should be the same as in locmem.F90 except am_nbead is
      !  not known. amoeba may not work with hybrid remd.
      allocate( hybrid_coord(3*natom + mxvar), &
                hybrid_force(3*natom + mxvar + 40), &
                hybrid_refc(3*natom + mxvar), stat=ierror)
      REQUIRE( ierror == 0 )

      ! If an output file was specified for hybrid REMD stripped
      !  coords trajectory (-hybridtraj FILE), open it. Only masters
      !  will write to it.
      if (master) then
         if (remstripcoord(1:1) == ' ') then
           write(6,'(a)') &
              "HYBRID REMD: Hybrid stripped traj file will not be written."
         else
            !write(remstripcoord,'(a16,i3.3)') "hybrid.stripped.", nodeid
            hybridwritetraj=.true.
            write(6,'(a,a)') &
               "HYBRID REMD: Opening hybrid stripped coord output: ",&
               remstripcoord
            if (master_master) &
               write(REMLOG_UNIT,'(a)') "# Writing stripped trajectories"
            call amopen(remstripcoord_unit,remstripcoord,'U','F','W')
            write(remstripcoord_unit,'(a80)') "Hybrid stripped coords"
         end if
      end if
   end if

   call mpi_barrier(commworld, ierror)

   ! What to do next depends on rem values
   if (rem == 3) then
      allocate(xtemp(3*natom + mxvar),    &
               ftemp(3*natom + mxvar+40), &
               stat=ierror)
      REQUIRE(ierror==0)
   else if (rem == 4) then
      if (master) then
         call mpi_allgather(solvph, 1, mpi_double_precision, &
                            statetable, 1, mpi_double_precision, &
                            remd_comm, ierror)
         ! Sort phs
         call sorttable(statetable)
         ! Check for duplicate pHs -- but table is sorted already so just check
         ! for adjacent matches
         do i=1,numreps-1
            j=i+1
            if (statetable(i) == statetable(j)) then
               write (6,*) "=============================="
               write (6,*) "pH VALUES OF 2 REPLICAS MATCH!"
               write (6,*) "pH= ", statetable(i)
               write (6,*) "NOT ALLOWED FOR pH REM"
               write (6,*) "=============================="
               call mexit(6,1)
            end if
         enddo
      end if ! master, making sorted statetable
   else if (rem == 5) then
      if (master) then
         call mpi_allgather(solve, 1, mpi_double_precision, &
                            statetable, 1, mpi_double_precision, &
                            remd_comm, ierror)
         ! Sort redoxes
         call sorttable(statetable)
         ! Check for duplicate Es -- but table is sorted already so just check
         ! for adjacent matches
         do i=1,numreps-1
            j=i+1
            if (statetable(i) == statetable(j)) then
               write (6,*) "=============================="
               write (6,*) "E VALUES OF 2 REPLICAS MATCH!"
               write (6,*) "E= ", statetable(i)
               write (6,*) "NOT ALLOWED FOR E REM"
               write (6,*) "=============================="
               call mexit(6,1)
            end if
         enddo
      end if ! master, making sorted statetable
   else
      ! ---=== MAKE SORTED TEMPERATURE TABLE ===---
      ! Now make the sorted temp0 table, used for exchanges. Only masters
      !  need this since they do the exchanges.
      if (master) then
         ! Master processes gather all of the target temperatures
         ! DAN ROE: Since mytargettemp isn't set yet gather temp0;
         !  for rem==2 would need to gather temp0les
         ! Could probably just set mytargettemp=temp0 or
         !  mytargettemp=temp0les and then gather mytargettemp
         call mpi_allgather(temp0, 1, mpi_double_precision, &
                            statetable, 1, mpi_double_precision, &
                            remd_comm, ierror)


         !  RXSGLD parameters
         if (isgld > 0) then
            ! build RXSGLD momentum guiding temperature table
            call mpi_allgather(sgfti, 1, mpi_double_precision, &
                               sgfttable, 1, mpi_double_precision, &
                               remd_comm, ierror)


            ! build RXSGLD force guiding factor table
            call mpi_allgather(sgffi, 1, mpi_double_precision, &
                               sgfftable, 1, mpi_double_precision, &
                               remd_comm, ierror)

            ! build RXSGLD effective force guiding factor table
            call mpi_allgather(psgldg, 1, mpi_double_precision, &
                               psgldgtable, 1, mpi_double_precision, &
                               remd_comm, ierror)

           ! Sort temperatures
           call sorttempsg(numreps, statetable, psgldgtable, sgfftable)
           ! Determine this replca's ID
           if (stagid == 0) then
             stagid = tempsglookup(numreps, temp0, &
               sgfti, sgffi, statetable, sgfttable, sgfftable)
           else
             sgfti = sgfttable(stagid)
             sgffi   = sgfftable(stagid)
             psgldg   = psgldgtable(stagid)
           end if
            ! Check to make sure the first replica has no guiding effect
           if(sgfttable(1) /= 0.0d0 .or. sgfftable(1) /= 0.0d0) then
                  write (6,*) "================================"
                  write (6,*) "No guiding effect on the first replica!"
                  write (6,*) "  sgft, sgff=: ", statetable(1), sgfttable(1), sgfftable(1)
                  write (6,*) "================================"
                  call mexit(6,1)
           end if

         else
            ! Sort temperatures
            call sorttable(statetable)

            ! Check for duplicate temperatures
            ! DAN ROE: Since table is already sorted, we can just check if
            !          any two adjacent temperatures match.
            do i = 1, numreps - 1
               if (statetable(i) == statetable(i+1)) then
                  write (6,*) "================================"
                  write (6,*) "TEMPERATURES OF 2 REPLICAS MATCH!"
                  write (6,*) "T= ",statetable(i)
                  write (6,*) "NOT ALLOWED FOR TEMPERATURE REMD"
                  write (6,*) "================================"
                  call mexit(6, 1)
               end if
            enddo
         end if
      end if ! master, making sorted statetable
   end if ! rem==3

   ! Set our replica indexes by looking up our spot in the sorted statetable
   if (master) then

      if (rem == 1 .or. rem == 2) then

         if (isgld > 0) then
            replica_indexes(1) = stagid
         else
            replica_indexes(1) = templookup(temp0, statetable)
         end if

      else if (rem == 3) then
         replica_indexes(1) = remd_rank + 1
      else if (rem == 4) then
         replica_indexes(1) = templookup(solvph, statetable)
      else if (rem == 5) then
         replica_indexes(1) = templookup(solve, statetable)
      end if

   end if
   ! Open and write out header to remlog file. Only overall master
   !  deals with the remlog.
   if (master_master) then
#ifdef VERBOSE_REMD
      write(6,'(a)') "| REMD: Initializing remlog."
#endif
      call amopen(REMLOG_UNIT, remlog, 'U', 'F', 'W')
      write(REMLOG_UNIT,'(a)') "# Replica Exchange log file"
      write(REMLOG_UNIT,'(a,i10)') "# numexchg is ",numexchg
      if (numwatkeep >= 0) then
         write(REMLOG_UNIT,'(a)') "# HYBRID REMD:"
         write(REMLOG_UNIT,'(a,i3,a)') "#   using igb= ",hybridgb,&
                                       " for exchange energies."
         write(REMLOG_UNIT,'(a,i10)') &
            "#   number of closest waters kept= ",numwatkeep
      end if
      if (rremd>0) write(REMLOG_UNIT,'(a,i3)') "# Reservoir REMD: ",rremd
      if (isgld > 0) then
        ! Write out RXSGLD filenames
        write(REMLOG_UNIT,'(a)') "# RXSGLD(RXSGMD) filenames:"
        write(REMLOG_UNIT,'(a,a)') "#   rxsgldlog= ",trim(remlog)
        write(REMLOG_UNIT,'(a,a)') "#   rxsgldtype= ",trim(remtype)
      else
        ! Write out REMD filenames
        write(REMLOG_UNIT,'(a)') "# REMD filenames:"
        write(REMLOG_UNIT,'(a,a)') "#   remlog= ",trim(remlog)
        write(REMLOG_UNIT,'(a,a)') "#   remtype= ",trim(remtype)
      end if
      if (rremd>0) then
         write(REMLOG_UNIT,'(a,a)') "#   saveene= ",trim(saveenefile)
         if (rremd==3) &
            write(REMLOG_UNIT,'(a,a)') "#   clusterinfo= ",trim(clusterinfofile)
         write(REMLOG_UNIT,'(a,a)') "#   reservoir= ",trim(reservoirname)
      end if
      ! Write out the column labels
      if (rem == 1 .or. rem == 2) then
         if (isgld > 0) then
         write(REMLOG_UNIT,'(a)')"# RXSGLD setup: stagid   temp0    sgft    sgff"
           do i=1,numreps
            write(REMLOG_UNIT,'(a,i4,f10.2,3f10.4)')"#      replica: ", &
                      i, statetable(i), sgfttable(i), sgfftable(i)
           enddo
            write(REMLOG_UNIT,'(a)') &
         "# Rep Stagid Vscale  SGscale Temp Templf Eptot Acceptance(i,i+1)"
         else
            write(REMLOG_UNIT,'(a)') &
            "# Rep#, Velocity Scaling, T, Eptot, Temp0, NewTemp0,&
            & Success rate (i,i+1), ResStruct#"
         end if
      else if (rem == 3) then
         write(REMLOG_UNIT, '(a)') &
           '# Rep#, Neibr#, Temp0, PotE(x_1), PotE(x_2), left_fe,&
           & right_fe, Success, Success rate (i,i+1)'
      else if (rem == 4) then
         write(REMLOG_UNIT, '(a)') &
            "# Rep#, N_prot, old_pH, new_pH, Success rate (i,i+1)"
      else if (rem == 5) then
         write(REMLOG_UNIT, '(a)') &
            "# Rep#, N_elec, old_E, new_E, Success rate (i,i+1)"
      else if (rem == 6) then
         write(REMLOG_UNIT, '(a)') &
            "# Rep#, Velocity Scaling, T, Eptot, Temp0, NewTemp0,&
            & Success rate (i,i+1), ResStruct#"
      end if
   end if ! master_master, remlog header

   return

end subroutine remd1d_setup

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Set up REMD run in multiple dimensions
subroutine multid_remd_setup(numexchg, numwatkeep, temp0, &
                             mxvar, natom, ig, solvph, solve, irest)

   use AmberNetcdf_mod, only : NC_readRestartIndices
   use sander_lib, only : strip, upper

   ! All go to's in this subroutine jump to the end to display an error message
   ! and quit

   implicit none

#  include "parallel.h"

   ! Formal arguments:
   integer, intent(in)     :: numexchg
   integer, intent(in)     :: numwatkeep
   _REAL_, intent(in out)  :: temp0
   integer, intent(in)     :: mxvar
   integer, intent(in)     :: natom
   integer, intent(in out) :: ig
   _REAL_, intent(in out)  :: solvph
   _REAL_, intent(in out)  :: solve
   integer, intent(in)     :: irest

   ! Parameter -- group size, this is because we can't have allocatables in
   ! namelists. This limits the number of dimensions and the number of replicas
   ! in each REMD group to 1000 (this should be major overkill)
   integer, parameter  :: GRPS = 1000

   ! Local variables
   !  &multirem namelist
   integer, dimension(GRPS,GRPS) :: group
   character(80)                 :: exch_type, desc
   integer                       :: repidxIn, crdidxIn

   ! Utility variables (counters, error markers, etc.)
   integer              :: i, j, idx
   integer              :: ifind
   integer              :: ierror
   integer              :: group_counter
   integer              :: add_fac
   integer, allocatable :: replica_assignments(:,:)
   integer, allocatable :: replica_indexes_buf(:,:)

   character(len=5)               :: extension
   character(len=80)              :: buf
   character(len=MAX_FN_LEN)      :: filename
   character(len=80), allocatable :: replica_desc(:)

   namelist / multirem /   group, exch_type, desc

! Variable descriptions:
!
!  group               : This is a set of integers that denote the replica #s
!                        involved in that set of communicating replicas
!  exch_type           : The type of exchange (TEMPERATURE)
!  temptype            : Token to declare an exchange in Temperature space
!  hamtype             : Token to declare an exchange in Hamiltonian space
!  i,j, group_counter  : counters
!  multirem            : namelist to extract information from
!  ifind               : for nmlsrc -- did we find our namelist? 0 = no, 1 = yes
!  alloc_failed        : did memory allocation fail?
!  replica_assignments : buffer for master_master when reading which group
!                        each replica belongs to in each dimension
!  replica_indexes_buf : buffer for reading replica ranks by master_master
!  extension           : file name extension (dimension #) for rem.log
!  filename            : full rem.log file name with ".extension" (above)
!  replica_desc        : The description string for each replica dimension

   ! Force an initial re-read of the NMR variables (e.g., temp0)

   initmodwt = .true.

   ! Initialize group_master_comm

   group_master_comm = mpi_comm_null
   remd_comm = mpi_comm_null

   ! Seed the random number generator

   ! Initialize the RNG (nodeid == repnum - 1).
   add_fac = mod(repnum-1, numreps) / 2 * 2 ! For compatibility with pmemd
   call amrset_gen(remd_gen, ig + add_fac)

   ! Open remlog and write some initial info. Then open and parse remd dimension
   ! file (only one thread does this)

   if (master_master) then

      ! Read the remd_dimension_file
      call amopen(REMD_DIM_UNIT, remd_dimension_file, 'O', 'F', 'R')
      remd_dimension = 0
      do
         call nmlsrc('multirem', REMD_DIM_UNIT, ifind)
         if (ifind == 0) exit
         ! nmlsrc backs up to the start of the namelist, so we have to eat this
         ! line of we want nmlsrc to go to the next &multirem namelist
         read(REMD_DIM_UNIT, '(a80)') buf
         remd_dimension = remd_dimension + 1
      end do

      ! We need at least 1 dimension
      if (remd_dimension == 0) go to 665

      ! Broadcast our remd_dimension so everyone knows in our commmaster
      call mpi_bcast(remd_dimension, 1, mpi_integer, 0, commmaster, ierror)
      ! Allocate our necessary data structures
      allocate(remd_types(remd_dimension),   &
               replica_desc(remd_dimension), &
               replica_assignments(remd_dimension, numgroup), &
               replica_indexes_buf(remd_dimension, numgroup), &
               index_list(numreps),          &
               stat=ierror)
      REQUIRE(ierror == 0)

      replica_assignments(:,:) = 0
      replica_indexes_buf(:,:) = 0

      ! Loop through all of the remd dimensions. At this point, REMD_DIM_UNIT
      ! has been rewound by subroutine nmlsrc
      do i = 1, remd_dimension

         call nmlsrc('multirem', REMD_DIM_UNIT, ifind)

         ! This should never happen because of the preliminary &multirem
         ! counting
         if (ifind == 0) then
            write(6, '(a)') 'ERROR: multid_remd_setup: Should not be here...'
            call mexit(6, 1)
         end if

         ! Initialize our namelist variables
         group(:,:) = 0
         exch_type = ' '
         desc = ' '

         ! Read that namelist
         read(REMD_DIM_UNIT, nml=multirem, err=666)

         ! Make the exch_type upper-case to make it case-insensitive
         call upper(exch_type)

         ! Store our replica description. It's just cosmetic. A form of
         ! self-documentation for the user which will also be dumped to the
         ! rem.log file for users' benefit. It has no impact on the simulation
         replica_desc(i) = desc

         ! Assign the remd_type
         select case(trim(exch_type))
            case ('TEMPERATURE')
               remd_types(i) = 1
            case ('TEMP')
               remd_types(i) = 1
            case ('HAMILTONIAN')
               remd_types(i) = 3
            case ('HREMD')
               remd_types(i) = 3
            case ('PH')
               remd_types(i) = 4
            case ('REDOX')
               remd_types(i) = 5
            case default
               write(6, '(2a)') 'ERROR: Unrecognized EXCH_TYPE ', &
                  trim(exch_type)
               call mexit(6, 1)
         end select

         do group_counter = 1, GRPS
            ! Jump out of the loop if this group is not assigned. That's the
            ! last group in this dimension
            if (group(group_counter, 1) == 0) exit

            ! Assign all of the replicas in this group
            do j = 1, GRPS

               idx = group(group_counter, j)

               ! bail out if this is 0 -- that means we've reached the last
               ! replica in this group
               if (idx == 0) exit

               ! Catch bad replica assignment
               if (idx <= 0 .or. idx > numgroup) go to 671

               ! Assign the replica ensuring each replica is assigned once and
               ! only once
               if (replica_assignments(i, idx) == 0) then
                  replica_assignments(i, idx) = group_counter
                  replica_indexes_buf(i, idx) = j
               else
                  go to 670
               end if

            end do ! j = 1, GRPS

         end do ! group_counter = 1, GRPS

         ! Make sure everyone was assigned here
         do j = 1, numgroup
            if (replica_assignments(i, j) == 0) go to 670
         end do

      end do ! i = 1, remd_dimension

      close(REMD_DIM_UNIT)

   else if (master) then ! master, but not master_master

      ! Receive the broadcast everywhere else
      call mpi_bcast(remd_dimension, 1, mpi_integer, 0, commmaster, ierror)
      allocate(remd_types(remd_dimension), &
               index_list(numreps),        &
               replica_assignments(remd_dimension, numgroup), &
               replica_indexes_buf(remd_dimension, numgroup), &
               stat=ierror)
      REQUIRE(ierror == 0)

   end if ! master_master

   ! Broadcast our remd_dimension for everyone (not just masters)
   call mpi_bcast(remd_dimension, 1, mpi_integer, 0, commsander, ierror)

   if (master) then

      ! Now allocate all data structures that everybody here needs
      allocate(group_num(remd_dimension),             &
               replica_indexes(remd_dimension),       &
               exchsuccess(remd_dimension, numgroup), &
               jumpright(remd_dimension),             &
               statetable(numgroup),                  &
               multid_print_data(numgroup),           &
               xtemp(3*natom + mxvar),                &
               ftemp(3*natom + mxvar+40),             &
               multid_print_data_buf(numgroup, numgroup),           &
               total_left_fe(remd_dimension, numgroup, numgroup),   &
               total_right_fe(remd_dimension, numgroup, numgroup),  &
               num_left_exchg(remd_dimension, numgroup, numgroup),  &
               num_right_exchg(remd_dimension, numgroup, numgroup), &
               stat=ierror)
      REQUIRE(ierror == 0)

      ! Broadcast remd_types to all masters
      call mpi_bcast(remd_types, remd_dimension, mpi_integer, 0, &
                     commmaster, ierror)
      ! Broadcast remd_types to all slaves
      call mpi_bcast(remd_types, remd_dimension, mpi_integer, 0, &
                     commsander, ierror)

      ! Initialize the data
      multid_print_data(:)       = NULL_REMLOG_DATA
      multid_print_data_buf(:,:) = NULL_REMLOG_DATA
      exchsuccess(:,:)           = 0
      jumpright(:)               = .false.
      total_left_fe(:,:,:)       = 0.d0
      total_right_fe(:,:,:)      = 0.d0
      num_left_exchg(:,:,:)      = 0
      num_right_exchg(:,:,:)     = 0

      ! Now set up all data structures for each dimension

      do i = 1, remd_dimension

         ! First scatter our placement in the group/replica ladders in each
         ! dimension
         call mpi_scatter(replica_assignments(i,:), 1, mpi_integer, &
                          group_num(i), 1, mpi_integer, &
                          0, commmaster, ierror)
         call mpi_scatter(replica_indexes_buf(i,:), 1, mpi_integer, &
                          replica_indexes(i), 1, mpi_integer, &
                          0, commmaster, ierror)

      end do ! i = 1, remd_dimension

      if (irest == 1) then
         ierror = NC_readRestartIndices(inpcrd, replica_indexes, group_num, &
                                        repidxIn, crdidxIn, remd_dimension)
         if (repidxIn .gt. -1 .and. crdidxIn .gt. -1) then
           remd_repidx = repidxIn
           remd_crdidx = crdidxIn
           write(6,'(2(a,i6))') '| Overall replica indices from restart: RepIdx=',&
                                remd_repidx, ' CrdIdx=', remd_crdidx
         else
           write(6,'(2(a,i6))') '| Initial overall replica indices: RepIdx=',&
                                remd_repidx, ' CrdIdx=', remd_crdidx
         endif
         if (ierror.ne.0) then ! TODO: mexit if -1?
            write(6,'(a)') '| Warning: Replica indices will NOT be used to &
                           &restart MultiD-REMD run.'
         else
            write(6,'(a)') '| Restarting MultiD-REMD run. This replica will use indices:'
            write(6,'(a,13i6)') '| ', (replica_indexes(i),i=1,remd_dimension)
         end if
      end if

      do i = 1, remd_dimension

         ! Now we set up our remd_comm and get our respective sizes

         remd_comm = mpi_comm_null
         call mpi_barrier(commmaster, ierror)
         call mpi_comm_split(commmaster, group_num(i), &
                             masterrank, remd_comm, ierror)
         call mpi_comm_size(remd_comm, remd_size, ierror)
         call mpi_comm_rank(remd_comm, remd_rank, ierror)

         ! Here we do any additional setup that's needed for any of the exchange
         ! types.

         select case(remd_types(i))

            case(1) ! TEMPERATURE

               hybridwritetraj = .false.
               if (numwatkeep >= 0) then
                  allocate(hybrid_coord(3*natom + mxvar), &
                           hybrid_force(3*natom + mxvar), &
                           hybrid_refc(3*natom + mxvar), stat=ierror)
                  REQUIRE(ierror == 0)
                  if (len_trim(remstripcoord) /= 0) &
                     write(6, '(a)') &
                        'HYBRID REMD: WARNING: stripped coord trajectory is not&
                        & written in multi-D REMD!'
               end if
               ! Collect all of the temps and determine our starting index
               call mpi_allgather(temp0, 1, mpi_double_precision, &
                                  statetable, 1, mpi_double_precision, &
                                  remd_comm, ierror)
               call sorttable(statetable)
               ! Check for duplicate temperatures
               do j = 1, remd_size - 1
                  if (statetable(j) == statetable(j+1)) then
                     write(6,*) 'ERROR: &
                        &TEMPERATURES OF 2 REPLICAS MATCH IN DIMENSION ', i
                     write(6, *) 'Temps ', j, ' and ', j+1, ' are ', &
                                 statetable(j), ' == ', statetable(j+1)
                     call mexit(6, 1)
                  end if
               end do

               replica_indexes(i) = templookup(temp0, statetable)
            ! END case(1) TEMPERATURE

            case(4) ! PH

               call mpi_allgather(solvph, 1, mpi_double_precision, &
                                  statetable, 1, mpi_double_precision, &
                                  remd_comm, ierror)
               call sorttable(statetable)
               ! Check for duplicate pHs
               do j = 1, remd_size - 1
                  if (statetable(j) == statetable(j+1)) then
                     write(6, *) 'ERROR: &
                        &pH VALUES OF 2 REPLICAS MATCH IN DIMENSION', i
                     write(6, *) 'pHs ', j, ' and ', j+1, ' are ', &
                                 statetable(j), ' == ', statetable(j+1)
                     call mexit(6, 1)
                  end if
               end do

               replica_indexes(i) = templookup(solvph, statetable)
            ! END case(4) PH

            case(5) ! REDOX

               call mpi_allgather(solve, 1, mpi_double_precision, &
                                  statetable, 1, mpi_double_precision, &
                                  remd_comm, ierror)
               call sorttable(statetable)
               ! Check for duplicate Es
               do j = 1, remd_size - 1
                  if (statetable(j) == statetable(j+1)) then
                     write(6, *) 'ERROR: &
                        &E VALUES OF 2 REPLICAS MATCH IN DIMENSION', i
                     write(6, *) 'Es ', j, ' and ', j+1, ' are ', &
                                 statetable(j), ' == ', statetable(j+1)
                     call mexit(6, 1)
                  end if
               end do

               replica_indexes(i) = templookup(solve, statetable)
            ! END case(5) REDOX

         end select
         ! Now we set up the rem.log files. Each dimension will get its own
         ! rem.log taking the -remlog <filename> and appending .1, .2, .3, ...
         ! etc. to it based on the dimension
         if (master_master) then

            write(extension, '(i5)') i
            call strip(extension)

            if (len_trim(remlog) + len_trim(extension) + 1 > MAX_FN_LEN) then
               write(6, '(2a)') 'ERROR: rem.log file name overflow. Increase &
                  &MAX_FN_LEN in file_io_dat.F90 and recompile'
                  call mexit(6, 1)
            end if

            filename = trim(remlog) // '.' // trim(extension)
            call amopen(REMLOG_UNIT, filename, owrite, 'F', 'W')

            write(REMLOG_UNIT, '(a)')       '# Replica Exchange log file'
            write(REMLOG_UNIT, '(a,i10)')   '# numexchg is ', numexchg
            write(REMLOG_UNIT, '(2(a,i4))') '# Dimension ', i, ' of ', &
                                             remd_dimension
            if (len_trim(replica_desc(i)) > 0) &
               write(REMLOG_UNIT, '(2a)')   '# Description: ', &
                                            trim(replica_desc(i))

            select case (remd_types(i))
               case(1)
                  write(REMLOG_UNIT, '(a)') '# exchange_type = TEMPERATURE'
                  write(REMLOG_UNIT, '(a)') '# REMD filenames:'
                  write(REMLOG_UNIT, '(2a)') '# remlog= ', trim(filename)
                  write(REMLOG_UNIT, '(2a)') '# remd dimension file= ', &
                                             trim(remd_dimension_file)
                  write(REMLOG_UNIT, '(a)') '# Rep#, Velocity Scaling, T, &
                      &Eptot, Temp0, NewTemp0, Success rate (i,i+1), ResStruct#'
               case(3)
                  write(REMLOG_UNIT, '(a)') '# exchange_type = HAMILTONIAN'
                  write(REMLOG_UNIT, '(a)') '# REMD filenames:'
                  write(REMLOG_UNIT, '(2a)') '# remlog= ', trim(filename)
                  write(REMLOG_UNIT, '(2a)') '# remd dimension file= ', &
                                        trim(remd_dimension_file)
                  write(REMLOG_UNIT, '(a)') '# Rep#, Neibr#, Temp0, PotE(x_1), &
                   &PotE(x_2), left_fe, right_fe, Success, Success rate (i,i+1)'
               case(4)
                  write(REMLOG_UNIT, '(a)') '# exchange_type = PH'
                  write(REMLOG_UNIT, '(a)') '# REMD filenames:'
                  write(REMLOG_UNIT, '(2a)') '# remlog= ', trim(filename)
                  write(REMLOG_UNIT, '(2a)') '# remd dimension file= ', &
                                             trim(remd_dimension_file)
                  write(REMLOG_UNIT, '(a)') '# Rep#, N_prot, old_pH, new_pH, &
                     &Success rate (i,i+1)'
               case(5)
                  write(REMLOG_UNIT, '(a)') '# exchange_type = REDOX'
                  write(REMLOG_UNIT, '(a)') '# REMD filenames:'
                  write(REMLOG_UNIT, '(2a)') '# remlog= ', trim(filename)
                  write(REMLOG_UNIT, '(2a)') '# remd dimension file= ', &
                                             trim(remd_dimension_file)
                  write(REMLOG_UNIT, '(a)') '# Rep#, N_elec, old_E, new_E, &
                     &Success rate (i,i+1)'
            end select

            close(REMLOG_UNIT)

         end if ! master_master

         ! Free the remd_comm
         call mpi_comm_free(remd_comm, ierror)

      end do ! i = 1, remd_dimension

   else   ! .NOT. master
      ! Allocate our data structures
      allocate(remd_types(remd_dimension), &
               xtemp(3*natom + mxvar),     &
               ftemp(3*natom + mxvar+40),  &
               stat=ierror)
      REQUIRE(ierror == 0)

      ! Receive remd_types
      call mpi_bcast(remd_types, remd_dimension, mpi_integer, 0, &
                     commsander, ierror)

      ! Do slave setup for each of the dimensions
      ! masters at the bcast call there
      do i = 1, remd_dimension

         if (remd_types(i) == 1 .and. rremd > 0) call load_reservoir_files()

         if (numwatkeep >= 0) then
            allocate(hybrid_coord(3*natom + mxvar), &
                     hybrid_force(3*natom + mxvar), &
                     hybrid_refc(3*natom + mxvar), stat=ierror)
            REQUIRE(ierror == 0)
         end if

      end do

   end if ! master

   ! Set the first value for next_rem_method
   next_rem_method = remd_types(1)

   ! Deallocate work arrays here
   if (allocated(replica_assignments)) deallocate(replica_assignments)
   if (allocated(replica_indexes_buf)) deallocate(replica_indexes_buf)
   if (allocated(replica_desc)) deallocate(replica_desc)

   return

! Error handlers

665   write(6, '(3a)') 'ERROR: Could not find &multirem in ', &
                       trim(remd_dimension_file), '!'
      call mexit(6, 1)

666   write(6, '(a,i3,3a)') 'ERROR: Could not read the ', i, &
                 'th &multirem namelist in ', trim(remd_dimension_file), '!'
      call mexit(6, 1)

670   write(6, '(a)') 'ERROR: Each replica must be assigned to one and &
                      &only one group in each dimension!'
      call mexit(6, 1)

671   write(6, '(a,i5)') 'ERROR: Bad replica assignment. Replicas must be &
                         &between 1 and ', numgroup
      call mexit(6, 1)

end subroutine multid_remd_setup

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Close REMD files, free memory
subroutine remd_cleanup()
   use AmberNetcdf_mod, only: NC_close

   implicit none
#  include "parallel.h"

   integer  :: ierror

   ! Close files
   if (master_master) then
      close(REMLOG_UNIT)
      close(remtype_unit)
   end if
   if (master .and. hybridwritetraj) close(remstripcoord_unit)
   if (master .and. rremd > 0) close(remin_unit)

   ! Deallocate multi-D REMD data structures
   if (allocated(jumpright)) deallocate(jumpright)
   if (allocated(remd_types)) deallocate(remd_types)
   if (allocated(group_num)) deallocate(group_num)
   if (allocated(multid_print_data)) deallocate(multid_print_data)
   if (allocated(multid_print_data_buf)) deallocate(multid_print_data_buf)
   if (allocated(replica_indexes)) deallocate(replica_indexes)
   if (allocated(index_list)) deallocate(index_list)

   ! Exchange success array and state table array
   if (allocated(exchsuccess)) deallocate(exchsuccess)
   if (allocated(statetable)) deallocate(statetable)

   ! Hamiltonian-REMD data structures
   if (allocated(xtemp)) deallocate(xtemp)
   if (allocated(ftemp)) deallocate(ftemp)
   if (allocated(total_left_fe)) deallocate(total_left_fe)
   if (allocated(total_right_fe)) deallocate(total_right_fe)
   if (allocated(num_left_exchg)) deallocate(num_left_exchg)
   if (allocated(num_right_exchg)) deallocate(num_right_exchg)

   ! Hybrid remd temp. storage
   if (allocated(hybrid_coord)) deallocate(hybrid_coord)
   if (allocated(hybrid_force)) deallocate(hybrid_force)
   if (allocated(hybrid_refc)) deallocate(hybrid_refc)

   ! Reservoir REMD storage
   if (allocated(saveene)) deallocate(saveene)
   if (reservoir_ncid.ne.-1) call NC_close(reservoir_ncid)

   ! DihedralCluster arrays
   if (allocated(clusternum)) deallocate(clusternum)
   if (allocated(clusterid)) deallocate(clusterid)
   if (allocated(clustersize)) deallocate(clustersize)
   if (allocated(dihclustat)) deallocate(dihclustat)
   if (allocated(currdihid)) deallocate(currdihid)
   if (allocated(dihclustmin)) deallocate(dihclustmin)
   if (allocated(dihcluststep)) deallocate(dihcluststep)

   ! RXSGLD
   if (allocated(sgfttable)) deallocate(sgfttable)
   if (allocated(sgfftable)) deallocate(sgfftable)

   ! Free communicators
   if (remd_comm /= mpi_comm_null) &
      call mpi_comm_free(remd_comm, ierror)
   if (group_master_comm /= mpi_comm_null) &
      call mpi_comm_free(group_master_comm, ierror)

   return

end subroutine remd_cleanup

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Perform actions necessary to exchange replicas in REMD
subroutine remd_exchange(rem_dim, rem_kind, x, v, amass, nr3, &
                         natom, nr, temp0)

   use sgld, only : isgld, rxsgld_scale

implicit none

#  include "parallel.h"

! Passed variables
   _REAL_, intent(in out) :: x(*)      ! Atomic coordinates
   _REAL_, intent(in out) :: v(*)      ! Atomic velocities
   _REAL_, intent(in)     :: amass(*)  ! Mass array (not inverted)
   integer, intent(in)    :: nr3       ! 3 * number of atoms
   integer, intent(in)    :: natom     ! number of atoms
   integer, intent(in)    :: nr        ! number of atoms (?)
   integer, intent(in)    :: rem_dim   ! Which dimension we are exchanging in
   integer, intent(in)    :: rem_kind  ! Which kind of remd are we doing?
   _REAL_, intent(in out) :: temp0     ! Current target temperature

   integer ierror

!-------------------
#  ifdef VERBOSE_REMD
   if (remd_master) &
     write(6,'(a,26("="),a,26("="))') '| ', "REMD EXCHANGE CALCULATION"
#  endif

   ! Set the partners array so we know who our partners are. NOTE, this also
   ! sets index_list
   if (master .and. .not. isgld > 0) call set_partners(rem_dim, remd_size)

   ! ---=== REMD EXCHANGE CALCULATION ===---
   ! Attempt Exchange with subrem() using energy from previous runmd call
   ! All procs call subrem so slave random # generators stay in sync
   call subrem(rem_dim)

   ! ---=== REMD COMMUNICATION ===---
   ! Broadcast my_remd_data
   call mpi_bcast(my_remd_data, SIZE_REMD_DATA, mpi_double_precision, &
                  0, commsander, ierror)

   if (isgld > 0) then
   ! ---=== RXSGLD property SCALING ===---
      if (my_remd_data%newtargettemp > 0.0) &
         temp0 = my_remd_data%newtargettemp
      call rxsgld_scale(stagid, nr, my_remd_data%myscaling, amass, x, v)
   else
      ! ---=== REMD VELOCITY SCALING ===---
      ! REMD: If an attempt is accepted, set the target temperature to
      !  the new one and rescale velocities.
#ifdef LES
      call remd_scale_velo(v, temp0, nr, nr3, rem_kind)
#else
      call remd_scale_velo(v, temp0, nr3, rem_kind)
#endif
   end if

   ! initmodwt forces modwt() in nmr.f to re-read values such as temp0 -
   ! otherwise they will be reset to the initial value after each
   ! exchange.
   initmodwt = .true.
#  ifdef VERBOSE_REMD
   if (remd_master) &
     write(6,'(a,24("="),a,24("="))') '| ', "END REMD EXCHANGE CALCULATION"
#  endif
   return

end subroutine remd_exchange

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ calculation of exchange probability for T-REMD, RXSGLD
subroutine subrem(rem_dim)

   use constants, only : TWO

   use sgld, only : isgld, tsgset, sgfti, sgffi, psgldg, epotlf, epotllf, &
                     templf,myscalsg, sgld_exchg, &
                    tempsglookup, stagidlookup

   implicit none

#  include "parallel.h"

   ! Passed arguments
   integer, intent (in) :: rem_dim     ! The dimension we are exchanging in

   ! Local variables
   _REAL_ delta, straw, o_scaling, o_scalsg, metrop
   _REAL_ pvterm
   _REAL_ l_temp0, r_temp0, o_temp0
   _REAL_ l_eptot, r_eptot, o_eptot
   _REAL_ o_stagid,o_psgldg,o_sgff,o_epotlf,o_epotllf,o_tlf
   _REAL_ dbeta,dbetamiu,depot,depotlf,depotllf
   _REAL_ alltlf(numreps), alleplf(numreps), allepllf(numreps)
   _REAL_ d_scalsg(numreps)
   integer allstagid(numreps)
   _REAL_ alltempi(numreps), alltemp0(numreps), alleptot(numreps)
   ! For use in remlog output
   _REAL_ d_scaling(numreps), d_o_temp0(numreps)
   ! DAN ROE: for recording RREMD structure output in log
   integer d_rremd_idx(numreps)
   ! index is the position of the replica's T in the sorted T list
   integer myindex, l_index, r_index, o_index
   ! repnum is the actual replica number
   integer my_repnum,l_repnum,r_repnum,o_repnum

   integer i, ierror, istatus(mpi_status_size)

   ! temporary for exchange success, will allgather from all nodes then
   ! add sum to actual exchsuccess array. also define an array for allgather sum
   integer texchsuccess(numreps), tempsuccess(numreps)
   _REAL_ exchfrac(numreps)

   logical exchange

! ---------------------
   ! Initialize variables
   texchsuccess(:) = 0

   delta = 0.0d0
   straw = 0.0d0
   metrop= 0.0d0
   pvterm= 0.0d0
   o_scaling = -0.1d0
   i = 0
   myindex = 0
   l_index = 0
   r_index = 0
   o_index = remd_rank

   ! Call amrand to get a random# which can eventually be used as a structure
   ! index in rremd. Do it here so the random # generators on all the
   ! threads stay in sync. Only the master of the highest T replica during
   ! jumpright will ever use this.
   ! Calling this even when we have no rremd will allow random # gen to stay
   ! in sync between remd and rremd runs.
   call amrand_gen(remd_gen, straw)

   ! Only the sander masters do the exchanges.
   if (master) then

      ! Gather current temperature, target temperature, and PE
      call mpi_allgather(my_remd_data%mytemp, 1, mpi_double_precision, &
                         alltempi, 1, mpi_double_precision, &
                         remd_comm, ierror)
      call mpi_allgather(my_remd_data%mytargettemp, 1, mpi_double_precision, &
                         alltemp0, 1, mpi_double_precision, &
                         remd_comm, ierror)
      call mpi_allgather(my_remd_data%myEptot, 1, mpi_double_precision, &
                         alleptot, 1, mpi_double_precision, &
                         remd_comm, ierror)

      ! alltemp0 and alleptot are indexed by replica #
      ! statetable is SORTED
      if (isgld > 0) then
         call mpi_allgather(Epotlf, 1, mpi_double_precision, &
                         alleplf, 1, mpi_double_precision, &
                         remd_comm, ierror)
         call mpi_allgather(epotllf, 1, mpi_double_precision, &
                         allepllf, 1, mpi_double_precision, &
                         remd_comm, ierror)
         call mpi_allgather(templf, 1, mpi_double_precision, &
                         alltlf, 1, mpi_double_precision, &
                         remd_comm, ierror)
         call mpi_allgather(stagid, 1, mpi_integer, &
                         allstagid, 1, mpi_integer, &
                         remd_comm, ierror)
      ! Find our position in the stage list
         myindex = stagid

      else

         myindex = replica_indexes(rem_dim)

      end if

      ! DAN ROE: Why do we look for both partners?

      ! Find neighbor Ts for exchange using Temperature list indexed
      !  by TEMPERATURE.
      ! Wrap around so replicas exchange in a circle,
      !  i.e. highest and lowest T may attempt exchange.
      l_index = myindex - 1
      if (l_index < 1) l_index = remd_size
      r_index = myindex + 1
      if (r_index > remd_size) r_index = 1

      if (isgld > 0) then
         my_repnum = stagidlookup(remd_size, myindex, allstagid)
         l_repnum = stagidlookup(remd_size, l_index, allstagid)
         r_repnum = stagidlookup(remd_size, r_index, allstagid)
         l_temp0 = statetable(l_index)
         r_temp0 = statetable(r_index)
         if(l_index /= allstagid(l_repnum) .or. &
            r_index /= allstagid(r_repnum) .or. &
            myindex /= allstagid(my_repnum))then
                  write (6,*) "================================"
                  write (6,*) "Stagid mismatch!"
                  write (6,*) "left : ",l_index , allstagid(l_repnum),l_repnum
                  write (6,*) "self : ",myindex , allstagid(my_repnum),my_repnum
                  write (6,*) "right: ",r_index , allstagid(r_repnum),r_repnum
                  write (6,*) "================================"
            call mexit(6,1)
         end if
      else ! NOT RXSGLD

         ! Now pull the replica number of the adjacent replicas (partners array)
         l_repnum = partners(1)
         r_repnum = partners(2)

         ! Now get their temperatures
         l_temp0 = alltemp0(l_repnum)
         r_temp0 = alltemp0(r_repnum)
      end if

      ! Get the energies of these neighbor replicas
      l_eptot = alleptot(l_repnum)
      r_eptot = alleptot(r_repnum)


      ! Set up partner replica information.
      ! By definition, only even replica #s will initiate the exchange.
      if(mod(myindex, 2) == 0) then
         ! this replica will calculate delta
         if (jumpright(rem_dim)) then
            ! Partner is to the right, i.e. at higher T
            ! DAN ROE: Integrate check for reservoir here?
            o_repnum = r_repnum
            o_index  = r_index
            o_eptot  = r_eptot
            o_temp0  = r_temp0
         else
            ! Partner is to the left, i.e. at lower T
            o_repnum = l_repnum
            o_index  = l_index
            o_eptot  = l_eptot
            o_temp0  = l_temp0
         end if
      else
         ! This replica will not calculate delta, right and left are
         !  switched compared to controlling replica.
         ! Only needs o_eptot for writing debug info
         if (jumpright(rem_dim)) then
            o_repnum = l_repnum
            o_index  = l_index
            o_eptot  = l_eptot
            o_temp0  = l_temp0
         else
            o_repnum = r_repnum
            o_index  = r_index
            o_eptot  = r_eptot
            o_temp0  = r_temp0
         end if
      end if  ! (even/odd rank check)

      ! Calculate pressure/volume correction if necessary
      if (use_pv) pvterm = pv_correction(my_remd_data%mytargettemp, o_temp0, o_repnum-1, remd_comm)

      ! Calculate the exchange probability if even # replica.
      if(mod(myindex, 2) == 0) then
         if (isgld > 0) then
            ! Replica exchange self-guided Langevin dynamics
            !  Works also for temperature-based replica exchange
            ! delta = D(beta)D(Ep)+D(beta*miu)*(D(Eplf)-D(Epllf))
            !  D(beta)=1/kT(neighbor)-1/kT
            !  D(betamiu)=(sgff(neighbor)-psgldg(neighbor))/kT(neighbor)-(sgff-psgldg)/kT
            !  D(Ep)=Eppt(neighbor)-Eppt
            !  D(Eplf)=Epotlf(neighbor)-Epotlf
            !  D(Epllf)=Epotllf(neighbor)-Epotllf
   
                 o_sgff=sgfftable(o_index)
                 o_psgldg=psgldgtable(o_index)
                 o_epotlf=alleplf(o_repnum)
                 o_epotllf=allepllf(o_repnum)
                 o_epotlf=alleplf(o_repnum)
                 dbeta = 503.01d0/ o_temp0-503.01d0 / my_remd_data%mytargettemp
                 dbetamiu=503.01d0*(o_sgff-o_psgldg)/ o_temp0-  &
                          503.01d0*(sgffi-psgldg) / my_remd_data%mytargettemp
                 depot=o_eptot- my_remd_data%myEptot 
                 depotlf=o_epotlf - epotlf 
                 depotllf=o_epotllf - epotllf 
             
                 delta = -dbeta * depot - dbetamiu*(depotlf-depotllf) 
         else
            ! Std REMD
            delta = (my_remd_data%myEptot - o_eptot) &
                  * (my_remd_data%mytargettemp - o_temp0) * 503.01d0 &
                  / (my_remd_data%mytargettemp * o_temp0)
         end if ! REMD exchange calculation

         metrop = exp(-delta + pvterm)

         ! Get random number between 0 and 1
         call amrand_gen(remd_gen, straw)

         ! Check for exchange
         exchange = straw < metrop

         ! If exchanged, set scaling, otherwise dont
         if (exchange) then
            ! Set velocity scaling
            my_remd_data%myscaling = sqrt(o_temp0 / my_remd_data%mytargettemp)
            o_scaling = 1.0d0 / my_remd_data%myscaling

            ! RXSGLD scaling factor
            if (isgld > 0) then
               myscalsg = sqrt(alltlf(o_repnum) / templf)
               o_scalsg = 1.0d0 / myscalsg
            end if

            ! Increment exchsuccess for lower replica #
            ! Use index, not repnum, since this is a property of temperatures
            if (jumpright(rem_dim)) then
               ! this is the lower temperature since this is the controlling rep
               ! NOTE THAT RANKS ARE 0 TO NUMREPS-1 BUT USE 1 TO NUMREPS FOR
               ! SUCCESS
               texchsuccess(myindex) = 1
            else
               ! other is the lower temperature
               texchsuccess(myindex-1) = 1
            end if
         else
            ! No exchange
            my_remd_data%myscaling = -1.0d0
            o_scaling = -1.0d0
            myscalsg = -1.0d0
            o_scalsg = -1.0d0
         end if ! exchange

! CARLOS: DO WE NEED BARRIER?

         call mpi_barrier(remd_comm, ierror)

         ! send the results to the partner
         ! NOTE THAT WE NEED TO SUBTRACT 1 FROM REPNUM SINCE SUBREM USES
         ! INDEX OF 1 TO NUMREPS WHILE THE MPI PROCESSES ARE 0 TO NUMREPS-1

         !write (6,*) "sending o_scaling of ",o_scaling,&
         !            "to replica ",o_repnum
         call mpi_send(o_scaling, 1, mpi_double_precision, &
                       o_repnum-1, 0, remd_comm, ierror)
         if (isgld > 0) &
            call mpi_send(o_scalsg, 1, mpi_double_precision, &
                          o_repnum-1, 501, remd_comm, ierror)
      else
         ! Not the replica controlling the exchange
         ! call rand to keep in sync with replicas that calculated Metropolis
         call amrand_gen(remd_gen, straw)
         ! DAN ROE: FileDebug

! CARLOS: NEED BARRIER?
         call mpi_barrier(remd_comm, ierror)
         ! receive the scaling data; if it is >0 we know exchange succeeded
         call mpi_recv(my_remd_data%myscaling, 1, mpi_double_precision, &
                       o_repnum-1, 0, remd_comm, istatus, ierror)
         if (isgld > 0) &
            call mpi_recv(myscalsg, 1, mpi_double_precision, &
                          o_repnum-1, 501, remd_comm, istatus, ierror)
      end if ! replica controlling the exchange (even #)

      ! toggle exchange direction

      call mpi_barrier(remd_comm, ierror)
      jumpright(rem_dim) = .not. jumpright(rem_dim)

      if(my_remd_data%myscaling < 0.0d0) then
         my_remd_data%newtargettemp = my_remd_data%mytargettemp
         exchange=.false.
      else
         my_remd_data%newtargettemp = o_temp0
         exchange=.true.
         if (isgld > 0) then
            stagid=allstagid(o_repnum)
            sgfti=sgfttable(stagid)
            sgffi=sgfftable(stagid)
            psgldg=psgldgtable(stagid)
            tsgset=o_temp0
            o_stagid = tempsglookup(remd_size, tsgset, sgfti, sgffi, &
                                    statetable, sgfttable, sgfftable)
            if(stagid /= o_stagid) &
               write(6,*) "Problem in exchange ID!", stagid, o_stagid
         end if
      end if

      if (exchange.and.o_index.gt.0) then
         ! Swap replica index and all of our group affiliations
         ! Do not do this if exchanging with the reservoir.
         call mpi_sendrecv_replace(replica_indexes, remd_dimension, &
                     mpi_integer, o_repnum-1, 22, o_repnum-1, 22, &
                     remd_comm, istatus, ierror)
         call mpi_sendrecv_replace(group_num, remd_dimension, mpi_integer, &
                     o_repnum-1, 23, o_repnum-1, 23, remd_comm, istatus, ierror)
         call mpi_sendrecv_replace(remd_repidx, 1, mpi_integer, &
                                   o_repnum-1, 24, o_repnum-1, 24, &
                                   remd_comm, istatus, ierror)
      end if

      ! REM: gather exchange log data

      call mpi_gather(my_remd_data%myscaling, 1, mpi_double_precision, &
                      d_scaling, 1, mpi_double_precision, &
                      0, remd_comm, ierror)
      if (isgld > 0) &
         call mpi_gather(myscalsg, 1, mpi_double_precision, &
                         d_scalsg, 1, mpi_double_precision, &
                         0, remd_comm, ierror)
      call mpi_gather(my_remd_data%newtargettemp, 1, mpi_double_precision, &
                      d_o_temp0, 1, mpi_double_precision, &
                      0, remd_comm, ierror)

       ! DAN ROE: RREMD: gather reservoir structure data
      call mpi_gather(rremd_idx, 1, MPI_INTEGER, &
                      d_rremd_idx, 1, MPI_INTEGER, &
                      0, remd_comm, ierror)

      call mpi_allreduce(texchsuccess, tempsuccess, remd_size,  MPI_INTEGER, &
                         MPI_SUM, remd_comm, ierror)

      ! add the current successes to overall array

      do i = 1, remd_size
         exchsuccess(rem_dim, i) = exchsuccess(rem_dim, i) + tempsuccess(i)

         ! multiple fraction by 2.0 since we attempt this pair every OTHER
         ! exchange attempt (alternating directions of exchange)

         exchfrac(i) = dble(exchsuccess(rem_dim, i)) &
                     / dble(mdloop) * dble(remd_dimension) * TWO
      enddo


      if (remd_master .and. rem > 0) then

      ! only overall master writes the log
      ! DAN ROE: added d_rremd_idx write

         write(unit = REMLOG_UNIT, fmt = '(a,i8)') '# exchange ', mdloop
         do i = 1, remd_size
            if (isgld > 0) then
               write(REMLOG_UNIT, '(i4,i4,2f8.4,2f8.2,e14.6,f8.4)') &
               i, & ! Replica #
               allstagid(i), & ! current temperature
               d_scaling(i), & ! scaling factor
               d_scalsg(i), & ! scaling factor
               alltempi(i), & ! current temperature
               alltlf(i), & ! current temperature
               alleptot(i), & ! current potential energy
               exchfrac(allstagid(i)) ! current exchange success fraction
            else
               write(REMLOG_UNIT, '(i2,f12.4,f10.2,f15.2,3f10.2,i8)') &
               i, & ! Replica #
               d_scaling(i), & ! scaling factor
               alltempi(i), & ! current temperature
               alleptot(i), & ! current potential energy
               alltemp0(i), & ! current target temperature
               d_o_temp0(i),& ! next target temperature
               exchfrac(index_list(i)), &
               d_rremd_idx(i) ! structure#
           end if
         end do
      else if (remd_master) then
         do i = 1, remd_size
            multid_print_data(i)%scaling     = d_scaling(i)
            multid_print_data(i)%real_temp   = alltempi(i)
            multid_print_data(i)%pot_ene_tot = alleptot(i)
            multid_print_data(i)%temp0       = alltemp0(i)
            multid_print_data(i)%new_temp0   = d_o_temp0(i)
            multid_print_data(i)%num_rep     = remd_size
            multid_print_data(i)%group_num   = group_num(rem_dim)
            multid_print_data(i)%success_ratio = &
               dble(exchsuccess(rem_dim, index_list(i))) / &
               dble(mdloop) * dble(remd_dimension) * 2
         end do
      end if

      ! DAN ROE: All sander masters write REMD exchange info
         write(6,'(a6,i10,a8,i1)') "| Exch= ", mdloop, " RREMD= ", rremd
         ! This Replica Information
         write(6,'(a16,a7,f6.2,2(a7,i2),a7,f10.2)') &
          "| Replica         "," Temp= ", my_remd_data%mytargettemp, &
          " Indx= ", myindex, " Rep#= ", repnum, " EPot= ", my_remd_data%myEptot
         write(6,'(a7,i6,a8,i6)') '| RepIdx=', remd_repidx, ' CrdIdx=', remd_crdidx
         ! Partner Information
         ! use .not.jumpright since jumpright has already been toggled
         ! Partner is normal replica
         write(6,'(a16,a7,f6.2,2(a7,i2),a7,f10.2)') &
            "| Partner         "," Temp= ",o_temp0," Indx= ",o_index, &
            " Rep#= ", o_repnum, " EPot= ", o_eptot
         ! Exchange calculation information
         if (mod(myindex,2)==0) then
            ! This replica controlled exchange and calculated metrop
            write(6,'(a8,E16.6,a8,E16.6,a12,f10.2)') &
               "| Metrop= ",metrop," delta= ",delta," o_scaling= ",o_scaling
         else
            ! This replica did not control exchange
            write(6,'(a)') "| Not controlling exchange."
         end if ! mod(myindex,2)
         ! Write random #, scaling, and success
         write(6,'(a8,E16.6,a12,f10.2,a10,L1)' ) &
            "| Rand=   ",straw," MyScaling= ", my_remd_data%myscaling, &
            " Success= ",exchange

   else  ! not part of remd_comm
      ! Keep random generator in sync
      call amrand_gen(remd_gen, straw)
   end if

   return
end subroutine subrem

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Scale velocities based on new temps after exchange
#ifdef LES
subroutine remd_scale_velo(v, temp0, nr, nr3, rem_kind)

   use les_data, only : temp0les, cnum
#else
subroutine remd_scale_velo(v, temp0, nr3, rem_kind)
#endif /* LES */

   implicit none
#  include "parallel.h"

   _REAL_, intent(in out) :: v(*)
   _REAL_, intent(in out) :: temp0
#ifdef LES
   integer, intent(in)    :: nr
#endif /* LES */
   integer, intent(in)    :: nr3
   integer, intent(in)    :: rem_kind

   integer i

!--------------------

   ! REMD: If an attempt is accepted, set the target temperature to
   !  the new one and rescale velocities.
   ! DAN ROE: Eventually take out debug info.
   if (rem_kind == 1 .or. rem_kind == 3) then
#ifdef VERBOSE_REMD
      if (master) &
         write (6,'(2(a,f6.2))') &
            "| REMD: checking to see if bath T has changed: ", &
            temp0, "->", my_remd_data%newtargettemp
#endif
      ! All processes set temperature. newtargettemp is set in subrem
      if (my_remd_data%newtargettemp > 0.0) &
         temp0 = my_remd_data%newtargettemp
      if (my_remd_data%myscaling > 0.0) then
         ! All processes scale velocities.
         ! DAN ROE: This could potentially be divided up as in runmd
         !  since when there are multiple threads per group each thread
         !  only ever knows about its own subset of velocities anyway.
#ifdef VERBOSE_REMD
         if (master) then
            write (6,'(a,f8.3,a,f8.3)') &
               "| REMD: scaling velocities by ", my_remd_data%myscaling,&
               " to match new bath T ", temp0
         end if
#endif
         do i = 1, nr3
            v(i) = v(i) * my_remd_data%myscaling
         enddo
      end if
#  ifdef LES
   elseif (rem_kind == 2) then
      if (my_remd_data%newtargettemp > 0.0) &
         temp0les = my_remd_data%newtargettemp
      if (my_remd_data%myscaling > 0.0) then
         do i = 1, nr
            if (cnum(i) > 0.0) then
               v(3*i-2) = v(3*i-2) * my_remd_data%myscaling
               v(3*i-1) = v(3*i-1) * my_remd_data%myscaling
               v(3*i)   = v(3*i)   * my_remd_data%myscaling
            end if
         enddo
      end if
#  endif  /* LES */
   end if ! rem_kind == 1, velocity scaling and bath T

   return

end subroutine remd_scale_velo


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! lookup temp in templist and return its index
integer function templookup(temp, templist)

   implicit none
#  include "parallel.h"

   _REAL_, intent(in) :: temp
   _REAL_, dimension(remd_size), intent(in) :: templist
   integer i

   templookup = -1
   do i = 1, remd_size
      if(abs(temp-templist(i)) < 1.0d-6) then
         templookup = i
         return
      end if
   end do
   REQUIRE( templookup > 0 )

end function templookup

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Determines partner replicas in a given dimension from everyone's ranks
subroutine set_partners(rem_dim, num_replicas)

   ! ALSO SETS index_list

   implicit none


   ! Passed Variables

   integer, intent(in) :: rem_dim
   integer, intent(in) :: num_replicas

   ! Local variables

   integer  :: lower_neibr
   integer  :: higher_neibr
   integer  :: my_idx
   integer  :: ierror

   integer  :: i ! counter

   ! Make sure only masters are here
   if (.not. master) return

   my_idx = replica_indexes(rem_dim)

   if (my_idx == 1) then
      lower_neibr = num_replicas
      higher_neibr = 2
   else if (my_idx == num_replicas) then
      lower_neibr = my_idx - 1
      higher_neibr = 1
   else
      lower_neibr = my_idx - 1
      higher_neibr = my_idx + 1
   end if

   ! Find out everyone's rank in this REMD dimension, then search for our
   ! neighbors.  NOTE that when we find our partners, they will be indexed
   ! starting from 1, NOT 0.

   call mpi_allgather(my_idx, 1, mpi_integer, index_list, 1, &
                      mpi_integer, remd_comm, ierror)

   partners(:) = -1

   do i = 1, num_replicas
      if (index_list(i) == lower_neibr)  partners(1) = i
      if (index_list(i) == higher_neibr) partners(2) = i
   end do

   if (partners(1) == -1 .or. partners(2) == -1) then
      write(6, '(a)') 'Creation of partners array failed! Bad REMD setup!'
      write(6, '(a)') 'Check that the numbering of replicas in the replica dimension file is'
      write(6, '(a)') 'consistent with the total number of replicas (and if present, that'
      write(6, '(a)') 'the replica indices in the input coordinates are as well).'
      do i = 1, num_replicas
        write(6, '(2(a,i6))') '    index_list(', i, ')= ', index_list(i)
      end do
      call mexit(6, 1)
   end if

   return

end subroutine set_partners

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ sort temp ascendingly
subroutine sorttable(temp)

   implicit none

#  include "parallel.h"

   _REAL_, dimension(remd_size), intent(inout) :: temp

   _REAL_ tempt
   integer i, j

   do i = 1, remd_size
      do j = i + 1, remd_size
         if(temp(j) < temp(i)) then
            tempt = temp(i)
            temp(i) = temp(j)
            temp(j) = tempt
         end if
      end do
   end do
end subroutine sorttable

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Load the specified structure from the reservoir into coords
subroutine load_reservoir_structure(x, v, nr3, natom)
   use netcdf
   use AmberNetcdf_mod, only: NC_error
   implicit none
#  include "parallel.h"

   _REAL_, dimension(*), intent(inout) :: x, v
   integer, intent(in) :: nr3, natom

   integer i, ierror
   character(len=80) line
   character(len=90) framename

!--------------------

   ! ----===== RREMD RESERVOIR LOADING =====----
   ! If we exchanged with the reservoir, swap the coordinates
   !  here (after the inpcrd have been read), Master process only.
   if (master) then
      ! Read restart file into coords. Need to set a filename for
      !  correct structure.
      write (6,'(a34)') "=========Reservoir Read==========="
      ! DAN ROE: integer size will need to correspond to maxreservoir
      ! Debug
      write (6,'(a,a)') "reservoirname=",trim(reservoirname)
      write (6,'(a,i6.6)') "rremd_idx=",rremd_idx
      if (reservoir_ncid.eq.-1) then
         ! DAN ROE: Should put a check in here to make sure framename
         !  doesn't get blown up
         write (framename,'(a,a1,i6.6)') &
            reservoirname(1:index(reservoirname," ")-1),".",rremd_idx
         write (6,'(a,a)') "Reservoir Filename= ",trim(framename)
         call amopen(reservoir_unit,framename,'O','F','R')
         ! Read title, # atoms
         read (reservoir_unit,"(a)") line
         write (6,"(a7,a)") "Title= ", line
         ! DAN ROE: should this read have a format?
         read (reservoir_unit,*) i
         if (i.ne.natom) then
            write (6,*) "restart file has wrong #atoms ",framename,i,natom
            backspace (reservoir_unit)
            read (reservoir_unit,"(a)") line
            write (6,"(a)") line
            call mexit(6,1)
         end if
         ! Read coordinates
         read (reservoir_unit,"(6(f12.7))") (x(i),i=1,nr3)
         ! Read Velocities
         if (reserv_velo) read (reservoir_unit,"(6(f12.7))") (v(i),i=1,nr3)
         close (unit=reservoir_unit)
         ! DAN ROE: Should be some error checking on reservoir read.
      else ! Netcdf reservoir
         ! Read coordinates
         if ( NC_error(nf90_get_var(reservoir_ncid, coordVID, x(1:nr3), &
                                    start = (/ 1, 1, rremd_idx /), &
                                    count = (/ 3, natom, 1 /)), &
                       'reading reservoir coordinates') ) call mexit(6,1)
         ! Read velocities
         if ( reserv_velo ) then
            if (NC_error(nf90_get_var(reservoir_ncid, velocityVID, v(1:nr3), &
                                      start = (/ 1, 1, rremd_idx /), &
                                      count = (/ 3, natom, 1 /)), &
                         'reading reservoir velocities')) call mexit(6,1)
         end if
      end if
      ! CARLOS: IMPORTANT! Scale the reservoir velocities!
      !  Scaling is set in subrem_reservoir
      ! If we have const P, need box change. Currently unsupported
      ! DAN ROE: checked for in mdread.
      !if (ifbox >= 1 .and. ntp > 0) then
      !   write (6,*) "const P not allowed for RREMD"
      !   call mexit (6,1)
      !end if
   end if ! master

   ! Now master needs to broadcast coords and velo for reservoir
   !  to all processes in its group.
   call mpi_bcast(x, nr3, mpi_double_precision, 0, commsander, ierror)
   if (reserv_velo) &
      call mpi_bcast(v, nr3, mpi_double_precision, 0, commsander, ierror)
   ! DAN ROE: Debug
   if (master) then
      write (6,'(a38)') "==========End Reservoir Read=========="
   end if

   ! ----===== END RREMD RESERVOIR LOADING =====----

   return

end subroutine load_reservoir_structure


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Get energy of stripped structure in hybrid remd
subroutine hybrid_remd_ene( &
              x,ix,ih,ipairs,qsetup,                        &
              numwatkeep,hybridgb,igb,ntr,nspm,t,temp0,     &
              ntb,cut,ener,do_list_update,nstep,onefac )
   use state
   implicit none
#  include "parallel.h"
#  include "../include/memory.h"

! sander.F90
   _REAL_ x(*)
   integer ix(*), ipairs(*)
   character(len=4) ih(*)
   logical qsetup
! md.h
   integer numwatkeep, hybridgb, igb, ntr, nspm
   _REAL_ t,temp0
! memory.h
!   integer natom, nres
! box.h
   integer ntb
   _REAL_ cut
! runmd.f
   _REAL_ onefac(*)
   type(state_rec) ::  ener
   logical do_list_update
   integer nstep

! nstep, and onefac only needed for printmd

! Temporary storage for natom,
!  ntb, and cut during stripped coord call to force.
! DAN ROE: Is dynamic allocation the way to go?
! NOTE: This allocation currently won't work with amoeba.
   integer hybrid_natom, ref_natom, hybrid_ntb, ier, i
   _REAL_ hybrid_cut

!--------------------

! This is a hybrid REMD run. Get energy of stripped system for next
!  exchange.
! DAN ROE: Note: hybrid code is placed here since it needs access to force()
! 1- First the current coordinates and forces are placed in a temporary
!    array. The original ntb, cut, and natom values are saved.
! 2- Then strip away all but numwatkeep waters from the temporary coords.
!    natom is changed to reflect the new system size.
! 3- Set igb, ntb, and cut to hybrid values, and call force using the
!    temporary coord (now stripped) and force arrays. The force call will use
!    GB based on the hybridgb variable. After we return from
!    force we should have the correct PE for the stripped system.
! DAN ROE: All threads need to do this because all need
!  to call force. All threads know about numwatkeep since it is a
!  namelist variable.

   if (master) write(6,'(17("="),a,i10,17("="))') &
      "HYBRID REMD: energy calc for exch ",mdloop+1

! DAN ROE: Debug
!         if (master) then
!            write (6,*) "Pre-strip coordinates: "
!            write (6,*) (xx(lcrd+i-1),i=1,10)
!            write (6,*) "Pre-force velocities:"
!            write (6,*) (xx(lvel+i-1),i=1,10)
!            write (6,*) "Pre-force forces:"
!            write (6,*) (xx(lforce+i-1),i=1,10)
!            call prntmd(nstep,t,ener,onefac,7,.false.)
!         end if

   ! 1- Store coords and forces in temp arrays,
   !    backup natom, ntb, and cut
   do i=1, natom
      hybrid_coord(3*i - 2) = x(lcrd+3*i - 3)
      hybrid_coord(3*i - 1) = x(lcrd+3*i - 2)
      hybrid_coord(3*i    ) = x(lcrd+3*i - 1)
      hybrid_force(3*i - 2) = x(lforce+3*i - 3)
      hybrid_force(3*i - 1) = x(lforce+3*i - 2)
      hybrid_force(3*i    ) = x(lforce+3*i - 1)
   enddo
   ! Store reference coordinates
   if ( ntr == 1 ) then
      do i=1, natom
         hybrid_refc(3*i - 2) = x(lcrdr+3*i - 3)
         hybrid_refc(3*i - 1) = x(lcrdr+3*i - 2)
         hybrid_refc(3*i    ) = x(lcrdr+3*i - 1)
      enddo
   end if
   hybrid_natom=natom
   hybrid_ntb=ntb
   hybrid_cut=cut

   ! 2- Strip waters. This will change the coords (hybrid_coord) and
   !    # of atoms (natom). The coords will be imaged.
   if (master) write(6,'(a)') "HYBRID REMD: Stripping waters"
   call stripwat(natom,hybrid_coord,ix(i02),ih(m02),x(lmass), &
                 ix(i70),nspm,nres,numwatkeep)
   ! Strip reference coordinates. Strip actual array since it is
   ! not passed into force. Will be restored from hybrid_refc after.
   if ( ntr == 1 ) then
      ref_natom = hybrid_natom
      call stripwat(ref_natom,x(lcrdr),ix(i02),ih(m02),x(lmass), &
                    ix(i70),nspm,nres,numwatkeep)
   end if

   if (master) then
      write(6,'(a,i8)') "HYBRID REMD: New natom= ",natom
      if (ntr == 1) &
         write(6,'(a,i8)') "HYBRID REMD: New ref natom= ",ref_natom
      if (hybridwritetraj) &
         call corpac(hybrid_coord,1,natom*3,remstripcoord_unit,.true.)
   end if
   call mpi_barrier(commworld, ier)

   ! 3- Call force to calculate energy using the stripped
   !    coordinates and GB (based on hybridgb).
   ! Make sure PBC off, and change the cutoff
   igb=hybridgb
   ntb=0
   cut=9801.0d0 !  = 99 * 99
   if (master) write(6,'(a)') "HYBRID REMD: Calling force."
   call force(x,ix,ih,ipairs,hybrid_coord,hybrid_force,ener,ener%vir, &
              x(l96),x(l97),x(l98),x(l99),qsetup, &
              do_list_update,nstep)
! DAN ROE: Debug
   if (master) then
!            write (6,*) "Post-force coordinates: "
!            write (6,*) (xx(lcrd+i-1),i=1,10)
!            write (6,*) "Post-force velocities:"
!            write (6,*) (xx(lvel+i-1),i=1,10)
!            write (6,*) "Post-force forces:"
!            write (6,*) (xx(lforce+i-1),i=1,10)
      call prntmd(nstep,t,ener,onefac,7,.false.)
   end if
   if (master) write(6,'(a,f13.4,a,f6.2)') &
      "HYBRID REMD: myEptot= ",ener%pot%tot," myTargetTemp= ",temp0

   ! 4- Restore original natom, igb, cut
   if (master) write(6,'(a)') "HYBRID REMD: Restoring..."
   natom=hybrid_natom
   igb=0
   ntb=hybrid_ntb
   cut=hybrid_cut
   ! Restore reference coordinates
   if ( ntr == 1 ) then
      do i=1, natom
         x(lcrdr+3*i - 3) = hybrid_refc(3*i - 2)
         x(lcrdr+3*i - 2) = hybrid_refc(3*i - 1)
         x(lcrdr+3*i - 1) = hybrid_refc(3*i    )
      enddo
   end if
   if (master) write(6,'(25("="),a,25("="))') &
      "END HYBRID REMD energy calc."

   return

end subroutine hybrid_remd_ene

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ strip water from structure for replica exchange with mixed solvent models
subroutine stripwat(nrp,x,ipres,lbres,amass,nsp,nspm,nres,numwatkeep)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!=====================================================================
! GMS
! ---
! This routine requires that the residues are organized so that waters
! come after any other residue, and there should be not residues other
! than water in the middle.
!
! Before calling this routine, make sure to have a copy of the original
! number of residues, coordinate and force arrays, since they will be
! overwritten here and will need to be restored later.
!
! It returns:
!    nrp --> New number of atoms (stripped system)
!    x   --> New coordinates (solute + closest waters)
! The way it works:
!
! 1. Find where the waters start, puts that residue # into 'firstwat'.
!    If there is any residue other than water after this point,
!    it nbombs the calculation.
!
! 2. Reimage the system so that waters surround the solute. Make sure
!    to image the waters correctly to the minimum distance.
!    (Remember that the hybrid remd code uses GB, which is not periodic.)
!
! 3. For every water, calculate the distance from it's oxygen to EVERY
!    non-water atom in the system, keeping only the shortest one.
!
!    * After this point, we have a list of all the waters, their
!      residue numbers and the distance to the closest non-water atom.
!
! 4. Sort the waters by distance.
!
!    * At this point, the routine uses a a "selection sort" algorithm,
!      which scales as N**2, N being the number of waters. If this
!      becomes a bottleneck, it can be changed later to "mergesort"
!      or a "binary tree" sort, which are N*logN.
!
! 5. Set a new logical array "keepwat(size)", which is a logical mask
!    indicating whether to keep this water or not. The first "numkeepwat"
!    waters are set to .true., all the rest is false.
!
!    * After this sort we have three arrays ordered by increasing distance
!      to the closest non-solute atom:
!
!      watresnum(size) --> The residue number of the water
!      closedist(size) --> The distance to the closest non-water residue
!      keepwat(size)   --> Logical. ".true." if we want to keep this water.
!
! 6. Now, sort those 3 arrays in the order by residue number instead
!
! 7. Get the coordinates for the atoms in the solute plus *only* the
!    waters we want to keep. That means looping through all residues and
!    copying the coordinates only of the non-water atoms, plus the atoms
!    from water residues marked with "keepwat = .true."
!
! 8. Copy those coordinates to the main coordinate array, update the "nrp"
!    variable with the number of atoms in this reduced system.
!
! 9. Return.
!=====================================================================

   implicit none
#  include "parallel.h"
! box.h needed for imaging water
#  include "box.h"

   integer, parameter :: size=10000

   integer nrp, ipres(*), nsp(*), nspm, nres, numwatkeep
   _REAL_  x(*), amass(*)
   character(len=4) lbres(*)

   _REAL_  x2(size), closedist(size), xcm(3)
   _REAL_  r2, xij, yij, zij, xi, yi, zi, tempr, aamass, tmassinv

   integer watresnum(size)
   integer firstwat, watpointer, newnrp, totwat, tempi
   integer i, i1, j, j1, k

   logical keepwat(size), templ

!----------------------------------------------------------
   ! find out where water starts
   firstwat = -1

      do i=1,nres
         if (lbres(i) == "WAT") then
            firstwat=i
            exit
         end if
      enddo

      REQUIRE( firstwat > 0 )

      ! next we need to reimage the system so that waters surround the solute
      ! we could use minimum image type calculation to get the
      ! closest distance but we still need the water imaged properly
      ! for the GB calculation (which is not periodic)
      ! follow the code for imaging from runmd's iwrap=1 code
      ! first center the system on the CM of the solute

      xcm(1) = 0.d0
      xcm(2) = 0.d0
      xcm(3) = 0.d0

      ! here tmassinv is only for non-water

      tmassinv=0.d0
      i = 0
      do k=1,firstwat-1
         do j =ipres(k),ipres(k+1)-1
            aamass = amass(j)
            xcm(1) = xcm(1) + x(i+1)*aamass
            xcm(2) = xcm(2) + x(i+2)*aamass
            xcm(3) = xcm(3) + x(i+3)*aamass
            i = i + 3
            tmassinv=tmassinv+aamass
         enddo
      end do

      tmassinv=1.d0/tmassinv

      xcm(1) = xcm(1) * tmassinv
      xcm(2) = xcm(2) * tmassinv
      xcm(3) = xcm(3) * tmassinv

      ! center all atoms, not just solute

      do i=1,nrp
         x(3*i-2) = x(3*i-2)-xcm(1)
         x(3*i-1) = x(3*i-1)-xcm(2)
         x(3*i)   = x(3*i)-xcm(3)
      enddo

      ! now re-image the box

      call wrap_molecules(nspm,nsp,x)
      if(ifbox == 2) call wrap_to(nspm,nsp,x,box)

      ! now start setting closest distance between water and solute

      ! in fact we do not need the distance, the distance squared is
      ! fine since we only want to sort them

      watpointer=0

      do i= firstwat,nres

         ! CHECK TO MAKE SURE IT IS WATER

         if (lbres(i).ne."WAT ") then
            if (master) then
               write (6,*) "solvent molecule is not water: ",lbres(i)
               write (6,*) "stopping water search"
            end if
            call mexit(6,1)
         end if

         watpointer=watpointer+1

         ! closedist(i) is the distance from the water to the
         ! closest solute atom

         closedist(watpointer)=9999
         watresnum(watpointer)=i

         !  do i1=ipres(i),ipres(i+1)-1

         ! for water, just take distance to oxygen (first atom in water)

         i1=ipres(i)
         xi = x(3*i1-2)
         yi = x(3*i1-1)
         zi = x(3*i1)

         ! loop over non-water residues/atoms

         do j=1,firstwat-1
            do j1=ipres(j),ipres(j+1)-1

               xij = xi - x(3*j1-2)
               yij = yi - x(3*j1-1)
               zij = zi - x(3*j1  )
               r2 = xij*xij + yij*yij + zij*zij

               if (r2.lt.closedist(watpointer)) closedist(watpointer)=r2
            enddo
         enddo
      enddo

      totwat=watpointer

      ! now we have a list of all waters - their residue number and distance
      ! to solute

      ! sort them by distance

      do i=1,totwat-1
         do j=i+1,totwat

           ! write (6,*) "working on ",i,j

            if (closedist(i).gt.closedist(j)) then
               tempr=closedist(i)
               tempi=watresnum(i)
               closedist(i)=closedist(j)
               watresnum(i)=watresnum(j)
               closedist(j)=tempr
               watresnum(j)=tempi

            end if
         enddo
      enddo

      ! now set save flags for closest numwatkeep

      do i=1,numwatkeep
         keepwat(i)=.true.
      enddo
      do i=numwatkeep+1,totwat
         keepwat(i)=.false.
      enddo

      ! now sort them back into the order by watresnum

      do i=firstwat,nres

         ! i is the residue number we are looking for

         ! i1 is the current water at this residue number
         ! this is in the 1 to totwat sense so we can pull those
         ! indices for swapping waters in list

         i1=i-firstwat+1

         ! look to see where this water is, and put it in place

         do j=1,totwat

            if (watresnum(j) == i) then

               ! found it, so swap them

               tempr=closedist(i1)
               tempi=watresnum(i1)
               templ=keepwat(i1)
               closedist(i1)=closedist(j)
               watresnum(i1)=watresnum(j)
               keepwat(i1)=keepwat(j)
               closedist(j)=tempr
               watresnum(j)=tempi
               keepwat(j)=templ

               ! get next i

               exit
            end if
         enddo
      enddo

      ! now go through and write the restart file. we need to write to temp
      ! array since we can't remove atoms while writing the file itself since
      ! there is more than 1 molecule per line

      newnrp=0
      do i=1,nres
         if (i.lt.firstwat.or.keepwat(i-firstwat+1)) then
            do j=ipres(i),ipres(i+1)-1
               newnrp=newnrp+1
               x2(3*newnrp-2) = x(3*j-2)
               x2(3*newnrp-1) = x(3*j-1)
               x2(3*newnrp)   = x(3*j)
            enddo
         end if
      enddo

      ! now copy this array to old one, reset nrp

      nrp=newnrp
      do i=1,nrp
        x(3*i-2) = x2(3*i-2)
        x(3*i-1) = x2(3*i-1)
        x(3*i)   = x2(3*i)
      enddo

   return

end subroutine stripwat

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ calculate cluster of current coordinates based on dihedral binning
subroutine calc_rremd_cluster(x)

!*********************************************************************
!               SUBROUTINE CALC_RREMD_CLUSTER
!*********************************************************************
! Used with rremd==3. Dihedral calc. based on Jcoupling routine from nmr.f
! Called from runmd(). Sets incluster. Only masters call this subroutine.

   use constants, only : pi
   implicit none

! Input:
!    X(I)  : Coordinate array.
   _REAL_, dimension(*), intent(in) :: x

! Internal variables
   integer i, j
!  I1,I2,I3: Atom pointers for this angle (3*(I-1), where I is absolute
!        I4  atom number.
   integer i1, i2, i3, i4
   _REAL_ Lx, Ly, Lz, Rx, Ry, Rz, Sx, Sy, Sz
   _REAL_ Lnorm, Rnorm, dih_value

!--------------------

   do i=1,nclustdih
      ! ----- Dihedral Calculation -----
      i1=3*(dihclustat(1,i)-1)
      i2=3*(dihclustat(2,i)-1)
      i3=3*(dihclustat(3,i)-1)
      i4=3*(dihclustat(4,i)-1)
      !write (6,*) "calculating dih for ",i1,i2,i3,i4

      Lx = ((x(i2+2)-x(i1+2))*(x(i3+3)-x(i2+3))) - ((x(i2+3)-x(i1+3))*(x(i3+2)-x(i2+2)))
      Ly = ((x(i2+3)-x(i1+3))*(x(i3+1)-x(i2+1))) - ((x(i2+1)-x(i1+1))*(x(i3+3)-x(i2+3)))
      Lz = ((x(i2+1)-x(i1+1))*(x(i3+2)-x(i2+2))) - ((x(i2+2)-x(i1+2))*(x(i3+1)-x(i2+1)))

      Rx = ((x(i4+2)-x(i3+2))*(x(i2+3)-x(i3+3))) - ((x(i4+3)-x(i3+3))*(x(i2+2)-x(i3+2)))
      Ry = ((x(i4+3)-x(i3+3))*(x(i2+1)-x(i3+1))) - ((x(i4+1)-x(i3+1))*(x(i2+3)-x(i3+3)))
      Rz = ((x(i4+1)-x(i3+1))*(x(i2+2)-x(i3+2))) - ((x(i4+2)-x(i3+2))*(x(i2+1)-x(i3+1)))

      Lnorm = sqrt(Lx*Lx + Ly*Ly + Lz*Lz)
      Rnorm = sqrt(Rx*Rx + Ry*Ry + Rz*Rz)

      Sx = (Ly*Rz) - (Lz*Ry)
      Sy = (Lz*Rx) - (Lx*Rz)
      Sz = (Lx*Ry) - (Ly*Rx)

      dih_value = (Lx*Rx + Ly*Ry + Lz*Rz) / (Lnorm * Rnorm)

      if ( dih_value > 1.d0 ) dih_value = 1.d0
      if ( dih_value < -1.d0 ) dih_value = -1.d0

      dih_value = acos( dih_value );

      if ( (Sx * (x(i3+1)-x(i2+1)) + Sy * (x(i3+2)-x(i2+2)) &
          + Sz * (x(i3+3)-x(i2+3))) < 0 ) &
         dih_value = -dih_value;

      ! now convert to degrees
      dih_value = dih_value * 180.d0/pi
      ! ----- End Dihedral Calculation ----

#ifdef RREMD_DEBUG
      write (6,'(a,f8.3)') "| RREMD_DEBUG: Torsion= ", dih_value
#endif
      ! dihedral is in -180 to 180 range. Subtract min value.
      dih_value = dih_value - dihclustmin(i)
      ! wrap values less than 0
      if (dih_value < 0) dih_value = dih_value + 360.d0
#ifdef RREMD_DEBUG
      write (6, '(a,f8.3)') "|            : ShiftedTorsion= ", dih_value
#endif
      ! calculate bin
      dih_value = dih_value / dihcluststep(i)
      currdihid(i) = INT( dih_value )
#ifdef RREMD_DEBUG
      write (6, '(a,i8)'  ) "|            : Bin= ", currdihid(i)
#endif
   enddo ! 1, nclustdih

   ! Set the default cluster to #0, which means "no match". clustersize(0) was
   !  set to 0 in multisander.F90. This means that "no match" structures will
   !  automatically be rejected.
   ! Eventually the reservoir could have "no match" structures added to it.
   incluster=0
   do i=1,nclust
      incluster = i
      do j=1,nclustdih
         ! As soon as a bin doesn't match, exit
         if (clusterid(i,j).ne.currdihid(j)) then
           incluster = 0
           exit ! exit j loop
         end if
      enddo
      if (incluster .ne. 0) then
        exit ! Exit i loop as soon as a cluster matches
      end if
   enddo
   write (6,'(a,i10)') &
      "RREMD: Current MD structure was assigned to cluster ",incluster

   return

end subroutine calc_rremd_cluster

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ load energy/cluster information for doing reservoir REMD
subroutine load_reservoir_files()
   use AmberNetcdf_mod, only: NC_checkTraj, NC_openRead, NC_setupReservoir,&
                              NC_readReservoir
   ! All threads call this subroutine
   implicit none

#  include "parallel.h"

   integer numatoms, iseed, ierror, i, nodeid, holder
! i1-4 for reading dihedral atom nums
   integer j, i1, i2, i3, i4
   _REAL_ minIn
   integer eptotVID, binsVID

! --------------------
   if (master) then
      ! Check if reservoirname is a Netcdf reservoir
      if (NC_checkTraj(reservoirname)) then
         if (NC_openRead(reservoirname, reservoir_ncid)) call mexit(6,1)
         if (NC_setupReservoir(reservoir_ncid, reservoirsize, restemp0, numatoms,&
                               coordVID, velocityVID, eptotVID, binsVID, iseed)) &
           call mexit(6,1)
         ! bins are required for rremd 3
         if (rremd.eq.3 .and. binsVID.eq.-1) then
            write(0,'(a)') 'FATAL: rremd=3 and netcdf reservoir does not contain cluster nums.'
            call mexit(6,1)
         end if
         reserv_velo = (velocityVID.ne.-1)
      else
         reservoir_ncid = -1
         ! Open saveene, which contains energies for each reservoir struct
         call amopen(remin_unit,saveenefile,'O','F','R')
         ! Read saveene file header line.
         ! Format:
         ! <# reservoir structures> <reservoir T> <#atoms> <random seed>
         !   <velocity flag [1 if reading velocity from reservoir]>
         ! Read # atoms - we dont seem to know natom yet.
         ! DAN ROE: Check if we need #atoms read here.
         ! Note: For cleanliness we could make the unit here the remlog
         !  unit since remlog is not used yet.
         read (remin_unit,*,err=5000) reservoirsize,restemp0,numatoms,iseed,holder
         reserv_velo = holder == 1
      end if
      ! Write the reservoir type and other info to remtype
      if (master_master) then
         write (remtype_unit,'(a,a)') &
            "RREMD: Info from saveene file ",saveenefile
         write (remtype_unit,'(a,i5)')   "  NumAtoms= ",numatoms
         write (remtype_unit,'(a,i5)')   "  ReservoirSize= ",reservoirsize
         write (remtype_unit,'(a,f6.2)') "  ReservoirTemp(K)= ",restemp0
         write (remtype_unit,'(a,i10)')  "  ReservoirRandomSeed= ",iseed
         if (reserv_velo) then
            write (remtype_unit,'(a)') "  Velocities will be read from reservoir"
         else
            write (remtype_unit,'(a)') "  Velocities will be assigned to structure after exchange"
         end if
        ! Print reservoir type information
         write (remtype_unit,'(a,i5)') "RREMD: Reservoir type ",rremd
         if (rremd == 1) then
            write (remtype_unit,'(a)') &
               "  Boltzmann weighted reservoir,&
               & exchange uses delta beta"
         elseif (rremd == 2) then
            write (remtype_unit,'(a)') &
               "  Non-Boltzmann 1/N weighted reservoir,&
               & exchange uses beta"
         elseif (rremd == 3) then
            write (remtype_unit,'(a)') &
               "  Non-Boltzmann weighted reservoir with&
               & defined weights"
            write (remtype_unit,'(a)') &
               "  (Currently only works via dihedral clustering)"
            write (remtype_unit,'(a)') &
               "  Exchange uses beta and weights."
            write (remtype_unit,'(a)') &
               "  FORCING WEIGHT OF 0 FOR STRUCTURES&
               & NOT IN RESERVOIR!"
         else
            write (remtype_unit,'(a,i5)') &
               "Unknown reservoir type: rremd=", rremd
            !call mexit(6,1)
         end if
      end if ! master_master remtype file write

      ! Exit if unknown reservoir type
      ! DAN ROE: Should this exit just be called above?
      if (rremd.lt.0 .or. rremd.gt.3) call mexit(6,1)

      ! Check reservoir size limits - non-netcdf reservoir only
      ! Currently the format for reservoir files allows a max integer
      ! width of 6. To accommodate a larger reservoir the format string
      ! in load_reservoir_structure would have to be changed.
      if (reservoir_ncid.eq.-1 .and. reservoirsize.gt.maxreservoir) then
         write (6,'(a,i6)') "RREMD ERROR: Reservoir size limit is currently ",&
                            maxreservoir
         write (6,'(a)') "To accommodate larger reservoir sizes edit the reservoir&
                         & file format string in load_reservoir_structure (remd.F90)"
         call mexit(6,1)
      end if

      ! Allocate memory for reservoir structure energies
      allocate( saveene(reservoirsize), stat=ierror)
      REQUIRE(ierror==0)
      ! Allocate memory to hold cluster numbers from DihedralCluster
      if (rremd == 3) then
        allocate( clusternum(reservoirsize), stat=ierror)
        REQUIRE(ierror==0)
      end if

      ! Read energy (and cluster number for rremd==3) for each structure
      ! All we really need to load are the energies; if the exchange
      !  is successful we can grab the coords (and velocity if necessary
      !  from the disk during the run.
      ! DAN ROE: For rremd<3, If cluster #s are present this will break.
      !  Should add a check. Maybe make cluster #s a separate file?
      if (reservoir_ncid.eq.-1) then
         do i=1,reservoirsize
            if (rremd==1.or.rremd==2) then
               ! 1 = Boltzmann weighted, read energies
               ! 2 = 1/N reservoir, read energies
               read(remin_unit,*,err=5000,iostat=ierror) saveene(i)
            elseif (rremd==3) then
               ! Weighted reservoir, read energies and corresponding cluster #s
               read(remin_unit,*,err=5000,iostat=ierror) saveene(i), clusternum(i)
               ! DAN ROE: Since I started numbering clusters at 0, add 1
               !  Change later to be consistent.
               !clusternum(i)=clusternum(i)+1
               if (clusternum(i)<1) then
                  write(6,'(a)') "Error: Cluster# < 1 not allowed."
                  goto 5000
               end if
            end if
            if (ierror<0.and.i<reservoirsize) then
               write(6,'(a)') "Error: EOF reached before all values read."
               goto 5000
            end if
         enddo
         close (unit=remin_unit)
      else ! Netcdf reservoir read
         if (NC_readReservoir(reservoir_ncid, reservoirsize, eptotVID, &
                              binsVID, saveene, clusternum))&
            call mexit(6,1)
      end if

      ! Read in the cluster info file for rremd==3
      ! DAN ROE: Put in error checking on reads
      if (rremd==3) then
         call amopen(remin_unit,clusterinfofile,'O','F','R')
         ! First read number of dihedrals to bin
         read (remin_unit,*,err=5020) nclustdih
         ! Allocate space for dihedral atoms, bins, min, and step
         allocate( dihclustat(4, nclustdih), &
                   dihclustnbin(nclustdih), &
                   currdihid(nclustdih), &
                   dihclustmin(nclustdih), &
                   dihcluststep(nclustdih), &
                   stat=ierror )
         REQUIRE(ierror==0)
         ! Now read atom #s and bins for each dihedral angle.
         ! Atom #s should start from 1
         do i=1,nclustdih
            ! Format: atom#1 atom#2 atom#3 atom#4 bins min
            read (remin_unit,"(5(i10,1x),f8.3)",err=5020) i1,i2,i3,i4,j,minIn
            dihclustat(1,i)=i1
            dihclustat(2,i)=i2
            dihclustat(3,i)=i3
            dihclustat(4,i)=i4
            dihclustnbin(i)=j
            dihclustmin(i)=minIn
            dihcluststep(i)=dble(j) ! Convert # of bins to float
            dihcluststep(i)=360.d0 / dihcluststep(i)
#ifdef RREMD_DEBUG
            write(6,'(5(i10,1x),2(f8.3,1x))') &
                  i1,i2,i3,i4,j,minIn,dihcluststep(i)
#endif
         enddo
         ! Read number of clusters
         read (remin_unit,*,err=5020) nclust
#ifdef RREMD_DEBUG
         write(6,'(a,i10)') "Nclust= ", nclust
#endif
         ! Allocate space for cluster IDs and sizes
         allocate( clusterid(nclust, nclustdih), &
                   clustersize(0:nclust), &
                   stat=ierror)
         REQUIRE(ierror==0)

         ! Read cluster weight and ID (bin values) for each cluster
         ! Set the cluster size for cluster 0, which corresponds to an MD
         !  structure not being present in the reservoir. The current code
         !  does not add this to a new cluster. NOTE: clustersize(0) is
         !  possible since clustersize is dimensioned from 0 above.
         ! This makes the reservoir exchange rigorous - structures not in
         !  the reservoir have a weight of zero, and therefore can't be
         !  exchanged.
         clustersize(0)=0
         do i=1,nclust
            ! Format: Cluster# Weight Bin1 ... BinNclustdih
            ! DAN ROE: This assumes clusters are in order!
            read (remin_unit,*,err=5020) &
                  j, clustersize(i), (clusterid(i,j),j=1,nclustdih)
#ifdef RREMD_DEBUG
            write(6,'(a,2(i10))') "Cluster# and size: ", i, clustersize(i)
#endif
         enddo

         close (unit=remin_unit)

         ! Overall master Write out information to remtype
         if (master_master) then
            write(remtype_unit,'(a,a)') "RREMD: clusterinfo file ",clusterinfofile
            write (remtype_unit,'(a,i10)') "  NumDihedrals= ", nclustdih
            do i=1,nclustdih
               write (remtype_unit,"(a14,i5,a7,4(1x,i5),a7,i5,2(a7,f8.3))") &
                  "    Dihedral #",i," atoms:",dihclustat(1,i), &
                  dihclustat(2,i),dihclustat(3,i),dihclustat(4,i), &
                  " bins: ",dihclustnbin(i)," step: ",dihcluststep(i), &
                  "  min: ",dihclustmin(i)
            enddo
            write (remtype_unit,'(a,i10)') "  NumClusters= ",nclust
            do i=1,nclust
               write (remtype_unit,'(a,2(i10,a))') &
                  "    Cluster ",i," has ",clustersize(i)," members"
               write (remtype_unit,*) "      Bins= ", &
                      (clusterid(i,j),j=1,nclustdih)
            enddo
         end if
         ! dihclustnbin no longer used, all necessary info is in dihcluststep
         if (allocated(dihclustnbin)) deallocate(dihclustnbin, stat=ierror)
         REQUIRE(ierror==0)
      end if ! rremd == 3
   end if ! master

   ! Set random # generator based on the seed in saveene; this is so
   !  we can do independent runs with different sequence of structures
   !  from the reservoir. Here we have no concept of the ig seed in
   !  md.in and we don't want to hardcode, so saveene is a convenient
   !  location. All threads do this.
   ! DAN ROE: FileDebug
   ! Send iseed to all threads since everyone needs to call amrset
   ! NOTE: repnum = nodeid + 1
   nodeid = repnum - 1
   call mpi_bcast(iseed,1,mpi_integer,0,commworld,ierror)

   call amrset(iseed + 17 * nodeid)
   ! Broadcast reserv_velo to all so we know it in sander() when we
   !  read reservoir.
   ! DAN ROE: Does every process need this or just masters?
   call mpi_bcast(reserv_velo,1,mpi_logical,0,commworld,ierror)

   return

5000  write (6,*) "RREMD: Error in reading saveene!"
      call mexit(6,1)

5020  write (6,*) "RREMD: Error in reading clusterinfo!"
      call mexit(6,1)

end subroutine load_reservoir_files

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Calculate pressure/volume correction to replica exchange delta.
function pv_correction(ourtemp, nbrtemp, neighbor, comm_rep_master)
  ! USE STATEMENTS
  use constants, only : KB
  ! ARGUMENTS
  implicit none
  _REAL_ pv_correction
  _REAL_, intent(in)  :: ourtemp, nbrtemp
  integer, intent(in) :: neighbor
  integer, intent(in) :: comm_rep_master
  ! INCLUDES
  ! LOCAL VARIABLES
  _REAL_, dimension(2) :: our, nbr
  integer ierr, istat(mpi_status_size)

  ! Pressure in position 1, volume in position 2.
  our(1) = my_pressure
  our(2) = my_volume
  call mpi_sendrecv(our, 2, mpi_double_precision, neighbor, 52, &
                    nbr, 2, mpi_double_precision, neighbor, 52, &
                    comm_rep_master, istat, ierr)
# ifdef VERBOSE_REMD
  write(6,'(3(a,f16.8))') '| REMD: OurPressure= ', our(1), ' atm, OurVolume= ', our(2), &
          ' Ang^3, OurTemp= ', ourtemp
  write(6,'(3(a,f16.8))') '| REMD: NbrPressure= ', nbr(1), ' atm, NbrVolume= ', nbr(2), &
          ' Ang^3, NbrTemp= ', nbrtemp
# endif
  pv_correction = (((1.d0/(KB*nbrtemp))*nbr(1)) - ((1.d0/(KB*ourtemp))*our(1))) &
                  * (nbr(2) - our(2))
  ! Need to convert from units of atm*Ang^3 to kcal/mol
  ! 1 atm*Ang^3 * 101325 Pa/atm * E-30 m^3/Ang^3 * J/Pa*m^3 * kcal/4184 J * 6.02E23/mol
  pv_correction = pv_correction * 1.4584E-5
# ifdef VERBOSE_REMD
  write(6,'(a,E16.8)') '| REMD: PvCorrection= ', pv_correction
# endif
end function pv_correction

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Broadcast new unit cell information in box1 and update unit cell.
subroutine remd_bcast_cell(box0, box1, commsander)
  ! USE STATEMENTS
  !use nblist, only : alpha, beta, gamma
  use nblist, only : fill_tranvec
  ! ARGUMENTS
  implicit none
  _REAL_, intent(in), dimension(:) :: box0       ! Original box
  _REAL_, intent(in), dimension(:) :: box1       ! New box
  integer, intent(in)              :: commsander ! COMM for broadcast
  ! INCLUDES
  ! LOCAL VARS
  integer ierr
  _REAL_ rmu(3)

  call mpi_bcast(box1, 3, mpi_double_precision, 0, commsander, ierr)
  !call fill_ucell(box(1), box(2), box(3), alpha, beta, gamma)
  ! Calculate scale from old box to new
  rmu(1:3) = box1(1:3) / box0(1:3)
  !write(6, '(a,6f8.3,3(x,E10.4))') 'DBG: Box {old, new, scale} = ', box0(1:3), box1(1:3), rmu(1:3)
  call redo_ucell( rmu )
  call fill_tranvec()
  ! TODO determine if grid sizes need to change?
end subroutine remd_bcast_cell

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Perform actions necessary to exchange replicas in REMD
subroutine hremd_exchange(rem_dim, x, ix, ih, ipairs, qsetup, do_list_update)

   use constants, only : TWO
   use state
   implicit none
  ! INCLUDES
#  include "box.h"
#  include "parallel.h"
#  include "../include/memory.h"

   ! sander.F90
   _REAL_  x(*)
   integer ix(*), ipairs(*)
   character(len=4) ih(*)
   logical :: qsetup
   logical :: do_list_update
   integer, intent(in) :: rem_dim

   ! local variables
   integer ierror, i, istatus(mpi_status_size)

   type(state_rec) :: expe

   _REAL_ delta, straw, metrop
   _REAL_ l_eptot, r_eptot, o_eptot
   _REAL_ oex_eptot, myexeptot
   _REAL_ l_temp0, r_temp0, o_temp0
   _REAL_ alleptot(numreps), allexeptot(numreps)
   _REAL_ alltempi(numreps), alltemp0(numreps)
   _REAL_ all_l_fe(numreps), all_r_fe(numreps)
   _REAL_ l_fe, r_fe
   _REAL_ original_box(3)
   _REAL_ new_box(3)
   _REAL_ pvterm

   integer rep, neighbor_rep ! counter/holder for remlog printing data
   character success_char ! character holder

   ! repnum is the actual replica number
   ! in h-rem, repnum is index
   integer myrepnum, l_repnum, r_repnum, o_repnum
   integer r_ex, l_ex ! right exchange and left exchange

   ! temporary for exchange success, will allgather from all nodes then
   ! add sum to actual exchsuccess array. also define an array for allgather sum
   integer texchsuccess(numreps),tempsuccess(numreps)
   integer all_l_ex(numreps), all_r_ex(numreps)
   _REAL_ exchfrac(numreps)

   logical exchange

! ---------------------

   ! Initialize variables
   texchsuccess(:) = 0
   all_r_fe(:) = 0.d0
   all_l_fe(:) = 0.d0
   all_l_ex(:) = 0
   all_r_ex(:) = 0
   l_fe = 0.0d0
   r_fe = 0.0d0
   l_ex = 0
   r_ex = 0
   myrepnum = 0
   l_repnum = 0
   r_repnum = 0
   o_temp0 = 0.d0
   o_eptot = 0.d0

   delta = 0.0d0
   straw = 0.0d0
   metrop= 0.0d0
   pvterm= 0.0d0

   if (master) then
      ! Set the partner array. NOTE: also sets index_list
      call set_partners(rem_dim, remd_size)
      myrepnum = replica_indexes(rem_dim)
      l_repnum = partners(1)
      r_repnum = partners(2)
   end if

   ! ---===  h-rem-fep ===---
   ! 1st: collect current PE, find out exchange neighbors
   ! 2nd: send and receive structures between neighbor pairs, find out new PE
   !      & we don't want to rescale velocities, we want to keep velocity
   !      & unchanged before and after exchange attempts
   ! 3rd: find out FE, determine whether exchange is accepted or not
   !      Only the sander masters do the exchanges.
   ! 4th: if exchange succeed, exchange coordinates and bcast
   !      & to non sander-masters.


   ! Step 1: collect current PE, find out exchange neighbors
   !         seems can be done by all procs.
   !
   ! first, bcast myrepnum from sandermaster to non-sandermasters.
   ! Only sandermasters know remd_rank, myrepnums are
   ! based on remd_ranks. So bcast myrepnum to non-sandermasters

   call mpi_bcast(myrepnum, 1, mpi_integer, 0, commsander, ierror)

   if (master) then
#ifdef VERBOSE_REMD
      write(6, '(a,24("="),a,24("="))') '| ',' H-REMD EXCHANGE CALCULATION '
      write(6, '(a,i10,a,i1)') '| Exch= ', mdloop, ' RREMD= ', 0
#endif
      ! Gather current PE, alleptot are indexed by replica number
      ! Gather current temp0, alltemp0 are indexed by replica number
      ! Gather current tempi, alltempi are indexed by replica number
      call mpi_allgather(my_remd_data%myEptot, 1, mpi_double_precision, &
                         alleptot, 1, mpi_double_precision, &
                         remd_comm, ierror)
      call mpi_allgather(my_remd_data%mytargettemp, 1, mpi_double_precision, &
                         alltemp0, 1, mpi_double_precision, &
                         remd_comm, ierror)
      call mpi_allgather(my_remd_data%mytemp, 1, mpi_double_precision, &
                         alltempi, 1, mpi_double_precision, &
                         remd_comm, ierror)

      ! Get the energies and temp0 of these neighbor replicas
      l_eptot = alleptot(l_repnum)
      r_eptot = alleptot(r_repnum)
      l_temp0 = alltemp0(l_repnum)
      r_temp0 = alltemp0(r_repnum)

      ! Set up partner replica information.
      ! By definition, only even replica #s will initiate the exchange.
      if(mod(myrepnum, 2) == 0) then
         ! this replica will calculate delta
         if(jumpright(rem_dim)) then
            ! Partner is to the right, i.e. larger in repnum
            o_repnum = r_repnum
            o_eptot  = r_eptot
            o_temp0  = r_temp0
         else
            ! Partner is to the left, i.e. smaller in repnum
            o_repnum = l_repnum
            o_eptot  = l_eptot
            o_temp0  = l_temp0
         end if
      else
         ! This replica will not calculate delta, right and left are
         ! switched compared to controlling replica.
         if(jumpright(rem_dim)) then
            o_repnum = l_repnum
            o_eptot  = l_eptot
            o_temp0  = l_temp0
         else
            o_repnum = r_repnum
            o_eptot  = r_eptot
            o_temp0  = r_temp0
         end if
      end if  ! (even/odd rank check)
   end if ! master in step 1

   ! Step 1 is finished here so far, sander masters know myrepnum and o_repnum
   ! however, sander non-masters know only myrepnum

   ! =====================================================================
   ! Step 2: send and receive coordinates;
   !         & calculate and collect new PE

   ! =====================================================
   ! Start exchanging coordinates and calculate PE
   ! unlike step 1, step 2 will be done by all processes
   ! =====================================================

   if (master) then
      !scaling factor for the velocities in the case of accepted exchange
      my_remd_data%myscaling = sqrt(my_remd_data%mytargettemp / o_temp0)
      ! Exchange of coordinates can be done with single sendrecv call. Store
      ! other replica's coordinate array, x(lcrd), in temporary coordinate
      ! array, xtemp
      call mpi_sendrecv(x(lcrd), 3*natom, mpi_double_precision, o_repnum-1, &
                        10, xtemp, 3*natom, mpi_double_precision, o_repnum-1, &
                        10, remd_comm, istatus, ierror)
      ! Exchange unit cells if necessary TODO are there other things that need exchanging
      if (exchange_ucell) then
        original_box(1:3) = box(1:3)
        call mpi_sendrecv(box,     3, mpi_double_precision, o_repnum-1, 11, &
                          new_box, 3, mpi_double_precision, o_repnum-1, 11, &
                          remd_comm, istatus, ierror)
        if (use_pv) pvterm = pv_correction(my_remd_data%mytargettemp, o_temp0, o_repnum-1, remd_comm)
      endif
   end if

   ! Broadcast our new temporary coordinate array
   call mpi_bcast(xtemp, 3*natom, mpi_double_precision, 0, commsander, ierror)
   ! Update the unit cell if necessary
   if (exchange_ucell) call remd_bcast_cell(box, new_box, commsander)
   !write(6,'(a,3f8.3)') 'DBG: New box: ', box(1:3)
   ! Call force on temporary coordinate array to get new energy
   call force(x, ix, ih, ipairs, xtemp, ftemp, expe, expe%vir, x(l96), &
              x(l97), x(l98), x(l99), qsetup, do_list_update,-1)

   ! =======================================================
   ! Step 3: calculate free energy difference, delta
   !         find out whether exchange or not
   !         This step is done by master
   !
   ! Step 4: also done by sander masters. so it is
   !         slightly easier to append step 4 to step 3
   ! ========================================================
   if (master) then
      call mpi_allgather(expe%pot%tot,1, mpi_double_precision, &
                         allexeptot, 1, mpi_double_precision, &
                         remd_comm, ierror)

      oex_eptot = allexeptot(o_repnum)
      myexeptot = allexeptot(myrepnum)

      ! Calculate the exchange probability if even # replica.
      if (mod(myrepnum,2)==0) then
         ! The below equation is generalized for each replica at a different
         ! temperature, but NOT exchanging that temperature (just coords).
         delta = 503.01d0 * (myexeptot - my_remd_data%myEptot) / &
                             my_remd_data%mytargettemp + 503.01d0 * &
                            (oex_eptot - o_eptot) / o_temp0
         metrop = exp(-delta + pvterm)

         if (jumpright(rem_dim)) then
            r_ex = 1
            r_fe = exp((my_remd_data%myEptot - oex_eptot) * 503.01d0 / &
                   my_remd_data%mytargettemp)
         else
            l_ex = 1
            l_fe = exp((my_remd_data%myEptot - oex_eptot) * 503.01d0 / &
                   my_remd_data%mytargettemp)
         end if

         ! Get random number between 0 and 1
         call amrand_gen(remd_gen, straw)

         ! Check for exchange
         exchange = straw < metrop

         ! send the results to the partner
         ! NOTE THAT WE NEED TO SUBTRACT 1 FROM REPNUM SINCE SUBREM USES
         ! INDEX OF 1 TO NUMREPS WHILE THE MPI PROCESSES ARE 0 TO NUMREPS-1
         call mpi_send(exchange, 1, mpi_logical, o_repnum-1, 1, &
                       remd_comm, ierror)

         ! Set replica coordinates based on exchange acceptance
         if (exchange) then
            ! Replace my coordinates with neighbor coordinates. Force isn't
            ! necessary right now because force() is called again immediately
            ! in the next step. This force call can be discarded, though, if
            ! forces are backed up. Note that x(lcrd) effectively indexes from 0
            ! whereas xtemp indexes from 1.
            do i=1, natom
               x(lcrd+3*i - 3) = xtemp(3*i - 2)
               x(lcrd+3*i - 2) = xtemp(3*i - 1)
               x(lcrd+3*i - 1) = xtemp(3*i)
               x(lforce+3*i - 3) = ftemp(3*i - 2)
               x(lforce+3*i - 2) = ftemp(3*i - 1)
               x(lforce+3*i - 1) = ftemp(3*i)
            end do

            ! Exchange the velocities and scale them for their new target
            ! temperature. (DSD 09/12)
            ! Use the ftemp array to hold velocities rather than allocate a
            ! temporary velocity array as well
            call mpi_sendrecv(x(lvel), 3*natom, mpi_double_precision,o_repnum-1, &
                        10, ftemp, 3*natom, mpi_double_precision, o_repnum-1, &
                        10, remd_comm, istatus, ierror)
            do i=1, natom
               x(lvel+3*i - 3) = ftemp(3*i - 2) * my_remd_data%myscaling
               x(lvel+3*i - 2) = ftemp(3*i - 1) * my_remd_data%myscaling
               x(lvel+3*i - 1) = ftemp(3*i) * my_remd_data%myscaling
            end do

            ! Exchange coordinate indices (even replica call)
            call mpi_sendrecv_replace(remd_crdidx, 1, mpi_integer, &
                                      o_repnum-1, 24, o_repnum-1, 24, &
                                      remd_comm, istatus, ierror)

            ! Increment exchsuccess for lower replica #
            ! Use index, not repnum, since this is a property of temperatures
            if (jumpright(rem_dim)) then
               ! this is the lower replica since this is the controlling rep
               ! NOTE THAT RANKS ARE 0 TO NUMREPS-1 BUT USE 1 TO NUMREPS FOR
               ! SUCCESS
               texchsuccess(myrepnum) = 1
            else
               ! other is the lower replica
               texchsuccess(o_repnum) = 1
            end if
         end if ! exchange

      else  ! (mod(myindex,2)==0)
         ! Not the replica controlling the exchange
         if (jumpright(rem_dim)) then
            l_ex = 1
            l_fe = exp((my_remd_data%myEptot - oex_eptot) * 503.01d0 / &
                    my_remd_data%mytargettemp)
         else
            r_ex = 1
            r_fe = exp((my_remd_data%myEptot - oex_eptot) * 503.01d0 / &
                   my_remd_data%mytargettemp)
         end if
         ! call rand to keep in sync with replicas that calculated Metropolis
         call amrand_gen(remd_gen, straw)

         call mpi_recv(exchange, 1, mpi_logical, &
                       o_repnum-1, 1, remd_comm, istatus, ierror)
         if (exchange) then
            do i=1, natom
               x(lcrd+3*i - 3) = xtemp(3*i - 2)
               x(lcrd+3*i - 2) = xtemp(3*i - 1)
               x(lcrd+3*i - 1) = xtemp(3*i)
               x(lforce+3*i - 3) = ftemp(3*i - 2)
               x(lforce+3*i - 2) = ftemp(3*i - 1)
               x(lforce+3*i - 1) = ftemp(3*i)
            end do

            ! Swap velocities using the ftemp array to hold
            call mpi_sendrecv(x(lvel), 3*natom, mpi_double_precision,o_repnum-1,&
                        10, ftemp , 3*natom, mpi_double_precision, o_repnum-1,&
                        10, remd_comm, istatus, ierror)
            do i=1, natom
               x(lvel+3*i - 3) = ftemp(3*i - 2) * my_remd_data%myscaling
               x(lvel+3*i - 2) = ftemp(3*i - 1) * my_remd_data%myscaling
               x(lvel+3*i - 1) = ftemp(3*i) * my_remd_data%myscaling
            end do

            ! Exchange coordinate indices (odd replica call)
            call mpi_sendrecv_replace(remd_crdidx, 1, mpi_integer, &
                                      o_repnum-1, 24, o_repnum-1, 24, &
                                      remd_comm, istatus, ierror)

         end if
      end if ! replica controlling the exchange (even #)

      ! toggle exchange direction
      jumpright(rem_dim) = .not. jumpright(rem_dim)

      ! REM: gather exchange log data
      call mpi_allreduce(texchsuccess, tempsuccess, remd_size,  &
                         MPI_INTEGER, mpi_sum, remd_comm, ierror)

      ! add the current successes to overall array
      do i = 1, remd_size
         exchsuccess(rem_dim, i) = exchsuccess(rem_dim, i) + tempsuccess(i)

         ! multiple fraction by 2.0 since we attempt this pair every OTHER
         ! exchange attempt (alternating directions of exchange)

         exchfrac(i) = dble(exchsuccess(rem_dim, i)) &
                     / dble(mdloop) * dble(remd_dimension) * TWO
      enddo

#ifdef VERBOSE_REMD
      write(6, '(a,f16.6)') '| My Eptot_1:       ', my_remd_data%myeptot
      write(6, '(a,f16.6)') '| My Eptot_2:       ', myexeptot
      write(6, '(a,f16.6)') '| Neighbor Eptot_1: ', o_eptot
      write(6, '(a,f16.6)') '| Neighbor Eptot_2: ', oex_eptot

      if (mod(myrepnum,2) == 0) then
         write(6, '(3(a,f16.6))') '| Delta= ', delta, ' Metrop= ', metrop, &
                                 ' Random #= ', straw
      else
         write(6, '(a)') '| Not controlling exchange.'
      end if

      if (exchange) then
         write(6, '(a)') '| Exchange Succeeded!'
      else
         write(6, '(a)') '| Exchange Failed!'
      end if
      write(6,'(a,26("="),a,26("="))') '| ', " END H-REMD CALCULATION"
#endif

   else
      ! call rand to keep rand generator in sync with masters
      call amrand_gen(remd_gen, straw)
   end if ! master

   call mpi_bcast(x(lcrd), 3*natom, mpi_double_precision, 0, commsander, ierror)
   call mpi_bcast(x(lforce), 3*natom, mpi_double_precision, 0, &
                  commsander, ierror)
   call mpi_bcast(x(lvel), 3*natom, mpi_double_precision, 0, commsander, ierror)
   if (exchange_ucell) then
     ! Broadcast exchange results so we know if unit cell needs to be restored.
     call mpi_bcast(exchange, 1, mpi_logical, 0, commsander, ierror)
     if (.not.exchange) then
       ! Restore original unit cell
       call remd_bcast_cell(box, original_box, commsander)
       !write(6,'(a,3f8.3)') 'DBG: Original box: ', box(1:3)
     endif
   endif
   ! initmodwt forces modwt() in nmr.f to re-read values such as temp0 -
   ! otherwise they will be reset to the initial value after each
   ! exchange.
   initmodwt = .true.

   ! Collect the information that we need to write to the remlog. Only masters
   ! do this part

   if (master) then
      ! Gather the r_fe and l_fe values from every thread

      call mpi_allgather(r_fe, 1, mpi_double_precision, &
                         all_r_fe, 1, mpi_double_precision, &
                         remd_comm, ierror)
      call mpi_allgather(l_fe, 1, mpi_double_precision, &
                         all_l_fe, 1, mpi_double_precision, &
                         remd_comm, ierror)
      call mpi_allgather(r_ex, 1, mpi_integer, all_r_ex, 1, mpi_integer, &
                         remd_comm, ierror)
      call mpi_allgather(l_ex, 1, mpi_integer, all_l_ex, 1, mpi_integer, &
                         remd_comm, ierror)

      ! Now print out the remlog data
      if (remd_master .and. rem > 0) then

         write(REMLOG_UNIT, '(a,i8)') '# exchange ', mdloop

         do rep = 1, remd_size

            ! Update FEP running sums and # of exchanges
            total_left_fe(1,1,rep) = total_left_fe(1,1,rep) + all_l_fe(rep)
            total_right_fe(1,1,rep) = total_right_fe(1,1,rep) + all_r_fe(rep)
            num_left_exchg(1,1,rep) = num_left_exchg(1,1,rep) + all_l_ex(rep)
            num_right_exchg(1,1,rep) = num_right_exchg(1,1,rep) + all_r_ex(rep)

            ! Calculate the left and the right free energy. We can overwrite
            ! what it was, since we already gathered it and is no longer
            ! necessary
            if (num_left_exchg(1,1,rep) > 0) then
               l_fe = alltemp0(rep) / 503.01d0 * &
                      log(total_left_fe(1,1,rep) / num_left_exchg(1,1,rep))
            else
               l_fe = 0.d0
            end if

            if (num_right_exchg(1,1,rep) > 0) then
               r_fe = alltemp0(rep) / 503.01d0 * &
                      log(total_right_fe(1,1,rep) / num_right_exchg(1,1,rep))
            else
               r_fe = 0.d0
            end if

            ! Figure out who our neighbor is. Note that jumpright has been
            ! toggled, so it's reversed of what it was above.
            if (jumpright(rem_dim)) then
               if (mod(rep, 2) == 0) then
                  neighbor_rep = rep - 1
               else
                  neighbor_rep = rep + 1
               end if
            else
               if (mod(rep, 2) == 0) then
                  neighbor_rep = rep + 1
               else
                  neighbor_rep = rep - 1
               end if
            end if
            ! Handle "boundary conditions"
            if (neighbor_rep > remd_size) neighbor_rep = 1
            if (neighbor_rep == 0) neighbor_rep = remd_size
            ! We succeeded either if our replica or our neighbor registers a
            ! success in the tempsuccess array
            if (tempsuccess(rep) == 1 .or. tempsuccess(neighbor_rep) == 1) then
               success_char = 'T'
            else
               success_char = 'F'
            end if

            write(REMLOG_UNIT, '(2i6,5f10.2,4x,a,2x,f10.2)') rep, &
               neighbor_rep, &
               alltemp0(rep), &
               alleptot(rep), &
               allexeptot(rep), &
               l_fe, &
               r_fe, &
               success_char, &
               exchfrac(index_list(rep))
         end do

      else if (remd_master) then
         ! Use texchsuccess to gather o_repnum's
         call mpi_gather(o_repnum, 1, mpi_integer, texchsuccess, 1, &
                         MPI_INTEGER, 0, remd_comm, ierror)
         do i = 1, remd_size
            neighbor_rep = texchsuccess(i)
            if (tempsuccess(i) == 1 .or. tempsuccess(neighbor_rep) == 1) then
               multid_print_data(i)%success = 1.d0
            else
               multid_print_data(i)%success = 0.d0
            end if
            multid_print_data(i)%neighbor_rep = neighbor_rep
            multid_print_data(i)%temp0        = alltemp0(i)
            multid_print_data(i)%pot_ene_tot  = alleptot(i)
            multid_print_data(i)%nei_pot_ene  = allexeptot(i)
            multid_print_data(i)%left_fe      = all_l_fe(i)
            multid_print_data(i)%right_fe     = all_r_fe(i)
            multid_print_data(i)%num_rep      = remd_size
            multid_print_data(i)%group_num    = group_num(rem_dim)
            multid_print_data(i)%repnum       = index_list(i)
            multid_print_data(i)%success_ratio = exchfrac(index_list(i))
            multid_print_data(i)%left_exchg   = all_l_ex(i)
            multid_print_data(i)%right_exchg  = all_r_ex(i)
         end do
      end if

      if (master .and. .not. remd_master) &
         call mpi_gather(o_repnum, 1, mpi_integer, texchsuccess, 1, &
                         MPI_INTEGER, 0, remd_comm, ierror)

      if (rem < 0) then
         ! Make sure everyone in our communicator trades their neighbor replica
         call mpi_gather(o_repnum, 1, mpi_integer, texchsuccess, 1, &
                         MPI_INTEGER, 0, remd_comm, ierror)
      end if
   end if ! master

   return

end subroutine hremd_exchange

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ calculation of exchange probability for Reservoir REMD
subroutine subrem_reservoir(rem_dim, reservoir_exchange_step)

   use constants, only : TWO

   implicit none

#  include "parallel.h"

   ! Passed arguments
   integer, intent (in) :: rem_dim     ! The dimension we are exchanging in
   integer, intent (in) :: reservoir_exchange_step ! Reservoir Exchange Frequency

   ! Local variables
   _REAL_ delta, straw, o_scaling, metrop
   _REAL_ pvterm
   _REAL_ l_temp0, r_temp0, o_temp0
   _REAL_ l_eptot, r_eptot, o_eptot
   _REAL_ alltempi(numreps), alltemp0(numreps), alleptot(numreps)
   ! For use in remlog output
   _REAL_ d_scaling(numreps), d_o_temp0(numreps)
   ! DAN ROE: for recording RREMD structure output in log
   integer d_rremd_idx(numreps)
   ! DAN ROE: cluster sizes for rremd==3
   integer myclustersize, o_clustersize
   ! index is the position of the replica's T in the sorted T list
   integer myindex, l_index, r_index, o_index
   ! repnum is the actual replica number
   integer l_repnum,r_repnum,o_repnum

   integer i, ierror, istatus(mpi_status_size)

   ! temporary for exchange success, will allgather from all nodes then
   ! add sum to actual exchsuccess array. also define an array for allgather sum
   integer texchsuccess(numreps), tempsuccess(numreps)
   _REAL_ exchfrac(numreps)

   logical exchange

! ---------------------
   ! Initialize variables
   texchsuccess(:) = 0

   delta = 0.0d0
   straw = 0.0d0
   metrop= 0.0d0
   pvterm= 0.0d0
   o_scaling = -0.1d0
   i = 0
   myindex = 0
   l_index = 0
   r_index = 0
   o_index = remd_rank
   myclustersize=0
   o_clustersize=0

   ! Call amrand to get a random# which can eventually be used as a structure
   ! index in rremd. Do it here so the random # generators on all the
   ! threads stay in sync. Only the master of the highest T replica during
   ! jumpright will ever use this.
   ! Calling this even when we have no rremd will allow random # gen to stay
   ! in sync between remd and rremd runs.
   call amrand_gen(remd_gen, straw)

   ! Only the sander masters do the exchanges.
   if (master) then

      ! Gather current temperature, target temperature, and PE
      call mpi_allgather(my_remd_data%mytemp, 1, mpi_double_precision, &
                         alltempi, 1, mpi_double_precision, &
                         remd_comm, ierror)
      call mpi_allgather(my_remd_data%mytargettemp, 1, mpi_double_precision, &
                         alltemp0, 1, mpi_double_precision, &
                         remd_comm, ierror)
      call mpi_allgather(my_remd_data%myEptot, 1, mpi_double_precision, &
                         alleptot, 1, mpi_double_precision, &
                         remd_comm, ierror)

      ! alltemp0 and alleptot are indexed by replica #
      ! statetable is SORTED

      myindex = replica_indexes(rem_dim)

      ! DAN ROE: Why do we look for both partners?

      ! Find neighbor Ts for exchange using Temperature list indexed
      !  by TEMPERATURE.
      ! Wrap around so replicas exchange in a circle,
      !  i.e. highest and lowest T may attempt exchange.
      l_index = myindex - 1
      if (l_index < 1) l_index = remd_size
      r_index = myindex + 1
      if (r_index > remd_size) r_index = 1

      ! Now pull the replica number of the adjacent replicas (partners array)
      l_repnum = partners(1)
      r_repnum = partners(2)

      ! Now get their temperatures
      l_temp0 = alltemp0(l_repnum)
      r_temp0 = alltemp0(r_repnum)

      ! Get the energies of these neighbor replicas
      l_eptot = alleptot(l_repnum)
      r_eptot = alleptot(r_repnum)


      ! Set up partner replica information.
      ! By definition, only even replica #s will initiate the exchange.
      if(mod(myindex, 2) == 0) then
         ! this replica will calculate delta
         if (jumpright(rem_dim)) then
            ! Partner is to the right, i.e. at higher T
            ! DAN ROE: Integrate check for reservoir here?
            o_repnum = r_repnum
            o_index  = r_index
            o_eptot  = r_eptot
            o_temp0  = r_temp0
         else
            ! Partner is to the left, i.e. at lower T
            o_repnum = l_repnum
            o_index  = l_index
            o_eptot  = l_eptot
            o_temp0  = l_temp0
         end if
      else
         ! This replica will not calculate delta, right and left are
         !  switched compared to controlling replica.
         ! Only needs o_eptot for writing debug info
         if (jumpright(rem_dim)) then
            o_repnum = l_repnum
            o_index  = l_index
            o_eptot  = l_eptot
            o_temp0  = l_temp0
         else
            o_repnum = r_repnum
            o_index  = r_index
            o_eptot  = r_eptot
            o_temp0  = r_temp0
         end if
      end if  ! (even/odd rank check)

      ! RREMD: If jumpright and this replica has highest T, then attempt an
      ! exchange with a random structure in the reservoir.
      ! Read coordinates from the file when we are in sander.
      ! Initialize rremd_idx, so that it is -1 for any replica other than the
      !  one that will use the reservoir coordinates.
      ! Change the partners for RREMD: lowest does not exchange, and highest
      !  uses the reservoir.
      rremd_idx = -1
      ! Highest replica will exchange with reservoir only if
      ! mod(mdloop, reservoir_exchange_step) is 1.
      ! Because the number of replicas should be even, the highest replica will always
      ! have an even index. The even replicas control the exchange (jumpright(rem_dim)) when
      ! the exchange number is odd. Therefore for the highest replica to exchange with
      ! the reservoir, the exchange number has to be odd and hence the modulus should
      ! be equal to 1 and not zero. Since we attempt exchange with the reservoir every
      ! alternate exchange, the reservoir_exchange_step has to be a
      ! multiple of 2 (condition checked for in mdread2.F90).
      if(mod(mdloop, reservoir_exchange_step) .eq. 1) then
         if (jumpright(rem_dim)) then
            if (myindex == 1) then
               ! Lowest T, will not do exchange or even wait for results from
               !  highest T
               o_index=-999
            else if (myindex == remd_size) then
               ! Highest T, will exchange with reservoir
               o_temp0 = restemp0
               o_index = -900
               ! Get random structure index
               ! amrand has already been called at the start of the routine so
               !  that all threads stay in sync.
               rremd_idx = int(straw * reservoirsize) + 1
               ! DAN ROE: Random number should not be outside reservoir size!
               ! By definition it wont be unless something breaks in amrand()
               if (rremd_idx.gt.reservoirsize) rremd_idx = reservoirsize
               ! Get PE, read in previously from saveene
               o_eptot = saveene(rremd_idx)
            end if
         end if ! jumpright(rem_dim)
      end if ! reservoir_exchange_step

      ! Calculate pressure/volume correction if necessary
      if (use_pv) pvterm = pv_correction(my_remd_data%mytargettemp, o_temp0, o_repnum-1, remd_comm)

      ! Calculate the exchange probability if even # replica.
      if(mod(myindex, 2) == 0) then
         ! RREMD: Here we calculate delta beta term for normal exchanges but
         !  only beta for exchanges with non-Boltzmann weighting (rremd>1)
         ! Reservoir exchanges are for top replica and jumpright only!
         ! Note that 503.01 is 1/kB in internal amber units.
         ! Exchange with reservoir only on reservoir_exchange_step
         if (rremd==2 .and. myindex==remd_size .and. jumpright(rem_dim) &
                  .and. mod(mdloop, reservoir_exchange_step) .eq. 1) then
            ! No delta beta, assign weight of 1/N to each structure in
            !  reservoir. The exchange criterion becomes:
            !  metrop = exp[-beta_replica * (E_reservoir-E_replica)]
            ! NB 1/N delta:
            delta = (o_eptot - my_remd_data%myEptot) * 503.01d0 &
                  / my_remd_data%mytargettemp
         ! Exchange with reservoir only on reservoir_exchange_step
         else if (rremd==3 .and. myindex==remd_size .and. jumpright(rem_dim) &
                  .and. mod(mdloop, reservoir_exchange_step) .eq. 1) then
            ! No delta beta, weight of each structure in the reservoir is
            !  user defined and has been read previously.
            ! runmd has determined which cluster the MD structure is in
            !  (incluster).
            ! Weights of each structure are obtained by clustersize(clusternum(i))
            !  or for the MD structure clustersize(incluster)
            ! Note that MD structures not present in the reservoir are given a
            !  weight of zero currently. In principle we could add to the reservoir,
            !  but that is for a future release.
            ! This makes the exchange criterion:
            !  metrop =  w_rep / w_res * exp[-beta_rep * (E_res-E_rep)]
            ! NB M/N delta:
            delta = (o_eptot - my_remd_data%myEptot) * 503.01d0 &
                  / my_remd_data%mytargettemp
            myclustersize = clustersize(incluster)
            o_clustersize = clustersize(clusternum(rremd_idx))
         else
            ! Std REMD, or RREMD with Boltzmann weighted reservoir, or NB RREMD
            !  but not exchanging with the reservoir this time.
            delta = (my_remd_data%myEptot - o_eptot) &
                  * (my_remd_data%mytargettemp - o_temp0) * 503.01d0 &
                  / (my_remd_data%mytargettemp * o_temp0)
         end if ! REMD exchange calculation

         metrop = exp(-delta + pvterm)
         if (rremd==3 .and. myindex==remd_size .and. jumpright(rem_dim) &
                 & .and. mod(mdloop, reservoir_exchange_step) .eq. 1) then
            metrop = metrop * myclustersize / o_clustersize
         end if

         ! Get random number between 0 and 1
         call amrand_gen(remd_gen, straw)

         ! Check for exchange
         exchange = straw < metrop

         ! If exchanged, set scaling, otherwise dont
         if (exchange) then
            ! Set velocity scaling
            my_remd_data%myscaling = sqrt(o_temp0 / my_remd_data%mytargettemp)
            o_scaling = 1.0d0 / my_remd_data%myscaling

            ! If RREMD and exchanging with the reservoir we are changing
            !  structure, not temperature, so myscaling is inverted.
            ! Dont scale if we aren't reading velocities in the reservoir.
            if(myindex == remd_size .and. jumpright(rem_dim) &
                    .and. mod(mdloop, reservoir_exchange_step) .eq. 1) then
               if (reserv_velo) then
                  my_remd_data%myscaling = o_scaling
               else
                  my_remd_data%myscaling = 1.d0
               end if
            end if
            ! Increment exchsuccess for lower replica #
            ! Use index, not repnum, since this is a property of temperatures
            if (jumpright(rem_dim)) then
               ! this is the lower temperature since this is the controlling rep
               ! NOTE THAT RANKS ARE 0 TO NUMREPS-1 BUT USE 1 TO NUMREPS FOR
               ! SUCCESS
               texchsuccess(myindex) = 1
            else
               ! other is the lower temperature
               texchsuccess(myindex-1) = 1
            end if
         else
            ! No exchange
            my_remd_data%myscaling = -1.0d0
            o_scaling = -1.0d0
         end if ! exchange

! CARLOS: DO WE NEED BARRIER?

         call mpi_barrier(remd_comm, ierror)

         ! send the results to the partner
         ! NOTE THAT WE NEED TO SUBTRACT 1 FROM REPNUM SINCE SUBREM USES
         ! INDEX OF 1 TO NUMREPS WHILE THE MPI PROCESSES ARE 0 TO NUMREPS-1

         ! RREMD: Dont communicate if we exchanged with reservoir
         if (o_index > 0) then
            !write (6,*) "sending o_scaling of ",o_scaling,&
            !            "to replica ",o_repnum
            call mpi_send(o_scaling, 1, mpi_double_precision, &
                          o_repnum-1, 0, remd_comm, ierror)
         end if
      else
         ! Not the replica controlling the exchange
         ! call rand to keep in sync with replicas that calculated Metropolis
         call amrand_gen(remd_gen, straw)
         ! DAN ROE: FileDebug

! CARLOS: NEED BARRIER?
         call mpi_barrier(remd_comm, ierror)
         ! receive the scaling data; if it is >0 we know exchange succeeded
         ! RREMD: Don't communicate if our partner is exchanging with
         !  the reservoir
         if (o_index > 0) then
            call mpi_recv(my_remd_data%myscaling, 1, mpi_double_precision, &
                          o_repnum-1, 0, remd_comm, istatus, ierror)
         else
            ! RREMD and lowest T replica, partner attempted exchange with reservoir
            my_remd_data%myscaling = -1.d0
         end if
      end if ! replica controlling the exchange (even #)

      ! toggle exchange direction

      call mpi_barrier(remd_comm, ierror)
      jumpright(rem_dim) = .not. jumpright(rem_dim)

      if(my_remd_data%myscaling < 0.0d0) then
         my_remd_data%newtargettemp = my_remd_data%mytargettemp
         exchange=.false.
         ! If RREMD, didnt exchange with the bath so get rid of the rremd_idx
         !  that we set. Ok to do for all replicas.
         ! DAN ROE: DO this later
         !rremd_idx=-1
      else
         ! RREMD: don't change temp0 if we exchanged with the reservoir.
         !  Instead, change coordinates to those of the new structure.
         ! HOW? Write over restart?
         ! OK to leave my_scaling since technically the new structure needs
         !  scaling.
         if (o_index > 0) then
            my_remd_data%newtargettemp = o_temp0
            exchange=.true.
         else ! RREMD highest replica
            my_remd_data%newtargettemp = my_remd_data%mytargettemp
            exchange = .true.
         end if
      end if

      if (exchange.and.o_index.gt.0) then
         ! Swap replica index and all of our group affiliations
         ! Do not do this if exchanging with the reservoir.
         call mpi_sendrecv_replace(replica_indexes, remd_dimension, &
                     mpi_integer, o_repnum-1, 22, o_repnum-1, 22, &
                     remd_comm, istatus, ierror)
         call mpi_sendrecv_replace(group_num, remd_dimension, mpi_integer, &
                     o_repnum-1, 23, o_repnum-1, 23, remd_comm, istatus, ierror)
         call mpi_sendrecv_replace(remd_repidx, 1, mpi_integer, &
                                   o_repnum-1, 24, o_repnum-1, 24, &
                                   remd_comm, istatus, ierror)
      end if

      ! REM: gather exchange log data

      call mpi_gather(my_remd_data%myscaling, 1, mpi_double_precision, &
                      d_scaling, 1, mpi_double_precision, &
                      0, remd_comm, ierror)
      call mpi_gather(my_remd_data%newtargettemp, 1, mpi_double_precision, &
                      d_o_temp0, 1, mpi_double_precision, &
                      0, remd_comm, ierror)

      ! DAN ROE: RREMD: gather reservoir structure data
      call mpi_gather(rremd_idx, 1, MPI_INTEGER, &
                      d_rremd_idx, 1, MPI_INTEGER, &
                      0, remd_comm, ierror)

      call mpi_allreduce(texchsuccess, tempsuccess, remd_size,  MPI_INTEGER, &
                         MPI_SUM, remd_comm, ierror)

      ! add the current successes to overall array

      do i = 1, remd_size
         exchsuccess(rem_dim, i) = exchsuccess(rem_dim, i) + tempsuccess(i)

         ! multiple fraction by 2.0 since we attempt this pair every OTHER
         ! exchange attempt (alternating directions of exchange)

         exchfrac(i) = dble(exchsuccess(rem_dim, i)) &
                     / dble(mdloop) * dble(remd_dimension) * TWO
      enddo


      if (remd_master .and. rem > 0) then

      ! only overall master writes the log
      ! DAN ROE: added d_rremd_idx write

         write(unit = REMLOG_UNIT, fmt = '(a,i8)') '# exchange ', mdloop
         do i = 1, remd_size
            write(REMLOG_UNIT, '(i2,f12.4,f10.2,f15.2,3f10.2,i8)') &
            i, & ! Replica #
            d_scaling(i), & ! scaling factor
            alltempi(i), & ! current temperature
            alleptot(i), & ! current potential energy
            alltemp0(i), & ! current target temperature
            d_o_temp0(i),& ! next target temperature
            exchfrac(index_list(i)), &
            d_rremd_idx(i) ! structure# attempted from  reservoir
         end do
      else if (remd_master) then
         do i = 1, remd_size
            multid_print_data(i)%scaling     = d_scaling(i)
            multid_print_data(i)%real_temp   = alltempi(i)
            multid_print_data(i)%pot_ene_tot = alleptot(i)
            multid_print_data(i)%temp0       = alltemp0(i)
            multid_print_data(i)%new_temp0   = d_o_temp0(i)
            multid_print_data(i)%num_rep     = remd_size
            multid_print_data(i)%group_num   = group_num(rem_dim)
            multid_print_data(i)%success_ratio = &
               dble(exchsuccess(rem_dim, index_list(i))) / &
               dble(mdloop) * dble(remd_dimension) * 2
         end do
      end if

      ! DAN ROE: All sander masters write REMD exchange info
#ifndef VERBOSE_REMD
      if (rremd > 0) then
#endif /* NOTE ifNdef */
         write(6,'(a6,i10,a8,i1)') "| Exch= ", mdloop, " RREMD= ", rremd
         ! This Replica Information
         write(6,'(a16,a7,f6.2,2(a7,i2),a7,f10.2)') &
          "| Replica         "," Temp= ", my_remd_data%mytargettemp, &
          " Indx= ", myindex, " Rep#= ", repnum, " EPot= ", my_remd_data%myEptot
         write(6,'(a7,i6,a8,i6)') '| RepIdx=', remd_repidx, ' CrdIdx=', remd_crdidx
         ! Partner Information
         ! use .not.jumpright since jumpright has already been toggled
         if (rremd>0 .and. .not. jumpright(rem_dim) &
                     .and. myindex == remd_size &
                     .and. mod(mdloop, reservoir_exchange_step) .eq. 1) then
            ! Partner is Reservoir
            write(6,'(a16,a7,f6.2,a10,i8,a7,f10.2)') &
               "Reservoir       ", " Temp= ", restemp0, &
               " Struct#= ", rremd_idx, " EPot= ", o_eptot
            if (exchange) then
               write(6,'(a20)') "ReservoirExchange= T"
            else
               write(6,'(a20)') "ReservoirExchange= F"
            end if
            if (rremd==3) then
               ! Reservoir has weights
               write(6,'(a11,i10)') "mycluster= ",incluster
               write(6,'(a15,i10,a16,i10)') &
                  "myclustersize= ", myclustersize,&
                  " o_clustersize= ", o_clustersize
            end if
         else if (rremd > 0 .and. .not. jumpright(rem_dim) &
                            .and. myindex == 1 &
                            .and. mod(mdloop, reservoir_exchange_step) .eq. 1) then
            ! Lowest T, partner would be highest T but that is exchanging
            ! with the reservoir so lowest T has no partner.
            write(6,'(a)') &
               "No partner, highest T exchanging w/ Reservoir."
         else
            ! Partner is normal replica
            write(6,'(a16,a7,f6.2,2(a7,i2),a7,f10.2)') &
               "| Partner         "," Temp= ",o_temp0," Indx= ",o_index, &
               " Rep#= ", o_repnum, " EPot= ", o_eptot
         end if ! not jumpright and myindex == remd_size
         ! Exchange calculation information
         if (mod(myindex,2)==0) then
            ! This replica controlled exchange and calculated metrop
            write(6,'(a8,E16.6,a8,E16.6,a12,f10.2)') &
               "| Metrop= ",metrop," delta= ",delta," o_scaling= ",o_scaling
         else
            ! This replica did not control exchange
            write(6,'(a)') "| Not controlling exchange."
         end if ! mod(myindex,2)
         ! Write random #, scaling, and success
         write(6,'(a8,E16.6,a12,f10.2,a10,L1)' ) &
            "| Rand=   ",straw," MyScaling= ", my_remd_data%myscaling, &
            " Success= ",exchange
#ifndef VERBOSE_REMD
      end if ! rremd>0, RREMD information writeout
#endif /* NOTE ifNdef */

      ! If no exchange occurred, reset rremd_idx
      if (.not. exchange) rremd_idx = -1

   else  ! not part of remd_comm
      ! Keep random generator in sync
      call amrand_gen(remd_gen, straw)
   end if

   return

end subroutine subrem_reservoir

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Perform actions necessary to exchange replicas in Reservoir REMD
subroutine reservoir_remd_exchange(rem_dim, rem_kind, x, v, amass, nr3, &
                         natom, nr, temp0, reservoir_exchange_step)

   use sgld, only : isgld, rxsgld_scale
   
implicit none

#  include "parallel.h"

! Passed variables
   _REAL_, intent(in out) :: x(*)      ! Atomic coordinates
   _REAL_, intent(in out) :: v(*)      ! Atomic velocities
   _REAL_, intent(in)     :: amass(*)  ! Mass array (not inverted)
   integer, intent(in)    :: nr3       ! 3 * number of atoms
   integer, intent(in)    :: natom     ! number of atoms
   integer, intent(in)    :: nr        ! number of atoms (?)
   integer, intent(in)    :: rem_dim   ! Which dimension we are exchanging in
   integer, intent(in)    :: rem_kind  ! Which kind of remd are we doing?
   _REAL_, intent(in out) :: temp0     ! Current target temperature
   integer, intent(in)    :: reservoir_exchange_step ! Reservoir Exchange Frequency

   integer ierror

!-------------------
#  ifdef VERBOSE_REMD
   if (remd_master) &
     write(6,'(a,26("="),a,26("="))') '| ', "RREMD EXCHANGE CALCULATION"
#  endif
   ! ---=== RREMD DIHEDRAL CLUSTERING CALC.===---
   ! For rremd==3 calculate cluster of current structure.
   if (master .and. rremd == 3) then
      call calc_rremd_cluster(x)
   end if

   ! Set the partners array so we know who our partners are. NOTE, this also
   ! sets index_list
   if (master .and. .not. isgld > 0) call set_partners(rem_dim, remd_size)

   ! ---=== REMD EXCHANGE CALCULATION ===---
   ! Attempt Exchange with subrem_reservoir() using energy from previous runmd call
   ! All procs call subrem so slave random # generators stay in sync
   call subrem_reservoir(rem_dim, reservoir_exchange_step)

   ! ---=== REMD COMMUNICATION ===---
   ! Broadcast reservoir structure index. All threads needs to know
   !  this to decide if they have to receive coords or not.
   if (rremd>0) then
      call mpi_bcast(rremd_idx, 1, mpi_integer, 0, commsander, ierror)
   else
     ! The above condition should be true. Otherwise there is 
     ! something terribly wrong here.
     write(6,*) 'ERROR: rremd has to be greater than 0 &
            &if rem == 6'
     call mexit(6,1)
   end if
   ! Broadcast my_remd_data
   call mpi_bcast(my_remd_data, SIZE_REMD_DATA, mpi_double_precision, &
                  0, commsander, ierror)


   ! ---=== RREMD RESERVOIR LOADING ===---
   ! If we exchanged with the reservoir, swap the coordinates
   !  here (after the inpcrd have been read).
   if (rremd_idx>0) &
      call load_reservoir_structure(x, v, nr3, natom)


   if (isgld > 0) then
   ! ---=== RXSGLD property SCALING ===---
      if (my_remd_data%newtargettemp > 0.0) &
         temp0 = my_remd_data%newtargettemp
      call rxsgld_scale(stagid, nr, my_remd_data%myscaling, amass, x, v)
   else
      ! ---=== REMD VELOCITY SCALING ===---
      ! REMD: If an attempt is accepted, set the target temperature to
      !  the new one and rescale velocities.
#ifdef LES
      call remd_scale_velo(v, temp0, nr, nr3, rem_kind)
#else
      call remd_scale_velo(v, temp0, nr3, rem_kind)
#endif
   end if

   ! initmodwt forces modwt() in nmr.f to re-read values such as temp0 -
   ! otherwise they will be reset to the initial value after each
   ! exchange.
   initmodwt = .true.
#  ifdef VERBOSE_REMD
   if (remd_master) &
     write(6,'(a,24("="),a,24("="))') '| ', "END RREMD EXCHANGE CALCULATION"
#  endif
   return

end subroutine reservoir_remd_exchange

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Manages exchange attempts in multiple dimensions
subroutine multid_remd_exchange(x, ix, ih, ipairs, qsetup, &
                                do_list_update, temp0, solvph, solve, &
                                reservoir_exchange_step) 

   use sander_lib, only : strip

   implicit none
#  include "../include/memory.h"
#  include "parallel.h"

   ! Formal arguments
   _REAL_, intent (in out)  :: x(*)     ! Global real array (NOT just coords)
   _REAL_, intent (in out)  :: temp0    ! Target temperature
   _REAL_, intent (in out)  :: solvph   ! Target solvent pH
   _REAL_, intent (in out)  :: solve    ! Target solvent Redox potential
   integer, intent (in out) :: ix(*)    ! Global integer array
   character(4), intent(in) :: ih(*)    ! Global hollerith array
   integer, intent (in out) :: ipairs(*)! Pairlist (?)
   integer, intent (in)     :: reservoir_exchange_step! Reservoir Exchange Frequency

   logical                  :: qsetup
   logical                  :: do_list_update

   ! Local variables

   integer :: my_dim       ! The dimension we are exchanging in
   integer :: i, j         ! Loop counters
   integer :: rep, gid     ! rep and group index
   integer :: ierror       ! Error code
   logical :: group_master ! Master of the group
   character :: torf_char  ! store if we are true or false
   _REAL_  :: r_fe, l_fe   ! right and left free energies for REFEP calculation

   character(len=MAX_FN_LEN) :: filename ! Name of remlog we need to write to
   character(len=MAX_FN_LEN) :: buf      ! Filename buffer

   ! Determine my_dim, starting with the lowest dimension. We want my_dim to
   ! cycle between 1 - remd_dimension, starting at 1

   my_dim = mod(mdloop - 1, remd_dimension) + 1

   ! Nobody is group master yet, this will be set later

   group_master = .false.

   ! Set up REMD communicators, and free them if they're formed

   if (master) then

      if (remd_comm /= mpi_comm_null) &
         call mpi_comm_free(remd_comm, ierror)
      remd_comm = mpi_comm_null

      call mpi_comm_split(commmaster, group_num(my_dim), &
                          masterrank, remd_comm, ierror)
      call mpi_comm_size(remd_comm, remd_size, ierror)
      call mpi_comm_rank(remd_comm, remd_rank, ierror)

      remd_master = remd_rank == 0

      if (group_master_comm /= mpi_comm_null) &
         call mpi_comm_free(group_master_comm, ierror)
      group_master_comm = mpi_comm_null

      if (remd_master) then
         call mpi_comm_split(commmaster, 0, masterrank, &
                             group_master_comm, ierror)
         call mpi_comm_rank(group_master_comm, group_master_rank, ierror)
         call mpi_comm_size(group_master_comm, group_master_size, ierror)
         group_master = group_master_rank == 0
      else
         ! If we are not remd_master, we'll call mpi_comm_split with
         ! MPI_UNDEFINEDs
         call mpi_comm_split(commmaster, MPI_UNDEFINED, MPI_UNDEFINED, &
                             group_master_comm, ierror)
      end if

      multid_print_data(:) = NULL_REMLOG_DATA
      multid_print_data_buf(:,:) = NULL_REMLOG_DATA

   end if ! master

   select case (remd_types(my_dim))

      case (1)
         if(rremd>0) then
            call reservoir_remd_exchange(my_dim, 1, x(lcrd), x(lvel), x(lmass), &
                            natom*3, natom, natom, temp0, reservoir_exchange_step)
         else
            call remd_exchange(my_dim, 1, x(lcrd), x(lvel), x(lmass), &
                            natom*3, natom, natom, temp0)
         end if
      case (3)
         call hremd_exchange(my_dim, x, ix, ih, ipairs, qsetup, do_list_update)

   end select

   ! Now collect all of the multid_print_data's for all of the various
   ! remd_masters to the group_master

   if (remd_master) then
      call mpi_allgather(multid_print_data, SIZE_REMLOG_DATA*numreps, &
                         mpi_double_precision, multid_print_data_buf, &
                         SIZE_REMLOG_DATA*numreps, mpi_double_precision, &
                         group_master_comm, ierror)

      ! Update the free energy data structures if this dimension is H-REMD

      if (remd_types(my_dim) == 3) then

         do i = 1, group_master_size

            gid = int(multid_print_data_buf(1,i)%group_num)

            do j = 1, int(multid_print_data_buf(1, i)%num_rep)

               rep = int(multid_print_data_buf(j,i)%repnum)
               num_left_exchg(my_dim, gid, rep) = &
                  num_left_exchg(my_dim, gid, rep) + &
                  int(multid_print_data_buf(j,i)%left_exchg)
               num_right_exchg(my_dim, gid, rep) = &
                  num_right_exchg(my_dim, gid, rep) + &
                  int(multid_print_data_buf(j,i)%right_exchg)
               total_left_fe(my_dim, gid, rep) = &
                  total_left_fe(my_dim, gid, rep) + &
                  multid_print_data_buf(j,i)%left_fe
               total_right_fe(my_dim, gid, rep) = &
                  total_right_fe(my_dim, gid, rep) + &
                  multid_print_data_buf(j,i)%right_fe

            end do

         end do ! i = 1, group_master_size

      end if ! remd_types(my_dim) == 3

   end if

   ! Now that we have all of the data on group_master, have *that* process
   ! open up the appropriate remlog and write out the data to that file. We
   ! call the intrinsic "open" so we can append to this old file (it should
   ! have already been opened via amopen in the setup routine).  This has
   ! the side-effect of flushing this rem.log every time we reach this point

   if (group_master) then

      write(buf, '(i5)') my_dim
      call strip(buf)
      filename = trim(remlog) // '.' // trim(buf)

      open(unit=REMLOG_UNIT, file=filename, status='OLD', position='APPEND')

      ! Modify your particular REMD's exchange method printout here:

      select case (remd_types(my_dim))

         case (1) ! TEMPERATURE

            do i = 1, group_master_size

               gid = int(multid_print_data_buf(1,i)%group_num)

               write(REMLOG_UNIT, '(2(a,i8))') &
                  '# exchange ', mdloop, ' REMD group ', gid

               do j = 1, int(multid_print_data_buf(1, i)%num_rep)
                  write(REMLOG_UNIT, '(i2,6f10.2,i8)') &
                     j, &
                     multid_print_data_buf(j,i)%scaling,       &
                     multid_print_data_buf(j,i)%real_temp,     &
                     multid_print_data_buf(j,i)%pot_ene_tot,   &
                     multid_print_data_buf(j,i)%temp0,         &
                     multid_print_data_buf(j,i)%new_temp0,     &
                     multid_print_data_buf(j,i)%success_ratio, &
                     int(multid_print_data_buf(j,i)%struct_num)
               end do

            end do

         case (3) ! HAMILTONIAN

            do i = 1, group_master_size

               gid = int(multid_print_data_buf(1,i)%group_num)

               write(REMLOG_UNIT, '(2(a,i8))') &
                  '# exchange ', mdloop, ' REMD group ', gid

               do j = 1, int(multid_print_data_buf(1, i)%num_rep)

                  rep = int(multid_print_data_buf(j,i)%repnum)

                  if (num_right_exchg(my_dim, i, rep) .eq. 0) then
                     r_fe = 0.d0
                  else
                     r_fe = multid_print_data_buf(j,i)%temp0 / 503.01d0 * &
                            log(total_right_fe(my_dim, gid, rep) / &
                            num_right_exchg(my_dim, gid, rep))
                  end if

                  if (num_left_exchg(my_dim, i, rep) .eq. 0) then
                     l_fe = 0.d0
                  else
                     l_fe = multid_print_data_buf(j,i)%temp0 / 503.01d0 * &
                            log(total_left_fe(my_dim, gid, rep) / &
                            num_left_exchg(my_dim, gid, rep))
                  end if

                  if (multid_print_data_buf(j,i)%success > 0.d0) then
                     torf_char = 'T'
                  else
                     torf_char = 'F'
                  end if

                  write(REMLOG_UNIT, '(2i6,5f10.2,4x,a,2x,f10.2)') j, &
                     int(multid_print_data_buf(j,i)%neighbor_rep),    &
                     multid_print_data_buf(j,i)%temp0,                &
                     multid_print_data_buf(j,i)%pot_ene_tot,          &
                     multid_print_data_buf(j,i)%nei_pot_ene,          &
                     l_fe,      &
                     r_fe,      &
                     torf_char, &
                     multid_print_data_buf(j,i)%success_ratio

               end do

            end do

         case (4) ! pH

            do i = 1, group_master_size

               gid = int(multid_print_data_buf(1,i)%group_num)

               write(REMLOG_UNIT, '(2(a,i8))') &
                  '# exchange ', mdloop, ' REMD group ', gid

               do j = 1, int(multid_print_data_buf(1,i)%num_rep)
                  write(REMLOG_UNIT, '(i6,x,i7,x,2f7.3,x,f8.4)') j, &
                     nint(multid_print_data_buf(j,i)%nprot),        &
                     multid_print_data_buf(j,i)%my_ph,              &
                     multid_print_data_buf(j,i)%nei_ph,             &
                     multid_print_data_buf(j,i)%success_ratio
               end do

            end do

         case (5) ! REDOX

            do i = 1, group_master_size

               gid = int(multid_print_data_buf(1,i)%group_num)

               write(REMLOG_UNIT, '(2(a,i8))') &
                  '# exchange ', mdloop, ' REMD group ', gid

               do j = 1, int(multid_print_data_buf(1,i)%num_rep)
                  write(REMLOG_UNIT, '(i6,x,i7,x,2f7.3,x,f8.4)') j, &
                     nint(multid_print_data_buf(j,i)%nelec),        &
                     multid_print_data_buf(j,i)%my_e,               &
                     multid_print_data_buf(j,i)%nei_e,              &
                     multid_print_data_buf(j,i)%success_ratio
               end do

            end do

      end select

      close(REMLOG_UNIT)

   end if ! group_master

   ! Set the next remd method
   next_rem_method = remd_types(mod(mdloop, remd_dimension) + 1)

   return

end subroutine multid_remd_exchange
#endif /* MPI */

end module remd
