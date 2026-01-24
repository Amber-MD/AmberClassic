integer commworld, commsander
integer worldrank
integer worldsize
integer numtasks,mytaskid

! Following are not used
integer commmaster, comm_lpimd
integer sanderrank, masterrank, lpimd_rank
integer sandersize, mastersize, lpimd_size
integer groupmaster, worldmaster
logical ng_sequential
integer iparpt
integer iparpt3,rcvcnt,rcvcnt3
logical mpi_orig
#undef  MPI_MAX_PROCESSORS
#define MPI_MAX_PROCESSORS 256
dimension iparpt(0:MPI_MAX_PROCESSORS)
dimension iparpt3(0:MPI_MAX_PROCESSORS)
dimension rcvcnt(0:MPI_MAX_PROCESSORS)
dimension rcvcnt3(0:MPI_MAX_PROCESSORS)
integer notdone

! Common block is for compatability:
common/parallel/numtasks,mytaskid,notdone, &
      iparpt,iparpt3,rcvcnt,rcvcnt3,mpi_orig
common/parallel_multi/commworld, commsander, commmaster, comm_lpimd, &
      worldrank, sanderrank, masterrank, lpimd_rank, &
      worldsize, sandersize, mastersize, lpimd_size, &
      groupmaster, worldmaster, ng_sequential
