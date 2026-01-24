!+ Specification and control of Amber's Input/Output
!+ R. Luo added pqr support

! File names
character(len=512) groupbuffer
character(len=80) mdin, mdout, inpcrd, parm, restrt, &
      refc, mdvel, mden, mdcrd, mdinfo, nmr, mincor, &
      vecs, radii, freqe,redir(8), &
      rstdip,mddip,inpdip,groups,gpes, &
      cpin, cpout, cprestrt, newtop, & ! B.Aguilar: newtop for mdar6
      pqr &
!#ifdef MMTSB
!      ,mmtsb_setup_file &
!#endif
!cnt N.TAKADA: For MDM
!#ifdef MDM_PDB
!      ,mdpdb &
!#endif
!cnt N.TAKADA: For MDM
; ! line terminator for free-form version

character owrite
common /files/ groupbuffer, mdin, mdout, inpcrd, parm, restrt, &
      refc, mdvel, mden, mdcrd, mdinfo, nmr, mincor, &
      vecs, radii, freqe, owrite, &
      rstdip,mddip,inpdip,groups,gpes, &
      cpin, cpout, cprestrt, newtop, & !B.Aguilar for mdar6 
      pqr &
!#ifdef MMTSB
!      ,mmtsb_setup_file &
!#endif
!cnt N.TAKADA: For MDM
!#ifdef MDM_PDB
!      ,mdpdb &
!#endif
!cnt N.TAKADA: For MDM
; ! line terminator for free-form version

! put this in a seperate common block to stop the compiler from
! complaining about misalignment
integer numgroup
common/nmgrp/ numgroup


! File units
! An I/O Unit resource manager does not exist.
integer     MDCRD_UNIT
integer     MDEN_UNIT
integer     MDINFO_UNIT
integer     MDVEL_UNIT
parameter ( MDINFO_UNIT =  7 )
parameter ( MDCRD_UNIT  = 12 )
parameter ( MDEN_UNIT   = 15 )
parameter ( MDVEL_UNIT  = 13 )
integer, parameter :: CNSTPH_UNIT = 18, CPOUT_UNIT = 19
integer, parameter :: PQR_UNIT = 20

! 18 was picked because CNSTPH uses it; conflicts are not expected.
!integer     MMTSB_UNIT
!parameter ( MMTSB_UNIT = 18 )


! File related controls and options
character(len=80) title,title1
common/runhed/ title, title1

logical mdin_pb &
!#ifdef MDM_MD
!      ,mdin_mdm &
!#endif
!#ifdef MDM_PDB
!      ,mdin_pdb &
!#endif
; ! line terminator for free-form version

common/mdin_flags/mdin_pb &
!#ifdef MDM_MD
!      ,mdin_mdm &
!#endif
!#ifdef MDM_PDB
!      ,mdin_pdb &
!#endif
; ! line terminator for free-form version

integer BC_HULP  ! size in integers of common HULP
parameter ( BC_HULP = 9 )

integer     ntpr,ntwr,ntwx,ntwv,ntwe,ntpp,ioutfm,ntwprt,ntave &
!#ifdef MDM_MD
!      ,mdm_nstep &
!#endif
; ! line terminator for free-form version
common/hulp/ntpr,ntwr,ntwx,ntwv,ntwe,ntpp,ioutfm,ntwprt,ntave &
!#ifdef MDM_MD
!      ! this variable is last in the common and does not increment BC_HULP
!      ,mdm_nstep &
!#endif
; ! line terminator for free-form version

!      NMRRDR : Contains information about input/output file redirection
!               REDIR and IREDIR contain information regarding
!               LISTIN, LISTOUT, READNMR, NOESY, SHIFTS, DUMPAVE,
!               PCSHIFT and DIPOLE respectively. If IREDIR(I) > 0,
!               then that input/output has been redirected.

integer iredir(8)
common/nmrrdr/redir,iredir
