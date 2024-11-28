
  ! Call DEGCNT to get the actual number of degrees of freedom for the
  ! solute and solvent. This call returns the correct numbers for belly
  ! simulations and simulations with separate solute/solvent scaling -- dap
  ! "IDUMAR" is dummy array. Used since this routine was also used w/ GIBBS.
  
  call degcnt(ibelly, nr, ix(ibellygp), nsolut, nbonh, nbona, 0, ix(iibh), &
              ix(ijbh), ix(iiba), ix(ijba), idumar, idumar, ntc, idumar, 0, &
              0, 0, idumar, rndfp, rndfs)

  ! RNDFP = # degrees of freedom for solute
  ! RNDFS = # degrees of freedom for solvent
  ! RNDF = total number of degrees of freedom.
#  ifdef MPI
  if (mdloop .eq. 0 .and. master) then
#  else
  if (master) then
#  endif /* MPI */
    write (6,'(a,f8.0)') "|  # of SOLUTE  degrees of freedom (RNDFP): ", &
        rndfp
    write (6,'(a,f8.0)') "|  # of SOLVENT degrees of freedom (RNDFS): ", &
        rndfs
  end if

  ! Modify RNDFP to reflect NDFMIN (set in mdread) and num_noshake
  rndfp = rndfp - ndfmin + num_noshake
  rndf = rndfp + rndfs
#  ifdef MPI
  if (mdloop .eq. 0 .and. master) then
#  else
  if (master) then
#  endif /* MPI */
    write (6,'(a,f8.0,a,i6,a,f8.0)') "|  NDFMIN = ", rndfp, &
      "     NUM_NOSHAKE = ", num_noshake, "     CORRECTED RNDFP = ", rndfp
    write (6,'(a,f8.0)') "|  TOTAL # of degrees of freedom (RNDF) = ", rndf
  end if
  ! Correct the degree of freedom count for extra points.  Note
  ! that rndfp, rndfs are uncorrected in an extra points context!
  call fix_degree_count(rndf)  
  if (isgld > 0) then

    ! Number of degrees of freedom in the SGLD part
    if (isgsta == 1) then
      sgsta_rndfp = 0
    else
      call degcnt(ibelly, nr, ix(ibellygp), isgsta-1, nbonh, nbona, 0, &
                  ix(iibh), ix(ijbh), ix(iiba), ix(ijba), idumar, idumar, &
                  ntc, idumar, 0, 0, 0, idumar, sgsta_rndfp, ignore_solvent)
    end if
    if (isgend == nr) then
      sgend_rndfp = rndf
    else
      call degcnt(ibelly, nr, ix(ibellygp), isgend, nbonh, nbona, 0, &
                  ix(iibh), ix(ijbh), ix(iiba), ix(ijba), idumar, &
                  idumar, ntc, idumar, 0, 0, 0, idumar, sgend_rndfp, &
                  ignore_solvent)
    end if

    ! Note that the solute ndf outputs above from degcnt are uncorrected
    ! also, ndfmin is not always being handled.
    call sg_fix_degree_count(sgsta_rndfp, sgend_rndfp, ndfmin, rndf)
  end if

#ifdef MPI
  ! Soft-core TI: count additional degrees of freedom.
  if (ifsc /= 0) then
    call sc_degrees_o_freedom(ndfmin)
  end if
#endif /* MPI */
  ! End of degrees of freedom setup.
