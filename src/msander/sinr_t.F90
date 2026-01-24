#include "../include/dprec.fh"
!-----------------------------------------------------------------------------
! Stochastic Isokinetic Nose-Hoover RESPA (SINR) Integrator Implementation
! By Jason Deckman, Mark Tuckerman, David Case
! January 2016
!-----------------------------------------------------------------------------
module sinr_t
  implicit none

  private

  type sinr
     sequence
   ! Stochastic Isokinetic Nose-Hoover chain RESPA integrator (SIN(R))

     integer :: L,nsys,nres,natom,istep,nrespa
#ifdef MPI
     integer :: mpirank,mpicomm,nproc,idummy
#endif

     double precision :: kb,errtol,wj(3)
     double precision :: dt,dt2,T,kbT,sigma
     double precision :: egt,sigsqe2gt
     double precision :: Q1,Q2,tau,gamma_n
     double precision :: lbeta,LLp1

   ! Arrays: auxilliary thermostat variables
     double precision, allocatable :: v1(:,:)
     double precision, allocatable :: v2(:,:)

  end type sinr

public sinr, sinr_init, init_sinr_vels, iLndt, iLvdt, iLudt, sinr_cleanup, &
       sinr_temp, sinr_write_vels, sinr_read_vels

#ifdef MPI
public sinr_mpi_init, sinr_test
#endif

contains

#ifdef MPI
   subroutine sinr_init(natom,nkija,dt,boltz,temp0,gamma_ln,tau,sd,mpic)
#else
   subroutine sinr_init(natom,nkija,dt,boltz,temp0,gamma_ln,tau,sd)
#endif
#ifdef MPI
    use mpi
#endif
    implicit none
     type(sinr), intent(inout) :: sd

     integer :: natom,nkija
#ifdef MPI
     integer :: ierr,mpic
#endif

     double precision :: dt,boltz,temp0,gamma_ln,tau

   ! MPI initialization

#ifdef MPI
     sd%mpicomm = mpic
     sd%mpirank = 0
     sd%nproc = 1
     call mpi_comm_rank(sd%mpicomm,sd%mpirank,ierr)
#endif

   ! Initialize SINR vairables

     sd%natom = natom
     sd%istep = 0
     sd%errtol = 0.00001

     sd%dt = dt
     sd%dt2 = dt/2.0d0

     sd%L = nkija
     sd%gamma_n = gamma_ln
     sd%kb = 2.0d0*boltz
     sd%kbT = sd%kb*temp0
     sd%T = temp0

     sd%lbeta = sd%L*sd%kbT
     sd%LLp1 = dble(sd%L)/dble(sd%L+1)

   ! Inititialize thermostat masses

     sd%Q1 = sd%kbT*tau*tau
     sd%Q2 = sd%Q1

   ! Stochastic operator initialization

     sd%sigma = sqrt(2.0d0*sd%kbT*gamma_ln/sd%Q2)
     sd%egt = exp(-1.0d0*gamma_ln*dt)
     sd%sigsqe2gt = sd%sigma*sqrt((1.0d0 - exp(-2.0d0*gamma_ln*dt))/(2.0d0*gamma_ln))

   ! Array initialization

     allocate(sd%v1(nkija,3*natom), sd%v2(nkija,3*natom))

   ! Suzuki-Yoshida factorisation operators initialization
   ! P. 3587 Tuckerman et. al., Mol. Phys. 111 (2013)

     sd%nsys = 3
     sd%nres = 2

     sd%wj(1) = (1.0d0/(2.0d0 - 2.0d0**(1.0d0/3.0d0)))*sd%dt/dble(sd%nres)
     sd%wj(3) = sd%wj(1)
     sd%wj(2) = -2.0d0**(1.0d0/3.0d0)*sd%wj(1)

   ! Write out important variables:

#ifdef MPI
     if(sd%mpirank==0) then
#endif
        write(6,*) ""
        write(6,'(A74)') "| Using Stochastic Isokinetic Nose-Hoover RESPA (SINR) integrator (ntt=10)"
        write(6,'(A74)') "| ------------------------------------------------------------------------"
        write(6,*) ""
        write(6,'(A71)') "| NOTE: Only the coordinates are canonical while the velocites are NOT."
        write(6,'(A71)') "| The reported temperature will thus appear anomolous, being about half"
        write(6,'(A69)') "| the desired simulation temperature for 1 thermostat DOF (nkija = 1)"
        write(6,'(A73)') "| and will approach but not exceed the set simulation temperature, temp0."
        write(6,'(A68)') "| However the coordinates are canonical and represent configurations"
        write(6,'(A74)') "| sampled from a Boltzman distribution at the specfied temperature, temp0."
        write(6,'(A73)') "| See SINR related references in the AMBER manual for a full explanation."

        write(6,*) ""

        write(6,"(A51,I2)") "| Number of SINR thermostat chain variables (DOF): ", sd%L
        write(6,"(A22,F12.3)") "| Thermostat mass Q1: ", sd%Q1
        write(6,*)
#ifdef MPI
     endif
#endif

   end subroutine sinr_init

   subroutine init_sinr_vels(v,m,sd)
     use random
#ifdef MPI
     use mpi
#endif
     implicit none
     type(sinr), intent(inout) :: sd

     integer :: ind,i,j,k

     double precision :: v(3*sd%natom),m(sd%natom)
     double precision :: vsq

     ind = 0

#ifdef MPI
     if(sd%mpirank==0) then
#endif
        do i=1,sd%natom
           do j=1,3

              call gauss(0.0d0,1.0d0,v(ind+j))

              do k=1,sd%L
                 call gauss(0.0d0,1.0d0,sd%v1(k,ind+j))
                 call gauss(0.0d0,1.0d0,sd%v2(k,ind+j))
                 sd%v2(k,ind+1) = sd%v2(k,ind+1)*sqrt(sd%kbT/sd%Q2)
              enddo

              vsq = v(ind+j)*v(ind+j)

              do k=1,sd%L
                 vsq = vsq + sd%v1(k,ind+j)*sd%v1(k,ind+j)
              enddo

              vsq = sqrt(vsq)

              v(ind+j) = v(ind+j)/vsq

              do k=1,sd%L
                 sd%v1(k,ind+j) = sd%v1(k,ind+j)/vsq
              enddo

              v(ind+j) = v(ind+j)/(sqrt(m(i)/(dble(sd%L)*sd%kbT)))

              do k=1,sd%L
                 sd%v1(k,ind+j) = sd%v1(k,ind+j)/(sqrt(sd%Q1/(dble(sd%L+1)*sd%kbT)))
              enddo

           enddo
           ind=ind+3
        enddo
#ifdef MPI
     endif
#endif
   end subroutine init_sinr_vels

#ifdef MPI
   subroutine sinr_mpi_init(v,ntask,sd)
     use mpi
     implicit none

     type(sinr), intent(inout) :: sd

     integer :: ntask,ierr
     double precision :: v(3*sd%natom)

     call mpi_bcast(v,3*sd%natom,MPI_DOUBLE_PRECISION, &
                    0,sd%mpicomm,ierr)

     call mpi_bcast(sd%v1,sd%L*3*sd%natom,MPI_DOUBLE_PRECISION, &
                    0,sd%mpicomm,ierr)

     call mpi_bcast(sd%v2,sd%L*3*sd%natom,MPI_DOUBLE_PRECISION, &
                    0,sd%mpicomm,ierr)

     sd%nproc = ntask

   end subroutine sinr_mpi_init
#endif

#ifdef MPI
   subroutine sinr_test(sd)
     implicit none
     type(sinr), intent(in) :: sd

     integer i,j,unit

     unit = 1173

     do i=1,sd%natom
        write(unit+sd%mpirank,*) sd%v1(1,3*(i-1)+1), sd%v1(1,3*(i-1)+2), &
                                 sd%v1(1,3*(i-1)+3)
     enddo

   end subroutine
#endif

   subroutine iLndt(v,m,istart,iend,sd)
     implicit none

     type(sinr), intent(inout) :: sd

     integer istart,iend
     integer :: ind,i,j,k,ii,jj

     double precision :: v(3*sd%natom),m(sd%natom)
     double precision :: Htnres,q1v1sq,v2kw(sd%L)

     ind = 3*(istart-1)

     do ii=istart,iend !,sd%natom
        do jj=1,3

           do i=1,sd%nsys
              do j=1,sd%nres

               ! Outer translation operator LN,2 eq. 66

                 do k=1,sd%L
                    sd%v2(k,ind+jj) = sd%v2(k,ind+jj) + 0.25d0*sd%wj(i)*(sd%Q1*(sd%v1(k,ind+jj)**2) - sd%kbT)/sd%Q2
                 enddo

               ! Inner translation operator LN,1 eq. 66, solutions with eq. 70

                 q1v1sq = 0.0d0

                 do k=1,sd%L
                    v2kw(k) = exp(-1.0d0*sd%v2(k,ind+jj)*sd%wj(i))
                    q1v1sq = q1v1sq + sd%Q1*sd%v1(k,ind+jj)*sd%v1(k,ind+jj)*v2kw(k)
                 enddo

                 Htnres = sqrt(sd%lbeta/(m(ii)*v(ind+jj)**2 + sd%LLp1*q1v1sq))

                 v(ind+jj) = v(ind+jj)*Htnres

                 do k=1,sd%L
                    sd%v1(k,ind+jj) = sd%v1(k,ind+jj)*Htnres*sqrt(v2kw(k))
                 enddo

               ! LN,2 again

                 do k=1,sd%L
                    sd%v2(k,ind+jj) = sd%v2(k,ind+jj) + 0.25d0*sd%wj(i)*(sd%Q1*(sd%v1(k,ind+jj)**2) - sd%kbT)/sd%Q2
                 enddo

              enddo
            enddo

         enddo
         ind=ind+3
     enddo

   end subroutine iLndt

   subroutine iLvdt(v,f,m,istart,iend,sd)
     implicit none
     type(sinr), intent(inout) :: sd
     integer istart,iend
     integer ind,i,j,k
     double precision :: v(3*sd%natom),f(3*sd%natom),m(sd%natom)
     double precision :: sinhbt,coshbt,fs,sdot,sbt,p,mlbeta,smlbeta

     ind = 3*(istart-1)

     do i=istart,iend !sd%natom
        do j=1,3

           ! Updated equations from 10.1021/acs.jpca.9b02771 and 10.1021/acs.jpcb.2c02262

           p = m(i)*v(ind+j)
           mlbeta = m(i)*(sd%lbeta)
           smlbeta = sqrt(mlbeta)

           sbt = (f(ind+j)/smlbeta)*(sd%dt2)
           sinhbt = sinh(sbt)
           coshbt = cosh(sbt)

           fs = smlbeta*sinhbt + p*(coshbt - 1.0d0)
           sdot = (p/smlbeta)*sinhbt + coshbt
 
           v(ind+j) = (v(ind+j) + (fs/m(i)))/sdot

           ! Update thermostats

           do k=1,sd%L
              sd%v1(k,ind+j) = sd%v1(k,ind+j)/sdot
           enddo

        enddo
        ind=ind+3
     enddo

   end subroutine iLvdt

   subroutine iLudt(q,v,istart,iend,sd)
     use random
     implicit none
     type(sinr), intent(inout) :: sd
     integer istart,iend
     integer ind,i,j,k
     double precision :: q(3*sd%natom),v(3*sd%natom),rnum

     ind = 3*(istart-1)

     do i=istart,iend
        do j=1,3

           q(ind+j) = q(ind+j) + sd%dt*v(ind+j)

           do k=1,sd%L
              call gauss(0.0d0,1.0d0,rnum)
              sd%v2(k,ind+j) = sd%v2(k,ind+j)*sd%egt + rnum*sd%sigsqe2gt
           enddo
        enddo
        ind=ind+3
     enddo
   end subroutine iLudt

   subroutine sinr_temp(v,m,istart,iend,sd)
#ifdef MPI
     use mpi
#endif
     implicit none
     type(sinr), intent(in) :: sd

     integer :: istart,iend,ind,i,j,k
     double precision :: v(3*sd%natom),m(sd%natom)
     double precision :: ke, ketot, ct
#ifdef MPI
     double precision :: vt(3*sd%natom), vtt(3*sd%natom)
     double precision :: v1t(sd%L,3*sd%natom), v1tt(sd%L,3*sd%natom)
     integer :: ierr

     if(sd%nproc .gt. 1) then

        vt=0.0d0
        vtt=0.0d0
        v1t=0.0
        v1tt=0.0

        vt(istart:iend) = v(istart:iend)

        do k=1,sd%L
          v1t(k,istart:iend) = sd%v1(k,istart:iend)
        enddo

        call mpi_reduce(vt,vtt,3*sd%natom,MPI_DOUBLE_PRECISION,mpi_sum,0, &
                   sd%mpicomm,ierr)

        call mpi_reduce(v1t,v1tt,sd%L*3*sd%natom,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                  sd%mpicomm,ierr)

     else if(sd%mpirank==0) then
         vtt = v
         v1tt = sd%v1
     endif

     if(sd%mpirank == 0) then

        ketot = 0.0d0
        ct = 0.0d0

        ind = 0

        do i=1,sd%natom
           do j=1,3
              ke = m(i)*vtt(ind+j)**2
              do k=1,sd%L
                 ke = ke + sd%LLp1*sd%Q1*(v1tt(k,ind+j)**2)
              enddo
              ketot = ketot + ke
            ! ct = ct + ke/(sd%L+sd%kb)
           enddo
           ind=ind+3
        enddo

        ! write(6,*) ""
        ! write(6,"(A,F12.3)") "| Total SINR kinetic energy: ", ketot
        ! write(6,"(A,F12.5)") "| SINR kinetic energy / N: ", ketot/sd%natom
        ! write(6,"(A,F12.5)") "| Isokinetic temperature: ", ketot/sd%natom/(sd%L+sd%kb)
        !  write(6,*) ""
    endif
#else
        ketot = 0.0d0
        ct = 0.0d0

        ind = 0

        do i=1,sd%natom
           do j=1,3
              ke = m(i)*v(ind+j)**2
              do k=1,sd%L
                 ke = ke + sd%LLp1*sd%Q1*(sd%v1(k,ind+j)**2)
              enddo
              ketot = ketot + ke
              ct = ct + ke/(sd%L+sd%kb)
           enddo
           ind=ind+3
        enddo

          ! write(6,*) ""
          ! write(6,"(A,F12.3)") "| Total SINR kinetic energy: ", ketot
          ! Write(6,"(A,F12.5)") "| SINR kinetic energy / N: ", ketot/sd%natom
          ! write(6,"(A,F12.5)") "| Isokinetic temperature: ", ct/sd%natom
          ! write(6,*) ""
#endif
   end subroutine sinr_temp

   subroutine sinr_write_vels(v,nsteps,istart,iend,tau,sd)
#ifdef MPI
     use mpi
#endif
     implicit none
     type(sinr), intent(in) :: sd
     integer :: ind,i,k,istart,iend,nsteps
     double precision :: tau,v(3*sd%natom)
#ifdef MPI
     double precision :: blksize,vt(3*sd%natom), vtt(3*sd%natom)
     double precision :: v1t(sd%L,3*sd%natom), v1tt(sd%L,3*sd%natom)
     double precision :: v2t(sd%L,3*sd%natom), v2tt(sd%L,3*sd%natom)
     integer :: ierr

      if(sd%nproc .gt. 1) then

        vt=0.0d0
        vtt=0.0d0

        v1t=0.0
        v1tt=0.0

        v2t=0.0
        v2tt=0.0

        vt(istart:iend) = v(istart:iend)

        do k=1,sd%L
          v1t(k,istart:iend) = sd%v1(k,istart:iend)
          v2t(k,istart:iend) = sd%v2(k,istart:iend)
        enddo

        call mpi_reduce(vt,vtt,3*sd%natom,MPI_DOUBLE_PRECISION,mpi_sum,0, &
                   sd%mpicomm,ierr)

        call mpi_reduce(v1t,v1tt,sd%L*3*sd%natom,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                  sd%mpicomm,ierr)

        call mpi_reduce(v2t,v2tt,sd%L*3*sd%natom,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                  sd%mpicomm,ierr)

     else if(sd%mpirank==0) then
        vtt = v
        v1tt = sd%v1
        v2tt = sd%v2
     endif

     if(sd%mpirank == 0) then

        open(579,file="sinrvels.rst")

        ind=0

        write(579,*) "NATOMS: ", sd%natom, "      STEPS: ", nsteps

        do i=1,sd%natom
           write(579,*) vtt(ind+1), vtt(ind+2), vtt(ind+3)
           ind=ind+3
        enddo

        write(579,*) "V1", sd%L, tau

        do k=1,sd%L
           ind=0
           do i=1,sd%natom
              write(579,*) v1tt(k,ind+1), v1tt(k,ind+2), v1tt(k,ind+3)
              ind=ind+3
           enddo
        enddo

        write(579,*) "V2", sd%L

        ind=0

        do k=1,sd%L
           ind=0
           do i=1,sd%natom
              write(579,*) v2tt(k,ind+1), v2tt(k,ind+2), v2tt(k,ind+3)
              ind=ind+3
           enddo
        enddo

        call flush(579)
        close(579)

     endif
#else
     open(579,file="sinrvels.rst")

     ind=0

     write(579,*) "NATOMS: ", sd%natom, "      STEPS: ", nsteps

     do i=1,sd%natom
        write(579,*) v(ind+1), v(ind+2), v(ind+3)
        ind=ind+3
     enddo

     write(579,*) "V1", sd%L, tau

     do k=1,sd%L
        ind=0
        do i=1,sd%natom
           write(579,*) sd%v1(k,ind+1), sd%v1(k,ind+2), sd%v1(k,ind+3)
           ind=ind+3
        enddo
     enddo

     write(579,*) "V2", sd%L

     ind=0

     do k=1,sd%L
        ind=0
        do i=1,sd%natom
           write(579,*) sd%v2(k,ind+1), sd%v2(k,ind+2),sd%v2(k,ind+3)
           ind=ind+3
        enddo
     enddo

     call flush(579)
     close(579)
#endif

   end subroutine sinr_write_vels

   subroutine sinr_read_vels(v,tau,sd)
     implicit none
     type(sinr), intent(inout) :: sd
     integer :: ind,i,k,nL
     logical :: exst
     double precision :: tau,taus,v(3*sd%natom)
     character(len=80) :: dummy

#ifdef MPI
     if(sd%mpirank==0) then
#endif
        INQUIRE(FILE='sinrvels.rst',EXIST=exst)

        if(.not.exst) then
           write(6,'(A67)') "| Cannot restart run: SINR restart file 'sinrvels.rst' not present."
           write(6,*) ""
           call flush(6)
           call mexit(6, 1)
        endif

        open(579,file="sinrvels.rst")

        ind=0

        read(579,*) dummy

        do i=1,sd%natom
           read(579,*) v(ind+1), v(ind+2), v(ind+3)
           ind=ind+3
        enddo

        read(579,*) dummy, nL, taus

        if(nL .ne. sd%L) then
            write(6,'(A72)') "| Number of specified chains does not match number in SINR restart file."
            write(6,'(A52,I2,A13)') "| Change 'nkija' to value saved in file, which is: '", nL, "' to restart."
            write(6,'(A10)') "| Exiting."
            write(6,*) ""
            call flush(6)
            call mexit(6, 1)
        endif

        if(taus .ne. tau) then
            write(6,'(A74)') "| Cannot restart: saved thermostat mass scale parameter 'sinrtau' does not"
            write(6,'(A78)') "| match that specified in input file - change sinrtau to match that in restart"
            write(6,'(A18,F12.5)') "| file, which is: ", taus
            write(6,'(A10)') "| Exiting."
            write(6,*) ""
            call flush(6)
            call mexit(6, 1)
        endif

        do k=1,sd%L
           ind=0
           do i=1,sd%natom
              read(579,*) sd%v1(k,ind+1), sd%v1(k,ind+2), sd%v1(k,ind+3)
              ind=ind+3
           enddo
        enddo

        read(579,*) dummy

        ind=0

        do k=1,sd%L
           ind=0
           do i=1,sd%natom
              read(579,*) sd%v2(k,ind+1), sd%v2(k,ind+2),sd%v2(k,ind+3)
              ind=ind+3
           enddo
        enddo

        close(579)

#ifdef MPI
     endif
#endif
   end subroutine sinr_read_vels

   subroutine sinr_cleanup(sd)
     implicit none
     type(sinr) :: sd

     deallocate(sd%v1,sd%v2)

   end subroutine sinr_cleanup

end module sinr_t
