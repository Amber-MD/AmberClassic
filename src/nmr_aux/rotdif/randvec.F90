      program randvec
      implicit none

      integer npts,seed
      double precision theta,phi,x,y,z
      double precision random
      double precision pi
      parameter (pi=3.141592653589793d0)
      integer i,iarg,last_arg
      character(len=20) arg

!     defaults:
      npts = 1000
      seed = 80531

!     Process command-line arguments:
      iarg = 0
      last_arg = command_argument_count()
      do while (iarg < last_arg)
         iarg = iarg + 1
         call getarg(iarg,arg)
         if( arg == '-n') then
            iarg = iarg + 1
            call getarg(iarg,arg)
            read(arg,*) npts
            write(6,*) 'setting npts to ',npts
         else if( arg == '-s') then
            iarg = iarg + 1
            call getarg(iarg,arg)
            read(arg,*) seed
            write(6,*) 'setting seed to ',seed
         else 
            write(6,*) 'bad argument: ', arg
            call exit(1)
         end if
      end do

!     generate a set of points uniformly covering surface of unit sphere

!     u,v: uniformly distributed random deviates 0 < u, v < 1
!     azimuthal angle phi: rho(phi)=1/2*pi -> phi=2*pi*u 0 < phi < 2*pi
!     polar angle theta: rho(theta)=sin(theta) 
!                        -> theta=arccos(1-v) 0 < theta < 0.5*pi
!     only a single hemisphere is required
!     rho(theta)=0.5d0*sin(theta) -> theta=arccos(1-2*v) 0 < theta < pi
!     if both hemispheres are required
!     see Numerical Recipes sec. 7.2
!     http://mathworld.wolfram.com/SpherePoi.html
      
      do i=1,npts
         phi=2d0*pi*random(seed)
         theta=dacos(1d0-random(seed))
         z=dcos(theta)
         x=dsin(theta)*dcos(phi)
         y=dsin(theta)*dsin(phi)
         write(6,'(i6,2x,3(f15.8,2x))') i,x,y,z
      end do

      stop
      end program randvec


      double precision function random(seed)
      implicit none
      real*4 mbig,mseed,mz,fac
      parameter(mbig=4000000.,mseed=1618033.,mz=0,fac=1./mbig)
      real*4 ma(55),mj,mk
      integer iff,i,ii,k,inext,seed,inextp
      save iff,inext,inextp,ma
      data iff /0/
      if(seed.lt.0.or.iff.eq.0) then
        iff=1
        mj=mseed-iabs(seed)
        mj=amod(mj,mbig)
        ma(55)=mj
        mk=1
        do i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.mz) mk=mk+mbig
          mj=ma(ii)
        enddo
        do k=1,4
          do i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.mz) ma(i)=ma(i)+mbig
          enddo
        enddo
        inext=0
        inextp=31
        seed=iabs(seed)
      endif
      inext=inext+1
      if(inext.eq.56) inext=1
      inextp=inextp+1
      if(inextp.eq.56) inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz) mj=mj+mbig
      ma(inext)=mj
      random=mj*fac
      return
      end function random
