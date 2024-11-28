      program RMScorr
      implicit none

!     given A) an initial normalized vector B) a series of rotation
!     matrices determined at each time step from RMS fitting the MD
!     snapshot to a reference structure:
!     1) rotate the vector using the given matrices
!     2) compute the time correlation function for the vector
!     3) iteratively estimate the area under the correlation function curve

      integer maxdat,ndat
      parameter (maxdat=100000)
      integer ncorr,itotframes
      real*8 x(maxdat),y(maxdat),z(maxdat),tdat(maxdat)
      real*8 p2(maxdat),p1(maxdat)
      integer ivec,ios
      real*8 x0,y0,z0,rotmat(3,3)
      character(len=80) rotmats,vecs
      common/ dat/ tdat,p2,ndat

      real*8 tfac

      integer i,j,irot,l,itmax
      real*8 ti,tf,delmin,d0,deff

!     ncorr:  length of time correlation function

      read(5,*) ncorr,tfac
      read(5,*) ti,tf
      read(5,*) itmax,delmin,d0,l
      read(5,*) rotmats
      read(5,*) vecs

      do i=1,ncorr+1
         tdat(i) = tfac*dble(i-1)
      end do
      ndat=ncorr+1

      open( unit=3,file=rotmats,status='OLD',iostat=ios)
      open( unit=2,file=vecs,status='OLD',iostat=ios)

      do j=1,999999
         read(2,*,end=99) ivec,x0,y0,z0
         x(1)=x0
         y(1)=y0
         z(1)=z0
         rewind(3)
         do i=1,999999
            if( i > maxdat ) stop 'too many rotation matrices'
            read(3,*,end=98) irot,rotmat(1,1),rotmat(2,1),rotmat(3,1), &
                        rotmat(1,2),rotmat(2,2),rotmat(3,2), &
                        rotmat(1,3),rotmat(2,3),rotmat(3,3)

            x(i+1)=rotmat(1,1)*x0+rotmat(2,1)*y0+rotmat(3,1)*z0
            y(i+1)=rotmat(1,2)*x0+rotmat(2,2)*y0+rotmat(3,2)*z0
            z(i+1)=rotmat(1,3)*x0+rotmat(2,3)*y0+rotmat(3,3)*z0

         end do
   98    itotframes = i-1

         call compute_corr(x,y,z,ncorr,itotframes,maxdat,p2,p1)
!        do i=1,ndat
!           write(6,*) tdat(i),p2(i)
!        end do

         call dlocint(ti,tf,itmax,delmin,d0,l,deff)
         write(6,'(i6,f15.8)') ivec,deff

      end do
   99 stop
      end


      subroutine compute_corr(x,y,z,ncorr,itotframes,maxdat,p2,p1)
      implicit none

!     x,y,z: arrays contain coordinates of (unnormalized) vector 
!     ncorr: maximum length to compute time correlation functions -
!      units of 'frames'
!     itotframes:  total number of frames provided

      integer ncorr,itotframes,maxdat
      real*8 x(maxdat),y(maxdat),z(maxdat)
      real*8 p2(maxdat),p1(maxdat)
   
      integer jmax,j,k
      real*8 mag,xj,yj,zj,xk,yk,zk,dot

      integer i

      do i=1,ncorr+1
         p2(i)=0d0
         p1(i)=0d0
      end do

!     i loop:  each value of i is a value of delay (correlation
!     function argument)
      do i=0,ncorr
         jmax=itotframes+1-i
         do j=1,jmax
            mag=dsqrt(x(j)*x(j)+y(j)*y(j)+z(j)*z(j))
            xj=x(j)/mag
            yj=y(j)/mag
            zj=z(j)/mag
            k=j+i
            mag=dsqrt(x(k)*x(k)+y(k)*y(k)+z(k)*z(k))
            xk=x(k)/mag
            yk=y(k)/mag
            zk=z(k)/mag
            dot=xj*xk+yj*yk+zj*zk
            p2(i+1)=p2(i+1)+1.5d0*dot*dot-0.5d0
            p1(i+1)=p1(i+1)+dot
         end do
         p2(i+1)=p2(i+1)/dble(jmax)
         p1(i+1)=p1(i+1)/dble(jmax)
      end do

      return
      end
      subroutine dlocint(ti,tf,itmax,delmin,d0,l,deff)
      implicit none

!     computes effect diffusion constant for a vector using its
!     correlation function as input
!     starting with definition 6*D=integral[0,inf;C(t)] 
!     integrates C(t) from ti -> tf yielding F(ti,tf)
!     iteratively solves the equation 
!     D(i+1)=[exp(6*D(i)*ti)-exp(6*D(i)*tf)]/[6*F(ti,tf)]
!     (numerator obtained by integrating exp(6*D*t) from ti -> tf)

!     modified so that itsolv now solves
!     F(ti,tf;C(t)]=integral[ti,tf;C(t)]
!     D(i+1)={exp[l*(l+1)*D(i)*ti]-exp[l*(l+1)*D(i)*tf)]}/[l*(l+1)*F(ti,tf)]

      real*8 ti,tf,dydx1,dydxn
      integer ndat,nmax
      parameter (nmax=100000)
      integer l
      real*8 tdat(nmax),ctdat(nmax),deriv2(nmax)
      real*8 sumct
      common/ dat/ tdat,ctdat,ndat
      common/deriv/ deriv2,dydx1,dydxn
      
      integer itmax,info
      real*8 delmin,d0
      real*8 di,deff

      integer i


!     ti,tf: integration limits
!     dydx1,dydxn:  (estimated) first derivatives of function to be
!                   interpolated (C(t)) by spline/splint; accurate
!                   estimate not needed
!     ndat:  # C(t) data points to be used for interpolation
!     itmax:  maximum number of iterations in subroutine itsolv
!     delmin:  convergence criterion used in subroutine itsolv;
!              maximum accepted fractional change in successive 
!              iterations
!     d0: initial guess for diffusion constant; accurate estimate not
!         needed
!     info:  =1; write out data on convergence of iterative solver
!     l: order of Legendre polynomial in the correlation function
!        <P(l)>

      dydx1 = 0.d0
      dydxn = 0.d0

      call intct(ti,tf,sumct)
      di=d0
      call itsolv(itmax,delmin,l,di,ti,tf,sumct,deff,info)

      return
      end subroutine dlocint

      subroutine intct(ti,tf,sumct)
      implicit none

!     Integrates the MD generated C(t) from ti->tf.
!     Calls subroutine qromb found in 4.3 of Numerical Recipes
!     Interpolates values of C(t) using modified versions of subroutines
!     spline and splint found in 3.3 of Numerical Recipes
!     C(t) from MD simulation provides the table of data needed for 
!     interpolation (to be read in in main program and placed in 
!     common block)

      real*8 ti,tf,sumct,ct

      integer nmax,ndat
      parameter (nmax=100000)
      real*8 tdat(nmax),ctdat(nmax)
      real*8 dydx1,dydxn,deriv2(nmax)
      common/ dat/ tdat,ctdat,ndat
      common/ deriv/ deriv2,dydx1,dydxn

      external ct

!     subroutine spline computes 2nd derivatives needed by splint
      call spline(tdat,ctdat,ndat,dydx1,dydxn,deriv2)
      call qromb(ct,ti,tf,sumct)

      return
      end subroutine intct

      subroutine qromb(func,a,b,ss)
      implicit none

!     Romberg integration routine found in 4.3 of Numerical Recipes

      integer jmax,jmaxp,k,km
      real*8 a,b,func,ss,eps
      external func
      parameter (eps=1d-6,jmax=20,jmaxp=jmax+1,k=5,km=k-1)

      integer j
      real*8 dss,h(jmaxp),s(jmaxp)

!     uses polint,trapzd

      h(1)=1d0
      do j=1,jmax
         call trapzd(func,a,b,s(j),j)
         if(j>=k)then
           call polint(h(j-km),s(j-km),k,0d0,ss,dss)
           if(dabs(dss)<=eps*dabs(ss)) return
         end if
         s(j+1)=s(j)
         h(j+1)=0.25d0*h(j)
      end do
      stop 'too many steps in qromb'

      end subroutine qromb

      subroutine trapzd(func,a,b,s,n)
      implicit none

!     trapezoid integration routine found in 4.2 of Numerical Recipes

      integer n
      real*8 a,b,s,func
      external func

      integer it,j
      real*8 del,sum,tnm,x

      if(n==1)then
        s=0.5d0*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5d0*del
        sum=0d0
        do j=1,it
           sum=sum+func(x)
           x=x+del
        end do
        s=0.5d0*(s+(b-a)*sum/tnm)
      end if

      return
      end subroutine trapzd

      subroutine polint(xa,ya,n,x,y,dy)
      implicit none

!     polynomial interpolation routine found in 3.1 of Numerical Recipes
!     modified to redimension data tables xa, ya
!     previously, only c, d were of dimension nmax

      integer n,nmax
      parameter (nmax=100000)
      real*8 dy,x,y,xa(nmax),ya(nmax)
!     real*8 dy,x,y,xa(n),ya(n)
!     parameter(nmax=100000)

      integer i,m,ns
      real*8 den,dif,dift,ho,hp,w,c(nmax),d(nmax)

      ns=1
      dif=dabs(x-xa(1))
      do i=1,n
         dift=dabs(x-xa(1))
         if(dift<dif)then
           ns=i
           dif=dift
         end if
         c(i)=ya(i)
         d(i)=ya(i)
      end do
      y=ya(ns)
      ns=ns-1
      do m=1,n-1
         do i=1,n-m
            ho=xa(i)-x
            hp=xa(i+m)-x
            w=c(i+1)-d(i)
            den=ho-hp
            if(den==0d0) stop 'failure in polint'
            den=w/den
            d(i)=hp*den
            c(i)=ho*den
         end do
         if(2*ns<n-m)then
           dy=c(ns+1)
         else
           dy=d(ns)
           ns=ns-1
         end if
         y=y+dy
      end do

      return
      end subroutine polint

      real*8 function ct(t)
      implicit none

!     evaluates correlation function C(t) at any value t by 
!     interpolating between simulated data from MD
!     using cubic spline interpolation
      
      real*8 t

      integer nmax,ndat
      parameter (nmax=100000)
      real*8 tdat(nmax),ctdat(nmax)
      real*8 dydx1,dydxn,deriv2(nmax)
      common/ dat/ tdat,ctdat,ndat
      common/ deriv/ deriv2,dydx1,dydxn

      call splint(tdat,ctdat,deriv2,ndat,t,ct)

      return
      end function ct

      subroutine splint(xa,ya,y2a,n,x,y)
      implicit none

!     cubic spline interpolation routine found in 3.3 of Numerical
!     Recipes
!     uses output y2a from spline
!     modified to redimension tables xa, y2a, ya
!     previously, xa, y2a, ya were of dimension n read in as argument

      integer n,nmax
      parameter (nmax=100000)
      real*8 x,y,xa(nmax),y2a(nmax),ya(nmax)
!     integer n
!     real*8 x,y,xa(n),y2a(n),ya(n)

      integer k,khi,klo
      real*8 a,b,h

      klo=1
      khi=n
1     if(khi-klo>1)then
        k=(khi+klo)/2
        if(xa(k)>x)then
          khi=k
        else
          klo=k
        end if
        go to 1
      end if
      h=xa(khi)-xa(klo)
      if(h==0d0) stop 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+  &
        ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6d0

      return
      end subroutine splint


      subroutine spline(x,y,n,yp1,ypn,y2)
      implicit none
!
!     given arrays x, y defining function y(x), and first derivatives
!     at x(1) and x(n) yp1, ypn, returns array containing second
!     derivatives of interpolating function needed by subroutine splint
!     which produces a cubic spline interpolation of y
!     from 3.3 in Numerical Recipes
!     note:  spline is called only once for a given function y(x);
!     interpolation routine splint is called for each desired value
!     of x using output y2
!     modified to redimension x, y, y2
!     previously, only u of dimension nmax

      integer n,nmax
      parameter (nmax=100000)
      real*8 yp1,ypn,x(nmax),y(nmax),y2(nmax)
!     real*8 yp1,ypn,x(n),y(n),y2(n)
!     parameter (nmax=500)

      integer i,k
      real*8 p,qn,sig,un,u(nmax)

      if(yp1>0.99d30)then
        y2(1)=0d0
        u(1)=0d0
      else
        y2(1)=-0.5d0
        u(1)=(3d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      end if
      do i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2d0
         y2(i)=(sig-1d0)/p
         u(i)=(6d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/  &
              (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      end do
      if(ypn>0.99d30)then
        qn=0d0
        un=0d0
      else
        qn=0.5d0
        un=(3d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      end if
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1d0)
      do k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+u(k)
      end do

      return
      end subroutine spline

      subroutine itsolv(itmax,delmin,l,d0,ti,tf,f,d,info)
      implicit none

!     solves the equation 6*D=[exp(-6*D*ti)-exp(-6*D*tf)]/F(ti,tf) iteratively, 
!     by putting 6*D(i+1)=[exp(-6*D(i)*ti)-exp(-6*D(i)*tf)]/F(ti,tf)
!     where F(ti,tf) is input (integral[dt*C(t)] from ti->tf) 


      integer itmax,info,l
      real*8 delmin,d0,ti,tf,f,d
      real*8 del
      integer i
      real*8 fac
      
      fac=dble(l*(l+1))
      i=1
      del=1d10
      do while ((i<=itmax).and.(del>delmin))
         d=(dexp(-fac*d0*ti)-dexp(-fac*d0*tf))
         d=d/(fac*f)
         del=dabs((d-d0)/d0)
!        if(info==1)then
!          write(3,'(i6,2x,3(e15.8,2x))') i,d0,d,del
!        end if
         d0=d
         i=i+1
      end do
      if((i>itmax).and.(del>delmin))then
         write(6,*) 'warning, itsolv did not converge'
         write(6,*) '# iterations=',i,'fractional change=',del
      else
!        write(6,*) 'converged: # iterations=',i
      end if

      return
      end subroutine itsolv
