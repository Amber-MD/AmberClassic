c=======================================================================
c
      subroutine sortf_as (num,x, isort)
c----------------------------------------------------------------------c
c     ascending order                                                  c
c----------------------------------------------------------------------c
      parameter (maxarr=100000)
      dimension  x(*), isort(*)
      dimension  xs(maxarr)
c
      if (num.ge.maxarr) then
        write(6,*) 
     +      'sortf_as.f: sort array too small: redimension/recompile'
        write(6,*) '  size needed is ', num
        stop
      endif
      isort(1) = 1
      if (num.eq.1) return
c
c
      do 10 i=2,num
         isort(i) = i
   10 continue
c
      xs(1) = x(1)
      do 100 i=2,num
         j = i - 1
   55    if (x(i).ge.xs(j)) then
            xs(j+1)    = x(i)
            isort(j+1) = i
            go to 100
         else
            xs(j+1)    = xs(j)
            isort(j+1) = isort(j)
            if (j.eq.1) then
               xs(1)    = x(i)
               isort(1) = i
               go to 100
            else
               j = j - 1
               go to 55
            endif
         endif
  100 continue
c
      return
      end
