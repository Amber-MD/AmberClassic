subroutine addsugarpsatoms(maxatom,maxres,natom,nres, &
                           igraph,ipres,lbres,ib,jb,nbond,c, &
                           ftype,fhybrid)
   
   ! Two pseudoatoms are inserted into the sugar rings of nucleic acids
   !   resulting into a 7 membered ring with 1 dof.
   !   The pseudoatoms are inserted between C1'-O4' and C4'-O4'. 
   !   The former corresponding bonds are removed, respectively. 
   ! Due to the increased ring size hc between C1' and C5' of the same 
   !   sugar ring are explicitely forbidden in forFIRSTteth.f
   ! Note that only standard nucleosides are considered.
   !
   ! Simone Fulle - 26.02.2008
   
   implicit none
   
   double precision, intent(inout) :: c(*)
   double precision :: x1, y1, z1, x2, y2, z2, xd, yd, zd, &
                       vecfac
   
   character(len=4), intent(inout) :: igraph(*), lbres(*)
   
   integer, intent(in) :: maxatom, maxres
   integer, intent(inout) :: natom, nres, ipres(*), fhybrid(*), nbond, &
                             ib(*), jb(*)
   
   ! --- fixed number of pseudo atoms applied in sugar rings *2---
   integer, parameter :: pseudonum = 1
   
   integer :: i, j, k, m, n, nat, nb, nrs, resid, atom
   
   character(len=1), intent(inout) :: ftype(*)
   
   ! --- Factor to calc coords of tether atoms ---
   vecfac = 1.0 / dble(pseudonum + 1)
   
   nrs = nres
   nat = natom
   nb = nbond
   
   ! --- Loop over all residues ---
   do resid=1,nrs
      ! -- in the case of a nucleic acid
      if (lbres(resid)(1:1) .eq. 'R' .or. &
          lbres(resid)(1:1) .eq. 'D') then
         ! --- Loop over all atoms i in residue r ---
         do i=ipres(resid), (ipres(resid+1)-1)
            if (igraph(i).eq.'O4''') then
               atom = i
               x1 = c(3*i-2)
               y1 = c(3*i-1)
               z1 = c(3*i)
            end if
         end do
         
         ! -- neighbor of C4' has index "atom-2" and of C1' has index "atom+1" --
         do n=1,2
            if (n .eq. 1) then
               j = atom-2
            else
               j = atom+1
            end if
            
            x2 = c(3*j-2)
            y2 = c(3*j-1)
            z2 = c(3*j)
            
            xd = x2 - x1
            yd = y2 - y1
            zd = z2 - z1
            
            nres = nres + 1
            lbres(nres) = 'BMH'
            ipres(nres) = nat + 1
            
            ! --- Update atom and bond information ---                  
            xd = xd * vecfac
            yd = yd * vecfac
            zd = zd * vecfac
            
            nat = nat + 1
            if (nat .gt. maxatom) then 
               write(6,*) 'Too many atoms: ', nat
               stop
            end if
            k = pseudonum
            
            c(3*nat-2) = x1 + k * xd
            c(3*nat-1) = y1 + k * yd
            c(3*nat  ) = z1 + k * zd
            
            igraph(nat) = 'X'
            ftype(nat) = 'N'
            fhybrid(nat) = 0
            nbond = nbond + 1
            if (nbond .gt. maxatom) then 
               write(6,*) 'Too many bonds: ', nbond
               stop
            end if
            
            ! -- insert bond between pseudoatom = nat and O4' = atom --         
            ib(nbond) = (atom-1)*3 
            jb(nbond) = (nat-1)*3
            
            ! -- instead of bond between C4' and O4 --
            ! -- insert either bond between pseudoatom = nat and C4' = atom-2 --
            if(n.eq.1) then 
               do m=1,nb
                  if (ib(m) .eq. ((atom-3)*3) .and. &
                      jb(m) .eq. ((atom-1)*3)) then
                     !ib(m) = (atom-3)*3   
                     jb(m) = (nat-1)*3
                  end if
               end do
            ! -- ... or bond between pseudoatom = nat and C1' = atom+1 --
            else
               do m=1,nb
                  if(ib(m) .eq. ((atom-1)*3) .and. &
                     jb(m) .eq. ((atom)*3)) then
                     ib(m) = (nat-1)*3   
                  end if
               end do           
            end if
            
            natom = nat
            if ((nres+1) .gt. maxres) then
               write(6,*) 'Too many residues: ', nres+1
               stop
            end if
            ipres(nres+1) = natom+1
         end do
      end if
   end do
   
   return
   
end subroutine addsugarpsatoms

