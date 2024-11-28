!<compile=optimized>

module memory_module

implicit none

#include "memory.h"

! Remove lastrst and lastist from the common block, since they're not in the
! sander memory.h common block. Eventually I think everything should be moved
! into the module (as in sander), but this is OK as a first step. Now,
!
! "use memory_module"
!
! has the same effect as
!
! #include "../include/memory.h"
!
! and should be used instead.  It also has the upside of being able to use the
! keyword "only".

integer, save :: lastrst, lastist

end module memory_module
