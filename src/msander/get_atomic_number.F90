! <compile=optimized>
#include "../include/dprec.fh"
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! ----------------------
  ! Identify atom elements
  ! ----------------------
  subroutine get_atomic_number(atom_name,atom_mass,atomic_number, errorFlag)
    !
    ! Assign atomic number based upon the first letter of the atom symbol
    ! which has been read in from the topology file and based upon the mass.
    ! The assumption here is that an atom name matches the first letter of the
    ! element. If it doesn't then this routine will need to be modified.

    use UtilitiesModule, only : Upcase

    implicit none

    character(len=4), intent(in)  :: atom_name
    _REAL_,           intent(in)  :: atom_mass
    logical, optional, intent(out) :: errorFlag
    integer,          intent(out) :: atomic_number

    logical::localErrorFlag
    localErrorFlag=.false.

    ! Lanthanides are not supported.
    ! Actinides are not supported.

    if( Upcase(atom_name(1:1)) .eq. 'A' ) then
       if(atom_mass > 24.0d0 .and. atom_mass <= 28.0d0) then
          atomic_number =  13 !Aluminium
       elseif(atom_mass > 35.0d0 .and. atom_mass <= 40.0d0) then
          atomic_number =  18 !Argon
       elseif(atom_mass > 73.0d0 .and. atom_mass <= 77.0d0) then
          atomic_number =  33 !Arsenic
       elseif(atom_mass > 106.0d0 .and. atom_mass <= 109.0d0) then
          atomic_number =  47 !Silver
       elseif(atom_mass > 195.0d0 .and. atom_mass <= 199.0d0) then
          atomic_number =  79 !Gold
       elseif(atom_mass > 208.0d0 .and. atom_mass <= 212.0d0) then
          atomic_number =  85 !Astatine
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'B' ) then
       if(atom_mass > 8.0d0 .and. atom_mass <= 10.0d0) then
          atomic_number =  4 !Beryllium
       elseif(atom_mass > 10.0d0 .and. atom_mass <= 12.0d0) then
          atomic_number =  5 !Boron
       elseif(atom_mass > 77.0d0 .and. atom_mass <= 81.0d0) then
          atomic_number =  35 !Bromine
       elseif(atom_mass > 135.0d0 .and. atom_mass <= 139.0d0) then
          atomic_number =  56 !Barium
       elseif(atom_mass > 207.0d0 .and. atom_mass <= 211.0d0) then
          atomic_number =  83 !Bismuth
       else
          write(6,*) 'Unable to correctly identify element ', atom_name
          call mexit(6,1)
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'C' ) then
       if(atom_mass > 10.0d0 .and. atom_mass <= 14.0d0) then
          atomic_number =  6 !Carbon
       elseif(atom_mass > 33.0d0 .and. atom_mass <= 37.0d0) then
          atomic_number =  17 !Chlorine
       elseif(atom_mass > 38.0d0 .and. atom_mass <= 42.0d0) then
          atomic_number =  20 !Calcium
       elseif(atom_mass > 50.0d0 .and. atom_mass <= 54.0d0) then
          atomic_number =  24 !Chromium
       elseif(atom_mass > 57.0d0 .and. atom_mass <= 61.0d0) then
          atomic_number =  27 !Cobalt
       elseif(atom_mass > 61.0d0 .and. atom_mass <= 65.0d0) then
          atomic_number =  29 !Copper
       elseif(atom_mass > 110.0d0 .and. atom_mass <= 114.0d0) then
          atomic_number =  48 !Cadmium
       elseif(atom_mass > 131.0d0 .and. atom_mass <= 135.0d0) then
          atomic_number =  55 !Cesium
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'D' ) then
       localErrorFlag=.true.

    elseif( Upcase(atom_name(1:1)) .eq. 'E' ) then
       localErrorFlag=.true.

    elseif( Upcase(atom_name(1:1)) .eq. 'F' ) then
       if(atom_mass > 17.0d0 .and. atom_mass <= 21.0d0) then
          atomic_number =  9 !Fluorine
       elseif(atom_mass > 54.0d0 .and. atom_mass <= 58.0d0) then
          atomic_number =  26 !Iron
       elseif(atom_mass > 218.0d0 .and. atom_mass <= 228.0d0) then
          atomic_number =  87 !Francium
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'G' ) then
       if(atom_mass > 67.0d0 .and. atom_mass <= 71.0d0) then
          atomic_number =  31 !Gallium
       elseif(atom_mass > 71.0d0 .and. atom_mass <= 75.0d0) then
          atomic_number =  32 !Germanium
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'H' ) then
       if(atom_mass > 0.0d0 .and. atom_mass <= 2.0d0) then
          atomic_number =  1 !Hydrogen
       elseif(atom_mass > 3.0d0 .and. atom_mass <= 5.0d0) then
          atomic_number =  2 !Helium
       elseif(atom_mass > 176.0d0 .and. atom_mass <= 180.0d0) then
          atomic_number =  72 !Hafnium
       elseif(atom_mass > 198.0d0 .and. atom_mass <= 202.0d0) then
          atomic_number =  80 !Mercury
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'I' ) then
       if(atom_mass > 112.0d0 .and. atom_mass <= 116.0d0) then
          atomic_number = 49 !Indium
       elseif(atom_mass > 125.0d0 .and. atom_mass <= 129.0d0) then
          atomic_number =  53 !Iodine
       elseif(atom_mass > 190.0d0 .and. atom_mass <= 194.0d0) then
          atomic_number =  77 !Iridium
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'J' ) then
       localErrorFlag=.true.

    elseif( Upcase(atom_name(1:1)) .eq. 'K' ) then
       if(atom_mass > 37.0d0 .and. atom_mass <= 41.0d0) then
          atomic_number = 19 !Potassium
       elseif(atom_mass > 77.0d0 .and. atom_mass <= 86.0d0) then
          atomic_number = 36 !Krypton
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'L') then
       if(atom_mass > 6.0d0 .and. atom_mass <= 8.0d0) then
          atomic_number = 3 !Lithium
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'M' ) then
       if(atom_mass > 22.0d0 .and. atom_mass <= 26.0d0) then
          atomic_number = 12 !Magnesium
       elseif(atom_mass > 53.0d0 .and. atom_mass <= 57.0d0) then
          atomic_number = 25 !Manganese
       elseif(atom_mass > 94.0d0 .and. atom_mass <= 98.0d0) then
          atomic_number = 42 !Molybdenem
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'N') then
       if(atom_mass > 13.0d0 .and. atom_mass <= 15.0d0) then
          atomic_number = 7 !Nitrogen
       elseif(atom_mass > 19.0d0 .and. atom_mass <= 22.0d0) then
          atomic_number = 10 !Neon
       elseif(atom_mass > 22.1d0 .and. atom_mass <= 23.0d0) then
          atomic_number = 11 !Sodium
       elseif(atom_mass > 57.0d0 .and. atom_mass <= 61.0d0) then
          atomic_number = 28 !Nickel
       elseif(atom_mass > 95.0d0 .and. atom_mass <= 99.0d0) then
          atomic_number = 41 !Niobium
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'O' ) then
       if(atom_mass > 14.0d0 .and. atom_mass <= 18.0d0) then
          atomic_number = 8 !Oxygen
       elseif(atom_mass > 188.0d0 .and. atom_mass <= 192.0d0) then
          atomic_number = 76 !Osmium
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'P' ) then
       if(atom_mass > 29.0d0 .and. atom_mass <= 33.0d0) then
          atomic_number = 15 !Phosphorus
       elseif(atom_mass > 104.0d0 .and. atom_mass <= 108.0d0) then
          atomic_number = 46 !Palladium
       elseif(atom_mass > 193.0d0 .and. atom_mass <= 197.0d0) then
          atomic_number = 78 !Platinum
       elseif(atom_mass > 205.0d0 .and. atom_mass <= 208.0d0) then
          atomic_number = 82 !Lead
       elseif(atom_mass > 208.0d0 .and. atom_mass <= 212.0d0) then
          atomic_number = 84 !Polonium
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'Q' ) then
       localErrorFlag=.true.

    elseif( Upcase(atom_name(1:1)) .eq. 'R' ) then
       if(atom_mass > 84.0d0 .and. atom_mass <= 88.0d0) then
          atomic_number = 37 !Rubidium
       elseif(atom_mass > 99.0d0 .and. atom_mass <= 102.0d0) then
          atomic_number = 44 !Ruthenium
       elseif(atom_mass > 102.0d0 .and. atom_mass <= 105.0d0) then
          atomic_number = 45 !Rhodium
       elseif(atom_mass > 184.0d0 .and. atom_mass <= 188.0d0) then
          atomic_number = 75 !Rhenium
       elseif(atom_mass > 210.0d0 .and. atom_mass <= 222.5d0) then
          atomic_number = 86 !Radon
       elseif(atom_mass > 223.0d0 .and. atom_mass <= 229.0d0) then
          atomic_number = 88 !Radium
       else
          localErrorFlag=.true.
       end if

    elseif( Upcase(atom_name(1:1)) .eq. 'S' ) then
       if(atom_mass > 26.0d0 .and. atom_mass <= 30.0d0) then
          atomic_number = 14 !Silicon
       elseif(atom_mass > 30.0d0 .and. atom_mass <= 34.0d0) then
          atomic_number = 16 !Sulphur
       elseif(atom_mass > 43.0d0 .and. atom_mass <= 47.0d0) then
          atomic_number = 21 !Scandium
       elseif(atom_mass > 77.0d0 .and. atom_mass <= 81.0d0) then
          atomic_number = 34 !Selenium
       elseif(atom_mass > 86.0d0 .and. atom_mass <= 89.0d0) then
          atomic_number = 38 !Strontium
       elseif(atom_mass > 116.0d0 .and. atom_mass <= 120.0d0) then
          atomic_number = 50 !Tin
       elseif(atom_mass > 120.0d0 .and. atom_mass <= 124.0d0) then
          atomic_number = 51 !Antimony
       else
          localErrorFlag=.true.
       end if

    elseif( Upcase(atom_name(1:1)) .eq. 'T' ) then
       if(atom_mass > 46.0d0 .and. atom_mass <= 50.0d0) then
          atomic_number = 22 !Titanium
       elseif(atom_mass > 96.0d0 .and. atom_mass <= 100.0d0) then
          atomic_number = 43 !Technetium
       elseif(atom_mass > 125.0d0 .and. atom_mass <= 130.0d0) then
          atomic_number = 52 !Tellurium
       elseif(atom_mass > 179.0d0 .and. atom_mass <= 183.0d0) then
          atomic_number = 73 !Tantalum
       elseif(atom_mass > 201.0d0 .and. atom_mass <= 206.0d0) then
          atomic_number = 81 !Thallium
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'U' ) then
       write(6,*) 'Unable to correctly identify element ', atom_name
       call mexit(6,1)

    elseif( Upcase(atom_name(1:1)) .eq. 'V' ) then
       if(atom_mass > 49.0d0 .and. atom_mass <= 53.0d0) then
          atomic_number = 23 !Vanadium
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'W' ) then
       if(atom_mass > 179.0d0 .and. atom_mass <= 183.0d0) then
          atomic_number = 74 !Tungsten
       else
          localErrorFlag=.true.
       endif

    elseif( Upcase(atom_name(1:1)) .eq. 'X' ) then
       if (atom_mass > 123.0d0 .and. atom_mass < 136.0d0) then
       else
          localErrorFlag=.true.
       end if

    elseif( Upcase(atom_name(1:1)) .eq. 'Z' ) then
       if(atom_mass > 61.0d0 .and. atom_mass <= 69.0d0) then
          atomic_number = 30 !Zinc
       elseif(atom_mass > 89.0d0 .and. atom_mass <= 93.0d0) then
          atomic_number = 40 !Zirconium
       else
          localErrorFlag=.true.
       end if

    elseif( Upcase(atom_name(1:1)) .eq. 'Z' ) then
       if(atom_mass > 89.0d0 .and. atom_mass <= 93.0d0) then
          atomic_number = 40 !Zirconium
       else
         localErrorFlag=.true.
       end if

    else
       localErrorFlag=.true.
    endif

    if (localErrorFlag) then
       if (present(errorFlag)) then
           ! not to force the program to terminate--left to the calling program
           ! what to do next---Taisung Lee (Rutgers, 2011)
           errorFlag=.true.
       else
           write(6,'(2a)') 'Unable to correctly identify element ', atom_name
           write(6,'(a)') 'Note: element guessing does not work with Hydrogen'
           write(6,'(a)') '      Mass Repartitioning if ATOMIC_NUMBER is not'
           write(6,'(a)') '      present in the topology file'
           call mexit(6,1)
       end if
    end if

  end subroutine get_atomic_number
