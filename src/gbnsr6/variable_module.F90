! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
#define REQUIRE(e) if(.not.(e)) call croak(__FILE__,__LINE__)
#include "pb_def.h"
#include "timer.h"

module variable_module

   implicit none

#  include "pb_constants.h"

   ! PBMD parameters

   _REAL_, parameter :: pbkb   = 1.3807D-23 / 1.6606D-27 / (1.00D+12)**2 * (1.00D+10)**2
   _REAL_, parameter :: fioni  = 6.0220D+23 / 1.00D+30
   _REAL_, parameter :: fiono  = ONE / fioni
   _REAL_, parameter :: eps0   = 8.8542D-12 / (1.6022D-19)**2 / (1.00D+10)**3 * (1.00D+12)**2 * 1.6606D-27
   _REAL_, parameter :: frcfac = 0.01D0 / 4.1840D0

   ! PBMD FD control variables

   logical :: outphi
   logical :: srsas
   logical :: scalerf
   logical :: outlvlset    !WMBS output total lvlset 
   logical :: outmlvlset   !WMBS output membrane lvlset

  
   integer :: phiform 
   integer :: saopt
   integer :: sasopt
   integer :: dbfopt
   integer :: eneopt
   integer :: npbopt
   integer :: solvopt
   integer :: frcopt
   integer :: intopt
   integer :: bcopt
   integer :: smoothopt
   integer :: xm
   integer :: ym
   integer :: zm
   integer :: xmym
   integer :: xmymzm
   integer :: nbuffer
   integer :: level
   integer :: nfocus
   integer :: fscale
   integer :: maxitn
   integer :: itn
   integer :: m, n
   integer :: savbcopt(MAXLEVEL)
   integer :: levelblock(MAXLEVEL)
   integer :: savxm(MAXLEVEL)
   integer :: savym(MAXLEVEL)
   integer :: savzm(MAXLEVEL)
   integer :: savxo(MAXLEVEL) !Origin offset
   integer :: savyo(MAXLEVEL) !
   integer :: savzo(MAXLEVEL)
   integer :: savxmym(MAXLEVEL)
   integer :: savxmymzm(MAXLEVEL)
   integer :: isurfchg           ! print surface charges
   integer :: membraneopt        !WMBS - membrane type 0 for none 1 for slab
   !integer :: mdielectricopt    !WMBS - place holder for future use
   integer :: poretype   !WMBS- membrane exclusion region flag 0 for none 1 for cylinder
   integer :: augtoltype         !WMBS - Aug IIM gmres Tolerance type
                                 !       0 for absolute, 1 for relative

   _REAL_ :: h
   _REAL_ :: gox
   _REAL_ :: goy
   _REAL_ :: goz
   _REAL_ :: fmiccg
   _REAL_ :: fmiccg2  
   _REAL_ :: accept
   _REAL_ :: laccept
   _REAL_ :: wsor
   _REAL_ :: lwsor
   _REAL_ :: norm
   _REAL_ :: inorm
   _REAL_ :: xmax
   _REAL_ :: xmin
   _REAL_ :: ymax
   _REAL_ :: ymin
   _REAL_ :: zmax
   _REAL_ :: zmin
   _REAL_ :: gxmax
   _REAL_ :: gxmin
   _REAL_ :: gymax
   _REAL_ :: gymin
   _REAL_ :: gzmax
   _REAL_ :: gzmin
   _REAL_ :: savxbox(MAXLEVEL)
   _REAL_ :: savybox(MAXLEVEL)
   _REAL_ :: savzbox(MAXLEVEL)
   _REAL_ :: cxbox(MAXLEVEL)
   _REAL_ :: cybox(MAXLEVEL)
   _REAL_ :: czbox(MAXLEVEL)
   _REAL_ :: savh(MAXLEVEL)
   _REAL_ :: savgox(MAXLEVEL)
   _REAL_ :: savgoy(MAXLEVEL)
   _REAL_ :: savgoz(MAXLEVEL)
   _REAL_ :: offx
   _REAL_ :: offy
   _REAL_ :: offz
   _REAL_ :: fillratio

   _REAL_ :: epsin
   _REAL_ :: epsout
   _REAL_ :: pbkappa
   _REAL_ :: istrng
   _REAL_ :: ivalence
   _REAL_ :: pbtemp
   _REAL_ :: totcrg
   _REAL_ :: totcrgp
   _REAL_ :: totcrgn

   _REAL_ :: pbgamma_int
   _REAL_ :: pbgamma_ext
   _REAL_ :: qef(3)
   _REAL_ :: ref(3)

   _REAL_ :: mthick     !WMBS- membrane thickness
   _REAL_ :: mctrdz     !WMBS- membrane z offset
   _REAL_ :: epsmemb    !WMBS- membrane dielectric constant
   _REAL_ :: poreradius !WMBS- radius for cylindrical region

   _REAL_ :: augctf     !WMBS- Aug IIM cluster radius (0 for h*h)
   _REAL_ :: augtol     !WMBS- Aug IIM gmres tolerance

   ! PBMD topology information

   integer              :: lastp
   integer              :: ngrdcrg

   integer, allocatable ::    icrd(:,:)
   integer, allocatable ::  grdcrg(:,:)
   _REAL_, allocatable :: qgrdcrg(:)
   _REAL_, allocatable ::    gcrd(:,:)
   _REAL_, allocatable ::    acrd(:,:)
   _REAL_, allocatable ::    acrg(:)
   _REAL_, allocatable ::    gcrg(:,:)
 
   ! PBMD nblist information

   integer              :: maxnbr
   integer              :: maxnba
   _REAL_              :: cutres, cutnb, cutfd, cutsa
 
   integer, allocatable ::   nshrt(:)
   integer, allocatable ::     nex(:)
   integer, allocatable ::     iex(:,:)
   integer, allocatable :: iprlong(:)
   integer, allocatable :: iprshrt(:)
   integer, allocatable ::  iar1pb(:,:)
   _REAL_, allocatable ::   cn1pb(:)
   _REAL_, allocatable ::   cn2pb(:)
   _REAL_, allocatable ::   cn3pb(:)

   ! PBMD cap water simulation information

   integer              :: mpopt
   integer              :: lmax
   integer              :: inatm
   integer              :: outwat
   integer              :: oution
   integer, allocatable :: outflag(:)
   integer, allocatable :: outflagorig(:)
   integer, allocatable :: mapout(:)
   integer, allocatable :: ibelly(:)
   _REAL_              :: sepbuf

   ! physical variables for energy and force calculations

   integer:: nbnd
   integer:: nbndx
   integer:: nbndy
   integer:: nbndz
   _REAL_, allocatable :: pos_crg(:,:,:)
   _REAL_, allocatable :: surf_crg(:,:)
   integer, allocatable :: ipos_crg(:,:,:)
   integer, allocatable :: crg_num(:)

   ! physical variable maps for numerical solutions

   _REAL_, allocatable ::     phi(:)
   _REAL_, allocatable ::      bv(:)
   _REAL_, allocatable ::   chgrd(:)
   _REAL_, allocatable ::    epsx(:)
   _REAL_, allocatable ::    epsy(:)
   _REAL_, allocatable ::    epsz(:)
   _REAL_, allocatable :: saltgrd(:)

   ! geometry maps for dielectric interface
 
   integer, allocatable ::  insas(:)
   integer, allocatable :: atmsas(:)
   _REAL_, allocatable ::  lvlset(:)
   _REAL_, allocatable ::      zv(:)

   ! physical variable maps for force calculations

   _REAL_, allocatable ::     cphi(:)
   integer, allocatable ::  iepsav(:,:)
   integer, allocatable :: iepsavx(:,:)
   integer, allocatable :: iepsavz(:,:)
   integer, allocatable :: iepsavy(:,:)
   _REAL_, allocatable ::   fedgex(:)
   _REAL_, allocatable ::   fedgey(:)
   _REAL_, allocatable ::   fedgez(:)

   ! saved phi array for pbmd
 
   _REAL_, allocatable :: xs(:)
   integer :: xsoffset

   ! ligand focusing options

   logical :: ligand
   character(len=256) ligandmask
   integer, allocatable :: liveflag(:)
   integer, allocatable :: realflag(:)
   integer :: ntrajmol
   _REAL_ :: buffer

   ! Multiple distributive fine grid geometry / Multiblock focusing

   logical :: multiblock                          !TRUE if specified multiblock
   logical :: firstleveldone                      !TRUE if level 1 is done once
   integer, allocatable :: blkxo(:)               !block origin in x dir
   integer, allocatable :: blkyo(:)               !block origin in y dir
   integer, allocatable :: blkzo(:)               !block origin in z dir
   integer, allocatable :: blkxlo(:)              !block lower bound in x
   integer, allocatable :: blkylo(:)              !block lower bound in y
   integer, allocatable :: blkzlo(:)              !block lower bound in z
   integer, allocatable :: blkxup(:)              !block upper bound in x
   integer, allocatable :: blkyup(:)              !block upper bound in y
   integer, allocatable :: blkzup(:)              !block upper bound in z
   integer, allocatable :: blknx(:)               !block grid length in x
   integer, allocatable :: blkny(:)               !block grid length in y
   integer, allocatable :: blknz(:)               !block grid length in z
   integer, allocatable :: blknxny(:)             !blknx . blkny
   integer, allocatable :: blknynz(:)             !blkny . blknz
   integer, allocatable :: blknxnz(:)             !blknx . blknz
   integer, allocatable :: blknxnynz(:)           !blknx . blkny . blknz
   integer              :: ngrdblkx               !# of grids per block in x
   integer              :: ngrdblky               !# of grids per block in y
   integer              :: ngrdblkz               !# of grids per block in z
   integer              :: xmblk                  !total blocks in x
   integer              :: ymblk                  !total blocks in y
   integer              :: zmblk                  !total blocks in z
   integer              :: ilower                 !
   integer              :: iupper                 !
   integer              :: jlower                 !
   integer              :: jupper                 !
   integer              :: klower                 !
   integer              :: kupper                 !
   integer              :: nfirstwpad(MAXBLOCK+1) !Described in the text
   integer              :: nfirst0pad(MAXBLOCK+1) !Described in the text
   integer              :: nfirst4s3f(MAXBLOCK+1) !Described in the text
   integer, allocatable :: lfocuswpad(:)          !Described in the text
   integer, allocatable :: lfocus0pad(:)          !Described in the text
   integer, allocatable :: lfocus4s3f(:)          !Described in the text
   integer, allocatable :: fineblockindex(:)      !for load balancing
   integer              :: iblkdebug              !
   _REAL_,  allocatable :: blkgox(:)              !block origin in x
   _REAL_,  allocatable :: blkgoy(:)              !block origin in y
   _REAL_,  allocatable :: blkgoz(:)              !block origin in z
   _REAL_,  allocatable :: coarsephi(:)           !saved 1st level solution
!  _REAL_, allocatable  :: spv(:)                 !Ion placement
   integer              :: saltout                !Ion placement
   logical              :: outsalt                !Ion placement
   _REAL_               :: stern                  !Ion Placement

end module variable_module
