#include "bc_memory.h"

integer       natom,nres,nbonh,nbona,ntheth,ntheta,nphih, &
      nphia,nnb,ntypes,nconp,maxmem,nwdvar,nparm, &
      natc,ibelly,natbel,ishake,nmxrs, &
      mxsub,natyp,npdec,i02,i04,i06,i08,i10, &
      iibh,ijbh,iicbh,iiba,ijba,iicba, &
      i24,i26,i28,i30,i32,i34,i36,i38,i40, &
      i42,i44,i46,i48,i50,i52,i54,i56,i58,ibellygp, &
      icnstrgp,i64,i65,i66,i68, &
      i70,i72,i74,i76,i78,i79,i80,i82,i84,i86, &
      icpstinf,icpresst,icptrsct, icpptcnt, &
      l05,l10,l15,lwinv,lpol,lcrd,lforce,l36,lvel,lvel2,l45,l50, &
      lcrdr,l60,l65,lmass,l75,l80,l85,l90,l95,l96,l97,l98,lfrctmp, &
      l105,l110,l115,l120,l125,l130,l135,l140,l145,l150, &
      l165,l170,l175,l180,l185,l186,l187,l188,l189,l190, &
      m02,m04,m06,m08,m10,m12,m14,m16,m18,i01, &
      lastr,lasti,lasth,lastpr,ncopy, &
      nbper,ngper,ndper,ifpert
!     lastpr,lastrst,lastist,
! ,lpolp, ncopy, &
!     imask1,imask2,numadjst,mxadjmsk

! 1     2      3      4      5       6       7      8      9      10
common/memory/ &
 natom ,nres  ,nbonh ,nbona ,ntheth ,ntheta ,nphih ,                       & ! 7
 nphia ,nnb   ,ntypes,nconp ,maxmem ,nwdvar ,nparm ,                      &  !14
 natc  ,ibelly,natbel,ishake,nmxrs  ,                                      & !19
 mxsub ,natyp ,npdec ,i02   ,i04    ,i06    ,i08   ,i10,                   & !27
 iibh  ,ijbh  ,iicbh ,iiba  ,ijba   ,iicba  ,                              & !33
 i24   ,i26   ,i28   ,i30   ,i32    ,i34    ,i36   ,i38   ,i40   ,         & !42
 i42   ,i44   ,i46   ,i48   ,i50    ,i52    ,i54   ,i56   ,i58   ,ibellygp,& !52
 icnstrgp,i64 ,i65   ,i66   ,i68    ,                                      & !57
 i70   ,i72   ,i74   ,i76   ,i78    ,i79    ,i80   ,i82   ,                & !65
 i84   ,i86   ,                                                            & !67
 l05   ,l10   ,l15   ,lwinv ,lpol   ,lcrd   ,lforce,l36   ,lvel  ,lvel2 ,  & !77
 l45   ,l50   ,                                                            & !79
 lcrdr ,l60   ,l65   ,lmass ,l75    ,l80    ,l85   ,l90   ,l95   ,l96   ,  & !89
 l97   ,l98   ,lfrctmp,                                                    & !92
 l105  ,l110  ,l115  ,l120  ,l125   ,l130   ,l135  ,l140  ,l145  ,l150  ,  & !102
 l165  ,l170  ,l175  ,l180  ,l185   ,l186   ,l187  ,l188  ,l189  ,l190  ,  & !112
 m02   ,m04   ,m06   ,m08   ,m10    ,m12    ,m14   ,m16   ,m18   ,i01   ,  & !122
 lastr ,lasti ,lasth ,lastpr,ncopy,                                       & !127
 nbper,ngper,ndper,ifpert                                                  !131 
!iifstwt,iifstwr,nrealb,nintb,nholb,npairb,lastr ,lasti ,lasth ,         & !129
!lastpr ,lpolp ,ncopy,          & !136
!imask1 ,imask2 ,numadjst,mxadjmsk                                         !140
