!     common block sizes:

integer bc_boxi,bc_boxr
parameter(bc_boxi=7)
parameter(bc_boxr=611)

! ... floats:

#  define _REAL_ double precision
_REAL_ box,cut,dielc,rad,wel,radhb,welhb, &     
      cutcap,xcap,ycap,zcap,fcap,rwell
common/boxr/box(3),cut,dielc, &               !5
      cutcap,xcap,ycap,zcap,fcap,rwell, &     !11
      rad(100),wel(100),radhb(200),welhb(200) !611

! ... integers:

integer ntb,ifbox,numpk,nbit,ifcap,natcap,isftrp

common/boxi/ntb,ifbox,numpk,nbit,ifcap,natcap,isftrp
