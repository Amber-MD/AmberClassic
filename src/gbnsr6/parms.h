
!    BC_PARMR is the number of reals in common RPARMS; BC_PARMI is the
!    number of ints in IPARMS.  (Change these if you change the sizes below).

#define BC_PARMR 24520
#define BC_PARMI 1200

_REAL_ rk(5000),req(5000),tk(900),teq(900),pk(1200),      &!13000
      pn(1200),phase(1200),cn1(1830),cn2(1830),solty(60), &!19120
      gamc(1200),gams(1200),fmn(1200),                    &!22720
      asol(200),bsol(200),hbcut(200),                     &!24520
      one_scee(1200) 
common/rparms/rk,req,tk,teq,pk, &
      pn,phase,cn1,cn2,solty, &
      gamc,gams,fmn, &
      asol,bsol,hbcut, &
      one_scee

integer ipn(1200)
common/iparms/ipn

#ifdef CHARMM
_REAL_ cn114,cn214,rkub,rub
common/p14/ cn114(1830),cn214(1830)
common/ub/rkub(900),rub(900)
#endif

! NPHB is the number of h-bond parameters. NIMPRP is the number of
! improper torsional parameters (NPTRA-NIMPRP is the number of regular
! torsional parameters).

integer       numbnd,numang,nptra,nphb,nimprp
common/prmlim/numbnd,numang,nptra,nphb,nimprp
