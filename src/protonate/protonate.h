c
c
c not used
c
ccc      integer prot_type, mol, srf
ccc      parameter  (prot_type = 1, mol = 1, srf = 2)
c
c NLINSZ = Maximum line length for filenames etc
c
      integer NLINSZ
      parameter (NLINSZ=1024)
c
c NR = maximum number of residues in the stdprt.dat file.
c NHSET = maximum number of sets of hydrogens for any residue:
c MAXAT = maximum number of atoms:
c
      integer NR, NHSET, MAXAT, NSSBOND
      parameter (NR=100)
      parameter (NHSET=25)
      parameter (NSSBOND=200)
      parameter (MAXAT=30000)
c
c protein common info
c
      integer nprot(NR), iprot(NR,NHSET),
     .        indat(NR,NHSET,4), nres(NR)
      common /proti/ nprot, iprot, indat, nres
c
      character*4 cdhname(NR,NHSET,3)
      character*4 cdatname(NR,50), cdatnamh(NR,NHSET)
      common /protc/ cdhname, cdatname, cdatnamh
c
c residue common info
c
      integer natom, ires(MAXAT)
      real    coord(3,MAXAT)
      logical isoldp(MAXAT), matched(MAXAT)
      common /crdi/ natom, ires, coord, isoldp, matched
c
      character*4 cnamat(MAXAT), catom(6)
      character*3 cresnam(MAXAT)
      character*1 chainid(MAXAT), crescod(MAXAT)
      character*6 ckey(MAXAT)
      common /crdc/ cnamat, catom, cresnam, chainid, ckey,
     .              crescod
c
c files common info
c
      integer nfilin,  nfilout, ndatin, nerrout,
     .        nlnkout, nedtout, nprmout, ndbgout
      common /filei/ nfilin,  nfilout, ndatin, nerrout,ndbgout
c
      character*(NLINSZ) datafile, infile, outfile, logfile, dbgfile
      common /filec/ datafile, infile, outfile, logfile, dbgfile
c
c temporary storage for writeout
c
      real    scoord
      integer icount, ilast, iresp, iresct
      logical atout, ipopt, keep, mismch, lrenum, lprint
c
      common /wouti/ scoord(3,6), icount, ilast, iresp, iresct,
     +               atout, ipopt, keep, mismch, lrenum, lprint
c
      character*1 chainp, cinsert
c
      common /woutc/ chainp, cinsert
c
c save the common blocks as statics
c
      save /proti/, /protc/, /crdi/, /crdc/, /wouti/, /woutc/,
     +     /filei/, /filec/
c
