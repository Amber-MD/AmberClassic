/*
 *     variables for Isotropic Periodic Sum (IPS) calculation
 *
 *     ips         IPS options: 1--for both ele and l-j
 *                              2--for ele only
 *                              3--for l-j only
 *     nnbips      Number of nonbonded atom pairs
 *     nnbipst     Provious Number of nonbonded atom pairs
 */
      static int ips,nnbipst,nnbips;
/*
 *   IPS parameters
 *   bipse*    Electrostatic
 *   bipsva*   Lennard-Jones repulsion
 *   bipsvc*   Lennard-Jones dispersion
 *   rips*     Radius of IPS local region 
 *   pips*0    Self IPS pair energies 
 *   pips*c    IPS system energy components 
 *   eipss*c   IPS system energies 
 *   virips*c  IPS system virials 
 *   vboxips   IPS local region volume
 */

#define NINT(x) ( (int) ((x) > 0 ? ((x) + 0.5) : ((x) - 0.5)) )

#define PI (3.1415926535897932384626433832795)
#define FOURPI (4.*3.1415926535897932384626433832795)
#define INVSQRT2 (1.0/1.4142135623730950488016887242097)

      static REAL_T  aipse,bipse0,bipse1,bipse2,bipse3,
       aipsva,bipsva0,bipsva1,bipsva2,bipsva3,
       aipsvc,bipsvc0,bipsvc1,bipsvc2,bipsvc3,
       rips,rips2,onerips6,onerips12,
       pipse0,pipsva0,pipsvc0,pipsec,pipsvac,pipsvcc,
       vboxips,eipssnb,eipssel,virips,ripsinv,rips2inv;

void ipssys( )
{

      int i,nsumit[500];
      int iti,itj,ic,niti,nitj;
      REAL_T  aij,cij,anij;
      REAL_T  aijsum,cijsum,cgsum,cgijsum, rips6;

/* IPS Radius: */

      rips=cut;
      rips2=cut*cut;
      rips6=rips2*rips2*rips2;
      onerips6=1./rips6;
      onerips12=onerips6*onerips6;
      ripsinv=1./sqrt(rips2);
      rips2inv = ripsinv*ripsinv;

/* Ele IPS parameters:  */

      bipse1=1.109466;
      bipse2=0.0708;
      bipse3=0.0079175;
      bipse0=1.0-3.0*bipse1-5.0*bipse2-7.0*bipse3;
      /* aipse=0.0*sqrt(2.0)-bipse0; */
      aipse=-bipse0;

/*  Dispersion IPS parameters:  */

      bipsvc1=-1.70554;
      bipsvc2=-0.469571;
      bipsvc3=0.0000;
      bipsvc0=1.0-(4.0*bipsvc1+5.0*bipsvc2+6.0*bipsvc3)/3.;
      aipsvc=8.0*0.4376630-bipsvc0;

/*  Repulsion IPS parameters: */

      bipsva1=2.664085;
      bipsva2=-0.611616;
      bipsva3=-0.776091;
      bipsva0=1.0-(7.0*bipsva1+8.0*bipsva2+9.0*bipsva3)/6.;
      aipsva=64.0*0.006353604-bipsva0;

/* Energy and force constants:  */

      pipsec=1.0+bipse0+bipse1+bipse2+bipse3;
      pipsvac=1.0+bipsva0+bipsva1+bipsva2+bipsva3;
      pipsvcc=1.0+bipsvc0+bipsvc1+bipsvc2+bipsvc3;
      pipsva0=bipsva0/64.0-pipsvac;
      pipsvc0=bipsvc0/8.-pipsvcc;
      pipse0=bipse0*INVSQRT2-pipsec;

      eipssnb=0.0;
      eipssel=0.0;

      for( i=0; i<prm->Ntypes; i++){
         nsumit[i] = 0;
      }
      cgsum=0.0;

      for( i=0; i<prm->Natom; i++){
        cgsum=cgsum+prm->Charges[i];
        iti=prm->Iac[i];
        if(iti > prm->Ntypes){
           fprintf( stderr, "problem with types: %5d %5d %5d\n",
              i, iti, prm->Ntypes );
           exit(1);
        }
        nsumit[iti]=nsumit[iti]+1;
      }
      cijsum=0.0;
      aijsum=0.0;
      nnbipst=0;

      /*  system energy is calculated based on all pairs: */

      for( iti=0; iti<prm->Ntypes; iti++ ){
        ic=iti*(iti-1)/2+iti;
        aij=prm->Cn1[ic];
        cij=prm->Cn2[ic];
        niti=nsumit[iti];
        anij=niti*niti*0.5;
        cijsum=cijsum+cij*anij;
        aijsum=aijsum+aij*anij;
        nnbipst=nnbipst+niti*niti;
        for( itj=iti+1; itj<prm->Ntypes; itj++){
          nitj=nsumit[itj];
          ic=itj*(itj-1)/2+iti;
          aij=prm->Cn1[ic];
          cij=prm->Cn2[ic];
          anij=niti*nitj;
          cijsum=cijsum+cij*anij;
          aijsum=aijsum+aij*anij;
          nnbipst=nnbipst+2*niti*nitj;
        }
      }
      cgijsum=cgsum*cgsum*0.5;
      if( cgijsum <  1.e-15 ) cgijsum = 0.;
      eipssnb=aijsum*(pipsvac+aipsva/64.0)*onerips12
                -cijsum*(pipsvcc+aipsvc*0.125)*onerips6;
      eipssel=cgijsum*(pipsec+aipse*INVSQRT2)*ripsinv;
      if(!teips)eipssel=0.0;
      if(!tvips)eipssnb=0.0;

      /* Calculate volume virials: */

      virips=-(eipssnb+eipssel);
      vboxips=FOURPI*rips2*rips/3.;
      if(get_mytaskid() == 0){
        fprintf( nabout, " IPS Radius: %6.2f\n", rips );
        if(teips){
          fprintf( nabout, " IPS parameters for electrostatic energy:\n");
          fprintf( nabout, 
                   "  aipse,bipse0,bipse1,bipse2,bipse3:\n%10.6f%10.6f%10.6f%10.6f%10.6f\n",
                   aipse,bipse0,bipse1,bipse2,bipse3 );
        }
        if(tvips){
          fprintf( nabout, " IPS parameters for L-J energy:\n" );
          fprintf( nabout, "  aipsvc,bipsvc0,bipsvc1,bipsvc2,bipsvc3:\n%10.6f%10.6f%10.6f%10.6f%10.6f\n",
                   aipsvc,bipsvc0,bipsvc1,bipsvc2,bipsvc3 );
          fprintf( nabout, "  aipsva,bipsva0,bipsva1,bipsva2,bipsva3:\n%10.6f%10.6f%10.6f%10.6f%10.6f\n",
                   aipsvc,bipsvc0,bipsvc1,bipsvc2,bipsvc3 );
        }
        fprintf (nabout, " eipssnb, eipssel=%20.7f%20.7f\n", eipssnb,eipssel );
        if(prm->IfBox)
          fprintf( nabout, "  IPS region volume=%20.7f\n\n", vboxips );
        else
          fprintf( nabout, "  Enclosing atom pairs/total pairs=%10d%10d\n\n",
                   nnbips,nnbipst );
      }
}

void ipsupdate( REAL_T volume)
{

/*
 *-----------------------------------------------------------------------
 *     Update parameters once IPS radius or the box size changed
 *-----------------------------------------------------------------------
 */
      REAL_T fips,change ;

      if(prm->IfBox){
        fips=vboxips/volume;
        /*fips=nnbips*1.0/nnbipst;*/
      }
      change=fips-1.0;
      change=change*change;

      /* Update system energies and forces: */
      if(change > 1.0e-8){
         virips=virips*fips;
         eipssnb=eipssnb*fips;
         eipssel=eipssel*fips;
      }
      vboxips=volume;
      nnbipst=nnbips;
}

void eexips(REAL_T *enb, REAL_T *eel, REAL_T *dx, REAL_T *x)
{
/*
 *----------------------------------------------------------------------
 *  3D IPS interaction between excluded atom pairs
 *  This routine must be called first to update IPS parameters when needed
 *      and to initalize electrostatic and vdw energies
 *
 *  by Xiongwu Wu  - 9/10/2004
 *----------------------------------------------------------------------
 */

      /* REAL_T vir[3][3];  */
      int i,j,k,ic,iti,itj;
      int i3,j3;
      REAL_T dxi, dyi, dzi,dij,dijx,dijy,dijz;
      REAL_T enbij,eelij;
      REAL_T  xi,yi,zi,xij,yij,zij;
      REAL_T  cgi,cgij,aij,cij;
      REAL_T  r2,u2,twou1,twou2,onetwou2,twou6,twou12;
      REAL_T  pe,peu,deu,pva,pvau,dvau,pvc,pvcu,dvcu;

/*
      ! check to see if volume or atom pair changed 
      CALL IPSUPDATE(NTB)
*/

      /* Setup constants for use in inner loops */
      *enb=eipssnb;
      *eel=eipssel;
      deu=0.0;
      dvau=0.0;
      dvcu=0.0;
      cgij=0.0;
      aij=0.0;
      cij=0.0;
      /* vir(1:3,1:3) = 0.;  */

      nnbips=0;
      for( i=0; i<prm->Natom; i++){
         nnbips++;
         if(teips){
            cgi=prm->Charges[i];
            cgij=cgi*cgi;
            eelij=0.5*cgij*pipse0*ripsinv;
            *eel+=eelij;
         }
         if(tvips){
            iti=prm->Iac[i];
            ic=iti*(iti-1)/2+iti;
            aij=prm->Cn1[ic];
            cij=prm->Cn2[ic];
            /* Atom i long-range reference and self-interaction */
            enbij=0.5*(aij*pipsva0*onerips6-cij*pipsvc0)*onerips6;
            *enb+=enbij;
         }
         i3=3*i;
         dxi=dx[i3];
         dyi=dx[i3+1];
         dzi=dx[i3+2];
         xi=x[i3];
         yi=x[i3+1];
         zi=x[i3+2];
         for( k=0; k<prm->Iblo[i]; k++ ){
           j = IexclAt[i][k]-1;
           if(j<0) continue;  /*needed?*/
           j3=3*j;
           xij=x[j3]-xi;
           yij=x[j3+1]-yi;
           zij=x[j3+2]-zi;
           r2=xij*xij+yij*yij+zij*zij;
           nnbips=nnbips+2;
           u2=r2*rips2inv;
           twou2=2.0-u2;
           onetwou2=1.0/twou2;
           if( onetwou2 < 0.0 ){
             printf( "negsqrt: %5d%5d %12.5f %12.5f\n", i+1,j+1,onetwou2,r2 );
             printf( " %12.5f %12.5f %12.5f   %12.5f %12.5f %12.5f\n",
                    xi,yi,zi, x[j3],x[j3+1],x[j3+2] );
             exit(1);
           }
           if(teips){
              cgij=cgi*prm->Charges[j];
              twou1=sqrt(onetwou2);
              pe=bipse0+u2*(bipse1+u2*(bipse2+u2*bipse3));
              peu=2.0*(bipse1+u2*(2.0*bipse2+3.0*bipse3*u2));
              deu=(peu+pe*onetwou2)*twou1;
              eelij=cgij*(pe*twou1-pipsec)*ripsinv;
              *eel+=eelij;
           }
           if(tvips){
              itj=prm->Iac[j];
              if(iti > itj)
                 ic=iti*(iti-1)/2+itj;
              else
                 ic=itj*(itj-1)/2+iti;
              aij=prm->Cn1[ic];
              cij=prm->Cn2[ic];
              twou6=onetwou2*onetwou2*onetwou2;
              twou12=twou6*twou6;

              /*  L-J r6 term:  */

              pvc=bipsvc0+u2*(bipsvc1+u2*(bipsvc2+u2*bipsvc3));
              pvcu=2.0*(bipsvc1+u2*(2.0*bipsvc2+3.0*bipsvc3*u2));
              dvcu=(pvcu+6.0*pvc*onetwou2)*twou6;
 
              /*  L-J r12 term:  */

              pva=bipsva0+u2*(bipsva1+u2*(bipsva2+u2*bipsva3));
              pvau=2.0*(bipsva1+u2*(2.0*bipsva2+3.0*bipsva3*u2));
              dvau=(pvau+12.0*pva*onetwou2)*twou12;
              enbij=aij*(pva*twou12-pipsvac)*onerips12
                  -cij*(pvc*twou6-pipsvcc)*onerips6;

              *enb+=enbij;
              dij=-(cgij*deu*ripsinv +(aij*dvau*onerips6 
                   -cij*dvcu)*onerips6)*rips2inv;
           } else {
              dij=-cgij*deu*ripsinv*rips2inv;
           }

           dijx=dij*xij;
           dijy=dij*yij;
           dijz=dij*zij;

           dxi+=dijx;
           dyi+=dijy;
           dzi+=dijz;
           dx[j3]  -= dijx;
           dx[j3+1]-= dijy;
           dx[j3+2]-= dijz;
/*
           vir[0][0] -= dijx*xij;
           vir[0][1] -= dijx*yij;
           vir[0][2] -= dijx*zij;
           vir[1][0] -= dijy*xij;
           vir[1][1] -= dijy*yij;
           vir[1][2] -= dijy*zij;
           vir[2][0] -= dijz*xij;
           vir[2][1] -= dijz*yij;
           vir[2][2] -= dijz*zij;
*/
         }
         dx[i3]  =   dxi;
         dx[i3+1] =  dyi;
         dx[i3+2] =  dzi;
      }
}

int nblist_box(REAL_T * x, int *npairs, int **pairlist, REAL_T cut)
{
   REAL_T dx, dy, dz, rrw, cut2, xi, yi, zi, boxxi, boxyi, boxzi;
   int tot_pair, i, j, ipair, jexcl, nexcl, jexcl_last, i3, j3;
   int skip, maxnb;

   cut2 = cut * cut;
   tot_pair = 0;
   boxxi = 1./prm->Box[0];
   boxyi = 1./prm->Box[1];
   boxzi = 1./prm->Box[2];

   maxnb = cut * cut * cut/ 1.25;
   for (i = 0; i < prm->Natom; i++) {
      if( !pairlist[i] ) pairlist[i] = ivector(0, maxnb );
   }

/*
    loop over all pairs of atoms:
*/
#ifdef OPENMP
#pragma omp parallel for schedule (static, 1) \
  private (i, i3, dx, dy, dz, j, j3, xi, yi, zi, rrw, \
       jexcl, jexcl_last, nexcl, skip, ipair)
#endif

   for (i = 0; i < prm->Natom; i++) {

      i3 = 3 * i;
      xi = x[i3];
      yi = x[i3 + 1];
      zi = x[i3 + 2];

      /* Perform some excluded list expansion. */

      jexcl = 0;
      jexcl_last = prm->Iblo[i] - 1;

      if (jexcl > jexcl_last) {
         nexcl = -1;
      } else {
         nexcl = IexclAt[i][jexcl] - 1;
      }
      ipair = 0;

      for (j = i + 1; j < prm->Natom; j++) {

         /* Perform some more excluded list expansion. */

         skip = 0;
         if (j == nexcl) {
            skip = 1;
            jexcl++;
            if (jexcl > jexcl_last) {
               nexcl = -1;
            } else {
               nexcl = IexclAt[i][jexcl] - 1;
            }
         }
         if (skip)
            continue;

         j3 = 3 * j;

#if 0
         dx = remainder(xi - x[j3 + 0], prm->Box[0]);
         dy = remainder(yi - x[j3 + 1], prm->Box[1]);
         dz = remainder(zi - x[j3 + 2], prm->Box[2]);
#elif 0
/* round requires _ISOC99_SOURCE or _GNU_SOURCE
http://osdir.com/ml/lib.glibc.bugs/2002-09/msg00059.html
*/
         dx = xi - x[j3 + 0];
         dy = yi - x[j3 + 1];
         dz = zi - x[j3 + 2];
   fprintf(nabout, "        %g  %g  %g \n", dx,dy,dz);
         dx -= round( dx*boxxi )*prm->Box[0];
         dy -= round( dy*boxyi )*prm->Box[1];
         dz -= round( dz*boxzi )*prm->Box[2];
   fprintf(nabout, "        %g  %g  %g \n", dx,dy,dz);
#else
         dx = xi - x[j3 + 0];
         dy = yi - x[j3 + 1];
         dz = zi - x[j3 + 2];
         dx -= NINT( dx*boxxi )*prm->Box[0];
         dy -= NINT( dy*boxyi )*prm->Box[1];
         dz -= NINT( dz*boxzi )*prm->Box[2];
#endif
         rrw = dx * dx + dy * dy + dz * dz;
         if (rrw < cut2) {
            pairlist[i][ipair++] = j;
         }
      }
      if( ipair > maxnb ){
         fprintf( stderr, "Too many pairs, max is %d\n", maxnb );
         exit(1);
      }
      npairs[i] = ipair;

   }
   
   for (i = 0; i < prm->Natom; i++){
      tot_pair += npairs[i];
   }

   if(get_mytaskid()==0){
     fprintf(nabout, "  nblist_box, number of atoms= %u\n", prm->Natom);
     fprintf(nabout, "  nblist_box, number of pairs= %u\n", tot_pair);
   }
   return (tot_pair);
}

int nbond_box(int *npairs, int **pairlist,
          REAL_T * x, REAL_T * f, REAL_T * enb, REAL_T * eel)
{
   int i, j, jn, ic, npr, iaci;
   REAL_T dumx, dumy, dumz, cgi, xw1, xw2, xw3, r2inv, df2, r6, f1,
       f2, r2;
   REAL_T df, fw1, fw2, fw3;
   REAL_T rinv;
   int threadnum, numthreads, foff;
   REAL_T elec, evdw;
   REAL_T boxxi, boxyi, boxzi;  /* inverse of box lengths  */
   REAL_T uips, uips2, twou2, twou, pipse, dpipse, dipse, eipse;

   evdw = 0.;
   elec = 0.;
   boxxi = 1./prm->Box[0];
   boxyi = 1./prm->Box[1];
   boxzi = 1./prm->Box[2];

#ifdef OPENMP
#pragma omp parallel reduction (+: elec, evdw) \
  private (i, npr, iaci, dumx, dumy, dumz, cgi, jn, j, xw1, xw2, xw3, \
           r2, r2inv, rinv, df2, ic, r6, f2, f1, df, fw1, fw2, fw3, \
           uips, uips2, twou2, twou, pipse, dpipse, dipse, eipse, \
           threadnum, numthreads, foff)
#endif
   {

#ifdef OPENMP
     threadnum = omp_get_thread_num();
     numthreads = omp_get_num_threads();
     foff = 3 * prm->Natom * threadnum;
#else
     threadnum = 0;
     numthreads = 1;
     foff = 0;
#endif

   for (i = threadnum; i < prm->Natom - 1; i += numthreads) {

      npr = npairs[i];
      if (npr > 0) {
         iaci = prm->Ntypes * (prm->Iac[i] - 1);
         dumx = 0.;
         dumy = 0.;
         dumz = 0.;
         cgi = prm->Charges[i];

         for (jn = 0; jn < npr; jn++) {
            j = pairlist[i][jn];
            xw1 = x[3 * i + 0] - x[3 * j + 0];
            xw2 = x[3 * i + 1] - x[3 * j + 1];
            xw3 = x[3 * i + 2] - x[3 * j + 2];
            xw1 -= NINT( xw1*boxxi )*prm->Box[0];
            xw2 -= NINT( xw2*boxyi )*prm->Box[1];
            xw3 -= NINT( xw3*boxzi )*prm->Box[2];
            r2 = xw1 * xw1 + xw2 * xw2 + xw3 * xw3;
            r2inv = 1. / r2;
            rinv = sqrt(r2inv);

         /*
          * -- use a ipse long range potential:
          *         eipse=e0+(b0+b1r^2+b2r^4+b3r^6)/sqrt(2-r^2)
          * compare: ele=1/r  fele=-1/r^2
          */

         df2 = -cgi* prm->Charges[j]*rinv;

      if( teips ){
         uips=ripsinv*r2*rinv;
         uips2=r2*rips2inv;
         twou2=1.0/(2.0-uips2);
         twou=sqrt(twou2);
         pipse = bipse0 + uips2*(bipse1 + uips2*(bipse2 + uips2*bipse3));
         dpipse = 2.0*bipse1 + uips2*(4.0*bipse2 + 6.0*uips2*bipse3);
         dipse = uips*(dpipse + pipse*twou2)*twou*rips2inv;
         eipse = -df2*uips*(pipse*twou - pipsec);
         elec += -df2 + eipse;
         df = df2 * (r2inv - dipse);
      } else {
         elec += -df2;
         df = df2 * r2inv;
      }

            ic = prm->Cno[iaci + prm->Iac[j] - 1];
            if (ic > 0) {
               ic--;
               r6 = r2inv * r2inv * r2inv;
               f2 = prm->Cn2[ic] * r6;
               f1 = prm->Cn1[ic] * r6 * r6;

/*
**						if( f1 > 500000. )
**							fprintf( stderr, "close contact: %d %d  %8.3f\n", 
**							i+1,j+1,1./rinv );
*/

               evdw += (f1 - f2);
               df += (6. * (f2 - f1 - f1)) * r2inv;
            }
            fw1 = xw1 * df;
            fw2 = xw2 * df;
            fw3 = xw3 * df;
            dumx = dumx + fw1;
            dumy = dumy + fw2;
            dumz = dumz + fw3;
            f[foff + 3*j + 0] -= fw1;
            f[foff + 3*j + 1] -= fw2;
            f[foff + 3*j + 2] -= fw3;
         }
         f[foff + 3*i + 0] += dumx;
         f[foff + 3*i + 1] += dumy;
         f[foff + 3*i + 2] += dumz;
      }
   }

   }
   *eel = elec;
   *enb = evdw;
   return (0);
}
