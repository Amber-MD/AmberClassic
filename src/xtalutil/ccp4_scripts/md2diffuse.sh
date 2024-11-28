#!/bin/bash

#  Overview of workflow to use MD snapshots to compute diffuse and Bragg
#    scattering.  Output is an mtz file (or formatted analog) containing
#    Bragg and diffuse scattering intensities.  An average electron density
#    map can also be created.

#  Basically, all of the input pdb-format files are fed to sfall to get
#    computed structure factors.  The remaining steps just process the
#    output mtz files, with an eye to minimizing disk storage.

#    (I generally remove files in the $pdbprefix and save_hkl folders
#    once I am done, since the input pdb files, and output mtz files,
#    can take up a *lot* of space for large problems.)

#  This script requires ccp4 programs to be in the PATH: sfall, f2mtz,
#    fft, cad, mtzdump.  It also requires a set of PDB-format files
#    from an MD (or other) ensemble;  Amber's cpptraj program is used
#    in the example to generate these, but other tools could replace that.

#  The expectation is that this file will be copied to your
#    working directory, and edited and run there.  Each major section is
#    enclosed in an "if" block, since one will often want to run things
#    step by step.  Keeping the edited version of this script around provdes
#    a record of what was done.  (In this master version, all the if blocks
#    start with "if false; then"; you need to change some or all of these
#    to "if true; then" to actually do anything.

#  Comments or problems: send email to dacase1@gmail.com, or (better)
#    create an issue at github.com/dacase/xtaltools.

#=============================================================================
#  Input variables: edit these to match your problem:

pdbprefix="../PDBdata/lys_ortho"     # pdbfiles will be called
                                      # $pdbprefix.$frame.pdb
dprefix="lys_ortho"                  # basename for final output files
numframes=2500                       # number of pdb-format snapshots

cell="CRYST1   30.490   56.385   73.824  90.00  90.00  90.00 P 21 21 21    4\n"
                                # should take this from the first pdb file
title="Diffuse/Bragg for lys_ortho"  # for mtz and map files
vf000=""                              # cell volume and number of electrons
grid=""                               # grid dimensions for final map
                                      # (see step 7 for vf000 and grid)

resolution=1.12                   # high resolution limit for output
let resolutionm=$resolution-0.1   # give a slightly better resolution to
                                  # sfall to esure that no hkl points get
                                  # omitted

#  (You will also need to edit Step 1 below if you plan to use cpptraj
#     from AmberTools to prepare the input PDB files.)

#=============================================================================
#  Step 1.  Run cpptraj to prepare PDB files

if false; then

cpptraj <<EOF
#  sample cpptraj script to create pdb snapshots for diffuse scattering analysis
# 
parm ../lys_ortho_1uc.parm7
reference ../md-1.x
trajin ../md9.nc
trajin ../md10.nc
trajin ../md11.nc
trajin ../md12.nc
trajin ../md13.nc
trajin ../md14.nc
trajin ../md15.nc
trajin ../md16.nc
trajin ../md17.nc
rms reference bb '@C,CA,N' norotate out fit.dat time 0.4
image byatom
trajout $pdbprefix.pdb pdb multi pdbv3 keepext sg "P 1"
go
EOF

fi

# Note that this might be a good place to run the bendfinder.com scripts on
# each of the above snapshots (if desired.)

#=============================================================================
# Step  2. Get a reference mtz file so each frame can have the same hkl values:

if false; then

sfall xyzin $pdbprefix.1.pdb hklout tmp.mtz memsize 10000000 <<EOF > sfall.log
mode sfcalc xyzin
symm P1
reso $resolution
NOSCALE
end
EOF

#  (note:  not sure if this step is really needed for current workflow....)
cad hklin1 tmp.mtz hklout frame1.mtz <<EOF
labin file 1 E1=FC E2=PHIC
labou file 1 E1=FP E2=SIGFP
ctypin file 1 E2=Q
EOF

#  (note: using mtzdump to take advantage of the LRESO feature....)
cat <<EOF > format_mtzdump.pl
#
#  take output from mtzdump and format to an rdb tab-delimited file
#
print "h\tk\tl\tq\n4N\t4N\t4N\t10N\n";

while (<STDIN>){
  \$lineno++;  next if \$lineno < 83;
  next if \$_ =~ /^ MTZDUMP/;
  next if \$_ =~ /^Times:/;
  next if \$_ =~ /^</;

  ( \$h, \$k, \$l, \$dinv2 ) = unpack("x a4 a4 a4 a9", \$_ );
  printf STDOUT "%d\t%d\t%d\t%10.5f\n", \$h, \$k, \$l, 6.2831853*sqrt(\$dinv2);
}
EOF

mtzdump hklin frame1.mtz <<EOF | perl ./format_mtzdump.pl > frame1.hkl
NREF -1
LRESO
GO
EOF

/bin/rm -f format_mtzdump.pl tmp.mtz

fi

#=============================================================================
#  Step 3.  On each frame, run sfall to get structure factors:

if false; then

mkdir -p save_hkl

#  script to add constant bfact to input snapshot PDB file:
cat <<EOF > modify_pdb
#!/usr/bin/perl -n
#
#   first, read in a pdb file line and unpack
#

if( \$_ =~ /^ATOM|^HETATM/){

#  expand to 80+ characters (really hokey version):
chop(\$_);
\$_ .= "                                                                    ";

( \$label, \$atno, \$atname, \$alt, \$resname, \$chainId, \$resno, \$iCode,
      \$x, \$y, \$z, \$occ, \$bfact, \$element, \$charge ) =
    unpack("a6 a5 x a4 a a3 x a a4 a x3 a8 a8 a8 a6 a6 x10 a2 a2",\$_);

#
#  do modifications necessary here:
#
\$bfact = 15.;

#
#  write back out
#
printf 
"%6s%5s %-4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
      \$label, \$atno,\$atname,\$alt, \$resname,\$chainId, \$resno, \$iCode,
        \$x,\$y,\$z, \$occ, \$bfact, \$element, \$charge;

} elsif ( \$_ =~ /^CRYST1/) {   # make sure we have a common record
    print "$cell";
} else {
   print;
}
EOF
chmod +x modify_pdb

#  script to save just fc,phic to intermediate files:
cat <<EOF > mtz2fcphic.c
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

//  dump rows of data from mtz file argv[1]
//  fc, phic written to stdout

int main( int argc, char** argv){

   int nread, nwrite, ncol, nrec, nbatch, head_loc;
   float *all, *fc, *phic;
   char head_string[4000];

   FILE* fp = fopen( argv[1], "r" );
   assert( fp );

   fseek( fp, 4, SEEK_SET );
   fread( &head_loc, sizeof(head_loc), 1, fp );
   fseek( fp, 4*head_loc, SEEK_SET );
   nread = fread( head_string, 1, 3999, fp );  /* slurp in much of the header */
   head_string[nread] = '\0';
   // fprintf( stderr, "|%s|\n", head_string );
   
   char* buf = strstr( head_string, "NCOL" );  assert( buf );
   sscanf( buf+4, "%d %d %d", &ncol, &nrec, &nbatch );
   // fprintf( stderr, "ncol, nrec, nbatch: %d,%d,%d\n", ncol, nrec, nbatch );

   //  allocate space for all the rows:
   all = malloc( sizeof(float) * ncol * nrec );
   assert( all );
   fc = malloc( sizeof(float) * nrec );
   assert( fc );
   phic = malloc( sizeof(float) * nrec );
   assert( phic );

   fseek( fp, 80, SEEK_SET );

   nread = fread( all, sizeof(float), nrec*ncol, fp );
   assert( nread == nrec*ncol );
   fclose( fp );

   for( int i=0; i<nrec; i++ ){
      fc[i] = all[ ncol*i + 5 ];
      phic[i] = all[ ncol*i + 6 ];
   }

   nwrite = fwrite( fc, sizeof(float), nrec, stdout );
   assert( nwrite == nrec );
   nwrite = fwrite( phic, sizeof(float), nrec, stdout );
   assert( nwrite == nrec );

}
EOF
gcc -std=gnu99  -o mtz2fcphic mtz2fcphic.c

#  Loop over input files:

for frame in {1..$numframes}; do

#  set all b-factors to 15:
./modify_pdb < ${pdbprefix}.$frame.pdb > $frame.pdb

#  run SFALL to get structure factors:
sfall xyzin $frame.pdb hklin frame1.mtz hklout $frame.pdb.mtz \
     memsize 10000000 <<EOF > $frame.sfall2.log
mode sfcalc xyzin hklin
symm P1
RESOLUTION $resolutionm
NOSCALE
end
EOF

#  for smaller, binary (fcphic) output:
./mtz2fcphic $frame.pdb.mtz > save_hkl/$frame.fcphic

/bin/rm $frame.pdb $frame.pdb.mtz $frame.sfall2.log

echo Done with frame $frame at `date` >> progress

done

/bin/rm modify_pdb mtz2fcphic.c mtz2fcphic

fi

#=============================================================================
#  Step 4.  Run the "diffuse1" program to compute block averages:

if false; then

cat <<EOF > diffuse1.c
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

//  Take input *.fcphic files in the save_hkl directory
//  input arguments: nfile1, nfile2, nskip, <ref.mtz>
//  intermediate file, with <freal> <fimag> <fsq>, to stdout

int main( int argc, char** argv )
{

   int head_loc;   // location in MTZ file of header
   char head_string[4000];
   int ncol, nrec, nbatch, nread, nwrite;
   int *h, *k, *l;
   char filename[250];
   double phi;

   // parse the input arguments:
   int nfile1 = atoi( argv[1] );
   int nfile2 = atoi( argv[2] );
   int nskip =  atoi( argv[3] );

   //  read the ref.mtz file to get hkl + number of relections:
   FILE* ref = fopen( argv[4], "r" ); assert( ref );
   fseek( ref, 4, SEEK_SET );
   fread( &head_loc, sizeof(head_loc), 1, ref );

   fseek( ref, 4*head_loc, SEEK_SET );
   nread = fread( head_string, 1, 3999, ref );
   head_string[nread] = '\0';
   // fprintf( stderr, "|%s|\n", head_string );
   char* buf = strstr( head_string, "NCOL" );  assert( buf );
   sscanf( buf+4, "%d %d %d", &ncol, &nrec, &nbatch );
   // fprintf( stderr, "ncol, nrec, nbatch: %d,%d,%d\n", ncol, nrec, nbatch );

   h = malloc( sizeof(int) * nrec );
   k = malloc( sizeof(int) * nrec );
   l = malloc( sizeof(int) * nrec );

   float row[ncol];  // holds data for each row
   fseek( ref, 80, SEEK_SET );

   for( int i=0; i<nrec; i++ ){
      fread( row, sizeof(row), 1, ref );
      h[i] = (int)row[0];
      k[i] = (int)row[1];
      l[i] = (int)row[2];
      // fprintf( stderr, "%d %d %d\n", h[i],k[i],l[i] );
   }
   fclose( ref );

   //  allocate arrays
   double* frsum = calloc( nrec, sizeof(double) ); assert( frsum );
   double* fisum = calloc( nrec, sizeof(double) ); assert( fisum );
   double* fsqsum = calloc( nrec, sizeof(double) ); assert( fsqsum );
   float* fc = malloc( sizeof(float) * nrec ); assert( fc );
   float* phic = malloc( sizeof(float) * nrec ); assert( phic );

   // Loop over the input files from sfall:

   int nfiles = 0;
   double DEG_TO_RAD=0.01745329252;
   FILE *fp;
   for( int ifile=nfile1; ifile<=nfile2; ifile+=nskip ){

      snprintf( filename, 250, "save_hkl/%d.fcphic", ifile );
      fprintf( stderr, "%s\n", filename );
      fp = fopen( filename, "r" ); assert( fp );
      nfiles += 1;
      nread = fread( fc, sizeof(float), nrec, fp );
      assert( nread == nrec );
      nread = fread( phic, sizeof(float), nrec, fp );
      assert( nread == nrec );
      fclose( fp );

      for( int i=0; i<nrec; i++ ){
         phi = DEG_TO_RAD*phic[i];
         frsum[i] += fc[i]*cos(phi);
         fisum[i] += fc[i]*sin(phi);
         fsqsum[i] += fc[i]*fc[i];
      }
   }
   fprintf( stderr, "Processed %d files\n", nfiles );

#if 0
   //  dump averages to formatted file:

   for( int i=0; i<nrec; i++ ){
      printf( "%5d\t%5d\t%5d\t%12.4f\t%12.4f\t%15.4f\n", 
         h[i],k[i],l[i],frsum[i]/nfiles, fisum[i]/nfiles, fsqsum[i]/nfiles );
   }
#else
   //  dump averages to unformatted file, to minimize loss of precision:

   for( int i=0; i<nrec; i++ ){
      frsum[i] /= nfiles;
      fisum[i] /= nfiles;
      fsqsum[i] /= nfiles;
   }

   nwrite = fwrite( frsum, sizeof(double), nrec, stdout );
   assert( nwrite == nrec );
   nwrite = fwrite( fisum, sizeof(double), nrec, stdout );
   assert( nwrite == nrec );
   nwrite = fwrite( fsqsum, sizeof(double), nrec, stdout );
   assert( nwrite == nrec );

#endif

}
EOF

gcc -std=gnu99  -O3 -o diffuse1 diffuse1.c -lm

# (This example just computes a single block average, over all frames.
# But the diffuse1 program could be run over subsets of the data, in order
#   look at block averages.)

./diffuse1 1 $numframes 1 frame1.mtz  > $dprefix.1.ihklb

/bin/rm -f diffuse1 diffuse1.c

fi

#=============================================================================
#  Step 5.  combine (if needed) several intermediate .ihkl files into a total:

if false; then

cat <<EOF > diffuse2.c
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

//  Take input intermediate files from diffuse1, sum all together,
//     correct for bfactor broadening
//  arguments:  ref.hkl, *.ihkl
//  formatted dhkl rdb file to stdout


int main( int argc, char** argv )
{

   char filename[250];
   char line[80];
   int hin, kin, lin, nread;
   int nfiles = argc - 2;

   //  read the ref.hkl file to get hkl + q:
   FILE* ref = fopen( argv[1], "r" ); assert( ref );

   //  first scan the file to get the number of records:
   int nrec = 0;
   while ( fgets( line, 80, ref )){ nrec++; }
   nrec -= 2;  // correct for header records 
   fprintf( stderr, "%s: nrec = %d, nfiles = %d\n", argv[1], nrec, nfiles );
   rewind( ref );

   int* h = malloc( sizeof(int) * nrec ); assert( h );
   int* k = malloc( sizeof(int) * nrec ); assert( k );
   int* l = malloc( sizeof(int) * nrec ); assert( l );
   double* q = malloc( sizeof(double) * nrec ); assert( q );

   fgets( line,  80, ref );  // skipping headers
   fgets( line,  80, ref );
   for( int i=0; i<nrec; i++ ){
      fgets( line,  80, ref );
      sscanf( line, "%d %d %d %lf", &h[i], &k[i], &l[i], &q[i] );
   }
   fclose( ref );

   //  allocate arrays
   double* frsum = calloc( nrec, sizeof(double) ); assert( frsum );
   double* fisum = calloc( nrec, sizeof(double) ); assert( fisum );
   double* fsqsum = calloc( nrec, sizeof(double) ); assert( fsqsum );
   double* fr = malloc( sizeof(double) * nrec ); assert( fr );
   double* fi = malloc( sizeof(double) * nrec ); assert( fi );
   double* fsq = malloc( sizeof(double) * nrec ); assert( fsq );

   // Loop over the intermediate *.ihkl files:

   FILE *fp;
   for( int ifile=2; ifile<argc; ifile++ ){

      fprintf( stderr, "%s\n", argv[ifile] );
      fp = fopen( argv[ifile], "r" ); assert( fp );
      nread = fread( fr, sizeof(double), nrec, fp );
      assert( nread == nrec );
      nread = fread( fi, sizeof(double), nrec, fp );
      assert( nread == nrec );
      nread = fread( fsq, sizeof(double), nrec, fp );
      assert( nread == nrec );

      for( int i=0; i<nrec; i++ ){
         frsum[i] += fr[i];
         fisum[i] += fi[i];
         fsqsum[i] += fsq[i];
      }
      fclose( fp );
   }

   //  dump averages to formatted file:

#if 1
   printf( "h\tk\tl\tq\tIdiff\tIBragg\tfav\tphiav\n" );
   printf( "4N\t4N\t4N\t10N\t15N\t15N\t12N\t12N\n" );
#else
   printf( "h\tk\tl\tdh\tdk\tdl\tq\tIdiff\tIBragg\tho\tko\tlo\n" );
   printf( "4N\t4N\t4N\t4N\t4N\t4N\t10N\t12N\t12N\t4N\t4N\t4N\n");
#endif

   double bfact = 15.;  // should this be an input parameter?
   double bs2over2, fav, phiav, fav2, fsqav, idiff;
   double RAD_TO_DEG = 57.29577951;
   int dh, dk, dl, hb, kb, lb, check;

   for( int i=0; i<nrec; i++ ){
      //  un-do B-factor. broadening:
      bs2over2 = exp(0.5*bfact*q[i]*q[i]/39.478417604357);
      fav2 = (frsum[i]*frsum[i] + fisum[i]*fisum[i])*bs2over2/(nfiles*nfiles);
      fav = sqrt( fav2 );
      phiav = RAD_TO_DEG*atan2(fisum[i],frsum[i]);
      fsqav = fsqsum[i]*bs2over2/nfiles;
      idiff = fsqav - fav2;

#if 1
      printf( "%d\t%d\t%d\t%10.5f\t%12.4f\t%12.4f\t%12.4f\t%10.4f\n", 
         h[i],k[i],l[i],q[i], idiff, fav2, fav, phiav );

#else  /* here for oversampled data: hardwire grids for now... */

       /* if one has multiple unit cells, and wishes to index the outputs
          according to the primary cell, turn this section on and enter
          the muliples in the h,k,l directions here:  */
       int hgrid = 7;
       int kgrid = 7;
       int lgrid = 7;

      //  find the central Bragg indices and offsets:
      hb = h[i]/hgrid;  dh = h[i]%hgrid; 
         if( dh > hgrid/2 ){ hb++; dh -= hgrid;}
         if( dh < -hgrid/2 ){ hb--; dh += hgrid; }
         // assert( hgrid*hb + dh == h[i] );

      kb = k[i]/kgrid;  dk = k[i]%kgrid; 
         if( dk > kgrid/2 ){ kb++; dk -= kgrid;}
         if( dk < -kgrid/2 ){ kb--; dk += kgrid; }
         // assert( kgrid*kb + dk == k[i] );

      lb = l[i]/lgrid;  dl = l[i]%lgrid; 
         if( dl > lgrid/2 ){ lb++; dl -= lgrid;}
         if( dl < -lgrid/2 ){ lb--; dl += lgrid; }
         // assert( lgrid*lb + dl == l[i] );

      printf( "%d\t%d\t%d\t%d\t%d\t%d\t%10.4f\t%12.4f\t%12.4f\t%d\t%d\t%d\n", 
         hb, kb, lb, dh, dk, dl, q[i], idiff, fav2, h[i], k[i], l[i] );
#endif

   }

}
EOF
gcc -std=gnu99  -O -o diffuse2 diffuse2.c -lm

#  (Again, the following line assumes a single "block" that includes all
#  frames.  Modify the following line if you created several ihklb files
#  in the previous step.)

./diffuse2 frame1.hkl $dprefix.1.ihklb > $dprefix.1.dhkl

/bin/rm -f diffuse2 diffuse2.c

fi
#=============================================================================
#  Step 6.  convert the .dhkl file above to mtz:

if false; then

f2mtz hklin $dprefix.dhkl hklout $dprefix.mtz <<EOF > makemap.log
CELL  $cell
SYMMETRY P1
LABOUT H K L Q IDIFFUSE FC PHIC
CTYPOUT H H H R J F P
SKIP 2
TITLE $title
END
EOF

fi

#=============================================================================
#  Step 7.  (Optional): Also get an average map file.  (You need to
#     know how to run the ccp4 "fft" program in order to set the VF000
#     and GRID variables.)

if false; then 

/bin/rm -f $dprefix.map
fft hklin $dprefix.mtz mapout $dprefix.map <<EOF >> makemap.log
LABIN F1=FC PHI=PHIC
VF000 $vf000
GRID $grid
EOF

fi


