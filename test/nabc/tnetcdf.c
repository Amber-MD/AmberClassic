//  test of netcdf capabilities.

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "AmberNetcdf.h"

int main( int argc, char *argv[] )
{
   struct AmberNetcdf ain, aout;
   double x[700];

  if( netcdfLoad( &ain, "trpcage.nc" )){
      fprintf( stderr, "Error: unable to load trpcage.nc trajectory\n" );
      exit(1);
  }
  netcdfDebug(&ain);
  netcdfInfo(&ain);

  netcdfGetFrame(&ain,0,x,NULL);
  printf("First 3 coords of frame 0: %lf %lf %lf\n",x[1],x[2],x[3]);
  netcdfCreate(&aout, "trpcageTest.nc", ain.ncatom, 0);
  netcdfDebug(&aout);
  netcdfInfo(&aout);
  while ( netcdfGetNextFrame(&ain,x,NULL) ) {
   printf("First 3 coords of frame %i: %lf %lf %lf\n",ain.currentFrame,
          x[1],x[2],x[3]);
    netcdfWriteNextFrame(&aout,x,NULL);
  }  
  netcdfClose(&ain);
  netcdfClose(&aout);

}
