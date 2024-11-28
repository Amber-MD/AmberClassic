// printneck.c : this will print the GB Neck lookup tables to make sure it's the
//               same as the Fortran code...

typedef double REAL_T; // needed by gbneck.h

#include "gbneck.h"
#include <stdio.h>

void print_array(const double[21][21], const char*);

int main() {

   int i;

   // Printing neckMaxPos
   print_array(neckMaxPos, "neckMaxPos");
   fprintf(stdout, "\n");
   print_array(neckMaxVal, "neckMaxVal");

   return 0;
}

void print_array(const double array[21][21], const char* name) {

   int i, j;

   for (i = 0; i <= 20; i++) {
      for (j = 0; j <= 20; j+=3)
         fprintf(stdout, "%s[%2i][%2i] = %10.8f; %s[%2i][%2i] = %10.8f; %s[%2i][%2i] = %10.8f\n",
                 name, i, j, array[i][j], name, i, j+1, array[i][j+1],
                 name, i, j+2, array[i][j+2]);
   }
}
