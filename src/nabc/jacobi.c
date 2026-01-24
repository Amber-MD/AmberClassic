#include <stdio.h>
#include <math.h>
#include "defreal.h"
#include "memutil.h"

int jacobi(REAL_T ** a, int n, REAL_T * d, REAL_T ** v, int *nrot)
{
   int i, j, iq, ip;
   REAL_T tresh, theta, tau, t, sm, s, h, g, c, *b, *z;

   b = vector(1, n);
   z = vector(1, n);

   for (ip = 1; ip <= n; ip++) {
      for (iq = 1; iq <= n; iq++)
         v[ip][iq] = 0.0;
      v[ip][ip] = 1.0;
   }

   for (ip = 1; ip <= n; ip++) {
      b[ip] = d[ip] = a[ip][ip];
      z[ip] = 0.0;
   }

   *nrot = 0;

   for (i = 1; i <= 50; i++) {

      sm = 0.0;
      for (ip = 1; ip <= n - 1; ip++) {
         for (iq = ip + 1; iq <= n; iq++)
            sm += fabs(a[ip][iq]);
      }
      if (sm == 0.0) {
         free_vector(z, 1, n);
         free_vector(b, 1, n);
         return (0);            /* all is well */
      }

      if (i < 4)
         tresh = 0.2 * sm / (n * n);
      else
         tresh = 0.0;

      for (ip = 1; ip <= n - 1; ip++) {
         for (iq = ip + 1; iq <= n; iq++) {
            g = 100.0 * fabs(a[ip][iq]);

            if (i > 4 &&
                fabs(d[ip]) + g == fabs(d[ip]) &&
                fabs(d[iq]) + g == fabs(d[iq]))
               a[ip][iq] = 0.0;
            else if (fabs(a[ip][iq]) > tresh) {
               h = d[iq] - d[ip];
               if (fabs(h) + g == fabs(h))
                  t = (a[ip][iq]) / h;
               else {
                  theta = 0.5 * h / (a[ip][iq]);
                  t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
                  if (theta < 0.0)
                     t = -t;
               }
               c = 1.0 / sqrt(1 + t * t);
               s = t * c;
               tau = s / (1 + c);
               h = t * a[ip][iq];
               z[ip] -= h;
               z[iq] += h;
               d[ip] -= h;
               d[iq] += h;
               a[ip][iq] = 0.0;
               for (j = 1; j <= ip - 1; j++) {
                  g = a[j][ip];
                  h = a[j][iq];
                  a[j][ip] = g - s * (h + g * tau);
                  a[j][iq] = h + s * (g - h * tau);
               }
               for (j = ip + 1; j <= iq - 1; j++) {
                  g = a[ip][j];
                  h = a[j][iq];
                  a[ip][j] = g - s * (h + g * tau);
                  a[j][iq] = h + s * (g - h * tau);
               }
               for (j = iq + 1; j <= n; j++) {
                  g = a[ip][j];
                  h = a[iq][j];
                  a[ip][j] = g - s * (h + g * tau);
                  a[iq][j] = h + s * (g - h * tau);
               }
               for (j = 1; j <= n; j++) {
                  g = v[j][ip];
                  h = v[j][iq];
                  v[j][ip] = g - s * (h + g * tau);
                  v[j][iq] = h + s * (g - h * tau);
               }
               ++(*nrot);
            }
         }
      }

      for (ip = 1; ip <= n; ip++) {
         b[ip] += z[ip];
         d[ip] = b[ip];
         z[ip] = 0.0;
      }

   }

   return (1);
}

void eigsrt(REAL_T * d, REAL_T ** v, int n)
{
   int k, j, i;
   REAL_T p;

   for (i = 1; i < n; i++) {
      p = d[k = i];
      for (j = i + 1; j <= n; j++) {
         if (d[j] >= p)
            p = d[k = j];
      }
      if (k != i) {
         d[k] = d[i];
         d[i] = p;
         for (j = 1; j <= n; j++) {
            p = v[j][i];
            v[j][i] = v[j][k];
            v[j][k] = p;
         }
      }
   }
}
