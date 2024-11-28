/* ******************************************************* */
/* com.c: center center-of-mas and remove velocity and rot */
/* funcions: com2zero(), com_vw2zero()                     */
/* ******************************************************* */
int com_psinv(double *, double *);

/* ******************************************************* */
/* com2zero(x): move center-of-mass to zero                */
/* ******************************************************* */
int com2zero(REAL_T *x, REAL_T *minv)
{
   int n, i;
   REAL_T xs, ys, zs, m_t;

   n = 3 * prm->Natom;

   xs = ys = zs = 0.0;

   /* total mass */
   for (i = 0; i < n/3; i++)
   {
      m_t += 1./minv[3*i];
   }

   /* center of mass */
   for (i = 0; i < n/3; i++)
   {
      xs += x[3*i  ] / minv[3*i  ];
      ys += x[3*i+1] / minv[3*i+1];
      zs += x[3*i+2] / minv[3*i+2];
   }
   xs /= m_t;
   ys /= m_t;
   zs /= m_t;

//   printf("com: %f, %f, %f\n", xs, ys, zs);

   /* move com to center */
   for (i = 0; i < n/3; i++)
   {
      x[3*i  ] -= xs;
      x[3*i+1] -= ys;
      x[3*i+2] -= zs;
   }

   return (0);
}


/* ******************************************************* */
/* com_vw2zero(x): remove com velocity and rotation        */
/* ******************************************************* */
int com_vw2zero(REAL_T *x, REAL_T *v,REAL_T *minv)
{
   int n,i;
   REAL_T Lx,Ly,Lz;
   REAL_T vsx,vsy,vsz;
   REAL_T xx,xy,xz,yy,yz,zz;
   REAL_T I[9], Iinv[9];
   REAL_T Omx,Omy,Omz;
   REAL_T m_t;
   #define SCFAC 1.e-3

   n = 3 * prm->Natom;

   vsx = vsy = vsz = 0.0;

   /* total mass */
   for(i=0;i<n/3;i++)
   {
      m_t += 1./minv[3*i];
   }

   /* com velocity */
   for (i = 0; i < n/3; i++)
   {
      vsx += v[3*i  ] / minv[3*i  ];
      vsy += v[3*i+1] / minv[3*i+1];
      vsz += v[3*i+2] / minv[3*i+2];
   }
   vsx /= m_t;
   vsy /= m_t;
   vsz /= m_t;

//   printf("com velocity: %f, %f, %f\n", vsx, vsy, vsz);
 
   /* remove com velocity */ 
   for (i = 0; i < n/3; i++)
   {
      v[3*i  ] -= vsx;
      v[3*i+1] -= vsy;
      v[3*i+2] -= vsz;
   }


   /* Angular momentum */
   Lx = Ly = Lz = 0.0;

   for(i = 0; i < n/3; i++)
   {
      Lx += SCFAC*(x[3*i+1]*v[3*i+2] - x[3*i+2]*v[3*i+1])/minv[3*i];
      Ly += SCFAC*(x[3*i+2]*v[3*i  ] - x[3*i  ]*v[3*i+2])/minv[3*i];
      Lz += SCFAC*(x[3*i  ]*v[3*i+1] - x[3*i+1]*v[3*i  ])/minv[3*i];
   }

   /* tensor of inertia */

   xx = xy = xz = yy = yz = zz = 0.0;

   for(i = 0; i < n/3; i++)
   {
      xx += x[3*i  ]*x[3*i  ]/minv[3*i];
      xy += x[3*i  ]*x[3*i+1]/minv[3*i];
      xz += x[3*i  ]*x[3*i+2]/minv[3*i];
      yy += x[3*i+1]*x[3*i+1]/minv[3*i];
      yz += x[3*i+1]*x[3*i+2]/minv[3*i];
      zz += x[3*i+2]*x[3*i+2]/minv[3*i];
   }
   I[0] = (yy+zz)*SCFAC;
   I[1] = -xy*SCFAC;
   I[2] = -xz*SCFAC;
   I[3] = -xy*SCFAC;
   I[4] = (xx+zz)*SCFAC;
   I[5] = -yz*SCFAC;
   I[6] = -xz*SCFAC;
   I[7] = -yz*SCFAC;
   I[8] = (xx+yy)*SCFAC;

   /* inverse of tensor */
   com_psinv(I, Iinv);

   /* rotation: angular velocity */
   Omx = Omy = Omz = 0.0;

   Omx = Iinv[0]*Lx + Iinv[1]*Ly + Iinv[2]*Lz;
   Omy = Iinv[3]*Lx + Iinv[4]*Ly + Iinv[5]*Lz;
   Omz = Iinv[6]*Lx + Iinv[7]*Ly + Iinv[8]*Lz;

//   printf("com rotation: %f, %f, %f\n", Omx, Omy, Omz);

   /* remove com rotation */
   for (i = 0; i < n/3; i++)
   {
      v[3*i  ] = v[3*i  ] - Omy*x[3*i+2] + Omz*x[3*i+1];
      v[3*i+1] = v[3*i+1] - Omz*x[3*i  ] + Omx*x[3*i+2];
      v[3*i+2] = v[3*i+2] - Omx*x[3*i+1] + Omy*x[3*i  ];
   }

   return (0);
}

/* brute force inverse of 3x3 matrix */
int com_psinv(double *a, double *ainv)
{ 

   REAL_T det, minor_det[9]; 

   minor_det[0] = a[4] * a[8] - a[5] * a[7];
   minor_det[1] = a[3] * a[8] - a[5] * a[6];
   minor_det[2] = a[3] * a[7] - a[4] * a[6];
   minor_det[3] = a[1] * a[8] - a[2] * a[7];
   minor_det[4] = a[0] * a[8] - a[2] * a[6];
   minor_det[5] = a[0] * a[7] - a[1] * a[6];
   minor_det[6] = a[1] * a[5] - a[2] * a[4];
   minor_det[7] = a[0] * a[5] - a[2] * a[3];
   minor_det[8] = a[0] * a[4] - a[1] * a[3];

   det = a[0]* minor_det[0] - a[1] * minor_det[1] + a[2] * minor_det[2];

   ainv[0] =  minor_det[0] / det;
   ainv[1] = -minor_det[3] / det;
   ainv[2] =  minor_det[6] / det;
   ainv[3] = -minor_det[1] / det;
   ainv[4] =  minor_det[4] / det;
   ainv[5] = -minor_det[7] / det;
   ainv[6] =  minor_det[2] / det;
   ainv[7] = -minor_det[5] / det;
   ainv[8] =  minor_det[8] / det;

   return (0);
}
