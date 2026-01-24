/*
** deform_energy.c
**
** DNA knowledge-based deformation energy calculations
** at the base-pair or base-pair step levels
**
** Thomas Gaillard <tgaillar@rci.rutgers.edu>
** David A. Case laboratory at Rutgers
**
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defreal.h"

/*
** Base-pair step parameters
** 
** Olson, W. K., Gorin, A. A., Lu, X. J., Hock, L. M. & Zhurkin, V. B. (1998).
** DNA sequence-dependent deformability deduced from protein-DNA crystal complexes.
** Proc. Natl. Acad. Sci. U.S.A. 95, 11163-11168.
*/

/*
** Protein-DNA step parameters average values
**       Shift     Slide      Rise      Tilt      Roll     Twist
**        0         1         2         3         4         5
*/

static float aver_step[4][4][6] = {
  {/*A*/
    {   -0.030,   -0.080,    3.270,   -1.400,    0.700,   35.100},/*AA*/
    {    0.000,   -0.590,    3.310,    0.000,    1.100,   29.300},/*AT*/
    {    0.090,   -0.250,    3.340,   -1.700,    4.500,   31.900},/*AG*/
    {    0.130,   -0.580,    3.360,   -0.100,    0.700,   31.500},/*AC*/
  },
  {/*T*/
    {    0.000,    0.050,    3.420,    0.000,    3.300,   37.800},/*TA*/
    {    0.030,   -0.080,    3.270,    1.400,    0.700,   35.100},/*TT*/
    {   -0.090,    0.530,    3.330,   -0.500,    4.700,   37.300},/*TG*/
    {    0.280,    0.090,    3.370,    1.500,    1.900,   36.300},/*TC*/
  },
  {/*G*/
    {   -0.280,    0.090,    3.370,   -1.500,    1.900,   36.300},/*GA*/
    {   -0.130,   -0.580,    3.360,    0.100,    0.700,   31.500},/*GT*/
    {    0.050,   -0.220,    3.420,   -0.100,    3.600,   32.900},/*GG*/
    {    0.000,   -0.380,    3.400,    0.000,    0.300,   33.600},/*GC*/
  },
  {/*C*/
    {    0.090,    0.530,    3.330,    0.500,    4.700,   37.300},/*CA*/
    {   -0.090,   -0.250,    3.340,    1.700,    4.500,   31.900},/*CT*/
    {    0.000,    0.410,    3.390,    0.000,    5.400,   36.100},/*CG*/
    {   -0.050,   -0.220,    3.420,    0.100,    3.600,   32.900} /*CC*/
  }
};

/*
** Protein-DNA step parameters force constants
**         Shift     Slide      Rise      Tilt      Roll     Twist
**          0         1         2         3         4         5
*/

static float fcst_step[4][4][6][6] = {
  {/*A*/
    {/*AA*/
      {    3.978,    1.559,    1.868,   -0.149,    0.017,    0.175},
      {    1.559,    6.160,    1.330,   -0.098,   -0.137,    0.022},
      {    1.868,    1.330,   21.748,   -0.288,    0.087,    0.072},
      {   -0.149,   -0.098,   -0.288,    0.100,   -0.001,    0.001},
      {    0.017,   -0.137,    0.087,   -0.001,    0.049,    0.031},
      {    0.175,    0.022,    0.072,    0.001,    0.031,    0.092}
    },
    {/*AT*/
      {    3.172,    0.000,    0.000,   -0.151,    0.000,    0.000},
      {    0.000,   10.694,    0.487,    0.000,   -0.060,   -0.090},
      {    0.000,    0.487,   25.547,    0.000,   -0.030,   -0.409},
      {   -0.151,    0.000,    0.000,    0.166,    0.000,    0.000},
      {    0.000,   -0.060,   -0.030,    0.000,    0.055,    0.029},
      {    0.000,   -0.090,   -0.409,    0.000,    0.029,    0.070}
    },
    {/*AG*/
      {    3.205,    1.380,    2.580,   -0.282,    0.003,    0.164},
      {    1.380,    7.187,    4.501,   -0.124,    0.100,    0.102},
      {    2.580,    4.501,   29.496,   -1.044,   -0.018,   -0.147},
      {   -0.282,   -0.124,   -1.044,    0.149,    0.026,    0.005},
      {    0.003,    0.100,   -0.018,    0.026,    0.096,    0.015},
      {    0.164,    0.102,   -0.147,    0.005,    0.015,    0.064}
    },
    {/*AC*/
      {    2.944,   -0.102,   -0.412,   -0.033,   -0.041,   -0.044},
      {   -0.102,    6.366,    2.510,   -0.010,    0.049,   -0.048},
      {   -0.412,    2.510,   23.860,    0.396,    0.245,   -0.276},
      {   -0.033,   -0.010,    0.396,    0.111,    0.014,    0.004},
      {   -0.041,    0.049,    0.245,    0.014,    0.080,    0.027},
      {   -0.044,   -0.048,   -0.276,    0.004,    0.027,    0.073}
    }
  },
  {/*T*/
    {/*TA*/
      {    3.860,    0.000,    0.000,   -0.176,    0.000,    0.000},
      {    0.000,    2.350,    0.234,    0.000,    0.061,   -0.067},
      {    0.000,    0.234,   21.914,    0.000,   -0.102,   -0.500},
      {   -0.176,    0.000,    0.000,    0.148,    0.000,    0.000},
      {    0.000,    0.061,   -0.102,    0.000,    0.029,    0.013},
      {    0.000,   -0.067,   -0.500,    0.000,    0.013,    0.052}
    },
    {/*TT*/
      {    3.978,   -1.559,   -1.868,   -0.149,   -0.017,   -0.175},
      {   -1.559,    6.160,    1.330,    0.098,   -0.137,    0.022},
      {   -1.868,    1.330,   21.748,    0.288,    0.087,    0.072},
      {   -0.149,    0.098,    0.288,    0.100,    0.001,   -0.001},
      {   -0.017,   -0.137,    0.087,    0.001,    0.049,    0.031},
      {   -0.175,    0.022,    0.072,   -0.001,    0.031,    0.092}
    },
    {/*TG*/
      {    3.733,   -0.871,   -1.747,   -0.037,   -0.018,    0.087},
      {   -0.871,    2.395,    2.369,   -0.062,    0.082,   -0.183},
      {   -1.747,    2.369,   18.235,   -0.189,    0.004,   -0.321},
      {   -0.037,   -0.062,   -0.189,    0.082,    0.001,    0.016},
      {   -0.018,    0.082,    0.004,    0.001,    0.048,    0.007},
      {    0.087,   -0.183,   -0.321,    0.016,    0.007,    0.043}
    },
    {/*TC*/
      {    6.542,   -1.750,   -0.148,   -0.132,    0.009,   -0.087},
      {   -1.750,    2.780,   -0.270,   -0.068,   -0.050,   -0.069},
      {   -0.148,   -0.270,   22.820,    0.502,    0.362,   -0.408},
      {   -0.132,   -0.068,    0.502,    0.087,    0.010,    0.000},
      {    0.009,   -0.050,    0.362,    0.010,    0.046,    0.011},
      {   -0.087,   -0.069,   -0.408,    0.000,    0.011,    0.071}
    }
  },
  {/*G*/
    {/*GA*/
      {    6.542,    1.750,    0.148,   -0.132,   -0.009,    0.087},
      {    1.750,    2.780,   -0.270,    0.068,   -0.050,   -0.069},
      {    0.148,   -0.270,   22.820,   -0.502,    0.362,   -0.408},
      {   -0.132,    0.068,   -0.502,    0.087,   -0.010,    0.000},
      {   -0.009,   -0.050,    0.362,   -0.010,    0.046,    0.011},
      {    0.087,   -0.069,   -0.408,    0.000,    0.011,    0.071}
    },
    {/*GT*/
      {    2.944,    0.102,    0.412,   -0.033,    0.041,    0.044},
      {    0.102,    6.366,    2.510,    0.010,    0.049,   -0.048},
      {    0.412,    2.510,   23.860,   -0.396,    0.245,   -0.276},
      {   -0.033,    0.010,   -0.396,    0.111,   -0.014,   -0.004},
      {    0.041,    0.049,    0.245,   -0.014,    0.080,    0.027},
      {    0.044,   -0.048,   -0.276,   -0.004,    0.027,    0.073}
    },
    {/*GG*/
      {    2.425,    0.309,    1.563,   -0.227,    0.119,    0.040},
      {    0.309,    3.542,    5.423,   -0.160,    0.143,   -0.061},
      {    1.563,    5.423,   30.312,   -0.851,    0.427,   -0.196},
      {   -0.227,   -0.160,   -0.851,    0.119,   -0.003,   -0.008},
      {    0.119,    0.143,    0.427,   -0.003,    0.064,   -0.001},
      {    0.040,   -0.061,   -0.196,   -0.008,   -0.001,    0.041}
    },
    {/*GC*/
      {    3.350,    0.000,    0.000,   -0.236,    0.000,    0.000},
      {    0.000,    6.244,    6.804,    0.000,    0.412,   -0.159},
      {    0.000,    6.804,   25.860,    0.000,    0.630,   -0.239},
      {   -0.236,    0.000,    0.000,    0.082,    0.000,    0.000},
      {    0.000,    0.412,    0.630,    0.000,    0.082,    0.005},
      {    0.000,   -0.159,   -0.239,    0.000,    0.005,    0.055}
    }
  },
  {/*C*/
    {/*CA*/
      {    3.733,    0.871,    1.747,   -0.037,    0.018,   -0.087},
      {    0.871,    2.395,    2.369,    0.062,    0.082,   -0.183},
      {    1.747,    2.369,   18.235,    0.189,    0.004,   -0.321},
      {   -0.037,    0.062,    0.189,    0.082,   -0.001,   -0.016},
      {    0.018,    0.082,    0.004,   -0.001,    0.048,    0.007},
      {   -0.087,   -0.183,   -0.321,   -0.016,    0.007,    0.043}
    },
    {/*CT*/
      {    3.205,   -1.380,   -2.580,   -0.282,   -0.003,   -0.164},
      {   -1.380,    7.187,    4.501,    0.124,    0.100,    0.102},
      {   -2.580,    4.501,   29.496,    1.044,   -0.018,   -0.147},
      {   -0.282,    0.124,    1.044,    0.149,   -0.026,   -0.005},
      {   -0.003,    0.100,   -0.018,   -0.026,    0.096,    0.015},
      {   -0.164,    0.102,   -0.147,   -0.005,    0.015,    0.064}
    },
    {/*CG*/
      {    1.586,    0.000,    0.000,   -0.137,    0.000,    0.000},
      {    0.000,    3.301,    0.191,    0.000,    0.005,    0.056},
      {    0.000,    0.191,   14.164,    0.000,   -0.020,   -0.162},
      {   -0.137,    0.000,    0.000,    0.068,    0.000,    0.000},
      {    0.000,    0.005,   -0.020,    0.000,    0.050,    0.024},
      {    0.000,    0.056,   -0.162,    0.000,    0.024,    0.047}
    },
    {/*CC*/
      {    2.425,   -0.309,   -1.563,   -0.227,   -0.119,   -0.040},
      {   -0.309,    3.542,    5.423,    0.160,    0.143,   -0.061},
      {   -1.563,    5.423,   30.312,    0.851,    0.427,   -0.196},
      {   -0.227,    0.160,    0.851,    0.119,    0.003,    0.008},
      {   -0.119,    0.143,    0.427,    0.003,    0.064,   -0.001},
      {   -0.040,   -0.061,   -0.196,    0.008,   -0.001,    0.041}
    }
  }
};

/*
** Base-pair parameters
** 
** Lankas F., Sponer J., Langowski J., Cheatham T. E. 3rd (2004).
** DNA deformability at the base pair level.
** J. Am. Chem. Soc. 126, 4124-4125.
*/

/*
** Base-pair parameters average values
**     Shear   Stretch   Stagger    Buckle Propeller   Opening
**       0         1         2         3         4         5
*/

static float aver_pair[4][6] = {
  {     0.01,     0.02,     0.01,    -1.60,   -10.30,     1.05},/*AT*/
  {    -0.01,     0.02,     0.01,     1.60,   -10.30,     1.05},/*TA*/
  {    -0.02,    -0.06,    -0.01,     0.10,    -7.82,    -0.03},/*GC*/
  {     0.02,    -0.06,    -0.01,    -0.10,    -7.82,    -0.03} /*CG*/
};

/*
** Base-pair parameters force constants
**       Shear   Stretch   Stagger    Buckle Propeller   Opening
**       0         1         2         3         4         5
*/

static float fcst_pair[4][6][6] = {
  {/*AT*/
    {   8.4574,  -0.1458,   0.0757,   0.0034,   0.0026,   0.0010},
    {  -0.1458,  42.2256,   0.1039,   0.0150,   0.0046,  -0.1313},
    {   0.0757,   0.1039,   4.0306,  -0.0026,  -0.0002,  -0.0320},
    {   0.0034,   0.0150,  -0.0026,   0.0066,   0.0001,  -0.0001},
    {   0.0026,   0.0046,  -0.0002,   0.0001,   0.0098,   0.0022},
    {   0.0010,  -0.1313,  -0.0320,  -0.0001,   0.0022,   0.0222}
  },
  {/*TA*/
    {   8.4574,   0.1458,  -0.0757,   0.0034,  -0.0026,  -0.0010},
    {   0.1458,  42.2256,   0.1039,  -0.0150,   0.0046,  -0.1313},
    {  -0.0757,   0.1039,   4.0306,   0.0026,  -0.0002,  -0.0320},
    {   0.0034,  -0.0150,   0.0026,   0.0066,  -0.0001,   0.0001},
    {  -0.0026,   0.0046,  -0.0002,  -0.0001,   0.0098,   0.0022},
    {  -0.0010,  -0.1313,  -0.0320,   0.0001,   0.0022,   0.0222}
  },
  {/*GC*/
    {   8.1138,  -1.0603,   0.1676,  -0.0014,   0.0030,   0.0276},
    {  -1.0603,  72.4022,  -1.6914,   0.0140,  -0.0998,  -1.2550},
    {   0.1676,  -1.6914,   5.9009,  -0.0197,   0.0638,   0.0813},
    {  -0.0014,   0.0140,  -0.0197,   0.0090,  -0.0004,  -0.0008},
    {   0.0030,  -0.0998,   0.0638,  -0.0004,   0.0105,   0.0057},
    {   0.0276,  -1.2550,   0.0813,  -0.0008,   0.0057,   0.0846}
  },
  {/*CG*/
    {   8.1138,   1.0603,  -0.1676,  -0.0014,  -0.0030,  -0.0276},
    {   1.0603,  72.4022,  -1.6914,  -0.0140,  -0.0998,  -1.2550},
    {  -0.1676,  -1.6914,   5.9009,   0.0197,   0.0638,   0.0813},
    {  -0.0014,  -0.0140,   0.0197,   0.0090,   0.0004,   0.0008},
    {  -0.0030,  -0.0998,   0.0638,   0.0004,   0.0105,   0.0057},
    {  -0.0276,  -1.2550,   0.0813,   0.0008,   0.0057,   0.0846}
  }
};

/*Index of nucleotide bases*/
static int baseindex(char base) {
  switch(base) {
    case 'A':return(0);
    case 'T':return(1);
    case 'G':return(2);
    case 'C':return(3);
    default:return(-1);
  }
}

/*Base-pair step level energy calculation*/
REAL_T step_ener(char *fname, int print) {

  FILE *fd;
  char line[100];
  REAL_T enertot = 0.0;
  float enerparamtot[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  int read = 0;
  int i;

  /*open 3DNA output file*/
  fd = fopen(fname,"r");
  if (fd == NULL) {
    fprintf(stderr,"error opening 3DNA output file %s!\n",fname);
    return -1;
  }

  if( print ){ /*print output header*/
    fprintf(stdout,"# DNA deformation energy calculations\n");
    fprintf(stdout,"# base-pair step level\n");
    fprintf(stdout,"#\n");
    fprintf(stdout,"# Energy by step\n");
    fprintf(stdout,"# num  step   Shift   Slide    Rise    Tilt    Roll   Twist   Total\n");
  }

  /*read 3DNA output file*/
  while (fgets(line,sizeof line,fd) != NULL) {
    
    /*stop reading*/
    if (read == 1 && strstr(line,"          ~~~") != NULL) {
      break;
    }
   
    /*read formatted line*/
    if (read == 1) {

      int k;
      char cb1, cb2, cb3, cb4;
      float param[6];
      int b1, b2;
      int i, j;
      float enerparam[6];
      float ener = 0.0;

      /*read parameters*/
      sscanf(line,"%4d %1c%1c/%1c%1c%10f%10f%10f%10f%10f%10f\n",&k,&cb1,&cb2,&cb3,&cb4,
        &param[0],&param[1],&param[2],&param[3],&param[4],&param[5]);
      b1 = baseindex(cb1);
      b2 = baseindex(cb2);

      /*calculate energy*/
      for (i=0; i<6; i++) {
        enerparam[i] = 0.5 * fcst_step[b1][b2][i][i] * (param[i] - aver_step[b1][b2][i]) * (param[i] - aver_step[b1][b2][i]);
        ener += enerparam[i];
        enerparamtot[i] += enerparam[i];
      }
      for (i=0; i<6; i++) {
        for(j=0; j<i; j++) {
          ener += fcst_step[b1][b2][i][j] * (param[i] - aver_step[b1][b2][i]) * (param[j] - aver_step[b1][b2][j]);
        }
      }
      enertot += ener;

      if( print ){ /*print step energy*/
        fprintf(stdout,"%5d %1c%1c/%1c%1c",k,cb1,cb2,cb3,cb4);
        for (i=0; i<6; i++) {
          fprintf(stdout,"%8.2f",enerparam[i]);
        }
        fprintf(stdout,"%8.2f\n",ener);
      }

    }

    /*begin reading*/
    if (strstr(line,"    step       Shift     Slide      Rise      Tilt      Roll     Twist") != NULL) {
      read = 1;
    }

  }
  fclose(fd);

  if( print ){ /*print total energy*/
    fprintf(stdout,"# Total step energy\n");
    fprintf(stdout,"#             Shift   Slide    Rise    Tilt    Roll   Twist   Total\n");
    fprintf(stdout,"#          ");
    for (i=0; i<6; i++) {
      fprintf(stdout,"%8.2f",enerparamtot[i]);
    }
    fprintf(stdout,"%8.2f\n",enertot);
  }

  return enertot;

}

/*Base-pair level energy calculation*/
REAL_T pair_ener(char *fname, int print) {

  FILE *fd;
  char line[100];
  REAL_T enertot = 0.0;
  float enerparamtot[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  int read = 0;
  int i;
  
  /*open 3DNA output file*/
  fd = fopen(fname,"r");
  if (fd == NULL) {
    fprintf(stderr,"error opening 3DNA output file %s\n",fname);
    return -1;
  }

  if( print ){ /*print output header*/
    fprintf(stdout,"# DNA deformation energy calculations\n");
    fprintf(stdout,"# base-pair level\n");
    fprintf(stdout,"#\n");
    fprintf(stdout,"# Energy by pair\n");
    fprintf(stdout,"# num pair  Shear  Stress Stagger  Buckle Propel. Opening   Total\n");
  }

  /*read 3DNA output file*/
  while (fgets(line,sizeof line,fd) != NULL) {
    
    /*stop reading*/
    if (read == 1 && strstr(line,"          ~~~") != NULL) {
      break;
    }
    
    /*read formatted line*/
    if (read == 1) {

      int k;
      char cb1, cb2;
      float param[6];
      int b1;
      int i, j;
      float enerparam[6];
      float ener;

      /*read parameters*/
      sscanf(line,"%5d %1c-%1c %10f%10f%10f%10f%10f%10f\n",&k,&cb1,&cb2,
        &param[0],&param[1],&param[2],&param[3],&param[4],&param[5]);
      b1 = baseindex(cb1);

      /*calculate energy*/
      ener = 0.0;
      for (i=0; i<6; i++) {
        enerparam[i] = 0.5 * fcst_pair[b1][i][i] * (param[i] - aver_pair[b1][i]) * (param[i] - aver_pair[b1][i]);
        ener += enerparam[i];
        enerparamtot[i] += enerparam[i];
      }
      for (i=0; i<6; i++) {
        for(j=0; j<i; j++) {
          ener += fcst_pair[b1][i][j] * (param[i] - aver_pair[b1][i]) * (param[j] - aver_pair[b1][j]);
        }
      }
      enertot += ener;

      if( print ){  /*print pair energy*/
        fprintf(stdout,"%5d %1c-%1c",k,cb1,cb2);
        for (i=0; i<6; i++) {
          fprintf(stdout,"%8.2f",enerparam[i]);
        }
        fprintf(stdout,"%8.2f\n",ener);
      }

    }

    /*begin reading*/
    if (strstr(line,"     bp        Shear    Stretch   Stagger    Buckle  Propeller  Opening") != NULL) {
      read = 1;
    }

  }
  fclose(fd);

  if( print ){  /*print total energy*/
    fprintf(stdout,"# Total pair energy\n");
    fprintf(stdout,"#           Shear  Stress Stagger  Buckle Propel. Opening   Total\n");
    fprintf(stdout,"#        ");
    for (i=0; i<6; i++) {
      fprintf(stdout,"%8.2f",enerparamtot[i]);
    }
    fprintf(stdout,"%8.2f\n",enertot);
  }

  return enertot;

}

