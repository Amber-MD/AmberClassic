#ifndef MatrixStructs
#define MatrixStructs

#ifndef PREP_API
#include "fftw3.h"
#endif

struct IMatrix {
  int row;
  int col;
  int* data;
  int** map;
};
typedef struct IMatrix imat;

struct DMatrix {
  int row;
  int col;
  int stride;
  double* data;
  double** map;
};
typedef struct DMatrix dmat;

struct CMatrix {
  int row;
  int col;
  char* data;
  char** map;
};
typedef struct CMatrix cmat;

#endif
