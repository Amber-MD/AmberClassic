#ifndef GridStructs
#define GridStructs

#ifndef PREP_API
#include "fftw3.h"

#include "MatrixDS.h"
#endif

struct IntegerBook {
  int pag;
  int row;
  int col;
  int* data;
  int*** map;
};
typedef struct IntegerBook ibook;

struct FloatBook {
  int pag;                  // Number of pages (bins in x-dimension)
  int row;                  // Number of rows (bins in y-dimension)
  int col;                  // Number of columns (bins in z-dimension)
  int isortho;              // Flag to indicate that the grid is orthogonal
  double orig[3];           // The Cartesian origin of the grid
  dmat U;                   // The transformation matrix that puts everything
                            //   in units of the grid length vectors
  dmat invU;                // The inverse of U
  dmat L;                   // The transformation matrix that puts everything
                            //   in units of the grid BIN length vectors
  dmat invL;                // The inverse of L
  float* data;              // Grid information
  float*** map;             // Three-dimensional map into data
};
typedef struct FloatBook fbook;

struct DoubleBook {
  int pag;
  int row;
  int col;
  int pfft;
  double* data;
  fftw_complex* fdata;
  double*** map;
  fftw_complex*** fmap;
};
typedef struct DoubleBook dbook;

struct CharacterBook {
  int pag;
  int row;
  int col;
  char* data;
  char*** map;
};
typedef struct CharacterBook cbook;

struct CompressedChargeBook {
  ibook qdata;
  int nexcp;
  int maxexcp;
  double tq;
  int* eidx;
  double* eval;
};
typedef struct CompressedChargeBook qbook;

#endif
