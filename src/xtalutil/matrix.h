#ifndef MATRIX_STRUCTS
#define MATRIX_STRUCTS

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

#ifndef MATRIX_FUNCS
#define MATRIX_FUNCS

imat CreateImat(int N, int M);

void DestroyImat(imat *A);

dmat CreateDmat(int N, int M);

void DestroyDmat(dmat *A);

cmat CreateCmat(int N, int M);

void DestroyCmat(cmat *A);

void DMatMult(dmat A, dmat B, dmat C);

void IMatMult(imat A, imat B, imat C);

void ttInv(dmat A, dmat invA);

void ttInvRow(dmat invA, double* pq, double* rs, int j, double deta);

void DMatAdd(dmat A, dmat B, double da, double db, dmat C);

void RotationMatrix(dmat *mat, double* vec, double angle);

void MatVecMult(dmat *A, double* B, double* C);

void AxbQRRxc(dmat A, double* b, int update_user);

void BackSub(dmat R, double* b);

void TRed2(dmat A, double* d, double* e);

void tqli(double* d, double* e, dmat z);

#endif
