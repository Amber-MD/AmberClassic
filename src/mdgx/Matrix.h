#ifndef MatrixHeadings
#define MatrixHeadings

#include "MatrixDS.h"

imat CreateImat(int N, int M);

void DestroyImat(imat *A);

imat ReallocImat(imat *A, int M, int N);

dmat CreateDmat(int N, int M, int alstride);

void DestroyDmat(dmat *A);

dmat ReallocDmat(dmat *A, int M, int N);

void CopyDmat(dmat *Ac, dmat *A, int Acex);

void AplusBmat(dmat *A, dmat *B, int pom);

void AddDmat(dmat *A, dmat *B, double Ascl, double Bscl);

void DMatVecMult(dmat *A, double* x, double* b);

cmat CreateCmat(int N, int M);

void DestroyCmat(cmat *A);

cmat ReallocCmat(cmat *A, int nrow, int ncol);

void DMatMult(dmat A, dmat B, dmat C);

void IMatMult(imat A, imat B, imat C);

void ttInv(dmat A, dmat invA);

void DMatAdd(dmat A, dmat B, double da, double db, dmat C);

void RotationMatrix(dmat *mat, double* vec, double angle);

void AxbQRRxc(dmat A, double* b, int update_user);

void BackSub(dmat R, double* b);

int jacobi(double** a, int n, double* d, double** v, int *nrot);

void eigsrt(double* d, double** v, int n);

int DetectInCmat(cmat *C, char* s, int *row, int *col, int initl, int direc);

int CompareCmat(cmat *C1, cmat *C2);

int DistinguishImat(imat *I1, imat *I2);

dmat RotateA2B(double* A, double* B);

#endif
