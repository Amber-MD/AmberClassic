#ifndef VECTOR_FUNCS
#define VECTOR_FUNCS

int* CpyIVec(int* V, int n);

double* CpyDVec(double* V, int n);

void ReflectIVec(int* C, int* V, int n);

void ReflectDVec(double* C, double* V, int n);

void SetIVec(int* V, int N, int s);

void SetDVec(double* V, int N, double s);

int* CountUp(int N);

int ISum(int* V, int n);

int DetectInIVec(int im, int* V, int n);

double DSum(double* V, int n);

double DAbsSum(double* V, int n);

void Normalize(double* V, int N);

double DAverage(double* V, int n);

double DStDev(double* V, int n);

int ISum(int* V, int n);

double pythag(double a, double b);

double DExtreme(double* V, int N, int iflag);

int IExtreme(int* V, int N, int iflag);

double DotP(double* V1, double* V2, int N);

double Dot3(double* V1, double* V2);

double DMag(double* V, int n);

void UnitVector3(double* V);

void IVecAdd(int* V, int N, int S);

void DVecAdd(double* V, int N, double S);

void IVec2VecAdd(int* V, int* W, int N);

void DVec2VecAdd(double* V, double* W, int N);

void IVecMult(int* V, int N, int S);

void DVecMult(double* V, int N, double S);

double TrimAcos(double acosarg);

void CrossP(double* p, double* q, double* cr);

void Project(double* r, double* q, double* p, int dim);

double Dihedral(double* alocs, double* blocs, double* clocs, double* dlocs);

double Pearson(double* vec_a, double* vec_b, int num_t);

double SumTrace3(double* A);

double* PascalTriangle(int n);

double VecRMSD(double* va, double* vb, int n);

void UnitNormal2Plane(double* unv, double* A, double* B, double* C);

#endif
