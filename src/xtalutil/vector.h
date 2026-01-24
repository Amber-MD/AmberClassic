#ifndef VECTOR_FUNCS
#define VECTOR_FUNCS

int* CpyIVec(int* V, int n);

int* CountUp(int N);

double DSum(double* V, int n);

void Normalize(double* V, int N);

double DAverage(double* V, int n);

double DStDev(double* V, int n);

int ISum(int* V, int n);

double pythag(double a, double b);

double* CpyDVec(double* V, int n);

double DExtreme(double* V, int N, int iflag);

int IExtreme(int* V, int N, int iflag);

void SetIVec(int* V, int N, int s);

double DotP(double* V1, double* V2, int N);

double DMag(double* V, int n);

void AddToIVec(int* V, int N, int S);

double TrimAcos(double acosarg);

void CrossP(double p[3], double q[3], double cr[3]);

void Project(double* r, double* q, double* p, int dim);

double Dihedral(double* alocs, double* blocs, double* clocs, double* dlocs);

double Pearson(double* vec_a, double* vec_b, int num_t);

#endif
