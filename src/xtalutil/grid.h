#ifndef GRID_STRUCTS
#define GRID_STRUCTS

struct boolGrid {
  long long int lx;
  long long int ly;
  long long int lz;
  long long int total_data;
  long long int cdatasize;
  double gx;
  double gy;
  double gz;
  double invgx;
  double invgy;
  double invgz;
  double ox;
  double oy;
  double oz;
  double mx;
  double my;
  double mz;
  unsigned char* data;
};
typedef struct boolGrid bgrid;

struct intGrid {
  long long int lx;
  long long int ly;
  long long int lz;
  long long int total_data;
  double gx;
  double gy;
  double gz;
  double invgx;
  double invgy;
  double invgz;
  double ox;
  double oy;
  double oz;
  int* data;
  int*** map;
};
typedef struct intGrid igrid;

#endif

#ifndef GRID_FUNCS
#define GRID_FUNCS

#include "matrix.h"

int BiColorGrid(bgrid tg, double* crds, int natm, dmat U, dmat invU,
		double R, int updateUser);

void CheckGridDims(int lx, int ly, int lz, double ox, double oy, double oz,
		   double mx, double my, double mz);

bgrid CreateBoolGrid(int lx, int ly, int lz, double ox, double oy, double oz,
		     double mx, double my, double mz);

void WriteBoolGrid(bgrid tg, int i, int j, int k, int V);

int ReadBoolGrid(bgrid tg, int i, int j, int k);

igrid CreateIgrid(long long int lx, long long int ly, long long int lz,
                  double gx, double gy, double gz, double ox, double oy,
                  double oz);

#endif
