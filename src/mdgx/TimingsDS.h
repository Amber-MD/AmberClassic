#ifndef TimingsStructs
#define TimingsStructs

#ifndef PREP_API
#include <sys/time.h>
#endif

struct ExecutionControl {
  struct timeval t0;
  struct timeval tC;
  struct timeval tti;
  struct timeval ttf;
  double bonds;
  double cellcomm;
  double nbInt;
  double nbDirAll;
  double nbBsp;
  double nbPtM;
  double nbFFT;
  double nbCnv;
  double nbMtP;
  double nbMtM;
  double nbRecAll;
  double Setup;
  double Integ;
  double Write;
  double Constraints;
  double Barostat;
  double Thermostat;
#ifdef MPI
  double mpiMeshPushWait;
  double mpiMeshPullWait;
  double mpiMeshPack;
#endif
};
typedef struct ExecutionControl execon;

#endif
