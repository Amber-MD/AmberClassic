#ifndef INC_AMBER_NETCDF_H
#define INC_AMBER_NETCDF_H

#include "AmberNetcdfDS.h"

// NOTE: To be NAB-useable, any functions added below must also be referenced
//       in nab/nabcode.h and nab/symbol.c
int netcdfDebug(struct AmberNetcdf*);
int netcdfLoad(struct AmberNetcdf*,char *);
int netcdfClose(struct AmberNetcdf*);
int netcdfWriteRestart(char *, int, double *, double *, double *, double,
		       double);
int netcdfCreate(struct AmberNetcdf*,char *, int, int);
int netcdfGetVelocity(struct AmberNetcdf*, int, double *);
int netcdfGetFrame(struct AmberNetcdf*, int, double*, double*);
int netcdfGetNextFrame(struct AmberNetcdf*,double*,double*);
int netcdfWriteFrame(struct AmberNetcdf*, int, double *, double *);
int netcdfWriteNextFrame(struct AmberNetcdf*, double *, double *);
int netcdfInfo(struct AmberNetcdf*);

#endif
