#ifndef ThermostatFunctions
#define ThermostatFunctions

#include "TrajectoryDS.h"
#include "TopologyDS.h"
#include "CrdManipDS.h"
#include "CellManipDS.h"

double KineticEnergy(cellgrid *CG, coord *crd, prmtop *tp, trajcon *tj);

double KineticEnergyTI(cellgrid *CG, coord *crd, prmtop *tp, trajcon *tj);

double SystemTemperature(cellgrid *CG, coord *crd, prmtop *tp, Energy *sysUV,
			 trajcon *tj, int updateKE);

double SystemTemperatureTI(cellgrid* CG, coord* crd, prmtop* tp, Energy* sysUV,
			   trajcon *tj, int updateKE);

void BerendsenThermostat(cellgrid *CG, coord *crd, trajcon *tj, Energy *sysUV);

void BerendsenThermostatTI(cellgrid* CG, coord* crd, trajcon *tj,
                           Energy* sysUV);

void AndersenThermostat(cellgrid *CG, coord *crd, prmtop *tp, trajcon *tj,
			double T);

void NHThermostat(cellgrid *CG, coord *crd, prmtop *tp, trajcon *tj,
		  Energy *sysUV, double ETAt, double *CHIt, double *CHItp1e,
		  double *CHItp1q, int barostat);

void PrepThermoBarostat(prmtop *tp, trajcon *tj);

void PrepLangevinThermostat(trajcon *tj);

#endif
