#ifndef pmeDirectStructs
#define pmeDirectStructs

struct pmeDirectControlData {
  int LRvdw;         // Flag to activate long-ranged van-der Waals correction
  double Ecut;       // The electrostatic non-bonded cutoff
  double Vcut;       // The van-der Waals non-bonded cutoff
  double Mcut;       // The maximum of Ecut and Vcut
  double MaxDens;    // The maximum expected density of atoms (determines the
                     //   storage size in direct space decomposition cells)
  double invMcut;    // The inverse of Mcut
  double invEcut;    // The inverse of Ecut
  double ewcoeff;    // The Ewald coefficient, 0.5/sigma
  double sigma;      // The Gaussian width for spreading charges in preparation
                     //   for Ewald calculations
  double Dtol;       // The direct sum tolerance, the point at which the
                     //   difference in the interactions of Gaussian and point
                     //   charges is so small as to be deemed negligible
  double lkpspc;     // The spacing of the direct-space lookup table
};
typedef struct pmeDirectControlData dircon;

#endif
