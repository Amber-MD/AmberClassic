#ifndef LoopBuildStructures
#define LoopBuildStructures

#include "ConfigSampDS.h"

struct LoopControl {
  int count;          // The number of structures to manage simultaneously
  int nbundle;        // The number of bundles to consider when attempting to
                      //   link each side of a very stretched bond
  int npath;          // The number of paths to try in each bundle
  int verbose;        // Verbosity level
  int nbelly;         // The number of bellymasks to consider when making
                      //   the movable part of system
  int allmove;        // Flag to indicate that the user is masochistic and has
                      //   made the whole system flexible
  int pathlength;     // Number of knots in each path
  int nlongbond;      // The number of long or stretched bonds in the system
  double maxstretch;  // Maximum stretch that a bond may make before it is
                      //   considered a "long bond" and subject to special
                      //   treatment (default 2.0A beyond equilibrium length)
  double lreconst;    // The length to which a bond can be stretched out of
                      //   shape before it invokes bond splinting and healing
                      //   (default 8.0A beyond the equilibrium length)
  double bdwidth;     // The width of the space covered by all paths in any
                      //  single bundle
  int* movable;       // The array of atom mobilities
  cmat belly;         // The belly masks dictating the atom mobilities
  track* paths;       // Paths along which dummy particles can guide the ends
                      //   of a long bond to come together
};
typedef struct LoopControl loopcon;

#endif
