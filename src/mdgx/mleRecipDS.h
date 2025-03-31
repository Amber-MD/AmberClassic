#ifndef mleRecipStructs
#define mleRecipStructs

#ifndef PREP_API
#include "MatrixDS.h"
#endif

struct GridToGridMap {
  int ng;
  int ggordr;
  dmat s;
  imat m;
};
typedef struct GridToGridMap g2gmap;

#endif
