#ifndef PARMEDIT_STRUCTS
#define PARMEDIT_STRUCTS

struct ParameterEditor {
  int addEP;               // Flag to indicate that any extra points added to
                           //   the original topology should be included in the
                           //   output (by default, this is set to 1, ON)
  int topchk;              // Flag to engage sanity checks on the new topology
  double qtol;             // Tolerance for residue charge sanity check
  char prmfile[MAXNAME];   // Name of the topology file to write
  char crdfile[MAXNAME];   // Name of the input coordinates file to write
  char prmtitle[MAXNAME];  // Title of the topology file to write
};
typedef struct ParameterEditor parmed;

#endif
