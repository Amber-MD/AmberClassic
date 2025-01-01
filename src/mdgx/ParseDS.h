#ifndef ParseDataStructures
#define ParseDataStructures

struct alpha_index {
  int pos;       // Position of the string in a character matrix
  int numeric;   // Integer obtained by condensing all digits from the string
  char* word;    // Pointer to the character string itself
};
typedef struct alpha_index WordNdex;

struct RegularExpression {
  int nparts;    // The number of parts that must be matched
  int matchhead; // Flag to indicate that the initial characters must match
  int matchtail; // Flag to indicate the final characters must match
  cmat* parts;   // The matrices of characters in each segment to match.
                 //   For instance, *Min[ABC]Reg??*[xyZ] has two parts.  The
                 //   first part MUST read "Min" followed by any of the letters
                 //   A, B, or C, followed by "Reg", followed by any two
                 //   characters as long as there are at least two.  The second
                 //   part must consist of one of the characters x, y, or Z AND
                 //   terminate the string.  In this case matchhead would
                 //   be zero to indicate that no, we don't have to match the
                 //   start of the string, and matchtail would be 1.
  imat noptions; // The number of options available for each character in
                 //   any of the parts
};
typedef struct RegularExpression regexp;

struct NamelistGroup {
  int count;     // The number of namelists encountered
  cmat titles;   // The titles of all namelists (e.g. &cntrl)
  cmat* nml;     // The array of namelists
};
typedef struct NamelistGroup nmlgroup;

#endif
