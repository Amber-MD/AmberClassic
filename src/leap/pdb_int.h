
# include       "pdb.h"

# define        STREQN(s1, s2, n)       (strncmp(s1, s2, n) == 0)

struct  pdb_format      {
        char    *scan_format;
        char    *print_format;
};

extern  struct  pdb_format      pdb_record_format[];

extern  void    pdb_sscanf(char *, char *, ...);
extern  void    pdb_sprintf(char *, char *, ...);
