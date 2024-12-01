
#define	E_INTERNAL_ERROR	"Internal nab program error, exiting.\n"
#define	E_UNEXPECTED_SYM_D	"Unexpected symbol %d.\n"
#define	E_NOMEM_FOR_S		"Unable to allocate space for %s.\n"
#define	E_REDEF_BUILTIN_S	"Builtin symbol %s can't be redefined.\n"
#define	E_REDEF_SYM_S		"Redefinition of symbol %s.\n"
#define	E_NO_DECL_S		"Symbol %s not declared.\n"
#define	E_UNKNOWN_ATTR_S	"Unknown attrbute %s.\n"
#define	E_NOSUCH_STRAND_S	"Strand %s not in molecule.\n"
#define	E_NOSUCH_RESIDUE_S	"Residue %s.\n"
#define	E_NOSUCH_END_S		"Strand ends are \"first\", \"last\" not %s.\n"
#define	E_NOSUCH_ATOM_S		"Atom %s.\n"
#define	E_CANT_OPEN_RESLIB_S	"Can't open residue library %s.\n"
#define	E_CANT_OPEN_S		"Can't open file %s.\n"
#define	E_BAD_RESLIB_HEADER_S	"Incorrect line in residue library header %s...\n"
#define	E_BAD_BNDFILE_HEADER_S	"Incorrect header line in bond file: %s...\n"
#define	E_BAD_BNDFILE_DATA_S	"Incorrect data line in bond file: %s...\n"
#define	E_LIGATE_BAD_ENDS_S	"end1/end2 in mergestr() must be be first/last or last/first, not %s\n"
#define	E_MISSING_FIELDS_S	"First use of struct tag %s requires field list\n"
#define	E_STRUCT_TAG_EXPECTED_S	"ID %s must be a struct tag\n"
#define	E_NOFIELDS_ALLOWED_S	"Only first decl of %s can be followed by fields\n"

void	errormsg( int, char* );
void	errormsg_s( int, char*, char* );
void	errormsg_2s( int, char*, char*, char* );
void	errormsg_d( int, char*, int );
