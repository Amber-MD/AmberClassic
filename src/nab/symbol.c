#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "nab.h"
#include "errormsg.h"
#include "symbol.h"

	/* Symbol table:	*/
	/* this MUST be ALPHABETIC on the 1st field:	*/

static	SYMREC_T	ssyms[] = {
	{ "EOF", T_INT, C_DEFINE, K_SCALAR, S_SYSTEM, CC_UNDEF,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "MAT_cube", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "MAT_cyclic", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "MAT_dihedral", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "MAT_fprint", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "MAT_fscan", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "MAT_getsyminfo", T_STRING, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "MAT_helix", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "MAT_ico", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "MAT_octa", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "MAT_orient", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "MAT_rotate", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "MAT_sprint", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "MAT_sscan", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "MAT_tetra", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "MAT_translate", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "NULL", T_NULL, C_NULL, K_SCALAR, S_SYSTEM, CC_UNDEF,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "acos", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "addresidue", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "addstrand", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "alignframe", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "andbounds", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "angle", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "anglep", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "argc", T_INT, C_VAR, K_SCALAR, S_SYSTEM, CC_UNDEF,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "argv", T_STRING, C_VAR, K_ARRAY, S_SYSTEM, CC_UNDEF,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "asin", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "atan", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "atan2", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "atof", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "atoi", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "axis2frame", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "bdna", T_MOLECULE, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "bonded_atoms", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "cap", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "ceil", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "circle", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "conjgrad", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "connectres", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "copymolecule", T_MOLECULE, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "cos", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "cosh", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "countmolatoms", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "countmolres", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "countmolstrands", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "countstrandresidues", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "date", T_STRING, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "db_viol", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "db_viol3", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "dg_helix", T_MOLECULE, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "dg_options", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "dist", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "distp", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "dt_to_bmat", T_BOUNDS, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "dt_to_prmtop", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "dumpatom", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "dumpbounds", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "dumpboundsviolations", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "dumpchiviolations", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "dumpmatrix", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "dumpmolecule", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "dumpresidue", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "embed", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "exit", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "exp", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "fabs", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "fclose", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "fd_helix", T_MOLECULE, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
	     FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "fflush", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "floor", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "fmod", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "fopen", T_FILE, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "fprintf", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "free", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "freemolecule", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "freeparm", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "freeresidue", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "fscanf", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_IO,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
#ifndef WIN32
	{ "ftime", T_STRING, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
#endif
	{ "gauss2", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "geodesics", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC, 
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "getchivol", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "getchivolp", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "getcif", T_MOLECULE, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "getcompcif", T_MOLECULE, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "getenv", T_STRING, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "getline", T_STRING, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "getmatrix", T_MATRIX, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "getpdb", T_MOLECULE, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "getpdb_prm", T_MOLECULE, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "getpdb_rlb", T_MOLECULE, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "getres", T_RESIDUE, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "getresidue", T_RESIDUE, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "getreslibkind", T_STRING, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "getresname", T_STRING, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "getseq_from_pdb", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "getxv", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "getxyz", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "getxyz_from_pdb", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "gsub", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "helixanal", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "index", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "length", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "link_na", T_MOLECULE, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "linkprot", T_MOLECULE, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "lmod", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "lmodC", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "lmod_opt_init", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "log", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "log10", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "md", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "mergestr", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "metrize", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "mm_options", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "mm_set_checkpoint", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "mme", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "mme2", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "mme4", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "mme_init", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "mme_rattle", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "mme_rism_max_memory", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "mme_timer", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "molsurf", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "mpierror", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "mpifinalize", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "mpiinit", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "mytaskid", T_INT, C_VAR, K_SCALAR, S_SYSTEM, CC_UNDEF,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "nabout", T_FILE, C_VAR, K_SCALAR, S_SYSTEM, CC_UNDEF,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "netcdfClose", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "netcdfCreate", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "netcdfDebug", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "netcdfGetFrame", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "netcdfGetNextFrame", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "netcdfGetVelocity", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
                FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "netcdfInfo", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "netcdfLoad", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "netcdfWriteFrame", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "netcdfWriteNextFrame", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
        { "netcdfWriteRestart", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
                FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "newbounds", T_BOUNDS, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "newmolecule", T_MOLECULE, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "newton", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "newtransform", T_MATRIX, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "nmode", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "numtasks", T_INT, C_VAR, K_SCALAR, S_SYSTEM, CC_UNDEF,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "orbounds", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "pair_ener", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "plane", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "pow", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "printf", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "putarc", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "putbnd", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "putcif", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "putdist", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "putmatrix", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "putpdb", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "putx", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "putxv", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "putxyz", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "rand2", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "readbinposfrm", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "readbinposhdr", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "readparm", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "rmsd", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "rot4", T_MATRIX, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "rot4p", T_MATRIX, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "rseed", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "safe_fopen", T_FILE, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "sasad", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "scanf", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_IO,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "setbounds", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "setboundsfromdb", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "setchiplane", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "setchivol", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "setframe", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "setframep", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "setmol_from_xyz", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "setmol_from_xyzw", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "setpoint", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "setreskind", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "setreslibkind", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "setseed", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "setxyz_from_mol", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "setxyzw_from_mol", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "showbounds", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "sin", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "sinh", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "split", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "split_n", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "sprintf", T_STRING, C_FUNC, K_SCALAR, S_SYSTEM, CC_IO,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "sqrt", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "sscanf", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_IO,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "stderr", T_FILE, C_DEFINE, K_SCALAR, S_SYSTEM, CC_UNDEF,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "stdin", T_FILE, C_DEFINE, K_SCALAR, S_SYSTEM, CC_UNDEF,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "stdout", T_FILE, C_DEFINE, K_SCALAR, S_SYSTEM, CC_UNDEF,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "step_ener", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "sub", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "substr", T_STRING, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "sugarpuckeranal", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "superimpose", T_MATRIX, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "system", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "tan", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "tanh", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "timeofday", T_STRING, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "torsion", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "torsionp", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "trans4", T_MATRIX, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "trans4p", T_MATRIX, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "transformmol", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "transformpts", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "transformres", T_RESIDUE, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "tsmooth", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "tss_init", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "tss_get", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "tss_next", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "tss_read", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "tss_set", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "tss_write", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "unlink", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "useboundsfrom", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "usemodeldist", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "wc_basepair", T_MOLECULE, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "wc_complement", T_STRING, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "wc_helix", T_MOLECULE, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "writebinposfrm", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "writebinposhdr", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_CC,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "xmin", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "xminC", T_FLOAT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 },
	{ "xmin_opt_init", T_INT, C_FUNC, K_SCALAR, S_SYSTEM, CC_NAB,
		FALSE, FALSE, FALSE, 0, NULL, NULL, NULL, NULL, NULL, 0 }
};
static	int	n_ssyms = sizeof( ssyms ) / sizeof( SYMREC_T );

int	n_gsyms = 0;
static	SYMREC_T	*gsyms = NULL;

static	int	n_lsyms = 0;
static	SYMREC_T	*lsyms = NULL;

static	int	n_usyms = 0;
static	int	u_init = FALSE;
static	SYMREC_T	*usyms = NULL;
static	SYMREC_T	*stag = NULL;

	/* Error Msg:		*/

static	char	e_msg[ 256 ];

	/* current function name	*/
char	c_funcname[ 256 ] = "main";

static	void	dsymtree( FILE *, char [], SYMREC_T * );
static	SYMREC_T	*nsym( char [], int, int, int, int, int );
static	SYMREC_T	*esym( SYMREC_T *, SYMREC_T * );
static	SYMREC_T	*fsym( SYMREC_T *, char [] );
static	SYMREC_T	*fssym( char [] );

SYMREC_T	*entersym( int scope, char name[],
	int type, int clas, int kind, int cconv )
{
	SYMREC_T	*sp;

	if( ( sp = nsym( name, type, clas, kind, scope, cconv ) ) == NULL )
		return( NULL );
	if( scope == S_USER ){
		usyms = esym( usyms, sp );
		u_init |= sp->s_type == T_STRING;
	}else if( scope == S_LOCAL )
		lsyms = esym( lsyms, sp );
	else
		gsyms = esym( gsyms, sp );
	return( sp );
}

SYMREC_T	*findsym( char name[] )
{
	SYMREC_T	*sp;

	if( ( sp = fsym( usyms, name )) ){
		return( sp );
	}else if( (sp = fsym( lsyms, name )) ){
		return( sp );
	}else if( (sp = fsym( gsyms, name )) ){
		return( sp );
	}else if( (sp =  fssym( name )) ){
		return( sp );
	}else{
		return( NULL );
	}
}

int	openusyms( SYMREC_T *st )
{
	int		ok = TRUE;

	stag = st;
	usyms = st->s_user;
	n_usyms = st->s_nusyms;
	u_init = FALSE;
	return( ok );
}

void	closeusyms( SYMREC_T *st )
{

	if( st != NULL ){
		st->s_init = u_init;
		st->s_user = usyms;
		st->s_nusyms = n_usyms;
	}
	u_init = FALSE;
	usyms = NULL;
	n_usyms = 0;
	stag = NULL;
}

void	freelsyms( void )
{

	lsyms = NULL;
	n_lsyms = 0;
}

void	dumpsyms( FILE *fp, int d_sys, int d_glb, int d_loc )
{
	int	s;

	if( d_sys ){
	fprintf( fp, "%3d System Symbols:\n", n_ssyms );
		for( s = 0; s < n_ssyms; s++ )
			dsym( fp, "", &ssyms[ s ] );
	}
	if( d_glb ){
		fprintf( fp, "%3d Global Symbols:\n", n_gsyms );
		dsymtree( fp, "", gsyms );
	}
	if( d_loc ){
		fprintf( fp, "%3d Local Symbols:\n", n_lsyms );
		dsymtree( fp, "", lsyms );
	}
}

void	dsym( FILE *fp, char offset[], SYMREC_T *sp )
{
	int	p;

	fprintf( fp, "%s\t%s\n", offset, sp->s_name );
	fprintf( fp, "%s\t\tType:\t", offset );
	switch( sp->s_type ){
	case T_UNDEF :
		fprintf( fp, "Undef\n"  );
		break;
	case T_INT :
		fprintf( fp, "Int\n"  );
		break;
	case T_SIZE_T :
		fprintf( fp, "Size_t\n" );
		break;
	case T_FLOAT :
		fprintf( fp, "Float\n"  );
		break;
	case T_STRING :
		fprintf( fp, "String\n"  );
		break;
	case T_FILE :
		fprintf( fp, "File\n"  );
		break;
	case T_POINT :
		fprintf( fp, "Point\n"  );
		break;
	case T_MATRIX :
		fprintf( fp, "Matrix\n"  );
		break;
	case T_RESIDUE :
		fprintf( fp, "Residue\n"  );
		break;
	case T_MOLECULE :
		fprintf( fp, "Molecule\n"  );
		break;
	case T_BOUNDS :
		fprintf( fp, "Bounds\n" );
		break;
	case T_USER :
		fprintf( fp, "User\n" );
		break;
	default :
		fprintf( fp, "?? (%d)\n", sp->s_type );
		break;
	}
	fprintf( fp, "%s\t\tClass:\t", offset );
	switch( sp->s_class ){
	case C_UNDEF :
		fprintf( fp, "Undef\n"  );
		break;
	case C_LIT :
		fprintf( fp, "Literal\n"  );
		break;
	case C_STRTAG :
		fprintf( fp, "StrTag\n" );
		break;
	case C_VAR :
		fprintf( fp, "Variable\n"  );
		break;
	case C_FUNC :
		fprintf( fp, "Function\n" );
		break;
	case C_EXPR :
		fprintf( fp, "Expression\n" );
		break;
	default :
		fprintf( fp, "?? (%d)\n", sp->s_class );
		break;
	}
	fprintf( fp, "%s\t\tKind:\t", offset );
	switch( sp->s_kind ){
	case K_UNDEF :
		fprintf( fp, "Undef\n"  );
		break;
	case K_SCALAR :
		fprintf( fp, "Scalar\n"  );
		break;
	case K_ARRAY :
		fprintf( fp, "Array\n"  );
		break;
	case K_DARRAY :
		fprintf( fp, "Dynamic Array\n"  );
		break;
	case K_HASHED :
		fprintf( fp, "Hashed\n"  );
		break;
	default :
		fprintf( fp, "?? (%d)\n", sp->s_kind );
		break;
	}
	fprintf( fp, "%s\t\tCconv:\t", offset );
	switch( sp->s_kind ){
	case CC_UNDEF :
		fprintf( fp, "Undef\n"  );
		break;
	case CC_NAB :
		fprintf( fp, "nab\n"  );
		break;
	case CC_CC :
		fprintf( fp, "C\n"  );
		break;
	case CC_FORTRAN :
		fprintf( fp, "Fortran\n"  );
		break;
	case CC_IO :
		fprintf( fp, "I/O\n"  );
		break;
	default :
		fprintf( fp, "?? (%d)\n", sp->s_cconv );
		break;
	}
	fprintf( fp, "%s\t\tisdyn  = %d\n", offset, sp->s_isdyn );
	fprintf( fp, "%s\t\tisparm = %d\n", offset, sp->s_isparm );
	fprintf( fp, "%s\t\tinit   = %d\n", offset, sp->s_init );
	fprintf( fp, "%s\t\tpcount = %d\n", offset, sp->s_pcount );
	if( sp->s_type == T_USER ){
		fprintf( fp, "%s\t\tUser  %s = {\n", offset,
			sp->s_uname ? sp->s_uname : "(NONE)" );
		dsymtree( fp, "\t\t", sp->s_user );
		fprintf( fp, "%s\t\t}\n", offset );
	}
	if( sp->s_parts ){
		fprintf( fp, "%s\t\tParts:  %d\n", offset, sp->s_pcount );
		for( p = 0; p < sp->s_pcount; p++ )
			dumpexpr( fp, sp->s_parts[ p ], 0 );
	}
	fprintf( fp, "%s\t}\n", offset );
}

static	void	dsymtree( FILE *fp, char offset[], SYMREC_T *sp )
{

	if( sp ){
		dsymtree( fp, offset, sp->s_left );
		dsym( fp, offset, sp );
		dsymtree( fp, offset, sp->s_right );
	}
}

static	SYMREC_T	*nsym( char name[],
	int type, int clas, int kind, int scope, int cconv )
{
	SYMREC_T	*sp;
	int		len;
	char		*snp;

	if( ( sp = ( SYMREC_T * )malloc( sizeof( SYMREC_T ) ) ) == NULL ){
		sprintf( e_msg, "SYMBOL %s", name );
		errormsg_s( TRUE, E_NOMEM_FOR_S, e_msg );
		return( NULL );
	}
	len = strlen( name ) + 1;
	if( ( snp = ( char * )malloc( len * sizeof(char) ) ) == NULL ){
		sprintf( e_msg, "SYMBOL name %s", name );
		errormsg_s( TRUE, E_NOMEM_FOR_S, e_msg );
		return( NULL );
	}else
		strcpy( snp, name );
	sp->s_name = snp;
	sp->s_type = type;
	sp->s_class = clas;
	sp->s_kind = kind;
	sp->s_scope = scope;
	sp->s_cconv = cconv;
	if( scope == S_GLOBAL )
		n_gsyms++;
	else if( scope == S_LOCAL )
		n_lsyms++;
	else
		n_usyms++;
	sp->s_isdyn = FALSE;
	sp->s_isparm = FALSE;
	sp->s_init = FALSE;
	sp->s_pcount = 0;
	sp->s_parts = NULL;
	sp->s_left = NULL;
	sp->s_right = NULL;
	sp->s_uname = NULL;
	sp->s_user = NULL;
	sp->s_nusyms = 0;
	return( sp );
}

static	SYMREC_T	*esym( SYMREC_T *root, SYMREC_T *nsp )
{
	int	cv;

	if( root == NULL )
		root = nsp;
	else if( ( cv = strcmp( root->s_name, nsp->s_name ) ) < 0 )
		root->s_right = esym( root->s_right, nsp );
	else if( cv > 0 )
		root->s_left = esym( root->s_left, nsp );

	return( root );
}

static	SYMREC_T	*fsym( SYMREC_T *root, char name[] )
{
	int	cv;

	if( root ){
		if( ( cv = strcmp( root->s_name, name ) ) == 0 ){
			return( root );
		}else if( cv < 0 )
			return( fsym( root->s_right, name ) );
		else
			return( fsym( root->s_left, name ) );
	}else
		return( NULL );
}

static	SYMREC_T	*fssym( char name[] )
{
	int	i, j, k;
	int	cv;
	SYMREC_T	*sp;

	for( i = 0, j = n_ssyms - 1; i <= j; ){
		k = ( i + j ) / 2;
		sp = &ssyms[ k ];
		if( ( cv = strcmp( sp->s_name, name ) ) == 0 )
			return( sp );
		else if( cv < 0 )
			i = k + 1;
		else
			j = k - 1;
	}
	return( NULL );
}
