#ifndef MOLUTIL_H
#define MOLUTIL_H

#include "nabtypes.h"

/* Functions defined in molutil.c in alphabetical order */

int addresidue (MOLECULE_T *mp, char *sname, RESIDUE_T *res);
int addstrand (MOLECULE_T *mp, char *sname);
REAL_T angle (MOLECULE_T *m, char *aex1, char *aex2, char *aex3);
REAL_T anglep (REAL_T *p1, REAL_T *p2, REAL_T *p3);
int cap (MOLECULE_T *mol, char *aex, int five, int three);
INT_T circle (REAL_T *p1, REAL_T *p2, REAL_T *p3, REAL_T *pc);
int connectres (MOLECULE_T *mol, char *sname, int ri, char *ainame, int rj, char *ajname);
EXTBOND_T *copyextbonds (RESIDUE_T *res);
RESIDUE_T *copyresidue (RESIDUE_T *);
int countmolatoms (MOLECULE_T *m, char *aex);
int countmolres (MOLECULE_T *m, char *aex);
int countmolstrands (MOLECULE_T *m, char *aex);
int countstrandresidues (MOLECULE_T *m, int strandnum);
REAL_T dist (MOLECULE_T *m, char *aex1, char *aex2);
REAL_T distp (REAL_T *pi, REAL_T *pj);
int freemolecule (MOLECULE_T *mol);
int freeparm (MOLECULE_T *mol);
int freeresidue (RESIDUE_T *res);
char *getresname (RESIDUE_T *res);
int mergestr (MOLECULE_T *mol1, char *strand1, char *end1, MOLECULE_T *mol2, char *strand2, char *end2);
int NAB_ainit (char **a, int s);
ATOM_T *NAB_anext (RESIDUE_T *res, ATOM_T *cap);
char **NAB_arc (ATOM_T *ap, char *key);
REAL_T *NAB_arf (ATOM_T *ap, char *key);
int *NAB_ari (ATOM_T *ap, char *key);
POINT_T *NAB_arp (ATOM_T *ap, char *key);
void NAB_initatom (ATOM_T *ap, int init_str);
void NAB_initres (RESIDUE_T *res, int init_str);
REF_MATRIX_T NAB_matcpy (REF_MATRIX_T mdst, REF_MATRIX_T msrc);
ATOM_T *NAB_mnext (MOLECULE_T *mol, ATOM_T *cap);
int *NAB_mri (MOLECULE_T *mol, char *key);
RESIDUE_T *NAB_rnext (MOLECULE_T *mol, RESIDUE_T *crp);
char **NAB_rrc (RESIDUE_T *res, char *key);
int *NAB_rri (RESIDUE_T *res, char *key);
MOLECULE_T *newmolecule (void);
REF_MATRIX_T newtransform (REAL_T dx, REAL_T dy, REAL_T dz, REAL_T rx, REAL_T ry, REAL_T rz);
REF_MATRIX_T rot4 (MOLECULE_T *mol, char *aex1, char *aex2, REAL_T angle);
REF_MATRIX_T rot4p (REAL_T *p1, REAL_T *p2, REAL_T angle);
int rt_errormsg_s (int, char *, char *);
int set_belly_mask (MOLECULE_T *m, char *aex, int *frozen);
int set_cons_mask (MOLECULE_T *m, char *aex, int *cons);
int setmol_from_xyz (MOLECULE_T **m, char **aex, REAL_T *xyz);
int setmol_from_xyzw (MOLECULE_T **m, char **aex, REAL_T *xyzw);
int setreskind (MOLECULE_T *m, char *aexp, char *rkind);
int setxyz_from_mol (MOLECULE_T **m, char **aex, POINT_T (*xyz));
int setxyzw_from_mol (MOLECULE_T **m, char **aex, REAL_T *xyzw);
REAL_T torsion (MOLECULE_T *mol, char *aei, char *aej, char *aek, char *ael);
REAL_T torsionp (REAL_T *pi, REAL_T *pj, REAL_T *pk, REAL_T *pl);
REF_MATRIX_T trans4 (MOLECULE_T *mol, char *aex1, char *aex2, REAL_T d);
REF_MATRIX_T trans4p (REAL_T *p1, REAL_T *p2, REAL_T d);
int transformmol (MATRIX_T mat, MOLECULE_T *mol, char *aexp);
int transformpts (MATRIX_T mat, POINT_T (*pts), int npts);
RESIDUE_T *transformres (MATRIX_T mat, RESIDUE_T *res, char *aexp);
void upd_molnumbers (MOLECULE_T *mp);
REF_MATRIX_T updtransform (MATRIX_T m1, MATRIX_T m2);

#endif

