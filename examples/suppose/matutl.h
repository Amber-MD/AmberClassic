/*  matutl.h    CCMATH mathematics library source code.
 *
 *  Copyright (C)  2000   Daniel A. Atkinson    All rights reserved.
 *  This code may be redistributed under the terms of the GNU library
 *  public license (LGPL). ( See the lgpl.license file for details.)
 * ------------------------------------------------------------------------
 */

/* Utility routines for real and complex matrix computations
 */
#ifndef CPX
struct complex {double re,im;};
typedef struct complex Cpx;
#define CPX 1
#endif
void trnm(double *a,int m);
void trncm(Cpx *s,int n);
void house(double *w,double *v,double *u,int m);
int qrevec(double *ev,double *v,double *d,int m);
int qreval(double *e,double *d,int m);
void hconj(Cpx *a,int m);
void chouse(Cpx *h,double *e,double *d,int m);
void chousv(Cpx *h,double *d, double *u,int m);
void qrecvc(double *e,Cpx *v,double *d,int m);
void cmmul(Cpx *a,Cpx *b,Cpx *c,int m);
void atou1(double *r,int i,int j);
void atovm(double *r,int i);
int qrbdu1(double *w,double *x,double *y,int k,double *z,int l);
void ldumat(double *x,double *y,int i,int k);
void ldvmat(double *x,double *y,int k);
int qrbdv(double *x,double *y,double *z,int i,double *w,int j);
int qrbdi(double *x,double *y,int l);
