#ifndef LEBEDEV_H
#define LEBEDEV_H

int available_table ( int rule );
int gen_oh ( int code, double a, double b, double v, double *x, 
  double *y, double *z, double *w );
void ld_by_order ( int order, double *x, double *y, double *z, double *w );
void ld0006 ( double *x, double *y, double *z, double *w );
void ld0014 ( double *x, double *y, double *z, double *w );
void ld0026 ( double *x, double *y, double *z, double *w );
void ld0038 ( double *x, double *y, double *z, double *w );
void ld0050 ( double *x, double *y, double *z, double *w );
void ld0074 ( double *x, double *y, double *z, double *w );
void ld0086 ( double *x, double *y, double *z, double *w );
void ld0110 ( double *x, double *y, double *z, double *w );
void ld0146 ( double *x, double *y, double *z, double *w );
void ld0170 ( double *x, double *y, double *z, double *w );
void ld0194 ( double *x, double *y, double *z, double *w );
void ld0230 ( double *x, double *y, double *z, double *w );
void ld0266 ( double *x, double *y, double *z, double *w );
void ld0302 ( double *x, double *y, double *z, double *w );
void ld0350 ( double *x, double *y, double *z, double *w );
void ld0434 ( double *x, double *y, double *z, double *w );
void ld0590 ( double *x, double *y, double *z, double *w );
void ld0770 ( double *x, double *y, double *z, double *w );
void ld0974 ( double *x, double *y, double *z, double *w );
void ld1202 ( double *x, double *y, double *z, double *w );
void ld1454 ( double *x, double *y, double *z, double *w );
void ld1730 ( double *x, double *y, double *z, double *w );
void ld2030 ( double *x, double *y, double *z, double *w );
void ld2354 ( double *x, double *y, double *z, double *w );
void ld2702 ( double *x, double *y, double *z, double *w );
void ld3074 ( double *x, double *y, double *z, double *w );
void ld3470 ( double *x, double *y, double *z, double *w );
void ld3890 ( double *x, double *y, double *z, double *w );
void ld4334 ( double *x, double *y, double *z, double *w );
void ld4802 ( double *x, double *y, double *z, double *w );
void ld5294 ( double *x, double *y, double *z, double *w );
void ld5810 ( double *x, double *y, double *z, double *w );
int order_table ( int rule );
int precision_table ( int rule );
void timestamp ( );
void xyz_to_tp ( double x, double y, double z, double *t, double *p );

#endif
