/*  eigen.c    CCMATH mathematics library source code.
 *
 *  Copyright (C)  2000   Daniel A. Atkinson    All rights reserved.
 *  This code may be redistributed under the terms of the GNU library
 *  public license (LGPL). ( See the lgpl.license file for details.)
 * ------------------------------------------------------------------------
 */
#include <stdlib.h>
#include "matutl.h"
#include "ccmath.h"

void eigen(double *a, double *ev, int *n)
{
	double *dp;
	int n_1;
	n_1 = (*n);

	dp = (double *) calloc(n_1, sizeof(double));
	housev(a, ev, dp, n_1);
	qrevec(ev, a, dp, n_1);
	trnm(a, n_1);
	free(dp);
}
