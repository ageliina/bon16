#ifndef LAF05_H
#define LAF05_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#define LAF05_N_Z  5000
#define LAF05_N_NH 23

struct laf05 {
    gsl_spline2d *spline;
	gsl_interp_accel *xacc;
	gsl_interp_accel *yacc;
};

extern double laf05_corr[LAF05_N_Z][LAF05_N_NH];

void laf05_init(void);

void laf05_free(void);

double laf05_get_corr(double log_N_H, double z);

#endif
