#ifndef UED14_H
#define UED14_H

#include <math.h>
#include <stdio.h>

/* Helper functions */
#define min(x,y) (((x) < (y)) ? (x) : (y))
#define max(x,y) (((x) > (y)) ? (x) : (y))
#define minmax(lo,mi,hi) min(max(mi,lo),hi)

/* Ueda+2014 parameters */
#define PSI_MIN 0.20
#define PSI_MAX 0.84
#define PSI_4375_0 0.43
#define BETA 0.24
#define EPSILON 1.7
#define A1 0.48
#define F_CTK 1.0

double ued14_get_Psi_4375(double z);
double ued14_get_Psi(double log_L_X, double z);
double ued14_get_f(double log_L_X, double z, double log_N_H);

#endif
