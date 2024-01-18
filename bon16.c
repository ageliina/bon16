/*
 * Implements the equations in Bon+16 for the bivariate distribution
 */

#include <stdio.h>
#include <math.h>

#include "dz.h"

#include "bon16.h"
#include "cap09.h"
#include "laf05.h"
#include "ued14.h"

/* For some reason pow10 is not defined on my system so define it here */
#define pow10(x) pow(10,(x))

double get_log_f_star(double log_M_star, double z, struct args_f_star args)
{
    double log_M_star_star = args.log_M_star_star + (z >= args.z1) * args.k_log_M_star_star * (z - args.z1);
    double log_x = log_M_star - log_M_star_star;
    return args.alpha * log_x - pow10(log_x) * log10(exp(1.));
}

double get_log_f_lambda_SAR(double log_lambda_SAR, double log_M_star, double z, struct args_f_lambda_SAR args)
{
    double log_lambda_SAR_star = args.log_lambda_SAR_star0 + args.k_lambda * (log_M_star - args.log_M_star0);
    double gamma1 = (z < args.z1) \
        ? args.gamma10 + args.k_gamma * (z - args.z0) \
        : args.gamma10 + args.k_gamma * (args.z1 - args.z0) + args.k_gamma1 * (z - args.z1);
    double x = pow10(log_lambda_SAR - log_lambda_SAR_star);
    return -1 * log10(pow(x, -gamma1) + pow(x, -args.gamma2));
}

double get_log_f_z(double z, struct args_f_z args)
{

    double xs[3] = { log10(1 + z), log10(1 + args.z0), log10(1 + args.z1) };
    double ps[3] = { args.p1, args.p2, args.p3 };

    double ret = NAN;
    switch ((z >= args.z0) + (z >= args.z1)) {
        case 0: ret = ps[0] * xs[0]; break;
        case 1: ret = ps[0] * xs[1] + ps[1] * (xs[0] - xs[1]); break;
        case 2: ret = ps[0] * xs[1] + ps[1] * (xs[2] - xs[1]) + ps[2] * (xs[0] - xs[2]); break;
    }
    return ret;
}

double get_log_Psi(double log_M_star, double log_lambda_SAR, double z, struct args_Psi args)
{
    return (
            args.log_Psi_star
            + get_log_f_lambda_SAR(log_lambda_SAR, log_M_star, z, args.args_l)
            + get_log_f_star(log_M_star, z, args.args_s)
            + get_log_f_z(z, args.args_z)
           );
}

double get_Phi_star(double x[], size_t dim, void *p)
{
    (void) (dim);
    struct args_Phi_star *args = (struct args_Phi_star *) p;
    return pow10(get_log_Psi(args->log_M_star, x[0], args->z, args->args));
}

double get_Phi_SAR(double x[], size_t dim, void *p)
{
    (void) (dim);
    struct args_Phi_SAR *args = (struct args_Phi_SAR *) p;
    return pow10(get_log_Psi(x[0], args->log_lambda_SAR, args->z, args->args));
}

double get_Phi_X(double x[], size_t dim, void *p)
{
    (void) (dim);
    struct args_Phi_X *args = (struct args_Phi_X *) p;
    double log_L_X = x[0] + x[1];
    double dlog_L_X = args->log_L_X_max - args->log_L_X_min;
    char is_in = (args->log_L_X_min <= log_L_X) * (log_L_X < args->log_L_X_max);
    return is_in ? pow10(get_log_Psi(x[0], x[1], args->z, args->args)) / dlog_L_X : 0;
}

double get_N_above_S(double x[], size_t dim, void *p)
{

    /*
     * NOTE: the "heavy" calculations are done by the lookup table, created by
     * 'create_dz.py'. It contains the values of 4 pi d_L(z)**2 and deltaV,
     * where d_L is the luminosity distance and deltaV the differential
     * comoving volume.
     */

    (void) (dim);
    struct args_N_above_S *args = (struct args_N_above_S *) p;
    double Nsum = 0;

    for (size_t i=0; i<NZ; i++) {

        double z = 0.5 * (DZ[i][0] + DZ[i][1]);
        double dV = 0.5 * (DZ[i][6] - DZ[i][5]);

        double l = x[0] + x[1];
        double l_lo = args->log_F_X_lo + 2 * log10(DZ[i][4]) + 50.077910954466375 - (2 - args->Gamma) * log10(1 + z);
        double l_hi = args->log_F_X_hi + 2 * log10(DZ[i][4]) + 50.077910954466375 - (2 - args->Gamma) * log10(1 + z);

        /* Check for the flux, log_L_X = log_lambda_SAR + log_M_star */
        if (l < l_lo)
            continue;
        if (l >= l_hi)
            continue;

        Nsum += pow10(get_log_Psi(x[0], x[1], z, args->args)) * dV;
    }

    /* Note: to deg-2 */
    return Nsum / 41252.96;
}

double get_N(double x[], size_t dim, void *p)
{
    /* Return Eq. 3 from Bon+16 */
    (void) (dim);
    struct args_Psi *args = (struct args_Psi *) p;

    /* Find out the zbin for the lookup table */
    size_t zbin = (size_t) (x[3] / 0.001);

    /* TODO: NH correction */
    //double log_f = x[1] + x[2] - 0.5 * (DZ[zbin][2] + DZ[zbin][3]);
    // NOTE: Magic number is log10(4 pi Mpc2_tp_cm2)
    double log_f = x[1] + x[2] - 2 * log10(DZ[zbin][4]) - 50.077910954466375 + laf05_get_corr(x[0], x[4]);

    double Psi = pow10(get_log_Psi(x[2], x[1], x[3], *args));
    double I = cap09_get_area(log_f) / cap09_get_area(99.0);
    double f = ued14_get_f(x[1] + x[2], x[3], x[0]);
    double dV_div_dz = (DZ[zbin][6] - DZ[zbin][5]) / 0.001;

    return Psi * I * f * dV_div_dz;
}
