#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>

#include "bon16.h"
#include "laf05.h"
#include "miy15.h"
#include "ued14.h"
#include "wea23.h"

#include "logn_logs.h"

#define N_CALLS 2000
#define SEED 20230530L

#define LOG_M_STAR_MIN 8.0
#define LOG_M_STAR_MAX 13.5
#define LOG_LAMBDA_SAR_MIN 32.0
#define LOG_LAMBDA_SAR_MAX 36.0
#define LOG_N_H_MIN 20.0
#define LOG_N_H_MAX 24.0
#define Z_MIN 0.01
#define Z_MAX 6.99


void integrate(gsl_monte_function *G, double xlo[], double xhi[], double *res, double *err)
{
    /* Init */
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc(T);
    gsl_rng_env_setup();
    gsl_rng_set(r, SEED);
    gsl_monte_plain_state *s = gsl_monte_plain_alloc(G->dim);

    /* Integrate */
    gsl_monte_plain_integrate(G, xlo, xhi, G->dim, N_CALLS, r, s, res, err);

    /* De-init */
    gsl_monte_plain_free(s);
    gsl_rng_free(r);
}

void routine_Phi_star(struct args_Psi args)
{
    /* Phi_star routine */
    double res, err;
    for (size_t i=0; i<N_WEA23; i++) {
        struct args_Phi_star args_Phi_star = { wea23[i][1], wea23[i][0], args };
        gsl_monte_function G = { &get_Phi_star, 1, &args_Phi_star };
        integrate(
                &G,
                (double[1]) {LOG_LAMBDA_SAR_MIN},
                (double[1]) {LOG_LAMBDA_SAR_MAX},
                &res, &err
                );
        printf("%e %e %e %e %e\n", wea23[i][0], wea23[i][1], res, err, err / res);
    }
}

void routine_Phi_X(struct args_Psi args)
{
    /* Phi_X routine */
    double res, err;
    for (size_t i=0; i<N_MIY15; i++) {
        struct args_Phi_X args_Phi_X = {
            miy15[i][1],
            miy15[i][2],
            miy15[i][0],
            args
        };
        gsl_monte_function G = { &get_Phi_X, 2, &args_Phi_X };
        integrate(
                &G,
                (double[2]) {LOG_M_STAR_MIN, LOG_LAMBDA_SAR_MIN},
                (double[2]) {LOG_M_STAR_MAX, LOG_LAMBDA_SAR_MAX},
                &res,
                &err
                );
        printf("%e %e %e %e %e %e\n", miy15[i][0], miy15[i][1], miy15[i][2], res, err, err / res);
    }
}

void routine_N_above_S(struct args_Psi args)
{

    /* logN-logS routine */
    double N = 0;
    double res, err;
    double Gamma_air21 = 1.9;
    //double Gamma_luo17 = 1.4;

    for (int i=N_LOGN_LOGS-1; i>=0; i--) {

        /* Convert to 2-10 keV using Luo+17 power-law Gamma */
        double F_X_lo = LOGN_LOGS[i+0][0];
        double F_X_hi = (i+1 < N_LOGN_LOGS) ? LOGN_LOGS[i+1][0] : 1e-10;
        struct args_N_above_S args_N_above_S = {
            log10(F_X_lo),
            log10(F_X_hi),
            Gamma_air21,
            //Gamma_luo17,
            args
        };

        gsl_monte_function G = { &get_N_above_S, 2, &args_N_above_S };
        integrate(
                &G,
                (double[2]) {LOG_M_STAR_MIN, LOG_LAMBDA_SAR_MIN},
                (double[2]) {LOG_M_STAR_MAX, LOG_LAMBDA_SAR_MAX},
                &res,
                &err
                );
        N += res;
        printf("%e %e %e %e\n", F_X_lo, N, 0., LOGN_LOGS[i][1]);
    }
}

void routine_N(struct args_Psi args)
{
    double res, err;

    laf05_init();
    gsl_monte_function G = { &get_N, 4, &args };
    integrate(
            &G,
            (double[4]) {LOG_N_H_MIN, LOG_LAMBDA_SAR_MIN, LOG_M_STAR_MIN, Z_MIN},
            (double[4]) {LOG_N_H_MAX, LOG_LAMBDA_SAR_MAX, LOG_M_STAR_MAX, Z_MAX},
            &res,
            &err
            );
    laf05_free();
    printf("%e %e\n", res, err);

}

int main(void)
{
    /* Setup the variables for the SMF/XLF */
    double theta[N_PARAMETERS];
    char buf[1000];

    while (fgets(buf, 1000, stdin)) {

        /* Parse the input arguments */
        char *token = strtok(buf, " ");
        size_t flag;
        for (size_t i=0; i<N_PARAMETERS+1; i++, token=strtok(NULL, " ")) {
            if (i == 0)
                sscanf(token, "%lu", &flag);
            else
                sscanf(token, "%lf", &theta[i-1]);
        }

        /* NOTE: order of theta is the same as in Bon+16 table */
        struct args_Psi args = {
            /* log_Psi_star */
            theta[0],
            /* f_lambda_SAR */
            {
                theta[4],   /* log_lambda_star0 */
                theta[5],   /* k_lambda */
                theta[6],   /* log_M_star0 */
                theta[7],   /* gamma_10 */
                theta[8],   /* k_gamma */
                theta[9],   /* k_gamma1 */
                theta[10],  /* gamma2 */
                theta[14],  /* z0 */
                theta[15]   /* z1 */
            },
            /* f_star */
            {
                theta[1],   /* log_M_star_star */
                theta[2],   /* alpha */
                theta[3],   /* k_log_M_star_star0 */
                theta[15]   /* z1 */
            },
            /* f_z */
            {
                theta[11],  /* p1 */
                theta[12],  /* p2 */
                theta[13],  /* p3 */
                theta[14],  /* z0 */
                theta[15]   /* z1 */
            }
        };

        if (flag & 1) routine_Phi_star  (args);
        if (flag & 2) routine_Phi_X     (args);
        if (flag & 4) routine_N_above_S (args);
        if (flag & 8) routine_N         (args);
    }

    return 0;
}
