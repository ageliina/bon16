#include "ued14.h"

double ued14_get_Psi_4375(double z)
{
    return PSI_4375_0 * pow(1 + min(z, 2.0), A1);
}

double ued14_get_Psi(double log_L_X, double z)
{
    double ret = ued14_get_Psi_4375(z) - BETA * (log_L_X - 43.75);
    return minmax(PSI_MIN, ret, PSI_MAX);
}

double ued14_get_f(double log_L_X, double z, double log_N_H)
{
    double Psi = ued14_get_Psi(log_L_X, z);
    double piecewise[2][5] =
    {
        {
            1.0 - (2.0 + EPSILON) / (1.0 + EPSILON) * Psi,
            1.0 / (1.0 + EPSILON) * Psi,
            1.0 / (1.0 + EPSILON) * Psi,
            EPSILON / (1.0 + EPSILON) * Psi,
            F_CTK / 2.0 * Psi,
        },
        {
            2.0 / 3.0 - (3.0 + 2.0 * EPSILON) / (3.0 + 3.0 * EPSILON) * Psi,
            1.0 / 3.0 - EPSILON / (3.0 + 3.0 * EPSILON) * Psi,
            1.0 / (1.0 + EPSILON) * Psi,
            EPSILON / (1.0 + EPSILON) * Psi,
            F_CTK / 2.0 * Psi,
        }
    };

    int idx1 = (int) Psi < (1.0 + EPSILON) / (3.0 + EPSILON);
    int idx2 = (int) log_N_H - 20.0;
    return piecewise[minmax(0, idx1, 1)][minmax(0, idx2, 4)];

}
