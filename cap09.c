#include <stdio.h>
#include <stdlib.h>
#include "cap09.h"

double cap09_logf_area[CAP09_N][4] = {
    {-12.8, 2.13},
    {-12.9, 2.13},
    {-13.0, 2.13},
    {-13.1, 2.13},
    {-13.2, 2.13},
    {-13.3, 2.13},
    {-13.4, 2.13},
    {-13.5, 2.13},
    {-13.6, 2.13},
    {-13.7, 2.13},
    {-13.8, 2.13},
    {-13.9, 2.09},
    {-14.0, 1.98},
    {-14.1, 1.77},
    {-14.2, 1.33},
    {-14.3, 0.67},
    {-14.4, 0.14},
    {-14.5, 0.01},
};

double cap09_get_area(double log_f)
{

    /* -inf < log_f < log_f_min */
    if (log_f < cap09_logf_area[CAP09_N - 1][0])
        return 0;

    /* log_f_max < log_f < inf */
    if (log_f >= cap09_logf_area[0][0])
        return cap09_logf_area[0][1];

    /* log_f_min < log_f < log_f_max */
    size_t idx = CAP09_N - 1 - ((size_t) ((log_f + 14.5) / 0.10));
    double x0 = cap09_logf_area[idx+1][0];
    double y0 = cap09_logf_area[idx+1][1];
    double x1 = cap09_logf_area[idx+0][0];
    double y1 = cap09_logf_area[idx+0][1];
    return (y1 - y0) / (x1 - x0) * (log_f - x0) + y0;
}
