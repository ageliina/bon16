#include <stdlib.h>
#include <stdio.h>
#include "laf05.h"

int main(void)
{

    laf05_init();
    size_t i=0, j=0;

    for (double z=0.000; z<5.501; z+=0.0005) {
        for (double n=19.0; n<26.0; n+=0.10) {
            //printf("%lu %lu %f %f %e %e\n", i, j, z, n, laf05_get_corr(n, z), laf05_corr[i][j]);
            printf("%lu %lu %f %f %e\n", i, j, z, n, laf05_get_corr(n, z));
            j++;
        }
        i++; j = 0;
    }

    laf05_free();
	return 0;
}
