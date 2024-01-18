#include <stdio.h>
#include "cap09.h"

int main(void)
{
    printf("%f\n", cap09_get_area(-10.0));
    printf("%f\n", cap09_get_area(-20.0));
    printf("%f\n", cap09_get_area(-15.0));
    printf("%f\n", cap09_get_area(-14.0));
    return 0;
}
