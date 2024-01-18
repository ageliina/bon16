#!/usr/bin/env python3
# encoding: utf-8
# Author: Akke Viitanen
# Email: akke.viitanen@helsinki.fi
# Date: 2023-06-02 14:45:40

"""
Create dz lookup table
"""

from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
import numpy as np


def main():

    # Helper function to write to file
    def my_print(*args, **kwargs):
        print(*args, **kwargs)
        print(*args, **kwargs, file=f)

    # Setup the cosmology and distances
    cosmo70 = FlatLambdaCDM(Om0=0.30, H0=70, Tcmb0=2.73)
    nz = 10000
    zbins = np.linspace(0.0, 10.0, nz + 1)

    # Write to file with C-formatting
    with open("dz.h", 'w') as f:
        my_print("#define NZ %dL" % nz)
        my_print("double DZ[%d][7] = {" % nz)
        my_print("/* Om0=0.30 H0=70 Tcmb0=2.73 */")
        my_print("/* All lengths in Mpc */")
        my_print("/* z1 z2 dc12 da12 dl12 dV1 dV2 */")
        for z1, z2 in zip(zbins[:-1], zbins[1:]):
            delim = ',' if z2 != zbins[-1] else ''
            my_print("    {%.10e, %.10e, %.10e, %.10e, %.10e, %.10e, %.10e}%s"
                     % (
                         z1,
                         z2,
                         cosmo70.comoving_distance(0.5 * (z1 + z2)).value,
                         cosmo70.angular_diameter_distance(0.5 * (z1 + z2)).value,
                         cosmo70.luminosity_distance(0.5 * (z1 + z2)).value,
                         cosmo70.comoving_volume(z1).value,
                         cosmo70.comoving_volume(z2).value,
                         delim
                     )
        )
        my_print("};")


if __name__ == "__main__":
    main()
