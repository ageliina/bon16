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


def convert_flux(
    S1,
    E1_min=2,
    E1_max=10,
    E2_min=2,
    E2_max=7,
    Gamma=1.9   # NOTE: aird+2021 below eq. 2
):
    """
    Convert flux S from bandpass E1 to bandpass E2 assuming a power-law
    spectrum with photon index Gamma
    """
    idx = (2 - Gamma)
    return S1 * np.true_divide(
        E2_max ** idx - E2_min ** idx,
        E1_max ** idx - E1_min ** idx
    )


def main():

    # Helper function to write to file
    def my_print(*args, **kwargs):
        print(*args, **kwargs)
        print(*args, **kwargs, file=f)

    luo2017 = np.loadtxt("../../data/luo2017/Default Dataset.csv")
    cap2009 = np.loadtxt("../../data/cappelluti2009/Default Dataset.csv")

    S = []
    N = []
    B = []
    for s, n in luo2017:
        s = convert_flux(s, 2.0, 7.0, 2.0, 10.0, 1.4)
        S.append(s)
        N.append(n)
        B.append(0)

    for s, n in cap2009:
        n /= (s / 1e-14) ** 1.5
        S.append(s)
        N.append(n)
        B.append(1)

    idx = np.argsort(S)
    S = np.array(S)[idx]
    N = np.array(N)[idx]
    B = np.array(B)[idx]

    # Write to file with C-formatting
    with open("logn_logs.h", 'w') as f:
        my_print("#define N_LOGN_LOGS %d" % len(S))
        my_print("double LOGN_LOGS[%d][3] = {" % len(S))
        my_print("/* combined logN-logS from Luo+2017 and Cappelluti+2009 */")
        my_print("/* S_2_10[erg/s/cm2] N(>S) */")
        for s, n, b in zip(S, N, B):
            delim = ',' if s != S[-1] else ''
            my_print("    {%e, %e, %d}%s" % (s, n, b, delim))
        my_print("};")


if __name__ == "__main__":
    main()
