#!/usr/bin/env python3
# encoding: utf-8
# Author: Akke Viitanen
# Email: akke.viitanen@helsinki.fi
# Date: 2023-05-11 18:18:25

"""
Fit Bon+16 to Miyaji+15 XLF and the stellar mass function
"""

import argparse
import glob
import math
import os
import random
import re
import subprocess
import sys
import time

import astropy as ap
import astropy.coordinates as c
import astropy.units as u
import fitsio
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import emcee

from smf import StellarMassFunctionCOSMOS2020, StellarMassFunctionBongiorno2016


MIY15 = np.loadtxt("data/miyaji2015/xlf_miyaji2015_table4.dat")
ZS = MIY15[:, 2]
SELECT = ZS >= 2.5
LXS = MIY15[SELECT, 5]
PHIS = MIY15[SELECT, -2]
DPHIS = MIY15[SELECT, -1]

MVEC = np.logspace(9.55, 12.45, 30)
SMF_GAL = StellarMassFunctionCOSMOS2020()
ZBINS = [z for z in SMF_GAL.get_zbins("total") if z[0] >= 2.5]
SMF_GAL = np.concatenate([
    SMF_GAL.get_stellar_mass_function(MVEC, np.mean(zbin), "all_best_fit")
    for zbin in ZBINS
])

#PARAMS = {
#    # name                  lower initial upper sigma
#    "log_Psi_star":         [-100, -6.86, 100, 1e-6],
#    "log_M_star_star":      [   0, 10.99, 100, 1e-6],
#    "k_log_M_star_star":    [-100, -0.25,   0, 1e-3],
#    "alpha":                [-100,  0.24, 100, 1e-6],
#    "log_lambda_SAR_star0": [   0, +33.8, 100, 1e-6],
#    "k_lambda":             [-100, -0.48, 100, 1e-6],
#    "log_M_star0":          [   0, 11.00, 100, 1e-6],
#    "gamma10":              [-100, -1.01, 100, 1e-6],
#    "k_gamma":              [-100, +0.58, 100, 1e-6],
#    "k_gamma1":             [-100, -0.06, 100, 1e-3],
#    "gamma2":               [-100, -3.72, 100, 1e-6],
#    "p1":                   [-100, +5.82, 100, 1e-6],
#    "p2":                   [-100, +2.36, 100, 1e-6],
#    "p3":                   [-100, -4.64, 100, 1e-3],
#    "z0":                   [   0, +1.10, 100, 1e-6],
#    "z1":                   [ 2.5, +2.59, 100, 1e-3],
##    "log_f":                [-100,  -2.0, 100, 0.01],   # fractional intrinsic scatter
#}

PARAMS = {
    # name                  lower initial upper sigma
    "log_Psi_star":         [-100, -6.86, 100, 1e-6],
    "log_M_star_star":      [   0, 10.99, 100, 1e-6],
    "k_log_M_star_star":    [-100, -0.25,   0, 1e-3],
    "alpha":                [-100,  0.24, 100, 1e-6],
    "log_lambda_SAR_star0": [   0, +33.8, 100, 1e-6],
    "k_lambda":             [-100, -0.48, 100, 1e-6],
    "log_M_star0":          [   0, 11.00, 100, 1e-6],
    "gamma10":              [-100, -1.01, 100, 1e-6],
    "k_gamma":              [-100, +0.58, 100, 1e-6],
    "k_gamma1":             [-100, +2.05, 100, 1e-3],
    "gamma2":               [-100, -3.72, 100, 1e-6],
    "p1":                   [-100, +5.82, 100, 1e-6],
    "p2":                   [-100, +2.36, 100, 1e-6],
    "p3":                   [-100, -1.92, 100, 1e-3],
    "z0":                   [   0, +1.10, 100, 1e-6],
    "z1":                   [ 2.5, +2.56, 100, 1e-3],
#    "log_f":                [-100,  -2.0, 100, 0.01],   # fractional intrinsic scatter
}



def log_prior(params):

    for k, (v_min, v, v_max, sigma) in params.items():
        if not (v_min <= v < v_max):
            return -np.inf

    return 0


from subprocess import Popen, PIPE, STDOUT

def log_likelihood(params):

    # Run Bongiorno+2016 external code to get Phi_star and Phi_X
    inp = (' '.join(["7"] + ["%f" % v[1] for v in params.values()])).encode()
    p = subprocess.run(
        ["src/bon16/main"],
        input=inp,
        capture_output=True,
    )

    ret = p.stdout[:-1].decode().split('\n')

    z_star = np.array([float(f.split()[0])    for f in ret[:360]])
    Phi_star = np.array([float(f.split()[-3]) for f in ret[:360]])
    Phi_star = Phi_star[z_star >= 2.5]

    z_X = np.array([float(f.split()[0])    for f in ret[360:-46]])
    Phi_X = np.array([float(f.split()[-3]) for f in ret[360:-46]])
    Phi_X = Phi_X[z_X >= 2.5]

    fx = np.array([float(f.split()[0]) for f in ret[-46:]])
    N_above_S = np.array([float(f.split()[-3]) for f in ret[-46:]])
    N_ABOVE_S = np.array([float(f.split()[-1]) for f in ret[-46:]])
    # De-select XMM-COSMOS limiting flux from Cap+09 Fig. 7
    N_above_S = N_above_S[fx < 3e-15]
    N_ABOVE_S = N_ABOVE_S[fx < 3e-15]

    #assert len(Phi_star) == 360, len(Phi_star)
    #assert len(Phi_X) == 196,    len(Phi_X)
    #assert len(N_above_S) == 46, len(N_above_S)

    plt.figure()
    plt.plot(SMF_GAL)
    plt.plot(Phi_star)

    plt.figure()
    plt.plot(PHIS)
    plt.plot(Phi_X)

    plt.figure()
    plt.plot(N_above_S)
    plt.plot(N_ABOVE_S)
    plt.show()

    if not np.all(np.isfinite(Phi_star)) or np.any(Phi_star > SMF_GAL):
        return -np.inf

    if not np.all(np.isfinite(Phi_X)):
        return -np.inf

    #sigma2 = DPHIS ** 2 + Phi_X ** 2 * np.exp(2 * params["log_f"][1])
    #assert np.all(sigma2 > 0)

    #return -0.5 * np.sum((Phi_X - PHIS) ** 2 / sigma2 + np.log(sigma2))
    #return -0.5 * np.sum((PHIS - Phi_X) ** 2 / sigma2 + np.log(2 * np.pi * sigma2))
    chisq_X = -0.5 * np.sum((Phi_X - PHIS) ** 2 / DPHIS ** 2)
    chisq_N = -0.5 * np.sum((N_above_S - N_ABOVE_S) ** 2 / (0.10 * N_ABOVE_S) ** 2)
    return chisq_X + chisq_N


def log_probability(theta, keys):

    # Set the parameters
    params = PARAMS
    for k, v in zip(keys, theta):
        params[k][1] = v

    lp = log_prior(params)
    if not np.isfinite(lp):
        return -np.inf
    ll = log_likelihood(params)

    #if np.isfinite(ll):
    chisq = -2 * ll / (len(PHIS) - 1)
    print(" ".join("%+6.4f" % t for t in theta), lp, "%7.2f" % ll, "%7.2f" % chisq)

    return lp + ll


def main():
    # Free parameters to fit for
    keys = "k_log_M_star_star", "k_gamma1", "p3", "z1"

    np.random.seed(20230529)
    ndim = len(keys)
    nwalkers = 2 * ndim

    pos = np.zeros((ndim, nwalkers))
    pos += np.array([PARAMS[k][1] for k in keys])[:, None]
    pos += np.random.randn(ndim, nwalkers) * np.array([PARAMS[k][-1] for k in keys])[:, None] * 10

    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=[keys])
    sampler.run_mcmc(pos.T, int(sys.argv[1]), progress=False)

    chain = sampler.get_chain()
    for i in range(nwalkers):
        np.savetxt("chain%d.txt" % i, chain[:, i, :])


if __name__ == "__main__":
    main()
