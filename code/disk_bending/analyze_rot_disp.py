#!/usr/bin/env python3
"""Fit the saved rigid-channel time series from test_rot_disp.

For each output/rotdisp_m<mode>_N<Nx>_seed<seed>.dat, fit

    a0(t) = A exp(-gamma t) cos(omega t + phi) + C

The damped-cosine fit is the authoritative frequency for damped modes: the
in-run periodogram is a peak estimate, and before mean subtraction it was
outright wrong for strongly damped modes (mode 4 read 0.06 instead of ~39).
The fitted gamma supersedes the crude first/last-fifth envelope estimate,
which is biased low by the offset C.

Prints one table row per file plus the razor-thin comparison columns.
Usage: python3 analyze_rot_disp.py [output_dir]
"""
import glob
import re
import sys

import numpy as np
from scipy.optimize import curve_fit

G = 4.300917270e-6           # kpc (km/s)^2 / Msun
SIGMA = 50.0 * 1.0e6         # Msun / kpc^2
LX = 16.0
SIGMA_X_REALISED = {1024: 38.05, 2048: None}  # None: print k^2 sx^2 with 40


def model(t, A, g, w, ph, C):
    return A * np.exp(-g * t) * np.cos(w * t + ph) + C


def fit_file(path):
    d = np.loadtxt(path)
    t, y = d[:, 0], d[:, 1]
    ym = y - y.mean()
    ws = np.linspace(0.5, 96, 3000)
    P = np.abs(np.exp(1j * np.outer(ws, t)) @ ym) ** 2
    w0 = ws[np.argmax(P)]
    p, cov = curve_fit(model, t, y, p0=[np.max(np.abs(ym)), 0.5, w0, 0.0,
                                        y.mean()], maxfev=40000)
    err = np.sqrt(np.diag(cov))
    if p[2] < 0:
        p[2] = -p[2]
    return w0, p, err


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else "output"
    files = sorted(glob.glob(f"{outdir}/rotdisp_m*_N*_seed*.dat"))
    if not files:
        print(f"no rotdisp_*.dat files under {outdir}/")
        return 1
    print(f"{'file':44s} {'k':>6s} {'omega_fit':>12s} {'gamma_fit':>12s} "
          f"{'w2/2piGSk':>10s}")
    for f in files:
        m = re.search(r"rotdisp_m(\d+)_N(\d+)_seed(\d+)\.dat", f)
        mode, nx = int(m.group(1)), int(m.group(2))
        k = 2.0 * np.pi * mode / LX
        try:
            w0, p, err = fit_file(f)
        except Exception as e:
            print(f"{f:44s} fit failed: {e}")
            continue
        wg2 = 2.0 * np.pi * G * SIGMA * k
        print(f"{f:44s} {k:6.3f} {p[2]:7.2f}+-{err[2]:4.2f} "
              f"{p[1]:7.2f}+-{err[1]:4.2f} {p[2]**2 / wg2:10.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
