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
With --plot, also renders figures/dispersion_relation.png: omega(k) against
the candidate relations on top, the measured Landau damping gamma(k) below.
Usage: python3 analyze_rot_disp.py [output_dir] [--plot]
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


H_SHEET = 0.2951        # kpc, measured scale height of the Sigma=50 sheet
OMEGA0 = 47.84          # km/s/kpc, midplane vertical frequency
SX = 38.05              # km/s, realised sigma_x at Nx = 1024

# Reference palette (dataviz skill): first three categorical slots, valid
# all-pairs; ink/chrome tokens for a light single-mode print figure.
C_BLUE, C_ORANGE, C_AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED, GRID, SURFACE = "#0b0b0b", "#898781", "#e1e0d9", "#fcfcfb"


def make_figure(rows, path="figures/dispersion_relation.png"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    kk = np.linspace(0.15, 2.15, 400)
    wg = np.sqrt(2.0 * np.pi * G * SIGMA * kk)
    wp2 = wg ** 2 - kk ** 2 * SX ** 2
    wthick = np.sqrt(wg ** 2 / (1.0 + kk * H_SHEET))

    fig, (ax, axg) = plt.subplots(
        2, 1, figsize=(6.0, 6.4), sharex=True, dpi=200,
        gridspec_kw={"height_ratios": [2.4, 1.0], "hspace": 0.08})
    for a in (ax, axg):
        a.set_facecolor(SURFACE)
        a.grid(True, color=GRID, linewidth=0.7)
        for s in a.spines.values():
            s.set_color(MUTED)
        a.tick_params(colors=MUTED, labelcolor=INK)
    fig.patch.set_facecolor(SURFACE)

    ax.plot(kk, wg, color=C_BLUE, lw=1.8, label="razor-thin, gravity only")
    ax.plot(kk, wthick, color=C_AQUA, lw=1.8,
            label="thickness-corrected  1/(1+kh)")
    kp = kk[wp2 > 0]
    ax.plot(kp, np.sqrt(wp2[wp2 > 0]), color=C_ORANGE, lw=1.8,
            label="with free-streaming pressure")
    ax.axhline(OMEGA0, color=MUTED, lw=1.0, ls=(0, (4, 3)))
    ax.text(0.17, OMEGA0 + 0.8, r"$\omega_0$ (vertical)", color=MUTED,
            fontsize=8)

    n1024 = [r for r in rows if r["nx"] == 1024]
    n2048 = [r for r in rows if r["nx"] != 1024]
    ax.errorbar([r["k"] for r in n1024], [r["w"] for r in n1024],
                yerr=[r["we"] for r in n1024], fmt="o", color=INK, ms=6,
                capsize=3, lw=1.2, zorder=5, label="measured (Nx = 1024)")
    if n2048:
        ax.errorbar([r["k"] for r in n2048], [r["w"] for r in n2048],
                    yerr=[r["we"] for r in n2048], fmt="o", mfc="none",
                    color=INK, ms=7, capsize=3, lw=1.2, zorder=5,
                    label="convergence check (Nx = 2048)")
    ax.set_ylabel(r"$\omega$  [km s$^{-1}$ kpc$^{-1}$]", color=INK)
    ax.set_ylim(0, 58)
    ax.legend(loc="lower right", fontsize=8, framealpha=0.95,
              edgecolor=GRID)

    axg.errorbar([r["k"] for r in n1024], [r["g"] for r in n1024],
                 yerr=[r["ge"] for r in n1024], fmt="o", color=INK, ms=6,
                 capsize=3, lw=1.2)
    axg.set_ylabel(r"$\gamma$  [unit$^{-1}$]", color=INK)
    axg.set_xlabel(r"$k$  [kpc$^{-1}$]", color=INK)
    axg.set_ylim(bottom=-0.5)

    ax.set_title("Stable bending branch, rotating box "
                 r"($Q = 2.45$, $\Sigma = 50\,M_\odot\,$pc$^{-2}$)",
                 color=INK, fontsize=10)
    fig.savefig(path, bbox_inches="tight")
    print(f"wrote {path}")


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    outdir = args[0] if args else "output"
    files = sorted(glob.glob(f"{outdir}/rotdisp_m*_N*_seed*.dat"))
    if not files:
        print(f"no rotdisp_*.dat files under {outdir}/")
        return 1
    rows = []
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
        rows.append({"mode": mode, "nx": nx, "k": k, "w": p[2],
                     "we": err[2], "g": abs(p[1]), "ge": err[1]})
    if "--plot" in sys.argv and rows:
        make_figure(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
