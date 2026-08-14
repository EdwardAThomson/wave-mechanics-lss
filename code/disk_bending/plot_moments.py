#!/usr/bin/env python3
"""Plot the vertical-moment history written by sheet_1d.

    python3 plot_moments.py output/kick_bulk [output/other ...]

Top panel: the corrugation <z> and the mean vertical velocity <v_z>/omega on
the same axis, so their quarter-cycle offset is visible directly. In 1D that
offset is the whole of the plan's §6 primary diagnostic; in 2D it becomes a
cross-correlation in x.

Middle: the oscillation envelope, which is where phase mixing shows up. An
undamped envelope means the perturbation is not winding, and in a harmonic
external potential it never will (Kohn's theorem).

Bottom: conservation and aliasing checks.
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_info(outdir):
    info = {}
    path = os.path.join(outdir, "run_info.dat")
    if os.path.exists(path):
        for line in open(path):
            p = line.split()
            if len(p) == 2:
                try:
                    info[p[0]] = float(p[1])
                except ValueError:
                    info[p[0]] = p[1]
    return info


def envelope(t, y):
    """Amplitude envelope from |analytic signal| via the Hilbert transform."""
    n = len(y)
    Y = np.fft.fft(y - y.mean())
    h = np.zeros(n)
    h[0] = 1
    if n % 2 == 0:
        h[n // 2] = 1
        h[1:n // 2] = 2
    else:
        h[1:(n + 1) // 2] = 2
    return np.abs(np.fft.ifft(Y * h))


def main():
    dirs = [a for a in sys.argv[1:] if not a.startswith("--")] or ["output/kick_bulk"]
    fig, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True)

    for outdir in dirs:
        d = np.loadtxt(os.path.join(outdir, "moments.dat"))
        info = read_info(outdir)
        omega = info.get("omega0", 75.0)
        label = os.path.basename(outdir.rstrip("/"))
        t = d[:, 1]
        zbar, vbar = d[:, 3], d[:, 4]
        energy, spill = d[:, 10], d[:, 11]

        # Only the first run goes in the top panel: overlaying several runs'
        # oscillations just produces a solid block of ink.
        if outdir == dirs[0]:
            axes[0].plot(t, zbar, lw=1.1,
                         label=rf"{label}: $\langle z\rangle$")
            axes[0].plot(t, vbar / omega, lw=1.1, alpha=0.8,
                         label=rf"{label}: $\langle v_z\rangle/\omega_z$")
        env = envelope(t, zbar)
        axes[1].plot(t, env, lw=1.3, label=label)
        axes[2].semilogy(t, np.maximum(np.abs(energy / energy[0] - 1), 1e-16),
                         lw=1.0, label=f"{label}: |dE/E|")
        axes[2].semilogy(t, np.maximum(spill, 1e-16), lw=1.0, ls="--",
                         label=f"{label}: spill")

    axes[0].set_ylabel("kpc")
    axes[0].legend(fontsize=8, ncol=2)
    axes[0].set_title("Corrugation and mean vertical velocity "
                      "(quarter-cycle offset = propagating oscillation)",
                      fontsize=10)
    axes[0].grid(alpha=0.3)

    axes[1].set_ylabel(r"envelope of $\langle z\rangle$ [kpc]")
    axes[1].set_title("Phase mixing: a decaying envelope is the spiral winding up. "
                      "Runs of equal\nperiod count end at different times because "
                      "$T_z$ differs. Ends are Hilbert edge artefacts.",
                      fontsize=9)
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3)

    axes[2].set_xlabel("t [Myr]")
    axes[2].set_ylabel("relative")
    axes[2].set_title("Conservation and velocity-space aliasing", fontsize=10)
    axes[2].legend(fontsize=7, ncol=2)
    axes[2].grid(alpha=0.3)

    fig.tight_layout()
    out = os.path.join(dirs[0], "moments.png")
    fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
