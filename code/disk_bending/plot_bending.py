#!/usr/bin/env python3
"""Plot the 2D bending-wave output written by sheet_2d.

    python3 plot_bending.py output/bend_2d

Panels:
  1. the corrugation <z>(x) at a few times
  2. the mode amplitude Re<z>_k(t) with its envelope, giving omega and the
     damping rate
  3. the <z> to <v_z> phase offset, the plan's §6 discriminator. Zero or +/-pi
     is a standing wave; +/- pi/2 is a propagating one. A displacement initial
     condition produces a standing wave by construction, so seeing 0 there is
     a check that the diagnostic works, not a null result.
"""
import sys
import os
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_info(outdir):
    info = {}
    p = os.path.join(outdir, "run_info.dat")
    if os.path.exists(p):
        for line in open(p):
            q = line.split()
            if len(q) == 2:
                try:
                    info[q[0]] = float(q[1])
                except ValueError:
                    info[q[0]] = q[1]
    return info


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else "output/bend_2d"
    info = read_info(outdir)
    d = np.loadtxt(os.path.join(outdir, "modes.dat"))
    t, zk_re, vk_re = d[:, 1], d[:, 2], d[:, 4]
    zabs, offset = d[:, 6], d[:, 8]

    fig, axes = plt.subplots(3, 1, figsize=(9, 9.5))

    # --- corrugation profiles ---
    files = sorted(glob.glob(os.path.join(outdir, "columns_*.dat")))
    pick = files[:: max(1, len(files) // 5)][:5]
    for f in pick:
        c = np.loadtxt(f)
        tm = float(open(f).readline().split("(")[1].split()[0])
        axes[0].plot(c[:, 0], c[:, 2], lw=1.2, label=f"{tm:.0f} Myr")
    axes[0].set_xlabel("x [kpc]")
    axes[0].set_ylabel(r"$\langle z\rangle$ [kpc]")
    axes[0].set_title("Corrugation profile", fontsize=10)
    axes[0].legend(fontsize=7, ncol=5)
    axes[0].grid(alpha=0.3)

    # --- mode amplitude and frequency ---
    axes[1].plot(t, zk_re, lw=1.0, label=r"Re$\langle z\rangle_k$")
    axes[1].plot(t, zabs, lw=1.2, color="k", alpha=0.6, label="envelope")
    axes[1].plot(t, -zabs, lw=1.2, color="k", alpha=0.6)
    # zero crossings give the measured frequency
    sign = np.sign(zk_re)
    cross = np.where(np.diff(sign) != 0)[0]
    if len(cross) > 2:
        tc = t[cross] + (t[cross + 1] - t[cross]) * (
            zk_re[cross] / (zk_re[cross] - zk_re[cross + 1]))
        half = np.mean(np.diff(tc)) / 977.79222  # Myr -> code
        w_meas = np.pi / half
        w_thin = info.get("omega_thin", 0.0)
        axes[1].set_title(
            rf"Mode {int(info.get('mode', 1))}: "
            rf"$\omega_{{\rm meas}}$ = {w_meas:.2f}, "
            rf"$\omega_{{\rm thin}}$ = {w_thin:.2f} km/s/kpc "
            rf"(ratio {w_meas/w_thin:.3f}), $kh$ = "
            rf"{info.get('k',0)*info.get('h_measured',1):.3f}", fontsize=10)
    axes[1].set_xlabel("t [Myr]")
    axes[1].set_ylabel(r"$\langle z\rangle_k$ [kpc]")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3)

    # --- phase offset ---
    axes[2].plot(t, offset / np.pi, lw=1.0)
    for y, lab in [(0.0, "standing"), (0.5, r"propagating (+$\pi/2$)"),
                   (-0.5, r"propagating ($-\pi/2$)")]:
        axes[2].axhline(y, ls="--", lw=0.8, color="grey")
        axes[2].text(t[-1], y, " " + lab, fontsize=7, va="center")
    axes[2].set_ylim(-1.05, 1.05)
    axes[2].set_xlabel("t [Myr]")
    axes[2].set_ylabel(r"$\arg\langle v_z\rangle_k - \arg\langle z\rangle_k$  [$\pi$]")
    axes[2].set_title(r"§6 discriminator: $\langle z\rangle$ to $\langle v_z\rangle$ "
                      "phase offset", fontsize=10)
    axes[2].grid(alpha=0.3)

    fig.tight_layout()
    out = os.path.join(outdir, "bending.png")
    fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
