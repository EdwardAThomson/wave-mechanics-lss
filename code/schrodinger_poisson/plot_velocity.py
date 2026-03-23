#!/usr/bin/env python3
"""Compare velocity fields from probability current vs phase gradient methods."""

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import re


def plot_velocity_comparison(output_dir, n_show=4):
    files = sorted(glob.glob(os.path.join(output_dir, "vel_*.dat")))
    if not files:
        print(f"No velocity files in {output_dir}/")
        return

    indices = np.linspace(0, len(files) - 1, min(n_show, len(files)), dtype=int)
    selected = [files[i] for i in indices]

    n = len(selected)
    fig, axes = plt.subplots(3, n, figsize=(4 * n, 10), squeeze=False)

    for col, fpath in enumerate(selected):
        with open(fpath) as f:
            header = f.readline()
        match = re.search(r'a = ([\d.e+-]+), z = ([\d.e+-]+)', header)
        a_val = float(match.group(1)) if match else 0
        z_val = float(match.group(2)) if match else 0

        data = np.loadtxt(fpath)
        Ng = int(np.sqrt(len(data)))
        rho = data[:, 2].reshape(Ng, Ng)
        vx_curr = data[:, 3].reshape(Ng, Ng)
        vy_curr = data[:, 4].reshape(Ng, Ng)
        vx_phase = data[:, 5].reshape(Ng, Ng)
        vy_phase = data[:, 6].reshape(Ng, Ng)

        speed_curr = np.sqrt(vx_curr**2 + vy_curr**2)
        speed_phase = np.sqrt(vx_phase**2 + vy_phase**2)

        # Row 0: density
        ax = axes[0][col]
        pcm = ax.imshow(np.log(2 + rho - 1).T, origin='lower',
                        extent=[0, 1, 0, 1], cmap='inferno', aspect='equal')
        ax.set_title(f"z={z_val:.1f} (a={a_val:.3f})", fontsize=10)
        if col == 0:
            ax.set_ylabel("Density log(2+d)")

        # Row 1: current method speed
        ax = axes[1][col]
        vmax = np.percentile(speed_curr, 99)
        pcm = ax.imshow(speed_curr.T, origin='lower', extent=[0, 1, 0, 1],
                        cmap='viridis', aspect='equal', vmin=0, vmax=vmax)
        if col == 0:
            ax.set_ylabel("|v| current method")
        plt.colorbar(pcm, ax=ax, shrink=0.7)

        # Row 2: phase method speed
        ax = axes[2][col]
        pcm = ax.imshow(speed_phase.T, origin='lower', extent=[0, 1, 0, 1],
                        cmap='viridis', aspect='equal', vmin=0, vmax=vmax)
        if col == 0:
            ax.set_ylabel("|v| phase method")
        plt.colorbar(pcm, ax=ax, shrink=0.7)

    plt.suptitle("Velocity comparison: current vs phase gradient", fontsize=13)
    plt.tight_layout()
    outfile = os.path.join(output_dir, "velocity_comparison.png")
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    print(f"Saved {outfile}")

    # --- Scatter plot: current vs phase at final time ---
    fpath = files[-1]
    with open(fpath) as f:
        header = f.readline()
    match = re.search(r'a = ([\d.e+-]+), z = ([\d.e+-]+)', header)
    a_val = float(match.group(1)) if match else 0
    z_val = float(match.group(2)) if match else 0

    data = np.loadtxt(fpath)
    rho = data[:, 2]
    vx_curr = data[:, 3]
    vx_phase = data[:, 5]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

    # Scatter: vx current vs vx phase
    ax1.scatter(vx_curr, vx_phase, s=1, alpha=0.3, c=np.log10(rho + 1e-10), cmap='inferno')
    lim = max(abs(vx_curr.max()), abs(vx_curr.min()),
              abs(vx_phase.max()), abs(vx_phase.min())) * 1.1
    ax1.plot([-lim, lim], [-lim, lim], 'k--', linewidth=0.5, alpha=0.5)
    ax1.set_xlabel("vx (current method)")
    ax1.set_ylabel("vx (phase method)")
    ax1.set_title(f"vx comparison at z={z_val:.1f}")
    ax1.set_xlim(-lim, lim)
    ax1.set_ylim(-lim, lim)
    ax1.set_aspect('equal')

    # Correlation coefficient
    mask = rho > np.percentile(rho, 10)  # exclude very low density
    corr = np.corrcoef(vx_curr[mask], vx_phase[mask])[0, 1]
    ax1.text(0.05, 0.95, f"r = {corr:.4f}", transform=ax1.transAxes,
             fontsize=11, va='top')

    # Histogram of difference
    diff = vx_curr - vx_phase
    ax2.hist(diff, bins=100, density=True, alpha=0.7)
    ax2.set_xlabel("vx_current - vx_phase")
    ax2.set_ylabel("PDF")
    ax2.set_title(f"Velocity method difference (rms = {np.std(diff):.2e})")
    ax2.axvline(0, color='k', linestyle='--', linewidth=0.5)

    plt.tight_layout()
    outfile = os.path.join(output_dir, "velocity_scatter.png")
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    print(f"Saved {outfile}")


if __name__ == "__main__":
    for d in sorted(glob.glob("output/sp_cosmo_*")):
        if os.path.isdir(d):
            plot_velocity_comparison(d)
    plt.show()
