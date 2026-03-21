#!/usr/bin/env python3
"""Plot tophat collapse: radial density profiles and z-slices."""

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import re


def plot_radial_profiles(output_dir):
    """Overlay radial density profiles at different times."""
    files = sorted(glob.glob(os.path.join(output_dir, "radial_*.dat")))
    if not files:
        print(f"No radial files in {output_dir}/")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    cmap = plt.cm.viridis
    n = len(files)

    for i, fpath in enumerate(files):
        with open(fpath) as f:
            header = f.readline()
        match = re.search(r'step = (\d+), time = ([\d.e+-]+)', header)
        step = int(match.group(1)) if match else 0
        time = float(match.group(2)) if match else 0.0

        data = np.loadtxt(fpath)
        r = data[:, 0]
        rho = data[:, 1]

        color = cmap(i / max(n - 1, 1))
        alpha = 0.3 + 0.7 * (i / max(n - 1, 1))

        # Only label a subset to keep legend readable
        label = f"t={time:.4f}" if i % max(n // 8, 1) == 0 or i == n - 1 else None

        ax1.plot(r, rho, color=color, alpha=alpha, linewidth=1.2, label=label)
        ax2.plot(r, rho, color=color, alpha=alpha, linewidth=1.2, label=label)

    ax1.set_xlabel("r (box units)")
    ax1.set_ylabel(r"$\langle|\psi|^2\rangle$ (spherical avg)")
    ax1.set_title("Radial density profile")
    ax1.legend(fontsize=8, loc='upper right')
    ax1.axvline(x=0.2, color='gray', linestyle='--', alpha=0.5, label='R_tophat')

    # Zoomed view near centre
    ax2.set_xlabel("r (box units)")
    ax2.set_ylabel(r"$\langle|\psi|^2\rangle$")
    ax2.set_title("Zoomed: central region")
    ax2.set_xlim(0, 0.3)
    ax2.axvline(x=0.2, color='gray', linestyle='--', alpha=0.5)

    plt.tight_layout()
    outfile = os.path.join(output_dir, "radial_profiles.png")
    plt.savefig(outfile, dpi=150)
    print(f"Saved {outfile}")


def plot_slices(output_dir, n_show=6):
    """Show selected z-slices."""
    files = sorted(glob.glob(os.path.join(output_dir, "slice_*.dat")))
    if not files:
        return

    # Select evenly spaced subset
    indices = np.linspace(0, len(files) - 1, min(n_show, len(files)), dtype=int)
    selected = [files[i] for i in indices]

    n = len(selected)
    fig, axes = plt.subplots(1, n, figsize=(3.5 * n, 3), squeeze=False)

    for i, fpath in enumerate(selected):
        with open(fpath) as f:
            header = f.readline()
        match = re.search(r'step = (\d+), time = ([\d.e+-]+)', header)
        step = int(match.group(1)) if match else 0
        time = float(match.group(2)) if match else 0.0

        data = np.loadtxt(fpath)
        N = int(np.sqrt(len(data)))
        rho = data[:, 2].reshape(N, N)

        ax = axes[0][i]
        pcm = ax.imshow(rho.T, origin='lower', extent=[0, 1, 0, 1],
                        cmap='inferno', aspect='equal')
        ax.set_title(f"t={time:.4f}", fontsize=10)
        ax.set_xlabel("x")
        if i == 0:
            ax.set_ylabel("y")
        plt.colorbar(pcm, ax=ax, shrink=0.8)

    plt.suptitle("Tophat collapse: z-slices", fontsize=13)
    plt.tight_layout()
    outfile = os.path.join(output_dir, "slices.png")
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    print(f"Saved {outfile}")


if __name__ == "__main__":
    d = "output/sp_3d_tophat_grav"
    if os.path.isdir(d):
        plot_radial_profiles(d)
        plot_slices(d)
    else:
        print(f"Directory {d} not found")
    plt.show()
