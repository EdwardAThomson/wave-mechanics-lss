#!/usr/bin/env python3
"""Plot z-slice output from sp_3d."""

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import re


def plot_slices(output_dir):
    files = sorted(glob.glob(os.path.join(output_dir, "slice_*.dat")))
    if not files:
        print(f"No slice files in {output_dir}/")
        return

    n = len(files)
    n_cols = min(4, n)
    n_rows = (n + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 3.5 * n_rows),
                             squeeze=False)

    for i, fpath in enumerate(files):
        with open(fpath) as f:
            header = f.readline()
        match = re.search(r'step = (\d+), time = ([\d.e+-]+)', header)
        step = int(match.group(1)) if match else 0
        time = float(match.group(2)) if match else 0.0

        data = np.loadtxt(fpath)
        x = data[:, 0]
        y = data[:, 1]
        rho = data[:, 2]

        N = int(np.sqrt(len(x)))
        X = x.reshape(N, N)
        Y = y.reshape(N, N)
        R = rho.reshape(N, N)

        ax = axes[i // n_cols][i % n_cols]
        pcm = ax.pcolormesh(X, Y, R, shading='auto', cmap='inferno')
        ax.set_title(f"t = {time:.4f} (step {step})", fontsize=10)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_aspect('equal')
        plt.colorbar(pcm, ax=ax, label=r"$|\psi|^2$", shrink=0.8)

    for i in range(n, n_rows * n_cols):
        axes[i // n_cols][i % n_cols].set_visible(False)

    plt.tight_layout()
    outfile = os.path.join(output_dir, "slices.png")
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    print(f"Saved {outfile}")


if __name__ == "__main__":
    for d in ["output/sp_3d", "output/sp_3d_grav"]:
        if os.path.isdir(d):
            plot_slices(d)
    plt.show()
