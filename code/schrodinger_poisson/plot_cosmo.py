#!/usr/bin/env python3
"""Plot density slices from cosmological S-P run."""

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import re


def plot_cosmo_slices(output_dir, n_show=6):
    files = sorted(glob.glob(os.path.join(output_dir, "slice_*.dat")))
    if not files:
        print(f"No files in {output_dir}/")
        return

    indices = np.linspace(0, len(files) - 1, min(n_show, len(files)), dtype=int)
    selected = [files[i] for i in indices]

    n = len(selected)
    fig, axes = plt.subplots(1, n, figsize=(3.5 * n, 3.2), squeeze=False)

    for i, fpath in enumerate(selected):
        with open(fpath) as f:
            header = f.readline()
        match = re.search(r'a = ([\d.e+-]+), z = ([\d.e+-]+)', header)
        a_val = float(match.group(1)) if match else 0.0
        z_val = float(match.group(2)) if match else 0.0

        data = np.loadtxt(fpath)
        N = int(np.sqrt(len(data)))
        # Use log(2+delta) for better contrast (thesis convention)
        log2rho = data[:, 3].reshape(N, N)

        ax = axes[0][i]
        pcm = ax.imshow(log2rho.T, origin='lower', extent=[0, 1, 0, 1],
                        cmap='inferno', aspect='equal')
        ax.set_title(f"z = {z_val:.1f}  (a = {a_val:.3f})", fontsize=10)
        ax.set_xlabel("x (box)")
        if i == 0:
            ax.set_ylabel("y (box)")
        plt.colorbar(pcm, ax=ax, label=r"$\log(2+\delta)$", shrink=0.8)

    plt.suptitle(f"S-P Cosmological Simulation ({os.path.basename(output_dir)})",
                 fontsize=12)
    plt.tight_layout()
    outfile = os.path.join(output_dir, "cosmo_slices.png")
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    print(f"Saved {outfile}")


if __name__ == "__main__":
    for d in sorted(glob.glob("output/sp_cosmo_*")):
        if os.path.isdir(d):
            plot_cosmo_slices(d)
    plt.show()
