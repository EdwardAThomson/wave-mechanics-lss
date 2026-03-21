#!/usr/bin/env python3
"""Plot output from sp_1d: density evolution with and without gravity."""

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import re
import sys


def load_outputs(output_dir):
    """Load all psi_*.dat files from a directory."""
    files = sorted(glob.glob(os.path.join(output_dir, "psi_*.dat")))
    if not files:
        return [], [], []

    datasets = []
    for fpath in files:
        with open(fpath) as f:
            header = f.readline()
        match = re.search(r'step = (\d+), time = ([\d.e+-]+)', header)
        step = int(match.group(1)) if match else 0
        time = float(match.group(2)) if match else 0.0
        data = np.loadtxt(fpath)
        datasets.append((step, time, data))

    return files, datasets


def plot_single(output_dir, title_prefix=""):
    """Plot density evolution from one directory."""
    files, datasets = load_outputs(output_dir)
    if not datasets:
        print(f"No output files found in {output_dir}/")
        return

    n = len(datasets)
    n_cols = min(4, n)
    n_rows = (n + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 3 * n_rows),
                             squeeze=False)

    for i, (step, time, data) in enumerate(datasets):
        ax = axes[i // n_cols][i % n_cols]
        ax.plot(data[:, 0], data[:, 1], 'b-', linewidth=0.8)
        ax.set_title(f"{title_prefix}t = {time:.4f}", fontsize=10)
        ax.set_xlabel("x")
        ax.set_ylabel(r"$|\psi|^2$")
        ax.set_xlim(0, 1)

    for i in range(n, n_rows * n_cols):
        axes[i // n_cols][i % n_cols].set_visible(False)

    plt.tight_layout()
    outfile = os.path.join(output_dir, "density_evolution.png")
    plt.savefig(outfile, dpi=150)
    print(f"Saved {outfile}")


def plot_comparison():
    """Compare free particle vs gravity side by side."""
    dir_free = "output/sp_1d"
    dir_grav = "output/sp_1d_grav"

    _, ds_free = load_outputs(dir_free)
    _, ds_grav = load_outputs(dir_grav)

    if not ds_free or not ds_grav:
        print("Need both output/sp_1d/ and output/sp_1d_grav/ for comparison")
        return

    # Match timesteps
    n = min(len(ds_free), len(ds_grav))

    fig, axes = plt.subplots(n, 2, figsize=(10, 2.5 * n), squeeze=False)
    fig.suptitle("Free particle (left) vs Gravity (right)", fontsize=14, y=1.01)

    for i in range(n):
        step_f, time_f, data_f = ds_free[i]
        step_g, time_g, data_g = ds_grav[i]

        # Free particle
        axes[i][0].plot(data_f[:, 0], data_f[:, 1], 'b-', linewidth=0.8)
        axes[i][0].set_ylabel(r"$|\psi|^2$")
        axes[i][0].set_xlim(0, 1)
        if i == 0:
            axes[i][0].set_title("Free particle", fontsize=11)
        axes[i][0].text(0.02, 0.95, f"t = {time_f:.4f}",
                        transform=axes[i][0].transAxes, va='top', fontsize=9)

        # Gravity
        axes[i][1].plot(data_g[:, 0], data_g[:, 1], 'r-', linewidth=0.8)
        axes[i][1].set_xlim(0, 1)
        if i == 0:
            axes[i][1].set_title("With gravity", fontsize=11)
        axes[i][1].text(0.02, 0.95, f"t = {time_g:.4f}",
                        transform=axes[i][1].transAxes, va='top', fontsize=9)

        # Also show V on secondary axis if column exists
        if data_g.shape[1] > 5:
            ax2 = axes[i][1].twinx()
            ax2.plot(data_g[:, 0], data_g[:, 5], 'k--', linewidth=0.5, alpha=0.5)
            ax2.set_ylabel("V", fontsize=8, color='gray')
            ax2.tick_params(labelsize=7, colors='gray')

    axes[-1][0].set_xlabel("x")
    axes[-1][1].set_xlabel("x")

    plt.tight_layout()
    outfile = "output/sp_1d_comparison.png"
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    print(f"Saved {outfile}")


if __name__ == "__main__":
    # Plot individual directories
    for d in ["output/sp_1d", "output/sp_1d_grav"]:
        if os.path.isdir(d):
            label = "Gravity: " if "grav" in d else "Free: "
            plot_single(d, title_prefix=label)

    # Comparison plot if both exist
    if os.path.isdir("output/sp_1d") and os.path.isdir("output/sp_1d_grav"):
        plot_comparison()

    plt.show()
