#!/usr/bin/env python3
"""
Plot 3D FPA results (z=p/2 slices).
Reproduces Figure 4.4 from Thomson (2011) Chapter 4.

Usage:
    python3 plot_fpa_3d.py
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import glob
import os


def load_slice(filename):
    """Load a z=p/2 slice file, returning D, x, y, delta arrays."""
    with open(filename) as f:
        line = f.readline()
        D = float(line.split("=")[1])

    data = np.loadtxt(filename)
    x = data[:, 0]
    y = data[:, 1]
    delta = data[:, 2]
    vx = data[:, 3]
    vy = data[:, 4]
    vz = data[:, 5]

    # Reshape to 2D grid
    N = int(np.sqrt(len(x)))
    X = x.reshape(N, N)
    Y = y.reshape(N, N)
    Delta = delta.reshape(N, N)

    return D, X, Y, Delta, vx.reshape(N, N), vy.reshape(N, N)


def plot_density_surface(datasets, output="fpa_3d_density.png"):
    """Plot 3D surface of δ vs x/p, y/p for each output time (Fig 4.4 style)."""
    fig = plt.figure(figsize=(12, 10))

    for i, (D, X, Y, Delta, vx, vy) in enumerate(datasets):
        if i >= 4:
            break
        ax = fig.add_subplot(2, 2, i + 1, projection="3d")
        ax.plot_surface(X, Y, Delta, cmap="viridis", alpha=0.8,
                        rstride=1, cstride=1, linewidth=0, antialiased=True)
        ax.set_xlabel("x/d")
        ax.set_ylabel("y/d")
        ax.set_zlabel("δ")
        ax.set_title(f"D = {D:.0f}")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)

    fig.suptitle(
        f"3D FPA Toy Model: Density contrast (z=p/2 slice, δ_a=0.01, N=64³)",
        fontsize=13,
    )
    plt.tight_layout()
    plt.savefig(output, dpi=150)
    print(f"Saved {output}")
    plt.close()


def plot_density_contour(datasets, output="fpa_3d_contour.png"):
    """Plot 2D contour of δ at z=0 for each output time."""
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    axes = axes.flatten()

    for i, (D, X, Y, Delta, vx, vy) in enumerate(datasets):
        if i >= 4:
            break
        ax = axes[i]
        c = ax.contourf(X, Y, Delta, levels=50, cmap="RdBu_r")
        fig.colorbar(c, ax=ax, label="δ")
        ax.set_xlabel("x/d")
        ax.set_ylabel("y/d")
        ax.set_title(f"D = {D:.0f}")
        ax.set_aspect("equal")

    fig.suptitle(
        f"3D FPA Toy Model: Density contrast contour (z=p/2 slice)",
        fontsize=13,
    )
    plt.tight_layout()
    plt.savefig(output, dpi=150)
    print(f"Saved {output}")
    plt.close()


def plot_mass_conservation(datasets, output="fpa_3d_mass.png"):
    """Check mass on the slice (note: only full 3D mass is conserved)."""
    fig, ax = plt.subplots(figsize=(8, 4))

    D_vals = []
    mass_vals = []
    for D, X, Y, Delta, vx, vy in datasets:
        mass = np.mean(1.0 + Delta)
        D_vals.append(D)
        mass_vals.append(mass)

    ax.plot(D_vals, mass_vals, "ko-", markersize=8)
    ax.set_xlabel("Growth factor D")
    ax.set_ylabel("Mean density <ρ>/<ρ₀> on slice")
    ax.set_title("Slice density (full 3D mass is conserved to ~10⁻¹³)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    print("  Slice mean density:")
    for D, m in zip(D_vals, mass_vals):
        print(f"    D = {D:8.0f}  slice_mean = {m:.6f}")

    plt.tight_layout()
    plt.savefig(output, dpi=150)
    print(f"Saved {output}")
    plt.close()


def main():
    files = sorted(glob.glob("fpa_3d_D*.dat"))
    if not files:
        print("No output files found. Run ./fpa_3d first.")
        return

    datasets = [load_slice(f) for f in files]
    datasets.sort(key=lambda d: d[0])  # sort by D

    print(f"Found {len(datasets)} output files:")
    for D, X, Y, Delta, vx, vy in datasets:
        print(f"  D = {D:8.0f}  delta_range = [{Delta.min():.6f}, {Delta.max():.6f}]")

    plot_density_surface(datasets)
    plot_density_contour(datasets)
    plot_mass_conservation(datasets)
    print("Done.")


if __name__ == "__main__":
    main()
