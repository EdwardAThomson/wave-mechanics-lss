#!/usr/bin/env python3
"""
Plot 1D FPA results.
Reproduces Figure 4.1 (density) and Figure 4.2 (velocity) from
Thomson (2011) Chapter 4.

Usage:
    python3 plot_fpa.py
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import os


def load_data(filename):
    """Load an FPA output file, returning x, delta, log(2+delta), velocity."""
    data = np.loadtxt(filename)
    # Parse D from header
    with open(filename) as f:
        line = f.readline()
        D = float(line.split("=")[1])
    return D, data[:, 0], data[:, 1], data[:, 2], data[:, 3]


def plot_density(datasets, output="fpa_1d_density.png"):
    """Plot log(2 + delta) vs x/p for all output times, matching Fig 4.1."""
    n = len(datasets)
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    axes = axes.flatten()

    for i, (D, x, delta, log_delta, vel) in enumerate(datasets):
        if i >= 4:
            break
        ax = axes[i]
        ax.plot(x, log_delta, "k-", linewidth=1.5)
        ax.set_xlabel("x/d")
        ax.set_ylabel("log(2 + δ)")
        ax.set_title(f"D = {D:.2f}")
        ax.set_xlim(0, 1)
        ax.set_ylim(0.294, 0.306)
        ax.grid(True, alpha=0.3)

    fig.suptitle(
        f"1D FPA: Density contrast (δ_a={0.01}, ν={1.0}, N={512})",
        fontsize=13,
    )
    plt.tight_layout()
    plt.savefig(output, dpi=150)
    print(f"Saved {output}")
    plt.close()


def plot_velocity(datasets, output="fpa_1d_velocity.png"):
    """Plot velocity vs x/p for all output times, matching Fig 4.2."""
    n = len(datasets)
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    axes = axes.flatten()

    for i, (D, x, delta, log_delta, vel) in enumerate(datasets):
        if i >= 4:
            break
        ax = axes[i]
        ax.plot(x, vel, "k-", linewidth=1.5)
        ax.set_xlabel("x/d")
        ax.set_ylabel("velocity")
        ax.set_title(f"D = {D:.2f}")
        ax.set_xlim(0, 1)
        ax.grid(True, alpha=0.3)

    fig.suptitle(
        f"1D FPA: Velocity field (δ_a={0.01}, ν={1.0}, N={512})",
        fontsize=13,
    )
    plt.tight_layout()
    plt.savefig(output, dpi=150)
    print(f"Saved {output}")
    plt.close()


def plot_mass_conservation(datasets, output="fpa_1d_mass.png"):
    """Plot total mass at each output time to verify conservation."""
    fig, ax = plt.subplots(figsize=(8, 4))

    D_vals = []
    mass_vals = []
    for D, x, delta, log_delta, vel in datasets:
        # mass = integral of (1+delta) dx, but since <delta>=0, mass = 1
        # More precisely: sum of |psi|^2 * dx = sum of (1+delta)*<rho> * dx
        rho = 10**log_delta - 2.0 + 1.0  # recover delta from log(2+delta), then rho=1+delta
        # Actually just use delta directly
        mass = np.mean(1.0 + delta)
        D_vals.append(D)
        mass_vals.append(mass)

    ax.plot(D_vals, mass_vals, "ko-", markersize=8)
    ax.axhline(y=1.0, color="r", linestyle="--", alpha=0.5, label="Expected")
    ax.set_xlabel("Growth factor D")
    ax.set_ylabel("Mean density <ρ>/<ρ₀>")
    ax.set_title("Mass conservation check")
    ax.set_ylim(0.999, 1.001)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Print actual values so user can see the precision
    print("  Mass conservation:")
    for D, m in zip(D_vals, mass_vals):
        print(f"    D = {D:8.2f}  mass = {m:.15f}  error = {abs(m - 1.0):.2e}")

    plt.tight_layout()
    plt.savefig(output, dpi=150)
    print(f"Saved {output}")
    plt.close()


def main():
    # Find all output files
    files = sorted(glob.glob("fpa_1d_D*.dat"))
    if not files:
        print("No output files found. Run ./fpa_1d first.")
        return

    datasets = [load_data(f) for f in files]
    datasets.sort(key=lambda d: d[0])  # sort by D

    print(f"Found {len(datasets)} output files:")
    for D, x, delta, log_delta, vel in datasets:
        print(f"  D = {D:8.2f}  delta_range = [{delta.min():.6f}, {delta.max():.6f}]")

    plot_density(datasets)
    plot_velocity(datasets)
    plot_mass_conservation(datasets)
    print("Done.")


if __name__ == "__main__":
    main()
