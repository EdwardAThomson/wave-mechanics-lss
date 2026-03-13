#!/usr/bin/env python3
"""
Diagnostic plots for 3D Zel'dovich/FPA with cosmological initial conditions.
Reproduces thesis Figures 4.5-4.7 style plots plus validation diagnostics.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.stats import norm
import struct
import glob
import os

OUTDIR = 'output/za'

plt.rcParams.update({
    'font.size': 12,
    'figure.figsize': (10, 8),
    'figure.dpi': 150,
})


def read_slice(filename):
    """Read a z-slice .dat file (x, y, delta)."""
    data = np.loadtxt(filename, comments='#')
    N = int(np.sqrt(len(data)))
    x = data[:, 0].reshape(N, N)
    y = data[:, 1].reshape(N, N)
    delta = data[:, 2].reshape(N, N)
    return x, y, delta


def read_density_3d(filename):
    """Read binary 3D density file."""
    with open(filename, 'rb') as f:
        D = struct.unpack('d', f.read(8))[0]
        N = struct.unpack('i', f.read(4))[0]
        L = struct.unpack('d', f.read(8))[0]
        delta = np.frombuffer(f.read(N**3 * 8), dtype=np.float64)
    return D, N, L, delta.reshape(N, N, N)


def read_histogram(filename):
    """Read histogram .dat file."""
    data = np.loadtxt(filename, comments='#')
    return data[:, 0], data[:, 1], data[:, 2]


def compute_power_spectrum(delta_3d, L, N):
    """Compute spherically-averaged power spectrum from 3D density field."""
    delta_k = np.fft.fftn(delta_3d)
    Pk_3d = np.abs(delta_k)**2 * (L / N)**3 / L**3

    kx = np.fft.fftfreq(N, d=L/N) * 2 * np.pi
    ky = np.fft.fftfreq(N, d=L/N) * 2 * np.pi
    kz = np.fft.fftfreq(N, d=L/N) * 2 * np.pi
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
    K = np.sqrt(KX**2 + KY**2 + KZ**2)

    k_fund = 2 * np.pi / L
    k_nyq = np.pi * N / L
    k_bins = np.arange(k_fund, k_nyq, k_fund)

    k_avg = []
    P_avg = []
    for i in range(len(k_bins) - 1):
        mask = (K >= k_bins[i]) & (K < k_bins[i+1])
        if np.sum(mask) > 0:
            k_avg.append(np.mean(K[mask]))
            P_avg.append(np.mean(Pk_3d[mask]))

    return np.array(k_avg), np.array(P_avg)


# ============================================================
# Figure 1: Surface plots of density at z=L/2 (Fig 4.5 style)
# ============================================================
slice_files = sorted(glob.glob(f'{OUTDIR}/za_slice_D*.dat'),
                     key=lambda f: int(f.split('_D')[1].replace('.dat', '')))
# Show subset matching thesis figures: D = 1, 15, 25, 34
thesis_D = {1, 15, 25, 34}
slice_subset = [f for f in slice_files if int(f.split('_D')[1].replace('.dat', '')) in thesis_D]
if slice_subset:
    n_panels = len(slice_subset)
    fig, axes = plt.subplots(1, n_panels, figsize=(5*n_panels, 5),
                             subplot_kw={'projection': '3d'})
    if n_panels == 1:
        axes = [axes]

    for ax, fname in zip(axes, slice_subset):
        x, y, delta = read_slice(fname)
        D = fname.split('_D')[1].replace('.dat', '')
        log_delta = np.log10(2.0 + delta)
        ax.plot_surface(x, y, log_delta, cmap='viridis',
                       linewidth=0, antialiased=True, alpha=0.9,
                       rstride=1, cstride=1)
        ax.set_title(f'D = {D}', fontsize=14)
        ax.set_xlabel('x/L')
        ax.set_ylabel('y/L')
        ax.set_zlabel(r'$\log_{10}(2+\delta)$')
        ax.view_init(elev=30, azim=-60)

    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/fig_density_slices.png', bbox_inches='tight')
    plt.close()
    print("Saved fig_density_slices.png")


# ============================================================
# Figure 2: Initial delta histogram (Fig 4.8 style)
# ============================================================
if os.path.exists(f'{OUTDIR}/delta_init_hist.dat'):
    centers, counts, normed = read_histogram(f'{OUTDIR}/delta_init_hist.dat')
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.bar(centers, normed, width=centers[1]-centers[0], alpha=0.7, label='Measured')

    mu = np.sum(centers * normed) * (centers[1] - centers[0])
    sig = np.sqrt(np.sum((centers - mu)**2 * normed) * (centers[1] - centers[0]))
    x_fit = np.linspace(centers[0], centers[-1], 200)
    ax.plot(x_fit, norm.pdf(x_fit, mu, sig), 'r-', lw=2,
            label=f'Gaussian: $\\mu$={mu:.4f}, $\\sigma$={sig:.4f}')

    ax.set_xlabel(r'$\delta$')
    ax.set_ylabel('PDF')
    ax.set_title('Initial Density Contrast Distribution')
    ax.legend()
    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/fig_delta_histogram.png', bbox_inches='tight')
    plt.close()
    print("Saved fig_delta_histogram.png")


# ============================================================
# Figure 3: Displacement component histograms (Fig 4.9 style)
# ============================================================
vel_files = [f'{OUTDIR}/vel_x_init_hist.dat', f'{OUTDIR}/vel_y_init_hist.dat', f'{OUTDIR}/vel_z_init_hist.dat']
vel_labels = [r'$s_x$', r'$s_y$', r'$s_z$']
if all(os.path.exists(f) for f in vel_files):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, fname, label in zip(axes, vel_files, vel_labels):
        centers, counts, normed = read_histogram(fname)
        ax.bar(centers, normed, width=centers[1]-centers[0], alpha=0.7)

        mu = np.sum(centers * normed) * (centers[1] - centers[0])
        sig = np.sqrt(np.sum((centers - mu)**2 * normed) * (centers[1] - centers[0]))
        x_fit = np.linspace(centers[0], centers[-1], 200)
        ax.plot(x_fit, norm.pdf(x_fit, mu, sig), 'r-', lw=2,
                label=f'$\\sigma$={sig:.2e}')

        ax.set_xlabel(label + ' (box units)')
        ax.set_ylabel('PDF')
        ax.legend()

    plt.suptitle('Initial Displacement Distributions')
    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/fig_velocity_histograms.png', bbox_inches='tight')
    plt.close()
    print("Saved fig_velocity_histograms.png")


# ============================================================
# Figure 4: Power spectrum P(k) at each D
# ============================================================
density_files = sorted(glob.glob(f'{OUTDIR}/za_density_D*.bin'),
                       key=lambda f: int(f.split('_D')[1].replace('.bin', '')))
if density_files:
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = plt.cm.viridis(np.linspace(0, 1, len(density_files)))

    for fname, color in zip(density_files, colors):
        D, N, L, delta_3d = read_density_3d(fname)
        k, Pk = compute_power_spectrum(delta_3d, L, N)
        ax.loglog(k, Pk, '-', color=color, label=f'D = {D:.0f}', lw=1.5)

    # Overlay BBKS theory line
    h_cosmo = 0.71
    Omega_m = 0.27
    Gamma = Omega_m * h_cosmo
    k_theory = np.logspace(-2, 1, 500)
    q = k_theory / Gamma
    T = np.log(1 + 2.34*q) / (2.34*q) * (1 + 3.89*q + (16.1*q)**2 + (5.46*q)**3 + (6.71*q)**4)**(-0.25)
    Pk_theory = k_theory * T**2

    # Match to D=1 data for normalisation
    if len(density_files) > 0:
        D0, N0, L0, d0 = read_density_3d(density_files[0])
        k0, Pk0 = compute_power_spectrum(d0, L0, N0)
        idx_match = np.argmin(np.abs(k0 - 0.5))
        idx_theory = np.argmin(np.abs(k_theory - 0.5))
        if Pk_theory[idx_theory] > 0 and len(k0) > idx_match:
            norm_factor = Pk0[idx_match] / Pk_theory[idx_theory]
            ax.loglog(k_theory, Pk_theory * norm_factor, 'k--', alpha=0.5,
                     label='BBKS (scaled to D=1)')

    ax.set_xlabel(r'$k$ [h/Mpc]')
    ax.set_ylabel(r'$P(k)$ [(Mpc/h)$^3$]')
    ax.set_title('Power Spectrum Evolution (Zel\'dovich Approx.)')
    ax.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/fig_power_spectrum.png', bbox_inches='tight')
    plt.close()
    print("Saved fig_power_spectrum.png")


# ============================================================
# Figure 5: 2D density colourmap slices
# ============================================================
if slice_files:
    n_panels = len(slice_files)
    fig, axes = plt.subplots(1, n_panels, figsize=(4*n_panels, 3.5))
    if n_panels == 1:
        axes = [axes]

    for ax, fname in zip(axes, slice_files):
        x, y, delta = read_slice(fname)
        D = fname.split('_D')[1].replace('.dat', '')
        log_delta = np.log10(2.0 + delta)
        im = ax.pcolormesh(x, y, log_delta, cmap='inferno', shading='auto')
        ax.set_title(f'D = {D}')
        ax.set_xlabel('x/L')
        ax.set_ylabel('y/L')
        ax.set_aspect('equal')
        plt.colorbar(im, ax=ax, label=r'$\log_{10}(2+\delta)$')

    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/fig_density_colourmap.png', bbox_inches='tight')
    plt.close()
    print("Saved fig_density_colourmap.png")


# ============================================================
# Figure 6: Growth rate check
# ============================================================
if len(density_files) >= 2:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    D_vals = []
    rms_vals = []
    for fname in density_files:
        D, N, L, delta_3d = read_density_3d(fname)
        rms = np.sqrt(np.mean(delta_3d**2))
        D_vals.append(D)
        rms_vals.append(rms)

    D_vals = np.array(D_vals)
    rms_vals = np.array(rms_vals)

    # Skip D=1 (delta=0) for plotting
    mask = np.array(rms_vals) > 0
    D_plot = D_vals[mask]
    rms_plot = rms_vals[mask]

    ax1.plot(D_plot, rms_plot, 'bo-', label=r'$\sigma_\delta$ (measured)')
    sigma_field = 0.112
    D_linear = np.linspace(2, D_vals[-1], 100)
    ax1.plot(D_linear, sigma_field * (D_linear - 1), 'r--', alpha=0.7,
            label=r'Linear ZA: $\sigma_0 \times (D-1)$')
    ax1.set_xlabel('D (growth factor)')
    ax1.set_ylabel(r'$\sigma_\delta$ (rms)')
    ax1.set_title('Growth of Density Fluctuations')
    ax1.legend()
    ax1.set_yscale('log')
    ax1.set_xscale('log')

    # Correlation with first non-zero density (D=5)
    ref_idx = next(i for i, r in enumerate(rms_vals) if r > 0)
    D_ref, N_ref, L_ref, d_ref = read_density_3d(density_files[ref_idx])
    corr_vals = []
    D_corr = []
    for i, fname in enumerate(density_files):
        if rms_vals[i] == 0:
            continue
        D, Ni, Li, di = read_density_3d(fname)
        norm = np.sqrt(np.sum(di**2) * np.sum(d_ref**2))
        corr = np.sum(di * d_ref) / norm if norm > 0 else 0
        corr_vals.append(corr)
        D_corr.append(D)

    ax2.plot(D_corr, corr_vals, 'go-')
    ax2.set_xlabel('D (growth factor)')
    ax2.set_ylabel(f'Correlation with D={D_ref:.0f}')
    ax2.set_title('Spatial Correlation (Structure Persistence)')
    ax2.axhline(y=0, color='k', linestyle=':', alpha=0.3)
    ax2.set_ylim(-0.1, 1.1)

    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/fig_growth_check.png', bbox_inches='tight')
    plt.close()
    print("Saved fig_growth_check.png")


print("\nAll plots generated.")
