#!/usr/bin/env python3
"""Replot all cosmological density slices with a fixed color scale.

Every plot uses log2(2 + delta) with vmin=0, vmax=8, cmap='inferno'.
This ensures a given color means the same density everywhere.

Reference:
  black  = empty void (delta = -1)
  dark   = underdense (delta < 0)
  orange = mean density (delta = 0, log2 = 1.0)
  yellow = overdense
  white  = strongly collapsed (delta > 100)
"""

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import h5py
import glob
import os

# Universal color scale
VMIN = 0.0   # log2(2 + delta) for delta = -1 (empty)
VMAX = 8.0   # log2(2 + delta) for delta = 254
CMAP = 'inferno'


def load_sp(fpath):
    with open(fpath, 'rb') as f:
        step = np.frombuffer(f.read(4), dtype=np.int32)[0]
        n = np.frombuffer(f.read(4), dtype=np.int32)[0]
        a = np.frombuffer(f.read(8), dtype=np.float64)[0]
        nu = np.frombuffer(f.read(8), dtype=np.float64)[0]
        psi = np.frombuffer(f.read(), dtype=np.complex128).reshape(n, n, n)
    rho = np.abs(psi)**2
    return n, step, a, rho / np.mean(rho) - 1.0


def load_gadget_tsc(fpath, N_grid):
    with h5py.File(fpath, 'r') as f:
        pos = f['PartType1/Coordinates'][:].astype(np.float64)
        a = f['Header'].attrs['Time']
        boxsize = f['Header'].attrs['BoxSize']
    pos_norm = pos / boxsize
    rho = np.zeros((N_grid, N_grid, N_grid))
    dx = 1.0 / N_grid
    gx = pos_norm[:, 0] / dx; gy = pos_norm[:, 1] / dx; gz = pos_norm[:, 2] / dx
    ix0 = np.round(gx).astype(int); iy0 = np.round(gy).astype(int); iz0 = np.round(gz).astype(int)
    fx = gx - ix0; fy = gy - iy0; fz = gz - iz0
    def tsc_w(d):
        ad = np.abs(d); w = np.zeros_like(d)
        m1 = ad < 0.5; w[m1] = 0.75 - d[m1]**2
        m2 = (ad >= 0.5) & (ad < 1.5); w[m2] = 0.5 * (1.5 - ad[m2])**2
        return w
    for dix in [-1, 0, 1]:
        wx = tsc_w(fx - dix); iix = (ix0 + dix) % N_grid
        for diy in [-1, 0, 1]:
            wy = tsc_w(fy - diy); iiy = (iy0 + diy) % N_grid
            for diz in [-1, 0, 1]:
                wz = tsc_w(fz - diz); iiz = (iz0 + diz) % N_grid
                np.add.at(rho, (iix, iiy, iiz), wx * wy * wz)
    return a, rho / np.mean(rho) - 1.0


def plot_slice(ax, delta, N, title):
    iz = N // 2
    im = ax.imshow(np.log2(2 + delta[:, :, iz]).T, origin='lower',
                   extent=[0, 1, 0, 1], cmap=CMAP, vmin=VMIN, vmax=VMAX)
    ax.set_title(title, fontsize=10)
    return im


# ================================================================
# 1. N=128 density slices (nu=1e-3)
# ================================================================
def replot_n128_slices():
    files = sorted(glob.glob('output/sp_cosmo_n128_nu_0.001000/psi_*.bin'))
    if not files:
        print("  No N=128 nu=1e-3 snapshots found, skipping")
        return

    indices = [0, 2, 4, 6, 8, 10] if len(files) >= 11 else list(range(min(6, len(files))))
    ncols = min(len(indices), 3)
    nrows = (len(indices) + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 4.5 * nrows))
    if nrows == 1: axes = [axes]
    if ncols == 1: axes = [[a] for a in axes]

    for col, fi in enumerate(indices):
        N, step, a, delta = load_sp(files[fi])
        row = col // ncols; c = col % ncols
        im = plot_slice(axes[row][c], delta, N, f'z={1/a-1:.1f} (step {step})')
        if c == 0: axes[row][c].set_ylabel('y')

    fig.colorbar(im, ax=axes, label='log₂(2 + δ)', shrink=0.6)
    fig.suptitle('N=128, ν=10⁻³: density slices (fixed scale)', fontsize=13)
    fig.tight_layout()
    fig.savefig('output/sp_cosmo_n128_nu_0.001000/density_slices.png', dpi=150)
    plt.close()
    print("  Saved N=128 density_slices.png")


# ================================================================
# 2. N=256 density slices (nu=1e-4) — already done but regenerate
# ================================================================
def replot_n256_slices():
    files = sorted(glob.glob('output/sp_cosmo_n256_nu_0.000100/psi_*.bin'))
    if not files:
        print("  No N=256 snapshots found, skipping")
        return

    indices = np.linspace(0, len(files) - 1, min(10, len(files)), dtype=int)
    ncols = 5
    nrows = (len(indices) + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows))

    for col, fi in enumerate(indices):
        N, step, a, delta = load_sp(files[fi])
        row = col // ncols; c = col % ncols
        im = plot_slice(axes[row][c], delta, N, f'z={1/a-1:.1f}')
        if c == 0: axes[row][c].set_ylabel('y')

    fig.colorbar(im, ax=axes, label='log₂(2 + δ)', shrink=0.6)
    fig.suptitle('N=256, ν=10⁻⁴: density slices (fixed scale)', fontsize=13)
    fig.tight_layout()
    fig.savefig('output/sp_cosmo_n256_nu_0.000100/density_slices_n256.png', dpi=150)
    plt.close()
    print("  Saved N=256 density_slices_n256.png")


# ================================================================
# 3. GADGET full comparison (matched ICs)
# ================================================================
def replot_gadget_comparison():
    gadget_dir = '/home/edward/Projects/PhD/gadget4/runs/wm-comparison/output'
    sp_dir = 'output/sp_cosmo_gadget_ic_nu_0.000300'
    gadget_snaps = sorted(glob.glob(f'{gadget_dir}/snapshot_00[1-9].hdf5'))
    sp_files = sorted(glob.glob(f'{sp_dir}/psi_*.bin'))

    if not gadget_snaps or not sp_files:
        print("  Missing GADGET or S-P snapshots, skipping")
        return

    # Match snapshots
    pairs = []
    for gsnap in gadget_snaps:
        with h5py.File(gsnap, 'r') as f:
            a_g = float(f['Header'].attrs['Time'])
        best_da = float('inf'); best_f = None
        for sf in sp_files:
            with open(sf, 'rb') as f:
                f.read(4); f.read(4)
                a_sp = np.frombuffer(f.read(8), dtype=np.float64)[0]
            if abs(a_sp - a_g) < best_da:
                best_da = abs(a_sp - a_g); best_f = sf
        pairs.append((gsnap, best_f, a_g))

    ncols = len(pairs)
    fig, axes = plt.subplots(2, ncols, figsize=(3.5 * ncols, 7), squeeze=False)

    for col, (gsnap, sp_file, a_g) in enumerate(pairs):
        z_g = 1.0 / a_g - 1.0
        N_g = 128
        a_gd, delta_g = load_gadget_tsc(gsnap, N_g)
        N_s, step, a_sp, delta_sp = load_sp(sp_file)

        corr = np.corrcoef(delta_g.ravel(), delta_sp.ravel())[0, 1]

        plot_slice(axes[0][col], delta_g, N_g, f'GADGET z={z_g:.1f}')
        if col == 0: axes[0][col].set_ylabel('y')
        im = plot_slice(axes[1][col], delta_sp, N_s, f'S-P r={corr:.3f}')
        if col == 0: axes[1][col].set_ylabel('y')

    fig.colorbar(im, ax=axes, label='log₂(2 + δ)', shrink=0.6)
    fig.suptitle('GADGET-4 vs S-P (matched ICs, N=128, fixed scale)', fontsize=13)
    fig.tight_layout()
    fig.savefig(f'{sp_dir}/gadget_full_comparison.png', dpi=150)
    plt.close()
    print("  Saved gadget_full_comparison.png")


# ================================================================
# 4. ZA comparison (N=128)
# ================================================================
def replot_za_comparison():
    sp_dir = 'output/sp_cosmo_n128_nu_0.001000'
    za_dir = '../fpa_3d_cosmo/output/za_N128'
    sp_files = sorted(glob.glob(f'{sp_dir}/psi_*.bin'))

    import struct
    za_D = [5, 15, 25, 34]
    Omega_m0 = 0.27; Omega_Lambda = 0.73; a_init = 1.0 / 51.0

    def hubble(a): return np.sqrt(Omega_m0 / a**3 + Omega_Lambda)
    def growth_factor(a):
        nsteps = 10000; da = a / nsteps; s = 0.0
        for i in range(nsteps):
            a0 = max(i * da, 1e-10); a1 = (i + 0.5) * da; a2 = (i + 1) * da
            s += (1 / (a0 * hubble(a0))**3 + 4 / (a1 * hubble(a1))**3 + 1 / (a2 * hubble(a2))**3) * da / 6
        return 2.5 * Omega_m0 * hubble(a) * s
    D_init = growth_factor(a_init)

    def D_to_a(D):
        lo, hi = a_init, 1.5
        for _ in range(100):
            mid = (lo + hi) / 2
            if growth_factor(mid) / D_init < D: lo = mid
            else: hi = mid
        return (lo + hi) / 2

    if not sp_files or not os.path.exists(za_dir):
        print("  Missing ZA or S-P data, skipping")
        return

    pairs = []
    for Dv in za_D:
        fname = f'{za_dir}/za_density_D{Dv}.bin'
        if not os.path.exists(fname): continue
        with open(fname, 'rb') as f:
            D = struct.unpack('d', f.read(8))[0]
            n = struct.unpack('i', f.read(4))[0]
            L = struct.unpack('d', f.read(8))[0]
            delta_za = np.frombuffer(f.read(), dtype=np.float64).reshape(n, n, n)
        a_za = D_to_a(Dv)
        best_da = float('inf'); best_f = None
        for sf in sp_files:
            with open(sf, 'rb') as f:
                f.read(4); f.read(4)
                a_sp = np.frombuffer(f.read(8), dtype=np.float64)[0]
            if abs(a_sp - a_za) < best_da:
                best_da = abs(a_sp - a_za); best_f = sf
        N_s, step, a_sp, delta_sp = load_sp(best_f)
        corr = np.corrcoef(delta_za.ravel(), delta_sp.ravel())[0, 1]
        pairs.append((Dv, a_za, delta_za, n, delta_sp, N_s, corr))

    if not pairs: return

    ncols = len(pairs)
    fig, axes = plt.subplots(2, ncols, figsize=(4.5 * ncols, 8), squeeze=False)

    for col, (Dv, a_za, dza, N_za, dsp, N_sp, corr) in enumerate(pairs):
        z_za = 1.0 / a_za - 1.0
        plot_slice(axes[0][col], dza, N_za, f'ZA z={z_za:.1f} (D={Dv})')
        if col == 0: axes[0][col].set_ylabel('y')
        im = plot_slice(axes[1][col], dsp, N_sp, f'S-P r={corr:.3f}')
        if col == 0: axes[1][col].set_ylabel('y')

    fig.colorbar(im, ax=axes, label='log₂(2 + δ)', shrink=0.6)
    fig.suptitle('ZA vs S-P (N=128, ν=10⁻³, fixed scale)', fontsize=13)
    fig.tight_layout()
    fig.savefig(f'{sp_dir}/za_comparison_n128.png', dpi=150)
    plt.close()
    print("  Saved za_comparison_n128.png")


# ================================================================
# Run all
# ================================================================
if __name__ == "__main__":
    print("Replotting all density slices with fixed color scale...")
    print(f"  vmin={VMIN}, vmax={VMAX}, cmap={CMAP}")
    print()

    print("1. N=128 density slices:")
    replot_n128_slices()

    print("2. N=256 density slices:")
    replot_n256_slices()

    print("3. GADGET comparison:")
    replot_gadget_comparison()

    print("4. ZA comparison:")
    replot_za_comparison()

    print("\nDone.")
