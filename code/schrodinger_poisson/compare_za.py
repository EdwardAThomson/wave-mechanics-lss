#!/usr/bin/env python3
"""Compare S-P density fields with Zel'dovich approximation at matched redshifts."""

import numpy as np
import matplotlib.pyplot as plt
import struct
import glob
import re
import os

N = 64
N3 = N * N * N

# ================================================================
# Cosmological parameters (must match both codes)
# ================================================================
Omega_m0 = 0.27
Omega_Lambda = 0.73
z_init = 50.0
a_init = 1.0 / (1.0 + z_init)


def hubble(a):
    return np.sqrt(Omega_m0 / a**3 + Omega_Lambda)


def growth_factor(a_target):
    nsteps = 10000
    da = a_target / nsteps
    integral = 0.0
    for i in range(nsteps):
        a0 = max(i * da, 1e-10)
        a1 = (i + 0.5) * da
        a2 = (i + 1) * da
        f0 = 1.0 / (a0 * hubble(a0))**3
        f1 = 1.0 / (a1 * hubble(a1))**3
        f2 = 1.0 / (a2 * hubble(a2))**3
        integral += (f0 + 4*f1 + f2) * da / 6.0
    return 2.5 * Omega_m0 * hubble(a_target) * integral


D_init = growth_factor(a_init)
D_ratio = growth_factor(1.0) / D_init


def D_to_a(D_target):
    """Find scale factor a such that D(a)/D(a_init) = D_target."""
    a_lo, a_hi = a_init, 1.5
    for _ in range(100):
        a_mid = (a_lo + a_hi) / 2
        D_mid = growth_factor(a_mid) / D_init
        if D_mid < D_target:
            a_lo = a_mid
        else:
            a_hi = a_mid
    return (a_lo + a_hi) / 2


# ================================================================
# Load ZA density (binary format from fpa_3d_cosmo)
# ================================================================
def load_za_density(D_val):
    fname = f"../fpa_3d_cosmo/output/za/za_density_D{D_val}.bin"
    if not os.path.exists(fname):
        return None, None, None
    with open(fname, 'rb') as f:
        D = struct.unpack('d', f.read(8))[0]
        n = struct.unpack('i', f.read(4))[0]
        L = struct.unpack('d', f.read(8))[0]
        delta = np.frombuffer(f.read(), dtype=np.float64).reshape(n, n, n)
    return D, D_to_a(D), delta


# ================================================================
# Load S-P density (binary psi format)
# ================================================================
def load_sp_density(psi_file):
    with open(psi_file, 'rb') as f:
        step = np.frombuffer(f.read(4), dtype=np.int32)[0]
        n = np.frombuffer(f.read(4), dtype=np.int32)[0]
        a = np.frombuffer(f.read(8), dtype=np.float64)[0]
        nu = np.frombuffer(f.read(8), dtype=np.float64)[0]
        psi = np.frombuffer(f.read(), dtype=np.complex128).reshape(n, n, n)
    rho = np.abs(psi)**2
    delta = rho / np.mean(rho) - 1.0
    return step, a, delta


def find_closest_sp_snapshot(target_a, sp_dir):
    """Find the S-P binary snapshot closest to the target scale factor."""
    files = sorted(glob.glob(os.path.join(sp_dir, "psi_*.bin")))
    best_file = None
    best_da = float('inf')
    for fpath in files:
        with open(fpath, 'rb') as f:
            _ = f.read(4)  # step
            _ = f.read(4)  # N
            a = np.frombuffer(f.read(8), dtype=np.float64)[0]
        if abs(a - target_a) < best_da:
            best_da = abs(a - target_a)
            best_file = fpath
    return best_file


# ================================================================
# Main comparison
# ================================================================
sp_dir = "output/sp_cosmo_nu_0.000100"

D_values = [5, 10, 15, 20, 25, 30, 34]

print(f"{'D':>4}  {'a_ZA':>7}  {'z_ZA':>6}  {'a_SP':>7}  {'corr':>8}  "
      f"{'rms_ZA':>8}  {'rms_SP':>8}")
print("-" * 65)

results = []
for Dv in D_values:
    D, a_za, delta_za = load_za_density(Dv)
    if delta_za is None:
        continue

    sp_file = find_closest_sp_snapshot(a_za, sp_dir)
    if sp_file is None:
        continue

    step_sp, a_sp, delta_sp = load_sp_density(sp_file)
    z_za = 1.0 / a_za - 1.0

    # Point-by-point correlation
    corr = np.corrcoef(delta_za.ravel(), delta_sp.ravel())[0, 1]
    rms_za = np.std(delta_za)
    rms_sp = np.std(delta_sp)

    print(f"{Dv:4d}  {a_za:7.4f}  {z_za:6.2f}  {a_sp:7.4f}  {corr:8.4f}  "
          f"{rms_za:8.4f}  {rms_sp:8.4f}")

    results.append((Dv, a_za, z_za, a_sp, corr, rms_za, rms_sp, delta_za, delta_sp))

# ================================================================
# Plot: side-by-side z-slices and scatter
# ================================================================
n_show = min(len(results), 4)
indices = np.linspace(0, len(results) - 1, n_show, dtype=int)

fig, axes = plt.subplots(3, n_show, figsize=(4 * n_show, 10), squeeze=False)

for col, i in enumerate(indices):
    Dv, a_za, z_za, a_sp, corr, rms_za, rms_sp, dza, dsp = results[i]
    iz = N // 2

    # Row 0: ZA slice
    vmax = max(np.percentile(np.log(2 + dza[:, :, iz]), 99),
               np.percentile(np.log(2 + dsp[:, :, iz]), 99))
    vmin = min(np.percentile(np.log(2 + dza[:, :, iz]), 1),
               np.percentile(np.log(2 + dsp[:, :, iz]), 1))

    axes[0][col].imshow(np.log(2 + dza[:, :, iz]).T, origin='lower',
                        extent=[0, 1, 0, 1], cmap='inferno',
                        vmin=vmin, vmax=vmax)
    axes[0][col].set_title(f"ZA  z={z_za:.1f} (D={Dv})", fontsize=10)
    if col == 0:
        axes[0][col].set_ylabel("y")

    # Row 1: S-P slice
    axes[1][col].imshow(np.log(2 + dsp[:, :, iz]).T, origin='lower',
                        extent=[0, 1, 0, 1], cmap='inferno',
                        vmin=vmin, vmax=vmax)
    axes[1][col].set_title(f"S-P  a={a_sp:.3f} (r={corr:.3f})", fontsize=10)
    if col == 0:
        axes[1][col].set_ylabel("y")

    # Row 2: scatter
    sub = np.random.choice(N3, min(5000, N3), replace=False)
    axes[2][col].scatter(dza.ravel()[sub], dsp.ravel()[sub],
                         s=1, alpha=0.3, color='steelblue')
    lim = max(abs(dza).max(), abs(dsp).max()) * 0.5
    axes[2][col].plot([-1, lim], [-1, lim], 'k--', linewidth=0.5)
    axes[2][col].set_xlim(-1, lim)
    axes[2][col].set_ylim(-1, lim)
    axes[2][col].set_xlabel("delta (ZA)")
    axes[2][col].text(0.05, 0.95, f"r = {corr:.4f}", transform=axes[2][col].transAxes,
                      va='top', fontsize=10)
    if col == 0:
        axes[2][col].set_ylabel("delta (S-P)")

plt.suptitle("Zel'dovich Approximation vs Schrödinger-Poisson", fontsize=13)
plt.tight_layout()
outfile = os.path.join(sp_dir, "za_comparison.png")
plt.savefig(outfile, dpi=150, bbox_inches='tight')
print(f"\nSaved {outfile}")
plt.show()
