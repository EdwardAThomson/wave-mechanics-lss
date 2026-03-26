#!/usr/bin/env python3
"""Statistical comparison between S-P and GADGET-4 density fields.

Computes: power spectrum P(k), density PDF, cross-power spectrum,
rms growth, and cumulative density distribution at matched redshifts.

Usage: python3 compare_statistics.py
"""

import numpy as np
import matplotlib.pyplot as plt
import h5py
import glob
import os

N = 128
N3 = N * N * N


# ================================================================
# Loading functions
# ================================================================

def load_sp(fpath):
    with open(fpath, 'rb') as f:
        step = np.frombuffer(f.read(4), dtype=np.int32)[0]
        n = np.frombuffer(f.read(4), dtype=np.int32)[0]
        a = np.frombuffer(f.read(8), dtype=np.float64)[0]
        nu = np.frombuffer(f.read(8), dtype=np.float64)[0]
        psi = np.frombuffer(f.read(), dtype=np.complex128).reshape(n, n, n)
    rho = np.abs(psi)**2
    return a, rho / np.mean(rho) - 1.0


def load_gadget_tsc(fpath):
    with h5py.File(fpath, 'r') as f:
        pos = f['PartType1/Coordinates'][:].astype(np.float64)
        a = f['Header'].attrs['Time']
        boxsize = f['Header'].attrs['BoxSize']
    pos_norm = pos / boxsize
    rho = np.zeros((N, N, N))
    dx = 1.0 / N
    gx = pos_norm[:, 0] / dx
    gy = pos_norm[:, 1] / dx
    gz = pos_norm[:, 2] / dx
    ix0 = np.round(gx).astype(int)
    iy0 = np.round(gy).astype(int)
    iz0 = np.round(gz).astype(int)
    fx = gx - ix0; fy = gy - iy0; fz = gz - iz0

    def tsc_w(d):
        ad = np.abs(d); w = np.zeros_like(d)
        m1 = ad < 0.5; w[m1] = 0.75 - d[m1]**2
        m2 = (ad >= 0.5) & (ad < 1.5); w[m2] = 0.5 * (1.5 - ad[m2])**2
        return w

    for dix in [-1, 0, 1]:
        wx = tsc_w(fx - dix); iix = (ix0 + dix) % N
        for diy in [-1, 0, 1]:
            wy = tsc_w(fy - diy); iiy = (iy0 + diy) % N
            for diz in [-1, 0, 1]:
                wz = tsc_w(fz - diz); iiz = (iz0 + diz) % N
                np.add.at(rho, (iix, iiy, iiz), wx * wy * wz)
    return a, rho / np.mean(rho) - 1.0


# ================================================================
# Statistical tools
# ================================================================

def power_spectrum(delta, L_box=32.0):
    """Compute spherically averaged power spectrum P(k).
    Returns k (h/Mpc), P(k) (Mpc/h)^3."""
    dk = np.fft.fftn(delta)
    pk = np.abs(dk)**2 / N3

    kx = np.fft.fftfreq(N, d=L_box / N)  # in h/Mpc
    KX, KY, KZ = np.meshgrid(kx, kx, kx, indexing='ij')
    K = np.sqrt(KX**2 + KY**2 + KZ**2)

    # Bin by |k|
    k_fund = 1.0 / L_box  # fundamental frequency
    k_nyq = N / (2.0 * L_box)
    k_bins = np.arange(k_fund, k_nyq, k_fund)
    k_centres = []
    pk_binned = []

    for i in range(len(k_bins) - 1):
        mask = (K >= k_bins[i]) & (K < k_bins[i + 1])
        if mask.sum() > 0:
            k_centres.append(0.5 * (k_bins[i] + k_bins[i + 1]))
            # P(k) = V * <|delta_k|^2> where V = L^3
            pk_binned.append(np.mean(pk[mask]) * L_box**3)

    return np.array(k_centres), np.array(pk_binned)


def cross_power_spectrum(delta1, delta2, L_box=32.0):
    """Compute cross-power spectrum and transfer function."""
    dk1 = np.fft.fftn(delta1)
    dk2 = np.fft.fftn(delta2)

    kx = np.fft.fftfreq(N, d=L_box / N)
    KX, KY, KZ = np.meshgrid(kx, kx, kx, indexing='ij')
    K = np.sqrt(KX**2 + KY**2 + KZ**2)

    pk_cross = np.real(dk1 * np.conj(dk2)) / N3
    pk_11 = np.abs(dk1)**2 / N3
    pk_22 = np.abs(dk2)**2 / N3

    k_fund = 1.0 / L_box
    k_nyq = N / (2.0 * L_box)
    k_bins = np.arange(k_fund, k_nyq, k_fund)
    k_centres = []
    r_k = []      # cross-correlation coefficient r(k)
    t_k = []      # transfer function sqrt(P_SP/P_G4)

    for i in range(len(k_bins) - 1):
        mask = (K >= k_bins[i]) & (K < k_bins[i + 1])
        if mask.sum() > 0:
            p12 = np.mean(pk_cross[mask])
            p11 = np.mean(pk_11[mask])
            p22 = np.mean(pk_22[mask])
            k_centres.append(0.5 * (k_bins[i] + k_bins[i + 1]))
            r_k.append(p12 / np.sqrt(p11 * p22) if p11 > 0 and p22 > 0 else 0)
            t_k.append(np.sqrt(p22 / p11) if p11 > 0 else 0)

    return np.array(k_centres), np.array(r_k), np.array(t_k)


def density_pdf(delta, nbins=100):
    """Compute density PDF (histogram of delta)."""
    # Use log(1+delta) for better dynamic range
    log_rho = np.log10(1 + delta.ravel())
    log_rho = log_rho[np.isfinite(log_rho)]
    hist, edges = np.histogram(log_rho, bins=nbins, density=True)
    centres = 0.5 * (edges[:-1] + edges[1:])
    return centres, hist


def cumulative_density(delta, nbins=200):
    """Fraction of volume above a given delta threshold."""
    sorted_d = np.sort(delta.ravel())
    n = len(sorted_d)
    thresholds = sorted_d[::max(1, n // nbins)]
    fractions = 1.0 - np.searchsorted(sorted_d, thresholds) / n
    return thresholds, fractions


# ================================================================
# Main
# ================================================================

def main():
    gadget_dir = '/home/edward/Projects/PhD/gadget4/runs/wm-comparison/output'
    sp_dir = 'output/sp_cosmo_gadget_ic_nu_0.000300'
    out_dir = os.path.join(sp_dir, 'statistics')
    os.makedirs(out_dir, exist_ok=True)

    gadget_snaps = sorted(glob.glob(f'{gadget_dir}/snapshot_00[1-9].hdf5'))
    sp_files = sorted(glob.glob(f'{sp_dir}/psi_*.bin'))

    if not gadget_snaps:
        print("No GADGET snapshots found!")
        return
    if not sp_files:
        print("No S-P snapshots found!")
        return

    # Match snapshots
    pairs = []
    for gsnap in gadget_snaps:
        with h5py.File(gsnap, 'r') as f:
            a_g = float(f['Header'].attrs['Time'])
        best_da = float('inf')
        best_f = None
        for sf in sp_files:
            with open(sf, 'rb') as f:
                f.read(4); f.read(4)
                a_sp = np.frombuffer(f.read(8), dtype=np.float64)[0]
            if abs(a_sp - a_g) < best_da:
                best_da = abs(a_sp - a_g)
                best_f = sf
        pairs.append((gsnap, best_f, a_g))

    # ============================================================
    # Figure 1: Power spectra at each redshift
    # ============================================================
    n_pairs = len(pairs)
    fig1, axes1 = plt.subplots(2, n_pairs, figsize=(4 * n_pairs, 8), squeeze=False)

    rms_z = []; rms_g4 = []; rms_sp = []; corr_all = []

    print(f"{'z':>6}  {'r_total':>8}  {'rms_G4':>8}  {'rms_SP':>8}")
    print("-" * 36)

    for col, (gsnap, sp_file, a_g) in enumerate(pairs):
        z_g = 1.0 / a_g - 1.0
        a_g_actual, delta_g = load_gadget_tsc(gsnap)
        a_sp, delta_sp = load_sp(sp_file)

        # Power spectra
        k_g, pk_g = power_spectrum(delta_g)
        k_s, pk_s = power_spectrum(delta_sp)

        # Cross-power
        k_x, r_k, t_k = cross_power_spectrum(delta_g, delta_sp)

        # Total correlation
        r_total = np.corrcoef(delta_g.ravel(), delta_sp.ravel())[0, 1]

        rms_z.append(z_g)
        rms_g4.append(np.std(delta_g))
        rms_sp.append(np.std(delta_sp))
        corr_all.append(r_total)

        print(f"{z_g:6.1f}  {r_total:8.4f}  {np.std(delta_g):8.4f}  {np.std(delta_sp):8.4f}")

        # Row 0: P(k)
        ax = axes1[0][col]
        ax.loglog(k_g, pk_g, 'b-', linewidth=1.5, label='GADGET-4')
        ax.loglog(k_s, pk_s, 'r-', linewidth=1.5, label='S-P')
        ax.set_xlabel('k (h/Mpc)')
        ax.set_title(f'z = {z_g:.1f}', fontsize=11)
        if col == 0:
            ax.set_ylabel('P(k) (Mpc/h)³')
            ax.legend(fontsize=8)

        # Row 1: cross-correlation r(k) and transfer function
        ax2 = axes1[1][col]
        ax2.semilogx(k_x, r_k, 'k-', linewidth=1.5, label='r(k)')
        ax2.semilogx(k_x, t_k, 'r--', linewidth=1, label='T(k) = √(P_SP/P_G4)')
        ax2.axhline(1.0, color='gray', linewidth=0.5, linestyle=':')
        ax2.axhline(0.0, color='gray', linewidth=0.5, linestyle=':')
        ax2.set_xlabel('k (h/Mpc)')
        ax2.set_ylim(-0.5, 3.0)
        if col == 0:
            ax2.set_ylabel('r(k) / T(k)')
            ax2.legend(fontsize=8)

    fig1.suptitle('Power Spectrum & Cross-correlation: GADGET-4 vs S-P', fontsize=13)
    fig1.tight_layout()
    fig1.savefig(os.path.join(out_dir, 'power_spectra.png'), dpi=150)
    print(f"\nSaved {out_dir}/power_spectra.png")

    # ============================================================
    # Figure 2: Density PDFs
    # ============================================================
    n_show = min(n_pairs, 4)
    idx_show = np.linspace(0, n_pairs - 1, n_show, dtype=int)
    fig2, axes2 = plt.subplots(1, n_show, figsize=(4 * n_show, 4), squeeze=False)

    for col_i, idx in enumerate(idx_show):
        gsnap, sp_file, a_g = pairs[idx]
        z_g = 1.0 / a_g - 1.0
        _, delta_g = load_gadget_tsc(gsnap)
        _, delta_sp = load_sp(sp_file)

        c_g, h_g = density_pdf(delta_g)
        c_s, h_s = density_pdf(delta_sp)

        ax = axes2[0][col_i]
        ax.semilogy(c_g, h_g, 'b-', linewidth=1.5, label='GADGET-4')
        ax.semilogy(c_s, h_s, 'r-', linewidth=1.5, label='S-P')
        ax.set_xlabel('log₁₀(1 + δ)')
        ax.set_title(f'z = {z_g:.1f}', fontsize=11)
        if col_i == 0:
            ax.set_ylabel('PDF')
            ax.legend(fontsize=8)

    fig2.suptitle('Density PDF: GADGET-4 vs S-P', fontsize=13)
    fig2.tight_layout()
    fig2.savefig(os.path.join(out_dir, 'density_pdfs.png'), dpi=150)
    print(f"Saved {out_dir}/density_pdfs.png")

    # ============================================================
    # Figure 3: Summary — rms growth + correlation vs redshift
    # ============================================================
    fig3, (ax3a, ax3b) = plt.subplots(1, 2, figsize=(11, 4.5))

    ax3a.plot(rms_z, rms_g4, 'bs-', markersize=5, label='GADGET-4')
    ax3a.plot(rms_z, rms_sp, 'ro-', markersize=5, label='S-P')
    ax3a.set_xlabel('Redshift z')
    ax3a.set_ylabel('rms(δ)')
    ax3a.set_title('Density fluctuation growth')
    ax3a.invert_xaxis()
    ax3a.legend()

    ax3b.plot(rms_z, corr_all, 'ko-', markersize=5)
    ax3b.set_xlabel('Redshift z')
    ax3b.set_ylabel('Correlation r')
    ax3b.set_title('Point-by-point correlation')
    ax3b.invert_xaxis()
    ax3b.set_ylim(0, 1)

    fig3.suptitle('Growth & Correlation: GADGET-4 vs S-P (matched ICs)', fontsize=13)
    fig3.tight_layout()
    fig3.savefig(os.path.join(out_dir, 'summary_stats.png'), dpi=150)
    print(f"Saved {out_dir}/summary_stats.png")

    plt.show()


if __name__ == "__main__":
    main()
