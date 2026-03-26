#!/usr/bin/env python3
"""Extract density and velocity potential from a GADGET-4 IC snapshot.

Reads particle positions and velocities, CIC-deposits onto an N^3 grid,
computes the velocity potential via div(v) in Fourier space, and saves
both fields as binary files for the S-P code to read.

Usage: python3 extract_gadget_ics.py <snapshot.hdf5> <output_dir>
"""

import numpy as np
import h5py
import sys
import os

def tsc_weight(d):
    """Triangular Shaped Cloud weight for offset d from cell centre.
    d should be in [-1.5, 1.5]. Returns weight for 3 cells: -1, 0, +1."""
    ad = np.abs(d)
    w = np.zeros_like(d)
    # |d| < 0.5: w = 0.75 - d^2
    mask1 = ad < 0.5
    w[mask1] = 0.75 - d[mask1]**2
    # 0.5 <= |d| < 1.5: w = 0.5 * (1.5 - |d|)^2
    mask2 = (ad >= 0.5) & (ad < 1.5)
    w[mask2] = 0.5 * (1.5 - ad[mask2])**2
    return w


def tsc_deposit(pos_norm, N):
    """TSC (Triangular Shaped Cloud) mass deposit onto an N^3 grid.

    TSC distributes each particle's mass across 27 cells (3^3) using
    quadratic weights. Smoother than CIC (8 cells, linear weights),
    giving a better representation of the continuous density field.

    pos_norm: (Npart, 3) positions normalised to [0, 1)
    Returns: (N, N, N) density field
    """
    dx = 1.0 / N
    result = np.zeros((N, N, N))

    # Grid position (cell-centred)
    gx = pos_norm[:, 0] / dx
    gy = pos_norm[:, 1] / dx
    gz = pos_norm[:, 2] / dx

    # Nearest cell
    ix0 = np.round(gx).astype(int)
    iy0 = np.round(gy).astype(int)
    iz0 = np.round(gz).astype(int)

    # Offset from cell centre
    fx = gx - ix0
    fy = gy - iy0
    fz = gz - iz0

    # Deposit to 27 surrounding cells
    for dix in [-1, 0, 1]:
        wx = tsc_weight(fx - dix)
        iix = (ix0 + dix) % N
        for diy in [-1, 0, 1]:
            wy = tsc_weight(fy - diy)
            iiy = (iy0 + diy) % N
            for diz in [-1, 0, 1]:
                wz = tsc_weight(fz - diz)
                iiz = (iz0 + diz) % N
                np.add.at(result, (iix, iiy, iiz), wx * wy * wz)

    return result


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <snapshot.hdf5> <output_dir>")
        sys.exit(1)

    snap_file = sys.argv[1]
    out_dir = sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)

    # Read GADGET snapshot (positions only — velocity derived from density)
    with h5py.File(snap_file, 'r') as f:
        pos = f['PartType1/Coordinates'][:].astype(np.float64)
        a = float(f['Header'].attrs['Time'])
        boxsize = float(f['Header'].attrs['BoxSize'])
        h = float(f['Parameters'].attrs['HubbleParam'])
        omega_m = float(f['Parameters'].attrs['Omega0'])

    npart = len(pos)
    N = round(npart ** (1.0 / 3.0))
    assert N * N * N == npart, f"Expected cubic particle count, got {npart}"

    z = 1.0 / a - 1.0
    print(f"Snapshot: a={a:.6f}, z={z:.1f}, N={N}^3={npart}")
    print(f"BoxSize={boxsize} Mpc/h, Omega_m={omega_m}, h={h}")

    # Normalise positions to [0, 1)
    pos_norm = pos / boxsize

    # TSC deposit: density
    print("\nTSC depositing density...")
    rho = tsc_deposit(pos_norm, N)
    rho_mean = np.mean(rho)
    delta = rho / rho_mean - 1.0
    print(f"delta: rms={np.std(delta):.6f}, range=[{delta.min():.4f}, {delta.max():.4f}]")

    # Compute velocity potential from density via Poisson equation
    # In the linear regime (z=50): v = -grad(phi_v), nabla^2(phi_v) = -delta
    # So phi_v_hat = -delta_hat / k^2
    print("Computing velocity potential from density (Poisson)...")
    delta_hat = np.fft.fftn(delta)

    kx = 2 * np.pi * np.fft.fftfreq(N, d=1.0 / N)  # box units: k = 2*pi*n
    KX, KY, KZ = np.meshgrid(kx, kx, kx, indexing='ij')
    K2 = KX**2 + KY**2 + KZ**2

    phi_v_hat = np.zeros_like(delta_hat)
    phi_v_hat[K2 > 0] = -delta_hat[K2 > 0] / K2[K2 > 0]
    phi_v = np.real(np.fft.ifftn(phi_v_hat))

    print(f"phi_v: rms={np.std(phi_v):.6e}, range=[{phi_v.min():.6e}, {phi_v.max():.6e}]")

    # Displacement field (for phase gradient check)
    # s = -grad(phi_v), so s_hat_j = -ik_j * phi_v_hat = ik_j * delta_hat / k^2
    sx_hat = np.zeros_like(delta_hat)
    sy_hat = np.zeros_like(delta_hat)
    sz_hat = np.zeros_like(delta_hat)
    sx_hat[K2 > 0] = 1j * KX[K2 > 0] * delta_hat[K2 > 0] / K2[K2 > 0]
    sy_hat[K2 > 0] = 1j * KY[K2 > 0] * delta_hat[K2 > 0] / K2[K2 > 0]
    sz_hat[K2 > 0] = 1j * KZ[K2 > 0] * delta_hat[K2 > 0] / K2[K2 > 0]
    sx = np.real(np.fft.ifftn(sx_hat))
    sy = np.real(np.fft.ifftn(sy_hat))
    sz = np.real(np.fft.ifftn(sz_hat))
    s_max = np.max(np.sqrt(sx**2 + sy**2 + sz**2))
    print(f"s_max = {s_max:.6f} (box units)")

    # Save as binary files
    # Format: N (int32), a (float64), then N^3 float64 values (row-major, z fastest)
    def save_field(filename, field, a_val):
        with open(filename, 'wb') as f:
            f.write(np.int32(N).tobytes())
            f.write(np.float64(a_val).tobytes())
            f.write(field.astype(np.float64).tobytes())
        print(f"  Saved {filename} ({os.path.getsize(filename) / 1e6:.1f} MB)")

    save_field(os.path.join(out_dir, "delta_gadget.bin"), delta, a)
    save_field(os.path.join(out_dir, "phi_v_gadget.bin"), phi_v, a)

    # Print phase budget for various nu values
    k_nyq = np.pi * N
    print(f"\nPhase budget (s_max={s_max:.6f}, k_Nyq={k_nyq:.1f}):")
    for nu in [3e-4, 5e-4, 1e-3, 3e-3, 1e-2]:
        pct = s_max / nu / k_nyq * 100
        print(f"  nu={nu:.0e}: {pct:.1f}% of Nyquist")

    print(f"\nDone. Files saved to {out_dir}/")


if __name__ == "__main__":
    main()
