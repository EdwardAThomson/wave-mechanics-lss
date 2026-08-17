#!/usr/bin/env python3
"""Plot the Husimi phase-space maps f(z, v_z) written by sheet_1d.

    python3 plot_phase_spiral.py output/kick_bulk [--absolute]

The vertical axis is v_z / omega_z, in kpc, so that an unperturbed orbit is a
circle and a wound-up spiral looks like a spiral rather than a sheared blob.
By default each map is divided by its own azimuthal average at fixed orbital
radius, which is how the Antoja snail is normally shown: the smooth
equilibrium divides out and the winding is what is left. --absolute shows the
raw distribution function instead.
"""
import sys
import os
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_info(outdir):
    info = {}
    path = os.path.join(outdir, "run_info.dat")
    if os.path.exists(path):
        for line in open(path):
            parts = line.split()
            if len(parts) == 2:
                try:
                    info[parts[0]] = float(parts[1])
                except ValueError:
                    info[parts[0]] = parts[1]
    return info


def load_husimi(path):
    with open(path) as fh:
        header = [fh.readline() for _ in range(3)]
    t_myr = float(header[0].split("(")[1].split()[0])
    nz = int(header[1].split("nz =")[1].split()[0])
    nv = int(header[1].split("nv =")[1].split()[0])
    d = np.loadtxt(path)
    z = d[:, 0].reshape(nz, nv)[:, 0]
    v = d[:, 1].reshape(nz, nv)[0, :]
    f = d[:, 2].reshape(nz, nv)
    return t_myr, z, v, f


def azimuthal_contrast(z, y, f, nb=90, mass_frac=0.99, z0=0.0, y0=0.0):
    """f / <f>_theta - 1, binned in orbital radius r = hypot(z, y).

    Restricted to the radius enclosing `mass_frac` of the distribution. Beyond
    that the distribution function is exponentially small and the relative
    contrast is dominated by nothing at all, which otherwise swamps the plot
    with meaningless colour in the empty corners.
    """
    Z, Y = np.meshgrid(z - z0, y - y0, indexing="ij")
    r = np.hypot(Z, Y)
    rmax = min(z.max() - abs(z0), y.max() - abs(y0))
    idx = np.clip((r / rmax * nb).astype(int), 0, nb - 1)

    # Radius enclosing mass_frac of the mass inside rmax.
    inside = r <= rmax
    order = np.argsort(r[inside])
    cum = np.cumsum(f[inside][order])
    if cum[-1] <= 0:
        return np.zeros_like(f), rmax
    r_cut = r[inside][order][np.searchsorted(cum, mass_frac * cum[-1])]
    r_cut = max(r_cut, 0.2 * rmax)

    out = np.zeros_like(f)
    for b in range(nb):
        m = idx == b
        if m.sum() < 8:
            continue
        bg = f[m].mean()
        if bg <= 0:
            continue
        out[m] = f[m] / bg - 1.0
    out[r > r_cut] = 0.0
    return out, r_cut


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else "output/kick_bulk"
    absolute = "--absolute" in sys.argv
    files = sorted(f for f in glob.glob(os.path.join(outdir, "husimi_*.dat"))
                   if "equilibrium" not in os.path.basename(f))
    if not files:
        sys.exit(f"no husimi_*.dat in {outdir}")

    info = read_info(outdir)
    omega = info.get("omega0", 75.0)
    comoving = "--fixed-frame" not in sys.argv

    # Centroid history, so the contrast can be taken about the moving centre of
    # mass. In a self-gravitating slab the bulk slosh is a coherent collective
    # mode that does not wind, and left in place it dominates the m = 1 signal
    # and hides the spiral entirely. The spiral lives in the residual.
    cen_t = cen_z = cen_v = None
    mom = os.path.join(outdir, "moments.dat")
    if comoving and os.path.exists(mom):
        m = np.loadtxt(mom)
        cen_t, cen_z, cen_v = m[:, 1], m[:, 3], m[:, 4]

    # Unperturbed equilibrium map, dumped by sheet_1d before the kick. Dividing
    # by it removes the equilibrium exactly, with no assumption that the
    # distribution is circular in (z, v/omega) -- which it is not, because the
    # potential is anharmonic.
    f_eq = None
    eq_path = os.path.join(outdir, "husimi_equilibrium.dat")
    if os.path.exists(eq_path):
        _, _, _, f_eq = load_husimi(eq_path)

    n = len(files)
    nrow = 1 if n <= 5 else int(np.ceil(n / np.ceil(n / 2.0)))
    ncol = int(np.ceil(n / nrow))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.1 * ncol, 3.3 * nrow),
                             squeeze=False)
    for k, path in enumerate(files):
        ax = axes[k // ncol][k % ncol]
        t, z, v, f = load_husimi(path)
        y = v / omega                      # kpc, so orbits are circles
        if absolute:
            img = np.log10(np.maximum(f, f.max() * 1e-4))
            kw = dict(cmap="magma")
            lim_r = min(z.max(), y.max())
            z0 = y0 = 0.0
        else:
            z0 = y0 = 0.0
            if cen_t is not None:
                z0 = float(np.interp(t, cen_t, cen_z))
                y0 = float(np.interp(t, cen_t, cen_v)) / omega
            if f_eq is not None and f_eq.shape == f.shape:
                floor = f_eq.max() * 3e-3
                img = np.where(f_eq > floor, f / np.maximum(f_eq, floor) - 1.0,
                               0.0)
                lim_r = min(z.max(), y.max())
                z0 = y0 = 0.0
            else:
                img, lim_r = azimuthal_contrast(z, y, f, z0=z0, y0=y0)
            nz_vals = np.abs(img[img != 0])
            lim = max(np.percentile(nz_vals, 99) if nz_vals.size else 1e-6, 1e-9)
            kw = dict(cmap="RdBu_r", vmin=-lim, vmax=lim)
        ax.imshow(img.T, origin="lower", extent=[z[0], z[-1], y[0], y[-1]],
                  aspect="equal", **kw)
        cx = z0 if not absolute else 0.0
        cy = y0 if not absolute else 0.0
        ax.set_xlim(cx - lim_r, cx + lim_r)
        ax.set_ylim(cy - lim_r, cy + lim_r)
        ax.set_title(f"t = {t:.0f} Myr", fontsize=9)
        ax.set_xlabel("z [kpc]", fontsize=8)
        if k % ncol == 0:
            ax.set_ylabel(r"$v_z/\omega_z$ [kpc]", fontsize=8)
        ax.tick_params(labelsize=7)
    for k in range(n, nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")

    kind = "absolute" if absolute else "contrast"
    sub = (f"$\\omega_z$ = {omega:.0f} km/s/kpc, "
           f"$T_z$ = {info.get('T_vert_Myr', 0):.0f} Myr, "
           f"{int(info.get('n_streams', 0))} streams, "
           f"$\\hbar_{{\\rm eff}}$ = {info.get('hbar', 0):.2f} kpc km/s")
    fig.suptitle(f"Vertical phase spiral, Husimi $f(z, v_z)$ ({kind})\n{sub}",
                 fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.93], h_pad=2.2)
    out = os.path.join(outdir, f"phase_spiral_{kind}.png")
    fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
