// Column moments for the 2D slab: the plan's §6 diagnostics, now with a real
// x axis to carry them.
//
//   Sigma(x)    = Integral rho dz                     surface density
//   <z>(x)      = Integral z rho dz / Sigma           the corrugation
//   <v_z>(x)    = Integral rho v_z dz / Sigma         mean vertical velocity
//
// The primary result in §6 is the phase relation between the last two. Here
// that is read off directly from their complex Fourier amplitudes at the
// excited wavenumber: a quarter-cycle (pi/2) offset is a propagating bending
// wave, 0 or pi a standing one.
#pragma once

#include <cmath>
#include <complex>
#include <cstdio>
#include <string>
#include <vector>

#include "evolve2d.h"
#include "grid2d.h"
#include "poisson_xz.h"
#include "units.h"

struct ColumnMoments {
    std::vector<double> Sigma;  // [Nx]
    std::vector<double> zbar;   // [Nx]
    std::vector<double> vbar;   // [Nx]
    double Sigma_total = 0.0;   // per unit y length
    double sigma_z_rms = 0.0;   // mass-weighted, about the local <v_z>
    double sigma_x_rms = 0.0;   // in-plane, from the x kinetic energy
    double ke = 0.0, pe = 0.0, energy = 0.0;
    double edge_mass = 0.0;
};

// Light version: the three column profiles only, one vertical derivative per
// stream and no Poisson solve. This is what the time series calls at every
// sample, so it must not carry the energy machinery.
inline ColumnMoments compute_columns_2d(const Evolver2D& ev,
                                        const SlabState2D& st) {
    const Grid2D& g = ev.grid();
    const double hb = ev.hbar();

    std::vector<double> rho;
    st.density(rho);
    std::vector<double> rho_v(g.size(), 0.0);
    std::vector<Complex> dpsi;
    for (int s = 0; s < st.n_streams(); ++s) {
        ev.ddz(st.psi[s], dpsi);
        for (int i = 0; i < g.size(); ++i) {
            rho_v[i] += hb * std::imag(std::conj(st.psi[s][i]) * dpsi[i]);
        }
    }

    ColumnMoments m;
    m.Sigma.assign(g.Nx, 0.0);
    m.zbar.assign(g.Nx, 0.0);
    m.vbar.assign(g.Nx, 0.0);
    for (int i = 0; i < g.Nx; ++i) {
        double s0 = 0.0, s1 = 0.0, sv = 0.0;
        for (int j = 0; j < g.Nz; ++j) {
            const double r = rho[g.idx(i, j)];
            s0 += r;
            s1 += r * g.z(j);
            sv += rho_v[g.idx(i, j)];
        }
        m.Sigma[i] = s0 * g.dz;
        m.zbar[i] = (s0 > 0.0) ? s1 / s0 : 0.0;
        m.vbar[i] = (s0 > 0.0) ? sv / s0 : 0.0;
    }
    for (int i = 0; i < g.Nx; ++i) m.Sigma_total += m.Sigma[i] * g.dx;
    return m;
}

inline ColumnMoments compute_moments_2d(const Evolver2D& ev,
                                        const SlabState2D& st,
                                        PoissonXZ& poisson) {
    const Grid2D& g = ev.grid();
    const double hb = ev.hbar();

    std::vector<double> rho;
    st.density(rho);

    std::vector<double> rho_v(g.size(), 0.0);
    std::vector<Complex> dpsi;
    for (int s = 0; s < st.n_streams(); ++s) {
        ev.ddz(st.psi[s], dpsi);
        for (int i = 0; i < g.size(); ++i) {
            rho_v[i] += hb * std::imag(std::conj(st.psi[s][i]) * dpsi[i]);
        }
    }

    ColumnMoments m;
    m.Sigma.assign(g.Nx, 0.0);
    m.zbar.assign(g.Nx, 0.0);
    m.vbar.assign(g.Nx, 0.0);

    double edge = 0.0, tot = 0.0;
    for (int i = 0; i < g.Nx; ++i) {
        double s0 = 0.0, s1 = 0.0, sv = 0.0;
        for (int j = 0; j < g.Nz; ++j) {
            const double r = rho[g.idx(i, j)];
            s0 += r;
            s1 += r * g.z(j);
            sv += rho_v[g.idx(i, j)];
            tot += r;
            if (std::fabs(g.z(j)) > 0.8 * 0.5 * g.Lz) edge += r;
        }
        m.Sigma[i] = s0 * g.dz;
        m.zbar[i] = (s0 > 0.0) ? s1 / s0 : 0.0;
        m.vbar[i] = (s0 > 0.0) ? sv / s0 : 0.0;
    }
    m.Sigma_total = 0.0;
    for (int i = 0; i < g.Nx; ++i) m.Sigma_total += m.Sigma[i] * g.dx;
    m.edge_mass = (tot > 0.0) ? edge / tot : 0.0;

    // Energies. The kinetic term uses the evolver's own operator so that the
    // reported energy is the quantity the scheme conserves.
    m.ke = ev.kinetic_energy(st);
    std::vector<double> Phi;
    // Ungauged: the Phi(0) = 0 convention adds a profile-dependent constant
    // whose drift would masquerade as an energy conservation failure.
    poisson.solve(rho, Phi, /*gauge_zero=*/false);
    double w = 0.0;
    for (int i = 0; i < g.size(); ++i) w += rho[i] * Phi[i];
    m.pe = 0.5 * w * g.dx * g.dz;
    m.energy = m.ke + m.pe;

    // Split the kinetic energy into its x and z shares to track heating in
    // each separately: in-plane dispersion is what drives firehose, vertical
    // dispersion is what damps.
    double kx2 = 0.0, kz2 = 0.0;
    {
        std::vector<Complex> dz_psi;
        for (int s = 0; s < st.n_streams(); ++s) {
            ev.ddz(st.psi[s], dz_psi);
            for (int i = 0; i < g.size(); ++i) kz2 += std::norm(dz_psi[i]);
        }
        kz2 *= 0.5 * hb * hb * g.dx * g.dz;
        kx2 = m.ke - kz2;
    }
    m.sigma_z_rms = std::sqrt(std::max(0.0, 2.0 * kz2 / m.Sigma_total));
    m.sigma_x_rms = std::sqrt(std::max(0.0, 2.0 * kx2 / m.Sigma_total));
    return m;
}

// Complex Fourier amplitude of a real field at mode index n (so that
// f(x) ~ 2 Re[a exp(i k_n x)]).
inline Complex mode_amplitude(const Grid2D& g, const std::vector<double>& f,
                              int n) {
    Complex acc(0.0, 0.0);
    for (int i = 0; i < g.Nx; ++i) {
        const double ph = -g.kx(n) * g.x(i);
        acc += f[i] * Complex(std::cos(ph), std::sin(ph));
    }
    return acc / static_cast<double>(g.Nx);
}

inline void write_columns(const std::string& path, const Grid2D& g,
                          const ColumnMoments& m, double t) {
    FILE* f = std::fopen(path.c_str(), "w");
    std::fprintf(f, "# t = %.8e  (%.4f Myr)\n", t, t * units::TIME_MYR);
    std::fprintf(f, "# x[kpc]  Sigma[Msun/kpc^2]  <z>[kpc]  <vz>[km/s]\n");
    for (int i = 0; i < g.Nx; ++i) {
        std::fprintf(f, "%.8e %.8e %.8e %.8e\n", g.x(i), m.Sigma[i], m.zbar[i],
                     m.vbar[i]);
    }
    std::fclose(f);
}
