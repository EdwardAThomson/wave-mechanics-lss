// Vertical moments and conserved quantities.
//
// The moments here are the 1D specialisation of the plan's §6 list. In a slab
// with no lateral structure the corrugation <z> and the mean vertical velocity
// <v_z> are single numbers rather than maps over (x, y), but they carry the
// same information: <z> and <v_z> in quadrature is a vertical oscillation, and
// their decay is the phase mixing that produces the spiral.
#pragma once

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "evolve.h"
#include "grid.h"
#include "poisson_z.h"
#include "units.h"

struct Moments {
    double Sigma = 0.0;      // Msun/kpc^2
    double zbar = 0.0;       // kpc
    double vbar = 0.0;       // km/s
    double z_rms = 0.0;      // kpc, about the mean
    double v_rms = 0.0;      // km/s, about the mean
    double ke = 0.0;         // kinetic energy per unit area
    double pe_self = 0.0;    // self-gravity energy per unit area (gauge free)
    double pe_ext = 0.0;
    double energy = 0.0;
    double edge_mass = 0.0;  // fraction of Sigma beyond 0.8 * L/2
};

// Fraction of the mass sitting in the outer fifth of the box. The kinetic step
// is applied on a periodic FFT grid, which is only legitimate while psi has
// decayed to ~0 at the edges; once a kicked sheet swings out far enough to
// wrap around, the isolated-boundary Poisson solve and the periodic kinetic
// step disagree and energy conservation collapses. This is the tell-tale, and
// it is a constraint on the box height set by the KICKED orbit, not by the
// equilibrium sheet.
inline double edge_mass_fraction(const Grid1D& g, const SlabState& st) {
    std::vector<double> rho;
    st.density(rho);
    double tot = 0.0, out = 0.0;
    for (int j = 0; j < g.N; ++j) {
        tot += rho[j];
        if (std::fabs(g.z(j)) > 0.8 * 0.5 * g.L) out += rho[j];
    }
    return (tot > 0.0) ? out / tot : 0.0;
}

inline Moments compute_moments(const Evolver& ev, const SlabState& st,
                               const std::vector<double>& Phi_ext) {
    const Grid1D& g = ev.grid();
    const int N = g.N;
    const double hb = ev.hbar();

    std::vector<double> rho;
    st.density(rho);

    std::vector<double> rho_v(N, 0.0);
    std::vector<Complex> dpsi;
    for (int s = 0; s < st.n_streams(); ++s) {
        ev.ddz(st.psi[s], dpsi);
        for (int j = 0; j < N; ++j) {
            rho_v[j] += hb * std::imag(std::conj(st.psi[s][j]) * dpsi[j]);
        }
    }
    // Use the evolver's own kinetic operator so that the reported energy is
    // the quantity the scheme actually conserves.
    const double ke = ev.kinetic_energy(st);

    Moments m;
    m.Sigma = g.integrate(rho);

    double sz = 0.0, sv = 0.0;
    for (int j = 0; j < N; ++j) {
        sz += rho[j] * g.z(j);
        sv += rho_v[j];
    }
    m.zbar = sz * g.dz / m.Sigma;
    m.vbar = sv * g.dz / m.Sigma;

    double s2 = 0.0;
    for (int j = 0; j < N; ++j) {
        const double dz = g.z(j) - m.zbar;
        s2 += rho[j] * dz * dz;
    }
    m.z_rms = std::sqrt(s2 * g.dz / m.Sigma);

    // <v^2> from the kinetic energy: KE = (1/2) int rho <v^2> dz, which for the
    // Madelung decomposition already includes the quantum-pressure share.
    const double v2 = 2.0 * ke / m.Sigma;
    m.v_rms = std::sqrt(std::max(0.0, v2 - m.vbar * m.vbar));

    // Gauge-independent self-energy: W = pi G int int rho rho |z - z'|
    std::vector<double> Phi_raw;
    poisson_isolated_1d(g, rho, Phi_raw, nullptr, /*gauge_zero=*/false);
    double w = 0.0, we = 0.0;
    for (int j = 0; j < N; ++j) {
        w += rho[j] * Phi_raw[j];
        we += rho[j] * Phi_ext[j];
    }
    m.ke = ke;
    m.pe_self = 0.5 * w * g.dz;
    m.pe_ext = we * g.dz;
    m.energy = m.ke + m.pe_self + m.pe_ext;
    m.edge_mass = edge_mass_fraction(g, st);
    return m;
}

// Density, mean velocity and potential on the grid, for plotting.
inline void write_profile(const std::string& path, const Evolver& ev,
                          const SlabState& st,
                          const std::vector<double>& Phi_ext,
                          double t) {
    const Grid1D& g = ev.grid();
    const int N = g.N;
    const double hb = ev.hbar();

    std::vector<double> rho;
    st.density(rho);
    std::vector<double> rho_v(N, 0.0);
    std::vector<Complex> dpsi;
    for (int s = 0; s < st.n_streams(); ++s) {
        ev.ddz(st.psi[s], dpsi);
        for (int j = 0; j < N; ++j) {
            rho_v[j] += hb * std::imag(std::conj(st.psi[s][j]) * dpsi[j]);
        }
    }
    std::vector<double> Phi;
    poisson_isolated_1d(g, rho, Phi);

    FILE* f = std::fopen(path.c_str(), "w");
    std::fprintf(f, "# t = %.8e  (%.4f Myr)\n", t, t * units::TIME_MYR);
    std::fprintf(f, "# z[kpc]  rho[Msun/kpc^3]  vz[km/s]  Phi_self  Phi_ext\n");
    for (int j = 0; j < N; ++j) {
        const double v = (rho[j] > 0.0) ? rho_v[j] / rho[j] : 0.0;
        std::fprintf(f, "%.8e %.8e %.8e %.8e %.8e\n",
                     g.z(j), rho[j], v, Phi[j], Phi_ext[j]);
    }
    std::fclose(f);
}
