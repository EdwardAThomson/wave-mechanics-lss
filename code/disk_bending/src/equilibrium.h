// Warm initial conditions for a self-gravitating vertical slab.
//
// ----------------------------------------------------------------------------
// Why eigenstates rather than sampled streams
// ----------------------------------------------------------------------------
// A single psi = sqrt(rho) exp(iS/hbar) is a cold stream. A disk with
// sigma_z = 20 km/s is not, and the warmth is not decorative: it sets the
// firehose boundary and supplies the Landau damping, so a cold run gets the
// wrong dispersion relation and shows no damping.
//
// The obvious fix, sampling streams from f(z, v_z) and summing |psi_j|^2
// incoherently, has a trap. A stream of constant v_j is not stationary in the
// sheet potential: it simply falls. The streams that *are* stationary are
// constant-energy tori, and those need two branches v = +/- sqrt(2(E - Phi))
// per energy with a density ~ 1/|v| that diverges at the turning point, where
// WKB fails. Building those cleanly is fiddly and it is the most likely way to
// lose a week (the plan's own §3 warning).
//
// The eigenstates of the self-consistent vertical Hamiltonian dodge all of it.
// |u_n|^2 is *exactly* stationary by construction, so validation test 4 passes
// to solver precision rather than being a thing to debug. Occupying them with
// w_n proportional to exp(-E_n / sigma_z^2) reproduces the isothermal
// distribution function to O(hbar_eff^2), which is the same order at which the
// whole Schrodinger-Poisson correspondence holds anyway. This is the standard
// construction in the wave-dark-matter halo literature (Lin, Chiueh et al.;
// Yavetz, Li & Hui 2021) and it transfers directly to a slab.
//
// The same eigenstate library serves both of the plan's IC options:
//   Option A (multi-stream)  : evolve each sqrt(w_n) u_n separately, sum |.|^2
//   Option B (single psi)    : psi = sum_n sqrt(w_n) exp(i theta_n) u_n
// so the two can be compared like for like on one basis. See README.md.
// ----------------------------------------------------------------------------
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "eigen1d.h"
#include "grid.h"
#include "poisson_z.h"
#include "units.h"

struct SheetConfig {
    double Sigma = 50.0 * units::MSUN_PC2;  // stellar surface density, Msun/kpc^2
    double sigma_z = 20.0;                  // target vertical dispersion, km/s
    double hbar = 0.3;                      // hbar_eff, kpc km/s

    // Rigid external vertical potential, standing in for the halo and any
    // component not being evolved:
    //
    //   Phi_ext(z) = nu_ext^2 z_ext^2 [ sqrt(1 + z^2/z_ext^2) - 1 ]
    //
    // harmonic with frequency nu_ext near the midplane, flattening to linear
    // above z_ext. Two separate reasons it is not optional:
    //
    //  1. Without any external term, an isolated self-gravitating slab given a
    //     bulk vertical kick just translates forever. Nothing restores the
    //     centre of mass, because self-gravity is internal.
    //
    //  2. Without ANHARMONICITY, there is still no phase spiral. In a purely
    //     harmonic external potential the centre-of-mass mode decouples
    //     exactly from every internal degree of freedom and oscillates
    //     undamped forever (Kohn's theorem), and every orbit shares one
    //     frequency so the distribution function rotates rigidly instead of
    //     winding. Set z_ext very large to recover that degenerate case and
    //     watch the spiral disappear; it is a useful check, not a bug.
    double nu_ext = 0.0;   // km/s/kpc, midplane frequency of the rigid term
    double z_ext = 1.0;    // kpc, height at which it flattens

    int n_states = 320;   // eigenstates to build; see the truncation warning
    int orth_window = 40;
    int scf_iters = 60;
    double scf_mix = 0.5;
    double scf_tol = 1e-8;  // on |dPhi| / sigma_z^2
    bool verbose = true;
};

struct WarmSheet {
    Grid1D grid;
    SheetConfig cfg;

    std::vector<double> Phi_self;  // self-gravity, gauge Phi_self(0) = 0
    std::vector<double> Phi_ext;   // rigid external
    std::vector<double> rho;       // Msun/kpc^3

    std::vector<double> E;                     // eigenvalues, (km/s)^2
    std::vector<double> w;                     // weights, sum w = Sigma
    std::vector<std::vector<double>> u;        // int |u|^2 dz = 1

    // diagnostics
    int n_significant = 0;      // bound states actually retained
    int n_requested = 0;        // states asked for
    double well_depth = 0.0;    // (V_edge - E_0) / sigma_z^2
    double energy_cut = 0.0;    // (E_last - E_0) / sigma_z^2, the truncation
    double tail_fraction = 0.0; // weight of the last state, relative to ground
    double sigma_measured = 0.0;
    double h_measured = 0.0;
    double rho0 = 0.0;
    double omega0 = 0.0;        // midplane vertical frequency, km/s/kpc
    double edge_potential = 0.0;
    double scf_residual = 0.0;
    int scf_used = 0;
};

// Effective vertical frequency including the rigid external term.
inline double sheet_omega0(const WarmSheet& s) { return s.omega0; }

// Drop states with E - E_0 above E_cap and renormalise the occupations back
// to Sigma. The 2D MatchedFD runs need this: the finite-difference kinetic
// band saturates at KE_max = 2 hbar^2 / dz^2, and the bound-state cut alone
// (V at the box edge) can sit above that, in which case the top thermal
// states are band-edge eigenvectors with checkerboard components that
// register as vertical spectral spill. Cap the library at a fraction of the
// band instead. The discarded thermal tail lowers sigma_z slightly;
// sigma_measured is recomputed so the deficit is visible, not hidden.
inline void truncate_sheet_energy(WarmSheet& s, double E_cap) {
    const int M0 = static_cast<int>(s.u.size());
    int M = 0;
    while (M < M0 && s.E[M] - s.E[0] < E_cap) ++M;
    if (M < 1) M = 1;
    if (M == M0) return;
    s.u.resize(M);
    s.E.resize(M);
    s.w.resize(M);
    double wsum = 0.0;
    for (double w : s.w) wsum += w;
    for (double& w : s.w) w *= s.cfg.Sigma / wsum;
    s.n_significant = M;
    s.energy_cut = (s.E[M - 1] - s.E[0]) / (s.cfg.sigma_z * s.cfg.sigma_z);
    s.tail_fraction = s.w[M - 1] / s.w[0];

    const Grid1D& g = s.grid;
    double ke = 0.0;
    for (int m = 0; m < M; ++m) {
        double pe = 0.0;
        for (int j = 0; j < g.N; ++j) {
            pe += (s.Phi_self[j] + s.Phi_ext[j]) * s.u[m][j] * s.u[m][j];
        }
        ke += s.w[m] * (s.E[m] - pe * g.dz);
    }
    s.sigma_measured = std::sqrt(2.0 * ke / s.cfg.Sigma);
}

inline WarmSheet build_warm_sheet(const Grid1D& g, const SheetConfig& cfg) {
    WarmSheet s;
    s.grid = g;
    s.cfg = cfg;

    const int N = g.N;
    const double hb = cfg.hbar;
    const double sig2 = cfg.sigma_z * cfg.sigma_z;

    // Rigid external potential
    s.Phi_ext.assign(N, 0.0);
    if (cfg.nu_ext != 0.0) {
        const double ze = cfg.z_ext;
        const double a = cfg.nu_ext * cfg.nu_ext * ze * ze;
        for (int j = 0; j < N; ++j) {
            const double z = g.z(j);
            s.Phi_ext[j] = a * (std::sqrt(1.0 + z * z / (ze * ze)) - 1.0);
        }
    }

    // Start from the analytic isothermal sheet; the SCF loop corrects it for
    // both the external term and the O(hbar^2) quantum pressure.
    IsothermalSheet analytic(cfg.Sigma, cfg.sigma_z);
    s.Phi_self.assign(N, 0.0);
    for (int j = 0; j < N; ++j) s.Phi_self[j] = analytic.Phi(g.z(j));

    std::vector<double> d(N), e(N - 1), V(N), Phi_new(N);
    const double kin_diag = hb * hb / (g.dz * g.dz);
    const double kin_off = -0.5 * hb * hb / (g.dz * g.dz);
    for (int j = 0; j < N - 1; ++j) e[j] = kin_off;

    const int M = std::min(cfg.n_states, N - 2);
    int M_bound = M;
    SymTridiagEigen eig;

    for (int it = 0; it < cfg.scf_iters; ++it) {
        for (int j = 0; j < N; ++j) {
            V[j] = s.Phi_self[j] + s.Phi_ext[j];
            d[j] = kin_diag + V[j];
        }

        eig = sym_tridiag_lowest(d, e, M, cfg.orth_window);

        // Keep only states genuinely bound by the well. An isolated slab's own
        // potential rises only linearly, 2 pi G Sigma |z|, so a box of height L
        // confines to a depth of about 2 pi G Sigma L / 2 and no more.
        // Eigenstates above that are box states: they turn around on the wall
        // rather than in the disk, and occupying them injects mass that
        // bounces off the boundary instead of orbiting. That shows up as a
        // gross non-stationarity that looks like a solver bug and is not one.
        const double v_edge = std::min(V[0], V[N - 1]);
        M_bound = 0;
        while (M_bound < M && eig.values[M_bound] < v_edge) ++M_bound;
        if (M_bound < 1) M_bound = 1;

        // Occupations: f(E) ~ exp(-E / sigma_z^2)
        std::vector<double> w(M_bound);
        double wsum = 0.0;
        for (int m = 0; m < M_bound; ++m) {
            w[m] = std::exp(-(eig.values[m] - eig.values[0]) / sig2);
            wsum += w[m];
        }
        for (int m = 0; m < M_bound; ++m) w[m] *= cfg.Sigma / wsum;

        // Density from the occupied states. Eigenvectors have 2-norm 1, so
        // dividing by dz converts sum u^2 = 1 into int |u|^2 dz = 1.
        std::vector<double> rho(N, 0.0);
        for (int m = 0; m < M_bound; ++m) {
            const auto& v = eig.vectors[m];
            const double wm = w[m] / g.dz;
            for (int j = 0; j < N; ++j) rho[j] += wm * v[j] * v[j];
        }

        poisson_isolated_1d(g, rho, Phi_new);

        double err = 0.0;
        for (int j = 0; j < N; ++j) {
            err = std::max(err, std::fabs(Phi_new[j] - s.Phi_self[j]));
        }
        err /= sig2;

        for (int j = 0; j < N; ++j) {
            s.Phi_self[j] = (1.0 - cfg.scf_mix) * s.Phi_self[j] +
                            cfg.scf_mix * Phi_new[j];
        }

        s.rho = rho;
        s.w = w;
        s.E.assign(eig.values.begin(), eig.values.begin() + M_bound);
        s.scf_residual = err;
        s.scf_used = it + 1;

        if (cfg.verbose && (it % 10 == 0 || err < cfg.scf_tol)) {
            std::printf("  scf %3d   |dPhi|/sigma^2 = %.3e\n", it, err);
        }
        if (err < cfg.scf_tol) break;
    }

    s.u.assign(M_bound, std::vector<double>(N));
    const double inv_sqrt_dz = 1.0 / std::sqrt(g.dz);
    for (int m = 0; m < M_bound; ++m) {
        for (int j = 0; j < N; ++j) s.u[m][j] = eig.vectors[m][j] * inv_sqrt_dz;
    }

    // ---- diagnostics ----
    s.rho0 = 0.0;
    for (int j = 0; j < N; ++j) {
        if (std::fabs(g.z(j)) < 0.5 * g.dz * 2.0) s.rho0 = std::max(s.rho0, s.rho[j]);
    }
    s.h_measured = cfg.Sigma / (4.0 * s.rho0);
    s.omega0 = std::sqrt(4.0 * units::PI * units::G * s.rho0 +
                         cfg.nu_ext * cfg.nu_ext);

    // Mass-weighted <v_z^2> from the kinetic energy of the occupied states.
    // For a real eigenstate <v_z> = 0 but <v_z^2> is not; in the semiclassical
    // limit this converges to sigma_z^2 with an O(hbar^2) excess.
    // Take it from the Hamiltonian rather than from a finite-difference
    // gradient of u: KE_n = E_n - int V |u_n|^2 is exact for the discrete
    // operator, whereas a centred difference underestimates |u'|^2 by
    // sin^2(k dz)/(k dz)^2, which for the highest occupied states is a
    // double-digit error and shows up as a spurious sigma_z deficit.
    double ke = 0.0;
    for (int m = 0; m < M_bound; ++m) {
        const auto& un = s.u[m];
        double pe = 0.0;
        for (int j = 0; j < N; ++j) pe += V[j] * un[j] * un[j];
        pe *= g.dz;
        ke += s.w[m] * (s.E[m] - pe);
    }
    s.sigma_measured = std::sqrt(2.0 * ke / cfg.Sigma);

    // Weight bookkeeping. The occupations are renormalised to sum to Sigma, so
    // a cumulative-mass criterion trivially returns M and says nothing. What
    // matters is where the library was truncated in energy: everything above
    // E_trunc is missing, and since those states carry the high velocities,
    // truncating too low shows up directly as a sigma_z deficit.
    s.energy_cut = (s.E[M_bound - 1] - s.E[0]) / sig2;
    s.n_significant = M_bound;
    s.n_requested = M;
    s.tail_fraction = s.w[M_bound - 1] / s.w[0];

    s.edge_potential = std::min(s.Phi_self[0] + s.Phi_ext[0],
                                s.Phi_self[N - 1] + s.Phi_ext[N - 1]);
    s.well_depth = (s.edge_potential - s.E[0]) / sig2;

    if (cfg.verbose) {
        std::printf("\n  warm sheet built with %d of %d requested states, "
                    "truncated at E = %.1f sigma_z^2\n",
                    M_bound, M, s.energy_cut);
        std::printf("  well depth       = %.1f sigma_z^2 "
                    "(V at box edge, above the ground state)\n", s.well_depth);
        std::printf("  Sigma            = %.4e Msun/kpc^2\n", cfg.Sigma);
        std::printf("  sigma_z target   = %.4f km/s\n", cfg.sigma_z);
        std::printf("  sigma_z measured = %.4f km/s   (excess %+.2f%%)\n",
                    s.sigma_measured,
                    100.0 * (s.sigma_measured / cfg.sigma_z - 1.0));
        std::printf("  h (analytic)     = %.5f kpc\n", analytic.h);
        std::printf("  h (measured)     = %.5f kpc\n", s.h_measured);
        std::printf("  rho0             = %.4e Msun/kpc^3\n", s.rho0);
        std::printf("  omega0           = %.3f km/s/kpc   (T_vert = %.2f Myr)\n",
                    s.omega0, 2.0 * units::PI / s.omega0 * units::TIME_MYR);
        std::printf("  hbar_eff         = %.4f kpc km/s\n", hb);
        std::printf("  hbar*omega0      = %.4f (km/s)^2   = %.4f sigma_z^2\n",
                    hb * s.omega0, hb * s.omega0 / sig2);
        std::printf("  lambda_dB        = %.5f kpc  = %.3f h\n",
                    2.0 * units::PI * hb / cfg.sigma_z,
                    2.0 * units::PI * hb / cfg.sigma_z / s.h_measured);
        std::printf("  top-state weight = %.3e of ground state\n", s.tail_fraction);
        std::printf("  scf residual     = %.3e after %d iters\n",
                    s.scf_residual, s.scf_used);
        if (M_bound < M) {
            std::printf("\n  NOTE: %d of the %d requested states were unbound "
                        "and have been dropped.\n"
                        "  An isolated slab's own potential rises only as "
                        "2 pi G Sigma |z|, so the box\n"
                        "  confines to about %.1f sigma_z^2 and no deeper. To "
                        "go higher, raise L or\n"
                        "  add a rigid external term with nu_ext > 0 (a halo). "
                        "Occupying unbound\n"
                        "  states would put mass into wall-bouncing orbits and "
                        "wreck stationarity.\n",
                        M - M_bound, M, s.well_depth);
        }
        if (s.energy_cut < 8.0) {
            std::printf("\n  WARNING: the state library reaches only "
                        "%.1f sigma_z^2.\n"
                        "  The missing high-energy states carry the fast "
                        "particles, so sigma_z will\n"
                        "  come out low (currently %+.1f%%). Deepen the well, "
                        "then raise n_states\n"
                        "  until sigma_measured stops moving.\n",
                        s.energy_cut,
                        100.0 * (s.sigma_measured / cfg.sigma_z - 1.0));
        }
    }
    return s;
}
