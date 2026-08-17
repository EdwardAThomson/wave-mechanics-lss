// Rotating-frame machinery for the 2D (x, z) slab: Coriolis as a Landau-gauge
// vector potential, epicyclic Landau levels as the warm in-plane basis.
//
// ----------------------------------------------------------------------------
// The formulation (output/rotation_spike.md has the feasibility numbers)
// ----------------------------------------------------------------------------
// In a frame rotating rigidly at Omega, Coriolis enters the Schrodinger
// operator as a vector potential with uniform "field" kappa = 2 Omega. In
// Landau gauge A = (0, kappa x, 0), and with nothing depending on y (a bending
// wave with purely radial wavevector), p_y = hbar k_y is conserved and the
// problem separates into independent sectors, one per k_y:
//
//   H_ky = -(hbar^2/2)(d2/dx2 + d2/dz2) + (hbar k_y - kappa x)^2 / 2 + Phi(x,z)
//
// The Coriolis force is a LOCAL extra potential: a harmonic trap of frequency
// kappa about the guiding centre x_g = hbar k_y / kappa. No spectral
// cross-terms, no new operator structure. The multi-stream machinery carries
// it directly, one stream per guiding centre, with the trap applied in the
// ordinary potential step (see Evolver2D::set_traps).
//
// For a rigidly rotating disk the background radial gravity balances the
// centrifugal force at every radius, so the local box keeps ONLY the trap:
// there is no residual tidal term. (Shear would add one; that is Stage 2.)
//
// ----------------------------------------------------------------------------
// Warm in-plane initial conditions: thermally occupied epicyclic levels
// ----------------------------------------------------------------------------
// The eigenstates of the per-stream x-Hamiltonian are harmonic oscillator
// (Landau) levels phi_n(x - x_g) with E_n = hbar kappa (n + 1/2). Each level
// carries <v_x^2> = <v_y^2> = E_n (virial split of the 1D oscillator, with
// v_y = hbar k_y - kappa x riding the same coordinate), so thermal occupation
// w_n ~ exp(-E_n / sigma_x^2) reproduces an ISOTROPIC in-plane Schwarzschild
// distribution: sigma_y / sigma_x = kappa / 2 Omega = 1 for rigid rotation,
// which is the physically correct axis ratio for this box. This is the exact
// in-plane analogue of the vertical eigenstate construction in equilibrium.h.
//
//   psi_g(x, z) = sum_{m,n} sqrt(w^z_m w^epi_n) e^{i theta_{g,m,n}}
//                 u_m(z) phi_n(x - x_g)
//
// Mass bookkeeping is exact: the phi_n are orthonormal, so each stream's mass
// is Sigma Lx / n_streams independent of the random phases, and the expected
// density is rho_z(z) uniform in x once guiding centres are packed finer than
// the thermal envelope a_epi = sigma_x / kappa (a sum of Gaussians of width
// a_epi at spacing d has relative ripple ~ 2 exp(-2 pi^2 a_epi^2 / d^2), which
// at d = a_epi / 4 is zero to machine precision).
//
// STATIONARITY. phi_n are continuum eigenstates, not eigenvectors of the
// discrete operator. With the SPECTRAL in-plane kinetic step the only error is
// aliasing of the momentum tail beyond the grid ceiling, and the level
// truncation below keeps every occupied level's classical turnaround a safe
// margin under that ceiling, so the tail is exp(-large) and the stack is
// stationary for practical purposes. (Matched-FD in x would make it exact but
// distorts velocities by sin(k dx)/(k dx), double digits at sigma_x = 40, so
// spectral-in-x is the right choice; MatchedFD remains right for z.)
//
// TRUNCATION. The grid represents in-plane velocities only up to
// v = hbar pi / dx (the repo's headline constraint), so levels are truncated
// at E_max = (margin * v_ceiling)^2 / 2. The discarded Maxwellian tail weight
// is reported; at the production parameters it is a few per cent, and the
// realised sigma_x (also reported) is correspondingly a little low. This is
// the same truncation the vertical library already makes, applied in x.
//
// ----------------------------------------------------------------------------
// Periodic wrapping instead of an edge buffer
// ----------------------------------------------------------------------------
// The trap (hbar k_y - kappa x)^2 / 2 is not periodic in x. The spike proposed
// an edge buffer as mitigation; wrapping is strictly better and costs nothing:
// give every stream the trap 0.5 kappa^2 wrap(x - x_g)^2 with wrap() the
// minimal periodic distance. The trap then has a cusp at the stream's antipode
// where its amplitude is ~exp(-(Lx / 2 a_epi)^2 / 2) (about e^-39 at
// production parameters), the box is exactly translation invariant modulo the
// guiding-centre spacing, and no part of it needs excluding from fits. The
// rigorous theta-function Landau states remain available if this residual
// ever matters; a warning below fires when the box is too small for that
// argument to hold.
#pragma once

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "equilibrium.h"
#include "evolve2d.h"
#include "grid2d.h"
#include "slab_ic.h"  // SplitMix64
#include "units.h"

#ifdef _OPENMP
#include <omp.h>
#endif

struct RotConfig {
    double kappa = 44.2;    // epicyclic frequency, km/s/kpc; rigid box: 2 Omega
    double sigma_x = 40.0;  // target in-plane dispersion, km/s
    // Guiding centres. 0 = automatic at spacing a_epi / 4, which makes the
    // density ripple exactly zero in practice; fewer streams only raise the
    // interference speckle (at 1/sqrt(n)), not the systematic ripple, until
    // the spacing approaches a_epi itself.
    int n_streams = 0;
    // Levels are kept while their classical turnaround velocity is below
    // v_margin * (hbar pi / dx). The discarded tail is reported.
    double v_margin = 0.9;
    uint64_t seed = 20260817ULL;
    bool verbose = true;
};

struct RotatingSlab {
    SlabState2D st;
    std::vector<double> xg;      // guiding centre per stream, kpc
    double kappa = 0.0;
    int n_levels = 0;            // epicyclic levels occupied per stream
    double sigma_x_realised = 0.0;  // sqrt(sum w_n E_n) of the truncated stack
    double tail_dropped = 0.0;      // Maxwellian weight beyond the level cut
    double a_epi = 0.0;             // thermal epicyclic envelope, kpc
    double ell = 0.0;               // oscillator length sqrt(hbar/kappa), kpc
};

inline RotatingSlab make_warm_rotating(const WarmSheet& s, const Grid2D& g,
                                       const RotConfig& rc) {
    RotatingSlab r;
    r.kappa = rc.kappa;
    r.st.grid = g;
    r.st.hbar = s.cfg.hbar;

    const double hb = s.cfg.hbar;
    const double kap = rc.kappa;
    const double sig2 = rc.sigma_x * rc.sigma_x;
    r.ell = std::sqrt(hb / kap);

    // ---- level ladder, truncated at the velocity ceiling ----
    const double v_ceil = hb * units::PI / g.dx;
    const double v_trunc = rc.v_margin * v_ceil;
    const double E_max = 0.5 * v_trunc * v_trunc;
    const double quantum = hb * kap;  // level spacing, (km/s)^2
    int n_levels = static_cast<int>(std::floor(E_max / quantum + 0.5));
    if (n_levels < 1) n_levels = 1;
    r.n_levels = n_levels;

    std::vector<double> w(n_levels), E(n_levels);
    double wsum = 0.0;
    for (int n = 0; n < n_levels; ++n) {
        E[n] = quantum * (n + 0.5);
        w[n] = std::exp(-E[n] / sig2);
        wsum += w[n];
    }
    for (int n = 0; n < n_levels; ++n) w[n] /= wsum;
    // Geometric tail of the un-truncated ladder beyond the cut.
    const double ratio = std::exp(-quantum / sig2);
    r.tail_dropped = std::pow(ratio, n_levels);

    double Ebar = 0.0;
    for (int n = 0; n < n_levels; ++n) Ebar += w[n] * E[n];
    r.sigma_x_realised = std::sqrt(Ebar);
    r.a_epi = r.sigma_x_realised / kap;

    // ---- guiding centres ----
    int n_g = rc.n_streams;
    if (n_g <= 0) n_g = std::max(1, static_cast<int>(
                                        std::ceil(4.0 * g.Lx / r.a_epi)));
    r.xg.resize(n_g);
    for (int t = 0; t < n_g; ++t) r.xg[t] = (t + 0.5) * g.Lx / n_g;

    // Support half-width: top turning point plus a generous quantum tail.
    const double xi_max = std::sqrt(2.0 * n_levels) + 10.0;
    const double support = xi_max * r.ell;

    const int M = static_cast<int>(s.u.size());

    if (rc.verbose) {
        const double Q = kap * rc.sigma_x /
                         (3.36 * units::G * s.cfg.Sigma);
        std::printf("  rotating slab: kappa = %.2f km/s/kpc (Omega_box = %.2f), "
                    "Q = %.2f\n", kap, 0.5 * kap, Q);
        std::printf("  ell = %.4f kpc, a_epi = %.4f kpc, levels = %d "
                    "(E cut %.1f sigma_x^2, tail %.2e dropped)\n",
                    r.ell, r.a_epi, n_levels, E[n_levels - 1] / sig2,
                    r.tail_dropped);
        std::printf("  v ceiling = %.1f km/s, level cut at %.1f km/s\n",
                    v_ceil, v_trunc);
        std::printf("  sigma_x target %.2f, realised %.2f km/s (%+.2f%%)\n",
                    rc.sigma_x, r.sigma_x_realised,
                    100.0 * (r.sigma_x_realised / rc.sigma_x - 1.0));
        std::printf("  %d guiding centres, spacing %.3f a_epi; "
                    "stream support +-%.2f kpc vs Lx/2 = %.2f\n",
                    n_g, g.Lx / n_g / r.a_epi, support, 0.5 * g.Lx);
    }
    // The basis is EXACTLY zero outside +-support (the mask below), so wrap
    // self-overlap is only possible at all once the support window exceeds
    // half the box.
    if (support > 0.5 * g.Lx) {
        std::printf("  WARNING: stream support (+-%.2f kpc) exceeds Lx/2 = "
                    "%.2f: wrapped tails\n"
                    "  overlap themselves and the periodic-trap argument "
                    "degrades. Widen the box\n"
                    "  or lower the level cut (v_margin).\n",
                    support, 0.5 * g.Lx);
    }

    // ---- random phases, drawn up front for thread-order independence ----
    std::vector<double> theta(static_cast<size_t>(n_g) * M * n_levels);
    {
        SplitMix64 rng(rc.seed);
        for (auto& t : theta) t = 2.0 * units::PI * rng.uniform();
    }

    // Per-stream vertical amplitudes: mass_g = Sigma Lx / n_g exactly.
    std::vector<double> amp(M);
    for (int m = 0; m < M; ++m) {
        amp[m] = std::sqrt(s.w[m] * g.Lx / static_cast<double>(n_g));
    }
    std::vector<double> sqrtw(n_levels);
    for (int n = 0; n < n_levels; ++n) sqrtw[n] = std::sqrt(w[n]);

    r.st.psi.assign(n_g, std::vector<Complex>(g.size(), Complex(0.0, 0.0)));

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int t = 0; t < n_g; ++t) {
        // Hermite-function table phi_n(xi) / sqrt(ell) on this stream's grid,
        // by upward recurrence (stable here: going up in n at fixed xi moves
        // from the forbidden region toward the classical one, so the wanted
        // solution is the growing one). Rows are then renormalised on the
        // grid, which also absorbs discretisation drift in the analytic norm.
        std::vector<double> H(static_cast<size_t>(n_levels) * g.Nx, 0.0);
        std::vector<char> in_support(g.Nx, 0);
        const double inv_sqrt_ell = 1.0 / std::sqrt(r.ell);
        for (int i = 0; i < g.Nx; ++i) {
            const double u = wrap_dx(g.x(i) - r.xg[t], g.Lx);
            const double xi = u / r.ell;
            if (std::fabs(xi) > xi_max) continue;
            in_support[i] = 1;
            double pm1 = 0.0;
            double p = 0.7511255444649425 * std::exp(-0.5 * xi * xi) *
                       inv_sqrt_ell;  // pi^(-1/4) e^(-xi^2/2) / sqrt(ell)
            H[i] = p;
            for (int n = 1; n < n_levels; ++n) {
                const double pn = std::sqrt(2.0 / n) * xi * p -
                                  std::sqrt((n - 1.0) / n) * pm1;
                pm1 = p;
                p = pn;
                H[static_cast<size_t>(n) * g.Nx + i] = p;
            }
        }
        for (int n = 0; n < n_levels; ++n) {
            double s2 = 0.0;
            double* row = &H[static_cast<size_t>(n) * g.Nx];
            for (int i = 0; i < g.Nx; ++i) s2 += row[i] * row[i];
            s2 *= g.dx;
            const double inv = (s2 > 0.0) ? 1.0 / std::sqrt(s2) : 0.0;
            for (int i = 0; i < g.Nx; ++i) row[i] *= inv;
        }

        // In-plane amplitude B_m(x) per vertical state, then the outer
        // product against u_m(z).
        std::vector<Complex> B(g.Nx);
        auto& psi = r.st.psi[t];
        for (int m = 0; m < M; ++m) {
            const size_t off = (static_cast<size_t>(t) * M + m) *
                               static_cast<size_t>(n_levels);
            for (int i = 0; i < g.Nx; ++i) {
                if (!in_support[i]) { B[i] = Complex(0.0, 0.0); continue; }
                Complex acc(0.0, 0.0);
                for (int n = 0; n < n_levels; ++n) {
                    const double th = theta[off + n];
                    acc += sqrtw[n] * Complex(std::cos(th), std::sin(th)) *
                           H[static_cast<size_t>(n) * g.Nx + i];
                }
                B[i] = acc;
            }
            const auto& um = s.u[m];
            for (int i = 0; i < g.Nx; ++i) {
                if (!in_support[i]) continue;
                const Complex c = amp[m] * B[i];
                Complex* row = &psi[g.idx(i, 0)];
                for (int j = 0; j < g.Nz; ++j) row[j] += c * um[j];
            }
        }
    }
    return r;
}
