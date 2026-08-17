// Strang split-step evolution for the vertical slab.
//
// psi(t + dt) = exp(-i K dt/2) exp(-i V dt / hbar) exp(-i K dt/2) psi(t)
//
// The kinetic factor is applied spectrally on the periodic FFT grid. That is
// legitimate even though the slab is *not* periodic, because the kinetic
// operator is local: all that is required is that psi has decayed to ~0 at the
// box edges. Periodicity is only a problem for gravity, which is long-ranged,
// and gravity is handled by the isolated solver in poisson_z.h. Separating
// those two roles is what makes a slab tractable with cosmological machinery.
//
// Two kinetic operators are offered:
//   Spectral  : exp(-i hbar k^2 dt / 2)              -- exact dispersion
//   MatchedFD : exp(-i hbar k_fd^2 dt / 2), with
//               k_fd^2 = (4/dz^2) sin^2(k dz/2)      -- the exact propagator
//               for the three-point finite-difference Laplacian
//
// MatchedFD exists because the warm equilibrium in equilibrium.h is built from
// eigenvectors of exactly that finite-difference Hamiltonian. Evolving with the
// matching propagator makes the equilibrium stationary to solver precision
// instead of breathing at the few-per-cent level that the operator mismatch
// would otherwise produce. Its drawback is that the group velocity
// hbar sin(k dz)/dz turns over and returns to zero at the Nyquist wavenumber,
// so its aliasing behaviour is worse than spectral; the difference is reported
// by tests/test_sheet.cpp so the trade is visible rather than assumed.
#pragma once

#include <cmath>
#include <complex>
#include <memory>
#include <vector>

#include "fft.h"
#include "grid.h"
#include "poisson_z.h"
#include "units.h"

#ifdef _OPENMP
#include <omp.h>
#endif

enum class Kinetic { Spectral, MatchedFD };

// The wavefunction ensemble. In multi-stream mode (the plan's Option A) each
// stream is an occupied eigenstate carrying weight w_n; in single-psi mode
// (Option B) there is one stream holding the coherent superposition.
struct SlabState {
    Grid1D grid;
    double hbar = 0.3;
    std::vector<std::vector<Complex>> psi;  // each already carries sqrt(w)

    int n_streams() const { return static_cast<int>(psi.size()); }

    void density(std::vector<double>& rho) const {
        const int N = grid.N;
        rho.assign(N, 0.0);
        for (const auto& p : psi) {
            for (int j = 0; j < N; ++j) rho[j] += std::norm(p[j]);
        }
    }

    double surface_density() const {
        std::vector<double> r;
        density(r);
        return grid.integrate(r);
    }
};

class Evolver {
public:
    Evolver(const Grid1D& g, double hbar, Kinetic kin)
        : grid_(g), hbar_(hbar), kin_(kin) {
        int nthreads = 1;
#ifdef _OPENMP
        nthreads = omp_get_max_threads();
#endif
        // FFTW planning is not thread-safe, so build the per-thread plans
        // serially up front; fftw_execute on distinct plans/buffers then is.
        for (int t = 0; t < nthreads; ++t) {
            ffts_.emplace_back(new FFT1D(g.N));
        }
        // Squared wavenumber used by the kinetic propagator.
        k2_.resize(g.N);
        for (int j = 0; j < g.N; ++j) {
            const double k = g.k(j);
            if (kin_ == Kinetic::Spectral) {
                k2_[j] = k * k;
            } else {
                const double s = std::sin(0.5 * k * g.dz);
                k2_[j] = 4.0 * s * s / (g.dz * g.dz);
            }
        }
    }

    double hbar() const { return hbar_; }
    const Grid1D& grid() const { return grid_; }
    Kinetic kinetic_mode() const { return kin_; }

    // Largest velocity the grid can represent, hbar * k_max. For the spectral
    // operator k_max is the Nyquist wavenumber; for MatchedFD the group
    // velocity peaks at k dz = pi/2 and falls back to zero, so the usable
    // ceiling is hbar/dz, a factor pi lower.
    double v_max() const {
        return (kin_ == Kinetic::Spectral) ? hbar_ * grid_.k_nyquist()
                                           : hbar_ / grid_.dz;
    }

    void kinetic_half(SlabState& st, double dt) const {
        const int N = grid_.N;
        const int S = st.n_streams();
        // Half a step of exp(-i hbar k^2 tau / 2) with tau = dt/2.
        const double a = -0.25 * hbar_ * dt;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int s = 0; s < S; ++s) {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            FFT1D& f = *ffts_[tid];
            f.load(st.psi[s]);
            f.forward();
            for (int j = 0; j < N; ++j) {
                const double ph = a * k2_[j];
                f.set(j, f.at(j) * Complex(std::cos(ph), std::sin(ph)));
            }
            f.backward();
            f.store(st.psi[s]);
        }
    }

    void potential_full(SlabState& st, const std::vector<double>& V,
                        double dt) const {
        const int N = grid_.N;
        const int S = st.n_streams();
        const double a = -dt / hbar_;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int s = 0; s < S; ++s) {
            for (int j = 0; j < N; ++j) {
                const double ph = a * V[j];
                st.psi[s][j] *= Complex(std::cos(ph), std::sin(ph));
            }
        }
    }

    // Spectral d/dz, used for the momentum density and kinetic energy.
    void ddz(const std::vector<Complex>& in, std::vector<Complex>& out) const {
        FFT1D& f = *ffts_[0];
        f.load(in);
        f.forward();
        for (int j = 0; j < grid_.N; ++j) {
            f.set(j, f.at(j) * Complex(0.0, grid_.k(j)));
        }
        f.backward();
        out.resize(grid_.N);
        f.store(out);
    }

    // Translate every stream by dz0 (exact, spectral).
    void translate(SlabState& st, double dz0) const {
        const int N = grid_.N;
        for (int s = 0; s < st.n_streams(); ++s) {
            FFT1D& f = *ffts_[0];
            f.load(st.psi[s]);
            f.forward();
            for (int j = 0; j < N; ++j) {
                const double ph = -grid_.k(j) * dz0;
                f.set(j, f.at(j) * Complex(std::cos(ph), std::sin(ph)));
            }
            f.backward();
            f.store(st.psi[s]);
        }
    }

    // Kinetic energy per unit area, computed with the SAME operator the
    // evolution uses. Using the spectral gradient here while stepping with
    // MatchedFD leaves an O((k dz)^2) offset that masquerades as an energy
    // conservation failure, so the two must be kept consistent.
    double kinetic_energy(const SlabState& st) const {
        const int N = grid_.N;
        double e = 0.0;
        for (int s = 0; s < st.n_streams(); ++s) {
            FFT1D& f = *ffts_[0];
            f.load(st.psi[s]);
            f.forward();
            for (int j = 0; j < N; ++j) {
                e += 0.5 * hbar_ * hbar_ * k2_[j] * std::norm(f.at(j));
            }
        }
        // FFTW's forward transform is unnormalised. With A_k = FFT(psi)_k, the
        // inverse carries the 1/N, so Parseval reads
        //   sum_j |d psi/dz|_j^2 = (1/N) sum_k k^2 |A_k|^2,
        // and the integral over z brings one further factor of dz.
        return e * grid_.dz / static_cast<double>(N);
    }

    // Fraction of |psi_k|^2 lying above |k| = k_cut, summed over streams.
    // This is the aliasing tell-tale: if it grows during a run, velocity-space
    // wrap-around is contaminating the result.
    double spectral_spill(const SlabState& st, double k_cut) const {
        const int N = grid_.N;
        double above = 0.0, total = 0.0;
        for (int s = 0; s < st.n_streams(); ++s) {
            FFT1D& f = *ffts_[0];
            f.load(st.psi[s]);
            f.forward();
            for (int j = 0; j < N; ++j) {
                const double p = std::norm(f.at(j));
                total += p;
                if (std::fabs(grid_.k(j)) > k_cut) above += p;
            }
        }
        return (total > 0.0) ? above / total : 0.0;
    }

private:
    Grid1D grid_;
    double hbar_;
    Kinetic kin_;
    std::vector<double> k2_;
    mutable std::vector<std::unique_ptr<FFT1D>> ffts_;
};

// ---------------------------------------------------------------------------
// Impulsive kicks
// ---------------------------------------------------------------------------
// In 1D the bending / breathing split of the plan's §6 is exactly the parity of
// the imposed velocity field:
//   bulk     v_z = v0            odd  vertical response -> <z> oscillates (bending)
//   compress v_z = A z           even vertical response -> <z^2> oscillates (breathing)
// A real disk crossing excites both; their ratio is set by impact geometry.

inline void kick_bulk(SlabState& st, double v0) {
    const double hb = st.hbar;
    for (auto& p : st.psi) {
        for (int j = 0; j < st.grid.N; ++j) {
            const double ph = v0 * st.grid.z(j) / hb;
            p[j] *= Complex(std::cos(ph), std::sin(ph));
        }
    }
}

inline void kick_compress(SlabState& st, double A) {
    const double hb = st.hbar;
    for (auto& p : st.psi) {
        for (int j = 0; j < st.grid.N; ++j) {
            const double z = st.grid.z(j);
            const double ph = 0.5 * A * z * z / hb;
            p[j] *= Complex(std::cos(ph), std::sin(ph));
        }
    }
}

// ---------------------------------------------------------------------------
// Potential assembly and one full Strang step
// ---------------------------------------------------------------------------
struct PotentialWorkspace {
    std::vector<double> rho, Phi_self, V;
};

inline void assemble_potential(const Grid1D& g, const SlabState& st,
                               const std::vector<double>& Phi_ext,
                               PotentialWorkspace& ws) {
    st.density(ws.rho);
    poisson_isolated_1d(g, ws.rho, ws.Phi_self);
    ws.V.resize(g.N);
    for (int j = 0; j < g.N; ++j) ws.V[j] = ws.Phi_self[j] + Phi_ext[j];
}

inline void strang_step(const Evolver& ev, SlabState& st,
                        const std::vector<double>& Phi_ext,
                        PotentialWorkspace& ws, double dt) {
    ev.kinetic_half(st, dt);
    assemble_potential(ev.grid(), st, Phi_ext, ws);
    ev.potential_full(st, ws.V, dt);
    ev.kinetic_half(st, dt);
}
