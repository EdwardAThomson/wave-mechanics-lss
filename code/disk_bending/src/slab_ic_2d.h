// Warm initial conditions for the 2D (x, z) slab.
//
// The equilibrium is homogeneous in x and stratified in z, so the construction
// separates exactly and the whole Stage 0 eigenstate library is reused
// untouched:
//
//   psi_t(x, z) = (1/sqrt(T)) sum_n sqrt(w_n) u_n(z) A_{t,n}(x)
//   A_{t,n}(x)  = sum_kx sqrt(g(kx)) exp(i theta) exp(i kx x)
//
// with u_n, w_n straight out of build_warm_sheet, and g(kx) the in-plane
// velocity distribution. Since an in-plane plane wave carries v_x = hbar kx,
// a Maxwellian of dispersion sigma_x means g(kx) proportional to
// exp(-hbar^2 kx^2 / 2 sigma_x^2), normalised to sum to one.
//
// This is Option A vertically (exact stationarity in z) and Option B in plane
// (random phases), the latter being forced: there is no in-plane eigenbasis to
// occupy, because x is a continuum of degenerate plane waves. Extra streams
// therefore only average down the in-plane speckle, at 1/sqrt(T), and that is
// exactly the quantity correction 5 in the README predicts will additionally
// average by sqrt(N_x) once a coherent k_x mode is fitted.
#pragma once

#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>

#include "equilibrium.h"
#include "evolve2d.h"
#include "grid2d.h"
#include "slab_ic.h"  // SplitMix64

#ifdef _OPENMP
#include <omp.h>
#endif

// `n_streams` independent realisations, summed incoherently. Each carries the
// full vertical eigenstate library.
inline SlabState2D make_warm_2d(const WarmSheet& s, const Grid2D& g,
                                double sigma_x, int n_streams,
                                uint64_t seed = 20260814ULL) {
    SlabState2D st;
    st.grid = g;
    st.hbar = s.cfg.hbar;
    st.psi.assign(n_streams, std::vector<Complex>(g.size(), Complex(0.0, 0.0)));

    const int M = static_cast<int>(s.u.size());
    const double hb = s.cfg.hbar;

    // In-plane weight per mode: Maxwellian in v_x = hbar kx.
    std::vector<double> gk(g.Nx);
    double gsum = 0.0;
    for (int i = 0; i < g.Nx; ++i) {
        const double vx = hb * g.kx(i);
        gk[i] = std::exp(-0.5 * vx * vx / (sigma_x * sigma_x));
        gsum += gk[i];
    }
    for (int i = 0; i < g.Nx; ++i) gk[i] = std::sqrt(gk[i] / gsum);

    const double stream_norm = 1.0 / std::sqrt(static_cast<double>(n_streams));

    // Pre-draw all phases so the result does not depend on thread scheduling.
    std::vector<double> theta(static_cast<size_t>(n_streams) * M * g.Nx);
    {
        SplitMix64 rng(seed);
        for (auto& t : theta) t = 2.0 * units::PI * rng.uniform();
    }

    // FFTW planning is not thread safe, so build the per-thread length-Nx
    // plans serially before the parallel region.
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    std::vector<std::unique_ptr<FFT1D>> ffts;
    for (int i = 0; i < nthreads; ++i) ffts.emplace_back(new FFT1D(g.Nx));

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int t = 0; t < n_streams; ++t) {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        FFT1D& f = *ffts[tid];
        std::vector<Complex> A(g.Nx), Ax(g.Nx);
        for (int m = 0; m < M; ++m) {
            // A(x) is the inverse transform of this (stream, eigenstate)'s
            // random-phase in-plane spectrum. The half-cell offset between
            // g.x(i) and the raw DFT grid is a per-mode phase, and the phases
            // are random already, so it can be dropped.
            const size_t off = (static_cast<size_t>(t) * M + m) * g.Nx;
            for (int i = 0; i < g.Nx; ++i) {
                const double th = theta[off + i];
                A[i] = gk[i] * Complex(std::cos(th), std::sin(th));
            }
            f.load(A);
            f.backward();          // carries a 1/Nx we do not want here
            f.store(Ax);
            const double a = std::sqrt(s.w[m]) * stream_norm * g.Nx;
            for (int i = 0; i < g.Nx; ++i) {
                const Complex c = a * Ax[i];
                for (int j = 0; j < g.Nz; ++j) {
                    st.psi[t][g.idx(i, j)] += c * s.u[m][j];
                }
            }
        }
    }
    return st;
}
