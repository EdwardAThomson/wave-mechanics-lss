// Strang split-step on the 2D (x, z) slab.
//
// Same structure as evolve.h: the kinetic factor is applied spectrally on a
// grid that is genuinely periodic in x and merely "wide enough" in z, while
// gravity gets the isolated treatment in poisson_xz.h.
//
// One asymmetry is deliberate. The MatchedFD option uses the finite-difference
// wavenumber in z but the exact one in x. That is not a fudge: the warm
// equilibrium is built from eigenstates of the three-point vertical
// Hamiltonian times exact in-plane plane waves, and plane waves are
// eigenfunctions of both the exact and the FD in-plane operator, so matching is
// only needed in z. The resulting discrete Hamiltonian is separable and the
// initial condition is an exact eigenbasis of it.
#pragma once

#include <cmath>
#include <complex>
#include <memory>
#include <vector>

#include <fftw3.h>

#include "evolve.h"   // for the Kinetic enum
#include "grid2d.h"
#include "poisson_xz.h"
#include "units.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// RAII pair of 2D plans plus a z-only batched pair, sharing one buffer.
class FFT2D {
public:
    explicit FFT2D(const Grid2D& g) : g_(g) {
        buf_ = fftw_alloc_complex(g_.size());
        fwd2_ = fftw_plan_dft_2d(g_.Nx, g_.Nz, buf_, buf_, FFTW_FORWARD,
                                 FFTW_MEASURE);
        bwd2_ = fftw_plan_dft_2d(g_.Nx, g_.Nz, buf_, buf_, FFTW_BACKWARD,
                                 FFTW_MEASURE);
        int nz[1] = {g_.Nz};
        // Nx transforms of length Nz; z is contiguous so stride 1, dist Nz.
        fwdz_ = fftw_plan_many_dft(1, nz, g_.Nx, buf_, nullptr, 1, g_.Nz, buf_,
                                   nullptr, 1, g_.Nz, FFTW_FORWARD, FFTW_MEASURE);
        bwdz_ = fftw_plan_many_dft(1, nz, g_.Nx, buf_, nullptr, 1, g_.Nz, buf_,
                                   nullptr, 1, g_.Nz, FFTW_BACKWARD, FFTW_MEASURE);
    }
    ~FFT2D() {
        fftw_destroy_plan(fwd2_);
        fftw_destroy_plan(bwd2_);
        fftw_destroy_plan(fwdz_);
        fftw_destroy_plan(bwdz_);
        fftw_free(buf_);
    }
    FFT2D(const FFT2D&) = delete;
    FFT2D& operator=(const FFT2D&) = delete;

    Complex* data() { return reinterpret_cast<Complex*>(buf_); }
    const Complex* data() const { return reinterpret_cast<const Complex*>(buf_); }

    void load(const std::vector<Complex>& a) {
        std::copy(a.begin(), a.end(), data());
    }
    void store(std::vector<Complex>& a) const {
        std::copy(data(), data() + g_.size(), a.begin());
    }

    void forward2() { fftw_execute(fwd2_); }
    void backward2() {
        fftw_execute(bwd2_);
        scale(1.0 / static_cast<double>(g_.size()));
    }
    void forward_z() { fftw_execute(fwdz_); }
    void backward_z() {
        fftw_execute(bwdz_);
        scale(1.0 / static_cast<double>(g_.Nz));
    }

private:
    void scale(double s) {
        Complex* d = data();
        for (int i = 0; i < g_.size(); ++i) d[i] *= s;
    }
    Grid2D g_;
    fftw_complex* buf_;
    fftw_plan fwd2_, bwd2_, fwdz_, bwdz_;
};

struct SlabState2D {
    Grid2D grid;
    double hbar = 0.6;
    std::vector<std::vector<Complex>> psi;

    int n_streams() const { return static_cast<int>(psi.size()); }

    void density(std::vector<double>& rho) const {
        rho.assign(grid.size(), 0.0);
        for (const auto& p : psi) {
            for (int i = 0; i < grid.size(); ++i) rho[i] += std::norm(p[i]);
        }
    }
};

class Evolver2D {
public:
    Evolver2D(const Grid2D& g, double hbar, Kinetic kin)
        : g_(g), hbar_(hbar), kin_(kin) {
        int nthreads = 1;
#ifdef _OPENMP
        nthreads = omp_get_max_threads();
#endif
        for (int t = 0; t < nthreads; ++t) ffts_.emplace_back(new FFT2D(g_));

        k2_.resize(g_.size());
        for (int i = 0; i < g_.Nx; ++i) {
            const double kx = g_.kx(i);
            for (int j = 0; j < g_.Nz; ++j) {
                const double kz = g_.kz(j);
                double kz2;
                if (kin_ == Kinetic::Spectral) {
                    kz2 = kz * kz;
                } else {
                    const double s = std::sin(0.5 * kz * g_.dz);
                    kz2 = 4.0 * s * s / (g_.dz * g_.dz);
                }
                k2_[g_.idx(i, j)] = kx * kx + kz2;
            }
        }
    }

    const Grid2D& grid() const { return g_; }
    double hbar() const { return hbar_; }

    double v_ceiling_x() const { return hbar_ * g_.kx_nyquist(); }
    double v_ceiling_z() const {
        return (kin_ == Kinetic::Spectral) ? hbar_ * g_.kz_nyquist()
                                           : hbar_ / g_.dz;
    }

    void kinetic_half(SlabState2D& st, double dt) const {
        const int S = st.n_streams(), N = g_.size();
        const double a = -0.25 * hbar_ * dt;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int s = 0; s < S; ++s) {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            FFT2D& f = *ffts_[tid];
            f.load(st.psi[s]);
            f.forward2();
            Complex* d = f.data();
            for (int i = 0; i < N; ++i) {
                const double ph = a * k2_[i];
                d[i] *= Complex(std::cos(ph), std::sin(ph));
            }
            f.backward2();
            f.store(st.psi[s]);
        }
    }

    void potential_full(SlabState2D& st, const std::vector<double>& V,
                        double dt) const {
        const int S = st.n_streams(), N = g_.size();
        const double a = -dt / hbar_;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int s = 0; s < S; ++s) {
            for (int i = 0; i < N; ++i) {
                const double ph = a * V[i];
                st.psi[s][i] *= Complex(std::cos(ph), std::sin(ph));
            }
        }
    }

    // Spectral d/dz, one call per stream.
    void ddz(const std::vector<Complex>& in, std::vector<Complex>& out) const {
        FFT2D& f = *ffts_[0];
        f.load(in);
        f.forward_z();
        Complex* d = f.data();
        for (int i = 0; i < g_.Nx; ++i) {
            for (int j = 0; j < g_.Nz; ++j) {
                d[g_.idx(i, j)] *= Complex(0.0, g_.kz(j));
            }
        }
        f.backward_z();
        out.resize(g_.size());
        f.store(out);
    }

    // Shift each x column vertically by shift[i], exactly, via the translation
    // operator exp(-i p d / hbar). This is how a bending perturbation is
    // imposed as a pure displacement with zero initial velocity.
    void displace_columns(SlabState2D& st,
                          const std::vector<double>& shift) const {
        for (int s = 0; s < st.n_streams(); ++s) {
            FFT2D& f = *ffts_[0];
            f.load(st.psi[s]);
            f.forward_z();
            Complex* d = f.data();
            for (int i = 0; i < g_.Nx; ++i) {
                for (int j = 0; j < g_.Nz; ++j) {
                    const double ph = -g_.kz(j) * shift[i];
                    d[g_.idx(i, j)] *= Complex(std::cos(ph), std::sin(ph));
                }
            }
            f.backward_z();
            f.store(st.psi[s]);
        }
    }

    double kinetic_energy(const SlabState2D& st) const {
        const int N = g_.size();
        double e = 0.0;
        for (int s = 0; s < st.n_streams(); ++s) {
            FFT2D& f = *ffts_[0];
            f.load(st.psi[s]);
            f.forward2();
            const Complex* d = f.data();
            for (int i = 0; i < N; ++i) {
                e += 0.5 * hbar_ * hbar_ * k2_[i] * std::norm(d[i]);
            }
        }
        return e * g_.dx * g_.dz / static_cast<double>(N);
    }

    // Fraction of |psi_k|^2 above the given ceilings, separately per axis.
    void spectral_spill(const SlabState2D& st, double kx_cut, double kz_cut,
                        double& spill_x, double& spill_z) const {
        double ax = 0.0, az = 0.0, tot = 0.0;
        for (int s = 0; s < st.n_streams(); ++s) {
            FFT2D& f = *ffts_[0];
            f.load(st.psi[s]);
            f.forward2();
            const Complex* d = f.data();
            for (int i = 0; i < g_.Nx; ++i) {
                for (int j = 0; j < g_.Nz; ++j) {
                    const double p = std::norm(d[g_.idx(i, j)]);
                    tot += p;
                    if (std::fabs(g_.kx(i)) > kx_cut) ax += p;
                    if (std::fabs(g_.kz(j)) > kz_cut) az += p;
                }
            }
        }
        spill_x = (tot > 0.0) ? ax / tot : 0.0;
        spill_z = (tot > 0.0) ? az / tot : 0.0;
    }

private:
    Grid2D g_;
    double hbar_;
    Kinetic kin_;
    std::vector<double> k2_;
    mutable std::vector<std::unique_ptr<FFT2D>> ffts_;
};

// ---------------------------------------------------------------------------
// Bending perturbations
// ---------------------------------------------------------------------------
// Two ways to excite a bending mode of wavenumber k, and the difference
// matters.
//
// `displace_columns` above shifts the equilibrium vertically by h0 cos(kx),
// leaving the velocity field untouched. That is the clean initial condition
// for measuring a dispersion relation: <z> starts at its maximum and <v_z> at
// zero, so the subsequent oscillation frequency is read off directly.
//
// The velocity kick below multiplies by exp(i v0 cos(kx) z / hbar). Note that
// it necessarily also imparts an in-plane velocity: the Madelung velocity
// field is a gradient, so demanding dS/dz = v0 cos(kx) forces
// dS/dx = -v0 k sin(kx) z. That is not a numerical artefact and cannot be
// removed. Any perturbation derived from a potential has the same property,
// including a real perturber's tidal field, so this is the more physical of
// the two and the less convenient one to interpret. Its size relative to the
// vertical kick is of order k z, so a few tens of per cent here.
inline void kick_bending(SlabState2D& st, double v0, double k) {
    const double hb = st.hbar;
    for (auto& p : st.psi) {
        for (int i = 0; i < st.grid.Nx; ++i) {
            const double c = v0 * std::cos(k * st.grid.x(i));
            for (int j = 0; j < st.grid.Nz; ++j) {
                const double ph = c * st.grid.z(j) / hb;
                p[st.grid.idx(i, j)] *= Complex(std::cos(ph), std::sin(ph));
            }
        }
    }
}
