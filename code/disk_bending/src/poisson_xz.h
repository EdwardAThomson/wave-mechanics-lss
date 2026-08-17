// Hybrid Poisson solver for a slab: periodic in x, isolated in z.
//
// ----------------------------------------------------------------------------
// The method
// ----------------------------------------------------------------------------
// Fourier transform in x, and the 2D Poisson equation separates into one
// ordinary differential equation per in-plane wavenumber,
//
//     d^2 Phihat/dz^2 - kx^2 Phihat = 4 pi G rhohat,
//
// whose decaying-boundary Green's function is -exp(-|kx||z-z'|) / (2|kx|), so
//
//     Phihat(kx, z) = -(2 pi G / |kx|) Integral rhohat(kx, z') e^{-|kx||z-z'|} dz'.
//
// The key point is that this integral, like the |z - z'| kernel of the 1D case,
// collapses to two-term recursions. Splitting it at z' = z,
//
//     L_j = e^{-k dz} L_{j-1} + rho_j dz        (sweep up)
//     R_j = e^{-k dz} R_{j+1} + rho_j dz        (sweep down)
//     S_j = L_j + R_j - rho_j dz                (j counted twice)
//
// so each wavenumber costs O(Nz) and the whole solve is O(Nx log Nx * Nz),
// exact for the quadrature rule, with no zero-padding and no periodic images
// in z. This is cheaper in memory than doubling the box and, unlike the
// zero-pad trick, it is exact rather than approximate.
//
// The kx = 0 mode is the horizontally averaged slab. Its kernel is the |z - z'|
// one from poisson_z.h: taking k -> 0 above gives -(2 pi G / k) * M, a constant
// in z that is the (infinite) self-energy of an infinite sheet and drops out of
// the force, plus exactly the 1D isolated result. So the two solvers agree in
// the limit by construction, which is what test_poisson_xz checks first.
// ----------------------------------------------------------------------------
#pragma once

#include <cmath>
#include <complex>
#include <vector>

#include <fftw3.h>

#include "grid2d.h"
#include "units.h"

using Complex = std::complex<double>;

// Apply the isolated 1D vertical kernel Phi = 2 pi G Integral rho |z - z'| dz'
// to one complex column, in O(Nz). Gauge: Phi(z = 0) = 0.
// Set `gauge_zero` false for energy work: the Phi(0) = 0 convention adds a
// constant that depends on the density profile, so it drifts as the profile
// evolves and shows up as a spurious energy non-conservation.
inline void isolated_column_k0(const Grid2D& g, const Complex* rho,
                               Complex* Phi, bool gauge_zero = true) {
    const int N = g.Nz;
    std::vector<Complex> B(N + 1, Complex(0.0, 0.0)), T(N + 1, Complex(0.0, 0.0));
    for (int j = N - 1; j >= 0; --j) {
        const Complex m = rho[j] * g.dz;
        B[j] = B[j + 1] + m;
        T[j] = T[j + 1] + m * g.z(j);
    }
    Complex A(0.0, 0.0), S(0.0, 0.0);
    const double c = 2.0 * units::PI * units::G;
    for (int j = 0; j < N; ++j) {
        const Complex m = rho[j] * g.dz;
        const double zj = g.z(j);
        const Complex Bj = B[j] - m;
        const Complex Tj = T[j] - m * zj;
        Phi[j] = c * (zj * A - S + Tj - zj * Bj);
        A += m;
        S += m * zj;
    }
    if (!gauge_zero) return;
    Complex phi0 = (N % 2 == 0) ? 0.5 * (Phi[N / 2 - 1] + Phi[N / 2]) : Phi[N / 2];
    for (int j = 0; j < N; ++j) Phi[j] -= phi0;
}

// Apply the exponential kernel for one non-zero |kx|, in O(Nz).
inline void isolated_column_k(const Grid2D& g, double k, const Complex* rho,
                              Complex* Phi) {
    const int N = g.Nz;
    const double decay = std::exp(-k * g.dz);
    std::vector<Complex> R(N + 1, Complex(0.0, 0.0));
    for (int j = N - 1; j >= 0; --j) R[j] = decay * R[j + 1] + rho[j] * g.dz;

    const double pref = -2.0 * units::PI * units::G / k;
    Complex L(0.0, 0.0);
    for (int j = 0; j < N; ++j) {
        const Complex m = rho[j] * g.dz;
        L = decay * L + m;
        Phi[j] = pref * (L + R[j] - m);
    }
}

class PoissonXZ {
public:
    explicit PoissonXZ(const Grid2D& g) : g_(g) {
        buf_ = fftw_alloc_complex(g_.size());
        int n[1] = {g_.Nx};
        // One length-Nx transform per z row. Consecutive x samples are Nz
        // apart in memory (stride), consecutive rows are 1 apart (dist).
        fwd_ = fftw_plan_many_dft(1, n, g_.Nz, buf_, nullptr, g_.Nz, 1, buf_,
                                  nullptr, g_.Nz, 1, FFTW_FORWARD, FFTW_MEASURE);
        bwd_ = fftw_plan_many_dft(1, n, g_.Nz, buf_, nullptr, g_.Nz, 1, buf_,
                                  nullptr, g_.Nz, 1, FFTW_BACKWARD, FFTW_MEASURE);
    }
    ~PoissonXZ() {
        fftw_destroy_plan(fwd_);
        fftw_destroy_plan(bwd_);
        fftw_free(buf_);
    }
    PoissonXZ(const PoissonXZ&) = delete;
    PoissonXZ& operator=(const PoissonXZ&) = delete;

    void solve(const std::vector<double>& rho, std::vector<double>& Phi,
               bool gauge_zero = true) {
        const int Nx = g_.Nx, Nz = g_.Nz;
        for (int i = 0; i < Nx * Nz; ++i) {
            buf_[i][0] = rho[i];
            buf_[i][1] = 0.0;
        }
        fftw_execute(fwd_);

        std::vector<Complex> col(Nz), out(Nz);
        for (int i = 0; i < Nx; ++i) {
            Complex* base = reinterpret_cast<Complex*>(buf_ + i * Nz);
            const double k = std::fabs(g_.kx(i));
            if (k == 0.0) {
                isolated_column_k0(g_, base, out.data(), gauge_zero);
            } else {
                isolated_column_k(g_, k, base, out.data());
            }
            for (int j = 0; j < Nz; ++j) base[j] = out[j];
        }

        fftw_execute(bwd_);
        Phi.resize(Nx * Nz);
        const double s = 1.0 / static_cast<double>(Nx);
        for (int i = 0; i < Nx * Nz; ++i) Phi[i] = buf_[i][0] * s;
    }

private:
    Grid2D g_;
    fftw_complex* buf_;
    fftw_plan fwd_, bwd_;
};

// Direct O(Nz^2) evaluation of one column's kernel, for certifying the
// recursions. Not used in production.
inline void isolated_column_direct(const Grid2D& g, double k,
                                   const std::vector<Complex>& rho,
                                   std::vector<Complex>& Phi) {
    const int N = g.Nz;
    Phi.assign(N, Complex(0.0, 0.0));
    if (k == 0.0) {
        const double c = 2.0 * units::PI * units::G;
        for (int j = 0; j < N; ++j) {
            Complex s(0.0, 0.0);
            for (int i = 0; i < N; ++i) {
                s += rho[i] * g.dz * std::fabs(g.z(j) - g.z(i));
            }
            Phi[j] = c * s;
        }
        Complex phi0 =
            (N % 2 == 0) ? 0.5 * (Phi[N / 2 - 1] + Phi[N / 2]) : Phi[N / 2];
        for (int j = 0; j < N; ++j) Phi[j] -= phi0;
    } else {
        const double pref = -2.0 * units::PI * units::G / k;
        for (int j = 0; j < N; ++j) {
            Complex s(0.0, 0.0);
            for (int i = 0; i < N; ++i) {
                s += rho[i] * g.dz * std::exp(-k * std::fabs(g.z(j) - g.z(i)));
            }
            Phi[j] = pref * s;
        }
    }
}
