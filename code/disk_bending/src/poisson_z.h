// Vertical Poisson solvers for a slab.
//
// ----------------------------------------------------------------------------
// Why this file is not an FFT
// ----------------------------------------------------------------------------
// A disk slab is isolated in z: the correct boundary condition is
// dPhi/dz -> +/- 2 pi G Sigma as z -> +/- infinity, not periodicity. A periodic
// solve stacks infinite image sheets above and below, and because vertical
// gravity from a sheet does not fall off with distance, the images do not
// become negligible however tall you make the box. That is the one place a
// cosmological triply-periodic code cannot simply be reused.
//
// In *one* dimension the isolated solve has a closed form,
//
//     Phi(z) = 2 pi G Integral rho(z') |z - z'| dz'
//
// so no Green's-function zero-padding and no tridiagonal solve is needed at
// all. Expanding |z - z'| into left/right partial sums turns the O(N^2) double
// sum into two O(N) cumulative sweeps that reproduce it *exactly* (to the
// quadrature rule, with zero image contamination). That is what
// poisson_isolated_1d does.
//
// For the 2D (x, z) stage the generalisation is: FFT in x, then for each k_x
// solve d^2 Phihat/dz^2 - k_x^2 Phihat = 4 pi G rhohat with decaying BCs, which
// has the exact kernel exp(-|k_x||z - z'|) / (2|k_x|) and degenerates to the
// expression above at k_x = 0. That is cheaper in memory than zero-padding by
// a factor of two and it is exact rather than approximate, so it, not the
// zero-pad trick, is the recommended route.
// ----------------------------------------------------------------------------
#pragma once

#include <cmath>
#include <vector>

#include "fft.h"
#include "grid.h"
#include "units.h"

// Isolated (vacuum) vertical Poisson solve.
//
//   d^2 Phi / dz^2 = 4 pi G rho,     Phi'(+/- inf) = +/- 2 pi G Sigma
//
// Exact for the midpoint quadrature rule, O(N), no periodic images.
// If `gz` is non-null it receives dPhi/dz (so the downward force is -gz).
// The gauge is fixed by Phi(z = 0) = 0.
// Set `gauge_zero` false to keep the raw kernel value, which is what the
// gauge-independent self-energy W = pi G int int rho rho |z - z'| needs.
inline void poisson_isolated_1d(const Grid1D& g,
                                const std::vector<double>& rho,
                                std::vector<double>& Phi,
                                std::vector<double>* gz = nullptr,
                                bool gauge_zero = true) {
    const int N = g.N;
    Phi.assign(N, 0.0);

    // Right-hand partial sums: B_j = sum_{i>j} rho_i dz,  T_j = sum_{i>j} rho_i dz z_i
    std::vector<double> B(N + 1, 0.0), T(N + 1, 0.0);
    for (int j = N - 1; j >= 0; --j) {
        const double m = rho[j] * g.dz;
        B[j] = B[j + 1] + m;
        T[j] = T[j + 1] + m * g.z(j);
    }

    // Sweep left to right accumulating A_j = sum_{i<j} rho_i dz, S_j likewise.
    double A = 0.0, S = 0.0;
    for (int j = 0; j < N; ++j) {
        const double m = rho[j] * g.dz;
        const double zj = g.z(j);
        const double Bj = B[j] - m;   // strictly i > j
        const double Tj = T[j] - m * zj;

        // Phi_j / (2 pi G) = z_j A - S + T - z_j B
        Phi[j] = 2.0 * units::PI * units::G * (zj * A - S + Tj - zj * Bj);
        if (gz) (*gz)[j] = 2.0 * units::PI * units::G * (A - Bj);

        A += m;
        S += m * zj;
    }

    if (!gauge_zero) return;

    // Gauge: Phi(0) = 0. For even N, z = 0 sits on the face between the two
    // central cells, so average them.
    double phi0;
    if (N % 2 == 0) {
        phi0 = 0.5 * (Phi[N / 2 - 1] + Phi[N / 2]);
    } else {
        phi0 = Phi[N / 2];
    }
    for (int j = 0; j < N; ++j) Phi[j] -= phi0;
}

// Periodic vertical Poisson solve (spectral), used only by the Jeans test,
// where periodicity is the physically correct boundary condition.
//
//   d^2 Phi / dz^2 = 4 pi G (rho - rhobar)
inline void poisson_periodic_1d(const Grid1D& g,
                                const std::vector<double>& rho,
                                std::vector<double>& Phi,
                                FFT1D& fft) {
    const int N = g.N;
    Phi.assign(N, 0.0);

    double mean = 0.0;
    for (int j = 0; j < N; ++j) mean += rho[j];
    mean /= N;

    std::vector<double> src(N);
    for (int j = 0; j < N; ++j) src[j] = rho[j] - mean;

    fft.load(src);
    fft.forward();
    fft.set(0, Complex(0.0, 0.0));
    for (int j = 1; j < N; ++j) {
        const double k = g.k(j);
        fft.scale(j, -4.0 * units::PI * units::G / (k * k));
    }
    fft.backward();
    fft.store_real(Phi);
}

// Analytic isothermal (Spitzer) self-gravitating sheet, for validation.
//   rho(z) = rho0 sech^2(z / 2h),  h = sigma^2 / (2 pi G Sigma),  Sigma = 4 rho0 h
//   Phi(z) = 2 sigma^2 ln cosh(z / 2h)
struct IsothermalSheet {
    double Sigma;   // surface density, Msun/kpc^2
    double sigma;   // vertical velocity dispersion, km/s
    double h;       // scale height, kpc
    double rho0;    // midplane density, Msun/kpc^3

    IsothermalSheet(double Sigma_, double sigma_)
        : Sigma(Sigma_), sigma(sigma_) {
        h = sigma * sigma / (2.0 * units::PI * units::G * Sigma);
        rho0 = Sigma / (4.0 * h);
    }
    double rho(double z) const {
        const double c = std::cosh(z / (2.0 * h));
        return rho0 / (c * c);
    }
    double Phi(double z) const {
        return 2.0 * sigma * sigma * std::log(std::cosh(z / (2.0 * h)));
    }
    // Small-amplitude vertical oscillation frequency at the midplane, km/s/kpc
    double omega0() const {
        return std::sqrt(4.0 * units::PI * units::G * rho0);
    }
};
