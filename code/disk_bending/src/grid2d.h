// 2D (x, z) grid for the bending-wave stage.
//
// x is periodic (in-plane, the direction the bending wave propagates along),
// z is isolated (vertical). Storage is row-major with z contiguous, so the
// vertical recursions in poisson_xz.h walk memory sequentially and the x
// transform is a strided batch.
#pragma once

#include <cmath>
#include <vector>

#include "units.h"

struct Grid2D {
    int Nx = 0, Nz = 0;
    double Lx = 0.0, Lz = 0.0;
    double dx = 0.0, dz = 0.0;

    Grid2D() = default;
    Grid2D(int nx, int nz, double lx, double lz)
        : Nx(nx), Nz(nz), Lx(lx), Lz(lz), dx(lx / nx), dz(lz / nz) {}

    inline int size() const { return Nx * Nz; }
    inline int idx(int i, int j) const { return i * Nz + j; }

    inline double x(int i) const { return (i + 0.5) * dx; }
    inline double z(int j) const { return -0.5 * Lz + (j + 0.5) * dz; }

    inline double kx(int i) const {
        const int n = (i <= Nx / 2) ? i : i - Nx;
        return 2.0 * units::PI * n / Lx;
    }
    inline double kz(int j) const {
        const int n = (j <= Nz / 2) ? j : j - Nz;
        return 2.0 * units::PI * n / Lz;
    }

    inline double kx_nyquist() const { return units::PI / dx; }
    inline double kz_nyquist() const { return units::PI / dz; }

    // Fundamental in-plane wavenumber, the longest bending mode the box holds.
    inline double k_fundamental() const { return 2.0 * units::PI / Lx; }
};

// Minimal periodic distance in x: u mapped into [-L/2, L/2). Used by the
// rotating-frame trap and its initial conditions (rotation.h) so that both
// wrap identically.
inline double wrap_dx(double u, double L) {
    return u - L * std::round(u / L);
}
