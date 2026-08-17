// 1D vertical grid.
//
// Cell-centred: z_j = -L/2 + (j + 1/2) dz, j = 0 .. N-1.
// The grid is used both as an FFT grid (for the kinetic step, which is local
// and therefore genuinely periodic-safe once psi has decayed at the edges) and
// as a quadrature grid for the isolated-BC Poisson solve, which does NOT
// assume periodicity. Keeping those two roles separate is the whole trick for
// a slab geometry; see README.md.
#pragma once

#include <cmath>
#include <vector>

#include "units.h"

struct Grid1D {
    int N = 0;
    double L = 0.0;
    double dz = 0.0;

    Grid1D() = default;
    Grid1D(int n, double box) : N(n), L(box), dz(box / n) {}

    // Cell-centre coordinate
    inline double z(int j) const { return -0.5 * L + (j + 0.5) * dz; }

    // Signed FFT wavenumber for index j (rad/kpc)
    inline double k(int j) const {
        const int n = (j <= N / 2) ? j : j - N;
        return 2.0 * units::PI * n / L;
    }

    // Nyquist wavenumber
    inline double k_nyquist() const { return units::PI / dz; }

    // Midpoint-rule integral of a cell-centred field
    inline double integrate(const std::vector<double>& f) const {
        double s = 0.0;
        for (int j = 0; j < N; ++j) s += f[j];
        return s * dz;
    }
};
