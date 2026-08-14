// Lowest-M eigenpairs of a real symmetric tridiagonal matrix.
//
// Sturm-sequence bisection for the eigenvalues, then shifted inverse iteration
// with full reorthogonalisation for the eigenvectors. Self-contained so the
// build keeps FFTW as its only external dependency, matching the rest of the
// repository. Correctness is pinned by tests/test_harmonic.cpp, which checks
// the computed spectrum against the analytic ladder E_n = hbar omega (n + 1/2).
//
// The matrix is the discrete vertical Hamiltonian
//     H = -(hbar_eff^2 / 2) d^2/dz^2 + V(z)
// with the standard three-point Laplacian and Dirichlet ends. Dirichlet is
// harmless here because every state we occupy has decayed to ~0 long before
// the box edge; see README.md for why the *evolution* is still done on the
// periodic FFT grid.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

struct SymTridiagEigen {
    std::vector<double> values;               // ascending, size M
    std::vector<std::vector<double>> vectors; // size M, each length N, 2-norm 1
};

namespace tridiag_detail {

// Number of eigenvalues strictly less than x (Sturm sequence).
inline int sturm_count(const std::vector<double>& d,
                       const std::vector<double>& e,
                       double x, double tiny) {
    const int n = static_cast<int>(d.size());
    double q = d[0] - x;
    int cnt = (q < 0.0) ? 1 : 0;
    for (int i = 1; i < n; ++i) {
        if (std::fabs(q) < tiny) q = (q < 0.0) ? -tiny : tiny;
        q = d[i] - x - e[i - 1] * e[i - 1] / q;
        if (q < 0.0) ++cnt;
    }
    return cnt;
}

// LU factorisation of (T - lambda I) with partial pivoting, stored as three
// upper diagonals uu/uv/uw plus the multipliers and swap flags.
struct TridiagLU {
    std::vector<double> uu, uv, uw, mult;
    std::vector<char> swapped;
};

inline TridiagLU factorise(const std::vector<double>& d,
                           const std::vector<double>& e,
                           double lambda, double tiny) {
    const int n = static_cast<int>(d.size());
    TridiagLU lu;
    lu.uu.assign(n, 0.0);
    lu.uv.assign(n, 0.0);
    lu.uw.assign(n, 0.0);
    lu.mult.assign(n, 0.0);
    lu.swapped.assign(n, 0);

    for (int i = 0; i < n; ++i) {
        lu.uu[i] = d[i] - lambda;
        lu.uv[i] = (i < n - 1) ? e[i] : 0.0;
    }

    for (int i = 0; i < n - 1; ++i) {
        const double sub = e[i];  // entry of row i+1 in column i
        if (std::fabs(lu.uu[i]) >= std::fabs(sub)) {
            lu.swapped[i] = 0;
            double piv = lu.uu[i];
            if (std::fabs(piv) < tiny) piv = lu.uu[i] = (piv < 0.0 ? -tiny : tiny);
            const double m = sub / piv;
            lu.mult[i] = m;
            const double next_diag = lu.uu[i + 1];      // d[i+1] - lambda
            const double next_super = lu.uv[i + 1];     // e[i+1]
            lu.uu[i + 1] = next_diag - m * lu.uv[i];
            lu.uv[i + 1] = next_super - m * lu.uw[i];
            // lu.uw[i+1] remains 0
        } else {
            lu.swapped[i] = 1;
            const double old_uu = lu.uu[i], old_uv = lu.uv[i], old_uw = lu.uw[i];
            const double next_diag = lu.uu[i + 1];
            const double next_super = lu.uv[i + 1];
            // Pivot row becomes the (i+1)-th original row.
            lu.uu[i] = sub;
            lu.uv[i] = next_diag;
            lu.uw[i] = next_super;
            const double m = old_uu / sub;
            lu.mult[i] = m;
            lu.uu[i + 1] = old_uv - m * next_diag;
            lu.uv[i + 1] = old_uw - m * next_super;
            lu.uw[i + 1] = 0.0;
        }
    }
    if (std::fabs(lu.uu[n - 1]) < tiny) {
        lu.uu[n - 1] = (lu.uu[n - 1] < 0.0) ? -tiny : tiny;
    }
    return lu;
}

inline void lu_solve(const TridiagLU& lu, std::vector<double>& b, double tiny) {
    const int n = static_cast<int>(b.size());
    for (int i = 0; i < n - 1; ++i) {
        if (!lu.swapped[i]) {
            b[i + 1] -= lu.mult[i] * b[i];
        } else {
            const double t = b[i];
            b[i] = b[i + 1];
            b[i + 1] = t - lu.mult[i] * b[i];
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        double s = b[i];
        if (i + 1 < n) s -= lu.uv[i] * b[i + 1];
        if (i + 2 < n) s -= lu.uw[i] * b[i + 2];
        double piv = lu.uu[i];
        if (std::fabs(piv) < tiny) piv = (piv < 0.0) ? -tiny : tiny;
        b[i] = s / piv;
    }
}

}  // namespace tridiag_detail

// Lowest `M` eigenpairs of the symmetric tridiagonal matrix with diagonal `d`
// and off-diagonal `e` (e[i] couples i and i+1, length N-1).
// `orth_window` limits reorthogonalisation to that many preceding vectors; for
// a smooth 1D well the spectrum is non-degenerate and a window of a few tens is
// ample, which keeps the cost O(M * window * N) rather than O(M^2 N). Pass a
// value >= M to reorthogonalise fully.
inline SymTridiagEigen sym_tridiag_lowest(const std::vector<double>& d,
                                          const std::vector<double>& e,
                                          int M, int orth_window = 40) {
    const int n = static_cast<int>(d.size());
    M = std::min(M, n);
    SymTridiagEigen out;
    out.values.assign(M, 0.0);
    out.vectors.assign(M, std::vector<double>(n, 0.0));

    // Gershgorin bounds and a scale for the singularity guards.
    double lo = 1e300, hi = -1e300, anorm = 0.0;
    for (int i = 0; i < n; ++i) {
        const double l = (i > 0) ? std::fabs(e[i - 1]) : 0.0;
        const double r = (i < n - 1) ? std::fabs(e[i]) : 0.0;
        lo = std::min(lo, d[i] - l - r);
        hi = std::max(hi, d[i] + l + r);
        anorm = std::max(anorm, std::fabs(d[i]) + l + r);
    }
    const double eps = 2.220446049250313e-16;
    const double tiny = eps * anorm;
    const double tol = 4.0 * eps * std::max(anorm, 1.0);

    // --- eigenvalues by bisection ---
    // The eigenvalues come out ascending, so each one can start from the
    // previous one as a lower bracket.
    double prev = lo;
    for (int m = 0; m < M; ++m) {
        double a = prev, b = hi;
        for (int it = 0; it < 200 && (b - a) > tol; ++it) {
            const double mid = 0.5 * (a + b);
            if (tridiag_detail::sturm_count(d, e, mid, tiny) <= m) {
                a = mid;
            } else {
                b = mid;
            }
        }
        out.values[m] = 0.5 * (a + b);
        prev = a;
    }

    // --- eigenvectors by shifted inverse iteration ---
    std::vector<double> v(n);
    for (int m = 0; m < M; ++m) {
        // Perturb the shift off the eigenvalue so the factorisation is not
        // exactly singular, and so degenerate-looking pairs separate.
        double shift = out.values[m];
        const double nudge = std::max(1e3 * eps * std::max(anorm, 1.0),
                                      1e-11 * std::fabs(shift));
        shift += ((m % 2) ? nudge : -nudge);
        auto lu = tridiag_detail::factorise(d, e, shift, tiny);

        // Start from a deterministic but non-symmetric seed so that no state
        // is orthogonal to it by parity.
        for (int i = 0; i < n; ++i) {
            v[i] = 1.0 / std::sqrt(static_cast<double>(n)) *
                   (1.0 + 0.5 * std::sin(0.7 * (i + 1) + 1.3 * (m + 1)));
        }

        for (int it = 0; it < 3; ++it) {
            tridiag_detail::lu_solve(lu, v, tiny);
            // Reorthogonalise against the recent window of found vectors.
            const int p0 = (orth_window >= m) ? 0 : (m - orth_window);
            for (int p = p0; p < m; ++p) {
                const auto& u = out.vectors[p];
                double dot = 0.0;
                for (int i = 0; i < n; ++i) dot += u[i] * v[i];
                for (int i = 0; i < n; ++i) v[i] -= dot * u[i];
            }
            double nrm = 0.0;
            for (int i = 0; i < n; ++i) nrm += v[i] * v[i];
            nrm = std::sqrt(nrm);
            if (nrm < 1e-300) {  // degenerate seed; reseed and retry
                for (int i = 0; i < n; ++i) v[i] = std::sin((i + 1) * (m + 2) * 0.37);
                nrm = 0.0;
                for (int i = 0; i < n; ++i) nrm += v[i] * v[i];
                nrm = std::sqrt(nrm);
            }
            for (int i = 0; i < n; ++i) v[i] /= nrm;
        }

        // Fix the sign convention: largest-magnitude component positive.
        int imax = 0;
        for (int i = 1; i < n; ++i) {
            if (std::fabs(v[i]) > std::fabs(v[imax])) imax = i;
        }
        if (v[imax] < 0.0) {
            for (int i = 0; i < n; ++i) v[i] = -v[i];
        }
        out.vectors[m] = v;

        // Rayleigh quotient refinement of the eigenvalue (cheap, and it
        // removes the deliberate shift nudge).
        double num = 0.0;
        for (int i = 0; i < n; ++i) {
            double hv = d[i] * v[i];
            if (i > 0) hv += e[i - 1] * v[i - 1];
            if (i < n - 1) hv += e[i] * v[i + 1];
            num += v[i] * hv;
        }
        out.values[m] = num;
    }

    return out;
}

// Max |H u - lambda u| over the returned set, as a fraction of the spectral
// scale. Used by the tests to certify the solver.
inline double eigen_residual(const std::vector<double>& d,
                             const std::vector<double>& e,
                             const SymTridiagEigen& res) {
    const int n = static_cast<int>(d.size());
    double worst = 0.0;
    for (size_t m = 0; m < res.values.size(); ++m) {
        const auto& v = res.vectors[m];
        double r = 0.0;
        for (int i = 0; i < n; ++i) {
            double hv = d[i] * v[i];
            if (i > 0) hv += e[i - 1] * v[i - 1];
            if (i < n - 1) hv += e[i] * v[i + 1];
            const double dv = hv - res.values[m] * v[i];
            r += dv * dv;
        }
        worst = std::max(worst, std::sqrt(r));
    }
    return worst;
}
