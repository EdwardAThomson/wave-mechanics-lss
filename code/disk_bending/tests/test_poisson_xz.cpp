// Validation test 5a: the hybrid (x, z) Poisson solver.
//
// Four checks, in increasing order of how much machinery they exercise.
//
//   (a) the O(Nz) recursions against direct O(Nz^2) evaluation of the same
//       kernel: certifies the recursion algebra to machine precision
//   (b) against a closed-form solution, for the convergence rate
//   (c) the kx = 0 column against the standalone 1D isolated solver: the two
//       must agree, since the exponential kernel degenerates to |z - z'|
//   (d) the full solver end to end on a separable density, which is the only
//       one of the four that exercises the strided batch transform and the
//       FFT normalisation, and therefore the one most likely to catch a bug
#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>

#include "../src/grid2d.h"
#include "../src/poisson_xz.h"
#include "../src/poisson_z.h"
#include "../src/units.h"

int main() {
    bool pass = true;
    std::printf("=== Test 5a: hybrid (x, z) Poisson ===\n\n");

    const double Sigma = 50.0 * units::MSUN_PC2;
    const double sigma_z = 20.0;
    IsothermalSheet sheet(Sigma, sigma_z);

    // ---- (a) recursions vs direct kernel evaluation ----
    {
        std::printf("(a) O(Nz) recursion vs direct O(Nz^2) kernel\n");
        std::printf("%10s %14s %16s\n", "k[1/kpc]", "max|rel diff|", "kernel");
        Grid2D g(8, 512, 8.0, 6.0);
        std::vector<Complex> rho(g.Nz), got(g.Nz), ref;
        for (int j = 0; j < g.Nz; ++j) {
            rho[j] = Complex(sheet.rho(g.z(j)), 0.3 * sheet.rho(g.z(j)));
        }
        for (double k : {0.0, 0.1, 0.7854, 3.0, 20.0}) {
            if (k == 0.0) {
                isolated_column_k0(g, rho.data(), got.data());
            } else {
                isolated_column_k(g, k, rho.data(), got.data());
            }
            isolated_column_direct(g, k, rho, ref);
            double worst = 0.0, scale = 0.0;
            for (int j = 0; j < g.Nz; ++j) {
                worst = std::max(worst, std::abs(got[j] - ref[j]));
                scale = std::max(scale, std::abs(ref[j]));
            }
            std::printf("%10.4f %14.3e %16s\n", k, worst / scale,
                        k == 0.0 ? "|z-z'|" : "exp(-k|z-z'|)");
            if (worst / scale > 1e-12) pass = false;
        }
        std::printf("\n");
    }

    // ---- (b) closed form ----
    // For rho(z) = rho0 exp(-k|z|), the exponential kernel integrates to
    //     Phi(z) = -(2 pi G rho0 / k) exp(-k|z|) (1/k + |z|)
    {
        std::printf("(b) against the closed form for rho = rho0 exp(-k|z|)\n");
        std::printf("%8s %14s %10s\n", "Nz", "max rel err", "order");
        const double k = 1.5, rho0 = 4.0e7;
        double prev = 0.0, order = 0.0;
        for (int Nz : {256, 512, 1024, 2048}) {
            Grid2D g(8, Nz, 8.0, 20.0);  // tall, so the tails are captured
            std::vector<Complex> rho(Nz), got(Nz);
            for (int j = 0; j < Nz; ++j) {
                rho[j] = Complex(rho0 * std::exp(-k * std::fabs(g.z(j))), 0.0);
            }
            isolated_column_k(g, k, rho.data(), got.data());
            double worst = 0.0, scale = 0.0;
            for (int j = 0; j < Nz; ++j) {
                const double az = std::fabs(g.z(j));
                if (az > 6.0) continue;
                const double exact = -(2.0 * units::PI * units::G * rho0 / k) *
                                     std::exp(-k * az) * (1.0 / k + az);
                worst = std::max(worst, std::fabs(got[j].real() - exact));
                scale = std::max(scale, std::fabs(exact));
            }
            const double rel = worst / scale;
            if (prev > 0.0) order = std::log2(prev / rel);
            std::printf("%8d %14.3e %10.2f\n", Nz, rel, order);
            prev = rel;
        }
        std::printf("    (second order, limited by the kink in the source at\n"
                    "     z = 0, which is the quadrature rule and not the solve)\n\n");
        // The convergence rate is the meaningful assertion here; the absolute
        // level just tracks the kink and the chosen resolution.
        if (std::fabs(order - 2.0) > 0.15) pass = false;
        if (prev > 1e-4) pass = false;
    }

    // ---- (c) kx = 0 reduces to the 1D isolated solver ----
    {
        Grid2D g2(8, 1024, 8.0, 8.0);
        Grid1D g1(1024, 8.0);
        std::vector<Complex> rho(g2.Nz), got(g2.Nz);
        std::vector<double> rho1(g1.N), phi1;
        for (int j = 0; j < g2.Nz; ++j) {
            rho[j] = Complex(sheet.rho(g2.z(j)), 0.0);
            rho1[j] = sheet.rho(g1.z(j));
        }
        isolated_column_k0(g2, rho.data(), got.data());
        poisson_isolated_1d(g1, rho1, phi1);
        double worst = 0.0, scale = 0.0;
        for (int j = 0; j < g2.Nz; ++j) {
            worst = std::max(worst, std::fabs(got[j].real() - phi1[j]));
            scale = std::max(scale, std::fabs(phi1[j]));
        }
        std::printf("(c) kx = 0 column vs the standalone 1D solver: "
                    "max rel diff %.3e\n\n", worst / scale);
        if (worst / scale > 1e-13) pass = false;
    }

    // ---- (d) full solver, separable density ----
    // rho(x, z) = rho_z(z) [1 + eps cos(k1 x)] has the exact solution
    //   Phi = Phi_0(z) + eps Phi_1(z) cos(k1 x)
    // with Phi_0 the kx = 0 kernel and Phi_1 the exponential kernel at k1.
    {
        std::printf("(d) full solver on a separable density\n");
        std::printf("%8s %8s %16s\n", "Nx", "Nz", "max rel err");
        const double eps = 0.3;
        for (int Nx : {32, 64, 128}) {
            Grid2D g(Nx, 512, 8.0, 6.0);
            const double k1 = g.k_fundamental();
            std::vector<double> rho(g.size()), Phi;
            for (int i = 0; i < g.Nx; ++i) {
                const double mod = 1.0 + eps * std::cos(k1 * g.x(i));
                for (int j = 0; j < g.Nz; ++j) {
                    rho[g.idx(i, j)] = sheet.rho(g.z(j)) * mod;
                }
            }
            PoissonXZ solver(g);
            solver.solve(rho, Phi);

            std::vector<Complex> rz(g.Nz), p0(g.Nz), p1(g.Nz);
            for (int j = 0; j < g.Nz; ++j) rz[j] = Complex(sheet.rho(g.z(j)), 0.0);
            isolated_column_k0(g, rz.data(), p0.data());
            isolated_column_k(g, k1, rz.data(), p1.data());

            double worst = 0.0, scale = 0.0;
            for (int i = 0; i < g.Nx; ++i) {
                const double c = std::cos(k1 * g.x(i));
                for (int j = 0; j < g.Nz; ++j) {
                    if (std::fabs(g.z(j)) > 2.0) continue;
                    const double exact = p0[j].real() + eps * p1[j].real() * c;
                    worst = std::max(worst,
                                     std::fabs(Phi[g.idx(i, j)] - exact));
                    scale = std::max(scale, std::fabs(exact));
                }
            }
            std::printf("%8d %8d %16.3e\n", Nx, g.Nz, worst / scale);
            if (worst / scale > 1e-12) pass = false;
        }
        std::printf("    (exact to round-off: the separable case is reproduced\n"
                    "     mode by mode, so the strided batch transform and the\n"
                    "     FFT normalisation are both right)\n\n");
    }

    std::printf("%s\n", pass ? "PASS" : "FAIL");
    return pass ? 0 : 1;
}
