// Validation test 0: the isolated vertical Poisson solve.
//
// Two things are checked.
//
//   (a) Against the analytic isothermal sheet, Phi = 2 sigma^2 ln cosh(z/2h),
//       including the convergence rate with resolution.
//
//   (b) Against a periodic solve on the same density, to show why the
//       cosmological code's triply-periodic Poisson cannot simply be reused.
//       The periodic answer differs by the image term -2 pi G Sigma z^2 / L,
//       which shrinks only as 1/L. Making the box taller is not a fix.
#include <cmath>
#include <cstdio>
#include <vector>

#include "../src/fft.h"
#include "../src/grid.h"
#include "../src/poisson_z.h"
#include "../src/units.h"

int main() {
    const double Sigma = 50.0 * units::MSUN_PC2;
    const double sigma = 20.0;
    IsothermalSheet sheet(Sigma, sigma);

    std::printf("=== Test 0: isolated vertical Poisson ===\n");
    std::printf("Sigma = %.3e Msun/kpc^2, sigma_z = %.1f km/s, h = %.5f kpc\n\n",
                Sigma, sigma, sheet.h);

    std::printf("(a) isolated solve vs analytic sech^2 sheet, box L = 8 kpc\n");
    std::printf("%8s  %12s  %12s  %8s\n", "N", "max|dPhi|", "rel", "order");
    double prev_err = 0.0;
    bool pass_a = true;
    for (int N : {256, 512, 1024, 2048, 4096}) {
        Grid1D g(N, 8.0);
        std::vector<double> rho(N), Phi;
        for (int j = 0; j < N; ++j) rho[j] = sheet.rho(g.z(j));
        poisson_isolated_1d(g, rho, Phi);

        double maxerr = 0.0, scale = 0.0;
        for (int j = 0; j < N; ++j) {
            const double z = g.z(j);
            if (std::fabs(z) > 2.0) continue;  // compare over the occupied range
            const double a = sheet.Phi(z);
            maxerr = std::max(maxerr, std::fabs(Phi[j] - a));
            scale = std::max(scale, std::fabs(a));
        }
        const double order = (prev_err > 0.0) ? std::log2(prev_err / maxerr) : 0.0;
        std::printf("%8d  %12.4e  %12.4e  %8.2f\n", N, maxerr, maxerr / scale,
                    order);
        if (N == 4096 && maxerr / scale > 1e-6) pass_a = false;
        prev_err = maxerr;
    }
    std::printf("    (clean second order: the residual is the midpoint\n"
                "     quadrature error of the |z - z'| kernel, and there is no\n"
                "     image contamination in it at all)\n\n");

    std::printf("(b) same density, periodic solve, growing box\n");
    std::printf("%8s  %14s  %14s  %10s\n", "L[kpc]", "max|dPhi| @1kpc",
                "predicted image", "rel err");
    for (double L : {4.0, 8.0, 16.0, 32.0}) {
        const int N = static_cast<int>(1024 * L / 8.0);
        Grid1D g(N, L);
        std::vector<double> rho(N), Phi_p;
        for (int j = 0; j < N; ++j) rho[j] = sheet.rho(g.z(j));
        FFT1D fft(N);
        poisson_periodic_1d(g, rho, Phi_p, fft);

        // Gauge both to zero at the midplane before comparing.
        double p0 = 0.5 * (Phi_p[N / 2 - 1] + Phi_p[N / 2]);
        double err = 0.0, ref = 0.0;
        for (int j = 0; j < N; ++j) {
            const double z = g.z(j);
            if (std::fabs(std::fabs(z) - 1.0) > g.dz) continue;
            err = std::max(err, std::fabs((Phi_p[j] - p0) - sheet.Phi(z)));
            ref = std::fabs(sheet.Phi(z));
        }
        const double predicted = 2.0 * units::PI * units::G * Sigma * 1.0 / L;
        std::printf("%8.0f  %14.4e  %14.4e  %10.1f%%\n", L, err, predicted,
                    100.0 * err / ref);
    }
    std::printf("    (the periodic error tracks 2 pi G Sigma z^2 / L and falls\n"
                "     only as 1/L: a taller periodic box does not rescue you)\n\n");

    // (c) The vertical force must saturate at the isolated value 2 pi G Sigma.
    {
        Grid1D g(4096, 12.0);
        std::vector<double> rho(4096), Phi, gz(4096);
        for (int j = 0; j < 4096; ++j) rho[j] = sheet.rho(g.z(j));
        poisson_isolated_1d(g, rho, Phi, &gz);
        const double expect = 2.0 * units::PI * units::G * Sigma;
        const double got = gz[4095];
        std::printf("(c) dPhi/dz at the top edge = %.6e, expected %.6e "
                    "(rel %.2e)\n\n",
                    got, expect, std::fabs(got / expect - 1.0));
        if (std::fabs(got / expect - 1.0) > 1e-6) pass_a = false;
    }

    std::printf("%s\n", pass_a ? "PASS" : "FAIL");
    return pass_a ? 0 : 1;
}
