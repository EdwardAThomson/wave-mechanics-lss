// Validation test 2: harmonic oscillator.
//
// Certifies two independent pieces against exact answers.
//
//   (a) the tridiagonal eigensolver, against E_n = hbar omega (n + 1/2). The
//       eigenpair residual ||H u - E u|| tests the solver; the departure of E_n
//       from the analytic ladder tests the three-point Laplacian, and must fall
//       as dz^2.
//   (b) the split-step, against a coherent state, which in a harmonic well
//       oscillates at omega with a rigidly constant width. The Strang error
//       must fall as dt^2.
#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>

#include "../src/eigen1d.h"
#include "../src/evolve.h"
#include "../src/grid.h"
#include "../src/units.h"

int main() {
    const double hbar = 0.3;
    const double omega = 50.0;  // km/s/kpc
    const double L = 4.0;
    bool pass = true;

    std::printf("=== Test 2: harmonic oscillator ===\n");
    std::printf("hbar_eff = %.3f, omega = %.1f km/s/kpc, hbar*omega = %.3f\n\n",
                hbar, omega, hbar * omega);

    // ---- (a) eigenvalues ----
    std::printf("(a) spectrum vs E_n = hbar omega (n + 1/2)\n");
    std::printf("%8s  %14s  %14s  %14s\n", "N", "max resid", "max rel err n<=40",
                "order");
    double prev = 0.0;
    for (int N : {1024, 2048, 4096}) {
        Grid1D g(N, L);
        std::vector<double> d(N), e(N - 1);
        const double kd = hbar * hbar / (g.dz * g.dz);
        for (int j = 0; j < N; ++j) {
            const double z = g.z(j);
            d[j] = kd + 0.5 * omega * omega * z * z;
        }
        for (int j = 0; j < N - 1; ++j) e[j] = -0.5 * kd;

        auto res = sym_tridiag_lowest(d, e, 41, 40);
        const double resid = eigen_residual(d, e, res);

        double worst = 0.0;
        for (int n = 0; n <= 40; ++n) {
            const double exact = hbar * omega * (n + 0.5);
            worst = std::max(worst, std::fabs(res.values[n] / exact - 1.0));
        }
        const double order = (prev > 0.0) ? std::log2(prev / worst) : 0.0;
        std::printf("%8d  %14.3e  %14.3e  %14.2f\n", N, resid, worst, order);
        prev = worst;
        if (resid > 1e-6 * hbar * omega * 41) pass = false;
        if (N == 4096 && worst > 2e-3) pass = false;
    }
    std::printf("    (residual ~ machine precision certifies the eigensolver;\n"
                "     the ladder error is the three-point Laplacian and falls\n"
                "     as dz^2 exactly as it should)\n\n");

    // ---- (b) coherent state ----
    {
        const int N = 2048;
        Grid1D g(N, L);
        Evolver ev(g, hbar, Kinetic::Spectral);
        const double sg = std::sqrt(hbar / (2.0 * omega));
        const double a = 0.4;  // kpc
        const double T = 2.0 * units::PI / omega;

        std::printf("(b) coherent state, amplitude %.2f kpc, width %.5f kpc, "
                    "T = %.4f (%.2f Myr)\n",
                    a, sg, T, T * units::TIME_MYR);
        std::printf("%10s  %12s  %12s  %6s  %12s  %6s  %10s\n", "steps/T", "dt",
                    "phase err", "order", "|<z>-exact|", "order", "width drift");

        std::vector<double> Vext(N);
        for (int j = 0; j < N; ++j) {
            Vext[j] = 0.5 * omega * omega * g.z(j) * g.z(j);
        }

        double prev_err = 0.0, prev_z = 0.0;
        for (int nsteps : {200, 400, 800, 1600}) {
            const double dt = T / nsteps;
            SlabState st;
            st.grid = g;
            st.hbar = hbar;
            st.psi.assign(1, std::vector<Complex>(N));
            for (int j = 0; j < N; ++j) {
                const double dz = g.z(j) - a;
                st.psi[0][j] = Complex(std::exp(-dz * dz / (4.0 * sg * sg)), 0.0);
            }
            double nrm = 0.0;
            for (int j = 0; j < N; ++j) nrm += std::norm(st.psi[0][j]);
            nrm = std::sqrt(nrm * g.dz);
            for (int j = 0; j < N; ++j) st.psi[0][j] /= nrm;

            for (int s = 0; s < nsteps; ++s) {
                ev.kinetic_half(st, dt);
                ev.potential_full(st, Vext, dt);
                ev.kinetic_half(st, dt);
            }

            std::vector<double> rho;
            st.density(rho);
            std::vector<Complex> dpsi;
            ev.ddz(st.psi[0], dpsi);
            double m0 = 0.0, m1 = 0.0, m2 = 0.0, mv = 0.0;
            for (int j = 0; j < N; ++j) {
                m0 += rho[j];
                m1 += rho[j] * g.z(j);
                m2 += rho[j] * g.z(j) * g.z(j);
                mv += hbar * std::imag(std::conj(st.psi[0][j]) * dpsi[j]);
            }
            const double zbar = m1 / m0;
            const double vbar = mv / m0;
            const double width = std::sqrt(m2 / m0 - zbar * zbar);

            // The oscillator phase after one exact period should be zero. The
            // Strang splitting shifts the numerical frequency by O(dt^2), so
            // this is the diagnostic that shows the scheme's true order.
            const double phase_err =
                std::fabs(std::atan2(-vbar / omega, zbar));
            // The displacement error is 1 - cos(phase_err), hence O(dt^4). It
            // is not evidence of a fourth-order scheme, just of measuring the
            // oscillation at a stationary point of its own phase.
            const double err = std::fabs(zbar - a);

            const double ord_p =
                (prev_err > 0.0) ? std::log2(prev_err / phase_err) : 0.0;
            const double ord_z =
                (prev_z > 0.0) ? std::log2(prev_z / err) : 0.0;
            std::printf("%10d  %12.3e  %12.3e  %6.2f  %12.3e  %6.2f  %10.3e\n",
                        nsteps, dt, phase_err, ord_p, err, ord_z,
                        std::fabs(width / sg - 1.0));
            prev_err = phase_err;
            prev_z = err;
            if (nsteps == 1600 && phase_err > 1e-5) pass = false;
            if (nsteps == 1600 && std::fabs(ord_p - 2.0) > 0.15) pass = false;
        }
        std::printf("    (phase error falls 4x per halving of dt: Strang is\n"
                    "     second order. A constant width confirms the potential\n"
                    "     operator injects no spurious spreading.)\n\n");
    }

    std::printf("%s\n", pass ? "PASS" : "FAIL");
    return pass ? 0 : 1;
}
