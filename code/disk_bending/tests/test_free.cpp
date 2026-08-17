// Validation test 1: free particle.
//
// Checks FFT conventions, the sign and normalisation of the kinetic operator,
// and the phase accuracy of the split-step, with gravity switched off.
//
//   (a) plane wave: psi ~ exp(i k z) must acquire exactly exp(-i hbar k^2 t / 2)
//   (b) Gaussian wavepacket: sigma(t)^2 = sigma0^2 + (hbar t / 2 sigma0)^2
//   (c) MatchedFD reproduces the finite-difference dispersion, not the exact
//       one, which is the whole point of offering it
#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>

#include "../src/evolve.h"
#include "../src/grid.h"
#include "../src/units.h"

int main() {
    const double L = 4.0;
    const int N = 1024;
    const double hbar = 0.3;
    Grid1D g(N, L);
    std::vector<double> no_ext(N, 0.0);
    PotentialWorkspace ws;
    bool pass = true;

    std::printf("=== Test 1: free particle ===\n");
    std::printf("N = %d, L = %.1f kpc, hbar_eff = %.3f kpc km/s\n\n", N, L, hbar);

    // ---- (a) plane wave ----
    {
        Evolver ev(g, hbar, Kinetic::Spectral);
        std::printf("(a) plane-wave phase after t = 0.01 (%.2f Myr)\n",
                    0.01 * units::TIME_MYR);
        std::printf("%6s  %12s  %14s  %12s\n", "mode", "k[1/kpc]",
                    "phase err[rad]", "v=hbar k");
        for (int n : {1, 4, 16, 64, 200}) {
            const double k = 2.0 * units::PI * n / L;
            SlabState st;
            st.grid = g;
            st.hbar = hbar;
            st.psi.assign(1, std::vector<Complex>(N));
            for (int j = 0; j < N; ++j) {
                const double ph = k * g.z(j);
                st.psi[0][j] = Complex(std::cos(ph), std::sin(ph));
            }
            const double t = 0.01;
            const int steps = 100;
            const double dt = t / steps;
            for (int s = 0; s < steps; ++s) {
                ev.kinetic_half(st, dt);
                ev.kinetic_half(st, dt);
            }
            const double expect = -0.5 * hbar * k * k * t;
            double got = std::arg(st.psi[0][N / 3] /
                                  Complex(std::cos(k * g.z(N / 3)),
                                          std::sin(k * g.z(N / 3))));
            double err = std::fmod(got - expect, 2.0 * units::PI);
            while (err > units::PI) err -= 2.0 * units::PI;
            while (err < -units::PI) err += 2.0 * units::PI;
            std::printf("%6d  %12.4f  %14.3e  %12.3f\n", n, k, std::fabs(err),
                        hbar * k);
            if (std::fabs(err) > 1e-10) pass = false;
        }
        std::printf("    v_max on this grid = %.2f km/s\n\n", ev.v_max());
    }

    // ---- (b) Gaussian spreading ----
    {
        Evolver ev(g, hbar, Kinetic::Spectral);
        const double w0 = 0.10;
        std::printf("(b) Gaussian spreading, sigma0 = %.3f kpc\n", w0);
        std::printf("%10s  %14s  %14s  %10s\n", "t", "sigma_num", "sigma_exact",
                    "rel");
        SlabState st;
        st.grid = g;
        st.hbar = hbar;
        st.psi.assign(1, std::vector<Complex>(N));
        for (int j = 0; j < N; ++j) {
            const double z = g.z(j);
            st.psi[0][j] = Complex(std::exp(-z * z / (4.0 * w0 * w0)), 0.0);
        }
        // normalise to unit surface density
        double nrm = 0.0;
        for (int j = 0; j < N; ++j) nrm += std::norm(st.psi[0][j]);
        nrm = std::sqrt(nrm * g.dz);
        for (int j = 0; j < N; ++j) st.psi[0][j] /= nrm;

        const double dt = 2.0e-5;
        int done = 0;
        for (double t_target : {0.005, 0.010, 0.020}) {
            const int want = static_cast<int>(t_target / dt);
            for (; done < want; ++done) {
                ev.kinetic_half(st, dt);
                ev.kinetic_half(st, dt);
            }
            const double t = done * dt;
            std::vector<double> rho;
            st.density(rho);
            double m0 = 0.0, m2 = 0.0;
            for (int j = 0; j < N; ++j) {
                m0 += rho[j];
                m2 += rho[j] * g.z(j) * g.z(j);
            }
            const double sig = std::sqrt(m2 / m0);
            const double exact =
                std::sqrt(w0 * w0 + std::pow(hbar * t / (2.0 * w0), 2));
            const double rel = std::fabs(sig / exact - 1.0);
            std::printf("%10.4f  %14.6f  %14.6f  %10.2e\n", t, sig, exact, rel);
            if (rel > 1e-6) pass = false;
        }
        std::printf("\n");
    }

    // ---- (c) matched finite-difference dispersion ----
    {
        Evolver ev(g, hbar, Kinetic::MatchedFD);
        std::printf("(c) MatchedFD kinetic operator: measured vs FD dispersion\n");
        std::printf("%6s  %10s  %14s  %14s  %10s\n", "mode", "k dz",
                    "omega_measured", "omega_fd", "omega_exact");
        for (int n : {1, 16, 64, 200, 400}) {
            const double k = 2.0 * units::PI * n / L;
            SlabState st;
            st.grid = g;
            st.hbar = hbar;
            st.psi.assign(1, std::vector<Complex>(N));
            for (int j = 0; j < N; ++j) {
                const double ph = k * g.z(j);
                st.psi[0][j] = Complex(std::cos(ph), std::sin(ph));
            }
            // Short enough that even the Nyquist-adjacent mode accumulates
            // less than pi of phase, so arg() does not wrap.
            const double t = 2.0e-5;
            ev.kinetic_half(st, t);
            ev.kinetic_half(st, t);
            double got = std::arg(st.psi[0][N / 3] /
                                  Complex(std::cos(k * g.z(N / 3)),
                                          std::sin(k * g.z(N / 3))));
            const double omega_meas = -got / t;
            const double s = std::sin(0.5 * k * g.dz);
            const double omega_fd = 0.5 * hbar * 4.0 * s * s / (g.dz * g.dz);
            const double omega_ex = 0.5 * hbar * k * k;
            std::printf("%6d  %10.4f  %14.4f  %14.4f  %10.4f\n", n, k * g.dz,
                        omega_meas, omega_fd, omega_ex);
            if (std::fabs(omega_meas / omega_fd - 1.0) > 1e-9) pass = false;
        }
        std::printf("    v_max (MatchedFD) = %.2f km/s vs %.2f km/s spectral;\n"
                    "    the FD group velocity turns over at k dz = pi/2, so its\n"
                    "    usable velocity ceiling is lower by a factor pi.\n\n",
                    ev.v_max(), hbar * g.k_nyquist());
    }

    std::printf("%s\n", pass ? "PASS" : "FAIL");
    return pass ? 0 : 1;
}
