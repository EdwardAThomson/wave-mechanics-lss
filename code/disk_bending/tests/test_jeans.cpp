// Validation test 3: Jeans growth in a periodic box.
//
// This is the test that catches a sign error in the Poisson coupling, which is
// otherwise easy to get wrong and produces plausible-looking nonsense. It runs
// periodic on purpose: an unstable uniform medium is the one case where
// periodicity is the physically correct boundary condition, so it isolates the
// gravity coupling from the slab boundary work.
//
// Linearised Schrodinger-Poisson about a uniform background gives
//
//     omega^2 = hbar^2 k^4 / 4 - 4 pi G rho0
//
// so the growth rate is gamma = sqrt(4 pi G rho0 - hbar^2 k^4 / 4) below the
// quantum Jeans wavenumber and the mode oscillates above it. Both branches are
// checked, which also pins the sign of the quantum-pressure term.
#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>

#include "../src/evolve.h"
#include "../src/fft.h"
#include "../src/grid.h"
#include "../src/poisson_z.h"
#include "../src/units.h"

namespace {

// Amplitude of the k = 2 pi / L mode of the density contrast.
double mode_amplitude(const Grid1D& g, const SlabState& st) {
    std::vector<double> rho;
    st.density(rho);
    double mean = 0.0;
    for (int j = 0; j < g.N; ++j) mean += rho[j];
    mean /= g.N;
    double c = 0.0, s = 0.0;
    const double k = 2.0 * units::PI / g.L;
    for (int j = 0; j < g.N; ++j) {
        c += (rho[j] - mean) * std::cos(k * g.z(j));
        s += (rho[j] - mean) * std::sin(k * g.z(j));
    }
    return 2.0 * std::hypot(c, s) / (g.N * mean);
}

struct Result {
    double measured;
    double predicted;
};

Result run(double L, int N, double hbar, double rho0, double t_end, int nsteps,
           bool growing) {
    Grid1D g(N, L);
    Evolver ev(g, hbar, Kinetic::Spectral);
    FFT1D fft(N);

    const double k = 2.0 * units::PI / L;
    const double eps = 1e-4;

    SlabState st;
    st.grid = g;
    st.hbar = hbar;
    st.psi.assign(1, std::vector<Complex>(N));
    for (int j = 0; j < N; ++j) {
        const double amp = std::sqrt(rho0 * (1.0 + eps * std::cos(k * g.z(j))));
        st.psi[0][j] = Complex(amp, 0.0);
    }

    const double dt = t_end / nsteps;
    std::vector<double> rho, Phi;
    const double d0 = mode_amplitude(g, st);

    for (int s = 0; s < nsteps; ++s) {
        ev.kinetic_half(st, dt);
        st.density(rho);
        poisson_periodic_1d(g, rho, Phi, fft);
        ev.potential_full(st, Phi, dt);
        ev.kinetic_half(st, dt);
    }

    const double d1 = mode_amplitude(g, st);
    const double w2 = 0.25 * hbar * hbar * std::pow(k, 4) -
                      4.0 * units::PI * units::G * rho0;

    Result r;
    if (growing) {
        const double gamma = std::sqrt(-w2);
        // A static initial perturbation splits equally into the growing and
        // decaying branches, so delta(t) = delta0 cosh(gamma t).
        r.predicted = gamma;
        r.measured = std::acosh(d1 / d0) / t_end;
    } else {
        const double omega = std::sqrt(w2);
        r.predicted = omega;
        // delta(t) = delta0 cos(omega t)
        r.measured = std::acos(std::max(-1.0, std::min(1.0, d1 / d0))) / t_end;
    }
    return r;
}

}  // namespace

int main() {
    bool pass = true;
    std::printf("=== Test 3: Jeans growth (periodic) ===\n\n");

    // ---- growing branch: gravity dominates ----
    {
        const double rho0 = 4.0e7;  // Msun/kpc^3
        const double L = 2.0, hbar = 0.3;
        const double k = 2.0 * units::PI / L;
        const double t_end = 0.03;
        std::printf("(a) growing branch: rho0 = %.2e, L = %.1f kpc, "
                    "hbar = %.2f\n",
                    rho0, L, hbar);
        std::printf("    4 pi G rho0 = %.4e, hbar^2 k^4 / 4 = %.4e\n",
                    4.0 * units::PI * units::G * rho0,
                    0.25 * hbar * hbar * std::pow(k, 4));
        std::printf("%10s  %14s  %14s  %10s\n", "steps", "gamma_meas",
                    "gamma_pred", "rel");
        for (int nsteps : {2000, 4000, 8000}) {
            auto r = run(L, 512, hbar, rho0, t_end, nsteps, true);
            const double rel = std::fabs(r.measured / r.predicted - 1.0);
            std::printf("%10d  %14.4f  %14.4f  %10.2e\n", nsteps, r.measured,
                        r.predicted, rel);
            if (nsteps == 8000 && rel > 5e-3) pass = false;
        }
        std::printf("\n");
    }

    // ---- oscillating branch: quantum pressure dominates ----
    {
        const double rho0 = 4.0e7;
        const double L = 0.2, hbar = 3.0;
        const double k = 2.0 * units::PI / L;
        const double w2 = 0.25 * hbar * hbar * std::pow(k, 4) -
                          4.0 * units::PI * units::G * rho0;
        const double t_end = 0.25 * (2.0 * units::PI / std::sqrt(w2));
        std::printf("(b) oscillating branch: L = %.1f kpc, hbar = %.1f\n", L,
                    hbar);
        std::printf("    4 pi G rho0 = %.4e, hbar^2 k^4 / 4 = %.4e\n",
                    4.0 * units::PI * units::G * rho0,
                    0.25 * hbar * hbar * std::pow(k, 4));
        std::printf("%10s  %14s  %14s  %10s\n", "steps", "omega_meas",
                    "omega_pred", "rel");
        for (int nsteps : {4000, 8000, 16000}) {
            auto r = run(L, 512, hbar, rho0, t_end, nsteps, false);
            const double rel = std::fabs(r.measured / r.predicted - 1.0);
            std::printf("%10d  %14.4f  %14.4f  %10.2e\n", nsteps, r.measured,
                        r.predicted, rel);
            if (nsteps == 16000 && rel > 5e-3) pass = false;
        }
        std::printf("    (this branch has the opposite sign of the quantum term,\n"
                    "     so agreeing on both pins the sign of BOTH the gravity\n"
                    "     coupling and the kinetic operator)\n\n");
    }

    std::printf("%s\n", pass ? "PASS" : "FAIL");
    return pass ? 0 : 1;
}
