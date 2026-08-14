// Stage 1: the firehose branch, which is the part of test 5 that an unsheared
// box CAN reach.
//
// The in-plane Jeans scan established that a non-rotating sheet fragments for
// k < k_J = 2 pi G Sigma / sigma_x^2. The razor-thin bending relation is
//
//     omega^2 = 2 pi G Sigma k - k^2 sigma_x^2
//
// which is stable for exactly the same k < k_J and firehose-unstable above it.
// So in an unsheared box the two conditions are complementary: wherever
// bending waves would propagate stably, the sheet is gravitationally unstable
// in plane, and vice versa. The stable bending branch is unreachable without
// rotation, but the FIREHOSE branch is reachable, because it lives at k > k_J
// where the sheet is in-plane stable.
//
// That is what this measures. Choose sigma_x above the box's in-plane
// stability threshold, so every mode the box holds has k > k_J, displace the
// sheet, and watch. Razor-thin theory says every mode should grow at
//
//     gamma = sqrt(k^2 sigma_x^2 - 2 pi G Sigma k).
//
// Finite thickness is known to stabilise the firehose strongly, so the
// interesting question is how much of that predicted growth survives at
// k h of order unity.
#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>

#include "../src/diagnostics2d.h"
#include "../src/equilibrium.h"
#include "../src/evolve2d.h"
#include "../src/grid2d.h"
#include "../src/poisson_xz.h"
#include "../src/slab_ic_2d.h"
#include "../src/units.h"

int main() {
    const double Sigma_pc2 = 50.0, sigma_z = 20.0, hbar = 0.6;
    const double sigma_x = 45.0;  // above sqrt(G Sigma Lx) = 29.3, so the
                                  // sheet is in-plane stable and cannot
                                  // fragment during the measurement
    const double Sigma = Sigma_pc2 * units::MSUN_PC2;
    const int Nx = 1024, Nz = 512, n_streams = 4;
    const double Lx = 4.0, Lz = 5.0;
    const int nmodes = 4;

    Grid2D g(Nx, Nz, Lx, Lz);
    Grid1D g1(Nz, Lz);
    SheetConfig cfg;
    cfg.Sigma = Sigma;
    cfg.sigma_z = sigma_z;
    cfg.hbar = hbar;
    cfg.nu_ext = 0.0;
    cfg.n_states = 300;
    cfg.verbose = false;
    WarmSheet sheet = build_warm_sheet(g1, cfg);

    Evolver2D ev(g, hbar, Kinetic::Spectral);
    PoissonXZ poisson(g);
    SlabState2D st = make_warm_2d(sheet, g, sigma_x, n_streams, 60607ULL);

    const double h = sheet.h_measured;
    const double k_J = 2.0 * units::PI * units::G * Sigma / (sigma_x * sigma_x);
    std::printf("=== Stage 1: firehose branch ===\n");
    std::printf("grid %d x %d, Lx = %.1f kpc, sigma_x = %.1f km/s, "
                "h = %.4f kpc\n", Nx, Nz, Lx, sigma_x, h);
    std::printf("k_J = %.3f /kpc, k_min = %.3f /kpc: every mode in the box is\n"
                "above k_J, so the sheet is in-plane stable and every mode is\n"
                "firehose-unstable in razor-thin theory.\n",
                k_J, 2.0 * units::PI / Lx);
    std::printf("v ceiling x = %.1f (%.1f sigma_x)\n\n", ev.v_ceiling_x(),
                ev.v_ceiling_x() / sigma_x);

    // Small displacement so the response stays linear.
    const double h0 = 0.01;
    std::vector<double> shift(g.Nx, 0.0);
    for (int i = 0; i < g.Nx; ++i) {
        for (int n = 1; n <= nmodes; ++n) {
            shift[i] += h0 * std::cos(g.kx(n) * g.x(i));
        }
    }
    ev.displace_columns(st, shift);

    // Razor-thin growth is fast, so run only a fraction of a vertical period.
    const double T = 2.0 * units::PI / sheet.omega0;
    const double dt = T / 400.0;
    const int nsteps = 220;
    const int sample_every = 5;

    std::vector<double> tser;
    std::vector<std::vector<double>> amp(nmodes), re(nmodes);
    std::vector<double> rho, V;
    for (int step = 0; step <= nsteps; ++step) {
        if (step % sample_every == 0) {
            ColumnMoments m = compute_columns_2d(ev, st);
            tser.push_back(step * dt);
            for (int n = 1; n <= nmodes; ++n) {
                const Complex a = mode_amplitude(g, m.zbar, n);
                amp[n - 1].push_back(std::abs(a));
                re[n - 1].push_back(a.real());
            }
        }
        if (step == nsteps) break;
        ev.kinetic_half(st, dt);
        st.density(rho);
        poisson.solve(rho, V);
        ev.potential_full(st, V, dt);
        ev.kinetic_half(st, dt);
    }

    std::printf("%5s %9s %7s %13s %13s %9s %10s\n", "mode", "k[1/kpc]", "k h",
                "gamma_thin", "gamma_meas", "ratio", "sign flips");
    for (int n = 1; n <= nmodes; ++n) {
        const double k = g.kx(n);
        const double w2 = 2.0 * units::PI * units::G * Sigma * k -
                          k * k * sigma_x * sigma_x;
        const double g_thin = (w2 < 0.0) ? std::sqrt(-w2) : 0.0;

        // Log slope of |amplitude|, least squares.
        const int m = static_cast<int>(amp[n - 1].size());
        double sa = 0, sb = 0, sc = 0, sd = 0;
        int used = 0;
        for (int i = 0; i < m; ++i) {
            if (amp[n - 1][i] <= 0.0) continue;
            const double y = std::log(amp[n - 1][i]);
            sa += tser[i]; sb += y; sc += tser[i] * y; sd += tser[i] * tser[i];
            ++used;
        }
        const double g_meas =
            (used > 2) ? (used * sc - sa * sb) / (used * sd - sa * sa) : 0.0;

        // An oscillating (stable) mode changes the sign of Re(amplitude);
        // a growing one does not.
        int flips = 0;
        for (int i = 1; i < m; ++i) {
            if (re[n - 1][i - 1] * re[n - 1][i] < 0.0) ++flips;
        }
        std::printf("%5d %9.3f %7.3f %13.2f %13.2f %9.3f %10d\n", n, k, k * h,
                    g_thin, g_meas, (g_thin > 0.0) ? g_meas / g_thin : 0.0,
                    flips);
    }

    std::printf("\n  gamma_thin is the razor-thin firehose growth rate. A mode\n"
                "  that is actually growing shows gamma_meas > 0 and no sign\n"
                "  flips; one that finite thickness has stabilised oscillates\n"
                "  instead, flipping sign and returning gamma_meas near zero.\n");
    return 0;
}
