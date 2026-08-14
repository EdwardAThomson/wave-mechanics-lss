// Stage 1: the heating is in-plane gravitational instability, not granule
// relaxation.
//
// The control experiment (test_granule2) showed the heating needs a
// time-dependent self-consistent potential and is not aliasing. But it does
// not fall with stream count, which rules out the stochastic-relaxation
// picture: seed fluctuations fall as 1/sqrt(N_streams), so a noise-driven
// process would fall as 1/N_streams.
//
// What is independent of the seed amplitude is an INSTABILITY. It grows
// exponentially from whatever seed is there, so halving the seed only delays
// it logarithmically.
//
// And a non-rotating self-gravitating sheet IS in-plane unstable. For a
// razor-thin sheet the marginal wavenumber is
//
//     k_J = 2 pi G Sigma / sigma_x^2
//
// with everything at k < k_J unstable. There is no rotation and no shear in a
// 2D (x, z) box, so nothing supplies the Coriolis support that keeps a real
// disk at Toomre Q > 1. The box is stable only if even its longest mode is
// stable, k_min = 2 pi / Lx > k_J, that is
//
//     sigma_x^2 > 2 pi G Sigma / k_min = G Sigma Lx.
//
// This program tests that directly by scanning sigma_x across the predicted
// threshold. If the heating switches off where the box becomes stable, the
// diagnosis is confirmed and the consequence is important: a 2D box with no
// shear cannot represent a stable disk, so the plan's "skip shear entirely at
// this stage" does not hold.
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

namespace {

// Power in the surface-density contrast at in-plane mode n.
double sigma_mode_power(const Grid2D& g, const std::vector<double>& Sig, int n) {
    double mean = 0.0;
    for (int i = 0; i < g.Nx; ++i) mean += Sig[i];
    mean /= g.Nx;
    Complex acc(0.0, 0.0);
    for (int i = 0; i < g.Nx; ++i) {
        const double ph = -g.kx(n) * g.x(i);
        acc += (Sig[i] / mean - 1.0) * Complex(std::cos(ph), std::sin(ph));
    }
    return std::abs(acc) / g.Nx;
}

}  // namespace

int main() {
    const double Sigma_pc2 = 50.0, sigma_z = 20.0, hbar = 0.6;
    const double Sigma = Sigma_pc2 * units::MSUN_PC2;
    const int Nx = 1024, Nz = 512;
    const double Lx = 4.0, Lz = 5.0;
    const int n_streams = 4, n_per = 2, spp = 200;

    Grid2D g(Nx, Nz, Lx, Lz);
    Grid1D g1(Nz, Lz);

    const double k_min = 2.0 * units::PI / Lx;
    const double sigma_crit = std::sqrt(units::G * Sigma * Lx);

    std::printf("=== Stage 1: in-plane gravitational instability ===\n");
    std::printf("grid %d x %d, Lx = %.1f, Lz = %.1f kpc, %d streams\n", Nx, Nz,
                Lx, Lz, n_streams);
    std::printf("Sigma = %.3e Msun/kpc^2, k_min = 2 pi / Lx = %.3f /kpc\n",
                Sigma, k_min);
    std::printf("box is in-plane stable when sigma_x > sqrt(G Sigma Lx) = "
                "%.2f km/s\n", sigma_crit);
    std::printf("(razor-thin estimate; finite thickness weakens self-gravity "
                "and lowers it a little)\n\n");

    std::printf("%9s %9s %10s %10s %12s %12s %10s %10s\n", "sigma_x", "k_J",
                "unstable?", "n_unstab", "dsx/period", "dsz/period",
                "dP1/P1", "max spill");

    for (double sigma_x : {10.0, 20.0, 30.0, 45.0, 60.0}) {
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
        SlabState2D st = make_warm_2d(sheet, g, sigma_x, n_streams, 31337ULL);

        const double k_J = 2.0 * units::PI * units::G * Sigma / (sigma_x * sigma_x);
        const int n_unstable = static_cast<int>(std::floor(k_J / k_min));

        const double T = 2.0 * units::PI / sheet.omega0;
        const double dt = T / spp;

        ColumnMoments m0 = compute_columns_2d(ev, st);
        // Full moments for the dispersions.
        ColumnMoments f0 = compute_moments_2d(ev, st, poisson);
        const double P1_0 = sigma_mode_power(g, m0.Sigma, 1);

        std::vector<double> rho, V;
        for (int step = 0; step < n_per * spp; ++step) {
            ev.kinetic_half(st, dt);
            st.density(rho);
            poisson.solve(rho, V);
            ev.potential_full(st, V, dt);
            ev.kinetic_half(st, dt);
        }

        ColumnMoments m1 = compute_columns_2d(ev, st);
        ColumnMoments f1 = compute_moments_2d(ev, st, poisson);
        const double P1_1 = sigma_mode_power(g, m1.Sigma, 1);
        double spx = 0.0, spz = 0.0;
        ev.spectral_spill(st, 0.5 * g.kx_nyquist(), 0.5 * g.kz_nyquist(), spx,
                          spz);

        std::printf("%9.1f %9.3f %10s %10d %12.4f %12.4f %10.2f %10.2e\n",
                    sigma_x, k_J, (k_J > k_min) ? "yes" : "no", n_unstable,
                    (f1.sigma_x_rms - f0.sigma_x_rms) / n_per,
                    (f1.sigma_z_rms - f0.sigma_z_rms) / n_per,
                    (P1_0 > 0.0) ? P1_1 / P1_0 : 0.0, spx);
    }

    std::printf("\n  'n_unstab' counts the in-plane modes the box holds below\n"
                "  k_J. 'dP1/P1' is the growth of the longest surface-density\n"
                "  mode over the run: greater than one means it is growing.\n\n"
                "  If the heating switches off as sigma_x crosses the\n"
                "  threshold, the in-plane heating is gravitational\n"
                "  fragmentation of an unsupported sheet, not granule\n"
                "  relaxation, and the fix is rotation, not resolution.\n");
    return 0;
}
