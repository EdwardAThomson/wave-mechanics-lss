// Stage 1: granule heating, and how it scales with the number of streams.
//
// This is the risk flagged in the Stage 0 write-up. The in-plane part of a 2D
// initial condition can only be a random-phase superposition, so |psi|^2
// carries interference granules of scale lambda_dB = 2 pi hbar_eff / sigma_x,
// and those granules self-gravitate. In the fuzzy-dark-matter literature that
// drives spurious two-body-like relaxation. Here it would heat the disk and
// corrupt exactly the quantities Stage 1 exists to measure: the firehose
// threshold and the Landau damping rate both depend on the dispersion.
//
// The scaling matters more than the absolute rate. Summing T independent
// realisations gives a density contrast falling as 1/sqrt(T), and the heating
// is driven by the square of that contrast, so the prediction is
//
//     heating rate  proportional to  1 / N_streams
//
// If that holds, 2D is safe because streams are cheap there. It is also the
// sharpest possible warning about 3D, where a single wavefunction is the only
// affordable option and the heating is therefore at its maximum.
#include <cmath>
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
    const double Sigma_pc2 = 50.0, sigma_z = 20.0, sigma_x = 10.0, hbar = 0.6;
    // Headroom set so that the Maxwellian tail does not alias even after the
    // disk has heated: 12 sigma in each direction at the start.
    const int Nx = 512, Nz = 768;
    const double Lx = 8.0, Lz = 6.0;

    Grid2D g(Nx, Nz, Lx, Lz);
    Grid1D g1(Nz, Lz);
    SheetConfig cfg;
    cfg.Sigma = Sigma_pc2 * units::MSUN_PC2;
    cfg.sigma_z = sigma_z;
    cfg.hbar = hbar;
    cfg.nu_ext = 0.0;
    cfg.n_states = 300;
    cfg.verbose = false;
    WarmSheet sheet = build_warm_sheet(g1, cfg);

    Evolver2D ev(g, hbar, Kinetic::Spectral);
    PoissonXZ poisson(g);

    const double T = 2.0 * units::PI / sheet.omega0;
    const int n_per = 2, spp = 250;
    const int samples_per_period = 4;
    const double dt = T / spp;

    std::printf("=== Stage 1: granule heating vs stream count ===\n");
    std::printf("grid %d x %d, Lx = %.1f, Lz = %.1f kpc\n", Nx, Nz, Lx, Lz);
    std::printf("sigma_z = %.1f, sigma_x = %.1f km/s, hbar_eff = %.2f\n",
                sigma_z, sigma_x, hbar);
    std::printf("in-plane lambda_dB = %.4f kpc = %.2f h,  h = %.4f kpc\n",
                2.0 * units::PI * hbar / sigma_x,
                2.0 * units::PI * hbar / sigma_x / sheet.h_measured,
                sheet.h_measured);
    std::printf("v ceiling: x %.1f (%.1f sigma_x), z %.1f (%.1f sigma_z)\n",
                ev.v_ceiling_x(), ev.v_ceiling_x() / sigma_x, ev.v_ceiling_z(),
                ev.v_ceiling_z() / sigma_z);
    std::printf("evolving %d vertical periods (%.0f Myr) at %d steps each\n\n",
                n_per, n_per * T * units::TIME_MYR, spp);

    std::printf("%8s %12s %12s %12s %12s %10s %10s\n", "streams",
                "rms drho/rho", "dsx/period", "dsz/period", "rate*N", "dE/E",
                "max spill_x");

    double rate_ref = 0.0;
    for (int n_streams : {2, 8, 32}) {
        SlabState2D st = make_warm_2d(sheet, g, sigma_x, n_streams,
                                      4242001ULL + 97ULL * n_streams);

        // Initial density contrast, measured near the midplane where the
        // granules live.
        double contrast = 0.0;
        {
            std::vector<double> rho;
            st.density(rho);
            double mean = 0.0;
            int cnt = 0;
            for (int i = 0; i < Nx; ++i) {
                for (int j = 0; j < Nz; ++j) {
                    if (std::fabs(g.z(j)) > 0.5 * sheet.h_measured) continue;
                    mean += rho[g.idx(i, j)];
                    ++cnt;
                }
            }
            mean /= cnt;
            double s2 = 0.0;
            for (int i = 0; i < Nx; ++i) {
                for (int j = 0; j < Nz; ++j) {
                    if (std::fabs(g.z(j)) > 0.5 * sheet.h_measured) continue;
                    const double d = rho[g.idx(i, j)] / mean - 1.0;
                    s2 += d * d;
                }
            }
            contrast = std::sqrt(s2 / cnt);
        }

        ColumnMoments m0 = compute_moments_2d(ev, st, poisson);
        std::vector<double> ts, sx, sz;
        double dE = 0.0, spx = 0.0, spz = 0.0, spill_max = 0.0;
        const int chunk = spp / samples_per_period;
        for (int p = 0; p < n_per * samples_per_period; ++p) {
            for (int k = 0; k < chunk; ++k) {
                ev.kinetic_half(st, dt);
                std::vector<double> rho, V;
                st.density(rho);
                poisson.solve(rho, V);
                ev.potential_full(st, V, dt);
                ev.kinetic_half(st, dt);
            }
            ColumnMoments m = compute_moments_2d(ev, st, poisson);
            ev.spectral_spill(st, 0.5 * g.kx_nyquist(), 0.5 * g.kz_nyquist(),
                              spx, spz);
            spill_max = std::max(spill_max, spx);
            ts.push_back((p + 1) * T / samples_per_period);
            sx.push_back(m.sigma_x_rms);
            sz.push_back(m.sigma_z_rms);
            dE = std::fabs((m.energy - m0.energy) / m0.energy);
        }

        // Least-squares slope over the whole (short) run, in km/s per vertical
        // period. Kept short deliberately so the heating stays in its linear
        // phase and the velocity headroom is not eaten by the heating itself.
        auto slope = [&](const std::vector<double>& y) {
            const int n = (int)y.size();
            double a = 0, b = 0, c = 0, d2 = 0;
            for (int i = 0; i < n; ++i) {
                a += ts[i]; b += y[i]; c += ts[i] * y[i]; d2 += ts[i] * ts[i];
            }
            return (n * c - a * b) / (n * d2 - a * a);
        };
        const double dsx = slope(sx) * T;
        const double dsz = slope(sz) * T;
        spx = spill_max;
        if (rate_ref == 0.0) rate_ref = dsx * 2.0;  // normalise on the n=2 row
        std::printf("%8d %12.4f %12.4f %12.4f %12.4f %10.2e %10.2e\n",
                    n_streams, contrast, dsx, dsz, dsx * n_streams / rate_ref,
                    dE, spx);
    }

    std::printf("\n  'rms drho/rho' should fall as 1/sqrt(N_streams) and\n"
                "  'rate*N' should be flat if the heating goes as 1/N_streams.\n"
                "  A flat last column means 2D is safe (streams are cheap) and\n"
                "  that 3D, where a single wavefunction is the only affordable\n"
                "  option, sits at the worst point of this scaling.\n");
    return 0;
}
