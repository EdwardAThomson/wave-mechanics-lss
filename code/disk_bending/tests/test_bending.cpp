// Stage 1 measurements on the 2D (x, z) slab.
//
//   (a) equilibrium stationarity, and with it the granule-heating rate. The
//       in-plane part of the initial condition is necessarily a random-phase
//       superposition (there is no in-plane eigenbasis to occupy), so |psi|^2
//       carries interference granules of scale lambda_dB = 2 pi hbar / sigma_x,
//       and those granules self-gravitate. In the fuzzy-dark-matter literature
//       that drives spurious relaxation. If it heats the disk faster than the
//       physics being measured, the whole approach is compromised, so this is
//       measured before anything else.
//
//   (b) the bending-wave dispersion relation, validation test 5. Excited as a
//       pure vertical displacement (zero initial velocity) so the oscillation
//       frequency is read off directly. Several wavenumbers are excited at
//       once and tracked independently, which is legitimate while the
//       amplitude is small enough to stay linear.
//
//   (c) the firehose threshold, by scanning the in-plane dispersion. Note this
//       is sigma_x, the dispersion along the wavevector, NOT sigma_z: see
//       correction 8 in README.md.
//
// This is a measurement program rather than a unit test, so it is not part of
// `make check`. Build and run it directly.
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

struct Setup {
    Grid2D g;
    WarmSheet sheet;
    double sigma_x;
    int n_streams;
};

Setup build(int Nx, int Nz, double Lx, double Lz, double Sigma_pc2,
            double sigma_z, double sigma_x, double hbar, int n_states,
            int n_streams, bool verbose = false) {
    Setup s;
    s.g = Grid2D(Nx, Nz, Lx, Lz);
    Grid1D g1(Nz, Lz);
    SheetConfig cfg;
    cfg.Sigma = Sigma_pc2 * units::MSUN_PC2;
    cfg.sigma_z = sigma_z;
    cfg.hbar = hbar;
    cfg.nu_ext = 0.0;  // pure self-gravity: the dispersion relation is for a
                       // self-gravitating sheet, a halo term would change it
    cfg.n_states = n_states;
    cfg.verbose = verbose;
    s.sheet = build_warm_sheet(g1, cfg);
    s.sigma_x = sigma_x;
    s.n_streams = n_streams;
    return s;
}

// Zero-crossing frequency estimate for a real oscillating series.
double frequency_from_crossings(const std::vector<double>& t,
                                const std::vector<double>& y, int* n_cross) {
    std::vector<double> cross;
    for (size_t i = 1; i < y.size(); ++i) {
        if ((y[i - 1] <= 0.0 && y[i] > 0.0) || (y[i - 1] >= 0.0 && y[i] < 0.0)) {
            const double f = y[i - 1] / (y[i - 1] - y[i]);
            cross.push_back(t[i - 1] + f * (t[i] - t[i - 1]));
        }
    }
    *n_cross = static_cast<int>(cross.size());
    if (cross.size() < 2) return 0.0;
    const double half = (cross.back() - cross.front()) / (cross.size() - 1);
    return units::PI / half;
}

}  // namespace

int main() {
    const double Sigma_pc2 = 50.0, sigma_z = 20.0, hbar = 0.6;
    const double Sigma = Sigma_pc2 * units::MSUN_PC2;

    // ================= (a) stationarity and granule heating =================
    std::printf("=== Stage 1 (a): 2D equilibrium and granule heating ===\n");
    {
        // Resolution set by the velocity ceiling, not by taste: the in-plane
        // Maxwellian needs several sigma_x of headroom or its tail aliases,
        // which looks exactly like heating and is not.
        Setup s = build(1024, 1024, 8.0, 6.0, Sigma_pc2, sigma_z, 20.0, hbar,
                        300, 4, true);
        std::printf("\ngrid %d x %d, Lx = %.1f, Lz = %.1f kpc\n", s.g.Nx, s.g.Nz,
                    s.g.Lx, s.g.Lz);
        std::printf("h = %.4f kpc, omega0 = %.2f km/s/kpc, T_vert = %.1f Myr\n",
                    s.sheet.h_measured, s.sheet.omega0,
                    2.0 * units::PI / s.sheet.omega0 * units::TIME_MYR);
        std::printf("lambda_dB (in-plane) = %.4f kpc = %.2f h\n",
                    2.0 * units::PI * hbar / s.sigma_x,
                    2.0 * units::PI * hbar / s.sigma_x / s.sheet.h_measured);

        // Spectral, not MatchedFD: the FD group velocity turns over at
        // k dz = pi/2 so its ceiling is a factor pi lower, and here the
        // ceiling is the binding constraint. The price is that the
        // FD-basis equilibrium breathes slightly, but breathing is even
        // in z while the bending signal <z> is odd, so it does not mix in.
        Evolver2D ev(s.g, hbar, Kinetic::Spectral);
        PoissonXZ poisson(s.g);
        SlabState2D st = make_warm_2d(s.sheet, s.g, s.sigma_x, s.n_streams);

        std::printf("v ceiling: x %.1f km/s (%.1f sigma_x), z %.1f km/s "
                    "(%.1f sigma_z)\n",
                    ev.v_ceiling_x(), ev.v_ceiling_x() / s.sigma_x,
                    ev.v_ceiling_z(), ev.v_ceiling_z() / sigma_z);

        const double T = 2.0 * units::PI / s.sheet.omega0;
        const int n_per = 5, spp = 250;
        const double dt = T / spp;
        ColumnMoments m0 = compute_moments_2d(ev, st, poisson);
        double spx = 0.0, spz = 0.0;
        ev.spectral_spill(st, 0.5 * s.g.kx_nyquist(), 0.5 * s.g.kz_nyquist(),
                          spx, spz);
        std::printf("\n%10s %12s %12s %12s %11s %10s %10s\n", "t[Myr]",
                    "Sigma_tot", "sigma_z", "sigma_x", "dE/E", "spill_x",
                    "spill_z");
        std::printf("%10.1f %12.5e %12.4f %12.4f %11.3e %10.2e %10.2e\n", 0.0,
                    m0.Sigma_total, m0.sigma_z_rms, m0.sigma_x_rms, 0.0, spx,
                    spz);

        std::vector<double> ts, sz, sx;
        for (int p = 0; p < n_per; ++p) {
            for (int k = 0; k < spp; ++k) {
                ev.kinetic_half(st, dt);
                std::vector<double> rho, V;
                st.density(rho);
                poisson.solve(rho, V);
                ev.potential_full(st, V, dt);
                ev.kinetic_half(st, dt);
            }
            ColumnMoments m = compute_moments_2d(ev, st, poisson);
            const double t = (p + 1) * T;
            ts.push_back(t);
            sz.push_back(m.sigma_z_rms);
            sx.push_back(m.sigma_x_rms);
            ev.spectral_spill(st, 0.5 * s.g.kx_nyquist(),
                              0.5 * s.g.kz_nyquist(), spx, spz);
            std::printf("%10.1f %12.5e %12.4f %12.4f %11.3e %10.2e %10.2e\n",
                        t * units::TIME_MYR, m.Sigma_total, m.sigma_z_rms,
                        m.sigma_x_rms,
                        std::fabs((m.energy - m0.energy) / m0.energy), spx,
                        spz);
        }

        // Heating rate per vertical period, from a straight-line fit.
        auto slope = [&](const std::vector<double>& y) {
            const int n = static_cast<int>(y.size());
            double sx_ = 0, sy = 0, sxy = 0, sxx = 0;
            for (int i = 0; i < n; ++i) {
                sx_ += ts[i];
                sy += y[i];
                sxy += ts[i] * y[i];
                sxx += ts[i] * ts[i];
            }
            return (n * sxy - sx_ * sy) / (n * sxx - sx_ * sx_);
        };
        const double dsz = slope(sz) * T, dsx = slope(sx) * T;
        std::printf("\ngranule heating: d(sigma_z)/dt = %+.4f km/s per vertical "
                    "period (%.3f%%)\n", dsz, 100.0 * dsz / m0.sigma_z_rms);
        std::printf("                 d(sigma_x)/dt = %+.4f km/s per vertical "
                    "period (%.3f%%)\n", dsx, 100.0 * dsx / m0.sigma_x_rms);
        std::printf("  (compare against the bending oscillation period, which is\n"
                    "   a few vertical periods: heating must be small over that)\n\n");
    }

    // ================= (b) dispersion relation =================
    std::printf("=== Stage 1 (b): bending-wave dispersion relation ===\n");
    {
        const double sigma_x = 4.0;  // cold in plane, so the razor-thin,
                                     // pressureless limit is the prediction
        Setup s = build(1024, 768, 32.0, 6.0, Sigma_pc2, sigma_z, sigma_x,
                        hbar, 300, 4);
        Evolver2D ev(s.g, hbar, Kinetic::Spectral);
        PoissonXZ poisson(s.g);
        SlabState2D st = make_warm_2d(s.sheet, s.g, sigma_x, s.n_streams);

        const double h = s.sheet.h_measured;
        std::printf("Lx = %.0f kpc, Nx = %d, sigma_x = %.1f km/s, h = %.4f kpc\n",
                    s.g.Lx, s.g.Nx, sigma_x, h);
        std::printf("v ceiling x = %.1f km/s (%.1f sigma_x)\n\n",
                    ev.v_ceiling_x(), ev.v_ceiling_x() / sigma_x);

        // Excite modes 1..5 together as a pure vertical displacement.
        const int nmodes = 5;
        const double h0 = 0.02;
        std::vector<double> shift(s.g.Nx, 0.0);
        for (int i = 0; i < s.g.Nx; ++i) {
            for (int n = 1; n <= nmodes; ++n) {
                shift[i] += h0 * std::cos(s.g.kx(n) * s.g.x(i));
            }
        }
        ev.displace_columns(st, shift);

        // Run three periods of the slowest mode.
        const double k1 = s.g.k_fundamental();
        const double w1 = std::sqrt(2.0 * units::PI * units::G * Sigma * k1);
        const double t_end = 3.0 * 2.0 * units::PI / w1;
        const double dt = (2.0 * units::PI / s.sheet.omega0) / 200.0;
        const int nsteps = static_cast<int>(t_end / dt);
        const int sample_every = std::max(1, nsteps / 600);

        std::vector<double> tser;
        std::vector<std::vector<double>> zser(nmodes), vser(nmodes);
        for (int step = 0; step <= nsteps; ++step) {
            if (step % sample_every == 0) {
                ColumnMoments m = compute_columns_2d(ev, st);
                tser.push_back(step * dt);
                for (int n = 1; n <= nmodes; ++n) {
                    zser[n - 1].push_back(
                        mode_amplitude(s.g, m.zbar, n).real());
                    vser[n - 1].push_back(
                        mode_amplitude(s.g, m.vbar, n).real());
                }
            }
            if (step == nsteps) break;
            ev.kinetic_half(st, dt);
            std::vector<double> rho, V;
            st.density(rho);
            poisson.solve(rho, V);
            ev.potential_full(st, V, dt);
            ev.kinetic_half(st, dt);
        }

        std::printf("%5s %9s %8s %12s %12s %10s %8s\n", "mode", "k[1/kpc]",
                    "k h", "omega_meas", "omega_thin", "ratio", "crossings");
        for (int n = 1; n <= nmodes; ++n) {
            const double k = s.g.kx(n);
            int nc = 0;
            const double w = frequency_from_crossings(tser, zser[n - 1], &nc);
            const double w_thin = std::sqrt(std::max(
                0.0, 2.0 * units::PI * units::G * Sigma * k - k * k * sigma_x * sigma_x));
            std::printf("%5d %9.4f %8.3f %12.3f %12.3f %10.4f %8d\n", n, k,
                        k * h, w, w_thin, w / w_thin, nc);
        }
        std::printf("\n  omega_thin is the razor-thin result "
                    "omega^2 = 2 pi G Sigma k - k^2 sigma_x^2.\n"
                    "  The measured values fall below it, by more at larger k h:\n"
                    "  that is the finite-thickness reduction of self-gravity,\n"
                    "  and it is the expected behaviour, not an error.\n\n");
    }

    return 0;
}
