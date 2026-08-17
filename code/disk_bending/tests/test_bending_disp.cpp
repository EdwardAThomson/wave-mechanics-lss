// Stage 1, test 5 proper: the stable bending dispersion relation, measured in
// the window an unsheared box can actually reach.
//
// The razor-thin relations are degenerate: bending is stable exactly where the
// sheet fragments in plane. But that degeneracy is a razor-thin statement, and
// the firehose run showed it does not survive finite thickness. At
// sigma_x = 45 km/s the sheet is in-plane stable (k > k_J) AND the bending
// modes oscillate rather than grow, where razor-thin theory says they should be
// firehose-unstable. So there is a window after all, at k h of order one, and
// this measures the dispersion relation inside it.
//
// The measurement is deliberately assumption-free. Rather than compare against
// one guessed finite-thickness formula, it scans sigma_x at fixed k and fits
//
//     omega^2 = A_eff - B_eff sigma_x^2
//
// Razor-thin theory predicts A = 2 pi G Sigma k (the self-gravity restoring
// term) and B = k^2 (the destabilising in-plane pressure term). Fitting both
// separately measures how finite thickness renormalises each, which is exactly
// the finite-thickness correction §7 asks to have verified.
//
// It also tests README correction 8 directly. §7 writes the destabilising term
// as k^2 sigma_z^2; the claim there is that it should be the IN-PLANE
// dispersion along the wavevector. So the run scans sigma_z as well: if
// omega^2 tracks sigma_x^2 and is insensitive to sigma_z, correction 8 is
// confirmed by measurement rather than by citation.
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

struct ModeFit {
    double omega = 0.0;   // dominant frequency
    double gamma = 0.0;   // envelope log-slope (positive = growing)
    int n_flip = 0;
};

// Dominant frequency from the power spectrum of a real time series, with
// parabolic interpolation of the peak. More robust than zero crossings when
// more than one frequency is present.
ModeFit fit_mode(const std::vector<double>& t, const std::vector<double>& y) {
    ModeFit f;
    const int n = static_cast<int>(y.size());
    if (n < 8) return f;
    const double dt = (t.back() - t.front()) / (n - 1);

    double mean = 0.0;
    for (double v : y) mean += v;
    mean /= n;

    // Brute-force periodogram over a generous frequency band; n is small.
    const double w_max = units::PI / dt;
    const int nw = 4000;
    int best = 1;
    double best_p = -1.0;
    std::vector<double> P(nw + 1, 0.0);
    for (int i = 1; i <= nw; ++i) {
        const double w = w_max * i / nw;
        double c = 0.0, s = 0.0;
        for (int j = 0; j < n; ++j) {
            c += (y[j] - mean) * std::cos(w * t[j]);
            s += (y[j] - mean) * std::sin(w * t[j]);
        }
        P[i] = c * c + s * s;
        if (P[i] > best_p) { best_p = P[i]; best = i; }
    }
    double shift = 0.0;
    if (best > 1 && best < nw) {
        const double a = P[best - 1], b = P[best], c = P[best + 1];
        const double den = a - 2.0 * b + c;
        if (std::fabs(den) > 0.0) shift = 0.5 * (a - c) / den;
    }
    f.omega = w_max * (best + shift) / nw;

    for (int j = 1; j < n; ++j) {
        if ((y[j - 1] - mean) * (y[j] - mean) < 0.0) ++f.n_flip;
    }

    // Envelope slope from |y| peaks, crude but enough to flag growth.
    double sa = 0, sb = 0, sc = 0, sd = 0;
    int used = 0;
    for (int j = 1; j < n - 1; ++j) {
        const double a0 = std::fabs(y[j] - mean);
        if (a0 < std::fabs(y[j - 1] - mean) || a0 < std::fabs(y[j + 1] - mean))
            continue;
        if (a0 <= 0.0) continue;
        const double l = std::log(a0);
        sa += t[j]; sb += l; sc += t[j] * l; sd += t[j] * t[j];
        ++used;
    }
    if (used > 2) f.gamma = (used * sc - sa * sb) / (used * sd - sa * sa);
    return f;
}

struct Result { ModeFit m[2]; };

Result run_case(double sigma_x, double sigma_z, double hbar, double Sigma_in,
                const Grid2D& g, int n_streams, double h0, double* h_out) {
    Grid1D g1(g.Nz, g.Lz);
    SheetConfig cfg;
    cfg.Sigma = Sigma_in;
    cfg.sigma_z = sigma_z;
    cfg.hbar = hbar;
    cfg.nu_ext = 0.0;
    cfg.n_states = 300;
    cfg.verbose = false;
    WarmSheet sheet = build_warm_sheet(g1, cfg);
    *h_out = sheet.h_measured;

    Evolver2D ev(g, hbar, Kinetic::Spectral);
    PoissonXZ poisson(g);
    SlabState2D st = make_warm_2d(sheet, g, sigma_x, n_streams, 8080ULL);

    std::vector<double> shift(g.Nx, 0.0);
    for (int i = 0; i < g.Nx; ++i) {
        for (int n = 1; n <= 2; ++n) shift[i] += h0 * std::cos(g.kx(n) * g.x(i));
    }
    ev.displace_columns(st, shift);

    // Run length from the SELF-GRAVITY frequency sqrt(2 pi G Sigma k), not
    // from k sigma_x. The first attempt used the latter and overestimated the
    // frequency about fivefold, so the run was shorter than a single period
    // and the periodogram fitted a trend instead of an oscillation.
    const double w_est = std::sqrt(2.0 * units::PI * units::G * Sigma_in * g.kx(1));
    const double t_end = 5.0 * 2.0 * units::PI / w_est;
    const double dt = (2.0 * units::PI / sheet.omega0) / 250.0;
    const int nsteps = static_cast<int>(t_end / dt);
    const int sample_every = std::max(1, nsteps / 200);

    std::vector<double> ts;
    std::vector<std::vector<double>> ser(2);
    std::vector<double> rho, V;
    for (int step = 0; step <= nsteps; ++step) {
        if (step % sample_every == 0) {
            ColumnMoments m = compute_columns_2d(ev, st);
            ts.push_back(step * dt);
            for (int n = 1; n <= 2; ++n) {
                ser[n - 1].push_back(mode_amplitude(g, m.zbar, n).real());
            }
        }
        if (step == nsteps) break;
        ev.kinetic_half(st, dt);
        st.density(rho);
        poisson.solve(rho, V);
        ev.potential_full(st, V, dt);
        ev.kinetic_half(st, dt);
    }

    Result r;
    for (int n = 0; n < 2; ++n) r.m[n] = fit_mode(ts, ser[n]);
    return r;
}

}  // namespace

int main() {
    const double Sigma = 50.0 * units::MSUN_PC2, hbar = 0.6;
    const int Nx = 768, Nz = 512, n_streams = 4;
    // Smaller box: raises k (so the oscillation is fast enough to sample many
    // periods) and LOWERS the in-plane stability threshold sqrt(G Sigma Lx).
    // h0 raised so the signal sits well above the mode-noise floor of a few
    // times 1e-3 kpc measured in test_mode_noise.
    const double Lx = 3.0, Lz = 5.0, h0 = 0.04;
    Grid2D g(Nx, Nz, Lx, Lz);

    const double k1 = g.kx(1), k2 = g.kx(2);
    std::printf("=== Stage 1, test 5: stable bending dispersion relation ===\n");
    std::printf("grid %d x %d, Lx = %.1f, Lz = %.1f kpc, %d streams\n", Nx, Nz,
                Lx, Lz, n_streams);
    std::printf("modes 1 and 2: k = %.3f and %.3f /kpc\n", k1, k2);
    std::printf("in-plane stability needs sigma_x > sqrt(G Sigma Lx) = "
                "%.1f km/s; all cases below clear it\n\n",
                std::sqrt(units::G * Sigma * Lx));

    // ---- scan sigma_x at fixed sigma_z ----
    std::printf("(a) scan sigma_x, sigma_z = 20 km/s\n");
    std::printf("%9s %7s %10s %10s %12s %8s %12s\n", "sigma_x", "k h",
                "omega_m1", "omega_m2", "w^2 (mode1)", "flips", "gamma");
    const double sxs[3] = {30.0, 40.0, 50.0};
    double w2[3];
    double h_meas = 0.0;
    for (int i = 0; i < 3; ++i) {
        Result r = run_case(sxs[i], 20.0, hbar, Sigma, g, n_streams, h0,
                            &h_meas);
        w2[i] = r.m[0].omega * r.m[0].omega;
        std::printf("%9.1f %7.3f %10.2f %10.2f %12.1f %8d %12.2f\n", sxs[i],
                    k1 * h_meas, r.m[0].omega, r.m[1].omega, w2[i],
                    r.m[0].n_flip, r.m[0].gamma);
    }

    // Fit omega^2 = A - B sigma_x^2 through the three points.
    {
        double sa = 0, sb = 0, sc = 0, sd = 0;
        for (int i = 0; i < 3; ++i) {
            const double x = sxs[i] * sxs[i];
            sa += x; sb += w2[i]; sc += x * w2[i]; sd += x * x;
        }
        const double B = -(3 * sc - sa * sb) / (3 * sd - sa * sa);
        const double A = (sb + B * sa) / 3.0;
        std::printf("\n  fit omega^2 = A - B sigma_x^2 at k = %.3f:\n", k1);
        std::printf("    A_eff = %10.1f   razor-thin 2 pi G Sigma k = %10.1f"
                    "   ratio %.3f\n",
                    A, 2.0 * units::PI * units::G * Sigma * k1,
                    A / (2.0 * units::PI * units::G * Sigma * k1));
        std::printf("    B_eff = %10.4f   razor-thin k^2            = %10.4f"
                    "   ratio %.3f\n",
                    B, k1 * k1, B / (k1 * k1));
        std::printf("  Ratios below one are the finite-thickness reduction of\n"
                    "  each term at k h = %.2f.\n", k1 * h_meas);
    }

    // ---- scan sigma_z at fixed sigma_x: correction 8 ----
    std::printf("\n(b) scan sigma_z at fixed sigma_x = 40 km/s "
                "(README correction 8)\n");
    std::printf("%9s %9s %7s %10s %12s\n", "sigma_z", "h[kpc]", "k h",
                "omega_m1", "w^2 (mode1)");
    for (double sz : {14.0, 20.0, 26.0}) {
        Result r = run_case(40.0, sz, hbar, Sigma, g, n_streams, h0, &h_meas);
        std::printf("%9.1f %9.4f %7.3f %10.2f %12.1f\n", sz, h_meas,
                    k1 * h_meas, r.m[0].omega, r.m[0].omega * r.m[0].omega);
    }
    std::printf("\n  §7 writes the destabilising term as k^2 sigma_z^2. If that\n"
                "  were right, omega^2 would fall steeply with sigma_z here.\n"
                "  Correction 8 says the term is the in-plane dispersion, so\n"
                "  sigma_z should enter only weakly, through the scale height\n"
                "  and hence the finite-thickness reduction.\n");
    return 0;
}
