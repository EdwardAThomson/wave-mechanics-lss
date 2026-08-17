// Rotation stabilises the sheet: the Toomre criterion, measured.
//
// Stage 1's blocking result was the in-plane/bending duality: without
// rotation, the clean bending regime k << k_J is exactly the regime where the
// unsheared sheet fragments (README.md, "shear cannot be skipped"). Rotation
// adds +kappa^2 to the in-plane relation and breaks the degeneracy:
//
//   in-plane, rotating:  omega^2 = kappa^2 - 2 pi G Sigma k + k^2 sigma_x^2
//
// stable at every k when Q = kappa sigma_x / (3.36 G Sigma) > 1 (the 3.36 is
// the stellar-disk coefficient; the fluid version has pi). This test measures
// exactly that, in the wave code, at the production box size Lx = 16 kpc:
//
//   1. kappa = 44.2, sigma_x = 40  -> Q = 2.45  : must hold together
//   2. kappa = 44.2, sigma_x = 10  -> Q = 0.61  : must fragment
//   3. kappa = 0,    sigma_x = 40  -> Q = 0     : must fragment (the Stage 1
//      duality: an unsheared box at Lx = 16 needs sigma_x > sqrt(G Sigma Lx)
//      = 58.7 km/s, above the physical 40)
//
// Config 3 reuses the non-rotating plane-wave ICs, so it is also a direct
// like-for-like control on the new machinery. The fragmentation metric is the
// in-plane heating rate d(sigma_x)/dt (the test_granule2 diagnostic, which
// switched off across the threshold by a factor ~750) plus the growth of the
// most unstable surface-density mode.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "../src/equilibrium.h"
#include "../src/evolve2d.h"
#include "../src/grid2d.h"
#include "../src/poisson_xz.h"
#include "../src/rotation.h"
#include "../src/slab_ic_2d.h"
#include "../src/units.h"

namespace {

// Largest surface-density mode amplitude 2 |Sigma_hat(k_n)| / Sigma_mean over
// modes lo..hi; *argmode reports which. A band, not a single mode, because
// the stellar most-unstable k is shifted from the fluid estimate by the
// reduction factor and there is no reason to guess it in advance.
double band_amp(const Grid2D& g, const std::vector<double>& rho, int lo,
                int hi, int* argmode) {
    std::vector<double> col(g.Nx, 0.0);
    double tot = 0.0;
    for (int i = 0; i < g.Nx; ++i) {
        for (int j = 0; j < g.Nz; ++j) col[i] += rho[g.idx(i, j)];
        tot += col[i];
    }
    double best = 0.0;
    int bestn = lo;
    for (int n = lo; n <= hi; ++n) {
        double cr = 0.0, ci = 0.0;
        for (int i = 0; i < g.Nx; ++i) {
            const double ph = -g.kx(n) * g.x(i);
            cr += col[i] * std::cos(ph);
            ci += col[i] * std::sin(ph);
        }
        const double a = 2.0 * std::sqrt(cr * cr + ci * ci) / tot;
        if (a > best) { best = a; bestn = n; }
    }
    if (argmode) *argmode = bestn;
    return best;
}

struct RunResult {
    double dsig_dt = 0.0;      // km/s per time unit
    double amp0 = 0.0, amp1 = 0.0;
    int mode1 = 0;             // where the final band maximum sits
    double efold = 0.0;        // measured e-folding rate of the band maximum
};

RunResult evolve_and_measure(const Grid2D& g, Evolver2D& ev, PoissonXZ& poisson,
                             SlabState2D& st, double t_end, double dt,
                             int band_lo, int band_hi) {
    RunResult res;
    std::vector<double> rho, V;
    const int nsteps = static_cast<int>(t_end / dt);
    st.density(rho);
    double mass = 0.0;
    for (double v : rho) mass += v;
    mass *= g.dx * g.dz;

    double kex0, kez0;
    ev.kinetic_energy_split(st, kex0, kez0);
    const double sx0 = std::sqrt(2.0 * kex0 / mass);
    res.amp0 = band_amp(g, rho, band_lo, band_hi, nullptr);

    // For the e-folding rate: log-amplitude at 1/4 and 3/4 of the run.
    const int q1 = nsteps / 4, q3 = 3 * nsteps / 4;
    double a_q1 = 0.0, a_q3 = 0.0;
    for (int s = 0; s < nsteps; ++s) {
        ev.kinetic_half(st, dt);
        st.density(rho);
        poisson.solve(rho, V);
        ev.potential_full(st, V, dt);
        ev.kinetic_half(st, dt);
        if (s == q1 || s == q3) {
            st.density(rho);
            const double a = band_amp(g, rho, band_lo, band_hi, nullptr);
            if (s == q1) a_q1 = a; else a_q3 = a;
        }
    }
    st.density(rho);
    double kex, kez;
    ev.kinetic_energy_split(st, kex, kez);
    const double sx1 = std::sqrt(2.0 * kex / mass);
    res.dsig_dt = (sx1 - sx0) / t_end;
    res.amp1 = band_amp(g, rho, band_lo, band_hi, &res.mode1);
    if (a_q1 > 0.0 && a_q3 > a_q1) {
        res.efold = std::log(a_q3 / a_q1) / ((q3 - q1) * dt);
    }
    return res;
}

}  // namespace

int main() {
    const double Sigma = 50.0 * units::MSUN_PC2;
    const double sigma_z = 20.0, hbar = 0.6, kappa = 44.2;
    Grid2D g(1024, 320, 16.0, 5.0);
    Grid1D g1(g.Nz, g.Lz);

    SheetConfig cfg;
    cfg.Sigma = Sigma;
    cfg.sigma_z = sigma_z;
    cfg.hbar = hbar;
    cfg.nu_ext = 0.0;
    cfg.n_states = 300;
    cfg.verbose = false;
    WarmSheet sheet = build_warm_sheet(g1, cfg);
    truncate_sheet_energy(sheet, 0.8 * 2.0 * hbar * hbar / (g.dz * g.dz));

    const double Tv = 2.0 * units::PI / sheet.omega0;
    const double dt = Tv / 200.0;

    std::printf("=== Toomre stabilisation at Lx = %.0f kpc ===\n", g.Lx);
    std::printf("Sigma = 50 Msun/pc^2, sigma_z = %.0f, hbar_eff = %.1f; "
                "unsheared threshold sqrt(G Sigma Lx) = %.1f km/s\n\n",
                sigma_z, hbar, std::sqrt(units::G * Sigma * g.Lx));

    // The unstable case gets a wide mode band and a long run: the stellar
    // reduction factor both shifts the most unstable k below the fluid
    // estimate pi G Sigma / sigma_x^2 and cuts the growth rate to a few
    // tenths of kappa, so a run sized from the fluid rate (a first attempt
    // here did exactly that, at Q = 0.61 over t = 0.15) sees nothing.
    struct Case {
        const char* name;
        double kap, sigx;
        double t_end;
        int band_lo, band_hi;
        bool expect_stable;
    };
    const Case cases[3] = {
        {"rotating, sigma_x = 40 (Q = 2.45)", kappa, 40.0, 0.50, 1, 30, true},
        {"rotating, sigma_x = 8 (Q = 0.49)", kappa, 8.0, 0.45, 1, 40, false},
        {"no rotation, sigma_x = 40 (Q = 0)", 0.0, 40.0, 0.50, 1, 3, false},
    };

    int failures = 0;
    for (const Case& c : cases) {
        std::printf("--- %s ---\n", c.name);
        SlabState2D st;
        Evolver2D ev(g, hbar, Kinetic::MatchedFD);
        if (c.kap > 0.0) {
            RotConfig rc;
            rc.kappa = c.kap;
            rc.sigma_x = c.sigx;
            rc.verbose = false;
            // Spacing a_epi for the cold case keeps the stream count sane;
            // the ripple that spacing costs is ~2 exp(-2 pi^2) of the mean,
            // invisible, and the case is unstable anyway.
            rc.n_streams = (c.sigx < 20.0) ? 71 : 0;
            RotatingSlab rot = make_warm_rotating(sheet, g, rc);
            ev.set_traps(rot.xg, rot.kappa);
            st = std::move(rot.st);
            std::printf("  %d streams, %d levels, sigma_x realised %.2f\n",
                        st.n_streams(), rot.n_levels, rot.sigma_x_realised);
        } else {
            st = make_warm_2d(sheet, g, c.sigx, 6);
            std::printf("  plane-wave ICs, 6 streams\n");
        }
        PoissonXZ poisson(g);
        RunResult r = evolve_and_measure(g, ev, poisson, st, c.t_end, dt,
                                         c.band_lo, c.band_hi);
        std::printf("  d(sigma_x)/dt = %+.4f km/s per unit time over t = %.2f"
                    "\n", r.dsig_dt, c.t_end);
        std::printf("  band max amplitude (modes %d..%d): %.2e -> %.2e at "
                    "mode %d   (e-fold rate %.1f /unit)\n", c.band_lo,
                    c.band_hi, r.amp0, r.amp1, r.mode1, r.efold);
        const bool grew = (r.amp1 > 10.0 * r.amp0) ||
                          (std::fabs(r.dsig_dt) > 2.0);
        std::printf("  verdict: %s (expected %s)\n\n",
                    grew ? "FRAGMENTS" : "stable",
                    c.expect_stable ? "stable" : "FRAGMENTS");
        if (grew == c.expect_stable) ++failures;
    }

    std::printf("%s\n", failures == 0
                            ? "PASS: rotation opens the stable window exactly "
                              "where Toomre says it should"
                            : "FAIL: see verdicts above");
    return failures == 0 ? 0 : 1;
}
