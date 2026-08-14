// Stage 1: is the in-plane heating granule self-gravity, or is it numerical?
//
// The Stage 1 first pass measured a large heating rate that refused to fall
// with stream count, which contradicts the obvious argument (density contrast
// falls as 1/sqrt(N), driving forces fall likewise, heating should fall as
// 1/N). Two candidate explanations: the argument is wrong, or the runs were
// contaminated by velocity-space aliasing, since the heating eats its own
// headroom.
//
// This program separates them.
//
// PART 1, the decisive control. Evolve one identical initial condition three
// ways:
//
//   self-gravity   V = Poisson(rho) recomputed every step. Granules pull on
//                  each other. This is the physical case.
//   frozen smooth  V = Poisson(rho_smooth) computed once from the horizontally
//                  averaged equilibrium profile, then held fixed. The vertical
//                  structure is identical, but granule self-gravity is gone.
//   frozen initial V = Poisson(rho) at t = 0, granules included, held fixed.
//                  Static, so no stochastic forcing, but the granule field is
//                  still there.
//
// Under "frozen smooth" the basis functions u_n(z) exp(i kx x) are EXACT
// eigenstates of the evolution operator, so sigma_x and sigma_z must be
// constant to solver precision. Any drift there is numerical, full stop. If
// heating appears only in the self-gravity run, it is granule relaxation and
// it is real for this model.
//
// PART 2, the stream scan redone with the velocity headroom raised from 12
// sigma to 24 sigma, so that even a doubling of sigma_x leaves the Maxwellian
// tail resolved and the spill column stays negligible throughout.
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "../src/diagnostics2d.h"
#include "../src/equilibrium.h"
#include "../src/evolve2d.h"
#include "../src/grid2d.h"
#include "../src/poisson_xz.h"
#include "../src/slab_ic_2d.h"
#include "../src/units.h"

namespace {

enum class Gravity { SelfConsistent, FrozenSmooth, FrozenInitial };

struct Trace {
    std::vector<double> t, sx, sz, spill;
    double dE = 0.0;
};

Trace run(const WarmSheet& sheet, const Grid2D& g, Evolver2D& ev,
          PoissonXZ& poisson, SlabState2D st, Gravity mode, double sigma_x,
          int n_per, int spp, int samples_per_period) {
    (void)sigma_x;
    const double T = 2.0 * units::PI / sheet.omega0;
    const double dt = T / spp;

    // Frozen potentials, built once.
    std::vector<double> V_frozen;
    if (mode == Gravity::FrozenSmooth) {
        std::vector<double> rho_smooth(g.size());
        for (int i = 0; i < g.Nx; ++i) {
            for (int j = 0; j < g.Nz; ++j) rho_smooth[g.idx(i, j)] = sheet.rho[j];
        }
        poisson.solve(rho_smooth, V_frozen);
    } else if (mode == Gravity::FrozenInitial) {
        std::vector<double> rho0;
        st.density(rho0);
        poisson.solve(rho0, V_frozen);
    }

    Trace tr;
    ColumnMoments m0 = compute_moments_2d(ev, st, poisson);
    const int chunk = spp / samples_per_period;
    std::vector<double> rho, V;
    for (int p = 0; p < n_per * samples_per_period; ++p) {
        for (int k = 0; k < chunk; ++k) {
            ev.kinetic_half(st, dt);
            if (mode == Gravity::SelfConsistent) {
                st.density(rho);
                poisson.solve(rho, V);
                ev.potential_full(st, V, dt);
            } else {
                ev.potential_full(st, V_frozen, dt);
            }
            ev.kinetic_half(st, dt);
        }
        ColumnMoments m = compute_moments_2d(ev, st, poisson);
        double spx = 0.0, spz = 0.0;
        ev.spectral_spill(st, 0.5 * g.kx_nyquist(), 0.5 * g.kz_nyquist(), spx,
                          spz);
        tr.t.push_back((p + 1) * T / samples_per_period);
        tr.sx.push_back(m.sigma_x_rms);
        tr.sz.push_back(m.sigma_z_rms);
        tr.spill.push_back(spx);
        tr.dE = std::fabs((m.energy - m0.energy) / m0.energy);
    }
    tr.t.insert(tr.t.begin(), 0.0);
    tr.sx.insert(tr.sx.begin(), m0.sigma_x_rms);
    tr.sz.insert(tr.sz.begin(), m0.sigma_z_rms);
    tr.spill.insert(tr.spill.begin(), 0.0);
    return tr;
}

double slope_per_period(const Trace& tr, const std::vector<double>& y,
                        double T) {
    const int n = static_cast<int>(y.size());
    double a = 0, b = 0, c = 0, d = 0;
    for (int i = 0; i < n; ++i) {
        a += tr.t[i]; b += y[i]; c += tr.t[i] * y[i]; d += tr.t[i] * tr.t[i];
    }
    return (n * c - a * b) / (n * d - a * a) * T;
}

}  // namespace

int main() {
    const double Sigma_pc2 = 50.0, sigma_z = 20.0, sigma_x = 10.0, hbar = 0.6;
    // 24 sigma_x and 14 sigma_z of headroom at t = 0, so that even a doubling
    // of the dispersion leaves the tail resolved.
    const int Nx = 512, Nz = 768;
    const double Lx = 4.0, Lz = 5.0;
    const int n_per = 2, spp = 200, sper = 4;

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
    const double lam = 2.0 * units::PI * hbar / sigma_x;

    std::printf("=== Stage 1: is the heating granule self-gravity? ===\n");
    std::printf("grid %d x %d, Lx = %.1f, Lz = %.1f kpc\n", Nx, Nz, Lx, Lz);
    std::printf("sigma_z = %.1f, sigma_x = %.1f km/s, hbar_eff = %.2f, "
                "h = %.4f kpc\n", sigma_z, sigma_x, hbar, sheet.h_measured);
    std::printf("in-plane lambda_dB = %.4f kpc, so Lx holds %.1f granules\n",
                lam, Lx / lam);
    std::printf("v ceiling: x %.1f (%.1f sigma_x), z %.1f (%.1f sigma_z)\n",
                ev.v_ceiling_x(), ev.v_ceiling_x() / sigma_x, ev.v_ceiling_z(),
                ev.v_ceiling_z() / sigma_z);
    std::printf("evolving %d vertical periods (%.0f Myr)\n\n", n_per,
                n_per * T * units::TIME_MYR);

    // ---------------- PART 1: the control ----------------
    std::printf("PART 1: same initial condition, three gravity treatments "
                "(8 streams)\n");
    std::printf("%-16s %10s %10s %12s %12s %10s %10s\n", "gravity", "sx_start",
                "sx_end", "dsx/period", "dsz/period", "dE/E", "max spill");

    const int ns_ctrl = 8;
    struct Case { Gravity mode; const char* name; };
    const Case cases[3] = {{Gravity::SelfConsistent, "self-gravity"},
                           {Gravity::FrozenSmooth, "frozen smooth"},
                           {Gravity::FrozenInitial, "frozen initial"}};
    for (const auto& c : cases) {
        SlabState2D st = make_warm_2d(sheet, g, sigma_x, ns_ctrl, 909001ULL);
        Trace tr = run(sheet, g, ev, poisson, st, c.mode, sigma_x, n_per, spp,
                       sper);
        double smax = 0.0;
        for (double s : tr.spill) smax = std::max(smax, s);
        std::printf("%-16s %10.4f %10.4f %12.4f %12.4f %10.2e %10.2e\n", c.name,
                    tr.sx.front(), tr.sx.back(),
                    slope_per_period(tr, tr.sx, T),
                    slope_per_period(tr, tr.sz, T), tr.dE, smax);
    }
    std::printf("\n  Under 'frozen smooth' the basis functions are exact\n"
                "  eigenstates of the evolution operator, so sigma_x and\n"
                "  sigma_z MUST be flat. Any drift in that row is numerical and\n"
                "  sets the noise floor for the other two.\n\n");

    // ---------------- PART 2: stream scan at high headroom ----------------
    std::printf("PART 2: stream scan, self-gravity, 24 sigma_x of headroom\n");
    std::printf("%8s %12s %12s %12s %12s %10s %10s\n", "streams",
                "rms drho/rho", "dsx/period", "dsz/period", "rate x N", "dE/E",
                "max spill");
    double ref = 0.0;
    for (int ns : {2, 8, 32}) {
        SlabState2D st = make_warm_2d(sheet, g, sigma_x, ns, 515001ULL + 13 * ns);

        double contrast = 0.0;
        {
            std::vector<double> rho;
            st.density(rho);
            double mean = 0.0; int cnt = 0;
            for (int i = 0; i < Nx; ++i)
                for (int j = 0; j < Nz; ++j) {
                    if (std::fabs(g.z(j)) > 0.5 * sheet.h_measured) continue;
                    mean += rho[g.idx(i, j)]; ++cnt;
                }
            mean /= cnt;
            double s2 = 0.0;
            for (int i = 0; i < Nx; ++i)
                for (int j = 0; j < Nz; ++j) {
                    if (std::fabs(g.z(j)) > 0.5 * sheet.h_measured) continue;
                    const double d = rho[g.idx(i, j)] / mean - 1.0;
                    s2 += d * d;
                }
            contrast = std::sqrt(s2 / cnt);
        }

        Trace tr = run(sheet, g, ev, poisson, st, Gravity::SelfConsistent,
                       sigma_x, n_per, spp, sper);
        double smax = 0.0;
        for (double s : tr.spill) smax = std::max(smax, s);
        const double dsx = slope_per_period(tr, tr.sx, T);
        if (ref == 0.0) ref = dsx * 2.0;
        std::printf("%8d %12.4f %12.4f %12.4f %12.4f %10.2e %10.2e\n", ns,
                    contrast, dsx, slope_per_period(tr, tr.sz, T),
                    dsx * ns / ref, tr.dE, smax);
    }
    std::printf("\n  'rate x N' flat means heating falls as 1/N_streams.\n"
                "  Rising means it falls more slowly than that, and the\n"
                "  1/sqrt(N) contrast argument does not carry over to the\n"
                "  heating rate.\n");
    return 0;
}
