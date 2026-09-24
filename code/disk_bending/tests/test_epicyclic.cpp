// Validation for the rotating-frame (Landau) machinery in rotation.h.
//
// Part A, dynamics: a coherent state displaced from its guiding centre must
// oscillate at exactly the epicyclic frequency kappa, with no amplitude decay
// and conserved energy. The frequency error must converge at the Strang rate,
// order 2 in dt. Frequency is measured from zero crossings of <x> - x_g,
// where linear interpolation is third-order accurate (the cosine is locally
// odd about its zeros), so the splitting error dominates as it should.
//
// Part B, equilibrium: the full warm rotating stack (thermal epicyclic levels
// times the vertical eigenstate library) with self-gravity on must be
// stationary: uniform surface density, sigma_x and z_rms constant, kinetic
// energy in x equal to the trap energy (virial split of the Landau levels),
// total energy conserved, no spectral spill.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "../src/equilibrium.h"
#include "../src/evolve2d.h"
#include "../src/grid2d.h"
#include "../src/poisson_xz.h"
#include "../src/rotation.h"
#include "../src/units.h"

namespace {

double mean_x(const SlabState2D& st) {
    const Grid2D& g = st.grid;
    double m = 0.0, mx = 0.0;
    std::vector<double> rho;
    st.density(rho);
    for (int i = 0; i < g.Nx; ++i) {
        double col = 0.0;
        for (int j = 0; j < g.Nz; ++j) col += rho[g.idx(i, j)];
        m += col;
        mx += col * g.x(i);
    }
    return mx / m;
}

// kappa measured from the zero crossings of a sampled cosine.
double kappa_from_crossings(const std::vector<double>& t,
                            const std::vector<double>& y) {
    std::vector<double> tc;
    for (size_t i = 1; i < y.size(); ++i) {
        if ((y[i - 1] > 0.0) != (y[i] > 0.0)) {
            const double f = y[i - 1] / (y[i - 1] - y[i]);
            tc.push_back(t[i - 1] + f * (t[i] - t[i - 1]));
        }
    }
    if (tc.size() < 2) return 0.0;
    return units::PI * (tc.size() - 1) / (tc.back() - tc.front());
}

double run_coherent(double kappa, double dt_frac, double* energy_drift,
                    double* amp_err) {
    const double hbar = 0.6;
    Grid2D g(256, 64, 4.0, 4.0);
    const double ell = std::sqrt(hbar / kappa);
    const double xc = 0.5 * g.Lx, dx0 = 0.3, sz = 0.3;

    SlabState2D st;
    st.grid = g;
    st.hbar = hbar;
    st.psi.assign(1, std::vector<Complex>(g.size()));
    std::vector<double> fx(g.Nx), fz(g.Nz);
    double nx = 0.0, nz = 0.0;
    for (int i = 0; i < g.Nx; ++i) {
        const double u = g.x(i) - xc - dx0;
        fx[i] = std::exp(-0.25 * u * u / (0.5 * ell * ell));
        nx += fx[i] * fx[i];
    }
    for (int j = 0; j < g.Nz; ++j) {
        const double z = g.z(j);
        fz[j] = std::exp(-0.25 * z * z / (sz * sz));
        nz += fz[j] * fz[j];
    }
    nx = std::sqrt(nx * g.dx);
    nz = std::sqrt(nz * g.dz);
    for (int i = 0; i < g.Nx; ++i) {
        for (int j = 0; j < g.Nz; ++j) {
            st.psi[0][g.idx(i, j)] = Complex(fx[i] / nx * fz[j] / nz, 0.0);
        }
    }

    Evolver2D ev(g, hbar, Kinetic::Spectral);
    ev.set_traps({xc}, kappa);
    const std::vector<double> V(g.size(), 0.0);

    const double T = 2.0 * units::PI / kappa;
    const double dt = T * dt_frac;
    const int nsteps = static_cast<int>(std::lround(2.0 / dt_frac));

    const double E0 = ev.kinetic_energy(st) + ev.trap_energy(st);
    std::vector<double> ts, xs;
    double max_drift = 0.0;
    for (int s = 0; s <= nsteps; ++s) {
        ts.push_back(s * dt);
        xs.push_back(mean_x(st) - xc);
        if (s % 16 == 0) {
            const double E = ev.kinetic_energy(st) + ev.trap_energy(st);
            max_drift = std::max(max_drift, std::fabs(E / E0 - 1.0));
        }
        if (s == nsteps) break;
        ev.kinetic_half(st, dt);
        ev.potential_full(st, V, dt);
        ev.kinetic_half(st, dt);
    }
    *energy_drift = max_drift;
    // Amplitude at the final sample, which lands on a maximum: t_end = 2T.
    *amp_err = std::fabs(xs.back() / dx0 - 1.0);
    return kappa_from_crossings(ts, xs);
}

}  // namespace

int main() {
    std::printf("=== rotating-frame validation ===\n\n");
    int failures = 0;

    // ---------------- Part A: coherent epicyclic oscillation ---------------
    const double kappa = 44.2;
    std::printf("[A] coherent state, displaced 0.3 kpc from its guiding "
                "centre, no gravity\n");
    std::printf("%12s %14s %12s %14s %12s\n", "dt", "kappa_meas",
                "rel err", "energy drift", "amp err");
    double err[2], edrift[2], aerr;
    const double fracs[2] = {1.0 / 128.0, 1.0 / 256.0};
    for (int c = 0; c < 2; ++c) {
        const double km = run_coherent(kappa, fracs[c], &edrift[c], &aerr);
        err[c] = std::fabs(km / kappa - 1.0);
        std::printf("%12.6f %14.6f %12.3e %14.3e %12.3e\n",
                    fracs[c] * 2.0 * units::PI / kappa, km, err[c], edrift[c],
                    aerr);
        if (aerr > 1e-3) ++failures;
    }
    // The energy "drift" is the bounded Strang quasi-energy oscillation, so
    // the requirement is that it scales as dt^2, like the frequency error,
    // rather than sitting under an absolute cap at every dt.
    const double order = std::log2(err[0] / err[1]);
    const double eorder = std::log2(edrift[0] / edrift[1]);
    std::printf("  convergence order = %.2f (frequency), %.2f (energy "
                "oscillation); expect ~2\n", order, eorder);
    if (order < 1.7 || order > 2.5 || err[1] > 1e-3) ++failures;
    if (eorder < 1.7 || eorder > 2.5 || edrift[1] > 1e-3) ++failures;

    // ---------------- Part B: warm rotating equilibrium --------------------
    std::printf("\n[B] warm rotating stack + self-gravity, stationarity\n");
    const double Sigma = 50.0 * units::MSUN_PC2;
    Grid2D g(512, 320, 8.0, 5.0);   // dx = dz = 1/64 kpc, production spacing
    Grid1D g1(g.Nz, g.Lz);
    SheetConfig cfg;
    cfg.Sigma = Sigma;
    cfg.sigma_z = 20.0;
    cfg.hbar = 0.6;
    cfg.nu_ext = 0.0;
    cfg.n_states = 300;
    cfg.verbose = false;
    WarmSheet sheet = build_warm_sheet(g1, cfg);
    // Keep the vertical library clear of the FD kinetic band edge; see
    // truncate_sheet_energy.
    truncate_sheet_energy(sheet, 0.8 * 2.0 * cfg.hbar * cfg.hbar /
                                     (g.dz * g.dz));
    std::printf("  vertical library: %d states to %.1f sigma_z^2, "
                "sigma_z = %.2f km/s, h = %.4f kpc\n", sheet.n_significant,
                sheet.energy_cut, sheet.sigma_measured, sheet.h_measured);

    RotConfig rc;
    rc.kappa = kappa;
    rc.sigma_x = 40.0;
    rc.verbose = true;
    RotatingSlab rot = make_warm_rotating(sheet, g, rc);
    SlabState2D& st = rot.st;

    Evolver2D ev(g, cfg.hbar, Kinetic::MatchedFD);
    ev.set_traps(rot.xg, rot.kappa);
    PoissonXZ poisson(g);

    std::vector<double> rho, V, Ve;
    st.density(rho);
    double mass0 = 0.0;
    for (double v : rho) mass0 += v;
    mass0 *= g.dx * g.dz;
    const double mass_err = std::fabs(mass0 / (Sigma * g.Lx) - 1.0);
    std::printf("  mass vs Sigma Lx: rel err = %.3e\n", mass_err);
    if (mass_err > 1e-6) ++failures;

    auto surface_dev = [&](const std::vector<double>& r) {
        double lo = 1e300, hi = -1e300;
        for (int i = 0; i < g.Nx; ++i) {
            double col = 0.0;
            for (int j = 0; j < g.Nz; ++j) col += r[g.idx(i, j)];
            col *= g.dz;
            lo = std::min(lo, col);
            hi = std::max(hi, col);
        }
        return 0.5 * (hi - lo) / (0.5 * (hi + lo));
    };
    auto zrms = [&](const std::vector<double>& r) {
        double m = 0.0, mz2 = 0.0;
        for (int i = 0; i < g.Nx; ++i) {
            for (int j = 0; j < g.Nz; ++j) {
                m += r[g.idx(i, j)];
                mz2 += r[g.idx(i, j)] * g.z(j) * g.z(j);
            }
        }
        return std::sqrt(mz2 / m);
    };

    double kex0, kez0;
    ev.kinetic_energy_split(st, kex0, kez0);
    const double trap0 = ev.trap_energy(st);
    const double sx0 = std::sqrt(2.0 * kex0 / mass0);
    std::printf("  t=0: sigma_x from KE_x = %.3f km/s (stack says %.3f), "
                "KE_x / E_trap = %.5f\n", sx0, rot.sigma_x_realised,
                kex0 / trap0);
    std::printf("  t=0: surface density ripple+speckle = %.3e\n",
                surface_dev(rho));
    if (std::fabs(sx0 / rot.sigma_x_realised - 1.0) > 0.02) ++failures;
    if (std::fabs(kex0 / trap0 - 1.0) > 0.02) ++failures;

    poisson.solve(rho, Ve, false);
    double W0 = 0.0;
    for (int i = 0; i < g.size(); ++i) W0 += rho[i] * Ve[i];
    W0 *= 0.5 * g.dx * g.dz;
    const double E0 = kex0 + kez0 + trap0 + W0;

    const double Tv = 2.0 * units::PI / sheet.omega0;
    const double Te = 2.0 * units::PI / kappa;
    const double t_end = 2.0 * std::max(Tv, Te);
    const double dt = Tv / 200.0;
    const int nsteps = static_cast<int>(t_end / dt);
    const double z0 = zrms(rho);

    std::printf("  evolving %d steps to t = %.3f (%.1f vertical, %.1f "
                "epicyclic periods)\n", nsteps, t_end, t_end / Tv, t_end / Te);
    double max_e = 0.0, max_sx = 0.0, max_z = 0.0, max_dev = 0.0;
    for (int s = 0; s < nsteps; ++s) {
        ev.kinetic_half(st, dt);
        st.density(rho);
        poisson.solve(rho, V);
        ev.potential_full(st, V, dt);
        ev.kinetic_half(st, dt);
        if ((s + 1) % 40 == 0 || s == nsteps - 1) {
            st.density(rho);
            double kex, kez;
            ev.kinetic_energy_split(st, kex, kez);
            poisson.solve(rho, Ve, false);
            double W = 0.0;
            for (int i = 0; i < g.size(); ++i) W += rho[i] * Ve[i];
            W *= 0.5 * g.dx * g.dz;
            const double E = kex + kez + ev.trap_energy(st) + W;
            max_e = std::max(max_e, std::fabs(E / E0 - 1.0));
            max_sx = std::max(
                max_sx, std::fabs(std::sqrt(2.0 * kex / mass0) / sx0 - 1.0));
            max_z = std::max(max_z, std::fabs(zrms(rho) / z0 - 1.0));
            max_dev = std::max(max_dev, surface_dev(rho));
        }
    }
    double spx, spz;
    ev.spectral_spill(st, 0.95 * g.kx_nyquist(), 0.95 * g.kz_nyquist(), spx,
                      spz);
    std::printf("  drift over run: energy %.3e, sigma_x %.3e, z_rms %.3e\n",
                max_e, max_sx, max_z);
    std::printf("  surface ripple max %.3e, spill x %.2e z %.2e\n", max_dev,
                spx, spz);
    if (max_e > 3e-4) ++failures;
    if (max_sx > 0.01) ++failures;
    if (max_z > 0.01) ++failures;
    if (spx > 1e-6 || spz > 1e-6) ++failures;

    std::printf("\n%s\n", failures == 0 ? "PASS: rotating-frame machinery "
                                          "validated"
                                        : "FAIL: see numbers above");
    return failures == 0 ? 0 : 1;
}
