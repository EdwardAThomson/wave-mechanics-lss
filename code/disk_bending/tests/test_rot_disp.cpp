// The stable bending branch, measured in the rotating box.
//
// This is the measurement Stage 1 could not make. Without rotation the clean
// bending regime k << k_J is exactly the fragmenting regime (the duality in
// README.md), and the noise floor at the largest stable box (Lx = 3 kpc) left
// the rigid-channel time series noise-dominated (test_bending_disp,
// test_mode_decomp: three seeds scattered). Rotation at Q = 2.45 holds the
// sheet together at the physical sigma_x = 40 km/s in a 16 kpc box
// (test_toomre), which simultaneously:
//   - reaches k h ~ 0.1, where omega^2 -> 2 pi G Sigma k is a real prediction,
//   - lowers the corrugation noise floor by sqrt(16/3) = 2.3x, since the
//     mode-amplitude noise scales with the in-plane granule count.
//
// Candidate dispersion relations for the measured frequency:
//   razor-thin, gravity only : omega^2 = 2 pi G Sigma k
//   razor-thin, with pressure: omega^2 = 2 pi G Sigma k - k^2 sigma_x^2
// The pressure term assumes stars free-stream through the corrugation; with
// rotation they epicycle instead, and k a_epi ~ 0.34 here, so the term may be
// suppressed. Both ratios are printed; the run decides.
//
// Method as in test_mode_decomp: excite ONE mode as a pure displacement,
// project delta_rho(k1, z, t) onto the Gram-Schmidt orthogonalised bending
// family, periodogram the rigid channel, and repeat over three random-phase
// seeds. Clustered = measured; scattered = noise.
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "../src/equilibrium.h"
#include "../src/evolve2d.h"
#include "../src/grid2d.h"
#include "../src/poisson_xz.h"
#include "../src/rotation.h"
#include "../src/units.h"

namespace {

void density_profile_mode(const Grid2D& g, const std::vector<double>& rho,
                          int n, std::vector<Complex>& out) {
    out.assign(g.Nz, Complex(0.0, 0.0));
    for (int i = 0; i < g.Nx; ++i) {
        const double ph = -g.kx(n) * g.x(i);
        const Complex e(std::cos(ph), std::sin(ph));
        for (int j = 0; j < g.Nz; ++j) out[j] += rho[g.idx(i, j)] * e;
    }
    for (int j = 0; j < g.Nz; ++j) out[j] /= static_cast<double>(g.Nx);
}

double dotr(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}

Complex dotc(const std::vector<Complex>& a, const std::vector<double>& b) {
    Complex s(0.0, 0.0);
    for (size_t i = 0; i < b.size(); ++i) s += a[i] * b[i];
    return s;
}

std::vector<double> drho_from_xi(const Grid2D& g,
                                 const std::vector<double>& rho0, int power) {
    const int N = g.Nz;
    std::vector<double> f(N), out(N, 0.0);
    for (int j = 0; j < N; ++j) f[j] = rho0[j] * std::pow(g.z(j), power);
    for (int j = 1; j < N - 1; ++j) out[j] = -(f[j + 1] - f[j - 1]) / (2.0 * g.dz);
    return out;
}

double dominant_frequency(const std::vector<double>& t,
                          const std::vector<Complex>& y, double w_max,
                          double* frac) {
    const int n = static_cast<int>(y.size());
    const int nw = 1500;
    double best_p = -1.0, best_w = 0.0, tot = 0.0;
    for (int i = 1; i <= nw; ++i) {
        const double w = w_max * i / nw;
        Complex acc(0.0, 0.0);
        for (int j = 0; j < n; ++j) {
            const double ph = w * t[j];
            acc += y[j] * Complex(std::cos(ph), std::sin(ph));
        }
        const double p = std::norm(acc);
        tot += p;
        if (p > best_p) { best_p = p; best_w = w; }
    }
    if (frac) *frac = (tot > 0.0) ? best_p / tot : 0.0;
    return best_w;
}

}  // namespace

int main(int argc, char** argv) {
    // Usage: test_rot_disp [mode] [nseeds]. Defaults: mode 1, 3 seeds. The
    // seed spread at mode 1 came out 0.15%, so higher modes can afford a
    // single seed each when filling in the dispersion relation.
    const int mode = (argc > 1) ? std::atoi(argv[1]) : 1;
    const int nseeds = (argc > 2) ? std::min(3, std::atoi(argv[2])) : 3;
    const double Sigma = 50.0 * units::MSUN_PC2;
    const double sigma_z = 20.0, sigma_x = 40.0, hbar = 0.6, kappa = 44.2;
    Grid2D g(1024, 320, 16.0, 5.0);
    Grid1D g1(g.Nz, g.Lz);
    const double h0 = 0.05;
    const unsigned long long seeds[3] = {8080ULL, 9091ULL, 10102ULL};

    SheetConfig cfg;
    cfg.Sigma = Sigma;
    cfg.sigma_z = sigma_z;
    cfg.hbar = hbar;
    cfg.nu_ext = 0.0;
    cfg.n_states = 300;
    cfg.verbose = false;
    WarmSheet sheet = build_warm_sheet(g1, cfg);
    truncate_sheet_energy(sheet, 0.8 * 2.0 * hbar * hbar / (g.dz * g.dz));

    const double h = sheet.h_measured;
    const double w0 = sheet.omega0;
    const double k1 = g.kx(mode);
    const double w_grav = std::sqrt(2.0 * units::PI * units::G * Sigma * k1);

    std::printf("=== stable bending branch, rotating box ===\n");
    std::printf("grid %d x %d, Lx = %.0f, Lz = %.0f kpc, kappa = %.1f "
                "(Q = %.2f), h0 = %.3f kpc (%.0f%% of h)\n", g.Nx, g.Nz, g.Lx,
                g.Lz, kappa,
                kappa * sigma_x / (3.36 * units::G * Sigma), h0,
                100.0 * h0 / h);
    std::printf("h = %.4f kpc, mode %d: k = %.4f /kpc, k h = %.3f, "
                "omega_0 = %.2f\n", h, mode, k1, k1 * h, w0);

    // Bending family and one breathing control, as in test_mode_decomp.
    std::vector<std::vector<double>> phi;
    phi.push_back(drho_from_xi(g, sheet.rho, 0));
    phi.push_back(drho_from_xi(g, sheet.rho, 2));
    phi.push_back(drho_from_xi(g, sheet.rho, 4));
    for (size_t m = 0; m < phi.size(); ++m) {
        for (size_t p = 0; p < m; ++p) {
            const double c = dotr(phi[m], phi[p]) / dotr(phi[p], phi[p]);
            for (int j = 0; j < g.Nz; ++j) phi[m][j] -= c * phi[p][j];
        }
        const double nrm = std::sqrt(dotr(phi[m], phi[m]));
        for (int j = 0; j < g.Nz; ++j) phi[m][j] /= nrm;
    }

    std::vector<double> shift(g.Nx, 0.0);
    for (int i = 0; i < g.Nx; ++i) shift[i] = h0 * std::cos(k1 * g.x(i));

    double w_press = 0.0;  // set per run from the realised sigma_x
    std::printf("\n%8s %13s %12s %12s %12s %14s %12s\n", "seed", "omega_rigid",
                "w^2/2piGSk", "w^2/(g-p)", "t0 frac", "noise/h0", "energy");

    double wr_all[3] = {0.0, 0.0, 0.0};
    for (int sidx = 0; sidx < nseeds; ++sidx) {
        RotConfig rc;
        rc.kappa = kappa;
        rc.sigma_x = sigma_x;
        rc.seed = seeds[sidx];
        rc.verbose = false;
        RotatingSlab rot = make_warm_rotating(sheet, g, rc);
        SlabState2D st = std::move(rot.st);

        Evolver2D ev(g, hbar, Kinetic::MatchedFD);
        ev.set_traps(rot.xg, rot.kappa);
        PoissonXZ poisson(g);

        const double sxr = rot.sigma_x_realised;
        const double wp2 = w_grav * w_grav - k1 * k1 * sxr * sxr;
        w_press = std::sqrt(std::max(0.0, wp2));

        // The run's own noise floor: the rigid-channel amplitude carried by
        // the IC speckle in the mode about to be excited, over the amplitude
        // the h0 displacement then puts there.
        std::vector<double> rho, V, Ve;
        st.density(rho);
        std::vector<Complex> prof;
        density_profile_mode(g, rho, mode, prof);
        const double noise = std::abs(dotc(prof, phi[0]));
        ev.displace_columns(st, shift);
        st.density(rho);
        density_profile_mode(g, rho, mode, prof);
        const double sig0 = std::abs(dotc(prof, phi[0]));
        const double noise_frac = (sig0 > 0.0) ? noise / sig0 : 0.0;

        const double t_end = 8.0 * 2.0 * units::PI /
                             ((w_press > 0.0) ? w_press : w_grav);
        const double dt = (2.0 * units::PI / w0) / 250.0;
        const int nsteps = static_cast<int>(t_end / dt);
        const int sample_every = std::max(1, nsteps / 400);

        double kex, kez;
        ev.kinetic_energy_split(st, kex, kez);
        poisson.solve(rho, Ve, false);
        double W = 0.0;
        for (int i = 0; i < g.size(); ++i) W += rho[i] * Ve[i];
        W *= 0.5 * g.dx * g.dz;
        const double E0 = kex + kez + ev.trap_energy(st) + W;

        std::vector<double> ts;
        std::vector<std::vector<Complex>> a(3);
        for (int step = 0; step <= nsteps; ++step) {
            if (step % sample_every == 0) {
                st.density(rho);
                ts.push_back(step * dt);
                density_profile_mode(g, rho, mode, prof);
                for (int m = 0; m < 3; ++m) a[m].push_back(dotc(prof, phi[m]));
            }
            if (step == nsteps) break;
            ev.kinetic_half(st, dt);
            st.density(rho);
            poisson.solve(rho, V);
            ev.potential_full(st, V, dt);
            ev.kinetic_half(st, dt);
        }
        st.density(rho);
        ev.kinetic_energy_split(st, kex, kez);
        poisson.solve(rho, Ve, false);
        W = 0.0;
        for (int i = 0; i < g.size(); ++i) W += rho[i] * Ve[i];
        W *= 0.5 * g.dx * g.dz;
        const double e_drift =
            std::fabs((kex + kez + ev.trap_energy(st) + W) / E0 - 1.0);

        const double p1 = std::norm(a[0][0]);
        const double s2 = p1 + std::norm(a[1][0]) + std::norm(a[2][0]);
        const double w_max = 2.0 * w0;
        double frac = 0.0;
        const double wr = dominant_frequency(ts, a[0], w_max, &frac);
        wr_all[sidx] = wr;
        std::printf("%8llu %13.2f %12.4f %12.4f %12.4f %14.4f %12.2e\n",
                    seeds[sidx], wr, wr * wr / (w_grav * w_grav),
                    (wp2 > 0.0) ? wr * wr / wp2 : 0.0, p1 / s2, noise_frac,
                    e_drift);
    }

    double wmean = 0.0;
    for (int i = 0; i < nseeds; ++i) wmean += wr_all[i];
    wmean /= nseeds;
    double wsd = 0.0;
    for (int i = 0; i < nseeds; ++i) {
        wsd += (wr_all[i] - wmean) * (wr_all[i] - wmean);
    }
    wsd = std::sqrt(wsd / nseeds);
    std::printf("\npredictions: sqrt(2piGSk) = %.2f, with pressure = %.2f "
                "km/s/kpc; resolution %.2f\n", w_grav, w_press,
                (w_press > 0.0 ? w_press : w_grav) / 8.0);
    std::printf("measured: %.2f +- %.2f km/s/kpc over %d seed(s)\n", wmean,
                wsd, nseeds);
    std::printf("clustered (sd well under the seed spread of the old "
                "noise-dominated runs) means the\nstable bending branch is "
                "measured; scattered means it is still noise.\n");
    return 0;
}
