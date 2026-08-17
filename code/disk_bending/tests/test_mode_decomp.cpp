// Stage 1: which frequency is the bending mode?
//
// Tracking <z>(x) gave a frequency insensitive to both k and sigma_x, implying
// a self-gravity term twice the razor-thin value. Something other than the
// bending mode was dominating, and this identifies what.
//
// WHY <z> IS NOT ENOUGH. With a Lagrangian displacement field xi(z), the
// density perturbation is delta_rho = -(rho_0 xi)'. Integrating by parts,
//
//     <z> = Integral z delta_rho dz / Sigma = Integral rho_0 xi dz / Sigma,
//
// the mass-weighted mean displacement. A rigid shift is xi = const, but the
// bending FAMILY is every even xi, so xi = 1, z^2, z^4, ... all move the whole
// column the same way and all contribute to <z>. They are not orthogonal under
// that integral, so <z> mixes them, and since they are not degenerate in
// frequency the mixture beats at several frequencies at once.
//
// (Odd xi = z, z^3, ... give even delta_rho: those are the breathing family,
// and a pure displacement should not excite them. That is checked too.)
//
// WHY THE RAW PROFILE IS NOT ENOUGH EITHER. A first attempt here analysed
// delta_rho(k, z, t) point by point in z. Its t = 0 overlap with the rigid
// displacement shape came out 0.63 instead of 1, because at one z-point the
// in-plane speckle is comparable to the signal. The projection below fixes
// this as a side effect: integrating against a smooth basis function averages
// the broadband z-noise down by roughly the square root of the number of
// independent vertical patches.
//
// SO: project delta_rho onto an orthogonalised bending family, giving one
// clean time series per mode, and frequency-analyse each separately.
#include <algorithm>
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

// delta_rho for displacement field xi(z): -(rho_0 xi)'
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

int main() {
    const double Sigma = 50.0 * units::MSUN_PC2;
    const double sigma_z = 20.0, sigma_x = 40.0, hbar = 0.6;
    // Nx and Nz cut relative to the first pass and the run doubled instead:
    // the mode-amplitude noise is set by the number of in-plane granule
    // patches (Lx / lambda_dB), which is grid-independent, so refining dx
    // buys nothing here while a longer run buys frequency resolution.
    const int Nx = 512, Nz = 384, n_streams = 6;
    const double Lx = 3.0, Lz = 5.0, h0 = 0.05;
    // ONE mode at a time. Exciting three at once put the peak displacement at
    // 81% of the scale height, which is not linear, and the modes then couple.
    const int NMODE = 1;
    const unsigned long long seeds[3] = {8080ULL, 9091ULL, 10102ULL};

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

    const double h = sheet.h_measured;
    const double w0 = sheet.omega0;
    const double k1 = g.kx(1);
    const double w_sg1 = std::sqrt(2.0 * units::PI * units::G * Sigma * k1);

    std::printf("=== Stage 1: bending-family mode decomposition ===\n");
    std::printf("grid %d x %d, Lx = %.1f, Lz = %.1f kpc, %d streams, "
                "h0 = %.3f kpc (%.0f%% of h)\n",
                Nx, Nz, Lx, Lz, n_streams, h0, 100.0 * h0 / h);
    std::printf("h = %.4f kpc, omega_0 = %.2f km/s/kpc; exciting modes 1..%d\n\n",
                h, w0, NMODE);

    // Bending family (even xi): xi = 1, z^2, z^4. Plus one breathing member
    // (odd xi = z) as a control.
    std::vector<std::vector<double>> phi;
    phi.push_back(drho_from_xi(g, sheet.rho, 0));
    phi.push_back(drho_from_xi(g, sheet.rho, 2));
    phi.push_back(drho_from_xi(g, sheet.rho, 4));
    const std::vector<double> breathe = drho_from_xi(g, sheet.rho, 1);

    // Gram-Schmidt, plain L2 in z.
    for (size_t m = 0; m < phi.size(); ++m) {
        for (size_t p = 0; p < m; ++p) {
            const double c = dotr(phi[m], phi[p]) / dotr(phi[p], phi[p]);
            for (int j = 0; j < Nz; ++j) phi[m][j] -= c * phi[p][j];
        }
        const double nrm = std::sqrt(dotr(phi[m], phi[m]));
        for (int j = 0; j < Nz; ++j) phi[m][j] /= nrm;
    }
    std::vector<double> br = breathe;
    {
        const double nrm = std::sqrt(dotr(br, br));
        for (int j = 0; j < Nz; ++j) br[j] /= nrm;
    }

    std::vector<double> shift(g.Nx, 0.0);
    for (int i = 0; i < g.Nx; ++i) shift[i] = h0 * std::cos(g.kx(1) * g.x(i));

    std::printf("REPRODUCIBILITY: same physics, three random-phase seeds.\n"
                "If the extracted frequency is real it must cluster; if it is\n"
                "noise it will scatter.\n\n");
    std::printf("%8s %13s %12s %12s %14s\n", "seed", "omega_rigid",
                "sqrt(2piGSk)", "ratio w^2", "t=0 rigid frac");

    for (int sidx = 0; sidx < 3; ++sidx) {
    SlabState2D st = make_warm_2d(sheet, g, sigma_x, n_streams, seeds[sidx]);
    ev.displace_columns(st, shift);

    // Sixteen periods of the slowest mode, for frequency resolution.
    const double t_end = 8.0 * 2.0 * units::PI / w_sg1;
    const double dt = (2.0 * units::PI / w0) / 250.0;
    const int nsteps = static_cast<int>(t_end / dt);
    const int sample_every = std::max(1, nsteps / 300);

    std::vector<double> ts;
    // a[in-plane mode][vertical channel][t]
    std::vector<std::vector<std::vector<Complex>>> a(
        NMODE, std::vector<std::vector<Complex>>(3));
    std::vector<std::vector<Complex>> abr(NMODE);
    std::vector<double> rho, V;
    std::vector<double> resid_frac;
    for (int step = 0; step <= nsteps; ++step) {
        if (step % sample_every == 0) {
            st.density(rho);
            ts.push_back(step * dt);
            for (int n = 1; n <= NMODE; ++n) {
                std::vector<Complex> p;
                density_profile_mode(g, rho, n, p);
                double tot = 0.0;
                for (int j = 0; j < Nz; ++j) tot += std::norm(p[j]);
                double acc = 0.0;
                for (int m = 0; m < 3; ++m) {
                    const Complex c = dotc(p, phi[m]);
                    a[n - 1][m].push_back(c);
                    acc += std::norm(c);
                }
                abr[n - 1].push_back(dotc(p, br));
                if (n == 1) resid_frac.push_back((tot > 0.0) ? 1.0 - acc / tot : 0.0);
            }
        }
        if (step == nsteps) break;
        ev.kinetic_half(st, dt);
        st.density(rho);
        poisson.solve(rho, V);
        ev.potential_full(st, V, dt);
        ev.kinetic_half(st, dt);
    }
    const int nt = static_cast<int>(ts.size());
    std::printf("sampled %d profiles, t_end = %.4f, frequency resolution "
                "%.2f km/s/kpc\n\n", nt, t_end, 2.0 * units::PI / t_end);

    const double p1 = std::norm(a[0][0][0]);
    const double s2 = p1 + std::norm(a[0][1][0]) + std::norm(a[0][2][0]);

    const double w_max = 4.0 * std::max(w0, w_sg1);
    double frac = 0.0;
    const double wr = dominant_frequency(ts, a[0][0], w_max, &frac);
    std::printf("%8llu %13.2f %12.2f %12.4f %14.4f\n", seeds[sidx], wr, w_sg1,
                wr * wr / (w_sg1 * w_sg1), p1 / s2);
    }

    std::printf("\n  frequency resolution is %.2f km/s/kpc.\n", w_sg1 / 8.0);
    std::printf("  Clustered values mean the bending frequency is measured.\n"
                "  Scattered values mean the rigid-channel time series is\n"
                "  noise-dominated and no dispersion relation can be quoted\n"
                "  from this configuration.\n");
    return 0;
}
