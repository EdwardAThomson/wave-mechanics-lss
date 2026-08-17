// Stage 1: a warm self-gravitating slab in 2D (x, z), periodic in x and
// isolated in z, with a bending perturbation at a chosen wavenumber.
//
// This is where the plan's §6 primary result lives. The corrugation <z>(x) and
// the mean vertical velocity <v_z>(x) are now fields rather than single
// numbers, so their phase relation is measurable: a quarter-cycle offset is a
// propagating bending wave, zero or half a cycle a standing one. That is the
// same discriminator the Cepheid "Great Wave" analysis uses observationally.
//
// Two ways to excite the mode, and they are not equivalent:
//
//   --kick displace   shift the equilibrium by h0 cos(kx), velocities
//                     untouched. A standing wave, clean for reading off the
//                     dispersion relation.
//   --kick velocity   multiply by exp(i v0 cos(kx) z / hbar). Physical, but it
//                     necessarily also imparts an in-plane velocity, because
//                     the Madelung velocity field is a gradient. See the note
//                     in evolve2d.h.
//
// Usage: ./sheet_2d [options]; --help for the list.

#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "src/diagnostics2d.h"
#include "src/equilibrium.h"
#include "src/evolve2d.h"
#include "src/grid2d.h"
#include "src/poisson_xz.h"
#include "src/slab_ic_2d.h"
#include "src/units.h"

namespace {

struct Options {
    int Nx = 512, Nz = 768;
    double Lx = 16.0, Lz = 6.0;
    double Sigma_pc2 = 50.0;
    double sigma_z = 20.0;
    double sigma_x = 10.0;
    double hbar = 0.6;
    double nu_ext = 0.0, z_ext = 0.8;
    int n_states = 300;
    int n_streams = 4;

    std::string kick = "displace";
    int mode = 1;
    double h0 = 0.03;    // kpc, for displace
    double v0 = 5.0;     // km/s, for velocity
    double periods = 3.0;  // in units of the predicted bending period
    int spp = 200;         // steps per VERTICAL period
    int n_outputs = 60;
    std::string kinetic = "spectral";
    std::string outdir = "output/bend_2d";
};

void usage() {
    std::printf(
        "sheet_2d: warm 2D (x, z) slab, bending perturbation, column moments.\n\n"
        "  --Nx <int> --Nz <int>      grid                    (512, 768)\n"
        "  --Lx <kpc> --Lz <kpc>      box                     (16, 6)\n"
        "  --sigma <Msun/pc^2>        surface density         (50)\n"
        "  --sigma-z <km/s>           vertical dispersion     (20)\n"
        "  --sigma-x <km/s>           in-plane dispersion     (10)\n"
        "  --hbar <kpc km/s>          hbar_eff                (0.6)\n"
        "  --nu-ext <km/s/kpc>        rigid halo term         (0)\n"
        "  --states <int>             eigenstates             (300)\n"
        "  --streams <int>            in-plane realisations   (4)\n"
        "  --kick displace|velocity   perturbation type       (displace)\n"
        "  --mode <int>               in-plane mode number    (1)\n"
        "  --h0 <kpc>                 displacement amplitude  (0.03)\n"
        "  --v0 <km/s>                velocity kick amplitude (5)\n"
        "  --periods <float>          bending periods to run  (3)\n"
        "  --spp <int>                steps per vertical period (200)\n"
        "  --outputs <int>            column dumps            (60)\n"
        "  --kinetic spectral|matched                         (spectral)\n"
        "  --out <dir>                output directory\n");
}

bool parse(int argc, char** argv, Options& o) {
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto need = [&](const char* n) -> const char* {
            if (i + 1 >= argc) { std::fprintf(stderr, "missing value for %s\n", n); std::exit(2); }
            return argv[++i];
        };
        if (a == "--help" || a == "-h") { usage(); std::exit(0); }
        else if (a == "--Nx") o.Nx = std::atoi(need("--Nx"));
        else if (a == "--Nz") o.Nz = std::atoi(need("--Nz"));
        else if (a == "--Lx") o.Lx = std::atof(need("--Lx"));
        else if (a == "--Lz") o.Lz = std::atof(need("--Lz"));
        else if (a == "--sigma") o.Sigma_pc2 = std::atof(need("--sigma"));
        else if (a == "--sigma-z") o.sigma_z = std::atof(need("--sigma-z"));
        else if (a == "--sigma-x") o.sigma_x = std::atof(need("--sigma-x"));
        else if (a == "--hbar") o.hbar = std::atof(need("--hbar"));
        else if (a == "--nu-ext") o.nu_ext = std::atof(need("--nu-ext"));
        else if (a == "--z-ext") o.z_ext = std::atof(need("--z-ext"));
        else if (a == "--states") o.n_states = std::atoi(need("--states"));
        else if (a == "--streams") o.n_streams = std::atoi(need("--streams"));
        else if (a == "--kick") o.kick = need("--kick");
        else if (a == "--mode") o.mode = std::atoi(need("--mode"));
        else if (a == "--h0") o.h0 = std::atof(need("--h0"));
        else if (a == "--v0") o.v0 = std::atof(need("--v0"));
        else if (a == "--periods") o.periods = std::atof(need("--periods"));
        else if (a == "--spp") o.spp = std::atoi(need("--spp"));
        else if (a == "--outputs") o.n_outputs = std::atoi(need("--outputs"));
        else if (a == "--kinetic") o.kinetic = need("--kinetic");
        else if (a == "--out") o.outdir = need("--out");
        else { std::fprintf(stderr, "unknown option %s\n", a.c_str()); return false; }
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options o;
    if (!parse(argc, argv, o)) return 2;

    Grid2D g(o.Nx, o.Nz, o.Lx, o.Lz);
    Grid1D g1(o.Nz, o.Lz);
    SheetConfig cfg;
    cfg.Sigma = o.Sigma_pc2 * units::MSUN_PC2;
    cfg.sigma_z = o.sigma_z;
    cfg.hbar = o.hbar;
    cfg.nu_ext = o.nu_ext;
    cfg.z_ext = o.z_ext;
    cfg.n_states = o.n_states;

    std::printf("=== sheet_2d ===\n");
    std::printf("grid %d x %d, Lx = %.2f kpc, Lz = %.2f kpc "
                "(dx = %.5f, dz = %.5f)\n",
                g.Nx, g.Nz, g.Lx, g.Lz, g.dx, g.dz);

    WarmSheet sheet = build_warm_sheet(g1, cfg);

    const Kinetic kin =
        (o.kinetic == "matched") ? Kinetic::MatchedFD : Kinetic::Spectral;
    Evolver2D ev(g, o.hbar, kin);
    PoissonXZ poisson(g);
    SlabState2D st = make_warm_2d(sheet, g, o.sigma_x, o.n_streams);

    const double Sigma = cfg.Sigma;
    const double k = g.kx(o.mode);
    const double h = sheet.h_measured;
    // Razor-thin, in-plane-pressure-corrected prediction. Finite thickness
    // reduces the self-gravity term, so the measured frequency comes out below
    // this; the ratio is the finite-thickness correction.
    const double w2_thin =
        2.0 * units::PI * units::G * Sigma * k - k * k * o.sigma_x * o.sigma_x;
    const double w_thin = (w2_thin > 0.0) ? std::sqrt(w2_thin) : 0.0;

    std::printf("\nbending mode %d: k = %.4f /kpc, lambda = %.2f kpc, "
                "k h = %.3f\n", o.mode, k, 2.0 * units::PI / k, k * h);
    std::printf("razor-thin prediction: omega = %.3f km/s/kpc "
                "(T = %.1f Myr)%s\n",
                w_thin, w_thin > 0 ? 2.0 * units::PI / w_thin * units::TIME_MYR : 0.0,
                w2_thin <= 0.0 ? "   FIREHOSE UNSTABLE" : "");
    std::printf("firehose threshold at this k: sigma_x = %.2f km/s "
                "(currently %.2f)\n",
                std::sqrt(2.0 * units::PI * units::G * Sigma / k), o.sigma_x);

    std::printf("\nresolution audit\n");
    std::printf("  v ceiling x = %8.1f km/s = %5.1f sigma_x\n",
                ev.v_ceiling_x(), ev.v_ceiling_x() / o.sigma_x);
    std::printf("  v ceiling z = %8.1f km/s = %5.1f sigma_z\n",
                ev.v_ceiling_z(), ev.v_ceiling_z() / o.sigma_z);
    if (ev.v_ceiling_x() < 6.0 * o.sigma_x || ev.v_ceiling_z() < 6.0 * o.sigma_z) {
        std::printf("  WARNING: less than 6 sigma of velocity headroom. The\n"
                    "  Maxwellian tail will alias, and aliasing looks exactly\n"
                    "  like heating. Watch the spill columns.\n");
    }

    // ---- perturbation ----
    if (o.kick == "velocity") {
        kick_bending(st, o.v0, k);
        std::printf("\nperturbation: velocity kick v0 = %.2f km/s at mode %d\n",
                    o.v0, o.mode);
    } else {
        std::vector<double> shift(g.Nx);
        for (int i = 0; i < g.Nx; ++i) shift[i] = o.h0 * std::cos(k * g.x(i));
        ev.displace_columns(st, shift);
        std::printf("\nperturbation: displacement h0 = %.4f kpc at mode %d "
                    "(%.2f h)\n", o.h0, o.mode, o.h0 / h);
    }

    std::string mk = "mkdir -p " + o.outdir;
    if (std::system(mk.c_str()) != 0) return 1;

    const double T_vert = 2.0 * units::PI / sheet.omega0;
    const double dt = T_vert / o.spp;
    const double T_bend = (w_thin > 0.0) ? 2.0 * units::PI / w_thin : 10.0 * T_vert;
    const long nsteps = static_cast<long>(o.periods * T_bend / dt);
    const long out_every = std::max<long>(1, nsteps / o.n_outputs);

    {
        FILE* f = std::fopen((o.outdir + "/run_info.dat").c_str(), "w");
        std::fprintf(f, "Nx %d\nNz %d\nLx %.10e\nLz %.10e\n", g.Nx, g.Nz, g.Lx, g.Lz);
        std::fprintf(f, "Sigma %.10e\nsigma_z %.10e\nsigma_x %.10e\nhbar %.10e\n",
                     Sigma, o.sigma_z, o.sigma_x, o.hbar);
        std::fprintf(f, "h_measured %.10e\nomega0 %.10e\nT_vert_Myr %.10e\n",
                     h, sheet.omega0, T_vert * units::TIME_MYR);
        std::fprintf(f, "mode %d\nk %.10e\nomega_thin %.10e\n", o.mode, k, w_thin);
        std::fprintf(f, "n_streams %d\nn_states %zu\nkick %s\n", st.n_streams(),
                     sheet.u.size(), o.kick.c_str());
        std::fclose(f);
    }

    FILE* ser = std::fopen((o.outdir + "/modes.dat").c_str(), "w");
    std::fprintf(ser, "# t[code] t[Myr] Re(z_k) Im(z_k) Re(v_k) Im(v_k) "
                      "|z_k| |v_k| phase_offset[rad] Sigma_tot\n");

    std::printf("\nrunning %ld steps, dt = %.3e (%.4f Myr), "
                "%.1f bending periods\n", nsteps, dt, dt * units::TIME_MYR,
                o.periods);
    std::printf("%10s %10s %12s %12s %12s\n", "step", "t[Myr]", "|<z>_k|[kpc]",
                "|<v>_k|[km/s]", "offset/pi");

    int out_index = 0;
    std::vector<double> rho, V;
    for (long step = 0; step <= nsteps; ++step) {
        const double t = step * dt;
        ColumnMoments cm = compute_columns_2d(ev, st);
        const Complex zk = mode_amplitude(g, cm.zbar, o.mode);
        const Complex vk = mode_amplitude(g, cm.vbar, o.mode);
        double offset = std::arg(vk) - std::arg(zk);
        while (offset > units::PI) offset -= 2.0 * units::PI;
        while (offset < -units::PI) offset += 2.0 * units::PI;
        std::fprintf(ser, "%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n",
                     t, t * units::TIME_MYR, zk.real(), zk.imag(), vk.real(),
                     vk.imag(), std::abs(zk), std::abs(vk), offset,
                     cm.Sigma_total);

        if (step % out_every == 0 || step == nsteps) {
            char name[256];
            std::snprintf(name, sizeof(name), "%s/columns_%04d.dat",
                          o.outdir.c_str(), out_index++);
            write_columns(name, g, cm, t);
            std::printf("%10ld %10.2f %12.6f %12.5f %12.4f\n", step,
                        t * units::TIME_MYR, std::abs(zk), std::abs(vk),
                        offset / units::PI);
            std::fflush(ser);
        }
        if (step == nsteps) break;

        ev.kinetic_half(st, dt);
        st.density(rho);
        poisson.solve(rho, V);
        ev.potential_full(st, V, dt);
        ev.kinetic_half(st, dt);
    }
    std::fclose(ser);
    std::printf("\ndone. output in %s/\n", o.outdir.c_str());
    return 0;
}
