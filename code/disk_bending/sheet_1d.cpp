// Stage 0 of the disk bending-wave ladder: a self-gravitating vertical slab,
// warm, isolated in z, kicked impulsively, diagnosed in phase space.
//
// What this can and cannot do
// ---------------------------
// It cannot do bending waves. A slab with no lateral coordinate has no k_x, so
// there is no bending mode, no dispersion relation and no firehose here: those
// need Stage 1, the 2D (x, z) box. What it *can* do is the vertical phase
// spiral, which is a genuinely 1D object in (z, v_z) and is, by the plan's own
// reckoning, the strongest single result available from the project. It also
// carries the bending/breathing distinction in the form of the parity of the
// kick, and it is where the resolution trade between hbar_eff, grid size and
// stream count can be measured cheaply before it becomes expensive.
//
// One physics point that changes the setup. An isolated self-gravitating slab
// given a uniform vertical kick just translates forever: nothing restores it,
// because momentum is conserved and there is no external frame. A real disk is
// kicked relative to a halo that stays put. So a rigid external vertical
// potential is not optional decoration here, it is what makes a bulk kick phase
// mix into a spiral at all. Set --nu-ext 0 to see the degenerate case.
//
// Usage:
//   ./sheet_1d [options]
// Run with --help for the full list.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "src/diagnostics.h"
#include "src/equilibrium.h"
#include "src/evolve.h"
#include "src/grid.h"
#include "src/husimi.h"
#include "src/slab_ic.h"
#include "src/units.h"

namespace {

struct Options {
    int N = 2048;
    double L = 4.0;
    double Sigma_pc2 = 50.0;   // Msun/pc^2
    double sigma_z = 20.0;     // km/s
    double hbar = 0.3;         // kpc km/s
    double nu_ext = 40.0;      // km/s/kpc, rigid halo term
    double z_ext = 0.8;        // kpc, where the rigid term flattens
    int n_states = 320;

    std::string ic = "multistream";  // multistream | single | batched
    int n_streams = 16;              // for ic = batched

    std::string kick = "bulk";  // bulk | compress | none
    double v_kick = 12.0;       // km/s for bulk
    double a_kick = 30.0;       // km/s/kpc for compress

    double n_periods = 12.0;
    int steps_per_period = 800;
    int n_outputs = 60;
    std::string kinetic = "matched";  // matched | spectral
    std::string outdir = "output/sheet_1d";
    int husimi_every = 10;
};

void usage() {
    std::printf(
        "sheet_1d: warm self-gravitating vertical slab, impulsive kick,\n"
        "          Husimi phase-space diagnostics.\n\n"
        "  --N <int>            grid points in z            (1024)\n"
        "  --L <kpc>            box height                  (4.0)\n"
        "  --sigma <Msun/pc^2>  stellar surface density     (50)\n"
        "  --sigma-z <km/s>     vertical dispersion         (20)\n"
        "  --hbar <kpc km/s>    hbar_eff                    (0.3)\n"
        "  --nu-ext <km/s/kpc>  rigid external vertical frequency (40)\n"
        "  --z-ext <kpc>        height where the rigid term flattens (0.8);\n"
        "                       make it huge for a purely harmonic well, which\n"
        "                       by Kohn's theorem gives no spiral at all\n"
        "  --states <int>       eigenstates in the library  (320)\n"
        "  --ic <mode>          multistream | single | batched (multistream)\n"
        "  --streams <int>      streams when ic = batched   (16)\n"
        "  --kick <mode>        bulk | compress | none      (bulk)\n"
        "  --vkick <km/s>       bulk kick amplitude         (12)\n"
        "  --akick <km/s/kpc>   compress kick amplitude     (30)\n"
        "  --periods <float>    vertical periods to run     (12)\n"
        "  --spp <int>          steps per vertical period   (800)\n"
        "  --outputs <int>      profile dumps               (60)\n"
        "  --kinetic <mode>     matched | spectral          (matched)\n"
        "  --husimi-every <int> dump a Husimi map every Nth output (10)\n"
        "  --out <dir>          output directory\n");
}

bool parse(int argc, char** argv, Options& o) {
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto need = [&](const char* name) -> const char* {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "missing value for %s\n", name);
                std::exit(2);
            }
            return argv[++i];
        };
        if (a == "--help" || a == "-h") { usage(); std::exit(0); }
        else if (a == "--N") o.N = std::atoi(need("--N"));
        else if (a == "--L") o.L = std::atof(need("--L"));
        else if (a == "--sigma") o.Sigma_pc2 = std::atof(need("--sigma"));
        else if (a == "--sigma-z") o.sigma_z = std::atof(need("--sigma-z"));
        else if (a == "--hbar") o.hbar = std::atof(need("--hbar"));
        else if (a == "--nu-ext") o.nu_ext = std::atof(need("--nu-ext"));
        else if (a == "--z-ext") o.z_ext = std::atof(need("--z-ext"));
        else if (a == "--states") o.n_states = std::atoi(need("--states"));
        else if (a == "--ic") o.ic = need("--ic");
        else if (a == "--streams") o.n_streams = std::atoi(need("--streams"));
        else if (a == "--kick") o.kick = need("--kick");
        else if (a == "--vkick") o.v_kick = std::atof(need("--vkick"));
        else if (a == "--akick") o.a_kick = std::atof(need("--akick"));
        else if (a == "--periods") o.n_periods = std::atof(need("--periods"));
        else if (a == "--spp") o.steps_per_period = std::atoi(need("--spp"));
        else if (a == "--outputs") o.n_outputs = std::atoi(need("--outputs"));
        else if (a == "--kinetic") o.kinetic = need("--kinetic");
        else if (a == "--husimi-every") o.husimi_every = std::atoi(need("--husimi-every"));
        else if (a == "--out") o.outdir = need("--out");
        else { std::fprintf(stderr, "unknown option %s\n", a.c_str()); return false; }
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options o;
    if (!parse(argc, argv, o)) return 2;

    Grid1D g(o.N, o.L);
    SheetConfig cfg;
    cfg.Sigma = o.Sigma_pc2 * units::MSUN_PC2;
    cfg.sigma_z = o.sigma_z;
    cfg.hbar = o.hbar;
    cfg.nu_ext = o.nu_ext;
    cfg.z_ext = o.z_ext;
    cfg.n_states = o.n_states;

    std::printf("=== sheet_1d ===\n");
    std::printf("grid: N = %d, L = %.2f kpc, dz = %.5f kpc\n", g.N, g.L, g.dz);

    WarmSheet sheet = build_warm_sheet(g, cfg);

    const Kinetic kin =
        (o.kinetic == "spectral") ? Kinetic::Spectral : Kinetic::MatchedFD;
    Evolver ev(g, cfg.hbar, kin);

    SlabState st;
    if (o.ic == "single") st = make_single_psi(sheet);
    else if (o.ic == "batched") st = make_batched(sheet, o.n_streams);
    else st = make_multistream(sheet);

    // ---- resolution audit, before anything is evolved ----
    const double v_ceiling = ev.v_max();
    // Estimate the velocity the run will actually reach: 3.3 sigma_z holds
    // 99.9% of the mass, plus whichever kick is being applied. Using the
    // formal top of the state library instead would quote a speed carried by
    // exp(-12) of the mass and would warn on runs that are perfectly fine.
    const double v_kick_used = (o.kick == "bulk")       ? std::fabs(o.v_kick)
                               : (o.kick == "compress") ? std::fabs(o.a_kick) * 0.5 * o.L
                                                        : 0.0;
    const double v_expected = 3.3 * cfg.sigma_z + v_kick_used;
    std::printf("\nresolution audit\n");
    std::printf("  streams              = %d\n", st.n_streams());
    std::printf("  v ceiling (this run) = %.1f km/s  (%s)\n", v_ceiling,
                kin == Kinetic::Spectral ? "spectral" : "matched FD");
    std::printf("  v ceiling (spectral) = %.1f km/s\n",
                cfg.hbar * g.k_nyquist());
    std::printf("  v ceiling (matchedFD)= %.1f km/s   "
                "(lower by pi: the FD group\n"
                "                                      velocity turns over at "
                "k dz = pi/2)\n",
                cfg.hbar / g.dz);
    std::printf("  v expected in run    = %.1f km/s  "
                "(3.3 sigma_z + kick)\n", v_expected);
    std::printf("  headroom             = %.2fx\n", v_ceiling / v_expected);
    if (v_ceiling < 2.0 * v_expected) {
        std::printf("  WARNING: less than 2x headroom. This is the same failure\n"
                    "  mode as the cosmological runs: velocities above the\n"
                    "  ceiling alias back with the wrong sign instead of being\n"
                    "  resolved. Raise N, or raise hbar_eff, or use --kinetic\n"
                    "  spectral for a factor pi more ceiling. Watch the spill\n"
                    "  column below: if it grows, the run is contaminated.\n");
    }

    const double T = 2.0 * units::PI / sheet.omega0;
    const double dt = T / o.steps_per_period;
    const double dt_kin = 2.0 * g.dz * g.dz / (units::PI * cfg.hbar);
    double vmaxpot = 0.0;
    for (int j = 0; j < g.N; ++j) {
        vmaxpot = std::max(vmaxpot,
                           std::fabs(sheet.Phi_self[j] + sheet.Phi_ext[j]));
    }
    const double dt_pot = units::PI * cfg.hbar / vmaxpot;
    std::printf("  dt                   = %.4e (%.4f Myr)\n", dt,
                dt * units::TIME_MYR);
    std::printf("  dt kinetic bound     = %.4e  (Nyquist; conservative)\n", dt_kin);
    std::printf("  dt potential bound   = %.4e\n", dt_pot);
    if (dt > dt_pot) {
        std::printf("  WARNING: dt exceeds the potential phase bound.\n");
    }

    // ---- reference: the unperturbed equilibrium in phase space ----
    // Dumped BEFORE the kick. Dividing later maps by this is the only
    // assumption-free way to show the perturbation: an azimuthal average in
    // (z, v/omega) would additionally assume the equilibrium is circular
    // there, and it is not, because the potential is anharmonic and orbits at
    // larger amplitude have lower frequency. That mismatch shows up as a
    // spurious m = 2 pattern that swamps the spiral.
    std::string mkdir_cmd = "mkdir -p " + o.outdir;
    if (std::system(mkdir_cmd.c_str()) != 0) {
        std::fprintf(stderr, "could not create %s\n", o.outdir.c_str());
        return 1;
    }
    {
        HusimiConfig hc0;
        const double vw = 4.0 * cfg.sigma_z +
                          ((o.kick == "bulk") ? std::fabs(o.v_kick) : 0.0);
        hc0.v_half = std::min(vw, 0.9 * v_ceiling);
        hc0.z_half = std::min(hc0.v_half / sheet.omega0, 0.45 * o.L);
        hc0.nz = 256;
        HusimiMap eq = husimi_transform(ev, st, hc0, sheet.h_measured,
                                        cfg.sigma_z);
        write_husimi(o.outdir + "/husimi_equilibrium.dat", eq, 0.0);
    }

    // ---- kick ----
    if (o.kick == "bulk") {
        kick_bulk(st, o.v_kick);
        std::printf("\nkick: bulk, v0 = %.2f km/s (odd -> bending analogue)\n",
                    o.v_kick);
    } else if (o.kick == "compress") {
        kick_compress(st, o.a_kick);
        std::printf("\nkick: compress, dv/dz = %.2f km/s/kpc "
                    "(even -> breathing analogue)\n",
                    o.a_kick);
    } else {
        std::printf("\nkick: none (stationarity run)\n");
    }
    if (o.kick == "bulk" && cfg.nu_ext == 0.0) {
        std::printf("  NOTE: nu_ext = 0, so nothing restores a bulk kick. The\n"
                    "  slab will translate at constant velocity and no spiral\n"
                    "  will form. This is physics, not a bug.\n");
    }

    // ---- run ----
    const long total_steps = static_cast<long>(o.n_periods * o.steps_per_period);
    const long out_every = std::max<long>(1, total_steps / o.n_outputs);

    PotentialWorkspace ws;
    const double k_cut = 0.5 * g.k_nyquist();

    // Machine-readable run parameters, so the plotting scripts do not have to
    // re-derive omega_z by eye from the moment history.
    {
        FILE* info = std::fopen((o.outdir + "/run_info.dat").c_str(), "w");
        std::fprintf(info, "omega0 %.10e\n", sheet.omega0);
        std::fprintf(info, "h_measured %.10e\n", sheet.h_measured);
        std::fprintf(info, "sigma_z %.10e\n", cfg.sigma_z);
        std::fprintf(info, "sigma_measured %.10e\n", sheet.sigma_measured);
        std::fprintf(info, "hbar %.10e\n", cfg.hbar);
        std::fprintf(info, "Sigma %.10e\n", cfg.Sigma);
        std::fprintf(info, "nu_ext %.10e\n", cfg.nu_ext);
        std::fprintf(info, "z_ext %.10e\n", cfg.z_ext);
        std::fprintf(info, "n_streams %d\n", st.n_streams());
        std::fprintf(info, "n_states %zu\n", sheet.u.size());
        std::fprintf(info, "v_ceiling %.10e\n", v_ceiling);
        std::fprintf(info, "T_vert_Myr %.10e\n", T * units::TIME_MYR);
        std::fprintf(info, "kick %s\n", o.kick.c_str());
        std::fprintf(info, "v_kick %.10e\n", o.v_kick);
        std::fprintf(info, "a_kick %.10e\n", o.a_kick);
        std::fprintf(info, "ic %s\n", o.ic.c_str());
        std::fclose(info);
    }

    FILE* series = std::fopen((o.outdir + "/moments.dat").c_str(), "w");
    std::fprintf(series,
                 "# t[code] t[Myr] Sigma zbar[kpc] vbar[km/s] z_rms v_rms "
                 "KE PE_self PE_ext E spill\n");

    // Match the phase-space window to the populated region, and make it square
    // in (z, v/omega) so that unperturbed orbits are circles. A window much
    // larger than the distribution makes the azimuthal-average contrast
    // meaningless in the empty annuli, which is where spurious structure comes
    // from in these plots.
    HusimiConfig hc;
    const double v_window = 4.0 * cfg.sigma_z + v_kick_used;
    hc.v_half = std::min(v_window, 0.9 * v_ceiling);
    hc.z_half = std::min(hc.v_half / sheet.omega0, 0.45 * o.L);
    hc.nz = 256;
    std::printf("\nHusimi window: |z| < %.3f kpc, |v_z| < %.1f km/s "
                "(square in z, v/omega)\n", hc.z_half, hc.v_half);

    FILE* wind = std::fopen((o.outdir + "/winding.dat").c_str(), "w");
    std::fprintf(wind, "# t[Myr]  radius[kpc]  phase[rad]  m1_amplitude\n");
    FILE* wsum = std::fopen((o.outdir + "/winding_summary.dat").c_str(), "w");
    std::fprintf(wsum, "# t[code] t[Myr] turns\n");

    std::printf("\nrunning %ld steps (%.1f vertical periods, %.1f Myr)\n",
                total_steps, o.n_periods, o.n_periods * T * units::TIME_MYR);
    std::printf("%10s %10s %12s %12s %12s %12s %10s %10s\n", "step", "t[Myr]",
                "<z>[kpc]", "<v>[km/s]", "z_rms", "dE/E", "spill", "edge");

    Moments m0 = compute_moments(ev, st, sheet.Phi_ext);
    int out_index = 0;
    bool edge_warned = false;

    for (long step = 0; step <= total_steps; ++step) {
        if (step % out_every == 0 || step == total_steps) {
            const double t = step * dt;
            Moments m = compute_moments(ev, st, sheet.Phi_ext);
            const double spill = ev.spectral_spill(st, k_cut);
            std::fprintf(series,
                         "%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e "
                         "%.8e %.8e\n",
                         t, t * units::TIME_MYR, m.Sigma, m.zbar, m.vbar,
                         m.z_rms, m.v_rms, m.ke, m.pe_self, m.pe_ext, m.energy,
                         spill);
            std::fflush(series);

            char name[256];
            std::snprintf(name, sizeof(name), "%s/profile_%04d.dat",
                          o.outdir.c_str(), out_index);
            write_profile(name, ev, st, sheet.Phi_ext, t);

            if (out_index % o.husimi_every == 0) {
                HusimiMap hm = husimi_transform(ev, st, hc, sheet.h_measured,
                                                cfg.sigma_z);
                std::snprintf(name, sizeof(name), "%s/husimi_%04d.dat",
                              o.outdir.c_str(), out_index);
                write_husimi(name, hm, t);

                WindingProfile wp = phase_spiral_winding(hm, sheet.omega0);
                for (size_t i = 0; i < wp.radius.size(); ++i) {
                    std::fprintf(wind, "%.6e %.6e %.6e %.6e\n",
                                 t * units::TIME_MYR, wp.radius[i], wp.phase[i],
                                 wp.power[i]);
                }
                std::fprintf(wind, "\n");
                std::fprintf(wsum, "%.8e %.8e %.8e\n", t, t * units::TIME_MYR,
                             wp.turns);
                std::fflush(wind);
                std::fflush(wsum);
            }

            std::printf("%10ld %10.2f %12.5f %12.4f %12.5f %12.3e %10.2e "
                        "%10.2e\n",
                        step, t * units::TIME_MYR, m.zbar, m.vbar, m.z_rms,
                        std::fabs((m.energy - m0.energy) / m0.energy), spill,
                        m.edge_mass);
            if (m.edge_mass > 1e-6 && !edge_warned) {
                std::printf("  WARNING: %.2e of the mass has reached the outer "
                            "fifth of the box.\n"
                            "  The kinetic step is periodic, so once the "
                            "kicked sheet swings out this\n"
                            "  far it wraps around and stops agreeing with the "
                            "isolated-BC gravity.\n"
                            "  Raise --L. The box has to contain the KICKED "
                            "orbit, not just the\n"
                            "  equilibrium sheet.\n",
                            m.edge_mass);
                edge_warned = true;
            }
            ++out_index;
        }
        if (step == total_steps) break;
        strang_step(ev, st, sheet.Phi_ext, ws, dt);
    }

    std::fclose(series);
    std::fclose(wind);
    std::fclose(wsum);
    std::printf("\ndone. output in %s/\n", o.outdir.c_str());
    return 0;
}
