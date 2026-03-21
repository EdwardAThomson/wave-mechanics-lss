// 3D Schrodinger-Poisson solver
// Goldberg/Cayley method with periodic BCs and Strang splitting
//
// Usage:
//   ./sp_3d           Free particle: 3D Gaussian disperses
//   ./sp_3d 5000      Gravity on: Gaussian held by self-gravity
//
// The kinetic half-step uses splitting operators (thesis section 5.3.1):
//   exp(-iKdt/2) = exp(-iK_x*dt/2) * exp(-iK_y*dt/2) * exp(-iK_z*dt/2)
// since [K_x, K_y] = [K_x, K_z] = [K_y, K_z] = 0.
// Each 1D sweep uses the Goldberg recursion with double-pass periodic BCs.
//
// Reference: Thomson (2011) Chapter 5

#include <cmath>
#include <complex>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <chrono>
#include <fftw3.h>

using Complex = std::complex<double>;

// === Grid ===
constexpr int N = 64;
constexpr int N3 = N * N * N;
constexpr double BOX = 1.0;
constexpr double epsilon = BOX / N;
constexpr double PI = M_PI;

// === Timestep (adjustable per IC type in main) ===
constexpr double delta = 1e-4;
int N_STEPS = 200;
int OUTPUT_EVERY = 50;

// lambda = 4 * epsilon^2 / delta  (split-operator kinetic half-step)
const double LAMBDA = 4.0 * epsilon * epsilon / delta;

// === Initial condition parameters ===
constexpr double ic_cx = 0.5, ic_cy = 0.5, ic_cz = 0.5;  // centre

// Gaussian IC
constexpr double w0 = 0.08;               // width
constexpr double px0 = 0.0, py0 = 0.0, pz0 = 0.0;  // initial momentum

// Tophat IC
constexpr double R_tophat = 0.2;           // radius (box units)
constexpr double delta_tophat = 0.5;       // fractional overdensity
constexpr double sigma_smooth = 3.0;       // smoothing width in grid cells

// === Gravitational coupling (command line, 0 = free particle) ===
double G_EFF = 0.0;

// ================================================================
// 3D index (row-major, z fastest)
// ================================================================
inline int idx(int ix, int iy, int iz) {
    return (ix * N + iy) * N + iz;
}

inline int wrap(int j) {
    return ((j % N) + N) % N;
}

// ================================================================
// 1D Goldberg/Cayley sweep on a pencil of length N
//
// This is the same recursion as in sp_1d.cpp but operates on
// an arbitrary Complex array rather than the full grid.
// Called once per pencil, per dimension, per kinetic half-step.
// ================================================================
void goldberg_sweep(Complex* pencil) {
    const Complex i_lam(0.0, LAMBDA);
    const Complex alpha = 2.0 - i_lam;
    const Complex beta  = 2.0 + i_lam;

    Complex Omega[N], e[N], f[N];

    // Omega (periodic indices)
    for (int j = 0; j < N; j++) {
        Omega[j] = -pencil[wrap(j + 1)] + beta * pencil[j] - pencil[wrap(j - 1)];
    }

    // e: double forward pass
    e[0] = alpha;
    for (int j = 1; j < N; j++) e[j] = alpha - 1.0 / e[j - 1];
    e[0] = alpha - 1.0 / e[N - 1];
    for (int j = 1; j < N; j++) e[j] = alpha - 1.0 / e[j - 1];

    // f: double forward pass
    f[0] = Omega[0];
    for (int j = 1; j < N; j++) f[j] = Omega[j] + f[j - 1] / e[j - 1];
    f[0] = Omega[0] + f[N - 1] / e[N - 1];
    for (int j = 1; j < N; j++) f[j] = Omega[j] + f[j - 1] / e[j - 1];

    // psi: double backward pass
    pencil[N - 1] = (pencil[0] - f[N - 1]) / e[N - 1];
    for (int j = N - 2; j >= 0; j--) pencil[j] = (pencil[j + 1] - f[j]) / e[j];
    pencil[N - 1] = (pencil[0] - f[N - 1]) / e[N - 1];
    for (int j = N - 2; j >= 0; j--) pencil[j] = (pencil[j + 1] - f[j]) / e[j];
}

// ================================================================
// Kinetic half-step in 3D: sequential sweeps along x, y, z
//
// For each dimension, we extract every 1D "pencil" along that axis,
// apply the Goldberg recursion, and write back.
//
// Since [K_x, K_y, K_z] all commute, the order doesn't matter
// (thesis section 5.3.1, Eq 5.37).
// ================================================================
void kinetic_half_step_3d(std::vector<Complex>& psi) {
    Complex pencil[N];

    // x-sweep: fix (iy, iz), vary ix
    for (int iy = 0; iy < N; iy++) {
        for (int iz = 0; iz < N; iz++) {
            for (int ix = 0; ix < N; ix++) pencil[ix] = psi[idx(ix, iy, iz)];
            goldberg_sweep(pencil);
            for (int ix = 0; ix < N; ix++) psi[idx(ix, iy, iz)] = pencil[ix];
        }
    }

    // y-sweep: fix (ix, iz), vary iy
    for (int ix = 0; ix < N; ix++) {
        for (int iz = 0; iz < N; iz++) {
            for (int iy = 0; iy < N; iy++) pencil[iy] = psi[idx(ix, iy, iz)];
            goldberg_sweep(pencil);
            for (int iy = 0; iy < N; iy++) psi[idx(ix, iy, iz)] = pencil[iy];
        }
    }

    // z-sweep: fix (ix, iy), vary iz — contiguous in memory
    for (int ix = 0; ix < N; ix++) {
        for (int iy = 0; iy < N; iy++) {
            int base = idx(ix, iy, 0);
            for (int iz = 0; iz < N; iz++) pencil[iz] = psi[base + iz];
            goldberg_sweep(pencil);
            for (int iz = 0; iz < N; iz++) psi[base + iz] = pencil[iz];
        }
    }
}

// ================================================================
// 3D Poisson solve: nabla^2 V = |psi|^2 - <|psi|^2>
//
// Discrete Green's function (thesis section 5.5):
//   V_hat = eps^2 * rho_hat / (2*kappa - 6)
//   kappa = cos(2*pi*n/N) + cos(2*pi*m/N) + cos(2*pi*o/N)
// ================================================================
void poisson_solve_3d(const std::vector<Complex>& psi,
                      std::vector<double>& V,
                      fftw_complex* work_in,
                      fftw_complex* work_out,
                      fftw_plan& plan_fwd,
                      fftw_plan& plan_bwd) {
    // Mean-subtracted density
    double rho_mean = 0.0;
    for (int i = 0; i < N3; i++) rho_mean += std::norm(psi[i]);
    rho_mean /= N3;

    for (int i = 0; i < N3; i++) {
        work_in[i][0] = std::norm(psi[i]) - rho_mean;
        work_in[i][1] = 0.0;
    }

    // Forward FFT
    fftw_execute(plan_fwd);

    // Apply discrete Green's function
    work_out[0][0] = 0.0;
    work_out[0][1] = 0.0;
    for (int ix = 0; ix < N; ix++) {
        double cx = cos(2.0 * PI * ix / N);
        for (int iy = 0; iy < N; iy++) {
            double cy = cos(2.0 * PI * iy / N);
            for (int iz = 0; iz < N; iz++) {
                if (ix == 0 && iy == 0 && iz == 0) continue;
                double cz = cos(2.0 * PI * iz / N);
                double eigenvalue = 2.0 * (cx + cy + cz) - 6.0;
                double factor = epsilon * epsilon / eigenvalue;
                int i = idx(ix, iy, iz);
                work_out[i][0] *= factor;
                work_out[i][1] *= factor;
            }
        }
    }

    // Backward FFT
    fftw_execute(plan_bwd);

    for (int i = 0; i < N3; i++) {
        V[i] = work_in[i][0] / static_cast<double>(N3);
    }
}

// ================================================================
// Potential step: pointwise Cayley rotation (thesis Eq 5.36)
// ================================================================
void potential_step_3d(std::vector<Complex>& psi,
                       const std::vector<double>& V) {
    for (int i = 0; i < N3; i++) {
        double vdt2 = G_EFF * V[i] * delta / 2.0;
        Complex num(1.0, -vdt2);
        Complex den(1.0,  vdt2);
        psi[i] *= num / den;
    }
}

// ================================================================
// Output: z-slice at iz = N/2
// ================================================================
void write_slice(const std::string& filename,
                 const std::vector<Complex>& psi,
                 const std::vector<double>& V,
                 int step, double time) {
    std::ofstream file(filename);
    file << "# step = " << step << ", time = " << time
         << ", G_eff = " << G_EFF << ", N = " << N
         << ", iz = " << N / 2 << "\n";
    file << "# x  y  rho  V\n";
    int iz = N / 2;
    for (int ix = 0; ix < N; ix++) {
        double x = (ix + 0.5) * epsilon;
        for (int iy = 0; iy < N; iy++) {
            double y = (iy + 0.5) * epsilon;
            int i = idx(ix, iy, iz);
            file << x << "  " << y << "  "
                 << std::norm(psi[i]) << "  " << V[i] << "\n";
        }
        file << "\n";  // blank line between x-rows (gnuplot compatible)
    }
    file.close();
}

// ================================================================
// Output: spherically averaged radial density profile
// ================================================================
void write_radial_profile(const std::string& filename,
                          const std::vector<Complex>& psi,
                          int step, double time) {
    constexpr int N_BINS = N / 2;
    double bin_width = 0.5 / N_BINS;   // bins from r=0 to r=0.5

    std::vector<double> rho_sum(N_BINS, 0.0);
    std::vector<int> count(N_BINS, 0);

    for (int ix = 0; ix < N; ix++) {
        double dx = (ix + 0.5) * epsilon - ic_cx;
        if (dx > 0.5) dx -= 1.0;
        if (dx < -0.5) dx += 1.0;
        for (int iy = 0; iy < N; iy++) {
            double dy = (iy + 0.5) * epsilon - ic_cy;
            if (dy > 0.5) dy -= 1.0;
            if (dy < -0.5) dy += 1.0;
            for (int iz = 0; iz < N; iz++) {
                double dz = (iz + 0.5) * epsilon - ic_cz;
                if (dz > 0.5) dz -= 1.0;
                if (dz < -0.5) dz += 1.0;
                double r = sqrt(dx * dx + dy * dy + dz * dz);
                int bin = static_cast<int>(r / bin_width);
                if (bin < N_BINS) {
                    rho_sum[bin] += std::norm(psi[idx(ix, iy, iz)]);
                    count[bin]++;
                }
            }
        }
    }

    std::ofstream file(filename);
    file << "# step = " << step << ", time = " << time << "\n";
    file << "# r  rho_avg  count\n";
    for (int b = 0; b < N_BINS; b++) {
        double r = (b + 0.5) * bin_width;
        double rho_avg = (count[b] > 0) ? rho_sum[b] / count[b] : 0.0;
        file << r << "  " << rho_avg << "  " << count[b] << "\n";
    }
    file.close();
}

// ================================================================
// Main
// ================================================================
int main(int argc, char* argv[]) {
    // Usage: ./sp_3d [G_eff] [gaussian|tophat]
    if (argc >= 2) G_EFF = std::stod(argv[1]);
    std::string ic_type = "gaussian";
    if (argc >= 3) ic_type = argv[2];
    bool gravity = (G_EFF != 0.0);
    bool tophat = (ic_type == "tophat");

    // Adjust step count for tophat (need more outputs to see collapse)
    if (tophat) {
        N_STEPS = 500;
        OUTPUT_EVERY = 25;
    }

    std::cout << "3D Schrodinger-Poisson solver\n";
    std::cout << "N = " << N << "^3 = " << N3
              << ", epsilon = " << epsilon << "\n";
    std::cout << "delta = " << delta << ", lambda = " << LAMBDA << "\n";
    std::cout << "Total time = " << N_STEPS * delta << "\n";
    std::cout << "G_eff = " << G_EFF
              << (gravity ? " (gravity ON)" : " (free particle)") << "\n";
    if (tophat) {
        std::cout << "IC: Tophat R=" << R_tophat << " delta=" << delta_tophat
                  << " smooth=" << sigma_smooth << " cells\n\n";
    } else {
        std::cout << "IC: 3D Gaussian w=" << w0 << " at ("
                  << ic_cx << "," << ic_cy << "," << ic_cz << ")\n\n";
    }

    // --- Initialise wavefunction ---
    std::vector<Complex> psi(N3);

    if (tophat) {
        // Tophat: uniform overdensity sphere with smooth edges
        // rho_in = 1 + delta_tophat (inside), rho_out adjusted for mass conservation
        double V_sphere = (4.0 / 3.0) * PI * R_tophat * R_tophat * R_tophat;
        double rho_in = 1.0 + delta_tophat;
        double rho_out = (1.0 - rho_in * V_sphere) / (1.0 - V_sphere);
        double smooth_width = sigma_smooth * epsilon;

        std::cout << "Tophat: rho_in = " << rho_in << ", rho_out = " << rho_out
                  << ", V_sphere = " << V_sphere << "\n";

        for (int ix = 0; ix < N; ix++) {
            double dx = (ix + 0.5) * epsilon - ic_cx;
            if (dx > 0.5) dx -= 1.0;
            if (dx < -0.5) dx += 1.0;
            for (int iy = 0; iy < N; iy++) {
                double dy = (iy + 0.5) * epsilon - ic_cy;
                if (dy > 0.5) dy -= 1.0;
                if (dy < -0.5) dy += 1.0;
                for (int iz = 0; iz < N; iz++) {
                    double dz = (iz + 0.5) * epsilon - ic_cz;
                    if (dz > 0.5) dz -= 1.0;
                    if (dz < -0.5) dz += 1.0;
                    double r = sqrt(dx * dx + dy * dy + dz * dz);
                    // Smooth transition: tanh profile at r = R_tophat
                    double profile = 0.5 * (1.0 - tanh((r - R_tophat) / smooth_width));
                    double rho = rho_out + (rho_in - rho_out) * profile;
                    // psi = sqrt(rho), real (no initial velocity)
                    psi[idx(ix, iy, iz)] = Complex(sqrt(std::max(rho, 0.0)), 0.0);
                }
            }
        }
    } else {
        // Gaussian wavepacket
        for (int ix = 0; ix < N; ix++) {
            double x = (ix + 0.5) * epsilon;
            double dx = x - ic_cx;
            if (dx > 0.5) dx -= 1.0;
            if (dx < -0.5) dx += 1.0;
            for (int iy = 0; iy < N; iy++) {
                double y = (iy + 0.5) * epsilon;
                double dy = y - ic_cy;
                if (dy > 0.5) dy -= 1.0;
                if (dy < -0.5) dy += 1.0;
                for (int iz = 0; iz < N; iz++) {
                    double z = (iz + 0.5) * epsilon;
                    double dz = z - ic_cz;
                    if (dz > 0.5) dz -= 1.0;
                    if (dz < -0.5) dz += 1.0;
                    double r2 = dx * dx + dy * dy + dz * dz;
                    double amp = exp(-r2 / (4.0 * w0 * w0));
                    double phase = px0 * x + py0 * y + pz0 * z;
                    psi[idx(ix, iy, iz)] = amp * Complex(cos(phase), sin(phase));
                }
            }
        }
    }

    // Normalise: integral |psi|^2 dV = 1
    double dV = epsilon * epsilon * epsilon;
    double norm_sq = 0.0;
    for (int i = 0; i < N3; i++) norm_sq += std::norm(psi[i]);
    norm_sq *= dV;
    double norm_factor = sqrt(norm_sq);
    for (int i = 0; i < N3; i++) psi[i] /= norm_factor;

    double mass0 = 0.0;
    for (int i = 0; i < N3; i++) mass0 += std::norm(psi[i]);
    mass0 *= dV;
    std::cout << "Initial mass = " << mass0 << "\n\n";

    // --- FFTW setup ---
    std::vector<double> V(N3, 0.0);
    fftw_complex* work_in = nullptr;
    fftw_complex* work_out = nullptr;
    fftw_plan plan_fwd, plan_bwd;

    if (gravity) {
        work_in = fftw_alloc_complex(N3);
        work_out = fftw_alloc_complex(N3);
        plan_fwd = fftw_plan_dft_3d(N, N, N, work_in, work_out,
                                     FFTW_FORWARD, FFTW_ESTIMATE);
        plan_bwd = fftw_plan_dft_3d(N, N, N, work_out, work_in,
                                     FFTW_BACKWARD, FFTW_ESTIMATE);
        poisson_solve_3d(psi, V, work_in, work_out, plan_fwd, plan_bwd);
    }

    // --- Output setup ---
    std::string outdir = "output/sp_3d";
    if (tophat) outdir += "_tophat";
    if (gravity) outdir += "_grav";
    system(("mkdir -p " + outdir).c_str());
    write_slice(outdir + "/slice_000000.dat", psi, V, 0, 0.0);
    if (tophat) {
        write_radial_profile(outdir + "/radial_000000.dat", psi, 0, 0.0);
    }

    printf("%-8s  %-8s  %-18s  %-12s  %-8s\n",
           "step", "time", "mass", "delta_mass", "sec/step");
    printf("--------  --------  ------------------  ------------  --------\n");
    printf("%-8d  %-8.4f  %-18.15f  %+.2e  -\n", 0, 0.0, mass0, 0.0);

    auto t_start = std::chrono::high_resolution_clock::now();

    // --- Time loop: Strang splitting K/2 - V - K/2 ---
    for (int step = 1; step <= N_STEPS; step++) {
        kinetic_half_step_3d(psi);             // K/2

        if (gravity) {
            poisson_solve_3d(psi, V, work_in, work_out, plan_fwd, plan_bwd);
            potential_step_3d(psi, V);          // V
        }

        kinetic_half_step_3d(psi);             // K/2

        if (step % OUTPUT_EVERY == 0) {
            auto t_now = std::chrono::high_resolution_clock::now();
            double elapsed = std::chrono::duration<double>(t_now - t_start).count();
            double sec_per_step = elapsed / step;

            double mass = 0.0;
            for (int i = 0; i < N3; i++) mass += std::norm(psi[i]);
            mass *= dV;

            double time = step * delta;
            printf("%-8d  %-8.4f  %-18.15f  %+.2e  %.4f\n",
                   step, time, mass, mass - mass0, sec_per_step);

            char fname[64];
            snprintf(fname, sizeof(fname), "/slice_%06d.dat", step);
            write_slice(outdir + fname, psi, V, step, time);
            if (tophat) {
                snprintf(fname, sizeof(fname), "/radial_%06d.dat", step);
                write_radial_profile(outdir + fname, psi, step, time);
            }
        }
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double total = std::chrono::duration<double>(t_end - t_start).count();
    printf("\nTotal wall time: %.1f s (%.4f s/step)\n", total, total / N_STEPS);

    // --- Cleanup ---
    if (gravity) {
        fftw_destroy_plan(plan_fwd);
        fftw_destroy_plan(plan_bwd);
        fftw_free(work_in);
        fftw_free(work_out);
    }

    std::cout << "Output in " << outdir << "/\n";
    return 0;
}
