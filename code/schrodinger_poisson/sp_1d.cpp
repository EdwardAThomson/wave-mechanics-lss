// 1D Schrodinger-Poisson solver: Goldberg/Cayley method with periodic BCs
//
// Tests:
//   ./sp_1d           Free particle: Gaussian disperses through periodic BCs
//   ./sp_1d 5000      Gravity on: Gaussian collapses under self-gravity
//
// Validates: kinetic step (Goldberg recursion + double-pass periodic BCs),
//            Poisson solve (FFT with discrete Green's function),
//            potential step (diagonal Cayley), mass conservation.
//
// Reference: Thomson (2011) Chapter 5
//
// Equation: i dpsi/dt = -d^2psi/dx^2 + G_eff * V * psi   (hbar=1, m=1/2)
// Gravity:  d^2V/dx^2 = |psi|^2 - <|psi|^2>
// Box: [0,1) with periodic BCs.

#include <cmath>
#include <complex>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <fftw3.h>

using Complex = std::complex<double>;

// === Grid ===
constexpr int N = 512;
constexpr double BOX = 1.0;
constexpr double epsilon = BOX / N;        // grid spacing
constexpr double PI = M_PI;

// === Timestep ===
constexpr double delta = 1e-5;             // timestep
constexpr int N_STEPS = 20000;             // total steps (total time = 0.2)
constexpr int OUTPUT_EVERY = 2000;         // output interval

// Kinetic-only lambda for split-operator form (thesis Eq 5.38):
//   lambda = 4 * epsilon^2 / delta
// The factor of 4 (vs 2 in the combined K+V form) accounts for the
// half-step: each kinetic_half_step evolves for delta/2.
const double LAMBDA = 4.0 * epsilon * epsilon / delta;

// === Initial condition: Gaussian wavepacket ===
//   psi(x) = A * exp(-(x-x0)^2 / (4*w^2)) * exp(i*p0*x)
constexpr double x0 = 0.5;                // centre of Gaussian
constexpr double w0 = 0.05;               // width (std dev of |psi|^2)
constexpr double p0 = 0.0;                // initial momentum (carrier wave)

// === Gravitational coupling (set from command line, 0 = free particle) ===
double G_EFF = 0.0;

// ================================================================
// Periodic index
// ================================================================
inline int wrap(int j) {
    return ((j % N) + N) % N;
}

// ================================================================
// Kinetic half-step: Goldberg/Cayley with double-pass periodic BCs
//
// Solves: (1 + iK*delta/4) psi^{n+1} = (1 - iK*delta/4) psi^n
// where K = -d^2/dx^2 (finite difference).
//
// The auxiliary functions e, f, Omega are computed with double
// recursion to enforce periodicity. This is the thesis author's
// original contribution (thesis section 5.2.1).
//
// With splitting operators (K only, no V), the recursion simplifies:
//   e[j] = (2 - i*lambda) - 1/e[j-1]
// because V does not appear in the kinetic step.
// ================================================================
void kinetic_half_step(std::vector<Complex>& psi) {
    const Complex i_lam(0.0, LAMBDA);
    const Complex alpha = 2.0 - i_lam;    // coefficient for e recursion
    const Complex beta  = 2.0 + i_lam;    // coefficient for Omega (RHS)

    std::vector<Complex> Omega(N), e(N), f(N);

    // --- Omega: no recursion needed (thesis Eq 5.15, kinetic-only) ---
    // Omega_j = -psi_{j+1} + (i*lambda + 2)*psi_j - psi_{j-1}
    for (int j = 0; j < N; j++) {
        Omega[j] = -psi[wrap(j + 1)] + beta * psi[j] - psi[wrap(j - 1)];
    }

    // --- e: first forward pass (thesis Eqs 5.20-5.22) ---
    e[0] = alpha;
    for (int j = 1; j < N; j++) {
        e[j] = alpha - 1.0 / e[j - 1];
    }

    // --- e: second forward pass, periodic (thesis Eqs 5.23-5.25) ---
    // Seed with e[N-1] from first pass to wrap across boundary
    e[0] = alpha - 1.0 / e[N - 1];
    for (int j = 1; j < N; j++) {
        e[j] = alpha - 1.0 / e[j - 1];
    }

    // --- f: first forward pass (thesis Eqs 5.26-5.27) ---
    f[0] = Omega[0];
    for (int j = 1; j < N; j++) {
        f[j] = Omega[j] + f[j - 1] / e[j - 1];
    }

    // --- f: second forward pass, periodic (thesis Eqs 5.28-5.29) ---
    f[0] = Omega[0] + f[N - 1] / e[N - 1];
    for (int j = 1; j < N; j++) {
        f[j] = Omega[j] + f[j - 1] / e[j - 1];
    }

    // --- psi: first backward pass (thesis Eq 5.30) ---
    // psi_{j} = (psi_{j+1} - f_j) / e_j, from j=N-1 down to 0
    // Seed: psi[N-1] = (psi[0] - f[N-1]) / e[N-1]  (periodicity)
    // On first pass, psi[0] is still the old (time n) value.
    psi[N - 1] = (psi[0] - f[N - 1]) / e[N - 1];
    for (int j = N - 2; j >= 0; j--) {
        psi[j] = (psi[j + 1] - f[j]) / e[j];
    }

    // --- psi: second backward pass (periodic convergence) ---
    // Now psi[0] has been updated by the first pass, giving a better seed.
    psi[N - 1] = (psi[0] - f[N - 1]) / e[N - 1];
    for (int j = N - 2; j >= 0; j--) {
        psi[j] = (psi[j + 1] - f[j]) / e[j];
    }
}

// ================================================================
// Poisson solve: d^2V/dx^2 = |psi|^2 - <|psi|^2>
//
// Uses FFT with the discrete Green's function (thesis section 5.5).
// The discrete Laplacian eigenvalue for mode k is:
//   (2*cos(2*pi*k/N) - 2) / epsilon^2
// So V_hat[k] = epsilon^2 * rho_hat[k] / (2*cos(2*pi*k/N) - 2)
// with V_hat[0] = 0 (mean potential is arbitrary).
// ================================================================
void poisson_solve(const std::vector<Complex>& psi,
                   std::vector<double>& V,
                   fftw_complex* work_in,
                   fftw_complex* work_out,
                   fftw_plan& plan_fwd,
                   fftw_plan& plan_bwd) {
    // Compute mean-subtracted density
    double rho_mean = 0.0;
    for (int j = 0; j < N; j++) {
        rho_mean += std::norm(psi[j]);
    }
    rho_mean /= N;

    for (int j = 0; j < N; j++) {
        work_in[j][0] = std::norm(psi[j]) - rho_mean;
        work_in[j][1] = 0.0;
    }

    // Forward FFT
    fftw_execute(plan_fwd);

    // Apply discrete Green's function in k-space
    work_out[0][0] = 0.0;
    work_out[0][1] = 0.0;
    for (int k = 1; k < N; k++) {
        double eigenvalue = 2.0 * cos(2.0 * PI * k / N) - 2.0;
        double factor = epsilon * epsilon / eigenvalue;
        // work_out already contains rho_hat from the forward FFT
        work_out[k][0] *= factor;
        work_out[k][1] *= factor;
    }

    // Backward FFT -> V(x)
    // Swap pointers: use work_out as input to backward transform
    fftw_execute(plan_bwd);

    // FFTW normalisation: divide by N
    for (int j = 0; j < N; j++) {
        V[j] = work_in[j][0] / N;
    }
}

// ================================================================
// Potential step: diagonal Cayley transform (thesis Eq 5.36)
//
// psi *= (1 - i*G_eff*V*delta/2) / (1 + i*G_eff*V*delta/2)
//
// This is a pure phase rotation (|factor| = 1), so exactly unitary.
// The full timestep delta is used here (not delta/2), because the
// Strang splitting assigns the full step to V: exp(-iVdt).
// The delta/2 in the formula is from the Cayley approximation of
// the exponential, not from the splitting.
// ================================================================
void potential_step(std::vector<Complex>& psi,
                    const std::vector<double>& V) {
    for (int j = 0; j < N; j++) {
        double vdt2 = G_EFF * V[j] * delta / 2.0;
        Complex numerator(1.0, -vdt2);
        Complex denominator(1.0, vdt2);
        psi[j] *= numerator / denominator;
    }
}

// ================================================================
// Output: density, phase, real/imag parts, potential
// ================================================================
void write_output(const std::string& filename,
                  const std::vector<Complex>& psi,
                  const std::vector<double>& V,
                  int step, double time) {
    std::ofstream file(filename);
    file << "# step = " << step << ", time = " << time
         << ", G_eff = " << G_EFF << "\n";
    file << "# x  rho  phase  re  im  V\n";
    for (int j = 0; j < N; j++) {
        double x = (j + 0.5) * epsilon;
        double rho = std::norm(psi[j]);
        double phase = std::arg(psi[j]);
        file << x << "  " << rho << "  " << phase << "  "
             << psi[j].real() << "  " << psi[j].imag() << "  "
             << V[j] << "\n";
    }
    file.close();
}

// ================================================================
// Main
// ================================================================
int main(int argc, char* argv[]) {
    // Parse command line: ./sp_1d [G_eff]
    if (argc >= 2) {
        G_EFF = std::stod(argv[1]);
    }
    bool gravity = (G_EFF != 0.0);

    std::cout << "1D Schrodinger-Poisson solver (Goldberg/Cayley + periodic BCs)\n";
    std::cout << "N = " << N << ", epsilon = " << epsilon << "\n";
    std::cout << "delta = " << delta << ", lambda = " << LAMBDA << "\n";
    std::cout << "Total time = " << N_STEPS * delta << "\n";
    std::cout << "G_eff = " << G_EFF
              << (gravity ? " (gravity ON)" : " (free particle)") << "\n";
    std::cout << "IC: Gaussian x0=" << x0 << " w=" << w0 << " p0=" << p0 << "\n";

    if (!gravity) {
        double t_final = N_STEPS * delta;
        double sigma_final = sqrt(w0 * w0 + 4.0 * t_final * t_final / (w0 * w0));
        std::cout << "Expected sigma at t=" << t_final << ": " << sigma_final
                  << " (box = " << BOX << ")\n";
    }
    std::cout << "\n";

    // --- Initialise wavefunction ---
    std::vector<Complex> psi(N);
    for (int j = 0; j < N; j++) {
        double x = (j + 0.5) * epsilon;
        double dx = x - x0;
        // Nearest image for periodic box
        if (dx > 0.5 * BOX) dx -= BOX;
        if (dx < -0.5 * BOX) dx += BOX;
        double amp = exp(-dx * dx / (4.0 * w0 * w0));
        psi[j] = amp * Complex(cos(p0 * x), sin(p0 * x));
    }

    // Normalise: integral |psi|^2 dx = 1
    double norm_sq = 0.0;
    for (int j = 0; j < N; j++) norm_sq += std::norm(psi[j]);
    norm_sq *= epsilon;
    double norm_factor = sqrt(norm_sq);
    for (int j = 0; j < N; j++) psi[j] /= norm_factor;

    double mass0 = 0.0;
    for (int j = 0; j < N; j++) mass0 += std::norm(psi[j]);
    mass0 *= epsilon;

    // --- FFTW setup for Poisson solve ---
    std::vector<double> V(N, 0.0);
    fftw_complex* work_in = nullptr;
    fftw_complex* work_out = nullptr;
    fftw_plan plan_fwd, plan_bwd;

    if (gravity) {
        work_in = fftw_alloc_complex(N);
        work_out = fftw_alloc_complex(N);
        plan_fwd = fftw_plan_dft_1d(N, work_in, work_out,
                                     FFTW_FORWARD, FFTW_ESTIMATE);
        plan_bwd = fftw_plan_dft_1d(N, work_out, work_in,
                                     FFTW_BACKWARD, FFTW_ESTIMATE);

        // Compute initial potential
        poisson_solve(psi, V, work_in, work_out, plan_fwd, plan_bwd);
    }

    // --- Output setup ---
    std::string outdir = gravity ? "output/sp_1d_grav" : "output/sp_1d";
    system(("mkdir -p " + outdir).c_str());
    write_output(outdir + "/psi_000000.dat", psi, V, 0, 0.0);

    printf("%-8s  %-8s  %-18s  %-12s\n", "step", "time", "mass", "delta_mass");
    printf("--------  --------  ------------------  ------------\n");
    printf("%-8d  %-8.5f  %-18.15f  %+.2e\n", 0, 0.0, mass0, 0.0);

    // --- Time loop: Strang splitting K/2 - V - K/2 ---
    for (int step = 1; step <= N_STEPS; step++) {
        kinetic_half_step(psi);                // K/2

        if (gravity) {
            poisson_solve(psi, V, work_in, work_out, plan_fwd, plan_bwd);
            potential_step(psi, V);             // V
        }

        kinetic_half_step(psi);                // K/2

        if (step % OUTPUT_EVERY == 0) {
            double mass = 0.0;
            for (int j = 0; j < N; j++) mass += std::norm(psi[j]);
            mass *= epsilon;

            double time = step * delta;
            printf("%-8d  %-8.5f  %-18.15f  %+.2e\n",
                   step, time, mass, mass - mass0);

            char fname[64];
            snprintf(fname, sizeof(fname), "/psi_%06d.dat", step);
            write_output(outdir + fname, psi, V, step, time);
        }
    }

    // --- Cleanup ---
    if (gravity) {
        fftw_destroy_plan(plan_fwd);
        fftw_destroy_plan(plan_bwd);
        fftw_free(work_in);
        fftw_free(work_out);
    }

    std::cout << "\nDone. Output written to " << outdir << "/\n";
    return 0;
}
