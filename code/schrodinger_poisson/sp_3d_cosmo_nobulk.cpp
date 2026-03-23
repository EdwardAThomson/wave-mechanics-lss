// 3D Schrodinger-Poisson solver with cosmological initial conditions
//
// Combines the BBKS IC generator (from Chapter 4 FPA code) with the
// full S-P solver (Chapter 5) including LCDM expansion.
//
// Usage: ./sp_3d_cosmo [nu]
//   nu = effective quantum parameter (default 1e-4)
//   Smaller nu -> more classical, but phase must be grid-resolvable
//
// The IC pipeline:
//   1. Generate BBKS Gaussian random density field at z_init
//   2. Compute velocity potential from Poisson equation
//   3. Construct Madelung wavefunction: psi = sqrt(1+delta) * exp(-i*phi_v/nu)
//   4. Evolve with expanding S-P solver
//
// Reference: Thomson (2011) Chapters 4-5

#include <cmath>
#include <complex>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <random>
#include <algorithm>
#include <chrono>
#include <fftw3.h>

using Complex = std::complex<double>;

// === Grid ===
constexpr int N = 64;
constexpr int N3 = N * N * N;
constexpr double BOX = 1.0;
constexpr double epsilon = BOX / N;
constexpr double PI = M_PI;

// === Cosmological parameters ===
constexpr double h_cosmo = 0.71;
constexpr double Omega_m0 = 0.27;
constexpr double Omega_Lambda = 0.73;
constexpr double sigma_8_target = 0.81;
constexpr double n_s = 1.0;
constexpr double L_phys = 32.0;         // Mpc/h
constexpr double Gamma_bbks = Omega_m0 * h_cosmo;
constexpr double z_init = 50.0;
constexpr double a_init = 1.0 / (1.0 + z_init);
constexpr unsigned int SEED = 42;

// === Quantum / expansion parameters ===
double NU = 1e-4;                       // set from command line

// L0 = L_param(a=1): derived from nu so that the Madelung phase and
// the S-P kinetic term are consistent (see algorithm.md).
// At high z (matter-dominated): L0 ≈ 1/(2*nu*sqrt(a_init))
double L0;

// === Expansion functions (flat LCDM, thesis Eq 5.55) ===
double L_param(double a) { return L0 * sqrt(a); }
double E_func(double a) { return Omega_m0 + Omega_Lambda * a * a * a; }
double C_func(double a) { return 2.0 * L_param(a) * sqrt(E_func(a)); }
double lambda_of_a(double a, double d_lna) {
    return 4.0 * C_func(a) * epsilon * epsilon / d_lna;
}
double g_eff_of_a(double a) {
    return 3.0 * Omega_m0 * L_param(a) * L_param(a) / C_func(a);
}

// ================================================================
// 3D index
// ================================================================
inline int idx(int ix, int iy, int iz) {
    return (ix * N + iy) * N + iz;
}

inline int wrap(int j) {
    return ((j % N) + N) % N;
}

// Wavenumber for index along one axis (box units)
double wavenumber(int index) {
    int k = (index <= N / 2) ? index : index - N;
    return 2.0 * PI * k / BOX;
}

// ================================================================
// BBKS transfer function and power spectrum (from fpa_3d_cosmo)
// ================================================================
double transfer_function(double k_phys) {
    if (k_phys == 0.0) return 1.0;
    double q = k_phys / Gamma_bbks;
    double log_term = log(1.0 + 2.34 * q) / (2.34 * q);
    double poly = 1.0 + 3.89 * q + pow(16.1 * q, 2)
                  + pow(5.46 * q, 3) + pow(6.71 * q, 4);
    return log_term * pow(poly, -0.25);
}

double power_spectrum_unnorm(double k_phys) {
    if (k_phys == 0.0) return 0.0;
    double Tk = transfer_function(k_phys);
    return pow(k_phys, n_s) * Tk * Tk;
}

double tophat_window(double kR) {
    if (kR < 1e-6) return 1.0;
    return 3.0 * (sin(kR) - kR * cos(kR)) / (kR * kR * kR);
}

double compute_sigma8_sq_unnorm() {
    double R8 = 8.0;
    int n_int = 10000;
    double k_min = 1e-4, k_max = 1e2;
    double dlnk = (log(k_max) - log(k_min)) / n_int;
    double sigma2 = 0.0;
    for (int i = 0; i < n_int; i++) {
        double lnk = log(k_min) + (i + 0.5) * dlnk;
        double k = exp(lnk);
        double W = tophat_window(k * R8);
        sigma2 += k * k * power_spectrum_unnorm(k) * W * W * k * dlnk;
    }
    sigma2 /= (2.0 * PI * PI);
    return sigma2;
}

// ================================================================
// Growth factor (from fpa_3d_cosmo)
// ================================================================
double hubble(double a) {
    return sqrt(Omega_m0 / (a * a * a) + Omega_Lambda);
}

double growth_factor(double a) {
    int nsteps = 10000;
    double da = a / nsteps;
    double integral = 0.0;
    for (int i = 0; i < nsteps; i++) {
        double a0 = std::max((i + 0.0) * da, 1e-10);
        double a1 = (i + 0.5) * da;
        double a2 = (i + 1.0) * da;
        double f0 = 1.0 / pow(a0 * hubble(a0), 3);
        double f1 = 1.0 / pow(a1 * hubble(a1), 3);
        double f2 = 1.0 / pow(a2 * hubble(a2), 3);
        integral += (f0 + 4.0 * f1 + f2) * da / 6.0;
    }
    return 2.5 * Omega_m0 * hubble(a) * integral;
}

double compute_growth_ratio() {
    return growth_factor(1.0) / growth_factor(a_init);
}

// ================================================================
// Gaussian random field generation (from fpa_3d_cosmo)
// ================================================================
void generate_gaussian_field(fftw_complex* delta_hat, double D_ratio) {
    double sigma8_sq_raw = compute_sigma8_sq_unnorm();
    double A = (sigma_8_target * sigma_8_target) / sigma8_sq_raw;
    double A_init = A / (D_ratio * D_ratio);
    double L3 = L_phys * L_phys * L_phys;
    double N3d = static_cast<double>(N3);

    std::mt19937_64 rng(SEED);
    std::normal_distribution<double> gauss(0.0, 1.0);

    for (int i = 0; i < N3; i++) {
        delta_hat[i][0] = 0.0;
        delta_hat[i][1] = 0.0;
    }

    for (int ix = 0; ix < N; ix++) {
        double kx_box = wavenumber(ix);
        for (int iy = 0; iy < N; iy++) {
            double ky_box = wavenumber(iy);
            for (int iz = 0; iz <= N / 2; iz++) {
                double kz_box = wavenumber(iz);
                int i = idx(ix, iy, iz);
                double k2_box = kx_box * kx_box + ky_box * ky_box + kz_box * kz_box;
                if (k2_box == 0.0) continue;

                double k_phys = sqrt(k2_box) / L_phys;
                double Pk_phys = A_init * power_spectrum_unnorm(k_phys);
                double Pk_box = Pk_phys / L3;
                double amplitude = N3d * sqrt(Pk_box);

                int jx = (N - ix) % N;
                int jy = (N - iy) % N;
                int jz = (N - iz) % N;

                if (iz > 0 && iz < N / 2) {
                    double g1 = gauss(rng);
                    double g2 = gauss(rng);
                    delta_hat[i][0] = amplitude * g1 / sqrt(2.0);
                    delta_hat[i][1] = amplitude * g2 / sqrt(2.0);
                    int j = idx(jx, jy, jz);
                    delta_hat[j][0] =  delta_hat[i][0];
                    delta_hat[j][1] = -delta_hat[i][1];
                } else {
                    bool self_conjugate = (ix == jx && iy == jy);
                    if (self_conjugate) {
                        delta_hat[i][0] = amplitude * gauss(rng);
                        delta_hat[i][1] = 0.0;
                    } else if (ix * N + iy < jx * N + jy) {
                        double g1 = gauss(rng);
                        double g2 = gauss(rng);
                        delta_hat[i][0] = amplitude * g1 / sqrt(2.0);
                        delta_hat[i][1] = amplitude * g2 / sqrt(2.0);
                        int j = idx(jx, jy, iz);
                        delta_hat[j][0] =  delta_hat[i][0];
                        delta_hat[j][1] = -delta_hat[i][1];
                    }
                }
            }
        }
    }
}

// ================================================================
// Velocity potential: phi_v_hat = -delta_hat / k^2
// ================================================================
void compute_velocity_potential(const fftw_complex* delta_hat,
                                std::vector<double>& phi_v) {
    fftw_complex* phi_hat = fftw_alloc_complex(N3);
    fftw_complex* phi_real = fftw_alloc_complex(N3);

    for (int ix = 0; ix < N; ix++) {
        double kx = wavenumber(ix);
        for (int iy = 0; iy < N; iy++) {
            double ky = wavenumber(iy);
            for (int iz = 0; iz < N; iz++) {
                double kz = wavenumber(iz);
                int i = idx(ix, iy, iz);
                double k2 = kx * kx + ky * ky + kz * kz;
                if (k2 == 0.0) {
                    phi_hat[i][0] = 0.0;
                    phi_hat[i][1] = 0.0;
                } else {
                    phi_hat[i][0] = -delta_hat[i][0] / k2;
                    phi_hat[i][1] = -delta_hat[i][1] / k2;
                }
            }
        }
    }

    fftw_plan ifft = fftw_plan_dft_3d(N, N, N, phi_hat, phi_real,
                                       FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(ifft);
    fftw_destroy_plan(ifft);

    for (int i = 0; i < N3; i++) {
        phi_v[i] = phi_real[i][0] / static_cast<double>(N3);
    }

    fftw_free(phi_hat);
    fftw_free(phi_real);
}

// ================================================================
// Displacement field (for phase gradient check)
// ================================================================
void compute_displacement(const fftw_complex* delta_hat,
                          std::vector<double>& sx,
                          std::vector<double>& sy,
                          std::vector<double>& sz) {
    fftw_complex* s_hat = fftw_alloc_complex(N3);
    fftw_complex* s_real = fftw_alloc_complex(N3);
    for (int comp = 0; comp < 3; comp++) {
        for (int ix = 0; ix < N; ix++) {
            double kx = wavenumber(ix);
            for (int iy = 0; iy < N; iy++) {
                double ky = wavenumber(iy);
                for (int iz = 0; iz < N; iz++) {
                    double kz = wavenumber(iz);
                    int i = idx(ix, iy, iz);
                    double k2 = kx * kx + ky * ky + kz * kz;
                    if (k2 == 0.0) { s_hat[i][0] = 0.0; s_hat[i][1] = 0.0; continue; }
                    double kj = (comp == 0) ? kx : (comp == 1) ? ky : kz;
                    s_hat[i][0] = -kj * delta_hat[i][1] / k2;
                    s_hat[i][1] =  kj * delta_hat[i][0] / k2;
                }
            }
        }
        fftw_plan ifft = fftw_plan_dft_3d(N, N, N, s_hat, s_real,
                                           FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(ifft);
        fftw_destroy_plan(ifft);
        std::vector<double>& s = (comp == 0) ? sx : (comp == 1) ? sy : sz;
        for (int i = 0; i < N3; i++) s[i] = s_real[i][0] / static_cast<double>(N3);
    }
    fftw_free(s_hat);
    fftw_free(s_real);
}

// ================================================================
// Goldberg/Cayley 1D sweep (from sp_3d.cpp)
// ================================================================
void goldberg_sweep(Complex* pencil, double lambda) {
    const Complex i_lam(0.0, lambda);
    const Complex alpha = 2.0 - i_lam;
    const Complex beta  = 2.0 + i_lam;

    Complex Omega[N], e[N], f[N];

    for (int j = 0; j < N; j++) {
        Omega[j] = -pencil[wrap(j + 1)] + beta * pencil[j] - pencil[wrap(j - 1)];
    }

    e[0] = alpha;
    for (int j = 1; j < N; j++) e[j] = alpha - 1.0 / e[j - 1];
    e[0] = alpha - 1.0 / e[N - 1];
    for (int j = 1; j < N; j++) e[j] = alpha - 1.0 / e[j - 1];

    f[0] = Omega[0];
    for (int j = 1; j < N; j++) f[j] = Omega[j] + f[j - 1] / e[j - 1];
    f[0] = Omega[0] + f[N - 1] / e[N - 1];
    for (int j = 1; j < N; j++) f[j] = Omega[j] + f[j - 1] / e[j - 1];

    pencil[N - 1] = (pencil[0] - f[N - 1]) / e[N - 1];
    for (int j = N - 2; j >= 0; j--) pencil[j] = (pencil[j + 1] - f[j]) / e[j];
    pencil[N - 1] = (pencil[0] - f[N - 1]) / e[N - 1];
    for (int j = N - 2; j >= 0; j--) pencil[j] = (pencil[j + 1] - f[j]) / e[j];
}

// ================================================================
// 3D kinetic half-step: sweep along x, y, z
// ================================================================
void kinetic_half_step_3d(std::vector<Complex>& psi, double lambda) {
    Complex pencil[N];

    for (int iy = 0; iy < N; iy++)
        for (int iz = 0; iz < N; iz++) {
            for (int ix = 0; ix < N; ix++) pencil[ix] = psi[idx(ix, iy, iz)];
            goldberg_sweep(pencil, lambda);
            for (int ix = 0; ix < N; ix++) psi[idx(ix, iy, iz)] = pencil[ix];
        }

    for (int ix = 0; ix < N; ix++)
        for (int iz = 0; iz < N; iz++) {
            for (int iy = 0; iy < N; iy++) pencil[iy] = psi[idx(ix, iy, iz)];
            goldberg_sweep(pencil, lambda);
            for (int iy = 0; iy < N; iy++) psi[idx(ix, iy, iz)] = pencil[iy];
        }

    for (int ix = 0; ix < N; ix++)
        for (int iy = 0; iy < N; iy++) {
            int base = idx(ix, iy, 0);
            for (int iz = 0; iz < N; iz++) pencil[iz] = psi[base + iz];
            goldberg_sweep(pencil, lambda);
            for (int iz = 0; iz < N; iz++) psi[base + iz] = pencil[iz];
        }
}

// ================================================================
// Bulk flow removal: Galilean boost to zero mean current
//
// 1. Compute mean probability current <j_i> over the whole grid
// 2. Derive bulk velocity: v_bulk = <j> / <rho>
// 3. Apply phase correction: psi(x) *= exp(-i * v_bulk . x / nu)
//
// This is a unitary operation (pure phase shift) — mass is preserved
// exactly. Analogous to subtracting mean velocity in N-body codes.
// ================================================================
void remove_bulk_flow(std::vector<Complex>& psi, double nu) {
    double inv_2eps = 1.0 / (2.0 * epsilon);
    double mean_jx = 0.0, mean_jy = 0.0, mean_jz = 0.0;
    double mean_rho = 0.0;

    for (int ix = 0; ix < N; ix++) {
        int ixp = (ix + 1) % N;
        int ixm = (ix - 1 + N) % N;
        for (int iy = 0; iy < N; iy++) {
            int iyp = (iy + 1) % N;
            int iym = (iy - 1 + N) % N;
            for (int iz = 0; iz < N; iz++) {
                int izp = (iz + 1) % N;
                int izm = (iz - 1 + N) % N;
                int i = idx(ix, iy, iz);
                Complex psi_star = std::conj(psi[i]);
                double rho = std::norm(psi[i]);

                Complex dpsi_dx = (psi[idx(ixp, iy, iz)] - psi[idx(ixm, iy, iz)]) * inv_2eps;
                Complex dpsi_dy = (psi[idx(ix, iyp, iz)] - psi[idx(ix, iym, iz)]) * inv_2eps;
                Complex dpsi_dz = (psi[idx(ix, iy, izp)] - psi[idx(ix, iy, izm)]) * inv_2eps;

                mean_jx += (psi_star * dpsi_dx).imag();
                mean_jy += (psi_star * dpsi_dy).imag();
                mean_jz += (psi_star * dpsi_dz).imag();
                mean_rho += rho;
            }
        }
    }

    // Bulk velocity = nu * <j> / <rho>  (the nu from j = nu*Im(psi* grad psi))
    double vx_bulk = nu * mean_jx / mean_rho;
    double vy_bulk = nu * mean_jy / mean_rho;
    double vz_bulk = nu * mean_jz / mean_rho;

    // Apply Galilean boost: psi *= exp(-i * v_bulk . x / nu)
    for (int ix = 0; ix < N; ix++) {
        double x = (ix + 0.5) * epsilon;
        for (int iy = 0; iy < N; iy++) {
            double y = (iy + 0.5) * epsilon;
            for (int iz = 0; iz < N; iz++) {
                double z = (iz + 0.5) * epsilon;
                double phase = -(vx_bulk * x + vy_bulk * y + vz_bulk * z) / nu;
                psi[idx(ix, iy, iz)] *= Complex(cos(phase), sin(phase));
            }
        }
    }
}

// ================================================================
// 3D Poisson solve
// ================================================================
void poisson_solve_3d(const std::vector<Complex>& psi,
                      std::vector<double>& V,
                      fftw_complex* work_in, fftw_complex* work_out,
                      fftw_plan& plan_fwd, fftw_plan& plan_bwd) {
    double rho_mean = 0.0;
    for (int i = 0; i < N3; i++) rho_mean += std::norm(psi[i]);
    rho_mean /= N3;

    for (int i = 0; i < N3; i++) {
        work_in[i][0] = std::norm(psi[i]) - rho_mean;
        work_in[i][1] = 0.0;
    }

    fftw_execute(plan_fwd);

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

    fftw_execute(plan_bwd);

    for (int i = 0; i < N3; i++) {
        V[i] = work_in[i][0] / static_cast<double>(N3);
    }
}

// ================================================================
// Potential step
// ================================================================
void potential_step_3d(std::vector<Complex>& psi,
                       const std::vector<double>& V,
                       double g_eff, double dt) {
    for (int i = 0; i < N3; i++) {
        double vdt2 = g_eff * V[i] * dt / 2.0;
        Complex num(1.0, -vdt2);
        Complex den(1.0,  vdt2);
        psi[i] *= num / den;
    }
}

// ================================================================
// Velocity computation: two methods
//
// Method 1: Probability current (robust, well-defined everywhere)
//   j = nu * Im(psi* grad(psi))
//   v = j / |psi|^2
//
// Method 2: Phase gradient (classical Madelung velocity)
//   v = nu * grad(arg(psi))
//   Breaks down at nodal points where |psi| = 0
//
// Both use centered finite differences for the gradient.
// ================================================================

// Periodic centered difference: (f[j+1] - f[j-1]) / (2*epsilon)
// applied along one dimension of the 3D grid.

void compute_velocity_current(const std::vector<Complex>& psi, double nu,
                              std::vector<double>& vx,
                              std::vector<double>& vy,
                              std::vector<double>& vz) {
    double inv_2eps = 1.0 / (2.0 * epsilon);
    for (int ix = 0; ix < N; ix++) {
        int ixp = (ix + 1) % N;
        int ixm = (ix - 1 + N) % N;
        for (int iy = 0; iy < N; iy++) {
            int iyp = (iy + 1) % N;
            int iym = (iy - 1 + N) % N;
            for (int iz = 0; iz < N; iz++) {
                int izp = (iz + 1) % N;
                int izm = (iz - 1 + N) % N;
                int i = idx(ix, iy, iz);
                Complex psi_star = std::conj(psi[i]);
                double rho = std::norm(psi[i]);

                // grad(psi) via centered differences
                Complex dpsi_dx = (psi[idx(ixp, iy, iz)] - psi[idx(ixm, iy, iz)]) * inv_2eps;
                Complex dpsi_dy = (psi[idx(ix, iyp, iz)] - psi[idx(ix, iym, iz)]) * inv_2eps;
                Complex dpsi_dz = (psi[idx(ix, iy, izp)] - psi[idx(ix, iy, izm)]) * inv_2eps;

                // j = nu * Im(psi* grad(psi)), v = j / rho
                if (rho > 1e-30) {
                    vx[i] = nu * (psi_star * dpsi_dx).imag() / rho;
                    vy[i] = nu * (psi_star * dpsi_dy).imag() / rho;
                    vz[i] = nu * (psi_star * dpsi_dz).imag() / rho;
                } else {
                    vx[i] = 0.0;
                    vy[i] = 0.0;
                    vz[i] = 0.0;
                }
            }
        }
    }
}

void compute_velocity_phase(const std::vector<Complex>& psi, double nu,
                            std::vector<double>& vx,
                            std::vector<double>& vy,
                            std::vector<double>& vz) {
    double inv_2eps = 1.0 / (2.0 * epsilon);
    for (int ix = 0; ix < N; ix++) {
        int ixp = (ix + 1) % N;
        int ixm = (ix - 1 + N) % N;
        for (int iy = 0; iy < N; iy++) {
            int iyp = (iy + 1) % N;
            int iym = (iy - 1 + N) % N;
            for (int iz = 0; iz < N; iz++) {
                int izp = (iz + 1) % N;
                int izm = (iz - 1 + N) % N;
                int i = idx(ix, iy, iz);

                // Phase differences (automatic unwrapping via arg of ratio)
                // arg(psi[j+1]/psi[j-1]) gives the wrapped phase difference
                double rho = std::norm(psi[i]);
                if (rho > 1e-30) {
                    double dphase_dx = std::arg(psi[idx(ixp, iy, iz)] / psi[idx(ixm, iy, iz)]) * inv_2eps;
                    double dphase_dy = std::arg(psi[idx(ix, iyp, iz)] / psi[idx(ix, iym, iz)]) * inv_2eps;
                    double dphase_dz = std::arg(psi[idx(ix, iy, izp)] / psi[idx(ix, iy, izm)]) * inv_2eps;
                    vx[i] = nu * dphase_dx;
                    vy[i] = nu * dphase_dy;
                    vz[i] = nu * dphase_dz;
                } else {
                    vx[i] = 0.0;
                    vy[i] = 0.0;
                    vz[i] = 0.0;
                }
            }
        }
    }
}

// Write velocity comparison: z-slab averaged slice with both methods
void write_velocity_slice(const std::string& filename,
                          const std::vector<Complex>& psi,
                          const std::vector<double>& vx_curr,
                          const std::vector<double>& vy_curr,
                          const std::vector<double>& vx_phase,
                          const std::vector<double>& vy_phase,
                          int step, double a) {
    std::ofstream file(filename);
    file << "# step = " << step << ", a = " << a
         << ", z = " << 1.0 / a - 1.0 << "\n";
    file << "# x  y  rho  vx_current  vy_current  vx_phase  vy_phase\n";
    int iz = N / 2;
    for (int ix = 0; ix < N; ix++) {
        double x = (ix + 0.5) * epsilon;
        for (int iy = 0; iy < N; iy++) {
            double y = (iy + 0.5) * epsilon;
            int i = idx(ix, iy, iz);
            file << x << "  " << y << "  " << std::norm(psi[i]) << "  "
                 << vx_curr[i] << "  " << vy_curr[i] << "  "
                 << vx_phase[i] << "  " << vy_phase[i] << "\n";
        }
        file << "\n";
    }
    file.close();
}

// ================================================================
// Output: z-slice at iz = N/2
// ================================================================
constexpr int SLAB_THICK = 8;

void write_slice(const std::string& filename,
                 const std::vector<Complex>& psi,
                 int step, double a) {
    std::ofstream file(filename);
    file << "# step = " << step << ", a = " << a
         << ", z = " << 1.0 / a - 1.0 << "\n";
    file << "# x  y  rho  log2rho\n";
    int iz_start = N / 2 - SLAB_THICK / 2;
    for (int ix = 0; ix < N; ix++) {
        double x = (ix + 0.5) * epsilon;
        for (int iy = 0; iy < N; iy++) {
            double y = (iy + 0.5) * epsilon;
            double avg = 0.0;
            for (int diz = 0; diz < SLAB_THICK; diz++) {
                int iz = (iz_start + diz + N) % N;
                avg += std::norm(psi[idx(ix, iy, iz)]);
            }
            avg /= SLAB_THICK;
            file << x << "  " << y << "  " << avg << "  "
                 << log(2.0 + avg - 1.0) << "\n";
        }
        file << "\n";
    }
    file.close();
}

// ================================================================
// Binary wavefunction dump (for checkpointing and post-processing)
//
// Format: [header] [psi data]
//   header: step (int), N (int), a (double), nu (double)
//   data:   N^3 complex doubles (re,im pairs), row-major z-fastest
//
// Total size: 16 + N^3 * 16 bytes = ~4 MB at N=64
// ================================================================
void write_psi_binary(const std::string& filename,
                      const std::vector<Complex>& psi,
                      int step, double a) {
    std::ofstream file(filename, std::ios::binary);
    int n = N;
    file.write(reinterpret_cast<const char*>(&step), sizeof(int));
    file.write(reinterpret_cast<const char*>(&n), sizeof(int));
    file.write(reinterpret_cast<const char*>(&a), sizeof(double));
    file.write(reinterpret_cast<const char*>(&NU), sizeof(double));
    // Complex is guaranteed to be two contiguous doubles (re, im)
    file.write(reinterpret_cast<const char*>(psi.data()), N3 * sizeof(Complex));
    file.close();
}

// ================================================================
// Main
// ================================================================
int main(int argc, char* argv[]) {
    if (argc >= 2) NU = std::stod(argv[1]);

    // Derive L0 from nu: at high z, L0 ≈ 1/(2*nu*sqrt(a_init))
    L0 = 1.0 / (2.0 * NU * sqrt(a_init));

    // Timestep: ~400 steps from a_init to a=1
    double d_lna = log(1.0 / a_init) / 400;
    int N_STEPS = 400;
    int OUTPUT_EVERY = 40;

    // Output at scale factors matching thesis comparison points
    // a ≈ 0.33, 0.62, 0.78, 0.93 (thesis section 4.2.2)

    std::cout << "3D Schrodinger-Poisson with cosmological ICs\n";
    std::cout << "N = " << N << "^3, L_phys = " << L_phys << " Mpc/h\n";
    std::cout << "Cosmology: h=" << h_cosmo << " Omega_m=" << Omega_m0
              << " sigma_8=" << sigma_8_target << "\n";
    std::cout << "nu = " << NU << ", L0 = " << L0 << "\n";
    std::cout << "z_init = " << z_init << " (a=" << a_init << ")\n";
    std::cout << "d(ln a) = " << d_lna << ", N_steps = " << N_STEPS << "\n\n";

    // --- Generate cosmological ICs ---
    std::cout << "Computing growth factor...\n";
    double D_ratio = compute_growth_ratio();
    std::cout << "D(z=0)/D(z_init) = " << D_ratio << "\n";

    std::cout << "Generating Gaussian random field...\n";
    fftw_complex* delta_hat = fftw_alloc_complex(N3);
    generate_gaussian_field(delta_hat, D_ratio);

    // delta(x)
    fftw_complex* work_tmp = fftw_alloc_complex(N3);
    fftw_plan ifft_tmp = fftw_plan_dft_3d(N, N, N, delta_hat, work_tmp,
                                           FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(ifft_tmp);
    fftw_destroy_plan(ifft_tmp);
    std::vector<double> delta_init(N3);
    for (int i = 0; i < N3; i++) delta_init[i] = work_tmp[i][0] / static_cast<double>(N3);
    fftw_free(work_tmp);

    double delta_rms = 0.0;
    for (int i = 0; i < N3; i++) delta_rms += delta_init[i] * delta_init[i];
    delta_rms = sqrt(delta_rms / N3);
    std::cout << "Initial delta_rms = " << delta_rms << "\n";

    // Velocity potential
    std::cout << "Computing velocity potential...\n";
    std::vector<double> phi_v(N3);
    compute_velocity_potential(delta_hat, phi_v);

    // Phase gradient check
    std::vector<double> sx(N3), sy(N3), sz(N3);
    compute_displacement(delta_hat, sx, sy, sz);
    fftw_free(delta_hat);

    double s_max = 0.0;
    for (int i = 0; i < N3; i++) {
        double s = sqrt(sx[i] * sx[i] + sy[i] * sy[i] + sz[i] * sz[i]);
        if (s > s_max) s_max = s;
    }
    double phase_grad_max = s_max / NU;
    double k_nyq = PI * N;
    std::cout << "s_max = " << s_max << "\n";
    std::cout << "max|grad(phase)| / k_Nyq = " << phase_grad_max / k_nyq
              << (phase_grad_max < k_nyq ? " (RESOLVABLE)" : " (ALIASED!)") << "\n";

    // --- Construct Madelung wavefunction ---
    std::cout << "\nConstructing Madelung wavefunction (nu=" << NU << ")...\n";
    std::vector<Complex> psi(N3);
    for (int i = 0; i < N3; i++) {
        double amp = sqrt(std::max(1.0 + delta_init[i], 0.0));
        double phase = -phi_v[i] / NU;
        psi[i] = amp * Complex(cos(phase), sin(phase));
    }

    // Check psi_hat power distribution
    {
        fftw_complex* psi_k = fftw_alloc_complex(N3);
        fftw_complex* in = fftw_alloc_complex(N3);
        for (int i = 0; i < N3; i++) { in[i][0] = psi[i].real(); in[i][1] = psi[i].imag(); }
        fftw_plan fwd = fftw_plan_dft_3d(N, N, N, in, psi_k, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(fwd);
        fftw_destroy_plan(fwd);
        double power_low = 0.0, power_high = 0.0;
        double k_half = k_nyq / 2.0;
        for (int ix = 0; ix < N; ix++) {
            double kx = wavenumber(ix);
            for (int iy = 0; iy < N; iy++) {
                double ky = wavenumber(iy);
                for (int iz = 0; iz < N; iz++) {
                    double kz = wavenumber(iz);
                    int i = idx(ix, iy, iz);
                    double k = sqrt(kx*kx + ky*ky + kz*kz);
                    double pw = psi_k[i][0]*psi_k[i][0] + psi_k[i][1]*psi_k[i][1];
                    if (k < k_half) power_low += pw; else power_high += pw;
                }
            }
        }
        double pct_low = 100.0 * power_low / (power_low + power_high);
        std::cout << "psi_hat power: " << pct_low << "% in low-k (want >80%)\n";
        fftw_free(psi_k);
        fftw_free(in);
    }

    // Mass
    double dV = epsilon * epsilon * epsilon;
    double mass0 = 0.0;
    for (int i = 0; i < N3; i++) mass0 += std::norm(psi[i]);
    mass0 *= dV;
    std::cout << "Initial mass = " << mass0 << "\n\n";

    // --- FFTW setup for Poisson solve ---
    std::vector<double> V(N3, 0.0);
    fftw_complex* work_in = fftw_alloc_complex(N3);
    fftw_complex* work_out = fftw_alloc_complex(N3);
    fftw_plan plan_fwd = fftw_plan_dft_3d(N, N, N, work_in, work_out,
                                           FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan plan_bwd = fftw_plan_dft_3d(N, N, N, work_out, work_in,
                                           FFTW_BACKWARD, FFTW_ESTIMATE);

    // --- Velocity arrays ---
    std::vector<double> vx_curr(N3), vy_curr(N3), vz_curr(N3);
    std::vector<double> vx_phase(N3), vy_phase(N3), vz_phase(N3);

    // --- Output setup ---
    std::string outdir = "output/sp_cosmo_nobulk_nu_" + std::to_string(NU);
    system(("mkdir -p " + outdir).c_str());
    write_slice(outdir + "/slice_000000.dat", psi, 0, a_init);
    write_psi_binary(outdir + "/psi_000000.bin", psi, 0, a_init);

    // Initial velocity output
    compute_velocity_current(psi, NU, vx_curr, vy_curr, vz_curr);
    compute_velocity_phase(psi, NU, vx_phase, vy_phase, vz_phase);
    write_velocity_slice(outdir + "/vel_000000.dat", psi,
                         vx_curr, vy_curr, vx_phase, vy_phase, 0, a_init);

    printf("%-6s  %-7s  %-6s  %-18s  %-11s  %-8s  %-7s  %-6s\n",
           "step", "a", "z", "mass", "delta_mass", "g_eff", "lambda", "s/step");
    printf("------  -------  ------  ------------------  "
           "-----------  --------  -------  ------\n");

    double a = a_init;
    printf("%-6d  %-7.4f  %-6.2f  %-18.15f  %+.2e  %-8.1f  %-7.1f  -\n",
           0, a, 1.0/a - 1.0, mass0, 0.0, g_eff_of_a(a), lambda_of_a(a, d_lna));

    auto t_start = std::chrono::high_resolution_clock::now();

    // --- Time loop ---
    for (int step = 1; step <= N_STEPS; step++) {
        double lam = lambda_of_a(a, d_lna);
        double g_eff = g_eff_of_a(a);

        kinetic_half_step_3d(psi, lam);
        poisson_solve_3d(psi, V, work_in, work_out, plan_fwd, plan_bwd);
        potential_step_3d(psi, V, g_eff, d_lna);
        kinetic_half_step_3d(psi, lam);

        // Remove bulk flow (Galilean boost to zero mean current)
        remove_bulk_flow(psi, NU);

        a *= exp(d_lna);

        if (step % OUTPUT_EVERY == 0) {
            auto t_now = std::chrono::high_resolution_clock::now();
            double elapsed = std::chrono::duration<double>(t_now - t_start).count();

            double mass = 0.0;
            for (int i = 0; i < N3; i++) mass += std::norm(psi[i]);
            mass *= dV;

            printf("%-6d  %-7.4f  %-6.2f  %-18.15f  %+.2e  %-8.1f  %-7.1f  %.3f\n",
                   step, a, 1.0/a - 1.0, mass, mass - mass0,
                   g_eff_of_a(a), lambda_of_a(a, d_lna), elapsed / step);

            char fname[64];
            snprintf(fname, sizeof(fname), "/slice_%06d.dat", step);
            write_slice(outdir + fname, psi, step, a);

            snprintf(fname, sizeof(fname), "/psi_%06d.bin", step);
            write_psi_binary(outdir + fname, psi, step, a);

            // Velocity output
            compute_velocity_current(psi, NU, vx_curr, vy_curr, vz_curr);
            compute_velocity_phase(psi, NU, vx_phase, vy_phase, vz_phase);
            snprintf(fname, sizeof(fname), "/vel_%06d.dat", step);
            write_velocity_slice(outdir + fname, psi,
                                 vx_curr, vy_curr, vx_phase, vy_phase, step, a);
        }
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double total = std::chrono::duration<double>(t_end - t_start).count();
    printf("\nTotal: %.1f s (%.3f s/step)\n", total, total / N_STEPS);

    fftw_destroy_plan(plan_fwd);
    fftw_destroy_plan(plan_bwd);
    fftw_free(work_in);
    fftw_free(work_out);

    std::cout << "Output in " << outdir << "/\n";
    return 0;
}
