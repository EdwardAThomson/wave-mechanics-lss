// FPA nu sweep test: generate cosmological ICs, construct Madelung
// wavefunction with given nu, evolve via free-particle propagator,
// compare density to ZA reference.
//
// Usage: ./fpa_test_nu <nu_value>
// e.g.:  ./fpa_test_nu 1e-4

#include <cmath>
#include <complex>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <random>
#include <algorithm>
#include <fftw3.h>

using Complex = std::complex<double>;

// N is set at runtime via command line (default 64)
int N = 64;
int N3 = N * N * N;
constexpr double p = 1.0;
double dx = p / N;
constexpr double PI = M_PI;
constexpr double L_phys = 32.0;
constexpr double h_cosmo = 0.71;
constexpr double Omega_cdm = 0.27;
constexpr double Omega_Lambda = 0.73;
constexpr double Omega_b = 0.0;
constexpr double Omega_m = Omega_cdm + Omega_b;
constexpr double sigma_8_target = 0.81;
constexpr double n_s = 1.0;
constexpr double z_init = 50.0;
constexpr double Gamma_bbks = Omega_m * h_cosmo;
constexpr unsigned int SEED = 42;

inline int idx(int ix, int iy, int iz) {
    return (ix * N + iy) * N + iz;
}

double wavenumber(int index) {
    int k = (index <= N / 2) ? index : index - N;
    return 2.0 * PI * k / p;
}

// --- Growth factor ---
double hubble(double a) {
    return sqrt(Omega_m / (a * a * a) + Omega_Lambda);
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
    return 2.5 * Omega_m * hubble(a) * integral;
}

double compute_growth_ratio() {
    double a_init = 1.0 / (1.0 + z_init);
    return growth_factor(1.0) / growth_factor(a_init);
}

// --- BBKS ---
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

// --- Gaussian random field (identical to main code, same seed) ---
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

// --- Velocity potential from Poisson equation ---
// phi_v_hat = -delta_hat / k^2  (k != 0)
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

// --- FPA evolution ---
void fpa_evolve(const std::vector<Complex>& psi_init,
                std::vector<Complex>& psi_out, double D, double nu) {
    fftw_complex* in  = fftw_alloc_complex(N3);
    fftw_complex* out = fftw_alloc_complex(N3);

    for (int i = 0; i < N3; i++) {
        in[i][0] = psi_init[i].real();
        in[i][1] = psi_init[i].imag();
    }

    fftw_plan fwd = fftw_plan_dft_3d(N, N, N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(fwd);
    fftw_destroy_plan(fwd);

    // Apply free-particle propagator: exp(-i nu |k|^2 (D-1) / 2)
    for (int ix = 0; ix < N; ix++) {
        double kx = wavenumber(ix);
        for (int iy = 0; iy < N; iy++) {
            double ky = wavenumber(iy);
            for (int iz = 0; iz < N; iz++) {
                double kz = wavenumber(iz);
                int i = idx(ix, iy, iz);
                double k2 = kx * kx + ky * ky + kz * kz;
                double phase = -nu * k2 * (D - 1.0) / 2.0;
                Complex propagator(cos(phase), sin(phase));
                Complex psi_k(out[i][0], out[i][1]);
                psi_k *= propagator;
                out[i][0] = psi_k.real();
                out[i][1] = psi_k.imag();
            }
        }
    }

    fftw_plan bwd = fftw_plan_dft_3d(N, N, N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(bwd);
    fftw_destroy_plan(bwd);

    for (int i = 0; i < N3; i++) {
        psi_out[i] = Complex(in[i][0], in[i][1]) / static_cast<double>(N3);
    }

    fftw_free(in);
    fftw_free(out);
}

// --- CIC deposit (for ZA reference) ---
void cic_deposit(std::vector<double>& rho, double x, double y, double z) {
    x -= floor(x); y -= floor(y); z -= floor(z);
    double gx = x / dx - 0.5, gy = y / dx - 0.5, gz = z / dx - 0.5;
    int ix0 = (int)floor(gx), iy0 = (int)floor(gy), iz0 = (int)floor(gz);
    double fx = gx - ix0, fy = gy - iy0, fz = gz - iz0;
    for (int dix = 0; dix <= 1; dix++) {
        double wx = (dix == 0) ? (1.0 - fx) : fx;
        int iix = ((ix0 + dix) % N + N) % N;
        for (int diy = 0; diy <= 1; diy++) {
            double wy = (diy == 0) ? (1.0 - fy) : fy;
            int iiy = ((iy0 + diy) % N + N) % N;
            for (int diz = 0; diz <= 1; diz++) {
                double wz = (diz == 0) ? (1.0 - fz) : fz;
                int iiz = ((iz0 + diz) % N + N) % N;
                rho[idx(iix, iiy, iiz)] += wx * wy * wz;
            }
        }
    }
}

void compute_displacement(const fftw_complex* delta_hat,
                          std::vector<double>& sx, std::vector<double>& sy,
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
        fftw_plan ifft = fftw_plan_dft_3d(N, N, N, s_hat, s_real, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(ifft);
        fftw_destroy_plan(ifft);
        std::vector<double>& s = (comp == 0) ? sx : (comp == 1) ? sy : sz;
        for (int i = 0; i < N3; i++) s[i] = s_real[i][0] / static_cast<double>(N3);
    }
    fftw_free(s_hat);
    fftw_free(s_real);
}

void zel_dovich_density(const std::vector<double>& sx, const std::vector<double>& sy,
                        const std::vector<double>& sz, double D,
                        std::vector<double>& delta) {
    std::vector<double> rho(N3, 0.0);
    double disp_factor = D - 1.0;
    for (int ix = 0; ix < N; ix++) {
        double qx = (ix + 0.5) * dx;
        for (int iy = 0; iy < N; iy++) {
            double qy = (iy + 0.5) * dx;
            for (int iz = 0; iz < N; iz++) {
                double qz = (iz + 0.5) * dx;
                int i = idx(ix, iy, iz);
                cic_deposit(rho, qx + disp_factor * sx[i],
                                 qy + disp_factor * sy[i],
                                 qz + disp_factor * sz[i]);
            }
        }
    }
    double rho_mean = 0.0;
    for (int i = 0; i < N3; i++) rho_mean += rho[i];
    rho_mean /= N3;
    for (int i = 0; i < N3; i++) delta[i] = rho[i] / rho_mean - 1.0;
}

// Slab-averaged slice writer (same as main code)
constexpr int SLAB_THICK = 8;

void write_slice(const std::string& filename, double D,
                 const std::vector<double>& delta) {
    std::ofstream file(filename);
    file << "# D = " << D << "\n";
    file << "# z-slab average, thickness = " << SLAB_THICK << " cells\n";
    file << "# x  y  delta_avg\n";
    int iz_start = N / 2 - SLAB_THICK / 2;
    for (int ix = 0; ix < N; ix++) {
        double x = (ix + 0.5) * dx;
        for (int iy = 0; iy < N; iy++) {
            double y = (iy + 0.5) * dx;
            double avg = 0.0;
            for (int diz = 0; diz < SLAB_THICK; diz++) {
                int iz = (iz_start + diz + N) % N;
                avg += delta[idx(ix, iy, iz)];
            }
            avg /= SLAB_THICK;
            file << x << "  " << y << "  " << avg << "\n";
        }
        file << "\n";
    }
    file.close();
}

void write_density_3d(const std::string& filename, double D,
                      const std::vector<double>& delta) {
    std::ofstream file(filename, std::ios::binary);
    file.write(reinterpret_cast<const char*>(&D), sizeof(double));
    int n = N;
    file.write(reinterpret_cast<const char*>(&n), sizeof(int));
    double l = L_phys;
    file.write(reinterpret_cast<const char*>(&l), sizeof(double));
    file.write(reinterpret_cast<const char*>(delta.data()), N3 * sizeof(double));
    file.close();
}

// Format nu as a clean string for directory names (e.g. 1e-4 -> "1e-04")
std::string nu_label(double nu) {
    char buf[32];
    snprintf(buf, sizeof(buf), "%.0e", nu);
    return std::string(buf);
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <nu> [N=64]\n";
        return 1;
    }
    double nu = std::stod(argv[1]);
    if (argc >= 3) {
        N = std::stoi(argv[2]);
        N3 = N * N * N;
        dx = p / N;
    }

    // Output directory: output/fpa_N<N>_nu_<value>/
    std::string outdir = "output/fpa_N" + std::to_string(N) + "_nu_" + nu_label(nu);
    std::string za_outdir = "output/za_N" + std::to_string(N);
    // Create directories (system call, simple and portable enough)
    system(("mkdir -p " + outdir).c_str());
    system(("mkdir -p " + za_outdir).c_str());

    double D_ratio = compute_growth_ratio();
    double k_nyq = PI * N;

    // --- Generate ICs (identical to main code, same seed) ---
    fftw_complex* delta_hat = fftw_alloc_complex(N3);
    generate_gaussian_field(delta_hat, D_ratio);

    // delta(x)
    fftw_complex* work = fftw_alloc_complex(N3);
    fftw_plan ifft = fftw_plan_dft_3d(N, N, N, delta_hat, work, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(ifft);
    fftw_destroy_plan(ifft);
    std::vector<double> delta_init(N3);
    for (int i = 0; i < N3; i++) delta_init[i] = work[i][0] / static_cast<double>(N3);
    fftw_free(work);

    // Velocity potential phi_v(x)
    std::vector<double> phi_v(N3);
    compute_velocity_potential(delta_hat, phi_v);

    // Displacement field (for ZA reference and phase gradient check)
    std::vector<double> sx(N3), sy(N3), sz(N3);
    compute_displacement(delta_hat, sx, sy, sz);

    double s_max = 0.0;
    for (int i = 0; i < N3; i++) {
        double s = sqrt(sx[i]*sx[i] + sy[i]*sy[i] + sz[i]*sz[i]);
        if (s > s_max) s_max = s;
    }
    double phase_grad_max = s_max / nu;
    bool resolvable = phase_grad_max < k_nyq;

    // --- Construct Madelung wavefunction ---
    std::vector<Complex> psi_init(N3);
    for (int i = 0; i < N3; i++) {
        double amp = sqrt(std::max(1.0 + delta_init[i], 0.0));
        double phase = -phi_v[i] / nu;
        psi_init[i] = amp * Complex(cos(phase), sin(phase));
    }

    // psi_hat power distribution
    double power_low = 0.0, power_high = 0.0;
    {
        fftw_complex* psi_k = fftw_alloc_complex(N3);
        fftw_complex* in = fftw_alloc_complex(N3);
        for (int i = 0; i < N3; i++) {
            in[i][0] = psi_init[i].real();
            in[i][1] = psi_init[i].imag();
        }
        fftw_plan fwd = fftw_plan_dft_3d(N, N, N, in, psi_k, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(fwd);
        fftw_destroy_plan(fwd);
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
                    if (k < k_half) power_low += pw;
                    else power_high += pw;
                }
            }
        }
        fftw_free(psi_k);
        fftw_free(in);
    }
    double pct_low = 100.0 * power_low / (power_low + power_high);

    // --- Print summary to stdout ---
    printf("=== FPA nu sweep: nu = %g, N = %d ===\n", nu, N);
    printf("Output directory: %s\n", outdir.c_str());
    printf("Grid: N = %d^3 = %d, k_Nyquist = %.1f\n", N, N3, k_nyq);
    printf("s_max = %.6f (box units)\n", s_max);
    printf("max|grad(theta)| / k_Nyq = %.3f  -->  %s\n",
           phase_grad_max / k_nyq, resolvable ? "RESOLVABLE" : "ALIASED");
    printf("psi_hat power in low-k = %.1f%%\n\n", pct_low);

    // --- Evolve and write output ---
    std::vector<double> D_values = {1.0, 5.0, 15.0, 25.0, 34.0};
    std::vector<Complex> psi(N3);
    std::vector<double> delta_fpa(N3), delta_za(N3);

    // Summary file
    std::ofstream summary(outdir + "/summary.dat");
    summary << "# FPA nu sweep results\n";
    summary << "# nu = " << nu << "\n";
    summary << "# N = " << N << ", L_phys = " << L_phys << " Mpc/h\n";
    summary << "# k_Nyquist = " << k_nyq << "\n";
    summary << "# s_max = " << s_max << "\n";
    summary << "# max|grad(theta)|/k_Nyq = " << phase_grad_max / k_nyq << "\n";
    summary << "# Phase resolvable: " << (resolvable ? "YES" : "NO") << "\n";
    summary << "# psi_hat low-k power = " << pct_low << "%\n";
    summary << "#\n";
    summary << "# D   FPA_rms   ZA_rms   FPA_ZA_corr   mass_FPA   mass_ZA\n";

    printf("%-6s  %-12s %-12s %-12s %-12s %-12s\n",
           "D", "FPA_rms", "ZA_rms", "FPA-ZA corr", "mass_FPA", "mass_ZA");
    printf("------  ------------ ------------ ------------ ------------ ------------\n");

    for (double D : D_values) {
        // FPA evolution
        fpa_evolve(psi_init, psi, D, nu);
        double rho_mean_fpa = 0.0;
        for (int i = 0; i < N3; i++) {
            delta_fpa[i] = std::norm(psi[i]);
            rho_mean_fpa += delta_fpa[i];
        }
        rho_mean_fpa /= N3;
        for (int i = 0; i < N3; i++) delta_fpa[i] = delta_fpa[i] / rho_mean_fpa - 1.0;

        // ZA evolution
        zel_dovich_density(sx, sy, sz, D, delta_za);

        // Statistics
        double fpa_rms = 0.0, za_rms = 0.0;
        for (int i = 0; i < N3; i++) {
            fpa_rms += delta_fpa[i] * delta_fpa[i];
            za_rms += delta_za[i] * delta_za[i];
        }
        fpa_rms = sqrt(fpa_rms / N3);
        za_rms = sqrt(za_rms / N3);

        // Cross-correlation
        double corr_num = 0.0, norm1 = 0.0, norm2 = 0.0;
        for (int i = 0; i < N3; i++) {
            corr_num += delta_fpa[i] * delta_za[i];
            norm1 += delta_fpa[i] * delta_fpa[i];
            norm2 += delta_za[i] * delta_za[i];
        }
        double corr = corr_num / sqrt(norm1 * norm2);

        double mass_fpa = rho_mean_fpa;
        double mass_za = 0.0;
        for (int i = 0; i < N3; i++) mass_za += (1.0 + delta_za[i]);
        mass_za /= N3;

        printf("%-6.1f  %-12.6f %-12.6f %-12.6f %-12.10f %-12.10f\n",
               D, fpa_rms, za_rms, corr, mass_fpa, mass_za);
        summary << D << "  " << fpa_rms << "  " << za_rms << "  "
                << corr << "  " << mass_fpa << "  " << mass_za << "\n";

        // Write slice and 3D density for both FPA and ZA
        int Di = static_cast<int>(D);
        write_slice(outdir + "/fpa_slice_D" + std::to_string(Di) + ".dat", D, delta_fpa);
        write_density_3d(outdir + "/fpa_density_D" + std::to_string(Di) + ".bin", D, delta_fpa);
        write_slice(za_outdir + "/za_slice_D" + std::to_string(Di) + ".dat", D, delta_za);
        write_density_3d(za_outdir + "/za_density_D" + std::to_string(Di) + ".bin", D, delta_za);
    }

    summary.close();
    printf("\nFPA output written to: %s/\n", outdir.c_str());
    printf("ZA  output written to: %s/\n", za_outdir.c_str());

    fftw_free(delta_hat);
    return 0;
}
