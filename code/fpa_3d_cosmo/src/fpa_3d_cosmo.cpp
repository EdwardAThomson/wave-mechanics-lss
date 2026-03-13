// 3D Zel'dovich Approximation with Cosmological Initial Conditions
// Reproduces Figures 4.5-4.7 from Thomson (2011) Chapter 4.
//
// Uses BBKS transfer function (Eqs 3.35-3.38) to generate a Gaussian
// random density field at z_init=50, then evolves via the Zel'dovich
// approximation: x(q,D) = q - (D-1) * grad(phi_v(q))
//
// The ZA density is computed by particle displacement + CIC gridding,
// which is exact (no nu/phase-space aliasing artifacts).
//
// The FPA with Madelung transform is also available for comparison,
// but requires N >> 64 to avoid phase aliasing at the grid Nyquist
// frequency. At N=64, the Madelung phase exp(-i*phi_v/nu) oscillates
// faster than the grid can resolve, aliasing 99%+ of psi_hat power
// into unphysical high-k modes.
//
// Cosmology: h=0.71, Omega_cdm=0.27, Omega_Lambda=0.73, Omega_b=0,
//            sigma_8=0.81, n_s=1
// Grid: N=64^3, L=32 Mpc/h
// z_init=50 -> ICs are at D=1, z=0 corresponds to D ~ 39

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

// Grid parameters (box units: box = [0,1)^3)
constexpr int N = 64;
constexpr int N3 = N * N * N;
constexpr double p = 1.0;
constexpr double dx = p / N;
constexpr double PI = M_PI;

// Physical box size
constexpr double L_phys = 32.0;      // Mpc/h

// Cosmological parameters
constexpr double h_cosmo = 0.71;
constexpr double Omega_cdm = 0.27;
constexpr double Omega_Lambda = 0.73;
constexpr double Omega_b = 0.0;
constexpr double Omega_m = Omega_cdm + Omega_b;  // = 0.27
constexpr double sigma_8_target = 0.81;           // at z=0
constexpr double n_s = 1.0;
constexpr double z_init = 50.0;

// BBKS shape parameter (Eq 3.37)
constexpr double Gamma_bbks = Omega_m * h_cosmo;  // = 0.1917

// Random seed
constexpr unsigned int SEED = 42;

// 3D index: row-major (z fastest)
inline int idx(int ix, int iy, int iz) {
    return (ix * N + iy) * N + iz;
}

// Wavenumber for index along one axis (box units: k = 2*pi*n/p)
double wavenumber(int index) {
    int k = (index <= N / 2) ? index : index - N;
    return 2.0 * PI * k / p;
}

// ================================================================
// Growth factor computation
// ================================================================

double hubble(double a) {
    return sqrt(Omega_m / (a * a * a) + Omega_Lambda);
}

// Linear growth factor D(a) via Simpson's rule
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

// ================================================================
// BBKS transfer function and power spectrum
// ================================================================

double transfer_function(double k_phys) {
    if (k_phys == 0.0) return 1.0;
    double q = k_phys / Gamma_bbks;
    double log_term = log(1.0 + 2.34 * q) / (2.34 * q);
    double poly = 1.0 + 3.89 * q
                  + pow(16.1 * q, 2)
                  + pow(5.46 * q, 3)
                  + pow(6.71 * q, 4);
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
// Gaussian random field generation
// ================================================================

void generate_gaussian_field(fftw_complex* delta_hat, double D_ratio) {
    double sigma8_sq_raw = compute_sigma8_sq_unnorm();
    double A = (sigma_8_target * sigma_8_target) / sigma8_sq_raw;
    double A_init = A / (D_ratio * D_ratio);
    std::cout << "sigma_8 raw (unnormalised) = " << sqrt(sigma8_sq_raw) << "\n";
    std::cout << "Normalisation A (z=0) = " << A << "\n";
    std::cout << "Normalisation A (z_init) = " << A_init << "\n";

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
// Zel'dovich approximation via particle displacement + CIC
// ================================================================

// Cloud-In-Cell assignment: deposit unit mass from position (x,y,z)
// onto the density grid. Positions are in box units [0,1).
void cic_deposit(std::vector<double>& rho, double x, double y, double z) {
    // Wrap to [0,1)
    x -= floor(x);
    y -= floor(y);
    z -= floor(z);

    // Grid position (cell-centered)
    double gx = x / dx - 0.5;
    double gy = y / dx - 0.5;
    double gz = z / dx - 0.5;

    int ix0 = static_cast<int>(floor(gx));
    int iy0 = static_cast<int>(floor(gy));
    int iz0 = static_cast<int>(floor(gz));

    double fx = gx - ix0;
    double fy = gy - iy0;
    double fz = gz - iz0;

    // Deposit to 8 surrounding cells
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

// Compute Zel'dovich density at growth factor D.
// Displacement: x(q, D) = q - (D-1) * grad(phi_v(q))
// where phi_hat = -delta_hat / k^2.
void zel_dovich_density(const std::vector<double>& sx,
                        const std::vector<double>& sy,
                        const std::vector<double>& sz,
                        double D, std::vector<double>& delta) {
    // Initialise density to zero
    std::vector<double> rho(N3, 0.0);

    double disp_factor = D - 1.0;

    // Displace each particle and deposit via CIC
    for (int ix = 0; ix < N; ix++) {
        double qx = (ix + 0.5) * dx;  // Cell-centered Lagrangian position
        for (int iy = 0; iy < N; iy++) {
            double qy = (iy + 0.5) * dx;
            for (int iz = 0; iz < N; iz++) {
                double qz = (iz + 0.5) * dx;
                int i = idx(ix, iy, iz);

                // Eulerian position
                double ex = qx + disp_factor * sx[i];
                double ey = qy + disp_factor * sy[i];
                double ez = qz + disp_factor * sz[i];

                cic_deposit(rho, ex, ey, ez);
            }
        }
    }

    // Convert to density contrast: delta = rho / rho_mean - 1
    double rho_mean = 0.0;
    for (int i = 0; i < N3; i++) rho_mean += rho[i];
    rho_mean /= N3;

    for (int i = 0; i < N3; i++) {
        delta[i] = rho[i] / rho_mean - 1.0;
    }
}

// ================================================================
// Compute displacement field from Poisson equation
// s = -grad(phi_v) where phi_hat = -delta_hat / k^2
// In Fourier space: s_hat_i = -ik_i * phi_hat = ik_i * delta_hat / k^2
// ================================================================

void compute_displacement(const fftw_complex* delta_hat,
                          std::vector<double>& sx,
                          std::vector<double>& sy,
                          std::vector<double>& sz) {
    fftw_complex* s_hat = fftw_alloc_complex(N3);
    fftw_complex* s_real = fftw_alloc_complex(N3);

    // Compute each component
    for (int comp = 0; comp < 3; comp++) {
        for (int ix = 0; ix < N; ix++) {
            double kx = wavenumber(ix);
            for (int iy = 0; iy < N; iy++) {
                double ky = wavenumber(iy);
                for (int iz = 0; iz < N; iz++) {
                    double kz = wavenumber(iz);
                    int i = idx(ix, iy, iz);
                    double k2 = kx * kx + ky * ky + kz * kz;

                    if (k2 == 0.0) {
                        s_hat[i][0] = 0.0;
                        s_hat[i][1] = 0.0;
                        continue;
                    }

                    // s_hat_j = i * k_j * delta_hat / k^2
                    // = i * k_j / k^2 * (re + i*im)
                    // = (-k_j * im / k^2) + i * (k_j * re / k^2)
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
        for (int i = 0; i < N3; i++) {
            s[i] = s_real[i][0] / static_cast<double>(N3);
        }
    }

    fftw_free(s_hat);
    fftw_free(s_real);
}

// ================================================================
// Output
// ================================================================

// Slab thickness for z-averaging (in cells)
constexpr int SLAB_THICK = 8;  // 8 cells = 4 Mpc/h at N=64, L=32

// Write a z-slab-averaged density projection
void write_slice(const std::string& filename, double D,
                 const std::vector<double>& delta) {
    std::ofstream file(filename);
    file << "# D = " << D << "\n";
    file << "# z-slab average, thickness = " << SLAB_THICK << " cells centred at iz=N/2\n";
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

void write_histogram(const std::string& filename,
                     const std::vector<double>& data,
                     int nbins, const std::string& label) {
    double vmin = *std::min_element(data.begin(), data.end());
    double vmax = *std::max_element(data.begin(), data.end());
    double bin_width = (vmax - vmin) / nbins;

    std::vector<int> counts(nbins, 0);
    for (double v : data) {
        int bin = static_cast<int>((v - vmin) / bin_width);
        if (bin >= nbins) bin = nbins - 1;
        if (bin < 0) bin = 0;
        counts[bin]++;
    }

    std::ofstream file(filename);
    file << "# " << label << " histogram\n";
    file << "# bin_center  count  normalised\n";
    for (int i = 0; i < nbins; i++) {
        double center = vmin + (i + 0.5) * bin_width;
        double norm = static_cast<double>(counts[i]) / (data.size() * bin_width);
        file << center << "  " << counts[i] << "  " << norm << "\n";
    }
    file.close();
}

// ================================================================
// Main
// ================================================================

int main() {
    std::cout << "3D Zel'dovich Approximation with Cosmological ICs\n";
    std::cout << "N = " << N << "^3 = " << N3
              << ", L_phys = " << L_phys << " Mpc/h"
              << ", box = [0," << p << ")^3\n";
    std::cout << "Cosmology: h=" << h_cosmo << " Omega_m=" << Omega_m
              << " sigma_8=" << sigma_8_target << " n_s=" << n_s
              << " z_init=" << z_init << "\n";
    std::cout << "BBKS shape parameter Gamma = " << Gamma_bbks << "\n\n";

    // --- Compute growth factor ratio ---
    double D_ratio = compute_growth_ratio();
    double a_init = 1.0 / (1.0 + z_init);
    std::cout << "a_init = " << a_init << "\n";
    std::cout << "D(z=0)/D(z_init) = " << D_ratio << "\n";
    std::cout << "sigma_8 at z_init = " << sigma_8_target / D_ratio << "\n\n";

    // --- Generate Gaussian random density field ---
    std::cout << "Generating Gaussian random field at z_init...\n";
    fftw_complex* delta_hat = fftw_alloc_complex(N3);
    generate_gaussian_field(delta_hat, D_ratio);

    // IFFT to get delta(x)
    fftw_complex* work = fftw_alloc_complex(N3);
    fftw_plan ifft = fftw_plan_dft_3d(N, N, N, delta_hat, work,
                                       FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(ifft);
    fftw_destroy_plan(ifft);

    std::vector<double> delta_init(N3);
    double max_imag = 0.0;
    for (int i = 0; i < N3; i++) {
        delta_init[i] = work[i][0] / static_cast<double>(N3);
        max_imag = std::max(max_imag, std::abs(work[i][1] / N3));
    }
    fftw_free(work);
    std::cout << "Max imaginary part of delta(x): " << max_imag << "\n";

    // Statistics
    double delta_mean = 0.0, delta_rms = 0.0;
    double delta_min = delta_init[0], delta_max = delta_init[0];
    for (int i = 0; i < N3; i++) {
        delta_mean += delta_init[i];
        delta_rms += delta_init[i] * delta_init[i];
        if (delta_init[i] < delta_min) delta_min = delta_init[i];
        if (delta_init[i] > delta_max) delta_max = delta_init[i];
    }
    delta_mean /= N3;
    delta_rms = sqrt(delta_rms / N3);
    std::cout << "Initial delta: mean=" << delta_mean << " rms=" << delta_rms
              << " min=" << delta_min << " max=" << delta_max << "\n";

    write_histogram("output/za/delta_init_hist.dat", delta_init, 100, "initial delta");

    // --- Compute displacement field ---
    std::cout << "\nComputing displacement field from Poisson equation...\n";
    std::vector<double> sx(N3), sy(N3), sz(N3);
    compute_displacement(delta_hat, sx, sy, sz);
    fftw_free(delta_hat);

    // Displacement statistics
    double s_rms = 0.0, s_max = 0.0;
    for (int i = 0; i < N3; i++) {
        double s2 = sx[i]*sx[i] + sy[i]*sy[i] + sz[i]*sz[i];
        s_rms += s2;
        double s = sqrt(s2);
        if (s > s_max) s_max = s;
    }
    s_rms = sqrt(s_rms / (3.0 * N3));
    std::cout << "Displacement rms (per component, box units): " << s_rms << "\n";
    std::cout << "Displacement max (box units): " << s_max << "\n";

    // Velocity = displacement at D=2 (one unit of D)
    // Physical velocity = s * L_phys * H_0 * f(Omega) where f ~ Omega^0.55
    // In box units, v = s (since dx/dD = s in ZA)
    write_histogram("output/za/vel_x_init_hist.dat", sx, 100, "initial s_x (box units)");
    write_histogram("output/za/vel_y_init_hist.dat", sy, 100, "initial s_y (box units)");
    write_histogram("output/za/vel_z_init_hist.dat", sz, 100, "initial s_z (box units)");

    // --- Zel'dovich evolution ---
    std::cout << "\nZel'dovich approximation evolution...\n";
    std::vector<double> D_values = {1.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 34.0};
    std::vector<double> delta(N3);

    for (double D : D_values) {
        zel_dovich_density(sx, sy, sz, D, delta);

        double dmin = delta[0], dmax = delta[0];
        double drms = 0.0;
        for (int i = 0; i < N3; i++) {
            if (delta[i] < dmin) dmin = delta[i];
            if (delta[i] > dmax) dmax = delta[i];
            drms += delta[i] * delta[i];
        }
        drms = sqrt(drms / N3);

        // Mass conservation check: sum of rho should = N^3
        double mass = 0.0;
        for (int i = 0; i < N3; i++) mass += (1.0 + delta[i]);
        mass /= N3;

        printf("D = %5.1f  delta_rms = %.6f  delta = [%.4f, %.4f]  <rho>/rho_mean = %.15f\n",
               D, drms, dmin, dmax, mass);

        std::string slice_file = "output/za/za_slice_D" + std::to_string(static_cast<int>(D)) + ".dat";
        write_slice(slice_file, D, delta);

        std::string density_file = "output/za/za_density_D" + std::to_string(static_cast<int>(D)) + ".bin";
        write_density_3d(density_file, D, delta);
    }

    std::cout << "\nOutput files written to output/za/:\n";
    std::cout << "  za_slice_D*.dat       - z-slab averaged density slices\n";
    std::cout << "  za_density_D*.bin     - full 3D density (binary)\n";
    std::cout << "  delta_init_hist.dat   - initial delta histogram\n";
    std::cout << "  vel_*_init_hist.dat   - initial displacement histograms\n";

    return 0;
}
