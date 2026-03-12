// 3D Free Particle Approximation (FPA) - Toy Model
// Reproduces Figure 4.4 from Thomson (2011) Chapter 4.
//
// 3D extension of the FPA with sinusoidal initial conditions.
// The FPA sets V=0 in the Schrödinger equation, giving:
//   ψ(x,D) = IFFT[ ψ̂(k) * exp(-iν |k|² D / 2) ]
//
// Initial conditions (Eqs 4.18-4.20):
//   δ_init = -δ_a cos(2πx/p) cos(2πy/p) cos(2πz/p)
//   φ_v,init = -(p/2π)² δ_init   (Eq 4.19)
//   ψ_init = (1+δ)^{1/2} exp(-i φ_v / ν)
//
// Parameters: δ_a = 0.01, p = 1.0, N = 64 per dimension
// Note: ν must be small enough for classical (Zel'dovich-like) growth.
// The thesis states ν=1 for the 1D case but does not explicitly state ν
// for the 3D toy model. With ν=1 the FPA gives quantum-like oscillations
// rather than growth; with ν ≈ 1e-5 the results match Figure 4.4.

#include <cmath>
#include <complex>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <fftw3.h>

using Complex = std::complex<double>;

constexpr int N = 64;
constexpr int N3 = N * N * N;
constexpr double delta_a = 0.01;
constexpr double nu = 1e-5;  // Small ν → classical (Zel'dovich) regime
constexpr double p = 1.0;
constexpr double dx = p / N;
constexpr double PI = M_PI;

// 3D index: row-major (z fastest, then y, then x)
inline int idx(int ix, int iy, int iz) {
    return (ix * N + iy) * N + iz;
}

// Wavenumber for index along one axis
double wavenumber(int index) {
    int k = (index <= N / 2) ? index : index - N;
    return 2.0 * PI * k / p;
}

// Initial conditions (Eqs 4.18-4.20)
void initialise(std::vector<Complex>& psi) {
    double poisson_factor = (p / (2.0 * PI)) * (p / (2.0 * PI));
    for (int ix = 0; ix < N; ix++) {
        double x = ix * dx;
        for (int iy = 0; iy < N; iy++) {
            double y = iy * dx;
            for (int iz = 0; iz < N; iz++) {
                double z = iz * dx;
                double delta_i = -delta_a * cos(2.0 * PI * x / p)
                                          * cos(2.0 * PI * y / p)
                                          * cos(2.0 * PI * z / p);
                double phi_v_i = -poisson_factor * delta_i;
                double amplitude = sqrt(1.0 + delta_i);
                double phase = -phi_v_i / nu;
                psi[idx(ix, iy, iz)] = amplitude * Complex(cos(phase), sin(phase));
            }
        }
    }
}

// Evolve to growth factor D using 3D FFT
void evolve(const std::vector<Complex>& psi_init,
            std::vector<Complex>& psi_out, double D) {
    fftw_complex* in  = fftw_alloc_complex(N3);
    fftw_complex* out = fftw_alloc_complex(N3);

    // Copy input
    for (int i = 0; i < N3; i++) {
        in[i][0] = psi_init[i].real();
        in[i][1] = psi_init[i].imag();
    }

    // Forward 3D FFT
    fftw_plan fwd = fftw_plan_dft_3d(N, N, N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(fwd);
    fftw_destroy_plan(fwd);

    // Apply free-particle propagator in k-space
    // ψ̂(k, D) = ψ̂(k, 1) * exp(-i ν |k|² (D-1) / 2)
    for (int ix = 0; ix < N; ix++) {
        double kx = wavenumber(ix);
        for (int iy = 0; iy < N; iy++) {
            double ky = wavenumber(iy);
            for (int iz = 0; iz < N; iz++) {
                double kz = wavenumber(iz);
                int i = idx(ix, iy, iz);
                double k2 = kx * kx + ky * ky + kz * kz;
                double phase = -nu * k2 * D / 2.0;
                Complex propagator(cos(phase), sin(phase));
                Complex psi_k(out[i][0], out[i][1]);
                psi_k *= propagator;
                out[i][0] = psi_k.real();
                out[i][1] = psi_k.imag();
            }
        }
    }

    // Inverse 3D FFT
    fftw_plan bwd = fftw_plan_dft_3d(N, N, N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(bwd);
    fftw_destroy_plan(bwd);

    // Normalise
    for (int i = 0; i < N3; i++) {
        psi_out[i] = Complex(in[i][0], in[i][1]) / static_cast<double>(N3);
    }

    fftw_free(in);
    fftw_free(out);
}

// Compute density contrast δ = |ψ|² / <|ψ|²> - 1
void compute_density(const std::vector<Complex>& psi,
                     std::vector<double>& delta) {
    double rho_mean = 0.0;
    for (int i = 0; i < N3; i++) {
        double rho = std::norm(psi[i]);
        delta[i] = rho;
        rho_mean += rho;
    }
    rho_mean /= N3;
    for (int i = 0; i < N3; i++) {
        delta[i] = delta[i] / rho_mean - 1.0;
    }
}

// Compute velocity using probability current: v = ν * Im(∇ψ / ψ)
// Central differences with periodic boundaries
// Returns 3 components (vx, vy, vz) for each point
void compute_velocity(const std::vector<Complex>& psi,
                      std::vector<double>& vx,
                      std::vector<double>& vy,
                      std::vector<double>& vz) {
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

                Complex grad_x = (psi[idx(ixp, iy, iz)] - psi[idx(ixm, iy, iz)]) / (2.0 * dx);
                Complex grad_y = (psi[idx(ix, iyp, iz)] - psi[idx(ix, iym, iz)]) / (2.0 * dx);
                Complex grad_z = (psi[idx(ix, iy, izp)] - psi[idx(ix, iy, izm)]) / (2.0 * dx);

                Complex ratio_x = grad_x / psi[i];
                Complex ratio_y = grad_y / psi[i];
                Complex ratio_z = grad_z / psi[i];

                vx[i] = nu * ratio_x.imag();
                vy[i] = nu * ratio_y.imag();
                vz[i] = nu * ratio_z.imag();
            }
        }
    }
}

// Compute total mass
double total_mass(const std::vector<Complex>& psi) {
    double mass = 0.0;
    double dV = dx * dx * dx;
    for (int i = 0; i < N3; i++) {
        mass += std::norm(psi[i]) * dV;
    }
    return mass;
}

// Write a z=0 slice of density to file (for plotting like Fig 4.4)
void write_slice(const std::string& filename, double D,
                 const std::vector<double>& delta,
                 const std::vector<double>& vx,
                 const std::vector<double>& vy,
                 const std::vector<double>& vz) {
    std::ofstream file(filename);
    file << "# D = " << D << "\n";
    file << "# z-slice at iz = N/2 (z = p/2)\n";
    file << "# x/p  y/p  delta  vx  vy  vz\n";

    int iz = N / 2;  // z = p/2 slice (where the overdense peak is centred)
    for (int ix = 0; ix < N; ix++) {
        double x = ix * dx;
        for (int iy = 0; iy < N; iy++) {
            double y = iy * dx;
            int i = idx(ix, iy, iz);
            file << x / p << "  " << y / p << "  "
                 << delta[i] << "  "
                 << vx[i] << "  " << vy[i] << "  " << vz[i] << "\n";
        }
        file << "\n";  // blank line between x-rows for gnuplot compatibility
    }
    file.close();
}

int main() {
    // Output times matching thesis Figure 4.4
    // D = 0 is the initial time (no evolution)
    std::vector<double> D_values = {0.0, 10.0, 30.0, 50.0};

    // Initialise wavefunction at D = 0
    std::vector<Complex> psi_init(N3);
    initialise(psi_init);

    std::cout << "3D Free Particle Approximation (Toy Model)\n";
    std::cout << "N = " << N << "^3 = " << N3
              << ", delta_a = " << delta_a
              << ", nu = " << nu << ", p = " << p << "\n";
    std::cout << "Initial mass = " << total_mass(psi_init) << "\n\n";

    std::vector<Complex> psi(N3);
    std::vector<double> delta(N3);
    std::vector<double> vx(N3), vy(N3), vz(N3);

    for (double D : D_values) {
        evolve(psi_init, psi, D);
        compute_density(psi, delta);
        compute_velocity(psi, vx, vy, vz);

        double mass = total_mass(psi);

        // Statistics
        double delta_min = delta[0], delta_max = delta[0];
        for (int i = 1; i < N3; i++) {
            if (delta[i] < delta_min) delta_min = delta[i];
            if (delta[i] > delta_max) delta_max = delta[i];
        }

        printf("D = %3.0f  mass = %.15f  delta = [%.6f, %.6f]\n",
               D, mass, delta_min, delta_max);

        // Write z=0 slice
        std::string filename = "fpa_3d_D" + std::to_string(static_cast<int>(D)) + ".dat";
        write_slice(filename, D, delta, vx, vy, vz);
    }

    std::cout << "\nOutput files written: fpa_3d_D*.dat (z=0 slices)\n";

    return 0;
}
