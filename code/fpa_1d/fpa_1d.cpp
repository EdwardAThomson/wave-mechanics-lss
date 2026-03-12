// 1D Free Particle Approximation (FPA)
// Reproduces the results of Short & Coles (2006) and Chapter 4 of
// Thomson (2011) "Schrödinger Wave-mechanics and Large Scale Structure"
//
// The FPA sets V=0 in the Schrödinger equation, giving an analytic
// Fourier-space solution. The wavefunction at any growth factor D is:
//   ψ(x,D) = IFFT[ ψ̂(k) * exp(-iν k² (D-1) / 2) ]
// where ψ̂(k) = FFT[ψ(x, D=1)] and D=1 is the initial time.
//
// Density: ρ = |ψ|² = ψψ*
// Density contrast: δ = ρ/<ρ> - 1
// Velocity: v = (ν) * Im(∇ψ / ψ)  [probability current method]
//
// Parameters (matching Short's thesis):
//   δ_a = 0.01, ν = 1.0, p = 1.0, N = 512
//   Domain: x ∈ [0, 1)

#include <cmath>
#include <complex>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <fftw3.h>

using Complex = std::complex<double>;

constexpr int N = 512;
constexpr double delta_a = 0.01;
constexpr double nu = 1.0;
constexpr double p = 1.0;
constexpr double dx = p / N;
constexpr double PI = M_PI;

// Wavenumbers for a periodic domain of length p with N points.
// FFTW stores: k = 0, 1, 2, ..., N/2, -(N/2-1), ..., -1
// Physical wavenumber: k_phys = 2π * k_index / p
double wavenumber(int index) {
    int k = (index <= N / 2) ? index : index - N;
    return 2.0 * PI * k / p;
}

// Initial conditions (Eq 4.11)
void initialise(std::vector<Complex>& psi) {
    for (int i = 0; i < N; i++) {
        double x = i * dx;
        double delta_i = -delta_a * cos(2.0 * PI * x / p);
        double phi_v_i = -(p / (2.0 * PI)) * (p / (2.0 * PI)) * delta_i;
        double amplitude = sqrt(1.0 + delta_i);
        double phase = -phi_v_i / nu;
        psi[i] = amplitude * Complex(cos(phase), sin(phase));
    }
}

// Evolve to growth factor D using FFT (Eq 4.5)
void evolve(const std::vector<Complex>& psi_init,
            std::vector<Complex>& psi_out, double D) {
    // Allocate FFTW arrays
    fftw_complex* in  = fftw_alloc_complex(N);
    fftw_complex* out = fftw_alloc_complex(N);

    // Copy input
    for (int i = 0; i < N; i++) {
        in[i][0] = psi_init[i].real();
        in[i][1] = psi_init[i].imag();
    }

    // Forward FFT
    fftw_plan fwd = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(fwd);
    fftw_destroy_plan(fwd);

    // Apply free-particle propagator in k-space
    // ψ̂(k, D) = ψ̂(k, 1) * exp(-i ν k² (D-1) / 2)
    for (int i = 0; i < N; i++) {
        double k = wavenumber(i);
        double phase = -nu * k * k * (D - 1.0) / 2.0;
        Complex propagator(cos(phase), sin(phase));
        Complex psi_k(out[i][0], out[i][1]);
        psi_k *= propagator;
        out[i][0] = psi_k.real();
        out[i][1] = psi_k.imag();
    }

    // Inverse FFT
    fftw_plan bwd = fftw_plan_dft_1d(N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(bwd);
    fftw_destroy_plan(bwd);

    // FFTW inverse is unnormalised: divide by N
    for (int i = 0; i < N; i++) {
        psi_out[i] = Complex(in[i][0], in[i][1]) / static_cast<double>(N);
    }

    fftw_free(in);
    fftw_free(out);
}

// Compute density contrast δ = |ψ|² / <|ψ|²> - 1
void compute_density(const std::vector<Complex>& psi,
                     std::vector<double>& delta) {
    double rho_mean = 0.0;
    for (int i = 0; i < N; i++) {
        double rho = std::norm(psi[i]);  // |ψ|²
        delta[i] = rho;
        rho_mean += rho;
    }
    rho_mean /= N;
    for (int i = 0; i < N; i++) {
        delta[i] = delta[i] / rho_mean - 1.0;
    }
}

// Compute velocity using probability current: v = ν * Im(∇ψ / ψ)
// Using central differences for ∇ψ with periodic boundaries
void compute_velocity(const std::vector<Complex>& psi,
                      std::vector<double>& velocity) {
    for (int i = 0; i < N; i++) {
        int ip = (i + 1) % N;
        int im = (i - 1 + N) % N;
        Complex grad_psi = (psi[ip] - psi[im]) / (2.0 * dx);
        Complex ratio = grad_psi / psi[i];
        velocity[i] = nu * ratio.imag();
    }
}

// Compute total mass (should be conserved)
double total_mass(const std::vector<Complex>& psi) {
    double mass = 0.0;
    for (int i = 0; i < N; i++) {
        mass += std::norm(psi[i]) * dx;
    }
    return mass;
}

// Write output to file
void write_output(const std::string& filename, double D,
                  const std::vector<double>& delta,
                  const std::vector<double>& velocity) {
    std::ofstream file(filename);
    file << "# D = " << D << "\n";
    file << "# x/p  delta  log(2+delta)  velocity\n";
    for (int i = 0; i < N; i++) {
        double x = i * dx;
        file << x / p << "  "
             << delta[i] << "  "
             << log10(2.0 + delta[i]) << "  "
             << velocity[i] << "\n";
    }
    file.close();
}

int main() {
    // Output times matching thesis Figure 4.1
    std::vector<double> D_values = {1.0, 59.06, 117.16, 174.98};

    // Initialise wavefunction at D = 1
    std::vector<Complex> psi_init(N);
    initialise(psi_init);

    std::cout << "1D Free Particle Approximation\n";
    std::cout << "N = " << N << ", delta_a = " << delta_a
              << ", nu = " << nu << ", p = " << p << "\n";
    std::cout << "Initial mass = " << total_mass(psi_init) << "\n\n";

    std::vector<Complex> psi(N);
    std::vector<double> delta(N);
    std::vector<double> velocity(N);

    for (double D : D_values) {
        evolve(psi_init, psi, D);
        compute_density(psi, delta);
        compute_velocity(psi, velocity);

        double mass = total_mass(psi);

        // Compute some statistics
        double delta_min = delta[0], delta_max = delta[0];
        for (int i = 1; i < N; i++) {
            if (delta[i] < delta_min) delta_min = delta[i];
            if (delta[i] > delta_max) delta_max = delta[i];
        }

        std::cout << "D = " << D
                  << "  mass = " << mass
                  << "  delta_min = " << delta_min
                  << "  delta_max = " << delta_max << "\n";

        // Write to file
        std::string filename = "fpa_1d_D" + std::to_string(static_cast<int>(D)) + ".dat";
        write_output(filename, D, delta, velocity);
    }

    std::cout << "\nShell crossing expected at D = 101\n";
    std::cout << "Output files written: fpa_1d_D*.dat\n";

    return 0;
}
