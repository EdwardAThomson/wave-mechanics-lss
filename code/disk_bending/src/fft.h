// Thin RAII wrapper around a 1D complex-to-complex FFTW plan pair.
//
// Used for the kinetic half-step and for the periodic Poisson solve. The
// isolated-BC vertical Poisson solve deliberately does NOT use this: see
// poisson_z.h.
#pragma once

#include <complex>
#include <vector>

#include <fftw3.h>

using Complex = std::complex<double>;

class FFT1D {
public:
    explicit FFT1D(int n) : N_(n) {
        buf_ = fftw_alloc_complex(N_);
        fwd_ = fftw_plan_dft_1d(N_, buf_, buf_, FFTW_FORWARD, FFTW_MEASURE);
        bwd_ = fftw_plan_dft_1d(N_, buf_, buf_, FFTW_BACKWARD, FFTW_MEASURE);
    }
    ~FFT1D() {
        fftw_destroy_plan(fwd_);
        fftw_destroy_plan(bwd_);
        fftw_free(buf_);
    }
    FFT1D(const FFT1D&) = delete;
    FFT1D& operator=(const FFT1D&) = delete;

    int size() const { return N_; }

    void load(const std::vector<Complex>& a) {
        for (int j = 0; j < N_; ++j) {
            buf_[j][0] = a[j].real();
            buf_[j][1] = a[j].imag();
        }
    }
    void load(const std::vector<double>& a) {
        for (int j = 0; j < N_; ++j) {
            buf_[j][0] = a[j];
            buf_[j][1] = 0.0;
        }
    }
    void store(std::vector<Complex>& a) const {
        for (int j = 0; j < N_; ++j) a[j] = Complex(buf_[j][0], buf_[j][1]);
    }
    void store_real(std::vector<double>& a) const {
        for (int j = 0; j < N_; ++j) a[j] = buf_[j][0];
    }

    void forward() { fftw_execute(fwd_); }
    // Backward transform including the 1/N normalisation
    void backward() {
        fftw_execute(bwd_);
        const double s = 1.0 / static_cast<double>(N_);
        for (int j = 0; j < N_; ++j) {
            buf_[j][0] *= s;
            buf_[j][1] *= s;
        }
    }

    Complex at(int j) const { return Complex(buf_[j][0], buf_[j][1]); }
    void set(int j, Complex v) {
        buf_[j][0] = v.real();
        buf_[j][1] = v.imag();
    }
    void scale(int j, double s) {
        buf_[j][0] *= s;
        buf_[j][1] *= s;
    }

private:
    int N_;
    fftw_complex* buf_;
    fftw_plan fwd_, bwd_;
};
