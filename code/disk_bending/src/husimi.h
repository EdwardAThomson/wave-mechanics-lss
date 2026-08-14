// Husimi transform: |psi|^2 -> f(z, v_z).
//
// This is where Schrodinger-Poisson earns its keep against N-body. The
// distribution function comes out of a single wavefunction on a grid with no
// particle sampling at all, so the Gaia-style vertical phase spiral appears
// with zero shot noise, at whatever contrast the phase-space resolution
// hbar_eff allows.
//
//   F(z0, v) = (1 / 2 pi hbar) | <g_{z0,v} | psi> |^2
//   g_{z0,v}(z) = (2 pi s^2)^{-1/4} exp( -(z - z0)^2 / 4 s^2 + i v z / hbar )
//
// The coherent-state window has widths (s, hbar/2s) in (z, v), so their
// product is the irreducible hbar_eff/2 phase-space cell. Choosing
// s = sqrt(hbar_eff h / sigma_z) balances the two against the sheet's own
// aspect ratio, which is the default here.
//
// For a multi-stream ensemble the streams are summed incoherently, which is
// the correct statement that they are independent samples of the same f.
// The result is normalised so that the double integral over (z, v) is the
// surface density, sidestepping any convention slip in the prefactor.
#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "evolve.h"
#include "fft.h"
#include "grid.h"
#include "units.h"

#ifdef _OPENMP
#include <omp.h>
#endif

struct HusimiMap {
    int nz = 0, nv = 0;
    double zmin = 0.0, zmax = 0.0;
    double vmin = 0.0, vmax = 0.0;
    std::vector<double> z, v;
    std::vector<double> F;  // row-major [iz * nv + iv], Msun/kpc^2 per (kpc km/s)

    double at(int iz, int iv) const { return F[static_cast<size_t>(iz) * nv + iv]; }
    double& at(int iz, int iv) { return F[static_cast<size_t>(iz) * nv + iv]; }
};

struct HusimiConfig {
    double z_half = 1.0;   // output half-range in z, kpc
    double v_half = 80.0;  // output half-range in v_z, km/s
    int nz = 256;          // output samples in z
    double s = 0.0;        // coherent-state width in z, kpc; 0 = auto
};

inline HusimiMap husimi_transform(const Evolver& ev, const SlabState& st,
                                  const HusimiConfig& cfg,
                                  double h_scale, double sigma_scale) {
    const Grid1D& g = ev.grid();
    const int N = g.N;
    const double hb = ev.hbar();

    double s = cfg.s;
    if (s <= 0.0) s = std::sqrt(hb * h_scale / sigma_scale);

    // Velocity samples come from the native FFT wavenumbers, v = hbar k.
    std::vector<int> modes;
    for (int j = 0; j < N; ++j) {
        const double v = hb * g.k(j);
        if (std::fabs(v) <= cfg.v_half) modes.push_back(j);
    }
    // Sort by velocity so the output is monotonic.
    std::sort(modes.begin(), modes.end(), [&](int a, int b) {
        return g.k(a) < g.k(b);
    });

    HusimiMap map;
    map.nz = cfg.nz;
    map.nv = static_cast<int>(modes.size());
    map.zmin = -cfg.z_half;
    map.zmax = cfg.z_half;
    map.z.resize(map.nz);
    map.v.resize(map.nv);
    for (int i = 0; i < map.nz; ++i) {
        map.z[i] = map.zmin + (i + 0.5) * (map.zmax - map.zmin) / map.nz;
    }
    for (int i = 0; i < map.nv; ++i) map.v[i] = hb * g.k(modes[i]);
    map.vmin = map.v.front();
    map.vmax = map.v.back();
    map.F.assign(static_cast<size_t>(map.nz) * map.nv, 0.0);

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    std::vector<std::unique_ptr<FFT1D>> ffts;
    for (int t = 0; t < nthreads; ++t) ffts.emplace_back(new FFT1D(N));

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int iz = 0; iz < map.nz; ++iz) {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        FFT1D& f = *ffts[tid];
        std::vector<Complex> windowed(N);
        const double z0 = map.z[iz];
        for (int sidx = 0; sidx < st.n_streams(); ++sidx) {
            for (int j = 0; j < N; ++j) {
                const double dz = g.z(j) - z0;
                const double wgt = std::exp(-dz * dz / (4.0 * s * s));
                windowed[j] = st.psi[sidx][j] * wgt;
            }
            f.load(windowed);
            f.forward();
            for (int iv = 0; iv < map.nv; ++iv) {
                map.at(iz, iv) += std::norm(f.at(modes[iv]));
            }
        }
    }

    // Normalise so that the (z, v) integral returns the surface density that
    // actually lies inside the output window.
    const double dzo = (map.zmax - map.zmin) / map.nz;
    const double dvo = (map.nv > 1) ? (map.v[1] - map.v[0]) : 1.0;
    double tot = 0.0;
    for (double x : map.F) tot += x;
    tot *= dzo * dvo;

    std::vector<double> rho;
    st.density(rho);
    double mass_in = 0.0;
    for (int j = 0; j < N; ++j) {
        if (g.z(j) >= map.zmin && g.z(j) <= map.zmax) mass_in += rho[j];
    }
    mass_in *= g.dz;

    if (tot > 0.0) {
        const double scale = mass_in / tot;
        for (double& x : map.F) x *= scale;
    }
    return map;
}

inline void write_husimi(const std::string& path, const HusimiMap& m, double t) {
    FILE* f = std::fopen(path.c_str(), "w");
    std::fprintf(f, "# t = %.8e  (%.4f Myr)\n", t, t * units::TIME_MYR);
    std::fprintf(f, "# nz = %d  nv = %d\n", m.nz, m.nv);
    std::fprintf(f, "# z[kpc] v[km/s] f\n");
    for (int iz = 0; iz < m.nz; ++iz) {
        for (int iv = 0; iv < m.nv; ++iv) {
            std::fprintf(f, "%.6e %.6e %.6e\n", m.z[iz], m.v[iv], m.at(iz, iv));
        }
    }
    std::fclose(f);
}

// ---------------------------------------------------------------------------
// Phase-spiral winding
// ---------------------------------------------------------------------------
// The snail is an m = 1 spiral in the (z, v) plane. Rescale v by the vertical
// frequency so that unperturbed orbits are circles, bin in radius, and take
// the phase of the m = 1 Fourier component of F in azimuth. A spiral shows up
// as that phase advancing monotonically with radius; the total advance across
// the occupied range, in turns, is a direct clock on the time since the kick.
//
// KNOWN LIMITATION. Rescaling by a single frequency only circularises the
// orbits if the vertical potential is harmonic, and if it were harmonic there
// would be no spiral to measure (Kohn's theorem; see README correction 6). In
// a real anharmonic well the equilibrium is elongated in these coordinates and
// that mismatch appears as a spurious m = 2 pattern which dominates the m = 1
// signal, so the numbers this returns are not yet a usable winding rate. The
// plotting script sidesteps the issue entirely by dividing each map by the
// measured pre-kick equilibrium, which assumes nothing. Doing it properly here
// needs action-angle coordinates: radius should be the vertical action and
// azimuth the orbital phase. That is Stage 1 work.
struct WindingProfile {
    std::vector<double> radius;  // in kpc, using v/omega as the vertical axis
    std::vector<double> phase;   // radians, unwrapped
    std::vector<double> power;   // m=1 amplitude relative to m=0
    double turns = 0.0;          // total advance across the profile, in turns
};

inline WindingProfile phase_spiral_winding(const HusimiMap& m, double omega,
                                           int n_r = 40) {
    WindingProfile wp;
    const double rmax = std::min(std::fabs(m.zmax),
                                 std::fabs(m.vmax) / omega);
    std::vector<double> re(n_r, 0.0), im(n_r, 0.0), m0(n_r, 0.0);

    for (int iz = 0; iz < m.nz; ++iz) {
        for (int iv = 0; iv < m.nv; ++iv) {
            const double x = m.z[iz];
            const double y = m.v[iv] / omega;
            const double r = std::hypot(x, y);
            if (r >= rmax || r == 0.0) continue;
            const int ir = static_cast<int>(r / rmax * n_r);
            if (ir < 0 || ir >= n_r) continue;
            const double th = std::atan2(y, x);
            const double val = m.at(iz, iv);
            m0[ir] += val;
            re[ir] += val * std::cos(th);
            im[ir] += val * std::sin(th);
        }
    }

    double prev = 0.0;
    bool first = true;
    for (int ir = 0; ir < n_r; ++ir) {
        if (m0[ir] <= 0.0) continue;
        const double r = (ir + 0.5) * rmax / n_r;
        double ph = std::atan2(im[ir], re[ir]);
        const double amp = std::hypot(re[ir], im[ir]) / m0[ir];
        if (!first) {
            while (ph - prev > units::PI) ph -= 2.0 * units::PI;
            while (ph - prev < -units::PI) ph += 2.0 * units::PI;
        }
        prev = ph;
        first = false;
        wp.radius.push_back(r);
        wp.phase.push_back(ph);
        wp.power.push_back(amp);
    }
    if (wp.phase.size() > 1) {
        wp.turns = (wp.phase.back() - wp.phase.front()) / (2.0 * units::PI);
    }
    return wp;
}
