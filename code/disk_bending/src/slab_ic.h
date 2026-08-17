// Turning a warm eigenstate library into a wavefunction ensemble.
//
// Both of the plan's §3 options are built from the same basis, which is what
// makes them directly comparable:
//
//   Option A (multi-stream): one stream per occupied eigenstate, densities
//     summed incoherently. Cost scales with the number of occupied states, and
//     that number is not free: it is roughly sigma_z * h / hbar_eff, the
//     phase-space volume of the disk measured in hbar_eff cells. So the plan's
//     "memory is free in 1D, re-evaluate in 3D" is really a three-way trade
//     between hbar_eff, grid size, and stream count, not a two-sided one.
//
//   Option B (single psi): one wavefunction, the coherent sum with random
//     relative phases. Worth being precise about what this costs. It is NOT a
//     colder or less faithful distribution function: a single psi at a given
//     hbar_eff already carries the full warm f, because the Husimi transform of
//     the superposition reproduces it. What the extra streams buy is variance
//     reduction on the O(1) interference speckle in |psi|^2, not warmth. So
//     N_streams controls noise, not physics, and in 3D, where matching the DF's
//     phase-space volume would need ~(sigma h / hbar)^3 streams, Option B plus
//     Husimi smoothing is the only viable route, not merely the cheaper one.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "equilibrium.h"
#include "evolve.h"

// Deterministic 64-bit PRNG so runs are reproducible without <random> policy.
struct SplitMix64 {
    uint64_t x;
    explicit SplitMix64(uint64_t seed) : x(seed) {}
    uint64_t next() {
        uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
    double uniform() { return (next() >> 11) * (1.0 / 9007199254740992.0); }
};

// One stream per occupied eigenstate. `max_states` <= 0 keeps every state the
// equilibrium builder produced.
inline SlabState make_multistream(const WarmSheet& s, int max_states = -1) {
    SlabState st;
    st.grid = s.grid;
    st.hbar = s.cfg.hbar;
    const int M = (max_states > 0)
                      ? std::min<int>(max_states, static_cast<int>(s.u.size()))
                      : static_cast<int>(s.u.size());
    st.psi.assign(M, std::vector<Complex>(s.grid.N));
    for (int m = 0; m < M; ++m) {
        const double a = std::sqrt(s.w[m]);
        for (int j = 0; j < s.grid.N; ++j) {
            st.psi[m][j] = Complex(a * s.u[m][j], 0.0);
        }
    }
    return st;
}

// A single wavefunction: the coherent superposition with random phases.
inline SlabState make_single_psi(const WarmSheet& s, uint64_t seed = 12345) {
    SlabState st;
    st.grid = s.grid;
    st.hbar = s.cfg.hbar;
    st.psi.assign(1, std::vector<Complex>(s.grid.N, Complex(0.0, 0.0)));
    SplitMix64 rng(seed);
    const int M = static_cast<int>(s.u.size());
    for (int m = 0; m < M; ++m) {
        const double a = std::sqrt(s.w[m]);
        const double th = 2.0 * units::PI * rng.uniform();
        const Complex c(a * std::cos(th), a * std::sin(th));
        for (int j = 0; j < s.grid.N; ++j) st.psi[0][j] += c * s.u[m][j];
    }
    return st;
}

// Intermediate: `n_streams` independent random-phase superpositions, summed
// incoherently. n_streams = 1 is Option B, n_streams = M with one state each
// is Option A, and anything between trades memory against speckle as 1/sqrt(n).
inline SlabState make_batched(const WarmSheet& s, int n_streams,
                              uint64_t seed = 12345) {
    SlabState st;
    st.grid = s.grid;
    st.hbar = s.cfg.hbar;
    st.psi.assign(n_streams, std::vector<Complex>(s.grid.N, Complex(0.0, 0.0)));
    SplitMix64 rng(seed);
    const int M = static_cast<int>(s.u.size());
    const double norm = 1.0 / std::sqrt(static_cast<double>(n_streams));
    for (int t = 0; t < n_streams; ++t) {
        for (int m = 0; m < M; ++m) {
            const double a = std::sqrt(s.w[m]) * norm;
            const double th = 2.0 * units::PI * rng.uniform();
            const Complex c(a * std::cos(th), a * std::sin(th));
            for (int j = 0; j < s.grid.N; ++j) st.psi[t][j] += c * s.u[m][j];
        }
    }
    return st;
}
