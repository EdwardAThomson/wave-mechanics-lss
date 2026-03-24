# LinkedIn Draft — The Nu Problem

---

15 years ago, my PhD simulated the universe using quantum mechanics instead of particles. The code produced noisy, messy density fields that didn't match N-body simulations. I never figured out why.

Last week, I finally found the smoking gun.

The method solves the Schrödinger-Poisson system: dark matter as a continuous wavefunction instead of discrete particles. The velocity field is encoded in the phase of the wavefunction: v = ν × ∇(phase). That phase must be representable on the computational grid.

Here's the problem: the exponential of a smooth function is NOT smooth. When you construct exp(-iφ/ν) with φ/ν ~ 100 radians, the result has Fourier harmonics far beyond what the grid can resolve. They alias back as noise. The velocity field is corrupted before the simulation even starts.

My thesis tried to fix this by making ν smaller — reasoning that a more "classical" limit would suppress quantum interference. That's exactly backwards. Smaller ν means larger phase values, more harmonics, more aliasing. The noise I attributed to "quantum interference effects" was actually numerical aliasing from an unresolvable phase.

The fix is simple: choose ν and N so that ν × π × N >> v_max.

Finding this required systematically ruling out every other hypothesis:
❌ Cayley time-stepping accuracy — sub-stepped 100x, no change
❌ Finite-difference dispersion — replaced with spectral FFT, same freeze
❌ Periodic boundary wrapping — spectral method needs none, same result
❌ Initial condition construction — zero-velocity ICs also freeze
✅ The fundamental ν×N resolution limit — confirmed across 8 experiments

With the right parameters (ν=10⁻³ at N=128), the code produces recognisable cosmic web structure for the first time: filaments, nodes, and voids, matching the Zel'dovich approximation at early times and developing sharper nonlinear structures through self-consistent gravity.

The modern rewrite runs in 3 minutes on a laptop. The original took days on a cluster. Sometimes you just need fresh eyes, better tools, and the patience to test one thing at a time.

Code: github.com/EdwardAThomson/wave-mechanics-lss

#PhD #Cosmology #Simulation #WaveMechanics #Debugging #Physics #CodeRewrite

---

[Attach: density_slices.png or za_comparison_n128.png]
