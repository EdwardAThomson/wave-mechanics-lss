# Schrödinger-Poisson Solver: Project Summary

## What we built

A 3D Schrödinger-Poisson cosmological simulation code, rewritten from
scratch in modern C++. The code evolves a complex wavefunction ψ on a
periodic grid, where |ψ|² gives the dark matter density and the phase
gradient gives the velocity field. Gravity is self-consistent via a
Poisson solve at every timestep.

### Code features

- **Spectral kinetic step**: exact free-particle propagator via 3D FFT
  (replaced the original Goldberg finite-difference scheme)
- **Cayley potential step**: exactly unitary (mass conserving to ~10⁻¹³)
- **Strang splitting**: K/2 - V - K/2 operator splitting
- **LCDM expansion**: time-stepping in ln(a), scale-factor-dependent
  Hamiltonian coefficients (thesis Eq 5.55)
- **BBKS initial conditions**: Gaussian random field with Madelung
  wavefunction construction
- **GADGET-4 IC support**: reads particle positions via TSC deposit,
  derives velocity potential from density via Poisson equation
- **Checkpoint/resume**: binary wavefunction dumps for chunked runs
- **Velocity diagnostics**: both probability current and phase gradient
  methods
- **Runtime**: 3 minutes at N=128, ~55 minutes at N=256 (on a laptop)

### Evolution from the thesis code

| Aspect | Thesis (2011) | This rewrite (2026) |
|---|---|---|
| Language | Fortran | C++ |
| Kinetic step | Goldberg finite-difference | Spectral FFT |
| Periodic BCs | Double-recursion auxiliary functions | Inherent in FFT |
| Unitarity | Exact (Cayley) | Exact (Cayley + spectral) |
| Poisson solve | FFT with discrete Green's function | FFT with continuum k⁻² |
| Typical run | Hours on a cluster | Minutes on a laptop |

## The nu problem: a 15-year-old mystery solved

### Background

The original thesis code (2011) produced density fields that were
"messier" than GADGET-2, with rapid gridpoint-to-gridpoint variability.
The velocity field was particularly problematic. The thesis attributed
this to quantum interference effects and attempted to suppress them by
making ν (= ℏ/m) smaller, reasoning that a more "classical" limit would
reduce quantum noise.

### What we found

The relationship is **backwards**. Making ν smaller at fixed N makes the
problem catastrophically worse.

The wavefunction encodes velocity in its phase: v = ν∇(phase). The
Madelung construction exp(-iφ_v/ν) with φ_v/ν ~ 100 radians generates
Fourier harmonics far beyond the grid's Nyquist frequency. These alias
back as noise, corrupting the velocity field before the simulation even
starts.

The constraint is: **v_max = ν × π × N**. To simulate at smaller ν
(more classical), you need proportionally larger N.

### How we found it

Systematic elimination of hypotheses, each tested with a separate
code variant:

| Test | Hypothesis | Result |
|---|---|---|
| Sub-stepped Cayley | Time-stepping accuracy | No change ✗ |
| Spectral FFT kinetic | Finite-difference dispersion | Same freeze ✗ |
| Spectral FFT kinetic | Auxiliary function wrapping | Same freeze ✗ |
| Zero-velocity ICs | IC aliasing only | Also freezes ✗ |
| Small L0 (ν=0.01) | Different ν regime | Correct growth ✓ |
| Round-trip IC test | Velocity encoding fidelity | 8% loss at ν=10⁻⁴ ✓ |

The round-trip test was definitive: construct ψ from (δ, φ_v), then
immediately extract velocity back. At ν=10⁻⁴, 8% of the velocity
amplitude is lost to aliasing. At ν=0.01, <1% is lost.

### The thesis mistake

The thesis used ν = 10⁻⁷ to 10⁻⁸ at N=64. At these values,
φ_v/ν ~ 10⁵ to 10⁶ radians — essentially 100% of the Madelung phase
is aliased. The noise was not quantum interference; it was numerical
aliasing from an unresolvable phase.

### Resolution

| N | Best ν | Phase budget | Growth | Structure |
|---|---|---|---|---|
| 64 | 0.01 | 0.5% | Correct | Fuzzy DM regime |
| 128 | 10⁻³ | 3% | Correct | Cosmic web visible |
| 256 | 10⁻⁴ | 12% | ~25x of 39x | Good cosmic web |
| 512 | 3×10⁻⁵ | 23% | — | Not yet run |

## Comparison with GADGET-4

### Setup

Both codes evolved 128³ dark-matter-only LCDM simulations with matched
initial conditions (GADGET-4 NGENIC ICs, TSC-deposited onto grid for
S-P, Poisson-derived velocity potential for Madelung construction).

### Results

**Power spectrum P(k)**:

- z=19: good agreement at low k, S-P has slight excess at high k
- z<10: S-P develops progressively more small-scale power
- After TSC deconvolution, GADGET P(k) is boosted at high k, but the
  S-P excess remains

**Density PDF**:

- z=19: nearly identical Gaussian distributions
- z<5: S-P has broader high-density tail (sharper peaks) and sharper
  low-density cutoff (emptier voids)
- This reproduces the thesis finding that wave-mechanics gives "more
  lows and more highs" than N-body

**Point-by-point correlation** (matched ICs):

- z=19: r = 0.49 (large-scale structures match)
- z<5: r ~ 0.1 (small-scale disagreement dominates)

**rms(δ) growth**:

- Both track together until z~10
- S-P saturates at rms~6.5 (phase budget exhausted at ν=3×10⁻⁴, N=128)
- GADGET continues growing to rms~21 at z=0

### Physical differences

1. **Phase saturation**: the S-P wavefunction cannot represent
   velocities beyond v_max = ν×π×N. Kinetic energy that exceeds this
   is redistributed as density structure rather than bulk flow.

2. **Quantum pressure**: the Bohm potential evacuates voids faster
   than classical gravity. This is physical for a given ν but absent
   in GADGET.

3. **Force resolution**: GADGET uses adaptive timestepping and tree
   gravity to resolve forces at scales far below the grid spacing.
   The S-P code resolves only to the grid scale.

## Auxiliary investigations

### Bulk flow (5 experiments)

Tested whether the solver produces spurious bulk velocity drift:

- 3-pass auxiliary functions: no effect
- Reversed sweep order (z,y,x): no effect
- Per-step Galilean correction: helps early, worsens late
- IC-only correction: no effect
- Static universe: also has drift, worse than expanding

Conclusion: the bulk flow is cosmic variance from the finite-box
random realisation, not a code artefact. Expansion naturally damps it.

### Tophat collapse

Validated gravity with and without expansion:

- Static: collapse on the expected freefall timescale (t_ff ≈ 0.049)
- Expanding (EdS): gradual overdensity growth as expected
- Spherical symmetry maintained despite split-operator treatment
- Mass conserved to ~10⁻¹³

### Spectral vs Goldberg

The spectral (FFT) kinetic step is recommended over Goldberg:

- Exact dispersion at all k (no finite-difference error)
- Naturally periodic (no auxiliary function wrapping needed)
- 14x faster (16s vs 220s for 400 steps at N=128)
- Simpler code

The Goldberg method with double-recursion periodic BCs was an elegant
solution within the finite-difference framework, but the spectral method
is superior for periodic boxes.

## Files

### Core solvers

- `sp_1d.cpp` — 1D solver (Goldberg + periodic BCs)
- `sp_3d.cpp` — 3D solver (Goldberg, tophat tests, expansion)
- `sp_3d_cosmo.cpp` — Cosmological solver (Goldberg, BBKS ICs)
- `sp_3d_cosmo_spectral.cpp` — Spectral kinetic step
- `sp_3d_cosmo_production.cpp` — Production code (ν=0.01, N=64)
- `sp_3d_cosmo_n128.cpp` — N=128, ν=10⁻³
- `sp_3d_cosmo_n256_test.cpp` — N=256, ν=10⁻⁴, checkpoint/resume
- `sp_3d_cosmo_gadget_ic.cpp` — Reads GADGET ICs for matched comparison

### Analysis

- `compare_za.py` — ZA vs S-P comparison
- `compare_statistics.py` — P(k), PDF, cross-correlation, rms growth
- `extract_gadget_ics.py` — TSC deposit + Poisson velocity from HDF5

### Investigation documents

- `output/bulk_flow_investigation.md`
- `output/growth_rate_investigation.md`
- `output/investigation_narrative.md`
- `output/sp_cosmo_gadget_ic_nu_0.000300/statistics/comparison_report.md`

## Next steps

- Run N=256 with GADGET-matched ICs for a direct comparison at higher N
- Try N=512 (18 GB RAM, ~22 min) for even more classical regime
- Velocity comparison once GADGET velocity units are resolved
- Power spectrum ratio P_SP/P_G4 as 2D function of (k, z)
- Consider adaptive mesh refinement (as in Schive et al. 2014)
