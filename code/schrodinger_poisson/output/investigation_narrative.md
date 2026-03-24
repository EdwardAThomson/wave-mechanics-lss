# The Nu Problem: Diagnosing a 15-Year-Old Bug

## Background

In 2011, Edward Thomson's PhD thesis developed a wave-mechanical approach
to cosmological structure formation — solving the Schrödinger-Poisson
system instead of tracking discrete particles. The code produced density
fields that were "messier" than N-body (GADGET-2) results, with rapid
gridpoint-to-gridpoint variability. The thesis attributed this to
interference effects from the single coherent wavefunction and hoped
that using a smaller nu (the quantum parameter hbar/m) would suppress
the noise.

In 2026, we rewrote the code from scratch in modern C++ and
systematically investigated the growth rate problem.

## The symptom

Mode-by-mode Fourier analysis showed density fluctuations growing ~6x
in the first 100 timesteps then freezing, instead of following the
linear growth factor D(a) through ~39x of growth from z=50 to z=0.

## The investigation

### Attempt 1: Cayley approximation accuracy
**Hypothesis**: The Cayley transform (1-ix)/(1+ix) approximation of
exp(-ix) breaks down for large phase arguments.
**Test**: Split each potential step into 100 sub-steps.
**Result**: Identical to baseline. Cayley accuracy was already 98.8%.
**Verdict**: Not the cause.

### Attempt 2: Goldberg finite-difference dispersion
**Hypothesis**: The Goldberg finite-difference Laplacian has 22% error
at high k, causing incorrect phase propagation.
**Test**: Replaced the entire Goldberg kinetic step with an exact
spectral (FFT) propagator — no finite differences, no tridiagonal
recursion, no auxiliary function wrapping.
**Result**: Slightly better (8x vs 6x) but still frozen.
**Verdict**: Not the cause. Also discovered the spectral method is
14x faster.

### Attempt 3: Auxiliary function wrapping (periodic BCs)
The thesis author's original contribution was a double-recursion method
to enforce periodic boundary conditions in the Goldberg scheme. Could
this be introducing artefacts?
**Test**: The spectral FFT method (Attempt 2) is naturally periodic —
no auxiliary functions needed. Same freeze.
**Verdict**: Not the cause. The wrapping method works correctly.

### Attempt 4: Initial condition aliasing
**Hypothesis**: The Madelung wavefunction exp(-i*phi_v/nu) aliases
at small nu, corrupting the initial velocity field.
**Test**: Round-trip test — construct psi from (delta, phi_v), then
extract velocity back from psi.
**Result**: 8% velocity amplitude loss at nu=1e-4; <1% at nu=0.01.
The exponential of a band-limited function generates harmonics beyond
Nyquist: phi_v/nu ~ 100 radians produces harmonics at all multiples
of the original bandwidth, most of which alias.
**Verdict**: Contributing factor, but not the whole story.

### Attempt 5: Zero-velocity initial conditions
**Hypothesis**: If the IC aliasing is the sole cause, starting with
no velocity (psi = sqrt(1+delta), real) should give correct growth.
**Test**: Gravity builds velocity self-consistently from rest.
**Result**: Also freezes at ~5.5x.
**Verdict**: The phase saturation occurs during evolution too, not
just in the ICs.

### The breakthrough

Comparing nu=1e-4 (frozen) with nu=0.01 (correct growth) revealed
the fundamental constraint: **v_max = nu × pi × N**.

The wavefunction's velocity field is v = nu × grad(phase). For the
phase gradient to be representable on a grid of N points, it cannot
exceed the Nyquist frequency pi×N. This sets a hard maximum velocity.

At nu=1e-4 on N=64: v_max = 0.02. The initial velocity field already
uses 54% of this budget, leaving only ~2x headroom before saturation.
Growth freezes when velocities exceed the representable range.

## The thesis mistake

The thesis (2011) used nu = 10^-7 to 10^-8 at N=64. At these values:
- phi_v/nu ~ 10^5 to 10^6 radians
- Phase gradients vastly exceed Nyquist
- Essentially 100% of the Madelung phase is aliased garbage

The thesis noted the resulting noise and concluded: "Although we hope
to limit this effect by using a small nu we may not have completely
suppressed these effects."

**The relationship is backwards.** Smaller nu requires proportionally
larger N. Making nu smaller at fixed N makes the problem catastrophically
worse. The noise attributed to "quantum interference" was actually
**numerical aliasing from an unresolvable Madelung phase**.

## The fix

Choose nu and N so that the phase budget is adequate:
- nu=0.01 at N=64: correct growth, but large quantum Jeans effects
- nu=1e-3 at N=128: correct growth AND recognisable cosmic web

The spectral (FFT) kinetic step is recommended over the Goldberg
finite-difference method: exact dispersion, naturally periodic, 14x
faster, simpler code.

## Validation

At N=128 with nu=1e-3:
- k=1 mode tracks D(a) through the linear regime
- Density slices show filamentary cosmic web at z=0
- Side-by-side comparison with the Zel'dovich approximation shows
  matching large-scale structure, with S-P developing sharper nonlinear
  peaks (from self-consistent gravity)
- Interference fringes visible in voids (quantum signature at this nu)

## Implications

The wave-mechanical approach to structure formation works correctly
when the resolution constraint nu×N >> v_max/pi is satisfied. The
poor results in the original thesis were not a failure of the method
but of the parameter choice. Modern "fuzzy dark matter" codes (Schive
et al. 2014, Mocz et al. 2017) use much higher resolution for this
reason.

The systematic investigation — 8 experimental variants, each ruling
out one hypothesis — would not have been feasible in 2011 with the
original Fortran code. The modern C++ rewrite with the spectral method
runs in seconds, enabling rapid hypothesis testing.
