# Growth Rate Investigation

The S-P cosmological simulation shows density fluctuations that grow
~6-8x then freeze, rather than following D(a) through ~39x of growth
from z=50 to z=0. This document tracks the systematic diagnosis.

## The problem

Mode-by-mode Fourier analysis shows all modes grow ~4-8x in the first
~100 steps then saturate. Linear theory predicts ~39x total growth.

## Experiments

### Test 1: Small L0 / large nu (sp_3d_cosmo_smallL0.cpp)
L0 = 357 (nu = 0.01). Much more "quantum" regime.
**Result**: k=1 mode tracks D(a) correctly through the linear regime,
then overshoots at late times (62x vs 39x) — nonlinear collapse as
expected. k=3 and k=5 suppressed at early times by quantum Jeans
effect, then grow via nonlinear mode coupling.
**Conclusion**: The dimensionless equations and growth calibration ARE
correct. The correct growth rate is recovered when nu is large enough.

### Test 2: Sub-stepped Cayley potential (sp_3d_cosmo_substep.cpp)
L0 = 35707 (nu = 1e-4), 100 sub-steps per potential step.
**Result**: Growth IDENTICAL to baseline — still freezes at ~6x.
The Cayley phase argument (~0.38) was already 98.8% accurate.
**Conclusion**: Cayley accuracy is NOT the problem.

### Test 3: Spectral (FFT) kinetic step (sp_3d_cosmo_spectral.cpp)
L0 = 35707 (nu = 1e-4), exact spectral propagator replaces Goldberg.
No finite-difference dispersion error. Naturally periodic (no auxiliary
function wrapping needed). 14x faster than Goldberg.
**Result**: Slightly better (~8x vs ~6x) but STILL freezes.
**Conclusion**: The Goldberg finite-difference Laplacian, the auxiliary
function wrapping, and the pencil sweep structure are all innocent.
The spectral method is better but does not solve the fundamental issue.

## Diagnosis

The freeze occurs with:
- Goldberg finite-difference kinetic step (6x)
- Spectral/FFT exact kinetic step (8x)
- Sub-stepped Cayley potential (6x)

But NOT with:
- Small L0 / large nu (nu=0.01): correct growth ✓

The distinguishing factor is **nu**, not the numerical method. This
points to a physical rather than numerical limitation.

### Root cause: WKB saturation at small nu

At nu = 1e-4, the wavefunction is deep in the semiclassical (WKB)
regime. Nearly all the physics is encoded in the phase of psi, not
the amplitude. The velocity field is v = nu * grad(phase), so the
phase must carry very high-frequency oscillations to represent even
moderate velocities:

    max|grad(phase)| / k_Nyquist = 0.54   (at t=0, before any growth)

As structure grows, the velocity field strengthens. When the phase
gradient reaches the grid Nyquist frequency, the wavefunction can no
longer represent stronger velocities — the phase "wraps" at the grid
scale and further growth is impossible. The wavefunction saturates.

This is the fundamental resolution limit of the wave-mechanical method:
for a given nu and N, there is a maximum velocity (and therefore a
maximum density contrast) that can be represented. Growth freezes when
this limit is hit.

The relationship: v_max ~ nu * k_Nyquist = nu * pi * N. For nu=1e-4,
N=64: v_max ~ 0.02 (box units). With the initial velocity field already
using 54% of this budget, there is only ~2x headroom before saturation.
The observed ~6-8x growth before freezing is consistent with this.

For nu=0.01, N=64: v_max ~ 2.0 (box units), giving ~400x headroom.
This is why the large-nu run shows correct growth all the way to z=0.

### Implications

This is NOT a code bug. It is the known resolution constraint of the
wave-mechanical approach:

    nu * N ~ v_max / pi

To simulate at smaller nu (more classical), you need proportionally
larger N. The thesis (Chapter 5) discusses this trade-off: smaller nu
gives better semiclassical approximation but requires finer grids.

For N=64 with cosmological ICs at z_init=50:
- nu = 0.01: correct growth, but strong quantum Jeans suppression
  at k > 3 (scales < ~10 Mpc/h). Good for large-scale structure.
- nu = 1e-4: negligible quantum effects, but phase saturates before
  z=0. Would need N ~ 6400 to resolve.
- nu = 1e-3: a possible compromise — moderate quantum effects,
  phase budget allows ~20x growth. Worth testing.

### Spectral vs Goldberg

The spectral kinetic step (Test 3) is recommended over Goldberg:
- Exact dispersion at all k (no 22% error at k = 0.54*k_Nyq)
- Naturally periodic (no auxiliary function wrapping)
- 14x faster (16s vs 220s for 400 steps at N=64)
- Simpler code (no tridiagonal recursion)

The Goldberg method with double-recursion periodic BCs was an elegant
solution within the finite-difference framework, but the spectral
method is superior for periodic boxes. The Goldberg method would still
be preferred for non-periodic boundary conditions.

## Summary table

| File | nu | L0 | Method | k=1 growth | Result |
|---|---|---|---|---|---|
| sp_3d_cosmo.cpp | 1e-4 | 35707 | Goldberg | 6.3x | Frozen ✗ |
| sp_3d_cosmo_substep.cpp | 1e-4 | 35707 | Goldberg + 100 sub-steps | 6.3x | Frozen ✗ |
| sp_3d_cosmo_spectral.cpp | 1e-4 | 35707 | Spectral FFT | 7.9x | Frozen ✗ |
| sp_3d_cosmo_smallL0.cpp | 0.01 | 357 | Goldberg | 62x | Correct ✓ |
