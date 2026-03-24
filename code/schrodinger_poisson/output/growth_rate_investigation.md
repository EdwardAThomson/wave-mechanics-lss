# Growth Rate Investigation

The S-P cosmological simulation shows incorrect growth of density
fluctuations: modes grow ~4-7x at early times then freeze, rather than
following D(a) through to z=0. This document tracks experiments to
diagnose and fix the issue.

## The problem

Mode-by-mode analysis shows all Fourier modes grow ~4-7x in the first
~100 steps then saturate. Linear theory predicts ~39x total growth from
z=50 to z=0. The growth rate is too fast initially, then stops entirely.

## Results

### Test 1: Small L0 (sp_3d_cosmo_smallL0.cpp)
L0 = 357 (nu = 0.01). Phase argument ~0.004 per step.
**Result**: k=1 mode tracks D(a) correctly through the linear regime,
then overshoots at late times (62x vs 39x) — nonlinear collapse as
expected. k=3 and k=5 suppressed at early times by quantum pressure
(Jeans effect at this large nu), then grow via nonlinear mode coupling.
**Conclusion**: Growth rate IS correctly calibrated. The dimensionless
equations reproduce the standard cosmological growth equation.

### Test 2: Sub-stepped potential (sp_3d_cosmo_substep.cpp)
L0 = 35707 (nu = 1e-4), 100 sub-steps per potential step.
**Result**: Growth IDENTICAL to baseline — still freezes at ~6x.
**Key insight**: The Cayley phase argument at early times is only ~0.38,
which the Cayley approximation handles at 98.8% accuracy. Sub-stepping
improves this to 99.9999% but makes no difference because Cayley accuracy
was never the problem.

### Revised diagnosis

The Cayley approximation is NOT the issue. The growth freeze at nu=1e-4
must be caused by something else. The small L0 run "worked" because it
changed the physics (nu=0.01 gives very different quantum pressure
dynamics), not because it fixed a numerical issue.

Possible true causes:
1. **Kinetic step phase resolution**: At nu=1e-4, the Madelung phase
   phi_v/nu has gradients that approach the Nyquist frequency. The
   Goldberg finite-difference Laplacian may not correctly propagate
   these high-frequency phases, causing the kinetic step to effectively
   "freeze" the wavefunction.
2. **Lambda too large**: With lambda ~ 500-7000, the kinetic step makes
   very small changes per step (alpha = 2 - i*lambda ≈ -i*lambda). The
   tridiagonal recursion might lose precision when lambda >> 1.
3. **Strang splitting error**: The splitting K/2-V-K/2 with self-consistent
   V might have systematic errors when the kinetic and potential scales
   are very different (large L0 = large scale separation).

## Next steps

- Test with intermediate L0 values to find where the transition occurs
- Check whether lambda being very large causes precision loss in the
  Goldberg recursion (e.g., test with artificially smaller lambda)
- Try spectral (FFT) kinetic step instead of Goldberg finite-difference
  as a cross-check

## Experimental variants

| File | L0 | nu | Description | Result |
|---|---|---|---|---|
| sp_3d_cosmo.cpp | 35707 | 1e-4 | Baseline | Growth freezes at ~6x |
| sp_3d_cosmo_smallL0.cpp | 357 | 0.01 | Small L0 | Correct growth ✓ |
| sp_3d_cosmo_substep.cpp | 35707 | 1e-4 | 100 sub-steps | Same as baseline ✗ |
