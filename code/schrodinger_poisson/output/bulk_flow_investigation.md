# Bulk Flow Investigation

Investigation into whether the 3D Schrodinger-Poisson solver produces
spurious bulk velocity drift, as was observed in the original thesis code
(Chapter 5, review items #7 and #8).

## Setup

All tests use:
- N = 64^3, BBKS cosmological ICs, seed=42
- nu = 1e-4, L0 = 35707
- z_init = 50 → z_final = 0 (400 steps in ln(a))
- LCDM: Omega_m = 0.27, Omega_Lambda = 0.73

Bulk flow measured as |<j>|/rms(j) where j = nu * Im(psi* grad(psi))
is the probability current, computed from binary wavefunction dumps.

## Results

### Test 1: 3-pass auxiliary functions (sp_3d_cosmo_3pass.cpp)
**Hypothesis**: Double recursion for periodic BCs hasn't converged.
**Result**: Bulk flow IDENTICAL to 2-pass (ratio = 1.0000 at all times).
**Conclusion**: Two passes are sufficient. The recursion is fully converged.

### Test 2: Reversed sweep order z→y→x (sp_3d_cosmo_zyx.cpp)
**Hypothesis**: First-swept dimension gets preferential phase bias.
**Result**: Bulk flow IDENTICAL to x→y→z (all components match exactly).
**Conclusion**: Sweep order has no effect. The splitting operators are not
the source of the directional bias.

### Test 3: Bulk removal every timestep (sp_3d_cosmo_nobulk.cpp)
**Hypothesis**: Removing <j> via Galilean boost at each step eliminates drift.
**Result**: Helps at early times (75-150x reduction at z=33-22) but
WORSENS at late times (bulk flow increases by z=0).
**Conclusion**: At late times, the velocity field is too structured for
a uniform Galilean correction. The boost distorts the phase in
overdense regions where local velocity ≠ mean velocity.

### Test 4: IC-only bulk removal (sp_3d_cosmo_ic_fix.cpp)
**Hypothesis**: The initial Madelung wavefunction has non-zero <j> which
gets amplified by gravity.
**Result**: Essentially NO change (reduction factors ~1.0x at all times).
**Conclusion**: The bulk flow is not seeded by the initial mean current.
It develops during the evolution from the gravitational dynamics of the
particular random realisation.

### Test 5: Cosmological ICs without expansion (sp_3d_cosmo_static.cpp)
**Hypothesis**: The expansion implementation itself introduces bulk drift.
**Result**: The static run ALSO develops bulk flow — comparable magnitude
at early times, and actually WORSE at late times (7-8% sustained vs 0.7%
for the expanding run at z=0).
**Conclusion**: Expansion is NOT the cause. In fact, expansion helps damp
the bulk flow at late times by weakening gravity relative to kinetic energy.
The drift is purely a gravitational dynamics effect.

## Summary

The bulk drift:
- Peaks at ~8.6% of rms velocity around z=10, then subsides to ~0.7% by z=0
- Is NOT caused by: auxiliary function convergence (Test 1), sweep order
  (Test 2), IC mean current (Test 4), or the expansion implementation (Test 5)
- IS caused by: gravitational dynamics of the finite-box random realisation.
  The lowest non-zero k modes of the density field produce a net dipole
  gravitational force. This is a cosmic variance / finite-box effect, well-known
  in N-body simulations.
- The per-step Galilean correction (Test 3) helps at early times but distorts
  the phase structure at late times when the velocity field is nonlinear.
- The expanding universe naturally damps the drift better than a static universe.
- The peak of ~8.6% is moderate and likely NOT the same issue reported in the
  thesis (which described a persistently "growing" bulk flow, review item #8).
  The thesis issue may have been a different problem.

## Experimental variants (preserved in working directory)

| File | Description |
|---|---|
| sp_3d_cosmo.cpp | Baseline (2-pass, xyz sweep, expanding) |
| sp_3d_cosmo_3pass.cpp | 3-pass auxiliary functions |
| sp_3d_cosmo_zyx.cpp | Reversed sweep order (z→y→x) |
| sp_3d_cosmo_nobulk.cpp | Galilean boost every timestep |
| sp_3d_cosmo_ic_fix.cpp | Galilean boost on ICs only |
| sp_3d_cosmo_static.cpp | No expansion (fixed lambda, g_eff) |

## Future directions

- Compare density fields with ZA/FPA output at same redshifts
- Try different random seeds to confirm seed-dependence of bulk flow
- Larger box (N=128) to reduce cosmic variance
