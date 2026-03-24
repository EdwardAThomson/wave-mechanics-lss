# Growth Rate Investigation

The S-P cosmological simulation shows incorrect growth of density
fluctuations: modes grow ~4x at early times then freeze, rather than
following D(a) through to z=0. This document tracks experiments to
diagnose and fix the issue.

## The problem

Mode-by-mode analysis shows all Fourier modes grow ~4-7x in the first
~100 steps then saturate. Linear theory predicts ~39x total growth from
z=50 to z=0. The growth rate is too fast initially, then stops entirely.

## Root cause hypothesis

The Cayley approximation of exp(-i*g_eff*U*dt) breaks down when
g_eff*U*dt >> 1. With L0 = 35707:
- g_eff ranges from 3900 (z=50) to 14500 (z=0)
- U ~ delta_rms ~ 0.01-0.1
- dt = d(ln a) = 0.0098
- Phase argument: g_eff * U * dt ~ 0.4 to 14

For arguments >> 1, the Cayley transform gives the wrong phase
(though still unit magnitude). The potential step scrambles rather
than coherently accumulates phase, freezing the growth.

## Experimental variants

| File | L0 | nu_eff | g_eff(z=50) | Phase arg | Description |
|---|---|---|---|---|---|
| sp_3d_cosmo.cpp | 35707 | 1e-4 | 3897 | ~0.4 | Baseline (broken) |
| sp_3d_cosmo_smallL0.cpp | 357 | 0.01 | 39 | ~0.004 | 100x smaller L0 |
| sp_3d_cosmo_substep.cpp | 35707 | 1e-4 | 3897 | ~0.004 | Sub-step potential (100 sub-steps) |

## Results

### Test 1: Small L0 (sp_3d_cosmo_smallL0.cpp)
L0 = 357, nu_eff = 0.01. Phase argument ~0.004 per step.
Cayley approximation should be very accurate.
Trade-off: more "quantum" behavior (larger Jeans length).
**Status**: TODO

### Test 2: Sub-stepped potential (sp_3d_cosmo_substep.cpp)
Keep L0 = 35707 but split each potential step into 100 sub-steps.
Each sub-step has phase argument 100x smaller. Kinetic step unchanged.
**Status**: TODO

### Test 3: More timesteps (smaller d_lna)
Keep everything the same but use 10x more timesteps (d_lna/10).
Brute force: 4000 steps instead of 400.
**Status**: TODO
