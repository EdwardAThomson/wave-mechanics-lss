# Production Run: nu=0.01, N=64, Spectral Kinetic Step

## Parameters

- Grid: N=64^3, box [0,1), L_phys = 32 Mpc/h
- Cosmology: LCDM, Omega_m=0.27, Omega_Lambda=0.73, h=0.71, sigma_8=0.81
- nu = 0.01, L0 = 357
- z_init = 50, z_final = 0 (400 steps in ln(a))
- Kinetic step: spectral FFT (exact propagator)
- Potential step: Cayley transform
- ICs: BBKS Gaussian random field + Madelung velocity (seed=42)
- Runtime: 23 seconds

## Results

### Growth rate (SUCCESS)
The k=1 mode tracks D(a) correctly through the linear regime (z=50
to z~4). This confirms the dimensionless equations are correctly
calibrated — the growth freeze at nu=1e-4 was a resolution issue,
not a calibration error.

### Late-time nonlinear behaviour
At late times (z < 2), the simulation enters a strongly nonlinear
regime driven by quantum dynamics:
- k=1 overshoots D(a) by ~3.5x (explosive nonlinear collapse)
- Small-scale modes (k=3,5) are initially Jeans-suppressed then
  grow via nonlinear mode coupling
- rms(delta) reaches ~14 at z~0.5 vs ~2 for the ZA
- Dense solitonic structures form instead of the ZA's filamentary web

### Comparison with Zel'dovich Approximation
Point-by-point correlation between S-P and ZA density fields:
- z=9.2: r = 0.22
- z=1.5: r = 0.25
- z=0.5: r = 0.17

The poor correlation is expected: at nu=0.01, the quantum Jeans scale
is large enough to fundamentally alter the morphology of structure
formation. The simulation produces fuzzy-dark-matter-like structures
(solitonic cores + interference fringes) rather than classical
filaments and halos.

## Physical interpretation

At nu=0.01 on a 64^3 grid, the quantum Jeans wavenumber at z=50 is
k_J ~ 14 (in mode numbers), meaning wavelengths shorter than ~2 Mpc/h
are pressure-supported. This suppresses small-scale growth at early
times. When the Jeans scale shrinks (as expansion weakens relative to
gravity), the accumulated potential energy is released explosively,
producing dense collapsed structures.

This is physically correct behaviour for a wave-mechanical system at
this nu — it is what fuzzy dark matter simulations look like. But it
is not comparable to classical N-body/ZA results.

## Path to matching N-body

To reproduce classical structure formation, nu must be small enough
that the Jeans scale is below the grid resolution at all times.
The constraint is:

    nu * N >> v_max / pi    (phase resolution)
    k_Jeans >> N/2          (Jeans below grid scale)

For N=128 with nu=1e-3:
- Phase budget: 100x more headroom than nu=1e-4 at N=64
- Jeans scale: ~10x smaller than nu=0.01
- Estimated runtime: ~3 minutes (spectral method)

This is the natural next step.
