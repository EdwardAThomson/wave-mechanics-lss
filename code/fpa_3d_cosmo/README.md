# 3D FPA with Cosmological Initial Conditions

Reproduces results from Thomson (2011) Chapter 4, Figures 4.5-4.7.

## Overview

Two approaches to evolving cosmological initial conditions on a 3D grid:

1. **Zel'dovich Approximation (ZA)** — classical particle displacement + CIC gridding.
   Exact in the linear regime, serves as the reference solution.

2. **Free Particle Approximation (FPA)** — wave-mechanical evolution via the
   Schrödinger equation with V=0. Encodes density and velocity in a
   wavefunction via the Madelung transform: `psi = sqrt(1+delta) * exp(-i phi_v / nu)`.
   Recovers ZA-like growth in the semiclassical limit (small nu), then
   develops quantum interference fringes at late times / large D.

## Physics

Initial conditions use the BBKS transfer function (Eqs 3.35-3.38) with
cosmological parameters h=0.71, Omega_m=0.27, sigma_8=0.81, n_s=1.
The Gaussian random field is generated at z_init=50 (growth factor D=1)
and evolved to D=34 (approximately z=0).

### The nu parameter

The parameter nu (analogous to hbar/m) controls the classical-quantum
transition:

- **Small nu**: semiclassical regime, FPA matches ZA. But the Madelung
  phase `phi_v/nu` oscillates rapidly and can alias on the grid.
- **Large nu**: phase is well-resolved but quantum pressure suppresses
  structure growth and interference fringes appear earlier.

The grid-resolvability constraint is: `s_max / nu < k_Nyquist = pi*N`,
where s_max is the maximum displacement. At N=64 this gives nu > 5.4e-5;
at N=128 it gives nu > 3e-5.

At late times (large D), shell crossing produces counter-propagating
wavefunctions that interfere, creating fringe patterns in the density
field. This is a genuine physical prediction of the wave-mechanical
approach, not a numerical artefact.

## Directory structure

```
fpa_3d_cosmo/
├── Makefile
├── README.md
├── src/
│   ├── fpa_3d_cosmo.cpp      # ZA main code (BBKS ICs + CIC gridding)
│   └── fpa_test_nu.cpp       # FPA nu sweep (Madelung + free propagator)
├── scripts/
│   └── plot_fpa_3d_cosmo.py  # ZA diagnostic plots
└── output/
    ├── za/                    # ZA reference results + plots
    └── fpa_N<N>_nu_<val>/     # FPA results per (N, nu) combination
```

## Build and run

```bash
make all                      # builds fpa_3d_cosmo and fpa_test_nu
make run                      # runs ZA, writes to output/za/
make plot                     # generates ZA diagnostic plots in output/za/
./fpa_test_nu <nu> [N]        # FPA sweep, e.g. ./fpa_test_nu 5e-5 128
```

## Key results

| N   | nu   | Phase resolvable | FPA-ZA correlation (D=5) | Notes                        |
|-----|------|------------------|--------------------------|------------------------------|
| 64  | 1e-4 | Yes (0.54)       | 0.81                     | Best at N=64                 |
| 64  | 5e-5 | No (1.08)        | 0.90                     | Aliased but highest corr     |
| 128 | 5e-5 | Yes (0.61)       | 0.74                     | Clean, interference at D=34  |
| 128 | 3e-5 | No (1.01)        | 0.81                     | Marginal aliasing            |

Interference fringes become visible at D~30-34 for nu~5e-5, corresponding
to the onset of shell crossing in the ZA. This is the regime where the
wave-mechanical approach predicts qualitatively different behaviour from
classical dynamics.
