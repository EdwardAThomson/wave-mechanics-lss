# The stable bending branch from a rotating Schrödinger-Poisson slab

Result write-up for the rotation stage of `code/disk_bending/`. The
feasibility analysis that preceded it is `rotation_spike.md`; the formulation
and validation live in `src/rotation.h` and `tests/test_epicyclic.cpp`. This
is the physics.

---

## Summary

Rotation is added to the 2D (x, z) slab as a Landau-gauge vector potential:
one wavefunction stream per guiding centre, Coriolis reduced to a static
per-stream harmonic trap. With the sheet Toomre-stabilised at the physical
solar-neighbourhood parameters, the stable bending branch, which Stage 1
proved unreachable without rotation, is measured cleanly at the first
attempt:

**omega = 19.48 +- 0.03 km/s/kpc at k = 0.393 /kpc** (three random-phase
seeds; 0.15% spread). The same measurement without rotation gave 42.98,
28.23, 1.99 across three seeds: pure noise.

Three findings:

1. **Rotation removes both Stage 1 blockers at once, as predicted.** The
   fragmentation blocker and the noise blocker had a common cause (the box
   size limit), so Toomre support at Q = 2.45 simultaneously opened the
   k << k_J regime and dropped the corrugation noise floor to 1-2% of the
   signal.
2. **Two k points, and the free-streaming pressure term is decisively
   suppressed.** At mode 2 (k = 0.785, one seed) omega = 28.51 against
   gravity-only 32.58 and gravity-minus-full-pressure 12.97: the
   full-pressure candidate is excluded by a factor 4.8 in omega^2. Across
   both points omega^2 sits at 72-77% of razor-thin gravity-only, and the
   implied pressure strength falls steeply with k a_epi, the direction
   epicyclic averaging predicts.
3. **Razor-thin Q is not the stability boundary of a thick sheet.** A
   rotating sheet with sigma_z = 20 km/s refuses to fragment even at
   razor-thin Q = 0.49, because the would-be unstable wavelengths sit at
   k h ~ 1-3 where thickness dilutes in-plane self-gravity by roughly
   1/(1 + k h). Thinning the sheet to sigma_z = 10, which drives the
   corrected relation negative, brings fragmentation straight back at 11.8
   e-folds per time unit. All four stability verdicts land as predicted.

---

## The formulation

In a frame rotating rigidly at Omega, Coriolis enters the Schrödinger
operator as a vector potential with uniform field kappa = 2 Omega. In Landau
gauge A = (0, kappa x, 0), with nothing depending on y, p_y = hbar k_y is
conserved and the problem separates into independent sectors:

```
H_ky = -(hbar^2/2)(d2/dx2 + d2/dz2) + (hbar k_y - kappa x)^2 / 2 + Phi(x, z)
```

The Coriolis force is a local harmonic trap of frequency kappa about the
guiding centre x_g = hbar k_y / kappa. One stream per guiding centre; the
trap is one extra diagonal term in the existing potential step; nothing else
in the solver changes. For rigid rotation the background radial gravity
balances the centrifugal force, so the local box keeps only the trap.

Warmth in the plane comes from thermally occupying the epicyclic (Landau)
levels about each guiding centre, the exact in-plane analogue of the vertical
eigenstate construction: E_n = hbar kappa (n + 1/2), w_n ~ exp(-E_n /
sigma_x^2). Each level carries <v_x^2> = <v_y^2> = E_n, so the occupation
reproduces an isotropic in-plane Schwarzschild distribution, which is the
correct axis ratio sigma_y / sigma_x = kappa / 2 Omega = 1 for a rigidly
rotating box.

Two implementation points that improved on the spike:

- **The trap is wrapped periodically**, 0.5 kappa^2 wrap(x - x_g)^2, and the
  level basis has exactly compact support, so the box stays fully translation
  invariant and the spike's 2.7 kpc edge buffer is unnecessary. The trap cusp
  at each stream's antipode sits where that stream's amplitude is ~e^-39.
- **Flux quantisation never appears.** In the reduced (x, z) problem k_y is a
  per-stream parameter, so guiding centres can be placed freely; the
  N_phi ~ 3001 counting in the spike only constrains a notional Ly.

The velocity ceiling v = hbar pi / dx, the repo's headline constraint,
reappears as the level truncation: levels are kept while their classical
turnaround stays under 0.9 x ceiling. At the production grid that truncates
the Maxwellian at 2.7 sigma, drops 2.5% of the tail weight and realises
sigma_x = 38.05 for a target of 40 (-4.9%). All quoted predictions use the
realised value.

## The model

| | value |
|---|---|
| stellar surface density Sigma | 50 M⊙/pc² |
| dispersions sigma_x, sigma_z | 40, 20 km/s |
| epicyclic frequency kappa | 44.2 km/s/kpc (Omega_box = 22.1) |
| Toomre Q | 2.45 |
| hbar_eff | 0.6 kpc km/s |
| grid | 1024 x 320 over 16 x 5 kpc (dx = dz = 1/64 kpc) |
| streams | 75 guiding centres, 222 epicyclic levels each |
| vertical library | 130 states to 5.9 sigma_z^2 (band-edge capped) |
| scale height h | 0.295 kpc |
| excitation | pure displacement h0 = 0.05 kpc = 17% of h, one mode |
| run length | 8 periods of the slower candidate |

Kinetic operator: MatchedFD in z (the vertical eigensolve is the same
three-point operator, so the library is exactly stationary), spectral in x
(Hermite functions are continuum eigenstates and every occupied level keeps a
margin below the Nyquist velocity, so the stack is stationary to the aliasing
of an exp(-large) tail; matched-FD in x would instead distort velocities by
double digits at sigma_x = 40).

## Result: the dispersion relation, two points so far

Mode 1, three random-phase seeds:

| seed | omega (rigid channel) | t=0 bending purity | noise / signal | energy drift |
|---|---|---|---|---|
| 8080 | 19.52 | 0.9997 | 1.9% | 1.2e-6 |
| 9091 | 19.46 | 0.9987 | 1.2% | 1.5e-6 |
| 10102 | 19.46 | 0.9995 | 0.8% | 2.3e-6 |

Frequency resolution of one run is 2.19 km/s/kpc; the three seeds agree to
0.06, so mode 2 was run with a single seed (purity 0.9998, noise 3.9%,
energy drift 2.9e-6). Compare the identical protocol without rotation
(test_mode_decomp, Lx = 3 kpc, the largest box that did not fragment):
42.98, 28.23, 1.99.

Against the razor-thin candidates, using the realised sigma_x = 38.05:

| mode | k [/kpc] | k h | k a_epi | omega measured | gravity only | with pressure | omega^2 / 2 pi G Sigma k |
|---|---|---|---|---|---|---|---|
| 1 | 0.393 | 0.116 | 0.34 | 19.48 +- 0.03 | 23.03 | 17.53 | 0.718 |
| 2 | 0.785 | 0.232 | 0.68 | 28.51 | 32.58 | 12.97 | 0.766 |

Mode 2 is the discriminating one: the candidates differ by a factor 2.5
there, and the measurement excludes the full free-streaming pressure term by
a factor 4.8 in omega^2 while sitting 23% below gravity-only. Reading the
two points together:

- **The free-streaming pressure term is strongly suppressed at both k.**
  This is what epicyclic confinement does: the -k^2 sigma_x^2 term is
  derived for stars streaming through the corrugation, and here stars
  oscillate within a_epi = 0.86 kpc instead.
- **Decomposing the residual 23-28% deficit is not yet possible.** Taking
  the thickness dilution as 1/(1 + k h) leaves an implied pressure strength
  of 0.43 at k a_epi = 0.34 falling to 0.06 at 0.68. The direction matches
  epicyclic averaging (more of the wavelength sampled per orbit, more
  cancellation), but the numbers lean on the approximate thickness form, and
  omega^2 / 2 pi G Sigma k being nearly flat (0.72 vs 0.77) across a
  doubling of k says the naive thickness-plus-pressure decomposition is not
  yet pinned down. Modes 3 and 4 would separate the two, since thickness
  suppression grows with k h while the epicyclic factor saturates.

## Result: Toomre stabilisation, with a thickness lesson

Four cases (`test_toomre`), each discriminating something different. The
fragmentation metrics are the in-plane heating rate d(sigma_x)/dt and the
growth of the surface-density band maximum:

| case | d(sigma_x)/dt | band max amplitude | verdict |
|---|---|---|---|
| rotating, sigma_x = 40, sigma_z = 20 (Q = 2.45) | +0.008 | 9.2e-3 → 7.9e-3 | stable |
| rotating, sigma_x = 8, sigma_z = 20 (razor-thin Q = 0.49, kh ~ 1-3) | +0.070 | 1.5e-2 → 1.6e-2 | stable |
| rotating, sigma_x = 8, sigma_z = 10 (thin, corrected omega^2 < 0) | +11.4 | 3.2e-2 → 8.8e-1 at mode 5, 11.8 e-folds/unit | fragments |
| no rotation, sigma_x = 40, sigma_z = 20 (Q = 0) | +5.9 | 7.1e-3 → 4.1e-1 at mode 1, 8.2 e-folds/unit | fragments |

Row 1 against row 4 is the headline: identical physical parameters, rotation
on versus off, and the heating rate drops by a factor ~760 while the mode 1
amplitude that grows 58x without rotation stays flat. Row 3 fragments
exactly in the band the thickness-corrected relation destabilises (modes
4-8; the fluid estimate peaks near k ~ 5, and mode 5 wins). Row 2 is the
thickness result discussed next.

The first attempts at this test are worth recording because they falsified
two successive explanations. At razor-thin Q = 0.61 and 0.49 the rotating
sheet refused to fragment. The first suspect was the guiding-centre lattice
(2.6 streams per unstable wavelength); a rerun at 352 streams changed
nothing, killing that hypothesis. The real answer is physics: with
sigma_z = 20 the sheet has h = 0.295 kpc, the cold in-plane unstable band
sits at k h ~ 1-3, and with the thickness dilution 1/(1 + k h) the corrected
relation

```
omega^2 = kappa^2 - 2 pi G Sigma k / (1 + k h) + k^2 sigma_x^2
```

is positive at every k for this Sigma, kappa and h, even as sigma_x -> 0.
A rotating sheet this thick cannot fragment at any in-plane temperature:
finite thickness widens the stable window well beyond razor-thin Q, the
in-plane analogue of the firehose suppression measured in test_firehose.
The decisive fragmentation control therefore uses a genuinely thin sheet
(sigma_z = 10, h = 0.074), for which the corrected relation goes properly
negative around k ~ 5.

## Validation

The fast ladder (`make check`) gains test_epicyclic:

| check | result |
|---|---|
| coherent state oscillates at kappa | order 2.00 in dt, amplitude error 8e-7 |
| energy with trap term | bounded oscillation, order 2.00 in dt |
| warm rotating stack, mass | exact to 1e-14 (orthonormality) |
| virial split KE_x / E_trap | 1.001 |
| stationarity with gravity, 2 epicyclic periods | energy 1.3e-6, sigma_x 9e-4, z_rms 1.4e-3 |
| spectral spill | < 1e-8 both axes |

One genuine bug was found on the way and is now guarded: the vertical
library's bound-state cut (V at the box edge) can sit above the
finite-difference kinetic band edge 2 hbar^2 / dz^2, and the band-edge
eigenstates carry checkerboard components that appeared as 1.6e-3 vertical
spill. `truncate_sheet_energy` caps the library at 80% of the band.

## What is not claimed

- **No shear.** A rigidly rotating box gives Toomre support, the stable
  branch and a propagating wave, but pattern winding needs the full
  Goldreich-Lynden-Bell shearing treatment. That remains the next stage.
- **Two k points so far**, the second single-seed. The clean statements are
  the frequencies themselves, their 0.15% seed reproducibility at mode 1,
  and the factor-4.8 exclusion of the full free-streaming pressure term at
  mode 2. The pressure-strength numbers (0.43 and 0.06) are
  interpretations stacked on the approximate 1/(1 + k h) thickness factor,
  and the flatness of omega^2 / 2 pi G Sigma k across the two modes warns
  against taking that decomposition quantitatively yet.
- sigma_x is realised 5% below target by the level truncation at the
  velocity ceiling; predictions use the realised value, but the comparison
  inherits that calibration.

## Reproducing

```sh
cd code/disk_bending
make && make check          # includes test_epicyclic, ~5 min
make measure                # builds the measurement binaries

./test_toomre               # four stability cases, ~40 min
./test_rot_disp             # mode 1, three seeds, ~1 h
./test_rot_disp 2 1         # mode 2, one seed
```
