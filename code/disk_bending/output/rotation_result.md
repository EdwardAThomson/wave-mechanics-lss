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
2. **Five k points: the free-streaming pressure term is dead, and Landau
   damping switches on.** The measured branch tracks the thickness-corrected
   relation omega^2 = 2 pi G Sigma k / (1 + k h) at k >= 0.79 /kpc (the
   full-pressure candidate is excluded by a factor 5 in omega^2 at mode 2),
   while the fundamental sits below it, consistent with residual pressure
   where epicyclic suppression is weakest. The damping rate rises from
   gamma = 0.5 at the fundamental to a peak of 6.4 near k a_epi ~ 1.4.
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

## Result: the dispersion relation and its damping

Mode 1 was measured first with three random-phase seeds (periodogram values
19.52, 19.46, 19.46: a 0.15% spread, against 1.99-42.98 scatter in the
identical non-rotating protocol), which justified single-seed runs for the
rest. The authoritative numbers below are damped-cosine fits to the saved
rigid-channel series (`analyze_rot_disp.py`), all at seed 8080, Nx = 1024;
`figures/dispersion_relation.png` is the picture:

| mode | k [/kpc] | k h | k a_epi | omega | gamma | omega^2 / 2 pi G Sigma k | 1/(1+kh) |
|---|---|---|---|---|---|---|---|
| 1 | 0.393 | 0.116 | 0.34 | 19.45 +- 0.01 | 0.48 +- 0.01 | 0.713 | 0.896 |
| 2 | 0.785 | 0.232 | 0.68 | 29.08 +- 0.14 | 3.67 +- 0.13 | 0.797 | 0.812 |
| 3 | 1.178 | 0.348 | 1.01 | 35.21 +- 0.37 | 5.64 +- 0.34 | 0.779 | 0.742 |
| 4 | 1.571 | 0.464 | 1.35 | 38.95 +- 0.65 | 6.37 +- 0.59 | 0.715 | 0.683 |
| 5 | 1.963 | 0.580 | 1.69 | 41.54 +- 0.48 | 4.37 +- 0.45 | 0.650 | 0.633 |

Convergence check: mode 1 rerun at Nx = 2048, which doubles the velocity
ceiling and eliminates the epicyclic level truncation (sigma_x realised
moves from 38.05 to ~40): omega = 19.56 +- 0.01, a 0.5% shift. The
truncation calibration was not biasing the measurement.

A method note recorded because it bit once: the in-run periodogram is
blind to damping and initially reported mode 2 as 28.51 (and a strongly
damped mode as 0.06 before mean subtraction). The damped-cosine fits
supersede the periodogram values wherever they differ.

Reading the table:

- **The free-streaming pressure term is excluded.** At mode 2 the
  gravity-minus-pressure candidate predicts omega = 13.0; measured 29.08,
  a factor 5.0 in omega^2. Modes 3-5 sit where that relation predicts no
  wave at all, oscillating cleanly. Epicyclic confinement kills the term:
  stars oscillate within a_epi = 0.86 kpc instead of streaming through the
  corrugation.
- **Modes 2-5 track the thickness-corrected curve** omega^2 =
  2 pi G Sigma k / (1 + k h), sitting on it to within a few per cent
  (slightly above at modes 2-3). The fundamental sits 20% below its
  thickness value in omega^2, consistent with residual pressure at the
  smallest k a_epi = 0.34, where epicyclic suppression is weakest. That
  crossover, pressure residual at low k a_epi and none above ~0.7, is the
  epicyclic-averaging picture in one figure, subject to the caveat that
  1/(1 + k h) is itself an approximation.
- **Landau damping switches on with k**: gamma rises from 0.48 at the
  fundamental (quality factor ~20) through 3.7 and 5.6 to a peak of 6.4
  near k a_epi ~ 1.4, then turns over by mode 5. A quantitative comparison
  against kinetic theory for a thick rotating sheet is future work; the
  measured curve is the falsifiable input for it.

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
- **Five k points, one seed each past mode 1.** The clean statements are
  the fitted frequencies and damping rates, mode 1's 0.15% seed
  reproducibility, and the factor-5 exclusion of the full free-streaming
  pressure term. The thickness-vs-pressure decomposition leans on the
  approximate 1/(1 + k h) form, so "tracks the thickness-corrected curve"
  is a statement about that form, not about an exact profile calculation.
- The damped-cosine fit frequency is the damped oscillation frequency; at
  the damping peak (mode 4) the undamped frequency would be ~1.3% higher.
  No conclusion changes at that level.
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
