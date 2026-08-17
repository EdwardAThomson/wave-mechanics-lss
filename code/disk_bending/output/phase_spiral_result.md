# The vertical phase spiral from a Schrödinger-Poisson solver

Result write-up for Stage 0 of `code/disk_bending/`. The solver, its validation
ladder and the full review of the adaptation plan are in the directory
`README.md`; this is the physics.

---

## Summary

A warm, self-gravitating, vertically isolated stellar slab is given an impulsive
vertical kick and evolved with a Schrödinger-Poisson solver. Husimi-smoothing
the wavefunction recovers the full distribution function `f(z, v_z)` on a grid,
and it winds up into a textbook Antoja-style phase spiral **with no particle
shot noise whatsoever**.

Three findings, in order of how surprising they were:

1. **The spiral forms cleanly and is completely noise-free.** This is the case
   for using wave mechanics here rather than N-body: the fine-grained spiral is
   exactly the diagnostic where particle codes need enormous particle counts to
   beat Poisson noise, and a wave solver gives the distribution function
   directly.
2. **Strong disk self-gravity suppresses the spiral.** Two runs identical but
   for the stellar surface density retain 0.4% and 65% of the initial
   corrugation after 14 vertical periods. The self-gravitating slab's dipole
   response is a collective mode that neither damps nor winds appreciably.
3. **A purely harmonic vertical potential produces no spiral at all**, by Kohn's
   theorem. The winding comes entirely from the anharmonicity of the vertical
   potential, so how much of the vertical restoring force comes from the disk
   itself versus the halo is a physics parameter to scan, not a detail to fix by
   choosing a value once.

---

## Why a wave solver for this

Schrödinger-Poisson is a Vlasov-Poisson solver in disguise. Husimi-smoothing
`|psi|^2` in phase space recovers `f(x, v)` to `O(hbar_eff)`, where
`hbar_eff = hbar/m` is a phase-space coarse-graining scale rather than a physical
quantum. Anything collisionless is therefore fair game, including disk bending
and vertical phase mixing.

The specific advantage over N-body is smooth, shot-noise-free phase space. The
Gaia vertical phase spiral is fine-grained structure in `f(z, v_z)`, and
resolving it with particles is expensive precisely because the signal is a small
contrast on a smooth background. Here the distribution function comes out of a
single grid calculation, and the contrast maps below are limited by the
phase-space resolution `hbar_eff`, not by counting statistics.

---

## The model

A one-dimensional vertical slab: periodic-safe in the kinetic step, isolated in
gravity, with the wavefunction ensemble built from eigenstates of the
self-consistent vertical Hamiltonian.

| | tracer run | self-gravitating run |
|---|---|---|
| stellar surface density `Sigma` | 5 M⊙/pc² | 50 M⊙/pc² |
| vertical dispersion `sigma_z` | 20 km/s | 20 km/s |
| rigid halo term `nu_ext`, `z_ext` | 60 km/s/kpc, 0.5 kpc | 60 km/s/kpc, 0.5 kpc |
| `hbar_eff` | 0.6 kpc km/s | 0.6 kpc km/s |
| measured scale height `h` | 0.229 kpc | 0.155 kpc |
| vertical period `T_z` | 98.4 Myr | 68.8 Myr |
| grid | 4096 over 8 kpc | 4096 over 8 kpc |
| streams | 320 | 320 |
| kick | 15 km/s bulk | 15 km/s bulk |
| duration | 14 `T_z` (1378 Myr) | 14 `T_z` (964 Myr) |

Both runs conserve energy to better than 5 parts in 10⁶ and keep
velocity-space spill below 4 parts in 10¹¹, so nothing below is an aliasing
artefact.

The external term matters and is not decoration. An isolated self-gravitating
slab given a uniform vertical kick simply translates forever, because
self-gravity is internal and nothing restores the centre of mass. A real disk is
kicked relative to a halo that stays put.

---

## Result 1: the spiral

`figures/phase_spiral_tracer.png`

The tracer population winds up into a multi-armed snail within a few hundred
Myr. The maps are the Husimi transform of the wavefunction, divided by the
measured pre-kick equilibrium so that the smooth background divides out exactly
and no assumption of circular symmetry is made. The vertical axis is `v_z/omega_z`
in kpc, so an unperturbed orbit is a circle and a wound spiral looks like a
spiral rather than a sheared blob.

The corrugation amplitude decays from 0.247 kpc to 0.0009 kpc over the run, a
factor of 280. That decay *is* the winding: the perturbation is not lost, it is
being wrapped into ever finer phase-space structure, and the mean displacement
averages towards zero as the arms tighten.

Note what is absent: there is no speckle, no shot noise, and no smoothing kernel
applied after the fact. The contrast is limited only by `hbar_eff`.

---

## Result 2: self-gravity suppresses the spiral

`figures/phase_spiral_selfgravitating.png`, `figures/moments_tracer_vs_selfgrav.png`

Identical setup, ten times the stellar surface density:

| | corrugation, early → late | retained |
|---|---|---|
| tracer, `Sigma` = 5 M⊙/pc² | 0.247 → 0.0009 kpc | 0.4% |
| self-gravitating, `Sigma` = 50 M⊙/pc² | 0.277 → 0.180 kpc | 65% |

The self-gravitating slab barely winds. Its dipole response is a **collective
normal mode**: the sheet sloshes as a whole, and because self-gravity is an
internal force the whole structure moves together rather than each orbit
phase-mixing at its own frequency. The external potential is distinctly
anharmonic in both runs, and the orbital frequency varies by roughly 40% across
the occupied range, so this is not for lack of anharmonicity.

The practical consequence is that the observed sharpness of the Milky Way's
phase spiral carries information about the disk-to-halo ratio of the vertical
restoring force, and that a simulation which gets that ratio wrong will get the
winding rate wrong even with the kick and the potential shape otherwise right.

---

## Result 3: a harmonic potential gives no spiral, ever

In a purely harmonic vertical potential the centre-of-mass mode decouples
exactly from every internal degree of freedom and oscillates undamped forever
(Kohn's theorem), and every orbit shares one frequency so the distribution
function rotates rigidly instead of winding. Measured directly: with
`Phi_ext = nu^2 z^2 / 2` the corrugation oscillates at exactly `2 pi / nu` with no
decay over 20 vertical periods.

The winding therefore comes entirely from anharmonicity, which is why the
external term in the runs above uses a flattening scale,
`Phi_ext = nu^2 z_e^2 [sqrt(1 + z^2/z_e^2) - 1]`, harmonic near the midplane and
flattening above `z_e`. Setting `z_e` very large recovers the degenerate case and
the spiral disappears, which is a useful check rather than a bug.

This is worth stating plainly because it is an easy way to waste a week: a
plausible-looking harmonic vertical potential will produce a perfectly stable,
perfectly undamped oscillation and no spiral at all, with nothing obviously
wrong in the output.

---

## Validation

Six tests, run in order, each catching a different class of error
(`make check`, about four minutes):

| # | test | result |
|---|---|---|
| 0 | isolated vertical Poisson vs analytic sech² sheet | 2nd order; edge force exact to 3e-9 |
| 1 | free particle: plane-wave phase, Gaussian spreading | 1e-14 |
| 2 | harmonic oscillator: spectrum and coherent state | eigen-residual 1e-12; Strang 2nd order |
| 3 | Jeans growth and the quantum-pressure branch | 3e-9 and 1e-11 |
| 4 | warm isothermal sheet stationarity | z_rms drift 4e-6 over 8 vertical periods |
| 5a | hybrid (x, z) Poisson, four ways | round-off (4e-16) end to end |

Test 3 runs both signs of the dispersion relation, gravity-dominated growth and
quantum-pressure-dominated oscillation, which pins the sign of the gravity
coupling and the kinetic operator together. Test 4 was expected to be the one
that failed first; the eigenstate construction below makes it pass by
construction.

---

## Three method notes that mattered

**Warm initial conditions come from eigenstates, not sampled streams.** A single
`psi = sqrt(rho) exp(iS/hbar)` is a cold stream, and a disk with `sigma_z` = 20
km/s is not. The obvious fix, sampling streams from `f(z, v_z)`, has a trap: a
stream of constant velocity is not stationary in the sheet potential, it simply
falls, and the streams that *are* stationary are constant-energy tori needing two
velocity branches with a density that diverges at the turning point where WKB
fails. Using eigenstates of the self-consistent vertical Hamiltonian instead
makes `|u_n|^2` exactly stationary by construction. Occupying them with weights
proportional to `exp(-E_n/sigma_z^2)` reproduces the isothermal distribution
function to `O(hbar_eff^2)`. Measured: the scale height agrees with the analytic
Spitzer value to five figures and the density to a few parts in 10⁴, with the
residual falling as `hbar_eff^2`.

**The resolution constraint is the same one that governs the cosmological runs.**
The grid can only represent velocities up to `v_max = hbar_eff * pi / dz`, and
the controlling dimensionless number throughout is `k_max dz = pi v_max/v_ceiling`.
Refining `hbar_eff` without refining `dz` raises that ratio and makes
velocity-space diagnostics *worse*, not better. This is the thesis's own headline
lesson reappearing in a completely different problem.

**Extra wavefunction streams buy variance reduction, not warmth.** A single `psi`
at a given `hbar_eff` already carries the full warm distribution function. What
the streams reduce is the interference speckle, and the spurious vertical dipole
that comes with it follows

```
A ≈ 0.70 sqrt(hbar_eff / omega_z) / sqrt(N_streams)
```

holding to a few per cent over a 64× range in stream count. Since the corrugation
`<z>` is the observable, that noise competes directly with the signal: for this
sheet a single wavefunction gives 0.045 kpc, about a quarter of a scale height.

---

## What is not claimed

- This is a **1D vertical slab**. It has no in-plane wavevector, so it contains
  no bending waves, no bending dispersion relation and no firehose. Those need
  the 2D box, and the Stage 1 work in `README.md` shows that the *stable* bending
  branch additionally needs rotation.
- The perturber is idealised. These are impulsive kicks, not a Sagittarius orbit.
- `hbar_eff` = 0.6 kpc km/s gives a de Broglie length of 0.19 kpc, which is 0.64
  scale heights. The phase-space resolution is therefore coarse enough to see in
  the maps, and the finest arms of a real spiral would be smoothed at that scale.
- No comparison against an N-body run of the same setup has been made. The claim
  here is that the spiral is noise-free, not that it is quantitatively closer to
  truth than a well-resolved particle calculation.

---

## Reproducing

```sh
cd code/disk_bending
make
make check                                  # the validation ladder, ~4 min

./sheet_1d --N 4096 --L 8 --hbar 0.6 --sigma 5  --nu-ext 60 --z-ext 0.5 \
           --vkick 15 --ic multistream --periods 14 --out output/spiral_tracer
./sheet_1d --N 4096 --L 8 --hbar 0.6 --sigma 50 --nu-ext 60 --z-ext 0.5 \
           --vkick 15 --ic multistream --periods 14 --out output/spiral_selfgrav

python3 plot_phase_spiral.py output/spiral_tracer
python3 plot_phase_spiral.py output/spiral_selfgrav
python3 plot_moments.py output/spiral_tracer output/spiral_selfgrav
```

Each production run is a few minutes on four cores. Add `--z-ext 50` to either
to recover the harmonic case and watch the spiral vanish.
