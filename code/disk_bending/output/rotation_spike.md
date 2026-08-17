# Rotation feasibility spike

Stage 1 established that the stable bending branch is unreachable in an
unsheared box: the in-plane Jeans and bending dispersion relations are the same
expression with opposite sign, so wherever bending waves propagate stably the
sheet is fragmenting. Rotation is what breaks that degeneracy in a real disk.

Before committing to implementing it, this spike answers the question that
would have sunk it: **does putting a Coriolis term in a wave code on a periodic
box force an awkward configuration?**

**Answer: no.** Flux quantisation is not a binding constraint, the stream count
is the same order as the existing runs, and the target science regime is
comfortably reachable. The one genuine complication is the box boundary, and it
has a cheap mitigation.

Reproduce the numbers below with the calculation embedded in this file's
history, or re-derive them from the formulae given.

---

## The formulation

In the rotating frame, Coriolis enters the Schrödinger operator as a vector
potential. In Landau gauge `A = (0, kappa x, 0)`, and with nothing depending on
`y` (which is the case for a bending wave with a purely radial wavevector),
`p_y` is conserved and the problem separates into independent sectors:

```
H_ky = -(hbar^2/2)(d2/dx2 + d2/dz2) + (hbar k_y - kappa x)^2 / 2 + Phi(x, z)
```

The Coriolis force becomes a **local extra potential**: a harmonic trap about
the guiding centre `x_g = hbar k_y / kappa`. That is the whole point. No
spectral cross-terms, no new operator structure, and the existing multi-stream
machinery carries it directly, one stream per `k_y`.

A rigidly rotating box has `kappa = 2 Omega`, so to model the solar
neighbourhood's `kappa = sqrt(2) Omega = 44.2 km/s/kpc` without shear, set the
box rotation to `Omega_box = kappa/2 = 22.1 km/s/kpc`.

## 1. Flux quantisation is not binding

A uniform "magnetic field" on a torus requires the flux through it to be an
integer number of flux quanta:

```
N_phi = kappa Lx Ly / (2 pi hbar_eff)   must be an integer
```

| Lx = Ly [kpc] | hbar_eff | N_phi | fractional tweak to fix |
|---|---|---|---|
| 8 | 0.30 | 1500.5 | 3.2e-4 |
| 8 | 0.60 | 750.3 | 3.5e-4 |
| 16 | 0.30 | 6002.1 | 1.8e-5 |
| 16 | 0.60 | 3001.1 | 1.8e-5 |
| 24 | 0.60 | 6752.4 | 5.5e-5 |

`N_phi` runs to thousands, so the allowed values sit ~0.03% apart in any
parameter. Hitting an integer needs a change to `hbar_eff` (or `Ly`) of order
one part in 10⁴, which is far below any physical significance. **This was the
risk that could have killed the approach, and it does not.**

## 2. Stream count is the same order as the existing runs

The Landau degeneracy `N_phi` is large, but the number of streams actually
needed is set by physics rather than by the degeneracy. Guiding centres are
spaced `dx_g = 2 pi hbar / (kappa Ly)`, and each stream's density is smeared
over its epicyclic amplitude `a_epi = sigma_x / kappa`. A uniform sheet only
requires the guiding centres packed comfortably finer than that envelope, so
every m-th one suffices.

At `sigma_x = 40 km/s`, `a_epi = 0.905 kpc`. Packing guiding centres at
`a_epi/8`:

| Lx [kpc] | N_phi | streams needed | memory at 1024 × 512 |
|---|---|---|---|
| 8 | 750 | 71 | 0.60 GB |
| 16 | 3001 | 142 | 1.19 GB |
| 24 | 6752 | 213 | 1.79 GB |

O(100) streams, against the 320 the 1D production runs already use. Memory is
around a gigabyte. **Affordable.**

## 3. The target regime is reachable, and only with rotation

This is the payoff. With `Sigma = 50 Msun/pc^2`, `sigma_z = 20`, `sigma_x = 40
km/s`, the Toomre parameter is

```
Q = kappa sigma_x / (3.36 G Sigma) = 2.45
```

comfortably stable, and matching the real solar neighbourhood. Then:

| Lx [kpc] | k | k h | 2 pi G Sigma k | k² sigma_x² | ratio | sigma_x needed without rotation |
|---|---|---|---|---|---|---|
| 4 | 1.571 | 0.465 | 2122 | 3948 | 0.54 | 29.3 |
| 8 | 0.785 | 0.233 | 1061 | 987 | 1.08 | 41.5 |
| 16 | 0.393 | 0.116 | 531 | 247 | **2.15** | 58.7 |
| 24 | 0.262 | 0.078 | 354 | 110 | **3.23** | 71.8 |

A ratio above one means the self-gravity restoring term dominates, which is the
clean bending regime where `omega^2 -> 2 pi G Sigma k` and the razor-thin
prediction becomes an actual falsifiable comparison. It also has `k h` around
0.1, genuinely thin, so finite-thickness corrections are small rather than
dominant.

The last column is the whole argument in one number. An unsheared box at
`Lx = 16 kpc` would need `sigma_x > 58.7 km/s` to avoid fragmenting, above the
physical 40. That is precisely why Stage 1 could not reach this regime. With
rotation at `Q = 2.4` the sheet holds together at the physical dispersion and
the window opens.

**It also fixes the noise floor for free.** The corrugation mode-noise floor
falls as `sqrt(Lx / lambda_dB)`, and the reason Stage 1 was stuck at `Lx = 3 kpc`
was fragmentation. Going to 16 kpc is a 2.3× improvement in signal to noise on
exactly the measurement that was failing. The physics constraint and the noise
constraint had a common cause, and rotation removes both.

## 4. The one real catch, and its mitigation

The trap `(hbar k_y - kappa x)^2 / 2` is **not periodic in x**. Rigorously,
periodic boundaries with a uniform field require magnetic translation operators
and theta-function Landau states rather than simple Gaussians.

In practice each stream is exponentially localised within `a_epi = 0.91 kpc` of
its guiding centre. In a 16 kpc box an interior stream has relative amplitude of
order `exp(-39)` at the far edge, so wrapping is harmless: this is the same
argument that already justifies applying a periodic kinetic step to a
non-periodic slab in `z`.

Only streams whose guiding centre sits near the boundary are genuinely affected.
The cheap honest mitigation is a buffer of about `3 a_epi ≈ 2.7 kpc` at each
edge, measuring modes in the interior. At `Lx = 16 kpc` that leaves 10.6 kpc of
clean interior, still comfortably more than a bending wavelength. The rigorous
torus treatment stays available if the boundary error ever matters.

---

## Verdict and estimated work

Feasible, and worth doing. In rough order:

1. **Epicyclic in-plane basis** replacing plane waves in the initial conditions.
   Harmonic-oscillator states in `x` about each guiding centre, occupied to
   reproduce the Schwarzschild distribution in `v_x`. Roughly 40 epicyclic levels
   are needed at `sigma_x = 40` with `hbar_eff = 0.6`, since the epicyclic quantum
   is `hbar kappa = 26.5 (km/s)^2`.
2. **One extra local potential term per stream** in the evolver. Genuinely small:
   the term is diagonal in real space and slots into the existing potential step.
3. **Edge buffer and an interior-only mode fit** in the diagnostics.
4. Re-run the Stage 1 dispersion measurement at `Lx = 16 kpc`, where the regime
   is clean and the noise floor is 2.3× lower.

What this does **not** cover: shear. A rigidly rotating box gives Toomre support
and therefore the stable bending branch, Landau damping and a propagating wave
for the phase-offset diagnostic. Pattern winding still needs the full
Goldreich–Lynden-Bell shearing-wave treatment, which remains Stage 2.
