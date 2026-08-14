# disk_bending: Schrodinger-Poisson for galactic disk corrugations

Stage 0 of the plan to adapt the `wave-mechanics-lss` solver from cosmological
large-scale structure to vertical bending waves and phase spirals in a galactic
disk: a warm, self-gravitating, vertically isolated slab, kicked impulsively and
diagnosed in phase space.

Everything here is new code. It shares the physics and the split-step idea with
`code/schrodinger_poisson/`, but not the source, for reasons set out in
"What actually carries over" below.

---

## Status

The validation ladder passes end to end (`make check`, about four minutes):

| # | Test | Result |
|---|---|---|
| 0 | Isolated vertical Poisson vs analytic sech² sheet | clean 2nd order; edge force exact to 3e-9 |
| 1 | Free particle: plane-wave phase, Gaussian spreading | 1e-14 |
| 2 | Harmonic oscillator: spectrum and coherent state | eigen-residual 1e-12; ladder error ∝ dz²; Strang 2nd order |
| 3 | Jeans growth and the quantum-pressure branch | 3e-9 and 1e-11 |
| 4 | Warm isothermal sheet stationarity | z_rms drift 4e-6, dE/E 5e-11 over 8 vertical periods |
| 5a | Hybrid (x, z) Poisson, four ways | round-off (4e-16) end to end; 2nd order vs closed form |

Test 3 runs both signs of the dispersion relation (gravity-dominated growth and
quantum-pressure-dominated oscillation), which pins the sign of the gravity
coupling and the kinetic operator together. Test 4 is the one the plan expected
to fail first; the eigenstate construction below makes it pass by construction.

Demonstrated physics: a vertical phase spiral formed from an impulsive kick,
completely free of particle noise (`figures/phase_spiral_tracer.png`), and the
contrasting case where slab self-gravity suppresses it
(`figures/phase_spiral_selfgravitating.png`). Both are 320-stream runs on a 4096
grid, conserving energy to 3e-6 with velocity-space spill below 1e-10.

---

## Stage 1 status (2D, x and z)

Built and validated:

- **Hybrid Poisson solver** (`src/poisson_xz.h`), the §4 option-2 route. FFT in
  x, then an exact O(Nz) recursion per `k_x` against the exponential kernel.
  `test_poisson_xz` certifies it four ways, ending at round-off (4e-16) on a
  separable density, which is what pins the strided batch transform and the FFT
  normalisation.
- **2D evolver and separable warm ICs** (`src/evolve2d.h`, `src/slab_ic_2d.h`).
  The Stage 0 eigenstate library is reused untouched: vertical eigenstates times
  an in-plane random-phase Maxwellian.
- **Column diagnostics and the §6 phase offset** (`src/diagnostics2d.h`),
  driver `sheet_2d.cpp`.

Measured, and not yet clean:

- **Bending-wave dispersion relation** (test 5). Frequencies come out within
  tens of per cent of the razor-thin prediction but scatter non-monotonically
  with `k h` (ratios 0.84, 1.17, 0.86, 0.97, 0.74 for modes 1 to 5), which is
  not physical: finite thickness should reduce the frequency monotonically.
  There were only 5 to 11 zero crossings per mode, and the run heated as it
  went. Needs longer runs and controlled heating before it is a result.
- **Landau damping** (test 6) is not attempted, because it depends on the
  dispersion and the dispersion is being changed by the heating below.

### The finding that gates the rest: granule heating

The in-plane part of any 2D or 3D initial condition can only be a random-phase
superposition, so `|psi|^2` carries interference granules of scale
`lambda_dB = 2 pi hbar_eff / sigma_x` and those granules self-gravitate. The
effective quasi-particle mass is large: at `hbar_eff = 0.6`, `sigma_x = 10 km/s`
it is of order `rho lambda_dB h` ~ 5e6 Msun per unit length, comparable to a
giant molecular cloud, and the relaxation time comes out shorter than a vertical
period.

Measured over two vertical periods, at 12 sigma of velocity headroom in both
directions:

| streams | rms drho/rho | d(sigma_x) per vertical period | max spill |
|---|---|---|---|
| 2 | 0.696 | 17.1 km/s | 0.17 |
| 8 | 0.346 | 15.3 km/s | 0.08 |
| 32 | 0.174 | 15.2 km/s | 0.11 |

Two things to take from this, one solid and one open.

**Solid:** the heating is severe. Starting from `sigma_x = 10 km/s` it roughly
quadruples in two vertical periods. It is much milder at larger `sigma_x` (2.4
km/s per period at `sigma_x = 20`), consistent with a steep dependence, which is
what the quasi-particle-mass argument predicts, since `lambda_dB` and hence the
granule mass grows as `1/sigma_x`.

**Open, and it contradicts the natural guess:** the density contrast falls
exactly as `1/sqrt(N_streams)`, as it must, so the potential fluctuations that
drive the heating should weaken and the heating should fall as `1/N_streams`. It
does not. Over a 16x range in stream count the rate barely moves. Either the
heating is not driven by the shared potential fluctuation in the way that
argument assumes, or these runs are still contaminated: the spill column shows
8 to 17 per cent of the in-plane power above half-Nyquist by the end, because the
heating itself eats the velocity headroom. Both need ruling out before the number
means anything.

**Consequence either way:** granule heating, not resolution, is the binding
constraint on Stage 1, and it must be controlled before the firehose threshold
or the Landau damping rate can be trusted, since both depend on the dispersion
that is drifting. It is also the sharpest warning yet about Stage 2, where a
single wavefunction is the only affordable option.

### A correction to correction 5

Stage 0 ended by predicting that corrugation noise would average down by a
further `sqrt(N_x)` in 2D, because the per-column dipoles are independent. That
was measurable, so it was measured, and it is wrong. The per-column dipoles are
correlated over the in-plane de Broglie scale, so the mode-amplitude noise is
**independent of the grid**: across a 16x range in `N_x` at fixed stream count it
moves by less than a factor of 1.5.

```
Nx      streams   rms <z>(x)   rms mode amplitude
128     1         0.0573       5.5e-3
512     1         0.0541       3.7e-3
2048    1         0.0539       4.1e-3
```

The right factor is `sqrt(L_x / lambda_corr)` with `lambda_corr` of order the
in-plane de Broglie length, so a **wider box** buys signal to noise and a finer
grid buys none. The conclusion Stage 0 was reaching for survives: the mode-noise
floor is a few times 1e-3 kpc even with one stream, comfortably under a realistic
0.1 kpc corrugation, so Option B stays viable at higher dimension.

---

## Review of the plan

The plan is sound in outline and its instincts about where the difficulty lies
were right: the initial conditions and the vertical boundary condition are the
two hard parts, and the dimensional ladder is the correct strategy. What follows
are the corrections and additions that came out of building it.

### 1. What actually carries over from the existing solver is less than §4 claims

§4 lists "Strang split-step core: `exp(-iVΔt/2ħ) → FFT → exp(-iħk²Δt/2) → iFFT`"
under "reuse directly". The main solver does not do that. `sp_1d.cpp` and the
`sp_3d_cosmo*` family use the **Goldberg/Cayley finite-difference** kinetic step
with the thesis's double-recursion trick for periodic boundaries; FFTs appear
only in the Poisson solve. There is exactly one spectral variant,
`sp_3d_cosmo_spectral.cpp`, and that is the only file with the operator §4
describes.

This matters beyond bookkeeping. The Goldberg step is second-order in space and
carries a modified dispersion relation, which is precisely the thing you cannot
afford when the deliverable is a measured dispersion relation (§7 test 5). Plan
to use the spectral operator, and treat the existing code as a reference for
conventions rather than a library to link against. There is also no HDF5 I/O and
no checkpointing anywhere in the repository, so those are new work, not reuse.

### 2. In 1D the isolated vertical Poisson solve needs no Green's function at all

§4 ranks the options as zero-pad ×2 (preferred), hybrid FFT-plus-tridiagonal,
then tall periodic box. For the 1D slab all three are unnecessary:

```
Phi(z) = 2 pi G Integral rho(z') |z - z'| dz'
```

Expanding `|z - z'|` into left and right partial sums evaluates that double sum
*exactly* in two O(N) cumulative sweeps, with no zero-padding, no images, and no
approximation beyond the quadrature rule. That is `poisson_isolated_1d` in
`src/poisson_z.h`, and test 0 confirms clean second-order convergence with the
vertical force at the box edge matching 2πGΣ to 3e-9.

For Stage 1 and beyond, the recommendation should be **inverted**: the hybrid
(FFT in x, then the exact kernel `exp(-|k_x||z-z'|)/2|k_x|` in z, degenerating to
the expression above at `k_x = 0`) is both cheaper in memory and *exact*, whereas
zero-padding is neither. Prefer it.

Test 0 also quantifies why a periodic box is not an option, which is worth having
as a number rather than an argument. The periodic error tracks the image term
`2πGΣz²/L` and falls only as `1/L`:

| box height | error in Φ at z = 1 kpc |
|---|---|
| 4 kpc | 41% |
| 8 kpc | 21% |
| 16 kpc | 10% |
| 32 kpc | 5% |

### 3. The periodicity problem is a *gravity* problem, not a wavefunction problem

Worth stating explicitly because it removes an apparent conflict. The kinetic
operator is local, so applying it spectrally on a periodic FFT grid is legitimate
for a non-periodic slab provided ψ has decayed at the box edges. Only gravity,
being long-ranged, needs the isolated treatment. So the solver runs a periodic
FFT kinetic step alongside an isolated-BC Poisson solve, and that combination is
correct rather than inconsistent. See the header comment in `src/evolve.h`.

The caveat is that "ψ has decayed at the edges" must hold for the **kicked**
orbit, not just the equilibrium. When a kicked sheet swings out far enough to
wrap around, energy conservation collapses (observed: 125% error). The driver
now reports the fraction of mass in the outer fifth of the box and warns.

### 4. Build warm ICs from eigenstates, not sampled streams (§3)

§3 is right that this is the hardest part, and right that a single ψ is a cold
stream. But both proposed options have avoidable problems, and there is a third
that is strictly better in 1D.

The trap in Option A as written: a stream of constant `v_j` is *not* stationary
in the sheet potential, it simply falls. The streams that are stationary are
constant-energy tori, which need two branches `v = ±sqrt(2(E - Φ))` per energy
with density `~1/|v|` diverging at the turning point where WKB fails. That is
the week the plan warns about losing.

**Use the eigenstates of the self-consistent vertical Hamiltonian instead.**
`|u_n|²` is exactly stationary by construction, so test 4 passes to solver
precision rather than being something to debug. Occupying them with
`w_n ∝ exp(-E_n/σ_z²)` and iterating Poisson to self-consistency reproduces the
isothermal distribution function to `O(ħ_eff²)`, the same order at which the
whole Schrodinger-Poisson correspondence holds. This is the standard construction
in the wave-dark-matter halo literature (Yavetz, Li & Hui 2021, and Lin, Chiueh
et al. before it) and it transfers directly to a slab. It is `build_warm_sheet`
in `src/equilibrium.h`.

Measured: the eigenstate sheet reproduces the analytic Spitzer sheet with `h`
correct to five figures and the density to a few parts in 10⁴, and that residual
falls as `ħ_eff²` (4.1× per halving) rather than being a discretisation artefact.

The same library serves both of the plan's IC options on one basis, which is what
makes them comparable at all.

### 5. Option A vs Option B: the measured answer (§3, §8)

§8 calls this "an empirical question worth settling in 2D first". It is settleable
in 1D, and cheaper. Two separate things are going on, and the plan's framing
conflates them.

**Warmth is not what the streams buy.** A single ψ at a given `ħ_eff` already
carries the full warm distribution function; its Husimi transform recovers `f`
correctly. What extra streams buy is variance reduction on the interference
speckle. So §3's "carries genuine shot noise, partially forfeiting the main
advantage over N-body" is too pessimistic about the physics and, as it turns out,
not pessimistic enough about the observable.

**The observable is what decides it.** Speckle in ρ is cosmetic if you only plot
Husimi-smoothed quantities. The spurious *dipole* is not, because `⟨z⟩` is the
bending-wave observable. A random superposition has a random vertical dipole
which oscillates at ω_z forever rather than averaging away. Measured over
ensembles of seeds, with `ℓ_ħ = sqrt(ħ_eff/ω_z)` the coherent-state width:

```
A_spurious  ~  0.70 * sqrt(hbar_eff / omega_z) / sqrt(N_streams)
```

The 1/√N scaling holds to a few per cent over a 64× range in N and a 4× range in
ħ_eff. For the default sheet a single ψ gives `A = 0.045 kpc`, which is 0.25 h and
would swamp a real corrugation; measuring a 0.1 kpc corrugation to 10% needs
about 20 streams.

**The caveat that rescues Option B in 3D:** in 1D the dipole is a single number.
In 2D and 3D the corrugation `⟨z⟩(x, y)` is a field, and fitting a coherent `k_x`
mode across `N_x` columns averages this noise down by a further `sqrt(N_x)`. So
Option B is more viable at higher dimension than the 1D number alone suggests.
That is a Stage 1 measurement, and it is now a cheap one.

**Stream count is not free to choose.** §8's `× N_s` column treats `N_s` as a
knob. It is bounded below by the phase-space volume: representing the DF with one
stream per eigenstate needs `N_s ~ σ_z h / ħ_eff` in 1D, and the 3D analogue
scales as the cube. A 3D warm DF at useful `ħ_eff` would need order 10⁶ streams,
so the "× 64 streams" 3D row in §8 does not describe a warm distribution
function. In 3D, Option B plus Husimi smoothing is the only viable route, not
merely the cheaper one.

### 6. Two things that will silently produce no phase spiral

Neither is in the plan, and both cost real time to diagnose.

**A harmonic external potential gives no spiral, ever.** In a purely harmonic
vertical potential the centre-of-mass mode decouples exactly from every internal
degree of freedom and oscillates undamped forever (Kohn's theorem), and every
orbit shares one frequency so the distribution function rotates rigidly instead
of winding. Observed directly: with `Φ_ext = ½ν²z²` the corrugation oscillates at
exactly `2π/ν` with no decay over 20 vertical periods. The winding comes entirely
from the **anharmonicity** of the vertical potential, so the external term needs a
flattening scale (here `Φ_ext = ν²z_e²[sqrt(1 + z²/z_e²) - 1]`).

**An external potential is also needed at all.** An isolated self-gravitating
slab given a uniform vertical kick just translates forever, because self-gravity
is internal and nothing restores the centre of mass. Real disks are kicked
relative to a halo that stays put.

### 7. Strong slab self-gravity suppresses the spiral

This was the least expected result. Two runs identical but for the disk surface
density, both kicked at 15 km/s and both run 14 vertical periods:

| | corrugation amplitude, early → late | retained |
|---|---|---|
| tracer, Σ = 5 M⊙/pc² | 0.247 → 0.0009 kpc | 0.4% |
| self-gravitating, Σ = 50 M⊙/pc² | 0.277 → 0.180 kpc | 65% |

The tracer population phase-mixes completely and winds into a textbook snail
(`figures/phase_spiral_tracer.png`). The self-gravitating slab does not: its
dipole response stays **coherent**, a collective normal mode that neither damps
nor winds appreciably, even though the external potential is distinctly
anharmonic and the orbital frequency varies by 40% across the occupied range
(`figures/phase_spiral_selfgravitating.png`).

The practical consequence for the project: how much of the vertical restoring
force comes from the disk itself, versus the halo, controls whether a spiral
forms at all. That is a physics parameter to scan, not a detail to fix by
choosing a value once.

### 8. §7's dispersion relation has the wrong dispersion in it

§7 gives `ω² ≈ 2πGΣ|k| − k²σ_z²` with `k_crit ≈ 2πGΣ/σ_z²`. The destabilising
term in the standard razor-thin result (Toomre 1966; Binney & Tremaine §6.6) is
the **in-plane** dispersion, not the vertical one: firehose is driven by the
momentum flux of in-plane streaming along the bent sheet. `σ_z` enters through
finite thickness, and stabilises. Since `σ_R ≈ 2σ_z` in a real disk, using `σ_z`
misplaces `k_crit` by roughly a factor of four. The plan already flags this
formula as needing verification; this is the specific thing to verify.

A direct consequence: **a 1D vertical slab has no in-plane dispersion at all**, so
tests 5 and 6 genuinely cannot be done at Stage 0. §2's "every piece of
bending-wave physics except lateral propagation lives here" is too generous, since
without `k_x` there is no bending mode, no dispersion relation and no firehose.
What Stage 0 genuinely delivers is the phase spiral (which is intrinsically 1D in
`(z, v_z)`), the bending/breathing split as the parity of the kick, and the
resolution trades. That is still the right first step, but for the phase spiral,
which §6 independently calls the strongest single result available.

### 9. Arithmetic in §2

`3750³` complex doubles is 5.27e10 points × 16 bytes = **844 GB**, not 800 TB.
The conclusion is unchanged (with any realistic stream count it is hopeless, and
844 GB for a single ψ is already a large-cluster job), but the figure is off by a
factor of 1000. §8's table is correct.

### 10. The resolution constraint, restated

§1's two-sided constraint is right. Two refinements from measurement:

- The aliasing bound is `Δx ≲ π ħ_eff / v_max` with a spectral operator, so §1's
  `Δx ≲ ħ_eff / v_max` is conservative by π. That π is real headroom and worth
  having.
- It is **not** available with a finite-difference kinetic operator. The FD group
  velocity `ħ sin(k dz)/dz` peaks at `k dz = π/2` and returns to zero at Nyquist,
  so its usable ceiling is `ħ/dz`, lower by exactly π. Test 1(c) shows the
  turnover.

The controlling dimensionless number throughout is `k_max dz = π v_max / v_ceiling`,
the same ratio that governs the LSS runs. Test 4 shows it biting in a way worth
remembering: the density profile is insensitive to it (dominated by well-resolved
low-energy states), but `σ_z`, a kinetic-energy moment weighting the highest-k
states, is off by 1% at `k_max dz = 1.28` and fine at 0.64. **Refining `ħ_eff`
without refining `dz` raises `k_max dz` and makes velocity-space diagnostics
worse, not better.**

### 11. Smaller notes

- The state library must be truncated at the **bound** states. An isolated slab's
  own potential rises only as `2πGΣ|z|`, so a box of height L confines to about
  `πGΣL/σ_z²` in units of `σ_z²` and no deeper; occupying states above that puts
  mass into wall-bouncing orbits and destroys stationarity in a way that looks
  exactly like a solver bug. `build_warm_sheet` now drops them and says so.
- Truncating the library too low biases `σ_z` low, because the missing states
  carry the fast particles. At 4.4σ_z² the deficit is 4.6%; it converges by about
  8σ_z².
- The Strang error in a harmonic well appears as `O(dt⁴)` if measured as the
  displacement error after a whole period, because that is `1 - cos` of the
  `O(dt²)` frequency error. Measure the phase, not the displacement, or you will
  conclude the scheme is fourth order.

---

## Layout

```
Stage 0, 1D vertical slab
  src/units.h          kpc, km/s, Msun; G and time conversions
  src/grid.h           1D vertical grid, FFT wavenumbers
  src/fft.h            RAII wrapper over FFTW plan pairs
  src/poisson_z.h      exact isolated-BC vertical Poisson; analytic Spitzer sheet
  src/eigen1d.h        symmetric tridiagonal eigensolver (bisection + inverse iteration)
  src/equilibrium.h    self-consistent warm sheet from eigenstates
  src/evolve.h         Strang split-step, spectral and matched-FD kinetic operators
  src/slab_ic.h        Option A / Option B / batched wavefunction ensembles
  src/husimi.h         Husimi transform to f(z, v_z); winding profile
  src/diagnostics.h    vertical moments, energies, edge-mass guard
  sheet_1d.cpp         driver

Stage 1, 2D (x, z) box
  src/grid2d.h         2D grid, periodic in x and isolated in z
  src/poisson_xz.h     hybrid Poisson: FFT in x, O(Nz) exponential-kernel recursion
  src/evolve2d.h       2D split-step; bending kicks and column displacement
  src/slab_ic_2d.h     separable warm IC: vertical eigenstates x in-plane Maxwellian
  src/diagnostics2d.h  column moments Sigma(x), <z>(x), <v_z>(x); mode amplitudes
  sheet_2d.cpp         driver

  tests/               validation ladder (tests 0-4 plus the hybrid Poisson)
  tests/test_bending.cpp   Stage 1 measurements; `make measure`, not `make check`
  plot_phase_spiral.py Husimi maps, normalised by the pre-kick equilibrium
  plot_moments.py      corrugation history, envelope, conservation
  plot_bending.py      corrugation profile, mode frequency, §6 phase offset
```

## Building and running

Depends on FFTW3 only (`apt-get install libfftw3-dev`). The tridiagonal
eigensolver is self-contained deliberately, so no LAPACK is needed; it is
certified against the analytic harmonic ladder by `tests/test_harmonic.cpp`.

```sh
make            # drivers and tests
make check      # the validation ladder in order, ~4 min
make measure    # Stage 1 measurements, tens of minutes
./sheet_1d --help
./sheet_2d --help
```

The phase spiral, and the self-gravitating case that does not produce one:

```sh
./sheet_1d --N 4096 --L 8 --hbar 0.6 --sigma 5  --nu-ext 60 --z-ext 0.5 \
           --vkick 15 --periods 14 --out output/spiral_tracer
./sheet_1d --N 4096 --L 8 --hbar 0.6 --sigma 50 --nu-ext 60 --z-ext 0.5 \
           --vkick 15 --periods 14 --out output/spiral_selfgrav
python3 plot_phase_spiral.py output/spiral_tracer
python3 plot_moments.py output/spiral_tracer output/spiral_selfgrav
```

## Known limitations

- `phase_spiral_winding` in `src/husimi.h` measures the m=1 phase in polar
  coordinates on `(z, v_z/ω_0)`, which assumes unperturbed orbits are circles
  there. They are not, because the potential is anharmonic, and the resulting
  numbers are dominated by that mismatch rather than by winding. The **plots** are
  sound (they divide by the measured pre-kick equilibrium and so assume nothing),
  but a quantitative winding rate needs proper action-angle coordinates. That is
  Stage 1 work and it is the natural home for the §6 phase-offset diagnostic too.
- The perturber (§5) is not implemented; both stages use idealised kicks, which
  is what §9 steps 3 and 5 ask for.
- No shear, so no pattern winding: that is Stage 2 (§9 step 7).

## Suggested next step

Not Stage 2. Granule heating has to be brought under control first, because
every remaining Stage 1 deliverable (the dispersion relation, the firehose
threshold, Landau damping) depends on a dispersion that is currently drifting
faster than the physics being measured.

In order:

1. **Settle whether the heating is physical or numerical.** Rerun the stream
   scan with the velocity headroom held fixed *as the disk heats* (start with
   30 sigma, not 12) so the spill column stays negligible throughout. If the
   rate still refuses to fall with stream count, the `1/N_streams` argument is
   wrong and the mechanism needs identifying.
2. **Map the heating against `hbar_eff` and `sigma_x`.** The quasi-particle-mass
   argument predicts a steep dependence on both. If it holds, there is a usable
   window at smaller `hbar_eff`, and its cost in grid points is then known
   rather than guessed.
3. **Redo the dispersion relation** inside that window, with long enough runs
   for 20-plus zero crossings per mode, and only then attempt Landau damping.
4. Stage 2 becomes a decision point again once 1 to 3 give a heating rate, since
   a single wavefunction sits at the worst point of whatever scaling emerges.

The §9 ordering otherwise holds up. Step 1 was easier than expected
(correction 2), step 2 was done by a different route (correction 4), and step 4
is now built but blocked on the above rather than on missing machinery.
