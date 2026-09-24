# Paper outline: wave-mechanical simulation of disk vertical dynamics

Working document for paper 1 from `code/disk_bending/`. Everything cited as a
result here is measured and committed; the "remaining" table tracks what is
still wanted before submission.

## Framing and novelty claim

**Claim:** first application (to our knowledge, checked August 2026) of the
Schrödinger method to the stellar phase-space dynamics of a galactic disk:
the Gaia-style vertical phase spiral and the bending-wave dispersion
relation, computed with zero shot noise on desktop hardware.

Positioning against the three adjacent literatures:

1. **The Schrödinger method itself** (Widrow & Kaiser 1993; Davies & Widrow
   1997; Uhlemann, Kopp & Haugg 2014; Kopp, Vattis & Skordis 2017, PRD 96
   123532; Mocz et al.; Eberhardt et al.; Thomson 2011, the thesis this
   repository revives). All cosmological or generic Vlasov work; none
   applies it to a stellar disk's vertical structure.
2. **FDM and the disk** (disc heating by FDM granulation, MNRAS 530:129,
   2024; dark subhalos and the Gaia snail, ApJ 2025, arXiv:2412.02757; the
   dark-matter-wake snail, MNRAS 524:801, 2023). In all of these the SP
   field is the dark matter and the stars are particles. We invert that: the
   SP field IS the stellar distribution function, with hbar_eff a numerical
   coarse-graining scale, not a boson mass.
3. **Phase-spiral modelling practice** (Antoja et al. 2018; Laporte et al.
   2019; Bland-Hawthorn & Tepper-Garcia 2021; Widrow et al. 2012/2014
   galactoseismology; Hunt et al. 2022; Weinberg-style kinetic theory).
   N-body needs 1e8-1e9 particles because the spiral is percent-level
   contrast in f(z, v_z); our method delivers f directly, limited by
   hbar_eff rather than counting statistics.

The June 2026 WaveDM.jl framework (arXiv:2606.25026) is the nearest current
code: SP + N-body baryons for galaxy scales, but again with the wave sector
as dark matter. Cite it as evidence the tooling era is arriving; our angle
(stars as the wave field, controlled linear measurements, validation ladder)
stays distinct.

## Target and shape

MNRAS, standard article. Methods plus first results. Roughly 12-15 pages.
Author: Edward A. Thomson. (Acknowledgement of the 2011 thesis lineage and
Martin Hendry as PhD supervisor.)

Title candidates (pick late):
- "Wave-mechanical simulation of galactic disk corrugations and the vertical
  phase spiral"
- "The Gaia phase spiral without shot noise: a Schrödinger-Poisson approach
  to disk vertical dynamics"
- "Schrödinger-Poisson simulation of bending waves and phase spirals in
  stellar disks"

## Abstract skeleton (5 sentences)

1. SP is a Vlasov solver: Husimi-smoothing |psi|^2 recovers f to
   O(hbar_eff), so anything collisionless is fair game, including disk
   vertical dynamics.
2. We build a validated solver ladder (1D slab, 2D unsheared box, 2D
   rigidly rotating box) with warm equilibria constructed from eigenstates
   of the self-consistent Hamiltonian, exact isolated-boundary Poisson
   solves, and a Landau-gauge treatment of Coriolis.
3. Results: a noise-free Antoja-style phase spiral, whose winding requires
   anharmonicity (Kohn) and is suppressed by disk self-gravity (65% vs 0.4%
   corrugation retained); a no-go duality for unsheared boxes; Toomre
   stabilisation measured, widened by finite thickness beyond razor-thin Q.
4. The stable bending branch measured at [2-5] wavenumbers: omega
   reproducible to 0.15% across phase realisations, full free-streaming
   pressure term excluded at 4.8x, both corrections (thickness, epicyclic
   confinement) acting at their expected sign and size.
5. Cost: minutes to an hour on a desktop; the equivalent shot-noise-free
   N-body measurement needs 1e9-1e10 particles in the patch.

## Sections

1. **Introduction.** Gaia phase spiral; the shot-noise problem for
   fine-grained f; the Schrödinger method and its cosmological lineage; this
   paper's inversion (stars as the wave field).
2. **Method.**
   2.1 SP as coarse-grained Vlasov; hbar_eff and the velocity ceiling
       v_max = hbar_eff pi / dx (the resolution law, worth its own figure).
   2.2 Warm equilibria from self-consistent eigenstates (vertical), thermal
       epicyclic Landau levels (in-plane, rotating case); why sampled
       streams fail (turning-point divergence).
   2.3 Exact isolated Poisson solves: O(N) cumulative sweep (1D),
       exponential-kernel recursion per k_x (2D).
   2.4 Rotation: Landau gauge, per-stream wrapped traps, no edge buffer;
       stream lattice and level truncation.
   2.5 Validation ladder table (orders of convergence, conservation).
3. **The vertical phase spiral** (Stage 0 results, 3 figures already made).
4. **Bending waves without rotation: a no-go and the firehose** (duality,
   fragmentation control, firehose suppression).
5. **The rotating box: Toomre stabilisation and the stable bending branch**
   (four-case stability table; dispersion relation figure; damping if
   measured; thickness lesson).
6. **Discussion.** Comparison with N-body cost; limitations (no shear, local
   box stretched at 16 kpc, hbar_eff coarseness, no N-body cross-check yet);
   the FDM follow-up (hbar_eff scan makes the quantum term physical: paper
   2).

## Figure list

| # | figure | status |
|---|---|---|
| 1 | phase spiral, tracer (Husimi map sequence) | exists, `figures/phase_spiral_tracer.png` |
| 2 | phase spiral, self-gravitating + moments comparison | exists (2 files) |
| 3 | resolution law: spill/heating vs hbar_eff N (from Stage 0/1 data) | to assemble from existing runs |
| 4 | dispersion relation omega^2/2piGSk vs k with candidate curves | needs the modes 3-5 runs, then a plot script |
| 5 | Toomre four-case summary (band amplitude vs time) | needs re-run with series output or table-only |
| 6 | validation ladder table | exists in READMEs |

## Remaining before submission

| item | status |
|---|---|
| modes 3-5 dispersion points + damping envelopes | running now |
| Nx = 2048 convergence check (sigma_x calibration) | queued |
| plot script for the dispersion figure | todo |
| decision: N-body cross-check or honest absence statement | user's call |
| literature check | done (August 2026 searches; gap confirmed as far as searchable) |
| LaTeX skeleton | todo after numbers land |

## Paper 2 (deliberately excluded from paper 1)

The hbar_eff scan: measure the dispersion relation as hbar_eff grows until
the quantum term hbar^2 k^4 / 4 bends it away from the classical curve. For
stars hbar_eff is numerical; for an FDM component it is m = hbar/hbar_eff,
and the measurement becomes a physical prediction for wave-dark-matter
disks, connecting to the FDM granulation-heating and Gaia-snail literature.
The machinery (test_rot_disp plus a loop) already exists.
