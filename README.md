# Wave-Mechanics of Large Scale Structure

**Schrödinger Wave-mechanics and Large Scale Structure**

*PhD Thesis by Edward A. Thomson, University of Glasgow, 2011*

---

## About

This repository contains my 2011 PhD thesis, which developed a novel approach to cosmological structure formation simulations using the Schrödinger-Poisson system. Instead of tracking discrete particles (as N-body codes do), the method describes dark matter as a continuous complex wavefunction whose modulus squared gives the density field.

The code solved the coupled Schrödinger-Poisson equations in 3D with:
- Self-consistent gravity
- Expanding (cosmological) coordinates
- Periodic boundary conditions
- Mass and momentum conservation (via unitary, symplectic integration)

Results in the original thesis were compared with the N-body codes Hydra and GADGET-2. (The 2026 rewrite compares against GADGET-4 — see below.)

In 2026, I revisited the thesis chapter-by-chapter to assess it with fresh eyes. Independent reviews were conducted using Claude Opus 4.6 and OpenAI's GPT 5.4, then reconciled into a single set of notes. The original LaTeX source has been recovered, split into per-chapter files, and the review corrections have been applied to the source across all chapters and appendices (a handful of larger items — such as a consistent `I`/`we` pass and a rebuilt corrected PDF — remain outstanding). The simulation code has been rewritten from scratch in modern C++: the Free Particle Approximation (Chapter 4) and the full Schrödinger-Poisson solver (Chapter 5) are reimplemented and validated, now with an exact spectral (FFT) kinetic step in place of the original Goldberg finite-difference scheme, cosmological initial conditions (BBKS and GADGET-4), and a statistical comparison against GADGET-4 on matched initial conditions.

The rewrite also resolved a 15-year-old open question from the thesis. The original code produced density and velocity fields "messier" than N-body, which the thesis attributed to quantum interference and hoped to suppress by making ν (= ℏ/m) smaller. The rewrite shows the relationship is **backwards**: because velocity is encoded in the wavefunction phase (v = ν∇phase), the grid can only represent velocities up to v_max = ν·π·N. Shrinking ν at fixed resolution makes an unresolvable phase alias into noise. The thesis's ν ≈ 10⁻⁷–10⁻⁸ at N=64 aliased essentially the entire velocity field; correct growth is recovered at, e.g., ν=10⁻³ on N=128. See `code/schrodinger_poisson/output/investigation_narrative.md`.

## YouTube videos

* [Using 2 LLMs to review my PhD thesis](https://www.youtube.com/watch?v=LpO6d4BPOio)
* [Claude wrote my PhD code from scratch in C++](https://www.youtube.com/watch?v=J77JQkqD7NE)

## Repository Structure

```
wave-mechanics-lss/
├── README.md
├── project_roadmap.md            # Roadmap for thesis revival project
├── thesis/
│   ├── 2011thomsonphd.pdf        # Original thesis (250 pages)
│   └── latex/                    # Original LaTeX source and figures
│       ├── thesis.tex            # Master document
│       ├── ch01_introduction.tex # Per-chapter files (ch01–ch07)
│       ├── ...
│       ├── appendices/           # Appendix LaTeX files
│       ├── bkg_pix/              # Background/cosmology figures
│       ├── fpa_pix/              # Free Particle Approximation figures
│       ├── fpa_cos/              # FPA with cosmological ICs figures
│       └── full/                 # Full Schrödinger-Poisson figures
├── review/                       # Reconciled review notes (2026)
│   ├── notes_ch00_front_matter.md
│   ├── notes_ch01_introduction.md
│   ├── ...
│   └── sources/                  # Original separate reviews
│       ├── claude/
│       └── gpt/
├── code/                         # Modern C++ rewrite (core complete)
│   ├── fpa_1d/                   # 1D Free Particle Approximation (Ch 4)
│   ├── fpa_3d/                   # 3D FPA toy model (Ch 4)
│   ├── fpa_3d_cosmo/             # 3D FPA with cosmological ICs, ZA vs FPA (Ch 4)
│   ├── disk_bending/             # SP applied to galactic disk corrugations
│   │   ├── src/                  # Isolated-BC vertical Poisson, warm eigenstate
│   │   │                         #   ICs, split-step, Husimi phase-space transform
│   │   ├── sheet_1d.cpp          # 1D vertical slab: kick + phase spiral
│   │   ├── tests/                # Validation ladder (Poisson, free, harmonic,
│   │   │                         #   Jeans, warm sheet stationarity)
│   │   └── README.md             # Design notes and review of the adaptation plan
│   └── schrodinger_poisson/      # Full Schrödinger-Poisson solver (Ch 5)
│       ├── sp_1d.cpp             # 1D solver (Goldberg + periodic BCs)
│       ├── sp_3d.cpp             # 3D solver (Goldberg, tophat tests, expansion)
│       ├── sp_3d_cosmo*.cpp      # Cosmological solvers (spectral kinetic step,
│       │                         #   BBKS/GADGET-4 ICs, N=128/256, ν variants)
│       ├── extract_gadget_ics.py # TSC deposit + Poisson velocity from GADGET-4 HDF5
│       ├── compare_statistics.py # P(k), density PDF, cross-correlation, rms growth
│       ├── algorithm.md          # Algorithm reference document
│       └── output/               # Run outputs, GADGET-4 comparison, and the
│                                 #   investigation write-ups (nu, growth, bulk flow)
└── web/                          # Web version of thesis (planned)
```

## Review Summary

| Chapter | Pages | Items |
|---------|-------|-------|
| 1. Introduction | 35 | 36 |
| 2. Review of Numerical Simulations | 16 | 27 |
| 3. Wave-mechanics | 35 | 35 |
| 4. Free Particle Approximation | 24 | 45 |
| 5. Solving the full Schrödinger-Poisson system | 67 | 78 |
| 6. Conclusion | 13 | 47 |
| 7. Epilogue: Vorticity and spin | 32 | 68 |
| Appendix A: Translation | 6 | 14 |
| Appendix B: Mathematical appendix | 7 | 27 |
| **Total** | **235** | **377** |

Items include typos, grammatical errors, mathematical issues, clarity improvements, and structural comments. The reconciled notes identify shared findings and reviewer-specific additions. Original reviews are preserved under `review/sources/`.

## Code Rewrite (2026)

The solver has been rebuilt from scratch in C++. It evolves a complex wavefunction ψ on a periodic grid (|ψ|² = density, phase gradient = velocity) with self-consistent gravity via a Poisson solve each step.

Key differences from the 2011 code:

| Aspect | Thesis (2011) | Rewrite (2026) |
|--------|---------------|----------------|
| Language | Fortran | C++ |
| Kinetic step | Goldberg finite-difference | Spectral FFT (exact propagator) |
| Periodic BCs | Double-recursion auxiliary functions | Inherent in FFT |
| Unitarity | Exact (Cayley) | Exact (Cayley + spectral) |
| Poisson solve | FFT with discrete Green's function | FFT with continuum k⁻² |
| Typical run | Hours on a cluster | Minutes on a laptop (~3 min at N=128) |

**Headline finding — the ν resolution constraint.** The wavefunction can only represent velocities up to **v_max = ν·π·N**. The thesis's "quantum interference" noise was numerical aliasing from an unresolvable Madelung phase at too-small ν; correct linear growth is recovered once ν·N is large enough (e.g. ν=0.01 at N=64, ν=10⁻³ at N=128, ν=10⁻⁴ at N=256). The diagnosis came from a systematic elimination of hypotheses (sub-stepped Cayley, spectral vs Goldberg, zero-velocity ICs, a Madelung round-trip test), each ruling out one candidate cause.

**GADGET-4 comparison.** On matched initial conditions (128³, 32 Mpc/h box), the codes agree in the linear regime (z≈19) and diverge nonlinearly: the S-P code produces excess small-scale power and a broader density PDF — "more lows and more highs" — reproducing the thesis's original comparison against GADGET-2, but now with matched ICs confirming it is physical rather than a random-seed artefact.

Full write-ups are under `code/schrodinger_poisson/output/` (`project_summary.md`, `investigation_narrative.md`, `growth_rate_investigation.md`, and the GADGET-4 `comparison_report.md`).

## Context

Since 2011, the approach explored in this thesis has become an active research area under names including "fuzzy dark matter," "ultralight axion dark matter," and "wave dark matter." Notable subsequent work includes Schive et al. (2014), Mocz et al. (2017), and many others who developed production-quality Schrödinger-Poisson solvers. The ideas in this thesis — particularly the connection between the Schrödinger equation, Madelung fluid equations, and cosmological structure formation — remain relevant.

## Thesis Abstract

The thesis presents the Schrödinger equation as an alternative to the Euler-Poisson system for simulating the evolution of large scale structure in the Universe. Through the Madelung transform, the Schrödinger equation can be recast as fluid equations, providing a natural framework for describing collisionless dark matter as a continuous field rather than discrete particles. The Free Particle Approximation (FPA) is verified and extended to 3D cosmological initial conditions. A full 3D Schrödinger-Poisson solver is developed using the Cayley transform (Goldberg scheme) with Suzuki splitting operators, and results are compared with the N-body code GADGET-2. An epilogue explores speculative extensions: gravitoelectromagnetism for vorticity, and a novel Pauli-like equation for spinning objects derived via Clifford algebra.

## License

- Thesis text and review notes: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
- Code: MIT License

The original thesis was deposited with the University of Glasgow under their standard thesis deposit terms ([http://theses.gla.ac.uk/2976/](http://theses.gla.ac.uk/2976/)), which permit downloading for personal non-commercial research or study. As the copyright holder, the author has chosen to release this work here under the more permissive CC BY 4.0 license.

## Author

Edward A. Thomson
- Thesis supervised by Prof. Martin Hendry, University of Glasgow
- Original thesis submitted 2011, review conducted 2026
