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

Results were compared with the N-body codes Hydra and GADGET-2.

In 2026, I revisited the thesis chapter-by-chapter to assess it with fresh eyes. Independent reviews were conducted using Claude Opus 4.6 and OpenAI's GPT 5.4, then reconciled into a single set of notes. The original LaTeX source has been recovered and split into per-chapter files, ready for corrections. The simulation code is being rewritten from scratch in modern C++, with the Free Particle Approximation (Chapter 4) and the full Schrödinger-Poisson solver (Chapter 5) now reimplemented and validated.

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
├── code/                         # Modern C++ rewrite (in progress)
│   ├── fpa_1d/                   # 1D Free Particle Approximation (Ch 4)
│   ├── fpa_3d/                   # 3D FPA toy model (Ch 4)
│   ├── fpa_3d_cosmo/             # 3D FPA with cosmological ICs (Ch 4)
│   └── schrodinger_poisson/      # Full Schrödinger-Poisson solver (Ch 5)
│       ├── sp_1d.cpp             # 1D solver (Goldberg + periodic BCs)
│       ├── sp_3d.cpp             # 3D solver (splitting operators + expansion)
│       ├── algorithm.md          # Algorithm reference document
│       └── output/               # Test outputs (free particle, gravity, tophat)
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
