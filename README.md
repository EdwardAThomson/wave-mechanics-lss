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

In 2026, I revisited the thesis chapter-by-chapter to assess it with fresh eyes. The detailed review notes (377 items across 9 chapters/appendices) are included in this repository.

## Repository Structure

```
wave-mechanics-lss/
├── README.md
├── project_roadmap.md        # Roadmap for thesis revival project
├── thesis/
│   └── 2011thomsonphd.pdf    # Original thesis (250 pages)
├── review/                   # Chapter-by-chapter review notes (2026)
│   ├── notes_ch01_introduction.md
│   ├── notes_ch02_review_numerical_simulations.md
│   ├── notes_ch03_wave_mechanics.md
│   ├── notes_ch04_free_particle_approximation.md
│   ├── notes_ch05_schrodinger_poisson_system.md
│   ├── notes_ch06_conclusion.md
│   ├── notes_ch07_epilogue_vorticity_spin.md
│   ├── notes_appendix_a_translation.md
│   └── notes_appendix_b_mathematical.md
├── web/                      # Web version of thesis (planned)
└── code/                     # Modern code rewrite (planned)
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

Items include typos, grammatical errors, mathematical issues, clarity improvements, and structural comments. Each review file contains an equations check verifying the mathematical content.

## Context

Since 2011, the approach explored in this thesis has become an active research area under names including "fuzzy dark matter," "ultralight axion dark matter," and "wave dark matter." Notable subsequent work includes Schive et al. (2014), Mocz et al. (2017), and many others who developed production-quality Schrödinger-Poisson solvers. The ideas in this thesis — particularly the connection between the Schrödinger equation, Madelung fluid equations, and cosmological structure formation — remain relevant.

## Thesis Abstract

The thesis presents the Schrödinger equation as an alternative to the Euler-Poisson system for simulating the evolution of large scale structure in the Universe. Through the Madelung transform, the Schrödinger equation can be recast as fluid equations, providing a natural framework for describing collisionless dark matter as a continuous field rather than discrete particles. The Free Particle Approximation (FPA) is verified and extended to 3D cosmological initial conditions. A full 3D Schrödinger-Poisson solver is developed using the Cayley transform (Goldberg scheme) with Suzuki splitting operators, and results are compared with the N-body code GADGET-2. An epilogue explores speculative extensions: gravitoelectromagnetism for vorticity, and a novel Pauli-like equation for spinning objects derived via Clifford algebra.

## License

- Thesis text and review notes: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
- Code (when added): MIT License

The original thesis was deposited with the University of Glasgow under their standard thesis deposit terms ([http://theses.gla.ac.uk/2976/](http://theses.gla.ac.uk/2976/)), which permit downloading for personal non-commercial research or study. As the copyright holder, the author has chosen to release this work here under the more permissive CC BY 4.0 license.

## Author

Edward A. Thomson
- Thesis supervised by Prof. Martin Hendry, University of Glasgow
- Original thesis submitted 2011, review conducted 2026
