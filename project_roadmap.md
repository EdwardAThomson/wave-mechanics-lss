# Project Roadmap: Wave-Mechanics LSS — Thesis Revival

## Repository Name
**`wave-mechanics-lss`** (or `schrodinger-poisson-lss`)

GitHub: `github.com/edwardthomson/wave-mechanics-lss`

---

## Phase 1: GitHub Repository Setup
**Goal:** Get the thesis and review notes online as a public archive.

- [ ] Create the GitHub repository with a clear README
  - Thesis title, abstract, author, university, year (2011)
  - Brief context: what the thesis is about and why it's being revisited
  - Link to the original university repository if available
- [ ] Upload the original thesis PDF (as-is, uncorrected)
- [ ] Upload the chapter-by-chapter review notes (from `/chapters/notes_*.md`)
- [ ] Add a LICENSE (CC BY 4.0 or similar for the text; MIT for any future code)
- [ ] Organise repo structure:
  ```
  wave-mechanics-lss/
  ├── README.md
  ├── LICENSE
  ├── thesis/
  │   └── 2011thomsonphd.pdf
  ├── review/
  │   ├── notes_ch01_introduction.md
  │   ├── notes_ch02_review_numerical_simulations.md
  │   ├── notes_ch03_wave_mechanics.md
  │   ├── notes_ch04_free_particle_approximation.md
  │   ├── notes_ch05_schrodinger_poisson_system.md
  │   ├── notes_ch06_conclusion.md
  │   ├── notes_ch07_epilogue_vorticity_spin.md
  │   ├── notes_appendix_a_translation.md
  │   └── notes_appendix_b_mathematical.md
  ├── web/          (Phase 2)
  └── code/         (Phase 3)
  ```

**Estimated effort:** A few hours.

---

## Phase 2: Web Conversion — Corrected Edition
**Goal:** Convert the thesis into a web-hosted version with corrections applied.

### Phase 2a: Conversion Pipeline
- [ ] Choose a static site framework (e.g. Astro, Hugo, or plain HTML + KaTeX)
- [ ] Convert the PDF chapter-by-chapter to Markdown/HTML using AI assistance
- [ ] Render equations with KaTeX or MathJax
- [ ] Proof every equation against the original PDF (use review notes as checklist)

### Phase 2b: Apply Corrections
- [ ] Fix typos and grammatical errors (straightforward, ~100 items across all chapters)
- [ ] Fix mathematical errors (sign conventions, Eq 5.62 Poisson formulation, Eq 7.4 angular momentum, etc.)
- [ ] Add correction annotations where substantive changes are made
  - Format: visible footnotes or margin notes showing what changed and why
  - Preserve intellectual honesty — this is a correction, not a rewrite
- [ ] Add cross-references between chapters where the review notes identified missing links
- [ ] Do NOT rewrite prose style or restructure arguments — keep the original voice

### Phase 2c: Host on Personal Website
- [ ] Deploy to personal site (e.g. `edwardthomson.com/thesis/` or similar)
- [ ] Add navigation: table of contents, chapter links, search
- [ ] Include the review notes as a companion "commentary track"
- [ ] Add a preface explaining the 2011 → 2026 context

**Estimated effort:** A chapter per session, roughly 9–12 sessions for the full thesis. The equation proofing is the bottleneck.

---

## Phase 3: Code Rewrite
**Goal:** A modern, performant implementation of the Schrödinger-Poisson solver.

### Phase 3a: Design Decisions
- [ ] Choose language (candidates: Julia, Rust, C++, Python+CuPy)
  - Julia: closest to mathematical notation, good FFT libraries, easy prototyping
  - Rust: performance + safety, steeper learning curve
  - C++: established in computational physics, good MPI support
  - Python: fastest to prototype, can offload to GPU via CuPy/JAX
- [ ] Choose numerical method:
  - FFT-based split-operator (used by Schive et al. 2014, Mocz et al. 2017) — trivially parallelisable, simpler than Cayley
  - Cayley/Goldberg scheme (as in original thesis) — unitary, but hard to parallelise
  - Recommendation: start with FFT split-operator, it's the community standard now
- [ ] Decide scope: 1D prototype first, then 3D

### Phase 3b: Implementation (1D Prototype)
- [ ] Free Particle Approximation (reproduce Chapter 4 results)
- [ ] Full Schrödinger-Poisson with split-operator
- [ ] Verify against known analytical solutions (tophat collapse, two-body)
- [ ] Compare FPA vs full S-P

### Phase 3c: Implementation (3D)
- [ ] 3D FFT-based solver
- [ ] Cosmological initial conditions (BBKS or CLASS transfer function)
- [ ] Periodic boundary conditions (natural with FFT)
- [ ] Poisson solver (Fourier space, trivial with FFT)
- [ ] Output: density fields, velocity fields, power spectra
- [ ] Compare with N-body results (GADGET or similar)

### Phase 3d: Investigate Open Questions from Thesis
- [ ] Velocity discrepancy with N-body codes — is it a resolution effect or fundamental?
- [ ] Interference effects — how do they scale with resolution?
- [ ] ν parameter sensitivity — systematic convergence study
- [ ] Higher resolution runs (256³, 512³) now feasible on modern hardware

**Estimated effort:** Phase 3a-3b in a few weeks. Phase 3c-3d is an ongoing project.

---

## Phase 4: Blog / Outreach (Optional, Ongoing)
**Goal:** Document the journey and share insights.

- [ ] Blog post: "Revisiting my PhD thesis 15 years later" (Phase 1 launch)
- [ ] Blog post: "377 mistakes in my thesis" (review notes summary)
- [ ] Blog series: key physics concepts made accessible
  - The Madelung transform: quantum mechanics as fluid dynamics
  - Why simulate dark matter with waves?
  - What is fuzzy dark matter?
- [ ] Blog posts as code milestones are reached (Phase 3)

---

## Notes
- Phases 1 and 2 are sequential (need repo before web conversion)
- Phase 3 can run in parallel with Phase 2
- Phase 4 is opportunistic — write when there's something to say
- The review notes (377 items) are the foundation for Phase 2b corrections
- The web version is the "living" corrected edition; the PDF is the historical record
