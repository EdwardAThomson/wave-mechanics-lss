# Project Roadmap: Wave-Mechanics LSS — Thesis Revival

## Repository Name
**`wave-mechanics-lss`**

GitHub: `github.com/edwardthomson/wave-mechanics-lss`

**Status at a glance (updated 2026-07):**
- ✅ **Phase 1 — GitHub archive:** complete
- 🔬 **Phase 2 — Web conversion:** not started (LaTeX source recovered; corrections and hosting pending)
- ✅ **Phase 3 — Code rewrite:** core complete (1D + 3D FPA and full Schrödinger-Poisson reimplemented, validated, and cosmological runs done); a few stretch items remain
- 🔄 **Phase 4 — Blog / outreach:** underway (2 YouTube videos published; written posts pending)

---

## Phase 1: GitHub Repository Setup ✅ COMPLETE
**Goal:** Get the thesis and review notes online as a public archive.

- [x] Create the GitHub repository with a clear README
  - Thesis title, abstract, author, university, year (2011)
  - Context on the 2011 → 2026 revival and how the review was conducted
- [x] Upload the original thesis PDF (as-is, uncorrected) — `thesis/2011thomsonphd.pdf`
- [x] Recover the original LaTeX source and split it into per-chapter files — `thesis/latex/`
- [x] Upload the chapter-by-chapter review notes — now under `review/notes_*.md`
  - Reconciled from two independent LLM reviews (Claude Opus 4.6 + GPT 5.4), preserved under `review/sources/`
  - Per-chapter correction checklists added under `review/checklists/`
- [x] Add a LICENSE (CC BY 4.0 for text/notes; MIT for code)
- [x] Organise repo structure (see README for the current layout)

**Outcome:** The archive is live. Structure evolved from the original plan — review
notes live under `review/` (not `chapters/`), the LaTeX source was recovered (a
head start on Phase 2), and `code/` already holds the working rewrite.

---

## Phase 2: Web Conversion — Corrected Edition 🔬 NOT STARTED
**Goal:** Convert the thesis into a web-hosted version with corrections applied.

The recovered per-chapter LaTeX source (`thesis/latex/`) is the natural starting
point for this phase — it removes the need to reconstruct equations from the PDF.

### Phase 2a: Conversion Pipeline
- [ ] Choose a static site framework (e.g. Astro, Hugo, or plain HTML + KaTeX)
- [ ] Convert the per-chapter LaTeX to Markdown/HTML (source already recovered)
- [ ] Render equations with KaTeX or MathJax
- [ ] Proof every equation against the original PDF (use review notes as checklist)

### Phase 2b: Apply Corrections
- [ ] Fix typos and grammatical errors (~377 items catalogued across all chapters — see `review/`)
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
- [ ] Fold in results and figures from the Phase 3 rewrite (the corrected edition can now cite reproduced, higher-resolution results)

**Estimated effort:** A chapter per session, roughly 9–12 sessions for the full
thesis. Equation proofing remains the bottleneck, though the recovered LaTeX
source makes it lighter than converting from PDF.

---

## Phase 3: Code Rewrite ✅ CORE COMPLETE
**Goal:** A modern, performant implementation of the Schrödinger-Poisson solver.

Rewritten from scratch in modern **C++**. Lives under `code/`. See
`code/schrodinger_poisson/output/project_summary.md` for the full write-up.

### Phase 3a: Design Decisions ✅
- [x] Language chosen: **C++**
- [x] Numerical method chosen: **spectral (FFT) split-operator** for the kinetic step
  (exact free-particle propagator), with an **exactly-unitary Cayley potential step**
  and **Strang splitting** — the spectral method replaced the original thesis's
  Goldberg finite-difference scheme (~14× faster, naturally periodic, mass-conserving to ~10⁻¹³)
- [x] Scope decided: 1D prototype first, then 3D

### Phase 3b: Implementation (1D Prototype) ✅
- [x] Free Particle Approximation, 1D — reproduces Chapter 4 results (`code/fpa_1d/`)
- [x] Full Schrödinger-Poisson, 1D, with split-operator (`code/schrodinger_poisson/sp_1d.cpp`)
- [x] Verify against known cases (free particle, gravity, tophat collapse) — see `output/`
- [x] Compare FPA vs full S-P

### Phase 3c: Implementation (3D) ✅
- [x] 3D FFT-based solver (`code/schrodinger_poisson/sp_3d*.cpp`)
- [x] 3D FPA, including the cosmological-IC toy model (`code/fpa_3d/`, `code/fpa_3d_cosmo/`)
- [x] Cosmological initial conditions — **BBKS transfer function** (Gaussian random field
  + Madelung wavefunction) and **GADGET-4 IC ingestion** via TSC deposit
- [x] Periodic boundary conditions (inherent in FFT)
- [x] Poisson solver (Fourier space, continuum k⁻²)
- [x] Output: density fields, velocity fields, power spectra
- [x] Compare with N-body / reference — Zel'dovich approximation and GADGET IC comparisons,
  with density PDFs, power spectra, and summary statistics
- [x] Higher-resolution runs — N=128 (~3 min) and N=256 (~55 min) on a laptop
- [ ] **Stretch:** direct side-by-side against a live GADGET-2/GADGET-4 N-body run (beyond IC matching)
- [ ] **Stretch:** 512³ runs (need more than a laptop)

### Phase 3d: Investigate Open Questions from Thesis ✅ (headline result solved)
- [x] **The ν / velocity-noise mystery — solved.** The thesis attributed messy density and
  velocity fields to "quantum interference" and tried to suppress it by shrinking ν. The
  relationship is *backwards*: at fixed N, smaller ν makes an unresolvable Madelung phase
  alias into noise. Correct results require ν·N large enough to resolve the phase
  (e.g. ν=10⁻³ at N=128 gives correct linear growth and a recognisable cosmic web).
  See `output/investigation_narrative.md`, `growth_rate_investigation.md`, `bulk_flow_investigation.md`.
- [x] Interference effects — characterised as a resolution/aliasing effect (quantum fringes are physical only when adequately resolved)
- [x] ν parameter sensitivity — systematic scan (8 experimental variants, each ruling out a hypothesis)
- [ ] **Stretch:** convergence study pushed to 512³ for a resolution-independence statement

**Outcome:** The core rewrite achieves the phase's goal — the method is reproduced,
validated, and the central 15-year-old open question is resolved. Remaining items are
stretch goals requiring more compute than a laptop.

---

## Phase 3e: New Application — Galactic Disk Corrugations 🔄 STAGE 0 COMPLETE
**Goal:** Point the same Schrödinger-Poisson machinery at vertical bending waves and
the Gaia vertical phase spiral instead of cosmological structure. Lives in
`code/disk_bending/`; see its `README.md`.

Schrödinger-Poisson is a Vlasov-Poisson solver, so anything collisionless is fair game.
The draw is that Husimi-smoothing |ψ|² gives the full distribution function f(z, v_z) on
a grid with **no particle shot noise**, which is exactly where N-body struggles for the
fine-grained phase spiral.

- [x] Isolated (vacuum) vertical Poisson solver — exact and O(N) in 1D, no Green's-function padding needed
- [x] Warm initial conditions from eigenstates of the self-consistent vertical Hamiltonian
      (stationary by construction, avoiding the turning-point singularity of stream sampling)
- [x] Validation ladder: free particle, harmonic oscillator, Jeans growth (both branches), sheet stationarity
- [x] Husimi transform to f(z, v_z); vertical phase spiral produced from an impulsive kick
- [x] Settled the multi-stream vs single-ψ question quantitatively:
      spurious dipole A ≈ 0.70·√(ħ_eff/ω_z)/√N_streams
- [ ] Stage 1: 2D (x, z) box — bending-wave dispersion relation and Landau damping
- [ ] Stage 2: shearing box and/or a Sagittarius-like perturber

**Notable findings:** a purely harmonic vertical potential produces no phase spiral at
all (Kohn's theorem), and strong disk self-gravity keeps the dipole response coherent
so it neither damps nor winds. How much vertical restoring force comes from the disk
versus the halo is therefore a physics parameter to scan, not a detail.

---

## Phase 4: Blog / Outreach 🔄 UNDERWAY
**Goal:** Document the journey and share insights.

- [x] Video: [Using 2 LLMs to review my PhD thesis](https://www.youtube.com/watch?v=LpO6d4BPOio)
- [x] Video: [Claude wrote my PhD code from scratch in C++](https://www.youtube.com/watch?v=J77JQkqD7NE)
- [ ] Blog post: "Revisiting my PhD thesis 15 years later" (Phase 1 launch)
- [ ] Blog post: "377 mistakes in my thesis" (review notes summary)
- [ ] Blog post: "The ν mystery — a 15-year-old bug that wasn't quantum" (Phase 3d result)
- [ ] Blog series: key physics concepts made accessible
  - The Madelung transform: quantum mechanics as fluid dynamics
  - Why simulate dark matter with waves?
  - What is fuzzy dark matter?
- [ ] Blog posts as further code milestones are reached

---

## Notes
- Phase 1 is done; Phase 2 (web edition) is the main remaining thread.
- Phase 3 ran ahead of Phase 2 rather than after it — the rewrite and the ν investigation
  are essentially complete, so the corrected web edition can now cite reproduced results.
- The recovered LaTeX source (`thesis/latex/`) is the launch point for Phase 2 and removes
  the "convert from PDF" step from the original plan.
- The review notes (377 items) remain the foundation for Phase 2b corrections.
- The web version is the "living" corrected edition; the PDF is the historical record.
- Next up: kick off Phase 2a (framework choice + LaTeX→web pipeline) and start drafting the Phase 4 posts, which now have concrete results to point at.
