# wave-mechanics-lss

Revival of Edward A. Thomson's 2011 PhD thesis (University of Glasgow): simulating
cosmological large-scale structure with the Schrödinger-Poisson system instead of
N-body particles. A complex wavefunction ψ evolves on a periodic grid where |ψ|² is
the density and the phase gradient (v = ν∇phase) is the velocity. Three strands:
archive the original thesis, publish a corrected web edition, and a modern C++ rewrite
of the solver.

## Layout

- `thesis/` — original PDF (`2011thomsonphd.pdf`) and recovered LaTeX source, split per
  chapter under `thesis/latex/` (`ch01`–`ch07`, `appendices/`, figure dirs `*_pix`/`*_cos`).
- `review/` — reconciled 2026 review notes (377 items, per-chapter `notes_*.md`).
  Original separate Claude and GPT reviews are preserved under `review/sources/`.
- `code/` — the C++ rewrite. `fpa_1d`, `fpa_3d`, `fpa_3d_cosmo` (Free Particle
  Approximation, Ch 4) and `schrodinger_poisson` (full solver, Ch 5).
- `web/` — planned web edition of the thesis (Phase 2, not built yet).
- `project_roadmap.md` — the phased plan (archive → web edition → code rewrite → outreach).

## Building and running the code

Each code subdir has its own `Makefile`; build with `make` inside it. In
`code/schrodinger_poisson/` the solver has many variants (`sp_1d`, `sp_3d`,
`sp_3d_cosmo*`) each with a matching `.cpp`; Python scripts (`plot_*.py`,
`compare_statistics.py`, `extract_gadget_ics.py`) handle analysis and GADGET-4
comparison. Run outputs and the investigation write-ups live in
`code/schrodinger_poisson/output/` — read `investigation_narrative.md`,
`growth_rate_investigation.md`, and `comparison_report.md` there before touching the
solver, they carry the physics reasoning.

## The one thing to know about the physics

The wavefunction can only represent velocities up to **v_max = ν·π·N** (ν = ℏ/m, N =
grid size). The thesis's "quantum interference" noise was numerical aliasing from an
unresolvable Madelung phase at too-small ν, not a physical effect. Correct linear
growth needs ν·N large enough (ν=0.01 at N=64, ν=1e-3 at N=128, ν=1e-4 at N=256).
This resolution constraint is the headline result of the rewrite; keep it in mind for
any change to ν, grid size, or the kinetic step.

## Working conventions

- The PDF is the historical record and stays as-is. The LaTeX source and web edition
  are the corrected "living" version. When applying review corrections, fix the error
  and annotate it; do not rewrite the original prose voice or restructure arguments.
- Corrections are driven by `review/notes_*.md` — treat those as the checklist.
- Author: Edward A. Thomson; original supervisor Prof. Martin Hendry.
- Licensing: text and notes CC BY 4.0; code MIT.

## House style (from the author's global preferences)

No em dashes in any prose (commas, colons, parentheses, or separate sentences instead).
Lead with the answer, keep it concise, and be explicit about anything that failed, was
skipped, or is uncertain.
