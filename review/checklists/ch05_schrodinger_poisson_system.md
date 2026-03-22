# Chapter 5: Solving the Full Schrödinger-Poisson System — Checklist

## Copy-editing (Shared #11-12)
- [x] Repeated articles ("the the", "a an"), Widrow/Widow typos — fixed
- [x] "Through out" → "Throughout"
- [x] "colloquially" → "commonly"
- [x] "as sen in" → "as seen in"
- [x] "how implement" → "how to implement" (summary)
- [ ] Pronoun inconsistency (`I` vs `we`) — deferred: pervasive throughout

## Poisson-equation presentation — Major fix (Shared #4)
- [x] Added mean-density subtraction clarification for periodic box
- [x] Removed misleading standalone integral expression for the potential
- [x] Replaced `ψ_k ψ_k*` with FT of `|ψ|²`; noted convolution distinction
- [x] Added V̂₀ = 0 convention and discrete Green's function note

## Mathematical errors (Shared)
- [x] Eq 5.64: 1D Gaussian normalisation — fixed exponent from 1/2 to 1/4 (#5)
- [x] Eq 5.69: momentum expression — fixed |ψ*∇ψ|² to -iν ψ*∇ψ (#6)
- [x] Conservation discussion — removed claim that mass+momentum → energy conservation (#6)

## Major unresolved scientific issues (Shared)
- [x] Velocity discrepancy with GADGET — added dedicated "Known Issues" subsection (#7)
- [x] Spurious growing bulk flow — reframed as major limitation with possible causes (#8)

## Conceptual / substantive (Shared)
- [N/A] Chapter too long — structure preserved; tightened where possible (#10)
- [x] Clearer signalling of new vs background — GADGET switch motivated, Gaussian smoothing labelled as diagnostic
- [x] Explicit separation of solid results and unresolved issues — new subsection added

## Conceptual (GPT-only)
- [x] p113-114: basis discussion — clarified V diagonal in position space, K tridiagonal; explained two-step Poisson/Schrödinger distinction
- [x] p118-120: implicit vs semi-implicit — identified as Crank-Nicolson (implicit midpoint), unconditionally stable, second-order
- [ ] p149-150: top-hat ν inconsistency — deferred: needs re-reading of top-hat results section
- [x] p133: switch from Hydra to GADGET-2 — motivated explicitly
- [ ] p153-156: resolution vs framework issues — deferred
- [x] p155-164: "inevitable"/"fundamental" claims — softened to "characteristic feature" / "open question"
- [x] p158-166: Gaussian smoothing — stated as diagnostic post-processing, not physical
- [N/A] p174: no-vorticity claim — language already appropriately cautious

## Conceptual (Claude-only)
- [N/A] p113: Cayley attribution — text already attributes to Goldberg via Widrow & Kaiser
- [x] p118: convention m = 1/2 — justified as standard simplification, constants restored later
- [x] p119: reusing δ for timestep — added disambiguation note
- [x] p115: Cayley factors commute — justified (polynomials of same matrix)
- [ ] p116-117: Goldberg → Suzuki/Strang splitting — deferred: substantial rewrite
- [x] p116: ad hoc V/3 split — motivated as naive attempt, explained why it fails (no theoretical basis, K/V non-commutativity)
- [ ] p119-129: recursion relations — deferred: worked example would be valuable but large
- [x] p123: boundary recursion convergence — added note: tridiagonal structure gives finite memory, two passes close the periodic loop exactly
- [x] p126-127: dimensionless parameter L — explained physical meaning (box/de Broglie ratio)
- [x] p132-133: GADGET-2 over Hydra — stated directly
- [ ] p139-141: interference/dispersion timescale — deferred: needs quantitative analysis
- [x] p168: GADGET velocity scaling — noted in Known Issues as possible source
- [x] p169: phase built from gravitational vs velocity potential — noted in Known Issues
- [x] p174-177: dedicated "Known Issues and Limitations" subsection — added

## Strengths to preserve
- [x] Periodic-boundary treatment — preserved
- [x] Widrow-Kaiser correction — preserved
- [x] Honesty about limitations — strengthened with dedicated subsection

## Metadata
- [ ] PDF metadata fields blank (deferred: will address when rebuilding LaTeX)
