# Chapter 5: Solving the Full Schrödinger-Poisson System — Checklist

## Copy-editing (Shared #11-12)
- [ ] Repeated articles, Widrow/Widow typos, awkward phrasing, colloquial wording — full pass
- [ ] Pronoun inconsistency (`I` vs `we`) — make consistent throughout

## Poisson-equation presentation — Major fix (Shared #4)
- [ ] Add mean-density subtraction clarification for periodic box
- [ ] Remove misleading standalone integral expression for the potential
- [ ] Replace `ψ_k ψ_k*` with Fourier transform of `|ψ|²`

## Mathematical errors (Shared)
- [ ] Eq 5.64: 1D Gaussian normalisation likely wrong (#5)
- [ ] Eq 5.69: not a correct total-momentum expression — fix (#6)
- [ ] Conservation discussion overstated — don't imply mass+momentum → energy conservation (#6)

## Major unresolved scientific issues (Shared)
- [ ] Velocity discrepancy with GADGET — reframe as major limitation, not minor caveat (#7)
- [ ] Spurious growing bulk flow — treat as important problem (#8)

## Conceptual / substantive (Shared)
- [ ] Chapter too long in places — tighter organisation (#10)
- [ ] Clearer signalling of what is new vs what is background
- [ ] Explicit separation between solid results and unresolved issues

## Conceptual (GPT-only)
- [ ] p113-114: basis discussion appears backwards (position vs momentum space) — fix
- [ ] p118-120: distinguish truly implicit scheme from semi-implicit/quasi-implicit implementation
- [ ] p149-150: top-hat discussion internally inconsistent about role of ν — reconcile
- [ ] p133: switch from Hydra (Ch4) to GADGET-2 here — motivate explicitly
- [ ] p153-156: separate resolution issues from deeper framework/implementation issues
- [ ] p155-164: "inevitable" / "fundamental" claims too strong — qualify
- [ ] p158-166: state whether Gaussian smoothing is diagnostic or part of physical comparison
- [ ] p174: no-vorticity claim — frame threshold/diagnostic logic more cautiously

## Conceptual (Claude-only)
- [ ] p113: "Cayley method" attribution unclear (Widrow & Kaiser vs Goldberg) — clarify
- [ ] p118: convention `m = 1/2` — justify briefly
- [ ] p119: reusing δ for timestep vs density contrast — change notation
- [ ] p115: Cayley factors commute — justify (both functions of same matrix)
- [ ] p116-117: Goldberg → Suzuki/Strang splitting — highlight as key methodological step
- [ ] p116: ad hoc `V/3` split — motivate or treat as failed intermediate idea
- [ ] p119-129: recursion relations — consider worked example or sweep-direction illustration
- [ ] p123: two boundary-recursion passes suffice — add convergence comment
- [ ] p126-127: dimensionless parameter L — explain physical meaning
- [ ] p132-133: reason for choosing GADGET-2 over Hydra — state directly
- [ ] p139-141: interference/dispersion timescale — make more quantitative
- [ ] p168: GADGET velocity scaling — double-check as possible mismatch source
- [ ] p169: possible explanation — phase built from gravitational vs velocity potential
- [ ] p174-177: add dedicated "Known Issues and Limitations" subsection

## Strengths to preserve
- [ ] Periodic-boundary treatment — one of the strongest contributions
- [ ] Widrow-Kaiser correction — important contribution
- [ ] Honesty about limitations — maintain this throughout revisions

## Metadata
- [ ] PDF metadata fields blank
