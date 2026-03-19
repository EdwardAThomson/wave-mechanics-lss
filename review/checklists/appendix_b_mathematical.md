# Appendix B: Mathematical Appendix — Checklist

## Copy-editing (Shared #8)
- [ ] Opening sentence garbled — rewrite ("equivalence of the Schrödinger equation and the continuity equation")
- [ ] "imaginary party" typo
- [ ] B.2 final QED-style statement too abrupt

## Mathematical errors (Resolved disagreement #1)
- [ ] B.2 phase convention / velocity formula sign-inconsistency:
  - Eq B.19 defines `ψ_i = (1+δ_i)^{1/2} e^{-iφ_i/ν}` (negative phase)
  - Standard probability current → `v_i = -∇φ_i`
  - But Eq B.33 states `v_i = ∇φ_i`
  - Eqs B.18 and B.34-B.35 not sign-consistent
  - Intermediate algebra B.29-B.30 also has sign flip
  - **Needs correction, not just clarification**

## Conceptual (GPT-only)
- [ ] B.1 conceptually overstated — not just continuity equation, it's the full Madelung pair
- [ ] p230 / Eqs B.8-B.10: reusing `V` for Bernoulli-type quantity and external potential — confusing
- [ ] p230-231: derivation omits pressure and potential to isolate continuity — state more clearly
- [ ] p235: `L` used for both dimensionless parameter and box size — clarify

## Conceptual (Claude-only)
- [ ] B.3: add brief reminder of what χ, U, y, L mean (inherited from Ch5)
- [ ] B.3 ends too abruptly after Eq B.46 — connect back to Ch5
- [ ] Add short appendix-level introduction (one sentence per section stating what it supports)

## Shared
- [ ] Widrow & Kaiser mistake — state explicitly what the mistake is (#7)
- [ ] Better signposting — clearer what each section is proving, what notation is inherited (#6)

## Metadata
- [ ] PDF metadata fields blank
