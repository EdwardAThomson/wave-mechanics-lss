# Chapter 4: Free Particle Approximation — Checklist

## Copy-editing (Shared #10)
- [x] p88: repeated "the the" → "the"
- [x] p95: "V_e is myself..." → "V_e is my result"
- [x] p97: "straight forward" → "straightforward"
- [x] p105: "In should be obvious" → "It should be obvious"
- [x] p106: "presences" → "presence"
- [x] p107-108: orphaned sentence fragment — joined with em dash
- [x] p109: "frame work" → "framework"

## Copy-editing (Shared #11)
- [ ] Pronoun inconsistency (`I` vs `we`) — deferred: pervasive throughout, needs careful pass

## Copy-editing (Claude-only)
- [x] p107: comma splice comparing FPA to Zel'dovich approximation — repaired

## Equation / presentation errors (Shared)
- [x] Eq 4.3: `f` is not defined (#3) — added f = d ln D/d ln a ≈ Ω_m^{0.6}
- [x] Eq 4.5: sign/presentation problem — removed spurious leading minus sign (#4)
- [x] Ch3/Ch4 velocity-potential sign convention — fixed line 97 to v = -∇φ_v (#5)
- [x] Eq 4.17: RMS formula wrong — added outer square root, fixed subscript (#6)
- [x] Eq 4.21: initial-potential formula — replaced ambiguous integral with Fourier/Green-function form (#7)

## Conceptual / substantive (Shared)
- [x] Speed/performance claims overstated — qualified to quasi-linear regime (#8)
- [x] Odd turn-over in velocity comparison — explained as phase-wrapping aliasing at velocity extremes (#9)

## Conceptual (GPT-only)
- [x] p87: "time evolution is unitary (that is, it conserves energy)" — separated: unitary conserves mass; exact solution also conserves energy
- [x] p95 / Fig 4.3: caption used ρ(V_e,V_c) but shows RMS — fixed caption notation
- [x] p101: non-monotonic correlation — added explanation (large-scale mode dominance)
- [x] p87-88: physical argument for zero effective potential — strengthened with explicit factor and regime
- [x] p90: role of growth factor D as effective time variable — stated more directly
- [x] p92: large ν prevents collapse — context already clear in text
- [x] p96-97: 3D extension simplicity specific to free-particle case — already stated in text
- [x] p104-105: near-identical velocity dispersions — emphasised as isotropy consistency check
- [N/A] p106-107: vorticity discussion — language already tentative

## Conceptual (Claude-only)
- [x] p88: multiplicative factor — tied to parameter regime and epoch
- [N/A] p89: `ρ = ψψ*` vs density contrast δ notation — both used correctly in context
- [x] p90 / Eq 4.6: same as Ch2 linear growth equation — added cross-reference
- [x] p90: incomplete Beta function growth factor — noted numerical integration
- [x] p91: `1 ≫ δ_a > 0` → `0 < δ_a ≪ 1`
- [x] p92: shell crossing at D=101 — stated as analytic (D = 1/δ_a from Zel'dovich)
- [x] p99: "jump to any D" feature — added paragraph emphasising as distinctive FPA advantage
- [x] p100: Ω_{b,0} = 0 — noted as CDM-only limitation
- [x] p101: correlation coefficients — context already clear (post-initial comparisons)
- [N/A] p106: `ψ = (0,0)` — kept; intentional emphasis on both complex components being zero
- [N/A] p107 / Fig 4.9: Gaussian caption — adequate as-is

## Metadata
- [ ] PDF metadata fields blank (deferred: will address when rebuilding LaTeX)
