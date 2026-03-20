# Chapter 3: Wave-mechanics — Checklist

## Copy-editing (Shared #9)
- [x] p55: stray hyphen in "the- particle"
- [x] p57: "it is may not be wise"
- [x] p60: "a the fluid's pressure"
- [x] p72: awkward "sesqui-linear" phrasing — replaced with |ψ|² = ψψ*
- [x] p78: repeated "an"
- [x] p85: repeated "the"
- [x] p138: "can not known be exactly known" → "cannot be exactly known"

## Copy-editing (Claude-only)
- [x] p86: comma splice — "Such phenomena do not exist..., they are strictly forbidden." → semicolon

## Mathematical errors (Shared)
- [x] p53: density interpretation wrong — `|ψ|` should be `|ψ|²` (#3)
- [x] Eq 3.32: definitely wrong — `e^{-i(K+V)dt} ≠ e^{-i(V+K)dt}` is meaningless; rewritten for non-commuting split operators (#4)
- [x] Phase-angle method needs `atan2` qualification — `arctan(Im/Re)` is numerically unsafe (#5)

## Mathematical (Claude-only)
- [x] p56-57 / Eqs 3.4-3.6: Madelung derivation — clarified division by ψ (valid where ψ ≠ 0, vortices discussed later)
- [x] p56, p63: link phase conventions between Madelung (φ) and Bohm (S) explicitly — added α=R, φ=S/m identification
- [x] p75 / Eq 3.27: distinguish symmetrised form from standard Strang splitting — added note after equation
- [ ] p75: units on `ν = ℏ/m` in ultralight-DM discussion — check dimensions (deferred: needs careful unit analysis)
- [x] p82 / Eq 3.41: sign conventions — added explicit convention statement at start of Phase Unwrapping section

## Conceptual / substantive (Shared)
- [x] Chapter needs better signposting; abrupt transition to literature review — added bridging paragraph (#6)
- [N/A] Second-quantisation/QFT section longer than its role warrants — kept as-is; provides necessary context (#7)
- [x] Bohmian section needs clearer statement of relevance — contextualised Bohm's 10⁻¹³ m claim; linked conventions (#7)

## Conceptual (GPT-only)
- [N/A] p58: "Free Particle Schrödinger = Fluid − Pressure" — kept as-is; mathematically correct and qualified in next paragraph
- [x] p79: unitarity/symplecticity/momentum/energy conservation paragraph too compressed — rewritten with Noether, shadow Hamiltonian
- [x] p85-86: vorticity discussion overstates what Kelvin's theorem rules out — softened to "generation from irrotational conditions"
- [x] p62: quantum mechanics fails below 10⁻¹³ m — contextualised as Bohm's era; noted modern experiments confirm QM further
- [x] p60-61: `ℏ → 0` classical limit — "breaks down" replaced with "singular: removes highest-order derivative"
- [x] p53: `ℏ` vs `ν = ℏ/m` shifts too loosely — added explicit ν definition and role statement at chapter opening
- [N/A] p60: pressure term matters only near shell crossing — kept; attributed to Short's result, qualified by context
- [x] p75: mention standard symmetric/Strang splitting explicitly — added after Woo & Chiueh equation
- [x] p82: `φ`, `φ_v`, sign conventions doing too much work — added explicit convention statement

## Conceptual (Claude-only)
- [x] p65: literature timeline — added "2011 - Thomson - this thesis"
- [N/A] p66-67: periodic boundaries in Widrow & Kaiser uncertain — already states "it is unclear"; no change needed
- [N/A] Structure: literature-review subsections — kept as-is; chronological structure serves the narrative

## Metadata
- [ ] PDF metadata fields blank (deferred: will address when rebuilding LaTeX)
