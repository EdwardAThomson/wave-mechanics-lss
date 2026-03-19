# Chapter 3: Wave-mechanics — Checklist

## Copy-editing (Shared #9)
- [ ] p55: stray hyphen in "the- particle"
- [ ] p57: "it is may not be wise"
- [ ] p60: "a the fluid's pressure"
- [ ] p72: awkward "sesqui-linear" phrasing
- [ ] p78: repeated "an"
- [ ] p85: repeated "the"

## Copy-editing (Claude-only)
- [ ] p86: comma splice — "Such phenomena do not exist..., they are strictly forbidden."

## Mathematical errors (Shared)
- [ ] p53: density interpretation wrong — `|ψ|` should be `|ψ|²` (#3)
- [ ] Eq 3.32: definitely wrong — `e^{-i(K+V)dt} ≠ e^{-i(V+K)dt}` is meaningless; rewrite for non-commuting split operators (#4)
- [ ] Phase-angle method needs `atan2` qualification — `arctan(Im/Re)` is numerically unsafe (#5)

## Mathematical (Claude-only)
- [ ] p56-57 / Eqs 3.4-3.6: Madelung derivation — clarify division by ψ; check sign in Eq 3.6 vs Hamilton-Jacobi
- [ ] p56, p63: link phase conventions between Madelung (φ) and Bohm (S) explicitly
- [ ] p75 / Eq 3.27: distinguish symmetrised form from standard Strang splitting
- [ ] p75: units on `ν = ℏ/m` in ultralight-DM discussion — check dimensions
- [ ] p82 / Eq 3.41: check sign in `φ = -ν W⁻¹(arg ψ)` against Madelung convention

## Conceptual / substantive (Shared)
- [ ] Chapter needs better signposting; abrupt transition to literature review — add bridging paragraph (#6)
- [ ] Second-quantisation/QFT section longer than its role warrants — tighten (#7)
- [ ] Bohmian section needs clearer statement of relevance (#7)

## Conceptual (GPT-only)
- [ ] p58: "Free Particle Schrödinger = Fluid − Pressure" → "pressure-free" conclusion too strong
- [ ] p79: unitarity/symplecticity/momentum/energy conservation paragraph too compressed — rewrite
- [ ] p85-86: vorticity discussion overstates what Kelvin's theorem rules out — soften
- [ ] p62: quantum mechanics fails below 10⁻¹³ m — check or contextualise
- [ ] p60-61: `ℏ → 0` classical limit — "Schrödinger equation breaks down" is not careful
- [ ] p53: `ℏ` vs `ν = ℏ/m` shifts too loosely — frame explicitly from the start
- [ ] p60: pressure term matters only near shell crossing — too categorical
- [ ] p75: mention standard symmetric/Strang splitting explicitly
- [ ] p82: `φ`, `φ_v`, sign conventions doing too much work — add explicit convention statement

## Conceptual (Claude-only)
- [ ] p65: literature timeline — include year of this thesis for context
- [ ] p66-67: periodic boundaries in Widrow & Kaiser uncertain — state uncertainty explicitly
- [ ] Structure: literature-review subsections — consider compressing or tabulating

## Metadata
- [ ] PDF metadata fields blank
