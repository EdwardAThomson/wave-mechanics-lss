# Chapter 3: Wave-mechanics - Reconciled Review Notes

## Source Coverage

- GPT review: `gpt-review/notes_ch03_wave_mechanics.md`
- Claude review: `claude-review/notes_ch03_wave_mechanics.md`
- Overall result: strong agreement, with GPT focusing more on overstatement and Claude adding a few extra equation-level checks

## Shared Findings

The two reviews independently agree on the following main points.

1. **This is one of the strongest and most important chapters in the thesis.** Both reviews treat it as conceptually central and commend its scope, ambition, and literature coverage.

2. **The literature review is a major strength.** Both reviews say the Chapter 3 survey is valuable in its own right and well-organised, even if it is dense.

3. **The density interpretation on `p53` is wrong as written.** Both reviews independently flag that `ψ ψ^* = |ψ|^2`, not `|ψ|`.

4. **Equation 3.32 is definitely wrong as printed.** Both reviews note that
   `e^{-i(K+V)dt} != e^{-i(V+K)dt}`
   is meaningless because `K+V = V+K` as an algebraic sum. The intended statement is about non-commuting split operators.

5. **The phase-angle method needs quadrant care.** Both reviews independently say that using `arctan(Im/Re)` without an `atan2`-type qualification is numerically unsafe.

6. **The chapter is dense and needs better signposting.** Both reviews say the transition from the interpretation material to the literature review is somewhat abrupt and would benefit from a bridging paragraph.

7. **Some digressive material could be tightened.** Both reviews think the second-quantisation/QFT section is longer than its actual role in the thesis, and both say the Bohmian section needs a clearer statement of why it matters for the rest of the chapter.

8. **The probability-current velocity section is one of the strongest parts of the chapter.** Both reviews treat the practical velocity discussion as especially successful and clear.

9. **There are a number of straightforward copy-editing issues.** Both reviews independently flag items such as:
   `p55` stray hyphen in "the- particle",
   `p57` "it is may not be wise",
   `p60` "a the fluid's pressure",
   `p72` awkward "sesqui-linear" phrasing,
   `p78` repeated "an",
   `p85` repeated "the".

## GPT-Only Additions Worth Keeping

These points appear in the GPT review but not the Claude review and are worth preserving.

1. **`p58`**: The conclusion "Free Particle Schrödinger = Fluid - Pressure" followed by "Conclusively, the Schrödinger-Poisson system is pressure-free" is too strong. The chapter should distinguish more carefully between the exact Madelung form and regimes where the quantum-pressure term is neglected or becomes small.

2. **`p79`**: The paragraph linking unitarity, symplecticity, momentum conservation, and energy conservation is too compressed and should be rewritten more carefully.

3. **`p85-86`**: The vorticity discussion overstates what Kelvin's circulation theorem rules out. The wording around angular momentum, irrotational flow, shell crossing, and vortex diagnostics should be softened.

4. **`p62`**: The statement that quantum mechanics fails below `10^-13 m` or agrees less well there is historically loaded and needs checking or contextualisation.

5. **`p60-61`**: The classical-limit discussion around `ℏ → 0` is too blunt; "the Schrödinger equation breaks down" is not a careful way to describe the semiclassical limit.

6. **`p53`**: The discussion of `ℏ` versus `ν = ℏ/m` shifts between quantum and effective-classical interpretations too loosely and should be framed more explicitly from the start.

7. **`p60`**: The claim that the pressure term matters only near shell crossing is probably too categorical unless established more carefully.

8. **`p75`**: It would help to mention standard symmetric/Strang splitting explicitly when discussing operator splitting methods.

9. **`p82`**: The notation around `φ`, `φ_v`, and sign conventions is doing too much work and would benefit from one explicit convention statement.

10. **Metadata**: the PDF metadata fields are blank here as well.

## Claude-Only Additions Worth Keeping

These points appear in the Claude review but not the GPT review and are worth preserving.

1. **`p56-57` / Eqs. 3.4-3.6**: The Madelung derivation presentation is not fully transparent about what has been divided by `ψ`, and the sign structure in Eq. 3.6 should be checked against the standard Hamilton-Jacobi form.

2. **`p56` and `p63`**: The phase conventions in the Madelung and Bohm forms should be linked explicitly so the relation between `φ` and `S` is not left implicit.

3. **`p65`**: The literature timeline would benefit from including the year of this thesis so the reader can see where the work sits historically.

4. **`p66-67`**: If periodic boundaries in Widrow & Kaiser truly cannot be determined from the paper, that uncertainty should be stated more explicitly.

5. **`p75` / Eq. 3.27**: The split-operator expression should distinguish more clearly between the symmetrised form shown there and standard Strang splitting.

6. **`p75`**: The units attached to `ν = ℏ/m` in the ultralight-dark-matter discussion look dimensionally suspect and should be checked.

7. **`p82` / Eq. 3.41**: The sign in `φ = - ν W^{-1}(arg ψ)` should be checked carefully against the Madelung convention.

8. **`p86`**: "Such phenomena do not exist in current simulations, they are strictly forbidden." is a comma splice and should be broken up even apart from the conceptual qualification.

9. **Structure**: The literature-review subsections might be made more digestible by compressing them or tabulating key differences between authors/methods.

## Reconciled Priorities

1. Fix the clear mathematical statements first:
   `|ψ|` versus `|ψ|^2`,
   Eq. 3.32,
   the Madelung/Bohm sign conventions,
   and the practical phase-angle/quadrant issue.

2. Rework the pressure discussion so the chapter does not overclaim that the Schrödinger-Poisson system is intrinsically pressure-free.

3. Tighten the operator-splitting and conservation-law discussion so unitarity, symplecticity, norm conservation, and long-term energy behaviour are not blurred together.

4. Soften the vorticity/angular-momentum claims and make the role of singular points of `ψ` more precise.

5. Trim or better frame the digressive sections:
   second quantisation/QFT,
   Bohmian mechanics,
   and the more speculative interpretive comments.

## Bottom Line

Chapter 3 is a strong and important chapter, and both reviews treat it as one of the thesis's best sections. The reconciled verdict is that it mostly needs targeted mathematical and interpretive tightening rather than structural overhaul, with the most important fixes concentrated in the Madelung/pressure discussion, operator splitting, and the vorticity claims.
