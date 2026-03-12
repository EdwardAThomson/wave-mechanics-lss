# Chapter 3: Wave-mechanics - Review Notes

## Overall Assessment
This is the most intellectually ambitious chapter so far and arguably the most important — it lays the theoretical foundation for the entire thesis. It covers the Schrödinger-Poisson system, the Madelung transform, Bohmian mechanics, a comprehensive literature review (sections 3.2.1–3.2.10), initial conditions, velocity calculation methods, and vorticity. The writing is generally clear and the physical reasoning is sound. The literature review is thorough and well-organised. However, there are some mathematical issues in the Madelung derivation, a few typos, and some passages where the reasoning could be tightened.

---

## Typos and Grammatical Errors

1. **p53**: "|ψ| = ψψ* provides the amplitude of the wave which we interpret as density" — |ψ| = √(ψψ*), not ψψ*. The quantity ψψ* = |ψ|² is the density. Either say "|ψ|² = ψψ* provides... the density" or "|ψ| = √(ψψ*) provides the amplitude."

2. **p55**: "enumeration of the- particle number" — stray hyphen: "enumeration of the particle number"

3. **p55**: "shown to be reliable then this extra complication" → "reliable **then** this extra complication" should be "reliable, **then** this extra complication" or better: "reliable**,** this extra complication"

4. **p57**: "it is may not be wise to pick a reference frame" → remove "is": "it **may** not be wise"

5. **p60**: "break down as a the fluid's pressure" → remove "a": "break down as **the** fluid's pressure"

6. **p72**: "the sesqui-linear quantity ψψ*" — ψψ* is bilinear (or sesquilinear in the sense of complex conjugation), but calling it "sesqui-linear" without context may confuse. It's more standard to call ψψ* the "modulus squared" or "norm squared" of ψ.

7. **p75**: "Prima facie" is correctly used but spelled "Prima facie" — correct.

8. **p78**: "an an implicit method" → double "an": "**an** implicit method"

9. **p83**: "tan(φ) = b/a" — the standard definition of the argument of c = a + ib is tan(φ) = b/a, but this gives the correct angle only when a > 0. Should mention atan2 or note the quadrant ambiguity explicitly, as this is computationally important.

10. **p83**: "would should a singularity" → "**would show** a singularity" (typo)

11. **p85**: "the the probability current" → double "the"

12. **p86**: "Such phenomena do not exist in current simulations, they are strictly forbidden." — comma splice; should be a semicolon or two sentences.

---

## Clarity and Content Issues

### Section 3.1 — Interpretation

13. **p53**: The Schrödinger equation (3.1) uses ℏ but the text says "ℏ sets the limit of spatial resolution, it can be considered as the unit size of a grid cell in phase space (ΔxΔp ~ ℏ)." This is the uncertainty principle relation, but in the classical wave-mechanics context of this thesis, ℏ is not Planck's constant — it's an effective parameter ν = ℏ/m. The conflation of the quantum interpretation (phase space cell) with the classical one (diffusion coefficient) starts here but isn't fully resolved until p60-61. Consider flagging this ambiguity earlier.

14. **p53-54**: The discussion about second quantisation and QFT (p54-55) is interesting but quite long for what amounts to "we don't need QFT for this thesis." The reader may wonder why so much space is devoted to something explicitly declared out of scope.

15. **p55**: "Dark Matter particles, whatever they are, are quantum in nature, so working with a Schrödinger equation (albeit in the classical limit of a large occupation number) seems like a natural approach to take." — This is a key philosophical claim but is stated without much justification. Why is the Schrödinger equation more natural than, say, classical field equations? The argument rests on the Madelung equivalence, which comes later. Consider reordering or at least forward-referencing.

### Section 3.1.2 — Madelung / Hydrodynamic form

16. **p56**: Equation 3.3: ψ = αe^{iφ/ν} — the text defines α = √ρ, ν = ℏ/m. But the phase convention here differs from the Bohm form on p63 (Eq 3.11): ψ = R e^{iS/ℏ}. In the Madelung form, v = ∇φ, while in Bohm's form v = ∇S/m. These are consistent if φ = -S/m (with appropriate sign conventions), but the sign conventions should be made explicit to avoid confusion.

17. **p56-57**: Equations 3.4-3.5 — the derivation of LHS and RHS of the transformed Schrödinger equation. I note that in Eq 3.4:
   LHS = iν(∂ψ/∂t) = iνψ/α · ∂α/∂t - ∂φ/∂t
   This should have a ψ multiplying the whole expression, since LHS = iν(∂ψ/∂t) and we need to match the RHS which has ψ throughout. The derivation appears correct but the presentation could be clearer about what's been divided through by ψ and what hasn't.

18. **p57**: Equation 3.6: Vψ = -∂φ/∂t · ψ - ½ψ(∇φ)² — this is the Bernoulli equation multiplied by ψ. But the sign might need checking: from the standard Madelung derivation, the Hamilton-Jacobi equation gives ∂φ/∂t + ½(∇φ)² + V + Q = 0 where Q is the quantum pressure. The signs in 3.6 should be consistent with this.

19. **p58**: Equation 3.8: ∂α²/∂t + ∇·(α²∇φ) = 0 — correct continuity equation.

20. **p58-59**: The discussion about reconciling the pressure-free Schrödinger equation with the fluid equations (which have pressure) is one of the most important conceptual points in the thesis. The key insight — "Free Particle Schrödinger = Fluid − Pressure" — is well stated. However, the logical chain could be tighter. Currently it reads somewhat discursively; a more structured presentation (perhaps numbered steps) would help.

### Section 3.1.2.1 — Bohmian mechanics

21. **p61-63**: The Bohmian mechanics section is well-written as a philosophical/historical overview but its connection to the practical work of the thesis could be made more explicit. What does the reader *need* from this section? The answer is equation 3.14 (quantum potential = quantum pressure). Consider stating this upfront.

### Section 3.2 — Literature review

22. **p65**: The timeline is very useful. Minor note: it would benefit from including the year of your own thesis (2011) at the end to show where your work fits.

23. **p66-67**: The Widrow & Kaiser review is good. One note: "it is unclear if their code implements periodic boundary conditions" — this seems like something that could be determined from the paper. If it genuinely can't be determined, that's worth emphasising as a gap.

24. **p72**: Equation 3.20: ψ(r,t) = ψ₀(a/a₀)^{-3/2} e^{A(r,t)+iB(r,t)/ℏ} — the (a/a₀)^{-3/2} factor accounts for dilution due to expansion (mass conservation in comoving volume). Good.

25. **p75**: Equation 3.27: e^{-i(K+V)dt} ≈ ½[e^{-iKdt}e^{-iVdt} + e^{-iVdt}e^{-iKdt}] — this is the symmetrised first-order splitting (Strang-like but not quite Strang). A standard Strang splitting would be e^{-iVdt/2}e^{-iKdt}e^{-iVdt/2}. Worth noting the distinction and which is more accurate.

26. **p75**: "ν = ℏ/m ~ 10⁻¹⁵/10⁻²² is of order 10⁷ (eV·seconds)" — the units here are confusing. ℏ has units of J·s or eV·s, and m has units of eV/c² (in natural units) or kg. The ratio ℏ/m has units of m²/s (in SI) or length²/time. The stated "eV·seconds" doesn't seem dimensionally correct. Should be checked.

### Sections 3.3-3.5 — Solving, initial conditions, velocity, vorticity

27. **p78**: Equation 3.32: e^{-i(K+V)dt} ≠ e^{-i(V+K)dt} — but this is just rewriting the same expression (addition is commutative!). What you mean is e^{-iKdt}e^{-iVdt} ≠ e^{-iVdt}e^{-iKdt} (the product of exponentials of non-commuting operators). The equation as stated is trivially wrong — both sides are identical since K+V = V+K. This needs correction.

28. **p81**: Equations 3.35-3.38 (BBKS transfer function) — these are standard and appear correct. The notation switches between |k| (bold, underlined) and k somewhat freely.

29. **p82**: Equation 3.41: φ = -νW⁻¹(arg(ψ)) — the sign convention here should be consistent with the Madelung form (Eq 3.3) where v = ∇φ. If ψ = αe^{iφ/ν} then arg(ψ) = φ/ν, so φ = ν·arg(ψ). The negative sign and the unwrapping operator W⁻¹ need to be carefully justified.

30. **p84**: Equations 3.45-3.52 (probability current method) — this is a nice contribution. The derivation is clear and the final result (Eq 3.52) v = (ℏ/m)ℑ(∇ψ/ψ) is elegant and correct.

31. **p85-86**: Section 3.5 on vorticity singularities is conceptually important and sets up Chapter 7 well. The identification of ψ = 0 points as potential vortex locations is correct and physically significant.

---

## Equations Check

- **Eq 3.1**: iℏ ∂ψ/∂t = (-ℏ²/2m ∇² + mV)ψ — correct
- **Eq 3.2**: ∇²V = 4πGψψ* — correct (Poisson equation with density = |ψ|²)
- **Eq 3.3**: ψ = αe^{iφ/ν} — correct (Madelung ansatz)
- **Eq 3.8**: ∂α²/∂t + ∇·(α²∇φ) = 0 — correct (continuity)
- **Eq 3.9**: Pressure definition — form looks correct
- **Eq 3.11**: ψ = R e^{iS/ℏ} — correct (Bohm form)
- **Eq 3.12**: -∂ρ/∂t = ∇·(ρ∇S/m) — correct (continuity)
- **Eq 3.13**: ∂S/∂t = -[V + (∇S)²/2m - ℏ²∇²R/(2mR)] — correct (quantum HJ)
- **Eq 3.14**: Quantum potential = quantum pressure equivalence — correct
- **Eq 3.30-3.32**: Time evolution via exponential — correct in principle but see note #27 on 3.32
- **Eq 3.45**: Probability current — correct
- **Eq 3.52**: v = (ℏ/m)ℑ(∇ψ/ψ) — correct

---

## Structural Comments

32. This is the longest and most complex of the review chapters. The structure works well overall: theory → interpretation → literature → practical methods.

33. The literature review (3.2) is comprehensive and valuable as a standalone reference. Each subsection clearly describes what each group did and how it relates to this thesis.

34. The transition from the interpretive material (3.1) to the literature review (3.2) is slightly abrupt. A bridging sentence would help.

35. The chapter covers a lot of ground — interpretation, Madelung, Bohm, 10 literature reviews, solving methods, initial conditions, three velocity methods, and vorticity. While all of this material is relevant, the chapter is dense and a reader new to the field might struggle. Consider whether the literature review subsections could be more concisely presented (perhaps summarising key differences in a table).
