# Chapter 1: Introduction - Review Notes

## Overall Assessment
The chapter is well-structured and readable. It serves its purpose as a broad introduction to the concordance model and motivates the wave-mechanical approach to LSS. The narrative flows logically from LSS overview → concordance model → spacetime/GR → Friedmann equation → cosmic history → CMB → dark matter → dark energy → fluid equations.

However, there are a number of typos, grammatical issues, and some passages where clarity could be improved.

---

## Typos and Grammatical Errors

1. **p3**: "has suggests a way" → should be "has suggested a way"

2. **p3**: "Penrose Interpretation of quantum mechanics" — capitalisation inconsistent. Should this be "Penrose interpretation" (lowercase 'i') for consistency?

3. **p12**: "Einstein's principle of equivalence" — what is described (physical states being the same in different frames) is closer to the principle of *general covariance* than the equivalence principle (which is about the equivalence of gravitational and inertial mass / local equivalence to flat spacetime). This could confuse a knowledgeable reader.

4. **p12**: "the notion of absolute objectivity (measurement of invariant quantities) is applied, not to translations but, to rotations" — the comma placement is awkward. Suggest: "is applied not to translations but to rotations."

5. **p13**: "using the newly discovered idea of using the metrics *g* as potentials" — the word "using" appears twice. Also, should be "metric" (singular) or "metric components *g_μν*".

6. **p14**: "As Universe was thought to be static" → "As **the** Universe was thought to be static"

7. **p15**: "the only unfamiliar term here is *k* the spatial curvature" → missing comma or colon: "the only unfamiliar term here is *k*, the spatial curvature"

8. **p15**: "the two American physicists that discovered it independently of Friedmann" — Robertson and Walker were not both American. Howard Percy Robertson was American, but Arthur Geoffrey Walker was British.

9. **p15**: "energy and momentum conservation (*T_{μν;ν}* = 0)" — this is the covariant divergence of the energy-momentum tensor. The notation is correct but might benefit from a brief word: "the covariant conservation law."

10. **p17**: "that lead to the notion" → "that **led** to the notion" (past tense)

11. **p19**: "the concept of temperature is not valid in this regime as temperature is only defined for a system under Maxwell-Boltzmann statistics" — this is not quite right. Temperature is defined for any system in thermal equilibrium, including those described by Fermi-Dirac or Bose-Einstein statistics. The issue at the Big Bang singularity is that the concept of thermal equilibrium itself breaks down, not that Maxwell-Boltzmann statistics are required.

12. **p21**: "for the symmetry breaking required to produced the matter-antimatter asymmetry" → "to **produce**"

13. **p21**: "GUT is an acronym for Grand Unification Theory, this theory has proposes that" → "this theory **proposes** that" (remove "has")

14. **p22**: "the the Universe would be in causal contact" → double "the"

15. **p23**: "Also at this time, all free electrons are thought to have been captured by the protons and Helium nuclei" — "Helium" should be lowercase "helium" (element names are not capitalised in standard usage).

16. **p23**: "physical processes that occurred at Last Scattering.." → double period

17. **p24**: "Doppler shits in frequency" → "Doppler **shifts**" (significant typo!)

18. **p24**: "all futures perturbations" → "all **future** perturbations"

19. **p25**: "an adiabatic component and and isocurvature component" → double "and"

20. **p29**: Table of cosmological parameters: Ω_k range shows "0.0179 < Ω_k < 0.0081" — this range is inverted (the lower bound is larger than the upper bound). Should likely be "-0.0179 < Ω_k < 0.0081" or similar, with a negative sign.

21. **p30**: "faster than he expected" — the "he" is ambiguous; the preceding sentence discusses Rubin and Ford (two people). Should be "they expected".

---

## Clarity and Content Issues

22. **p2 (opening)**: "the lack of data has often meant that cosmology was a 'playground' for theoreticians" — this is a nice point but slightly undermined by the next sentence switching to past tense ("We had to make assumptions"). The tense shifts between past and present throughout the opening paragraphs could be smoothed.

23. **p3**: The paragraph about Penrose and quantum gravity feels somewhat tangential. The connection between the Schrödinger-Poisson system used in LSS simulations and Penrose's interpretation of quantum mechanics is interesting but is introduced without much motivation and then dropped. A sentence explicitly stating "this connection is not pursued in this thesis" (or the opposite) would help.

24. **p8**: "In wave-mechanics there is no 'number of particles' so mass resolution is limited by machine precision" — this is an important and interesting claim that deserves a bit more explanation. How exactly does machine precision limit mass resolution in a wave-mechanical code? Is it about the dynamic range of the wavefunction amplitude?

25. **p10**: "It is not be obvious that the Schrödinger formalism can be applied to Large Scale Structure" — typo: "not be obvious" → "not obvious" (remove "be"). Also, the sentence about performing numerical experiments on reflection and tunnelling and omitting them "for brevity" is tantalising — a footnote or forward reference to where these might be found would be helpful.

26. **p12**: The discussion of Einstein's principle of equivalence, coordinate dependence of inertial forces, and the distinction between rotational and translational invariance is somewhat compressed and may confuse readers. The logical flow from "frames of reference describe the same point" to "inertial forces are coordinate dependent" to "absolute objectivity applies to rotations not translations" is not obvious without more connective tissue.

27. **p13**: Equation 1.2 uses Einstein summation convention (no explicit sum) while equation 1.1 has an explicit sum. The text says for 1.1 "the summation is given explicitly (not by convention)" but then 1.2 silently switches. This is fine if deliberate, but worth noting for consistency.

28. **p14-15**: The scale factor is introduced on p14 as the ratio of physical to comoving distance, with r = a x. The notation uses underlined r and x. But equation 1.4 uses ξ² for the 3-space metric, which then gets rewritten in polar coords in 1.6. The connection between ξ and x (comoving coordinate) could be made more explicit.

29. **p15**: "the two American physicists" — as noted above, Walker was British, not American.

30. **p17**: The discussion about expansion of space being "contentious" (p17-18) is philosophically interesting but quite long for an introduction chapter. It somewhat disrupts the flow of the concordance model narrative. Consider whether this belongs in a footnote or can be condensed.

31. **p19**: The claim that temperature "is only defined for a system under Maxwell-Boltzmann statistics" is incorrect — see typo note #11. Temperature is well-defined for quantum statistical systems too. The valid point here is that near the singularity, the concept of equilibrium (and hence temperature) may not apply.

32. **p24**: Eq 1.11: P_k = <δ_k δ*_{k'}> = Ak^n — the text says "power in Fourier space P_k is assumed to be linearly proportional to the wavenumbers k." But P = Ak^n with n ≈ 1 means P is proportional to k (linear in k), while the spectrum is a *power law* in k. The word "linearly" here is technically correct only for n = 1 but might mislead readers into thinking it's always linear. Suggest clarifying that n ≈ 1 (Harrison-Zel'dovich) gives approximately linear proportionality.

33. **p28**: Equation 1.12: δT/T = -δz/(1+z) = δρ/ρ — the sign convention here may need checking. The Sachs-Wolfe effect gives δT/T = -(1/3)δΦ/c² for large scales, while the intrinsic (adiabatic) relation involves δT/T ∝ δρ/ρ. The equality as written conflates the relationship between redshift perturbations and density perturbations in a way that might not be precisely correct — it's a crude approximation as stated, which is acknowledged in the text, but worth double-checking the signs and factors.

34. **p34**: Equations 1.14 and 1.15 (continuity and Euler equations) — these are presented without derivation, which is fine for an intro, but the notation switches from the Boltzmann equation's (x, v_x, v_y, v_z) to (ρ, v). A brief sentence connecting these would help: "Taking moments of the collisionless Boltzmann equation yields..."

---

## Structural / Flow Comments

35. **Section 1.2.1 "A briefer history of time"** — clever title, but this subsection covers a lot of ground (Big Bang, inflation, CMB, Last Scattering) in a narrative style before section 1.2.2 introduces the mathematical framework. The overlap with the later more detailed subsections (CMB, Dark Matter, Dark Energy, etc.) means some material is covered twice. Consider whether the narrative overview and the detailed subsections could be better integrated to reduce repetition.

36. **Chapter ending (p35)**: The chapter ends rather abruptly with "The next chapter provides a review of the standard simulation techniques in LSS." A brief summary paragraph recapping the key points of Chapter 1 and how they set up the rest of the thesis would provide a smoother transition.

---

## Equations Check

- **Eq 1.1**: Q = Σ g_μν dx_μ dx_ν — correct (explicit sum form of metric)
- **Eq 1.2**: ds² = g_μν dx_μ dx_ν — correct (Einstein convention)
- **Eq 1.3**: G_νμ + Λ g_νμ = (8πG/c⁴) T_νμ — correct (Einstein field equations)
- **Eq 1.4**: ds² = c²dt² - a(t)²dξ² — correct (FRW metric, compact form)
- **Eq 1.5**: r = a x — correct (physical = scale factor × comoving)
- **Eq 1.6**: dξ² = dr²/(1-kr²) + r²[dθ² + sin²θ dφ²] — correct (3-metric in polar coords)
- **Eq 1.7**: ȧ² = (8πGρ/3)a² - kc² — correct (Friedmann equation)
- **Eq 1.8**: Ω = ρ/ρ_c = 8πGρ/(3H²) — correct
- **Eq 1.9**: H² = H₀² [Ω_r₀(1+z)⁴ + Ω_m₀(1+z)³ + Ω_k₀(1+z)² + Ω_Λ₀] — correct
- **Eq 1.10**: v = H₀d — correct (Hubble's law)
- **Eq 1.11**: P_k = <δ_k δ*_{k'}> = Ak^n — correct (power spectrum definition)
- **Eq 1.12**: δT/T = -δz/(1+z) = δρ/ρ — approximately correct but see note #33
- **Eq 1.13**: Boltzmann equation — correct
- **Eq 1.14**: ∂ρ/∂t + (1/a)∇·(ρv) = 0 — correct (continuity, comoving)
- **Eq 1.15**: ∂v/∂t + (1/a)v·∇v + Hv = -(1/a)∇Φ — correct (Euler, comoving)
