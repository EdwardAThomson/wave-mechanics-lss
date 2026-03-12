# Chapter 4: Free Particle Approximation - Review Notes

## Overall Assessment
A well-structured chapter (24 pages) that serves as a bridge between the theoretical framework of Chapter 3 and the full Schrödinger-Poisson system of Chapter 5. It covers the FPA theory, 1D verification against Short's results, 3D toy model extension, 3D cosmological comparison with the Hydra N-body code, consistency checks, vorticity detection, and evaluation. The progression from 1D → 3D toy → 3D cosmological is logical and builds confidence incrementally. The results are clearly presented with good figures and statistical comparisons. However, there are some typos, a few mathematical presentation issues, and some places where the writing could be tighter.

---

## Typos and Grammatical Errors

1. **p88**: "the the velocity potential φ_v" → double "the": "**the** velocity potential φ_v"

2. **p88**: "the multiplicative factor is close to one" — the factor is 3Ω_cdm/(2f²D). For this to be close to one requires specific conditions on Ω_cdm, f, and D. It would help to state explicitly for what parameter values and at what epoch this approximation holds well.

3. **p89**: "the density field (ρ = ψψ*)" — throughout the thesis, ρ sometimes means total density and sometimes means the density contrast field. Here ρ = ψψ* = |ψ|² is the density (not the density contrast δ). Consistent notation would help.

4. **p90**: "here H is the familiar Hubble parameter that is statement of expansion" → grammatically awkward: "that is **a** statement of **the** expansion" or better: "which describes the rate of expansion"

5. **p90**: "and a is the expansion factor" — this is clear but earlier (p88) the scale factor was introduced without being named. Minor inconsistency.

6. **p92**: "The Γ parameter in Short's thesis is not used explicitly here but it is the combination of the effective Planck constant and the period: Γ = ν/p²" — this is presented as an aside mid-paragraph. It would read better as a footnote or set off more clearly.

7. **p95**: "V_e is myself and V_c is Chris Short" — should be "V_e is **mine** and V_c is Chris Short's" or "V_e is from my code and V_c is from Short's code"

8. **p96**: "Zel'dovich" is spelled "Zel'dovich" here but "Zel'dovich" on p88. The apostrophe/accent varies throughout — should be consistent.

9. **p97**: "straight forward" → should be one word: "**straightforward**"

10. **p100**: "v 4.2.1" — appears to be a typo for section reference; likely should be "(v. 4.2.1)" or "(version 4.2.1)" referring to the Hydra version number. The formatting is unclear.

11. **p105**: "The slight skewing of the distributions is a cause of some concern but the overall behaviour of the simulations are at least believable" → subject-verb agreement: "the overall behaviour... **is** at least believable"

12. **p105**: "hence we do not suspect that something is terribly wrong" — this is rather colloquial for a thesis. Consider: "which suggests the simulation is behaving correctly"

13. **p106**: "the presences of vortices" → "the **presence** of vortices"

14. **p107**: "see the rest of the velocity vectors 'circling' around the wavefunctions null-point" — this sentence fragment appears to be an orphan, disconnected from any preceding text. It reads as if the beginning of the sentence was accidentally deleted. Should be "...we can see the rest of the velocity vectors..."

15. **p108**: "Schrödinger equation." — this appears as an orphan sentence fragment at the top of a paragraph. It looks like it should be joined to the preceding paragraph about the wavefunction construction artefact: "...is an artefact from the construction of the wavefunction **using the** Schrödinger equation."

16. **p109**: "frame work" → should be one word: "**framework**"

---

## Clarity and Content Issues

### Section 4.0 — FPA Theory

17. **p87-88**: The explanation of why V = Φ_g - φ_v = 0 is the core physical argument of the FPA but it's somewhat buried. The key insight — that in the linear regime the gravitational potential and velocity potential are approximately equal — deserves more prominent treatment. Consider a display equation showing V = Φ_g - (3Ω_cdm)/(2f²D) φ_v ≈ 0 with explicit conditions.

18. **p88**: Equation 4.3: V = Φ_g - (3Ω_c)/(2f²D) φ_v = 0 — the factor involves f (the growth rate f = d ln D/d ln a) but f is not defined here. It appears implicitly through the relationship between the potentials. Should be defined.

19. **p89**: Equation 4.4: v = -∇φ_v — the sign convention here (velocity is minus the gradient of the velocity potential) is the opposite of Chapter 3's convention where v = ∇φ (Eq 3.8). This needs reconciliation. If φ_v here corresponds to -φ in Chapter 3, this should be stated explicitly.

20. **p89**: Equation 4.5: ψ = -(1/(2π)³) ∫ ψ̂_init(k) e^(-iν(D-1)k²/2) e^(ik·x) d³k — the leading minus sign is unusual for an inverse Fourier transform. Standard convention gives a positive prefactor. Either the minus sign is a convention choice that should be explained, or it may be a typo.

21. **p89**: "In should be obvious" → "**It** should be obvious" (typo)

### Section 4.0.1 — Linear Growth Factor

22. **p90**: Equation 4.6: ∂²δ/∂t² + 2H ∂δ/∂t - 4πGρ̄δ = 0 — this is the standard linear perturbation equation. However, note that this is identical to Eq 2.6 from Chapter 2. A cross-reference would be helpful rather than re-deriving it.

23. **p90**: Equation 4.9 uses the incomplete Beta function β_α(5/6, 2/3). This is a relatively obscure special function. While it's fine to state the result, it would be helpful to note that for practical computation, numerical integration of Eq 4.7 may be simpler.

24. **p90**: "the integral for the growing mode in the case of a flat Universe (Ω_cdm + Ω_Λ = 1) can be expressed as:" — this should note that this is an Einstein-de Sitter generalisation to include Λ. For a pure Einstein-de Sitter universe (Ω_Λ = 0), D_+ = a, which is a useful limiting case to mention.

### Section 4.1 — 1D FPA

25. **p91**: Equation 4.11: φ_{v,init} = -(p/2π)² δ_i — this relation between velocity potential and density assumes the linear regime Poisson equation. Worth noting this explicitly: in the linear regime, ∇²φ ∝ δ, so for a single Fourier mode with wavenumber 2π/p, φ ∝ δ/k² = δ(p/2π)².

26. **p91**: "1 ≫ δ_a > 0" — this notation is non-standard. Usually written as "0 < δ_a ≪ 1" (i.e., δ_a is small and positive).

27. **p92**: The statement that shell crossing occurs at D = 101 — is this an analytical result or a numerical observation? If analytical, the derivation (or at least a reference to where it's derived) would strengthen the claim.

### Section 4.2 — 3D FPA

28. **p96**: Section 4.2 states "the generalization to three dimensions is a simple extension to the existing equations" — this is true for the FPA but the text should note more carefully that this simplicity is specific to the FPA (because all momentum operators commute). In the full S-P system, the 3D extension is more complex, as is previewed in section 4.2.1.1.

29. **p98-99**: Equation 4.21: The gravitational potential is computed as φ_{g,init} = ∫∫ 4πG δ_init dV — this double integral notation is unusual. If this is meant to represent solving the Poisson equation ∇²φ = 4πGδ (which gives φ_k ~ δ_k/k² in Fourier space), the integral form should be the Green's function solution: φ(x) = -G ∫ δ(x')/|x-x'| dV'. The double integral as written is unclear.

30. **p99**: The 3D computational algorithm bullet list is a nice addition. However, the step "Evolve wavefunction (ψ). Jump to any D" is the key distinguishing feature of the FPA (no time-stepping required) and could be emphasised more.

### Section 4.2.3 — Cosmological Test

31. **p100**: The cosmological parameters used (h = 0.71, Ω_cdm,0 = 0.27, Ω_Λ,0 = 0.73, Ω_b,0 = 0, σ_8,0 = 0.81) set Ω_b = 0, meaning baryons are ignored entirely. This is a reasonable simplification for a CDM-only comparison but should be noted as a limitation.

32. **p101**: The correlation coefficients r = 0.9795, 0.8336, 0.7871, 0.8109 — the non-monotonic behaviour (decreasing then increasing) is noted but not explained. The text says "it is not clear why the correlation steadily decreases but then improves for the last comparison." This honesty is good, but some speculation would be welcome. Could it be related to the Hydra code's adaptive timestep resolution, or to shell-crossing effects that the FPA handles differently?

33. **p101**: "The initial conditions are, of course, the same, hence the correlation parameter is 1.0" — but the initial correlation is not listed in the four values given (0.9795, 0.8336, 0.7871, 0.8109). These four correspond to t = 162, 531, 739, 936. The text could be clearer that these are all post-initial comparisons.

### Section 4.2.4 — Consistency Checks

34. **p104-105**: The velocity tables show σ_{v_x} = σ_{v_y} = σ_{v_z} = 43.5139 to 6 significant figures — this remarkable isotropy is reassuring. Worth highlighting this more explicitly as a strong consistency check.

### Section 4.2.5 — Vorticity

35. **p106**: "That is when v = ∇φ no longer makes sense. This occurs most obviously when ψ = (0,0), ψ ∈ ℂ" — writing ψ = (0,0) for ψ = 0 + 0i is unusual. Simply "ψ = 0" is clearer.

36. **p107**: Figure 4.9 caption: "The velocity components have a theoretical gaussian overplotted upon them to show how far they deviate from theory" — "Gaussian" should be capitalised (it's derived from a proper name). Also, the caption says the components deviate from theory, but the text on p103 says they are close to Gaussian. The caption should say "to show how **closely** they **match** theory" or similar.

### Section 4.3 — Conclusion

37. **p107**: "The FPA is a fast and efficient method for probing the quasi-linear regime of density perturbations, it proves to be a good match for the Zel'dovich approximation" — comma splice. Should be a semicolon or two sentences.

38. **p108**: "The benefit of the FPA is that it runs much faster than all known N-body codes" — this is a strong claim. It's true for the quasi-linear regime, but should be qualified: "for the quasi-linear regime, the FPA runs much faster than N-body codes." An N-body code that only evolves to the same quasi-linear epoch would also be fast.

39. **p108**: Figure 4.10 shows the velocity comparison with a "turn-over at each end" described as "unexpected" and "an artefact from the construction of the wavefunction." This is an important observation that deserves more investigation. Is this a boundary effect? A consequence of the phase wrapping? The dismissal as simply an "artefact" without further analysis is unsatisfying.

40. **p109**: Short's perturbative approach is mentioned (adding a perturbation term to the free particle Hamiltonian with V still set to zero). The distinction between this perturbative approach and the full non-perturbative approach of Chapter 5 is clearly drawn, which is good.

---

## Equations Check

- **Eq 4.1**: iν ∂ψ/∂D = (-ν²/2 ∇² + V)ψ — correct (Schrödinger equation with D as time variable)
- **Eq 4.2**: ∇²Φ_g = 4πGρ_{b,c}(|ψ|² - 1) — correct (Poisson equation with overdensity)
- **Eq 4.3**: V = Φ_g - 3Ω_c/(2f²D) φ_v = 0 — correct (FPA condition)
- **Eq 4.4**: ψ, δ, v definitions — correct in form, but sign convention for v needs cross-check with Ch.3
- **Eq 4.5**: Fourier evolution of free particle — correct in principle; leading minus sign needs checking
- **Eq 4.6**: ∂²δ/∂t² + 2H ∂δ/∂t - 4πGρ̄δ = 0 — correct (linear growth equation)
- **Eq 4.7**: D_+ ∝ H ∫ da/(aH)³ — correct
- **Eq 4.8**: D_- ∝ H — correct
- **Eq 4.9**: D_+ with incomplete Beta function — plausible, standard result for ΛCDM
- **Eq 4.10**: α = Ω_{Λ,0}a³/(Ω_{cdm,0} + Ω_{Λ,0}a³) — correct
- **Eq 4.11**: 1D initial conditions — correct
- **Eq 4.12-4.13**: v_init = ∇φ_{v,init}; v_{x,init} = (p/2π)δ_a sin(2πx/p) — correct
- **Eq 4.14**: X_rms definition — correct
- **Eq 4.15**: Mean definition — correct
- **Eq 4.16**: Correlation coefficient — correct (Pearson correlation)
- **Eq 4.17**: V_difference RMS — correct but the (V_diff)_rms formula has 1/n^{1/2} rather than the standard √(1/n Σ...). As written: (V_diff)_rms = (1/n^{1/2}) Σ(V_diff - <V_diff>)². This is missing the square root over the sum. Should be √(1/n Σ(V_diff - <V_diff>)²) to be consistent with Eq 4.14.
- **Eq 4.18-4.20**: 3D initial conditions — correct (natural extension of 1D)
- **Eq 4.21**: Initial wavefunction construction — correct in form, but the double integral for φ_g is unclear (see note #29)

---

## Structural Comments

41. The chapter has a clear and logical progression: theory → 1D verification → 3D toy model → 3D cosmological test → consistency checks → vorticity → evaluation → bridge to Chapter 5. This is effective.

42. The 1D section (4.1) serves primarily as verification of the implementation against Short's published results. This is methodologically sound but could be more concise — the reader may feel that extensive 1D results for a known method are less interesting than the new 3D results.

43. The consistency checks section (4.2.4) is a valuable addition. The Gaussian density distribution, Gaussian velocity components, and Maxwellian speed distribution all provide confidence in the initial conditions. This section could be placed before the cosmological comparison (4.2.3) since it validates the inputs.

44. The chapter ending (4.3.1) effectively bridges to Chapter 5 by clearly articulating the FPA's limitations (shell-crossing breakdown, inability to probe non-linear regime) and motivating the need for the full Schrödinger-Poisson solver.

45. Pronoun consistency: the chapter switches between "I" (p92, p95, p100) and "we" (p97, p101, p105) without clear reason. Should be consistent throughout.
