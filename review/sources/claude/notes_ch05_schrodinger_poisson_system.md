# Chapter 5: Solving the full Schrödinger-Poisson system - Review Notes

## Overall Assessment
This is the central chapter of the thesis (67 pages) and presents the main original contribution: a 3D wave-mechanics code that solves the full coupled Schrödinger-Poisson system with self-consistent gravity, expanding coordinates, periodic boundaries, and mass conservation. The chapter covers the numerical method (Cayley transform, Goldberg scheme, splitting operators), periodic boundary implementation (novel contribution), expansion (Einstein-de Sitter and general flat FLRW), the Poisson solver, initial conditions, and extensive results including 1D tests, 2D gravitational interaction, 3D tophat collapse, and full 3D cosmological simulations compared with GADGET-2.

The chapter is impressively thorough and honest about limitations — the velocity scaling discrepancy, the net bulk motion problem, and the interference effects are all openly discussed. However, the chapter is very long and could benefit from tighter organisation. There are several typos, some mathematical presentation issues, and the results section would benefit from clearer separation of what works well from what remains problematic.

---

## Typos and Grammatical Errors

1. **p112**: "This requirements are outlined" → "**These** requirements are outlined"

2. **p112**: "expression of a exponential" → "expression of **an** exponential"

3. **p113**: "colloquially called the Cayley method" — "colloquially" is slightly odd here; "commonly" would be better.

4. **p113**: "I decided to use the same method as Widrow & Kaiser, which is a finite difference scheme based upon using Cayley's decomposition of an exponential (colloquially called the Cayley method)." — This sentence attributes the Cayley method to Widrow & Kaiser but then the text on p117 says it's from Goldberg et al (1967). The attribution is unclear.

5. **p115**: "the space for skew-Hermitian matrices forms the Lie algebra u(n) of the Lie group U(n)" — correct mathematical statement, well put.

6. **p116**: "the the operator" → double "the" (if present in text on p116 about operators acting)

7. **p118**: "The Crank-Nicolson method is a 2nd order finite difference scheme that is implicit and unconditionally stable. The scheme employed here is implicit but conditionally stable" — this is a crucial distinction that could easily be missed. The reader might assume unconditional stability from the Crank-Nicolson reference. Consider highlighting this difference more prominently.

8. **p118**: "For this derivation ℏ = 1 and m = 1/2" — setting m = 1/2 is an unusual convention. While it simplifies the kinetic energy term to -∇², it may confuse readers who expect m = 1. Worth a brief justification.

9. **p119**: "where λ = 2ε²/δ" — the reuse of δ for the timestep (having already used it for density contrast and for Kronecker/Dirac delta throughout the thesis) is potentially confusing. Consider using Δt or τ for the timestep.

10. **p121**: "My attempts to reproduce Watanabe's 'adhesive operators'... has been unsuccessful" → "**have** been unsuccessful" (subject-verb agreement with plural "attempts")

11. **p124**: "Widow & Kaiser" → "**Widrow** & Kaiser" (typo — missing 'r')

12. **p127**: "Widow & Kaiser" appears again — same typo as above.

13. **p128**: "the the operator" — check for double "the" in the discussion of higher dimensions.

14. **p129**: "so the the operator closest to the wavefunction" → double "the"

15. **p130**: "straight forward" → "**straightforward**" (one word)

16. **p132**: "Through out" → "**Throughout**" (one word)

17. **p139**: "minorly apparent" → awkward phrasing; "slightly apparent" or "barely visible" would be better.

18. **p141**: "over lap" → "**overlap**" (one word)

19. **p145**: "As the two peaks move towards each other they are observed to squeeze themselves under gravity" — "squeeze themselves" is somewhat colloquial. Consider: "are observed to contract under their mutual gravitational attraction."

20. **p149**: "(Coles, (Coles 1997))" — double citation formatting error.

21. **p155**: "the gaussian shape of the density field skews over time with the peak leaning towards the higher mass end" — "Gaussian" should be capitalised (derived from proper name Gauss). This applies throughout the chapter.

22. **p158**: "there is a a an accepted standard" → "there is **an** accepted standard" (triple article error)

23. **p164**: "there is a a" — check for similar stuttering.

24. **p170**: "the table indicates that there is a net motion in the x direction, this is also true of the other components" — comma splice.

---

## Clarity and Content Issues

### Section 5.1 — Specific requirements

25. **p111-112**: The five requirements for a cosmological wave-mechanics code (3D, self-consistent gravity, expansion, periodic boundaries, mass conservation) are clearly stated and form a good framework for the chapter. However, requirement 5 ("mass conserving") could be strengthened to "conserves mass and momentum" since both are demonstrated later (Fig 5.22).

26. **p112**: Equation 5.2: ∇²V = 4πGψψ* — this is the Poisson equation with density = ψψ*. But compared with the FPA's Eq 4.2 which has ∇²Φ_g = 4πGρ_{b,c}(|ψ|² - 1), the subtraction of the mean density (-1 term) is missing here. The mean subtraction is important for periodic boundary conditions (the Poisson equation has no solution for a non-zero mean source in a periodic box). This should be noted.

### Section 5.1 — Matrix exponential and Cayley transform

27. **p114**: Equations 5.5-5.6: The Cayley method of writing ψ(t+dt) = (e^{i½Hdt})^{-1} e^{-i½Hdt} ψ(t) is elegant and ensures unitarity. The explanation is clear.

28. **p114-115**: The Cayley transform (Eq 5.7): A = (I+M)(I-M)^{-1} — the connection to the numerical scheme is well made. However, the statement "the order of the terms on the right hand side does not matter, as the 'top' and 'bottom' brackets commute with each other" on p115 needs qualification. (I+M) and (I-M)^{-1} commute because they are both functions of M. This is not obvious and could use a sentence of explanation.

29. **p115**: Equation 5.9: [x,p] = xp - px = iℏ — correct. However, this is the canonical commutation relation for quantum mechanics. In the numerical scheme, the non-commutativity comes from the discrete representation of operators in different bases (position vs momentum), not directly from the canonical commutation relation. The text could make this distinction clearer.

### Section 5.1 — Splitting operators

30. **p116**: Equation 5.12 (Goldberg scheme): ψ(t+dt) = [e^{-i(K_x+V_x)dt}[e^{-i(K_y+V_y)dt}[e^{-i(K_z+V_z)dt}ψ(t)]]] — the splitting of V into V_x, V_y, V_z (one-third each) is ad hoc. The potential V is a function of position and doesn't naturally split into directional components. This is acknowledged ("does not solve the problem of commutativity") but the motivation for splitting V equally into thirds could be explained better.

31. **p116-117**: The transition from the Goldberg scheme to Suzuki splitting operators is the key methodological insight. The text states the Goldberg scheme "breaks the unitarity" — it would strengthen the argument to quantify by how much (e.g., what level of mass non-conservation results).

32. **p116**: Equation 5.13: e^{a(P+X)} = [S_m(a/n)]^n + O(a^{m+1/n^m}) — this general formula for splitting operators is stated without much context. The relationship between the order m, the number of sub-steps n, and the error should be made clearer.

33. **p117**: Equation 5.15: f_2(P,X) = S(a) + O(a³) = e^{(a/2)P} e^{aX} e^{(a/2)P} — this is the standard Strang (second-order) splitting. It would help to explicitly state this is the one used in the final code.

### Section 5.2 — Numerical method

34. **p118**: Equation 5.16: ψ(x,t+dt) = (1 - ½iHdt)/(1 + ½iHdt) ψ(x,t) — this is the Cayley form of the time evolution. The relationship to Crank-Nicolson is noted correctly but the distinction between unconditional (CN) and conditional (this scheme) stability is important and somewhat buried.

35. **p119**: Equations 5.20-5.27: The recursion relations for the auxiliary functions e, f, and Ω are the core of the Goldberg method. The derivation is clear but dense. A worked example or diagram showing the sweep direction would help readers implement this.

36. **p119**: "the advanced potential cannot be known as the system has not evolved" — this is the key approximation: using V^n instead of V^{n+1}. The text correctly notes this requires small timesteps. However, it would be worth noting that this makes the scheme not truly implicit (it's semi-implicit or quasi-implicit).

### Section 5.2.2 — Periodic boundaries

37. **p121-123**: The periodic boundary implementation is described as a novel contribution. The core idea — re-iterating the recursion relations for e, f, Ω, and ψ using the last value to seed the first — is clever and clearly explained. Figure 5.1 provides good visual evidence.

38. **p121**: "My attempts to reproduce Watanabe's 'adhesive operators'... has been unsuccessful; however, I developed my own method" — the honesty is commendable, but it would be useful to state *why* the adhesive operator approach failed (was it a coding error? a misunderstanding of the method? a fundamental incompatibility?).

39. **p123**: "In practice we found that the recursion relations only have to be performed twice" — this is an important practical result. Is there a theoretical argument for why two iterations suffice? The convergence behaviour (exponential? linear?) would be useful to characterise.

### Section 5.2.3 — Expansion

40. **p124-126**: The analysis of Widrow & Kaiser's equations, discovering that their density is comoving rather than physical, is a genuine contribution. The detective work is well presented. However, the passage is quite long (3 pages) for what amounts to a factor-of-a³ correction. Consider condensing.

41. **p125**: Equations 5.40-5.45: The chain of reasoning showing <χχ*> = 1 (comoving) vs <ψψ*> = 1/(6πGt₀²) is correct and important. The presentation is clear.

42. **p126**: Equations 5.46-5.47: "here I provide them as they should have appeared" — identifying a typographical error in Widrow & Kaiser (coupled equations appearing as one) is a useful contribution. Cross-referencing Appendix B.3 for the derivation is good.

43. **p126**: "L = ma^{1/2}L²/(ℏt₀)" — this dimensionless parameter L (script L) is crucial but its physical meaning is not well explained. It's essentially the ratio of the box size to the de Broglie wavelength. A sentence on its physical significance would help.

44. **p127**: Section 5.2.3.2: The generalisation to flat models with non-zero cosmological constant is straightforward but important. Equation 5.54 is the key result — the generalised Schrödinger equation with the [Ω_{m₀} + (1-Ω_{m₀})a³]^{1/2} factor.

### Section 5.2.4 — Higher dimensions and splitting operators

45. **p128-129**: Equation 5.55: The naive Goldberg extension splits V/3 into each dimension. Equation 5.56 then uses Strang splitting e^{-iKdt/2} e^{-iVdt} e^{-iKdt/2}. Combined with [P_x, P_y] = 0, the kinetic energy splits exactly (Eq 5.57). This is the key technical result and is clearly presented.

46. **p129**: Equation 5.57: ψ(x,t+dt) = [(1-i(dt/2)K_x)/(1+i(dt/2)K_x)] × [(1-i(dt/2)K_y)/(1+i(dt/2)K_y)] × [(1-i(dt/2)K_z)/(1+i(dt/2)K_z)] ψ(x,t) — this shows that the kinetic energy can be split exactly because momentum operators in different directions commute. Good.

47. **p130**: Equation 5.58: e = 2 - iλ — the new form of e without the V term (since potential is handled separately by splitting) is a significant simplification. Worth highlighting that this makes the kinetic sweep independent of the potential.

48. **p130**: Equation 5.59: e^{-iVdt}ψ = (1 - iδV/2)(1 + iδV/2)^{-1}ψ — since V is diagonal in position space, this is exact (no approximation beyond the Cayley first-order expansion). The text says "straight forward" but this is actually a key advantage of the splitting approach.

### Section 5.3 — Poisson equation

49. **p130-131**: Equation 5.60: ∇²V = 4πGψψ* — again, the mean subtraction is missing. In a periodic box, one must solve ∇²V = 4πG(ψψ* - <ψψ*>). The k=0 mode of V is conventionally set to zero. This should be stated.

50. **p131**: Equation 5.62: V_k = 4πGψ_kψ*_k/k² — this is incorrect. The Fourier transform of ψψ* is the convolution (ψ_k * ψ*_k), not the product ψ_kψ*_k. What should be written is: ρ_k = FT[ψψ*], then V_k = 4πGρ_k/k². The equation as written conflates the density in Fourier space with the product of the wavefunction's Fourier transform and its conjugate.

51. **p131**: Equation 5.63: V_k = 4πGψ_kψ*_k/(2κ-3) with the cosine kernel — same issue as above regarding ψ_kψ*_k vs ρ_k. Also, this is the same discrete Poisson kernel as Eq 2.4 from Chapter 2. Good to cross-reference.

### Section 5.5 — Initial conditions

52. **p132**: Equation 5.64: ψ_i = (2πw²)^{-1/2} exp(-(x_i-x₀)²/(4w²) + ip₀x_i) — correct Gaussian wave-packet. The normalisation factor should be (2πw²)^{-1/4} for a properly normalised 1D Gaussian wavefunction (since |ψ|² must integrate to 1 and the Gaussian integral of e^{-x²/(2w²)} gives (2πw²)^{1/2}). As written with (2πw²)^{-1/2}, the density |ψ|² would have normalisation (2πw²)^{-1}, which integrates to (2πw²)^{-1} × (2πw²)^{1/2} = (2πw²)^{-1/2} ≠ 1.

53. **p132**: Equation 5.65: The two-Gaussian initial condition drops the normalisation factor. The text justifies this by saying total normalisation just needs to be constant over time, which is valid.

54. **p133**: The cosmological initial conditions use GADGET-2 (Springel 2005) as the comparison code rather than Hydra (used in Chapter 4). This switch is noted but the reason could be stated more explicitly. Is it simply that GADGET-2 is more widely used?

### Section 5.6 — Results

55. **p134**: The structure of the results section (periodic boundaries → 1D tests → 2D gravitational interaction → 3D tophat → cosmological simulation) is logical, building from simple to complex. Good.

56. **p135-137**: Section 5.6.1 on periodic boundaries: The explanation of the double recursion and Figure 5.1 showing convergence after two iterations is convincing. The 3D wave-packet crossing boundaries (Figure 5.2) is a strong visual demonstration.

57. **p139-140**: The discussion of interference/feedback effects from periodic boundaries is physically insightful. The key point — that for cosmological values of ν (~10⁻⁷), dispersion is slow enough that feedback is negligible during the simulation time — is reassuring but should be stated more quantitatively. What is the dispersion timescale compared to the simulation time?

58. **p141**: The observation that gravity suppresses interference effects is interesting and physically reasonable (gravity maintains coherence by keeping the wavefunction localised).

59. **p143-148**: Section 5.6.2 (2D two-body gravitational interaction): The sequence of figures (5.7, 5.8, 5.9) showing attraction, pass-through, and turn-around is compelling and demonstrates multi-streaming. This is one of the strongest results in the chapter. The double-slit-like interference pattern during the pass-through (timesteps 35-50) is a genuinely interesting result.

60. **p149-150**: Section 5.6.3 (tophat collapse): The dependence on ν is clearly demonstrated. The observation that larger ν causes faster collapse is consistent with the de Broglie wavelength interpretation (larger ν → larger quantum pressure → but also faster dispersion?). Actually, re-reading, the text says collapse happens faster with larger ν, but physically larger ν means more quantum pressure which should *resist* collapse. This seems contradictory and may need checking.

61. **p150**: "In tests where ν is much larger the over-density expands so fast that the wavefunction crosses over the boundaries" — this confirms that larger ν leads to faster dispersion, not faster collapse. The earlier statement about collapse happening "faster" with ν = 10⁻⁷ vs 10⁻⁸ may be because both values are small enough that quantum pressure is negligible, and the dynamics are dominated by the expansion/time-stepping. This subtlety should be clarified.

62. **p153**: "The simulations at resolutions of 32³ and 64³ prove to be unsatisfactory as the magnitude of δ appears to be too high" — this is an important limitation. The initial |δ| ~ 0.3 growing to δ ~ 6 in 100 timesteps is indeed too fast. The suggestion that it's related to the TSC smoothing of particle positions is plausible but not proven.

63. **p155**: "The highest peak of density (and over-density) never goes as high as it does in the Gadget outputs" — combined with "the lowest trough of density is far lower", this suggests the wave-mechanics code produces a wider dynamic range but with lower peak densities. This is a significant and clearly stated finding.

64. **p155**: "ρ_min ~ 10⁶ for wave-mechanics while Gadget only goes as low as ρ_min ~ 10⁹·⁵" — these are absolute densities (in code units). The factor of ~10³·⁵ difference in minimum density is large and unexplained.

65. **p158**: "there is a a an accepted standard" — beyond the typo, the discussion of Gaussian smoothing with σ = 8 Mpc is relevant. The connection to the σ₈ normalisation convention could be made more explicit.

### Section 5.6.4.2 — Conserved quantities

66. **p164-167**: The conservation of mass (to 7 significant figures) and momentum over 3500 timesteps (Figure 5.22) is an important result that validates the numerical scheme. This should perhaps be presented earlier or more prominently, as it addresses one of the five core requirements.

### Section 5.7 — Velocity & Vorticity

67. **p168**: Equation 5.70: v = (ℏ/m)ℑ(∇ψ/ψ) — correct, consistent with Eq 3.52.

68. **p168**: "The velocity data from the Gadget output files requires a simple scaling in order to calculate the real physical velocities: v_real = √a · v_code" — is this the standard GADGET velocity scaling? GADGET uses "peculiar velocities" = a^{1/2} × dx/dt. This is worth double-checking as an incorrect scaling could explain the velocity discrepancy.

69. **p169**: The velocity discrepancy is the most troubling open issue. The initial velocities from wave-mechanics don't match GADGET's, despite identical initial density fields. Two possible explanations are offered: (1) gravitational potential to velocity potential scaling, (2) velocity field calculation error. A third possibility not mentioned: the wave-mechanics velocity is calculated from the probability current of the constructed wavefunction, but the wavefunction was constructed using the gravitational potential (not velocity potential) as the phase. If the two potentials are not exactly proportional at the initial time (especially after TSC smoothing), the derived velocity will differ from GADGET's.

70. **p170**: The velocity table shows mean[V_x] growing from -3.7 to -165 km/s over the simulation — this is a significant net bulk motion and clearly unphysical. The text honestly acknowledges this as "worrying" and an unsolved problem. This is a substantial issue that undermines confidence in the velocity outputs.

71. **p170**: "as sen in section 4.2.5" → "as **seen** in section 4.2.5" (typo)

72. **p174**: "min[ψ_final] ~ 10" — this means ψ never gets close to zero, so no vorticity singularities form. This is surprising given that the FPA showed near-zero values at D = 30. The full S-P system with gravity may be keeping ψ away from zero, but this deserves more discussion.

### Section 5.8 — Summary

73. **p174-177**: The summary is honest and well-balanced, clearly stating both achievements (novel periodic boundary implementation, mass/momentum conservation, generalised Schrödinger equation) and limitations (velocity discrepancy, interference effects, messy density fields). The acknowledgement that "we cannot conclude whether our results are completely correct" is admirably candid for a thesis.

---

## Equations Check

- **Eq 5.1**: iℏ∂ψ/∂t = (-ℏ²/2m ∇² + mV)ψ — correct
- **Eq 5.2**: ∇²V = 4πGψψ* — correct in form but missing mean subtraction for periodic box
- **Eq 5.3**: ψ(x,t+dt) = e^{-iHdt}ψ = e^{-i(K+V)dt}ψ — correct
- **Eq 5.4**: e^M = 1 + M + M²/2! + ... — correct (matrix exponential)
- **Eq 5.5-5.6**: Cayley evolution — correct
- **Eq 5.7**: A = (I+M)(I-M)^{-1} — correct (Cayley transform)
- **Eq 5.8**: (e^{-M})^{-1}e^M = (1+M)(1-M)^{-1} — correct to first order
- **Eq 5.9**: [x,p] = iℏ — correct
- **Eq 5.10**: e^K e^V ≠ e^V e^K — correct (non-commuting operators)
- **Eq 5.12**: Goldberg 3D scheme — correct in form
- **Eq 5.13**: General splitting operator formula — correct
- **Eq 5.14**: First order splitting — correct
- **Eq 5.15**: Second order (Strang) splitting — correct
- **Eq 5.16**: Cayley time evolution — correct
- **Eq 5.17-5.18**: Discretised evolution — correct
- **Eq 5.19**: Centred difference — correct
- **Eq 5.20-5.29**: Recursion relations — appear correct (standard Goldberg scheme)
- **Eq 5.30-5.34**: Periodic boundary double recursion — appear correct
- **Eq 5.35-5.36**: S-P system in expanding coordinates — correct
- **Eq 5.37-5.39**: Dimensionless transforms — correct
- **Eq 5.40-5.45**: Comoving density argument — correct
- **Eq 5.46-5.47**: Corrected Widrow & Kaiser equations — plausible
- **Eq 5.48**: Evolution in Einstein-de Sitter — correct
- **Eq 5.49-5.54**: General flat FLRW equations — appear correct
- **Eq 5.55-5.57**: Higher-dimensional splitting — correct
- **Eq 5.58**: Simplified e without V — correct
- **Eq 5.59**: Potential update via Cayley — correct
- **Eq 5.60**: ∇²V = 4πGψψ* — missing mean subtraction
- **Eq 5.61**: V = ∫4πGψψ*dV — unclear notation (see note #50)
- **Eq 5.62**: V_k = 4πGψ_kψ*_k/k² — **incorrect**: should be V_k = 4πGρ_k/k² where ρ_k = FT[ψψ*]
- **Eq 5.63**: Discrete Poisson with cosine kernel — same ψ_kψ*_k issue
- **Eq 5.64**: Gaussian wave-packet — normalisation should be (2πw²)^{-1/4}
- **Eq 5.65**: Two-Gaussian initial condition — correct (unnormalised)
- **Eq 5.66-5.67**: Gaussian smoothing — correct
- **Eq 5.68**: M_total = ∫ψψ*dx³ — correct
- **Eq 5.69**: P_total = ∫|ψ*∇ψ|²dx³ — this doesn't look right. Momentum should be P = ∫ψ*(-iℏ∇)ψ dx³ or equivalently P = ℏ∫ℑ(ψ*∇ψ)dx³. The expression |ψ*∇ψ|² has the wrong dimensions and form for momentum. **Needs checking.**
- **Eq 5.70**: v = (ℏ/m)ℑ(∇ψ/ψ) — correct

---

## Structural Comments

74. At 67 pages, this is by far the longest chapter and arguably tries to cover too much ground. The mathematical method (sections 5.1-5.3), the expansion generalisation (5.2.3), and the results (5.6) could each be substantial chapters in their own right. Consider whether the expansion derivation (particularly the Widrow & Kaiser correction) could be placed in an appendix.

75. The results section (5.6) occupies approximately 35 of the 67 pages and contains many figures. While the thoroughness is admirable, some figures are similar enough that they could be consolidated. For example, Figures 5.4, 5.5, 5.6 (three 1D two-peak tests) could be presented more concisely.

76. The chapter's most novel contributions are: (1) the periodic boundary implementation via double recursion, (2) the correction to Widrow & Kaiser's equations, (3) the generalisation to non-EdS cosmologies, and (4) the 3D cosmological results. These deserve to be highlighted more prominently, perhaps with a clear statement at the start of each section saying "this is new."

77. The unresolved issues (velocity scaling, net bulk motion, interference effects, initial conditions sensitivity) are honestly presented but spread throughout the chapter. A dedicated "Known Issues and Limitations" subsection within 5.8 would help the reader understand the state of the art at a glance.

78. Pronoun inconsistency continues from Chapter 4: "I" and "we" are mixed throughout without clear convention.
