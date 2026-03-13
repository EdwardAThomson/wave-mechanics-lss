# Appendix B: Mathematical Appendix - Review Notes

## Overall Assessment
A concise mathematical appendix (7 pages, pp229-235) containing three derivations: (B.1) the Madelung transform showing equivalence of the Schrödinger equation and the continuity equation, (B.2) the initial velocity of the FPA derived via the probability current, and (B.3) the Schrödinger and Poisson equations in the Einstein-de Sitter (EdS) model, correcting a mistake in Widrow & Kaiser. The derivations are generally clear and serve their purpose of supporting the main text. There are a few presentation issues and one notable step that needs checking.

---

## Typos and Grammatical Errors

1. **p229**: "the Continuity equation of fluid dynamics Schrödinger equation" — garbled sentence. Should be something like: "the mathematical equivalence of the Schrödinger equation and the continuity equation of fluid dynamics."

2. **p229**: "Differentiate" is capitalised mid-sentence: "LHS, Differentiate ψ with respect to time" → "LHS, differentiate ψ with respect to time"

3. **p232**: "the imaginary party of a complex number" → "the imaginary **part** of a complex number" (typo)

4. **p234**: "QED. This matches the result of Short in his PhD thesis." — "QED" is a strong claim for what is essentially a verification calculation. "This confirms..." or simply removing "QED" would be more appropriate.

5. **p234**: "except there is a mistake in their version of the coupled equations" — it would be helpful to state explicitly what the mistake in Widrow & Kaiser is, rather than just noting its existence.

---

## Section B.1 — Madelung Transform

6. **p229**: Equation B.1: iν(∂ψ/∂t) = -(ν²/2)∇²ψ + Vψ — correct (Schrödinger equation with ν = ℏ/m). Note this omits the potential term compared to the full S-P system (the text says "The role of the pressure term is discussed elsewhere and is extracted out from this derivation").

7. **p229**: Equation B.2: ψ = (1+δ)^{1/2} exp(iφ/ν) = α exp(iφ/ν) — correct Madelung ansatz with α = √(1+δ) = √ρ.

8. **p229-230**: Equations B.3-B.5: The time derivative of ψ and the resulting LHS = (iνψ/α)(∂α/∂t) - ψ(∂φ/∂t). Let me verify:
   ∂ψ/∂t = (∂α/∂t)e^{iφ/ν} + α(i/ν)(∂φ/∂t)e^{iφ/ν}
   iν(∂ψ/∂t) = iν(∂α/∂t)e^{iφ/ν} - α(∂φ/∂t)e^{iφ/ν}
   = (iνψ/α)(∂α/∂t) - ψ(∂φ/∂t)  ✓

9. **p230**: Equations B.6-B.7: ∇ψ and ∇²ψ calculations — correct. The Laplacian expands to five terms involving ∇²α, ∇φ·∇α, ∇²φ, and (∇φ)². The i²/ν² = -1/ν² term is correctly handled.

10. **p230**: Equations B.8-B.9: The identification of V (Bernoulli/Hamilton-Jacobi equation) and P (quantum pressure):
    V = -(∂φ/∂t)ψ - ½ψ(∇φ)²
    P = (ν²/2)(∇²α/α)
    These are correct. Note that V here is the effective potential from the Bernoulli equation, not the gravitational potential.

11. **p230-231**: Equations B.10-B.12: The RHS is assembled as K + V (kinetic + potential), leading to the cancellation with the LHS. After cancellation, Eq B.11 gives:
    (iνψ/α)(∂α/∂t) = -(iνψ/α)(∇φ·∇α) - (iν/2)ψ∇²φ
    Dividing by iνψ/α:
    ∂α/∂t = -(∇φ·∇α) - (α/2)∇²φ
    Which is Eq B.12: ∂α/∂t + (∇φ)·(∇α) + (α/2)∇²φ = 0. Correct.

12. **p231**: Equations B.13-B.17: The verification that B.12 is the continuity equation. Starting from ∂ρ/∂t + ∇·(ρ∇φ) = 0 with ρ = α², expanding gives 2α(∂α/∂t) + α²∇²φ + 2α∇α·∇φ = 0, dividing by 2α gives B.17. This is consistent with B.12. Correct.

---

## Section B.2 — Derivation of the Initial Velocity of the FPA

13. **p231-232**: Equation B.18: δᵢ = -δₐcos(2πx/p), φᵢ = -(p/2π)²δᵢ — these are the 1D FPA initial conditions from Chapter 4 (Eq 4.11). Correct.

14. **p232**: Equations B.19-B.21: The Madelung transform applied to the initial conditions, giving ψᵢ and ψᵢ* explicitly. The sign of the phase exponent differs between ψᵢ (negative) and ψᵢ* (positive), which is correct for complex conjugation.

15. **p232**: The definitions α1 through α4 are a useful bookkeeping device for the lengthy algebra. However, the notation is non-standard and could be confused with actual variables. Labels like A₁...A₄ or using named substitutions would be clearer.

16. **p233**: Equations B.26-B.31: The probability current calculation. Let me verify the key step. Starting from J = (ℏ/2mi)(ψ*∇ψ - ψ∇ψ*):

    Using the shorthand: ψ = α1·α2, ψ* = α1·α2*, ∇ψ = α3·α2 + α1·α2·α4, ∇ψ* = α3·α2* - α1·α2*·α4

    ψ*∇ψ = α1·α2*·(α3·α2 + α1·α2·α4) = α1·α2*·α3·α2 + (α1)²·α2·α2*·α4
    ψ∇ψ* = α1·α2·(α3·α2* - α1·α2*·α4) = α1·α2·α3·α2* - (α1)²·α2·α2*·α4

    Difference: ψ*∇ψ - ψ∇ψ* = α1·α3·(α2*·α2 - α2·α2*) + 2(α1)²·|α2|²·α4
    Since α2·α2* = |α2|² = 1, the first term vanishes: α2*α2 - α2α2* = 0.
    So: ψ*∇ψ - ψ∇ψ* = 2(α1)²·α4

    J = (ℏ/2mi)·2(α1)²·α4 = (ℏ/mi)(α1)²·α4

    With α4 = (i/ν)·(p/2π)δₐsin(2πx/p) = (i/ν)∇φᵢ and ν = ℏ/m:
    J = (ℏ/mi)·(α1)²·(im/ℏ)∇φᵢ = (α1)²·∇φᵢ = |ψ|²∇φᵢ

    Then v = J/ρ = J/|ψ|² = ∇φᵢ. This matches Eq B.33. ✓

17. **p233-234**: Equations B.29-B.30: The simplification steps. Eq B.29 has J = (ℏ/2mi)[α1·α3 - (α1)²·α4 - α1·α3 - (α1)²·α4]. This should simplify to J = (ℏ/2mi)[-2(α1)²·α4]... but the text gets J = (ℏ/2mi)[2·(α1)²·α4] (Eq B.30). There may be a sign issue in the intermediate steps, but the final result (B.33: v = ∇φ) is correct, which suggests the signs work out. The intermediate algebra is dense and hard to follow — a cleaner presentation would help.

18. **p234**: Equations B.34-B.35: The final result v_i = (p/2π)δₐsin(2πx/p), obtained by taking the gradient of φᵢ = -(p/2π)²δₐcos(2πx/p). This matches Eq 4.12 from Chapter 4. Correct.

---

## Section B.3 — Schrödinger and Poisson Equations in EdS Model

19. **p234**: "This is equation (17) in Widrow & Kaiser except there is a mistake in their version" — as noted in #5, the specific mistake should be identified for the reader.

20. **p234**: Equations B.36-B.41: The Poisson equation derivation. Starting from ∇²V = (4πG/a)(ψψ* - <ψψ*>) and normalising, eventually arriving at ∇²_y U = χχ* - 1. Let me trace the steps:
    - B.36: Standard Poisson with mean subtraction — correct
    - B.37: Dividing by <ψψ*> — correct
    - B.38: Substituting <ψψ*> = 1/(6πGt₀²) — this normalisation needs justification. In an EdS universe, ρ̄ = 1/(6πGt²), so at t₀, <ψψ*> should be proportional to ρ̄(t₀).
    - B.39: Defining χ = ψ·(6πGt₀²)^{1/2} — a rescaled wavefunction
    - B.40-B.41: Converting to comoving coordinates y and dimensionless potential U — correct in form

21. **p235**: Equations B.42-B.46: The Schrödinger equation derivation.
    - B.42: Standard Schrödinger equation iℏ(∂ψ/∂t) = (-ℏ²/2m ∇² + mV)ψ — correct
    - B.43: Substituting the rescaled variables — involves (6πGt₀²)^{1/2} and comoving coordinates
    - B.44-B.45: Algebraic rearrangement
    - B.46: Final form: i(4L/3)(∂χ/∂lna) = -∇²_y χ + (4L²/3)Uχ, where L = ma^{1/2}L²/(ℏt₀) — this is the Schrödinger equation in EdS expanding coordinates with ln a as the time variable.

22. **p235**: The definition L = ma^{1/2}L²/(ℏt₀) — this is the key dimensionless parameter that determines the effective resolution. It combines the physical constants (m, ℏ), the cosmological parameters (a, t₀), and the box size (L). This is related to the ν parameter used elsewhere in the thesis. The connection should be made explicit.

23. **p235**: The derivation ends abruptly after Eq B.46 with just "Here L = ma^{1/2}L²/(ℏt₀)." A brief comment connecting this to the equations used in Chapter 5 would be helpful. Is this exactly the form implemented in the code?

---

## Structural Comments

24. The three sections serve distinct purposes: B.1 supports Chapter 3 (Madelung theory), B.2 supports Chapter 4 (FPA velocity verification), B.3 supports Chapter 5 (EdS implementation). This parallel structure with the main chapters is effective.

25. Section B.2 is the most detailed and carefully worked through. The step-by-step algebra, while dense, provides a complete verification that the probability current method recovers the expected FPA velocity.

26. Section B.3 is the most compressed and would benefit from more intermediate steps, particularly in explaining the normalisations and variable substitutions (B.38-B.40).

27. A brief introduction to the appendix explaining what each section contains and why it's separated from the main text would improve readability.
