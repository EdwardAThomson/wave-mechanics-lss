# Chapter 2: Review of Numerical Simulations - Review Notes

## Overall Assessment
A solid, concise review chapter (16 pages). It covers the necessary ground: N-body methods (direct summation, particle-mesh), symplectic integrators, fluid/perturbation theory, and motivates the wave-mechanical approach. Considerably cleaner than Chapter 1 in terms of typos. The symplectic manifold section (2.2) is ambitious for a review chapter but relevant given the thesis's emphasis on energy conservation.

---

## Typos and Grammatical Errors

1. **p38**: "that under go classical gravitational attraction" → "that **undergo** classical gravitational attraction" (one word)

2. **p38**: "scales as N × N" — more conventionally written as N² or O(N²), which you do use later on p43. Minor inconsistency.

3. **p39**: "mass resolution is  10⁹" — appears to be missing a symbol, likely "∼ 10⁹" or "≈ 10⁹". The tilde/approximately may have been lost in typesetting.

4. **p40**: "(Where N is the number of bodies)" — "Where" should not be capitalised mid-sentence.

5. **p41**: "computational practise" — in British English "practise" is the verb and "practice" is the noun. Here it's used as a noun ("computational practice"), so should be "practice".

6. **p43**: "Leapfrog integrator. While this method only gives positions accurate to third order (total error is O(h²))" — this is potentially confusing. The Leapfrog/Störmer-Verlet method is second-order accurate (local error O(h³), global error O(h²)). Saying "positions accurate to third order" and then "(total error is O(h²))" is mixing local and global error descriptions. Suggest clarifying: "second-order accurate (global error O(h²))".

7. **p43**: "Direction summation methods" → "**Direct** summation methods" (typo)

8. **p43**: "colloquially called P³M ('P-cubed-M')" — earlier on p37 you write "P³M" and "P-P" consistently, which is good.

9. **p46**: "if ω = 0 then either dq or dp must be zero. For real quantities this is always the case" — this is confusingly worded. The point is that for the canonical symplectic form on ℝ²ⁿ, non-degeneracy is automatic because the matrix representation of ω is the standard symplectic matrix J, which has det(J) = 1 ≠ 0. As written, "for real quantities this is always the case" seems to say ω is always non-degenerate, which is the opposite of ω = 0.

10. **p46**: "divergenceless (dω = 0)" — the closure condition dω = 0 is not the same as being "divergenceless". The term "divergenceless" typically refers to a vector field; "closed" is the standard term for a differential form with vanishing exterior derivative. Suggest simply saying "closed (dω = 0)".

11. **p47**: "all symplectic manifolds (therefore, groups) are simply connected" — this is not correct. Symplectic manifolds are not in general simply connected. For example, the 2-torus T² admits a symplectic structure but is not simply connected. The symplectic *group* Sp(2n,ℝ) is connected but not simply connected (its fundamental group is ℤ). This claim needs correction or significant qualification.

12. **p47**: The chain of reasoning "From Noether's theorem we know that the Poincaré group, of which the Abelian (symplectic) group is a subgroup..." — the symplectic group is not Abelian (except in the trivial 2D case), and the relationship between the Poincaré group and the symplectic group is not standard. The Poincaré group is the isometry group of Minkowski spacetime; the symplectic group acts on phase space. They are related through the broader structure of Hamiltonian mechanics but one is not a subgroup of the other in any standard sense. This paragraph needs reworking.

13. **p49**: "a function only of space and time (3 + 1D)" — minor: the parenthetical should be "(3+1)D" with the parenthesis around "3+1" to clarify it's 3 spatial + 1 temporal dimension, which you write correctly.

14. **p50**: "there is difficulty following fine structure in phase space" — should be "there is **a** difficulty" or "it has difficulty".

---

## Clarity and Content Issues

15. **p37**: "N = 1024³" for the Millennium Simulation — you give N = 10²⁴³ on this page and later describe "billion body" simulations. Worth noting explicitly that 1024³ ≈ 10¹⁰ ≈ 10 billion, so readers can see the connection.

16. **p38-39**: The discussion of superparticles and mass resolution is clear but somewhat repetitive — the concept of a "superparticle" is introduced on p38 and then again on p39. Could be tightened.

17. **p40**: The discussion of Sundman's regularization for N=3 and Wang's generalization for N>3 is interesting historical context but the relevance to cosmological N-body simulations (where N ~ 10⁹⁻¹⁰) is not explicitly drawn. A sentence connecting this back would help: something like "while these analytical results are of theoretical interest, cosmological simulations with N >> 3 must rely entirely on numerical methods."

18. **p42**: Equation 2.2 — the softened force equation uses |ε| << 1, but the absolute value bars seem unnecessary if ε is a length parameter (always positive). Also, the condition |ε| << 1 is dimensionful — should this be ε << L for some characteristic length L, or ε << mean interparticle separation?

19. **p44**: Equation 2.3: ∇²Φ = 4πGa²ρ̄δ — this is the Poisson equation in comoving coordinates. Worth noting that this uses the over-density δ = (ρ - ρ̄)/ρ̄ and that the source is the *perturbation* to the density, not the total density. This avoids the issue of a constant background producing a uniform potential.

20. **p44**: Equation 2.4: Φ_k = 4πGρ_k / (2κ - 3) — I'd like to double-check this. The standard discrete Poisson solver using a 3D cosine kernel gives a denominator involving the sum of cosines minus a constant, but the factor of 3 in "(2κ - 3)" implies κ is the sum of three cosine terms, each contributing 1 when the mode indices are zero. This seems correct in form but the specific numerical factors depend on the discretisation scheme. Worth cross-referencing with Chapter 5 where you say this is detailed further.

21. **p46-47**: Section 2.2 on symplectic manifolds is quite mathematically dense for a review chapter aimed at astrophysicists. The key message — "symplectic integrators conserve energy, which is important for cosmological simulations" — could be delivered more concisely. The current version attempts to derive the full mathematical machinery (differential forms, Liouville's theorem, Noether's theorem, Poincaré group) in about two pages, which is very compressed and leads to some of the errors noted above (#11, #12).

22. **p48**: The linearised fluid equations (Eq 2.6) drop the v·∇v term from the Euler equation (Eq 1.15 in Ch.1) without comment. It would help to note this is the linearisation step — dropping the non-linear advection term.

23. **p49**: "solutions for δ ∝ tⁿ are n = 2/3 (growing mode) or −1 (decaying mode)" — for an Einstein-de Sitter universe with a ∝ t^{2/3}, the growing mode of the density perturbation goes as δ ∝ a ∝ t^{2/3} and the decaying mode as δ ∝ t^{-1} ∝ a^{-3/2}. This checks out.

24. **p50**: The summaries of the N-body, Phase Space, and Fluid methods are presented in an indented block format that reads almost like bullet points. This is fine but slightly inconsistent with the rest of the chapter's prose style. Consider whether these could be formatted as a proper comparison table.

---

## Equations Check

- **Eq 2.1**: F_ij = -G Σ m_i m_j (r_i - r_j)/|r_i - r_j|³ — correct (Newton's gravitational force)
- **Eq 2.2**: r̈_i = -G Σ m_j (r_i - r_j)/(|r_i - r_j|² + ε²)^{3/2} — correct (softened acceleration)
- **Eq 2.3**: ∇²Φ = 4πGa²ρ̄δ — correct (comoving Poisson equation)
- **Eq 2.4**: Φ_k = 4πGρ_k/(2κ-3), κ = cos(2πn/L)+cos(2πm/L)+cos(2πo/L) — plausible discrete form, needs cross-check with Ch.5
- **Eq 2.5**: ω = Σ dqⁱ ∧ dp_i — correct (canonical symplectic form)
- **Eq 2.6**: Linearised fluid equations (continuity + Euler + Poisson) — correct

---

## Structural Comments

25. The chapter flows well: N-body overview → direct summation → particle-mesh → symplectic integrators → fluid/perturbation theory → motivation for wave-mechanics. The logical progression makes the case for wave-mechanics effectively.

26. Section 2.2 (symplectic) is the weakest section due to compressed mathematical treatment leading to errors. Consider either expanding it properly or reducing it to the key practical message about energy conservation and citing references for the mathematical details.

27. The chapter ending works better than Chapter 1's — it explicitly sets up Chapter 3.
