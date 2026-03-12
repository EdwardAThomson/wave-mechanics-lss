# Chapter 6: Conclusion - Review Notes

## Overall Assessment
A short chapter (13 pages) that serves two functions: summarising the thesis contributions and results (pp178-182), and presenting future work directions (pp183-190). The summary section is honest about both achievements and shortcomings, which is commendable. The future work section (6.1) is thoughtful and well-structured, covering multiple fluids, additional physics, periodic boundaries via adhesive operators, mesh refinement, and parallelisation. However, there are some grammatical issues, a few places where claims could be better qualified, and the equations in the future work section have some presentation issues.

---

## Typos and Grammatical Errors

1. **p178**: "previous inconsistencies and misconceptions were explored and resolved to the best of our ability" — "to the best of our ability" is somewhat hedging for a thesis conclusion. Consider "were explored and, where possible, resolved" or simply "were explored and resolved."

2. **p179**: "It is not a term that is 'added' in by the Schrödinger equation" — "added in by" is colloquial. Better: "added **to** the Schrödinger equation" or "introduced by the Schrödinger equation."

3. **p179**: "Johnston, considered a different approach from myself and Short" — comma after "Johnston" is incorrect (it separates subject from verb). Also "from myself" is non-standard; should be "from **me** and Short" or "from that of Short and **myself**."

4. **p180**: "colloquially we call it the Cayley method" — "colloquially" is not quite right here; this is more of an informal shorthand than colloquial language. "We refer to it as the Cayley method" would be cleaner.

5. **p180**: "this process is trickier than it may first appear" — fine, but note the switch from "I" (previous paragraph) to implied "we" here. Pronoun consistency issue throughout.

6. **p180**: "(3+1 d)" — should be "(3+1)D" or "(3+1)-dimensional" for consistency with standard notation.

7. **p181**: "Through out this work" — if this appears as two words, it should be one: "**Throughout** this work." (Note: may be correct in the PDF; hard to tell from typesetting.)

8. **p181**: "suffers less from discretization problems" — should be "**suffer** less" (plural subject: "they describe... and suffer less").

9. **p182**: "although not necessarily inconsistent, results for density and velocity" — the comma placement makes this hard to parse. Suggest: "the unexpected — although not necessarily inconsistent — results for density and velocity."

10. **p182**: "this was a problem for generating our initial conditions" — the antecedent of "this" is ambiguous. The coarse resolution was the problem, but the sentence structure makes it unclear. Suggest: "this coarse resolution posed problems for generating initial conditions."

11. **p182**: "gaussian filter" and "gaussian" (lower case) — should be capitalised: "**Gaussian** filter" (derived from a proper name). Appears twice on this page.

12. **p183**: "I see my thesis and all previous work as being the first generation of wave-mechanics codes" — this is a bold but interesting claim. However, "first generation" typically implies a group of related but independent efforts; the text could clarify whether this means "pioneering generation" or "first wave of development."

13. **p184**: "they are less conventional and hence far more speculative" — comma splice. Should be a semicolon or "and are hence."

14. **p186**: "straight forward" — should be one word: "**straightforward**" (if it appears as two words; same issue as noted in Chapter 4).

15. **p186**: "I don't believe" and "I don't think" — contractions are informal for a thesis. Consider "I do not believe."

16. **p190**: "the next second processor" — unclear phrasing. Should be "the second processor on the next iteration" or similar.

---

## Clarity and Content Issues

### Summary Section (pp178-182)

17. **p178**: "The quantum nature of dark matter particles is mostly unknown but not thought to be significant. However this does not affect the main outcome of this work." — This is an important caveat but it's somewhat buried. The argument is: even if DM is quantum, the high occupation number means the single-particle Schrödinger equation (classical limit) is adequate. This could be stated more precisely.

18. **p179**: The discussion of the pressure term is one of the clearest articulations in the thesis: "The free particle Schrödinger equation describes a fluid that has no internal pressure." This key insight deserves even more prominence — perhaps as a display quote or boxed statement.

19. **p179**: "Short discovered (under the limit of ℏ → 0) that the term is only important in regions of high-density" — this is the quantum pressure / Bohm potential term. The limit ℏ → 0 is the classical limit, and saying the term is "only important in high-density regions" seems backwards. In the ℏ → 0 limit, the quantum pressure term vanishes *everywhere*. For finite ℏ, it's important where ∇²√ρ/√ρ is large, which occurs at density *gradients* (not necessarily high density per se). Worth clarifying.

20. **p180**: "The benefit of the FPA is that it runs much faster than all other known methods of simulation" — same overly strong claim as noted in the Chapter 4 review (#38). Should be qualified: "for the quasi-linear regime."

21. **p180**: The five goals are clearly listed: (1) 3D coordinates, (2) self-consistent gravity, (3) expanding coordinates, (4) periodic boundaries, (5) mass conserving. This is a strong way to demonstrate achievement. However, momentum conservation is mentioned as "a bonus" — but it was presented as an important result in Chapter 5. It deserves more credit here.

22. **p180-181**: "Sabiu's run with Gadget took approximately the same time to complete as my Wave-mechanics code, of the order of a few hours; however, no rigorous speed testing was conducted in this thesis. Both codes were run on different machines so does not represent a fair comparison" — "so does not" should be "so **this** does not." More importantly, the honest caveat about lack of rigorous benchmarking is good, but it somewhat undermines the earlier claim about running "much faster than all other known methods."

23. **p181**: "we also believe that wave-mechanics should provide a better representation of the Universe (over N-body codes) because they describe a continuous density field and suffers less from discretization problems" — this conflates two different things. The wave-mechanics code is *also* discretised (on a grid). The advantage is that the density field is inherently continuous (|ψ|²) rather than a collection of point particles, so there's no shot noise. But the grid discretisation still limits spatial resolution. The claim should be more precise.

24. **p181-182**: The honest discussion of problems (interference effects, low resolution, velocity discrepancy) is appropriate for a conclusion. However, the statement "The messiness of the evolved density fields is potentially due to interference" identifies a fundamental issue — a single coherent wavefunction will always interfere with itself. This is not just a resolution issue but a conceptual limitation of the single-wavefunction approach. The text could acknowledge this more directly.

25. **p182**: "The choice of the ν parameter influences the dynamics in such a way that it can greatly affect the collapse speed; in some cases, collapse can be completely prevented" — this is an important admission. The ν parameter is supposed to be a resolution parameter (analogous to grid spacing), but if it affects the physics (collapse speed), that's more concerning than a simple resolution issue. This deserves more discussion about what the "correct" value of ν should be.

26. **p182**: "Initially, we incurred problems from running low resolution simulations where collapse was too fast in the initial timesteps" — this seems to contradict the previous statement that ν can *prevent* collapse. Both too-fast and too-slow collapse depending on parameters suggests the dynamics are quite sensitive to ν. This tension should be acknowledged.

### Future Work: Multiple Fluids (6.1.1, pp184-185)

27. **p184**: "Johnston has shown how to include another fluid of the same mass into a Schrödinger code" — the phrase "of the same mass" is important but could be expanded. Does this mean the same total mass, the same particle mass (hence same ν), or both? The next sentence clarifies somewhat but the initial statement is ambiguous.

28. **p185**: The multi-species algorithm outline (steps 1-3) is clear and well-presented. One concern: step 2 computes a common gravitational potential V̂₁...N from the combined density of all species. This requires computing |ψ₁|² + |ψ₂|² + ... + |ψ_N|² as the source for the Poisson equation. Worth stating this explicitly.

29. **p185**: "the addition of any term that does not commute with the kinetic energy is very likely to break the energy conservation of the code" — this is a strong and important statement. It effectively says that extending the physics is constrained by the algebraic structure of the splitting operators. Good to flag, but it would help to explain *why* non-commuting additions break energy conservation (it's because the Suzuki splitting is only exact for a specific decomposition structure).

### Future Work: Additional Physics (6.1.2, pp185-186)

30. **p185**: Equation 6.2: exp(K̂')ψ(t) = exp(K̂ + P̂)ψ(t) — defining K̂' = K̂ + P̂ is fine but the text says "I don't believe that it would be entirely correct to modify the kinetic operator to account for pressure." This raises the question: why not? If K and P commute (both being functions of ∇²), then exp(K+P) = exp(K)exp(P) exactly. The issue arises only if they don't commute. Worth clarifying under what conditions K and P commute.

31. **p186**: Equation 6.3: exp(K̂ + V̂ + P̂)ψ(t) → exp(½K̂)exp(½P̂)exp(V̂)exp(½P̂)exp(½K̂)ψ(t) — this is a Suzuki-type decomposition for three operators. However, this particular ordering assumes K and P should be split symmetrically around V, with P̂ half-steps inside K̂ half-steps. Is this the unique correct splitting? With three non-commuting operators, the Suzuki decomposition has multiple possible orderings. The choice may affect accuracy and stability.

### Future Work: Periodic Boundaries (6.1.3, pp186-188)

32. **p187**: Equation 6.4: ψ(r + R, t) = ψ(r, t)e^{iφ}, φ = k·R — this is Bloch's theorem. The text says "k is the Bloch wavenumber and R is the length of the lattice (although Watanabe calls this the unit vector, implying it should be the spacing between gridpoints)." This confusion about whether R is the lattice vector or the grid spacing is significant. In solid-state physics, R is a lattice translation vector (the full box size for a simple lattice). The text should resolve this ambiguity definitively rather than leaving it uncertain.

33. **p187-188**: The confusion about R (lattice size vs. grid spacing) propagates to the interpretation of the phase φ = k·R. If R is the box size L and we want strict periodicity, then φ = 2πn for integer n, making e^{iφ} = 1. If R is the grid spacing Δx, then φ = kΔx, which is a non-trivial phase. The physical meaning changes entirely. This needs clarification.

34. **p188**: Equations 6.5-6.6: The matrix formulation with e^{±iφ} corner elements is clear and well-presented. The key insight — that the tridiagonal structure is broken only by corner elements — is important for the discussion of efficient solvers.

35. **p188**: "∂²_x = ∂²_{x-td} + ∂²_{x-ad}" — decomposing the periodic Laplacian into tridiagonal plus adhesive (corner) contributions is a clever approach. The text notes this can be solved by treating the first and last lines separately. This is essentially the Sherman-Morrison formula for rank-2 updates to tridiagonal systems — worth mentioning by name for readers who want to implement it.

36. **p188**: "An alternative method for periodic boundaries would be the use of Fourier transforms as they are inherently periodic. However, our attempts to implement the kinetic energy operator via an FFT method have been unsuccessful." — This is a tantalising statement. *Why* was the FFT approach unsuccessful? Was it a unitarity issue? A stability issue? An accuracy issue? The reader deserves more than one sentence on this, as an FFT-based kinetic operator is a very natural approach (and is what many later codes use successfully).

### Future Work: Mesh Refinement (6.1.4, pp189)

37. **p189**: The discussion correctly distinguishes between static mesh refinement (pre-binning) and adaptive mesh refinement (AMR). The statement that AMR "is missing in the current wave-mechanics of LSS literature" positions this as an open problem clearly.

38. **p189**: "The potential energy calculation requires separate consideration" — this is an understatement. The Poisson equation on a non-uniform mesh is significantly more complex than on a uniform mesh, requiring multigrid methods or similar. This is arguably the harder problem and deserves more attention.

### Future Work: Parallelisation (6.1.5, pp189-190)

39. **p189-190**: The discussion of parallelisation difficulties is honest and identifies the core problem: the recursion relations in the Goldberg/Cayley scheme are inherently sequential. The staggering idea (overlapping e and f recursions) is a partial mitigation but not a true parallelisation.

40. **p190**: "Future research into this area may have to choose a different Schrödinger solver if such a solver is able to be parallelized while being reliable." — This is an important conclusion. It effectively admits that the Goldberg/Cayley method, while elegant and unitary, may be a dead end for large-scale simulations. The FFT-based split-operator method (which *is* trivially parallelisable) is the obvious alternative — and indeed is what subsequent work (Schive et al. 2014, etc.) adopted.

---

## Equations Check

- **Eq 6.1**: ψ(x, t+dt) = e^{-i(K+V)dt}ψ(x,t) = e^{-iKdt/2}e^{-iVdt}e^{-iKdt/2}ψ(x,t) — correct (Strang splitting)
- **Eq 6.2**: exp(K̂')ψ(t) = exp(K̂ + P̂)ψ(t) — correct (definition)
- **Eq 6.3**: exp(K̂ + V̂ + P̂)ψ(t) → exp(½K̂)exp(½P̂)exp(V̂)exp(½P̂)exp(½K̂)ψ(t) — this is a *specific choice* of Suzuki splitting for three operators. It's second-order accurate but assumes a particular nesting order. The splitting is not unique; one could equally put P̂ outside K̂. The choice should be justified (e.g., by computational convenience or commutativity properties).
- **Eq 6.4**: ψ(r+R, t) = ψ(r,t)e^{iφ}, φ = k·R — correct (Bloch's theorem)
- **Eq 6.5-6.6**: Matrix form with e^{±iφ} corners — correct (periodic tridiagonal system with Bloch phase)

---

## Structural Comments

41. The chapter has a natural two-part structure: retrospective (summary + honest assessment) and prospective (future work). Both parts are well-executed.

42. The summary section effectively recaps all major chapters without being overly repetitive. The honest acknowledgment of problems (interference, resolution, velocity discrepancy) adds credibility.

43. The future work section is one of the strongest in the thesis. Rather than vague hand-waving about "future improvements," it presents concrete ideas with mathematical formulation (Eqs 6.1-6.6) and practical analysis of feasibility. The discussion of adhesive operators and parallelisation challenges is particularly valuable.

44. **Missing from the conclusion**: The unresolved velocity scaling discrepancy with GADGET (discussed extensively in Chapter 5) is mentioned only briefly. Given that this was one of the most significant unresolved issues, it deserves a dedicated paragraph in the conclusion discussing possible causes and how future work might address it.

45. **Missing from the conclusion**: No discussion of convergence testing — how does the solution change as resolution increases? The brief mention of 128³ tests looking "better behaved" is insufficient. A systematic convergence study would be important for validating the method.

46. Pronoun inconsistency continues: "I" (pp180, 183, 185) vs. "we" (pp178, 180, 181, 182) with no clear pattern. The conclusion should be consistent throughout.

47. The future work section is somewhat long relative to the summary (8 pages vs. 5 pages). For a conclusion chapter, the retrospective analysis could be expanded (particularly regarding the unresolved issues) and the future work slightly condensed.
