# Chapter 5: Solving the full Schrödinger-Poisson system - Review Notes

## Overall Assessment
This is the core chapter of the thesis and the main original contribution. It is also the most methodologically interesting chapter so far: the numerical scheme, periodic-boundary treatment, cosmological generalisation, and the extensive test suite all add real value.

The chapter is strongest when it is explicit about implementation details and when it is honest about unresolved problems. The main weaknesses are a handful of important equation-level issues, several claims about conservation/stability that are stated too strongly, and a few places where the results section slides from evidence into interpretation a bit too quickly.

---

## Higher-Priority Technical Issues

1. **p113-114**: The discussion of the matrix structure looks wrong in basis terms. In position space, the potential operator is diagonal, while the kinetic finite-difference operator is banded. The text says the potential term "will contain off diagonal elements," which appears backwards unless another basis is intended and stated explicitly.

2. **p118**: The claim that the Cayley-based scheme used here is "implicit but conditionally stable" needs much clearer justification. The standard Cayley / Crank-Nicolson style treatment of the linear Schrödinger equation is usually introduced as unconditionally stable. If the use of `V^n` instead of a truly implicit potential update changes that, the chapter should say so explicitly and explain why.

3. **p119-120**: The scheme is described as implicit, but the advanced potential is replaced by the current-time potential. That makes the practical implementation at least partly semi-implicit / approximate. This is not necessarily wrong, but the distinction should be made clearer.

4. **p130 / Eq. 5.60**: The Poisson equation is written as `∇²V = 4πG ψψ*`. For a periodic box, the mean source term must be handled carefully. Elsewhere in the thesis the mean density is subtracted; here that subtraction is absent. This needs explicit clarification because a periodic Poisson problem with non-zero mean source is not well posed in the usual way.

5. **p130 / Eq. 5.61**: `V = ∫ 4πG ψψ* dV` is not a correct standalone integral solution of the Poisson equation. If the intent is a Green-function solution, that kernel is missing. If the actual numerical method is Fourier-space inversion, the equation should say that directly instead of presenting this misleading integral form.

6. **p130-131 / Eqs. 5.62-5.63**: The Fourier-space Poisson solution is written in terms of `ψ_k ψ_k*`. That is the wrong object if the source is `ρ = |ψ|^2` in real space. What is needed is the Fourier transform of the density field, not the product of the Fourier transforms of `ψ` and `ψ*` as written here.

7. **p132 / Eq. 5.64**: The normalisation of the 1D Gaussian wavefunction appears wrong. The prefactor written is not the usual one for a properly normalised Gaussian packet.

8. **p149-150**: The tophat-collapse discussion appears internally inconsistent about the role of `ν`. One passage says larger `ν` makes collapse happen faster, while nearby discussion says larger `ν` increases dispersion and can push the wavefunction across boundaries. That interpretation needs to be reconciled more carefully.

9. **p165-167**: The claim that mass and momentum conservation together imply energy conservation is too strong. The displayed checks show mass and the quantity called momentum staying flat, but that is not sufficient on its own to conclude energy conservation.

10. **p167 / Eq. 5.69**: The definition of `P_total` as `∫ |ψ* ∇ψ|^2 dx^3` does not look like the standard total momentum. This needs checking carefully, because the subsequent conservation claim depends on it.

11. **p168-170**: The velocity discrepancy is a major unresolved problem, not a minor caveat. The chapter is commendably honest about it, but from a review perspective this is one of the most serious issues in the thesis because it affects the physical interpretation of the cosmological results.

12. **p170**: The non-zero mean velocities growing over time imply a spurious bulk flow. Since this is explicitly unphysical, it should be treated as a major unresolved defect in the current implementation rather than just an open question for future work.

13. **p174-177**: The chapter summary is admirably candid, but it slightly overstates what has been established numerically. The mass-conservation result is strong; the cosmological validity of the velocity field is much less secure.

---

## Typos and Grammatical Issues

14. **p112**: "This requirements" should be "These requirements".

15. **p112**: "a exponential" should be "an exponential".

16. **p124 / p127**: "Widow & Kaiser" is a typo for "Widrow & Kaiser".

17. **p128**: "the the operator" contains a repeated "the".

18. **p132**: "Through out" should be "Throughout".

19. **p139**: "pass significantly pass the boundary" is ungrammatical.

20. **p149**: `(Coles, (Coles 1997))` contains a citation-formatting error.

21. **p157-158**: "there is a an accepted standard" has repeated articles.

22. **p170**: "as sen in section 4.2.5" should be "as seen in section 4.2.5".

23. **p174**: "Given sufficiently enough time" should be "Given enough time" or "Given sufficient time".

24. **p174**: "how implement a wave-mechanics code" should be "how to implement".

25. **p175**: "I’m less sure" is too conversational in register for a thesis chapter.

26. **p177**: "beyond all doubts" should normally be "beyond all doubt".

27. **Throughout**: PDF metadata fields such as `Title` and `Author` are blank here as well.

---

## Clarity / Interpretation Notes

28. **p111-112**: The five requirements are a good organising device. They work well and could be referenced more explicitly later when the results are summarised.

29. **p114-115**: The Cayley-transform discussion is interesting, but some of the group-theory framing is more abstract than the chapter needs. The practical numerical consequences are the important part here.

30. **p115**: The statement that the factors on the right-hand side of the Cayley transform commute should be justified briefly. They commute because they are both functions of the same matrix, which is not obvious to every reader.

31. **p116-117**: The move from the naive higher-dimensional Goldberg extension to Suzuki splitting operators is one of the key ideas of the chapter and deserves even more emphasis than it gets.

32. **p121-123**: The periodic-boundary construction is one of the strongest parts of the chapter. The argument is clear and Figure 5.1 supports it well.

33. **p124-126**: The analysis of the Widrow-Kaiser inconsistency is a genuine contribution, but this section is a little long. It could be tightened without losing substance.

34. **p127 / Eq. 5.54**: The generalisation to flat models with non-zero cosmological constant looks important and should perhaps be highlighted more strongly as one of the chapter’s main original results.

35. **p133**: The switch from Hydra in Chapter 4 to GADGET-2 here is reasonable, but the reason for that change could be stated more explicitly.

36. **p153-156**: The chapter is very honest about the unsatisfactory lower-resolution cosmological runs. That is good science. It would help to separate more clearly what is believed to be a resolution issue versus what may be a deeper issue with the wave-mechanics framework or its implementation.

37. **p155-164**: The density and overdensity comparisons are interesting, but the interpretation sometimes outruns the evidence. In particular, claims that certain behaviours are "inevitable" or "fundamental" to wave-mechanics should be softened unless they are demonstrated more directly.

38. **p158-166**: The Gaussian smoothing section is useful and clearly improves the visual comparability to GADGET. It would be good to say a little more explicitly whether the smoothing is being used as a diagnostic visual aid or as part of the physically relevant comparison.

39. **p168-170**: The discussion of the velocity discrepancy is one of the most important parts of the whole chapter. It is handled honestly, but it should probably be framed as a central limitation of the current results rather than a secondary anomaly.

40. **p174**: The statement that the full S-P results show no signs of vorticity because `ψ` does not get close enough to zero is interesting, but the threshold and diagnostic logic are still somewhat qualitative.

---

## Structural Notes

41. The chapter is long, but the internal progression mostly works: method first, then boundary conditions, expansion, Poisson solver, algorithm, initial conditions, and finally results.

42. The results section is strong because it builds from simple tests to the full cosmological run. That progression makes the chapter much easier to trust.

43. The most convincing original contributions in this chapter are:
   the periodic-boundary treatment,
   the careful implementation discussion,
   the cosmological generalisation beyond Einstein-de Sitter,
   and the explicit conservation checks.

44. The main residual risk after reading the chapter is that the density-side results look promising while the velocity-side results still have a serious unresolved defect. That asymmetry should be made very explicit whenever the chapter’s success is summarised.
