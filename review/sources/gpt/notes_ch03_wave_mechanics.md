# Chapter 3: Wave-mechanics - Review Notes

## Overall Assessment
This is the most conceptually important chapter so far. It lays out the interpretation of the Schrödinger-Poisson system, the Madelung and Bohm connections, the literature review, the initial-condition machinery, the velocity methods, and the vorticity discussion.

The chapter is strong in scope and ambition, but it also contains a handful of important technical statements that look wrong or overstated, especially in the Madelung/pressure discussion and in the numerical-operator discussion. There are also a number of smaller typos and phrasing issues.

---

## Higher-Priority Technical Issues

1. **p53**: "the norm of the function, `|ψ| = ψψ*`" is incorrect. `ψψ* = |ψ|^2`, not `|ψ|`. This should be fixed because it affects the basic density interpretation of the wavefunction.

2. **p58**: The conclusion "Free Particle Schrödinger = Fluid - Pressure" followed by "Conclusively, the Schrödinger-Poisson system is pressure-free" is too strong and, as written, misleading. In the standard Madelung picture, the Schrödinger equation leads to fluid-like equations with an additional quantum-pressure / quantum-potential term. The usual statement is that one recovers pressureless dust only when that term is neglected or becomes negligible, not that the Schrödinger-Poisson system is intrinsically pressure-free.

3. **p79 / Eq. 3.32**: The displayed non-commutation statement is wrong as written:
   `e^{-i(K+V)dt} != e^{-i(V+K)dt}`
   Since `K + V = V + K` as algebraic sums, both sides are the same expression. What is presumably intended is a statement about split operators, e.g. `e^{-iKdt}e^{-iVdt} != e^{-iVdt}e^{-iKdt}`.

4. **p79**: "The Schrödinger equation naturally conserves energy and momentum but this property is only preserved when a unitary solver is employed as such a method is inherently symplectic." This is too compressed and not reliably stated. Unitarity, symplecticity, norm conservation, and exact energy conservation are related but not interchangeable; this should be rewritten more carefully.

5. **p83 / Eqs. 3.43-3.44**: The phase-angle method uses `arctan(Im/Re)` in a way that ignores quadrant ambiguity unless extra care is taken. Computationally this should really refer to `atan2`-type logic, or at least explicitly discuss the quadrant issue. Since this section is about practical numerical methods, this matters.

6. **p85-86**: The vorticity discussion is too categorical in places. "Extrinsic angular momentum ... is forbidden by Kelvin's circulation theorem" overstates the result. Kelvin's theorem depends on restrictive assumptions, and the connection between cosmological vorticity, shell crossing, singular points of `ψ`, and numerical diagnostics needs a more careful wording.

7. **p62**: The statement that quantum mechanics has difficulty describing systems below `10^-13 m` and that experiments do not agree as well there is historically loaded and reads as if it were an accepted modern limitation. That needs checking or contextualising much more carefully.

8. **p60-61**: The discussion of the classical limit `ħ -> 0` is not stated carefully enough. Saying the Schrödinger equation "breaks down" in that limit is too blunt; the issue is subtler than that and is usually handled through asymptotics / WKB / semiclassical limits rather than literal substitution.

---

## Typos and Grammatical Issues

9. **p55**: "enumeration of the- particle number" contains a stray hyphen.

10. **p57**: "it is may not be wise" should be "it may not be wise".

11. **p60**: "break down as a the fluid's pressure" should be "break down as the fluid's pressure".

12. **p69**: "under going gravitational collapse" should be "undergoing gravitational collapse".

13. **p72**: "sesqui-linear quantity `ψψ*`" is unusual phrasing here and likely more confusing than helpful. "modulus squared" or "norm squared" would be clearer.

14. **p77**: "This paper includes has some similarities" should be "This paper has some similarities".

15. **p78**: "the first is an an implicit method" contains a repeated "an".

16. **p84**: "This latter point would should a singularity" should be "would show a singularity".

17. **p85**: "Appendix B.2 shows the consistency of the the probability current..." contains a repeated "the".

18. **Throughout**: PDF metadata fields such as `Title` and `Author` are blank here as well.

---

## Clarity / Interpretation Notes

19. **p53**: The role of `ħ` / `ν = ħ/m` is introduced somewhat ambiguously. The text shifts between quantum-mechanical and effective-classical interpretations. Since this ambiguity is central to the thesis, it would help to make the intended interpretation explicit earlier.

20. **p53**: "`p^2` appears directly in the Hamiltonian and is interpreted as the momentum" is imprecise. The momentum operator is `p`; `p^2` is the kinetic-energy contribution.

21. **p54-55**: The section on second quantisation and QFT is interesting, but it is longer than the role it finally plays in the thesis. It may be worth tightening unless the later chapters rely more strongly on it.

22. **p56**: "ansätze form" is awkward. "ansatz" or "ansatz form" would read better.

23. **p57-59**: The Madelung derivation would benefit from a clearer separation between exact transformed equations and later interpretive claims. At present, the algebra and the physical interpretation are interleaved in a way that makes the logic harder to follow.

24. **p60**: The claim that the pressure term is "only important in the domain of shell crossing" should probably be softened unless this is explicitly established for the approximations used.

25. **p64-77**: The literature review is comprehensive and valuable, but there are occasional places where interpretive comments are presented more confidently than the chapter’s own evidence supports. The Woo & Chiueh discussion in particular could use a more careful separation of what they claim versus what you infer.

26. **p75**: The discussion of operator splitting is useful, but it would help to mention the standard symmetric / Strang-style split explicitly when contrasting methods.

27. **p76**: The units attached to `ν = ħ/m` are not clearly handled in the discussion of extremely light dark matter. The dimensional analysis should be checked.

28. **p80-82**: The initial-conditions section is solid overall. The narrative would be slightly cleaner with one explicit sentence connecting the inflationary spectrum, transfer function, and random-phase realisation in the actual generator workflow.

29. **p82**: The notation around the phase potential is a little confusing because `ϕ`, `ϕ_v`, and sign conventions are doing a lot of work. A brief explicit convention statement would help.

30. **p83-85**: The probability current method is one of the strongest parts of the chapter. It is clear and practically motivated. If anything, this section could be given slightly more prominence relative to the Bohmian detour.

31. **p85**: "This is a statement of zero angular momentum" is too strong. Irrotational flow does not simply collapse to a statement about total angular momentum in that way.

---

## Structural Notes

32. The chapter covers a lot of ground and mostly holds together, but it is dense. A short bridging paragraph between the interpretation sections and the literature review would help.

33. The Bohmian mechanics section is interesting, but its payoff for the rest of the thesis is mainly the comparison with the quantum-potential term. That purpose could be stated more directly at the start of the section.

34. The literature review is genuinely useful as a standalone reference and is one of the strengths of the chapter.
