# Chapter 2: Review of Numerical Simulations - Review Notes

## Overall Assessment
This chapter gives a useful high-level survey of standard simulation methods and does a good job motivating why wave-mechanics was considered as an alternative.

The main concern is that the section on symplectic methods appears to contain several mathematically inaccurate or overstated claims. Outside that section, most issues are copy-editing, clarity, or occasional technical imprecision.

---

## Higher-Priority Technical Issues

1. **Section 2.2**: "The word symplectic is synonymous with the word Hamiltonian" is too strong and not really correct. Symplectic geometry is closely related to Hamiltonian mechanics, but the terms are not synonymous.

2. **Section 2.2**: "As the phase-space volume is conserved then the Hamiltonian (total energy of the system) is also conserved" is not generally true. Volume preservation and exact energy conservation are distinct properties. This should be corrected.

3. **Section 2.2**: "all symplectic manifolds ... are simply connected" is false. This is a substantive mathematical error.

4. **Section 2.2**: "the Poincaré group, of which the Abelian (symplectic) group is a subgroup" is not a coherent group-theoretic statement. This needs either a careful rewrite or removal.

5. **Section 2.2**: The discussion of symplectic integrators seems to overstate what they guarantee. In practice, symplectic integrators preserve the symplectic structure and usually have good long-term energy behaviour, but they do not in general exactly conserve the Hamiltonian.

6. **p40**: The claim that Sundman's regularization "assumes zero angular momentum, a fair assumption for gravity which is a radial force" is misleading. Gravitational systems can have non-zero angular momentum, and in fact most interesting orbital systems do.

---

## Typos and Grammatical Issues

7. **p36**: "thousands or millions times denser" should be "thousands or millions of times denser".

8. **p37**: "what they did was to investigate ... and has come to be known" should be reworked. The verb agreement is off.

9. **p37**: "one of the first billion-body (N = 10243 ) gravitational simulations" is unclear as written. This likely intends `1024^3`, not `10243`.

10. **p37**: "more particles allows for higher resolution" should be "more particles allow higher resolution" or similar.

11. **p38**: "under go classical gravitational attraction" should be "undergo classical gravitational attraction".

12. **p38**: "It is more accurate to say that the computational particle represents collection of real particles" needs an article: "represents a collection...".

13. **p38**: "only account for the gravitational force all other physics is essentially smoothed out" needs punctuation or rewriting.

14. **p39**: "Consequently, the minimum useful mass resolution is 10^9..." is awkwardly phrased. The parenthetical explanation is also hard to parse.

15. **p41**: "it is not good computational practise" should be "practice" in standard modern usage here.

16. **p41**: "To avoid these problems Aarseth introduced a regulator, that is a minimum distance (ϵ) for impact ϵ" is awkward and repetitive.

17. **p42**: "The softening parameter is non-derivable" is unusual phrasing. "Phenomenological" or "introduced ad hoc" would be clearer.

18. **p42**: "larger enough" should be "large enough".

19. **p43**: "Direction summation methods" should be "Direct summation methods".

20. **p46**: "The word symplectic is a construction by Weyl" reads awkwardly. Better would be "The term 'symplectic' was introduced by Weyl".

21. **p47**: "a differentiable symmetry" is not the usual phrasing in the Noether statement. The sentence should be checked for accuracy and wording.

22. **p48**: "short coming" should be "shortcoming".

23. **p50**: "As the N-body discrete particles to model a fluid" is missing a verb, likely "uses".

24. **p50**: "While phase-space methods deal with a continuous distribution function and avoid two-body relaxation they explicitly deal..." needs a comma after "relaxation".

25. **Throughout**: PDF metadata fields such as `Title` and `Author` are blank here as well.

---

## Clarity / Technical Precision Notes

26. **p37**: The sentence about GADGET-2 being "a mix of the direct summation and particle mesh methods" is broadly fine, but the wording around P3M / tree / SPH could be cleaner. As written, it conflates several algorithmic components.

27. **p39**: "the minimum useful mass resolution is 10^9" seems to refer to dynamic range rather than resolution in the usual sense. This deserves clearer wording.

28. **p40**: The discussion of Sundman and Wang is interesting, but it is slightly too compressed to be reliable without very careful phrasing. I would check the historical/mathematical details here.

29. **p44**: The sentence "the density ρ = ρ̄δ" is imprecise if `δ` is the density contrast in the standard sense, where usually `ρ = ρ̄(1 + δ)`. This should be checked against the notation used elsewhere in the thesis.

30. **p45**: "The dominant calculation is still the Particle-Mesh force routine" may be true in the intended context, but it reads more like a broad claim than a chapter-local statement. It may be worth softening unless specifically sourced.

31. **p46-47**: The symplectic discussion shifts from a reasonable motivation for leapfrog to a very abstract mathematical digression. It may be more effective to keep this section practical and directly tied to numerical integration.

32. **p48**: "Most perturbation modes will evolve in a linear manner" is a little too sweeping unless clearly restricted to early times / large scales.

33. **p50**: "He evolved the Euler equation for density and velocity fields (see equation 2.6)" is imprecise because equation 2.6 includes continuity and Poisson equations as well.

---

## Structural Notes

34. **Section 2.2** feels disproportionate relative to the rest of the chapter. The abstract mathematical treatment of symplectic manifolds is much more formal than the surrounding review, and currently it is also the least reliable section.

35. The chapter works best when it stays close to simulation practice. A shorter, more concrete explanation of why leapfrog is preferred in N-body work would likely improve the flow.
