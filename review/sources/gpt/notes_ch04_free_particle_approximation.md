# Chapter 4: Free Particle Approximation - Review Notes

## Overall Assessment
This chapter is well structured and methodologically strong. The progression from FPA theory to 1D verification, then 3D toy tests, then a 3D cosmological comparison against Hydra is clear and persuasive.

The main issues are not in the broad structure but in a few equation-level presentation problems, several over-strong numerical/physical claims, and a handful of obvious copy-editing errors.

---

## Higher-Priority Technical Issues

1. **p87**: "time evolution is unitary (that is, it conserves energy)" is too strong / not reliably stated. Unitarity guarantees norm conservation, but it should not be casually equated with exact energy conservation in that way.

2. **p88 / Eq. 4.3**: The FPA condition is central, but the quantity `f` in
   `V = Φ_g - 3Ω_cdm φ_v / (2f^2 D)`
   is not defined here. Since the argument of the whole approximation depends on this relation, `f` needs to be introduced explicitly.

3. **p89 / Eq. 4.5**: The Fourier-evolution formula appears with a leading minus sign in front of the inverse-transform prefactor. That sign is unusual and looks like a likely typo unless a specific convention is being used and explained elsewhere.

4. **p89 vs Chapter 3**: The velocity-potential sign convention appears inconsistent. Here the chapter uses `v = -∇ϕ_v`, while Chapter 3 used `v = ∇ϕ` in the Madelung discussion. If these are intentionally different conventions, that needs to be stated explicitly.

5. **p95 / Eq. 4.17**: The RMS formula for `V_diff` appears wrong as printed. It is missing the outer square root and instead has a prefactor `1 / n^(1/2)` multiplying the sum. This should be corrected to match the earlier RMS definition.

6. **p95 / Figure 4.3 caption**: The caption calls the plotted quantity `ρ(V_e, V_c)` while the text says it is the RMS deviation of the difference of the two velocities. `ρ` was introduced earlier as the correlation coefficient, so the caption appears inconsistent with the text.

7. **p99 / Eq. 4.21**: The displayed formula for the initial gravitational potential is unclear:
   `φ_g,init(x,y,z) = ∫∫ 4πG δ_init(x,y,z) dV`
   This is not a standard or clearly correct way to write the Poisson solution. If the actual computation is done in Fourier space via `φ_k ~ δ_k / k^2`, the equation should say that directly.

8. **p101**: The text says the correlation "steadily decreases but then improves for the last comparison," but the listed values are not monotonic in a simple way and no explanation is attempted. That is acceptable as honesty, but it leaves a noticeable gap in the interpretation of the main comparison metric.

9. **p108**: "The benefit of the FPA is that it runs much faster than all known N-body codes" is too broad. This should be qualified to the regime and comparison being made.

10. **p108**: The observed turn-over in Figure 4.10 is flagged as an artefact from wavefunction construction, but this is an important enough discrepancy that it deserves either a little more analysis or a more cautious statement.

---

## Typos and Grammatical Issues

11. **p88**: "the the velocity potential" contains a repeated "the".

12. **p88**: "Structure can not form" should normally be "cannot form".

13. **p89**: "In should be obvious" should be "It should be obvious".

14. **p91**: `1 >> δ_a > 0` would read more naturally as `0 < δ_a << 1`.

15. **p95**: "V_e is myself and V_c is Chris Short" should be something like "V_e is from my code and V_c is from Chris Short's code."

16. **p96**: "Zel'dolvich" is misspelled.

17. **p97**: "straight forward" should be "straightforward".

18. **p100**: "Hydra (v 4.2.1)" is unclear as written. It looks like either a typo or an awkward version reference.

19. **p105**: "the overall behaviour of the simulations are at least believable" should be "is at least believable".

20. **p105**: "hence we do not suspect that something is terribly wrong" is too colloquial for a thesis.

21. **p106**: "the presences of vortices" should be "the presence of vortices".

22. **p106-107**: "If we zoom into to look closer" contains a duplicated "to".

23. **p107**: There is an orphan sentence at the bottom of the page: "see the rest of the velocity vectors 'circling' around the wavefunctions null-point." It appears to have been detached from the preceding paragraph.

24. **p108**: "Schrödinger equation." also appears as an orphan fragment following the Figure 4.10 discussion.

25. **p109**: "frame work" should be "framework".

26. **Throughout**: PDF metadata fields such as `Title` and `Author` are blank here as well.

---

## Clarity / Interpretation Notes

27. **p87-88**: The explanation of why the effective potential can be set to zero is the key physical idea of the chapter, but it is somewhat buried in prose. It would benefit from a slightly more prominent presentation.

28. **p89**: "`|ψ|^2` defines the density field" is clear, but the notation around `ρ`, `δ`, and spatial averages could be kept slightly tighter to avoid switching too quickly between density and density contrast.

29. **p90**: The growth-factor subsection is useful, but the practical role of `D` in the rest of the chapter could be stated even more directly: it is the effective time variable for the FPA comparisons.

30. **p92**: The statement that large `ν` "prevents collapse and hence appears to prevent singularities from forming" is plausible in the intended sense, but it should be phrased more carefully because the chapter is moving between numerical robustness and physical interpretation.

31. **p96-97**: The claim that the 3D extension is simple is fair for the FPA, but it would help to say explicitly that this simplicity is specific to the free-particle case and does not carry over to the full Schrödinger-Poisson solver.

32. **p98-99**: The initial-condition construction is generally clear. The computational algorithm list is helpful and could perhaps be used a little earlier as a signpost for the section.

33. **p104-105**: The near-identical velocity dispersions in the three Cartesian components are a strong consistency check and could be highlighted more positively.

34. **p106-107**: The vorticity section is interesting, but the interpretation remains somewhat tentative. Since the chapter itself frames these vortical points as diagnostics, the prose should stay careful about what is established versus merely suggestive.

35. **Pronouns throughout**: The chapter shifts between `I` and `we` without an obvious reason. It would be better to choose one style and keep it consistent.

---

## Structural Notes

36. The overall chapter structure is effective and probably one of the stronger pieces of the thesis so far.

37. Section 4.2.4 (consistency checks) is valuable and strengthens the later Hydra comparison. It might even work slightly earlier as validation of the inputs before the main comparison.

38. Section 4.3.1 does a good job of motivating Chapter 5 by making the FPA limitations concrete.
