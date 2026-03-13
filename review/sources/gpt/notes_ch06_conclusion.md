# Chapter 6: Conclusion - Review Notes

## Overall Assessment
This is a good conclusion chapter overall. It does two jobs clearly: it summarises the thesis and it proposes concrete future directions rather than vague generalities.

The main review concern is not missing structure but overstatement. A few claims in the summary section are stronger than the evidence established in the earlier chapters, especially given the unresolved velocity/bulk-motion problem and the still-ambiguous interpretation of the pressure term.

---

## Higher-Priority Content Issues

1. **p178-179**: The conclusion repeats the Chapter 3 pressure-term interpretation very confidently:
   "The free particle Schrödinger equation describes a fluid that has no internal pressure."
   This is too strong given the earlier ambiguity. The thesis has not fully settled that point, and the conclusion should not state it more confidently than Chapter 3 actually justified.

2. **p179**: "Short discovered (under the limit of `ħ -> 0`) that the term is only important in regions of high-density" is not stated carefully enough. In the strict `ħ -> 0` limit the quantum-pressure term vanishes everywhere. For finite `ħ`, its importance is tied to gradients/curvature in the density field, not simply "high density" as such.

3. **p180**: "The benefit of the FPA is that it runs much faster than all other known methods of simulation" is too broad. Earlier chapters did not establish that level of benchmarking, and Chapter 5 explicitly says rigorous speed testing was not done.

4. **p180-181**: The conclusion says the main goals of the full code were achieved and that the simulations "prove to be robust and stable." That is only partly true. Mass conservation and several test problems are strong evidence, but the cosmological velocity field still has a serious unresolved defect. The summary should make that limitation more explicit.

5. **p180**: The comparison with GADGET runtime is not rigorous enough to support any general speed claim. The text itself admits the runs were on different machines and that no proper benchmarking was performed, so this should be framed more cautiously.

6. **p181**: The claim that wave-mechanics should provide a "better representation of the Universe" than N-body codes is too strong as written. The continuous density field is an advantage, but the wave-mechanics code is still discretised on a grid and still has unresolved issues in the present implementation.

7. **p181-182**: The conclusion does mention the messy density fields and possible interference, but it could be more direct that this may be a conceptual limitation of using a single coherent wavefunction for a classical cosmological system, not merely a low-resolution nuisance.

8. **p182**: The `ν` parameter is treated as both a numerical/resolution parameter and something that materially affects collapse behaviour. That tension is important and should be acknowledged more explicitly in the conclusion.

9. **p182**: The chapter says higher-resolution runs looked better behaved, but without a systematic convergence study this remains suggestive rather than decisive. The conclusion would be stronger if it stated that more carefully.

10. **p183-190**: The future-work section is strong, but it could more clearly distinguish between:
    straightforward engineering extensions,
    mathematically nontrivial extensions,
    and ideas that are still speculative.

---

## Typos and Grammatical Issues

11. **p179**: "Johnston, considered a different approach..." has an unnecessary comma after `Johnston`.

12. **p180**: "colloquially we call it the Cayley method" would read better as "we refer to it as the Cayley method."

13. **p180**: `(3+1 d)` should be `(3+1)D` or `(3+1)-dimensional`.

14. **p181**: "Through out this work" should be `Throughout this work`.

15. **p181**: "they describe a continuous density field and suffers less from discretization problems" should be `suffer less`.

16. **p182**: "the unexpected, although not necessarily inconsistent, results for density and velocity" is awkwardly punctuated and would read better with em-dash-style interruption or a recast sentence.

17. **p182**: `gaussian` / `gaussian filter` should be capitalised as `Gaussian`.

18. **p184**: "they are less conventional and hence far more speculative" would benefit from punctuation revision; currently it reads as a splice.

19. **p186**: Contractions such as `I don't` are a little informal for the tone of the thesis.

20. **p190**: "the next second processor" is unclear phrasing.

21. **Throughout**: PDF metadata fields such as `Title` and `Author` are blank here as well.

---

## Clarity / Interpretation Notes

22. **p178**: The opening summary paragraph is useful, but the sentence "The quantum nature of dark matter particles is mostly unknown but not thought to be significant" could use more careful phrasing. The thesis is really arguing that the present modelling choice does not rely on that quantum nature being dynamically important.

23. **p180**: The statement that momentum conservation comes "as a bonus" undersells one of the more persuasive validation checks from Chapter 5.

24. **p181**: The contrast with N-body methods is one of the conceptual strengths of the thesis, but the conclusion should avoid implying that wave-mechanics has already outperformed N-body in practice. The thesis has shown promise, not dominance.

25. **p183-184**: The future-work list is good and concrete. It helps that the ideas are tied to specific technical mechanisms rather than just broad themes.

26. **p184-185**: In the multiple-fluids discussion, it would help to say explicitly that the shared gravitational potential must be sourced by the combined density of all species.

27. **p185-186**: The discussion of adding pressure via an additional split operator is sensible, but the conclusion should be careful not to imply that this decomposition is obviously unique or straightforward for arbitrary new physics.

28. **p186-188**: The adhesive-operator discussion is interesting and technically valuable. The ambiguity around the Bloch/lattice notation could be made clearer if this section is kept.

29. **p189-190**: The parallelisation section is one of the better future-work discussions because it is candid about the possibility that the current solver may simply not be the right long-term route for very large simulations.

---

## Structural Notes

30. The chapter’s two-part structure works well: retrospective summary first, then future work.

31. The future-work section is arguably stronger than the summary section because it is more precise and more technically grounded.

32. The summary would benefit from one dedicated paragraph on the unresolved velocity-scaling / bulk-motion issue, since that is one of the most important open problems left by the thesis.

33. The conclusion would also be stronger with a brief explicit statement of what is solidly established versus what remains provisional.
