# Chapter 6: Conclusion - Reconciled Review Notes

## Source Coverage

- GPT review: `gpt-review/notes_ch06_conclusion.md`
- Claude review: `claude-review/notes_ch06_conclusion.md`
- Overall result: very strong agreement; the chapter is good overall, but the retrospective summary overstates some conclusions

## Shared Findings

The two reviews independently agree on the following main points.

1. **This is a solid conclusion chapter overall.** Both reviews think the chapter succeeds in doing two jobs:
   summarising the thesis,
   and outlining concrete future work.

2. **The future-work section is one of the strongest parts of the chapter.** Both reviews praise it for being specific, technically grounded, and more useful than generic future-work lists.

3. **The main problem is overstatement in the summary section.** Both reviews independently say that some retrospective claims go beyond what the earlier chapters actually established.

4. **The pressure-term summary is too confident.** Both reviews independently object to the conclusion treating the free-particle Schrödinger equation as straightforwardly pressure-free, given the earlier ambiguities in Chapter 3.

5. **The speed claims are too broad.** Both reviews agree that the conclusion overstates the runtime advantage of the FPA and the full code, especially given the lack of rigorous benchmarking.

6. **The summary of the full code is too confident given the unresolved velocity problem.** Both reviews say the conclusion should state more explicitly that the cosmological velocity field remains a major unresolved limitation.

7. **The comparison with N-body methods needs more care.** Both reviews think the conclusion slides too easily from "wave-mechanics has advantages" to "wave-mechanics gives a better representation of the Universe," which is stronger than the thesis has shown.

8. **The `ν` parameter tension should be acknowledged more directly.** Both reviews note the unresolved ambiguity between treating `ν` as a numerical/resolution parameter and observing that it materially changes collapse behaviour.

9. **The chapter would be stronger with a more explicit limitations paragraph.** Both reviews independently say the unresolved velocity/bulk-motion issue deserves its own dedicated treatment in the conclusion.

10. **There are multiple straightforward copy-editing issues.** Both reviews independently flag:
   the comma after `Johnston`,
   awkward use of "colloquially",
   `(3+1 d)`,
   `Through out`,
   singular/plural agreement problems,
   lower-case `Gaussian`,
   and pronoun inconsistency.

11. **Pronoun inconsistency remains a thesis-wide issue here too.** Both reviews note continued shifting between `I` and `we`.

## GPT-Only Additions Worth Keeping

These points appear in the GPT review but not the Claude review and are worth preserving.

1. **`p179`**: The statement about the `ℏ → 0` limit and the pressure term being important only in high-density regions is not carefully phrased. The issue is density gradients/curvature, not simply density itself.

2. **`p181-182`**: The messy density fields and interference patterns could be framed more directly as a possible conceptual limitation of using a single coherent wavefunction, not just as a resolution nuisance.

3. **`p182`**: The conclusion should acknowledge more explicitly that better behaviour at higher resolution is suggestive but not decisive without a systematic convergence study.

4. **`p183-190`**: The future-work section could benefit from a clearer distinction between:
   straightforward engineering extensions,
   mathematically difficult extensions,
   and genuinely speculative ideas.

5. **`p180`**: Momentum conservation is described as "a bonus," but the conclusion should give it more credit as one of the stronger validation checks.

6. **Metadata**: the PDF metadata fields are blank here as well.

## Claude-Only Additions Worth Keeping

These points appear in the Claude review but not the GPT review and are worth preserving.

1. **`p178`**: The opening statement that previous inconsistencies were resolved "to the best of our ability" is somewhat hedged and could be rewritten more cleanly.

2. **`p178`**: The argument about the quantum nature of dark matter versus the practical modelling choice should be stated more precisely.

3. **`p180`**: The five-goals recap works well and should perhaps be expanded to give momentum conservation more prominence.

4. **`p181`**: The discretisation comparison with N-body should be made more precise: wave-mechanics still uses a grid, so the advantage is not the absence of discretisation per se, but the continuous density representation and reduced shot noise.

5. **`p184-185`**: In the multiple-fluids discussion, the common gravitational potential should be stated explicitly as being sourced by the sum of all species densities.

6. **`p185-186`**: The pressure/split-operator future-work idea should be more careful about whether the proposed operator decomposition is unique or obviously correct.

7. **`p187-188`**: The Bloch/lattice-vector notation in the adhesive-operator discussion is ambiguous enough that the text should resolve it rather than leaving it open.

8. **`p188`**: The FFT alternative for periodic boundaries is mentioned but dismissed too quickly; the reader needs more explanation of why those attempts failed.

9. **`p189`**: The Poisson problem on non-uniform meshes is a major obstacle for AMR and deserves a more explicit statement.

10. **`p190`**: The parallelisation section rightly suggests that the current Goldberg/Cayley solver may not be the right long-term route for very large simulations; this is an important conclusion and worth keeping prominent.

11. **Missing from the conclusion**: a convergence-study discussion should be mentioned explicitly.

## Reconciled Priorities

1. Tone down the retrospective claims first:
   pressure-free interpretation,
   speed/performance claims,
   robustness/stability claims,
   and any suggestion that wave-mechanics has already outperformed N-body methods.

2. Add one explicit paragraph on the unresolved velocity issue:
   mismatch with GADGET,
   spurious bulk flow,
   and the fact that this remains the main limitation of the current implementation.

3. Make the `ν` issue more candid:
   if it changes collapse behaviour materially, then its role is not merely a harmless numerical convenience.

4. Keep the future-work section strong, but sharpen the classification of its ideas and clarify the gravitational coupling / operator-structure assumptions where needed.

5. Clean up the copy-editing and pronoun consistency.

## Bottom Line

Chapter 6 is a good conclusion and a strong future-work chapter. The reconciled verdict is that its main weakness is not structure but tone: the summary sometimes states the thesis's strongest claims more confidently than the earlier evidence supports. The future-work section is stronger, more precise, and in some ways more persuasive than the retrospective summary.
