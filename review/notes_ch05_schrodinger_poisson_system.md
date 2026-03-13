# Chapter 5: Solving the Full Schrödinger-Poisson System - Reconciled Review Notes

## Source Coverage

- GPT review: `gpt-review/notes_ch05_schrodinger_poisson_system.md`
- Claude review: `claude-review/notes_ch05_schrodinger_poisson_system.md`
- Overall result: extremely strong agreement on both the chapter's originality and its main unresolved defects

## Shared Findings

The two reviews independently agree on the following main points.

1. **This is the core technical chapter of the thesis and contains several genuine contributions.** Both reviews treat it as the central original work, especially for:
   periodic boundaries,
   the full 3D implementation,
   the cosmological generalisation,
   and the detailed numerical methodology.

2. **The periodic-boundary treatment is one of the strongest parts of the chapter.** Both reviews independently praise the double-recursion boundary method and treat it as one of the most convincing original contributions.

3. **The Widrow-Kaiser analysis/correction is a real contribution.** Both reviews agree that the comoving-density clarification and the corrected equations are important, even if that section could be shortened.

4. **The Poisson-equation presentation is a major problem area.** Both reviews independently flag:
   the missing mean-density subtraction for a periodic box,
   the misleading standalone integral expression for the potential,
   and the incorrect use of `ψ_k ψ_k^*` instead of the Fourier transform of `|ψ|^2`.

5. **Equation 5.64 likely has the wrong Gaussian normalisation.** Both reviews independently flag the 1D Gaussian packet normalisation as incorrect.

6. **The momentum/conservation discussion is overstated or wrong as written.** Both reviews treat the conservation claims as too strong, and both flag Eq. 5.69 as not looking like a correct total-momentum expression.

7. **The velocity discrepancy is one of the most serious unresolved issues in the thesis.** Both reviews agree that this is not a minor caveat. It is a major limitation because it undermines confidence in the physical interpretation of the cosmological velocity field.

8. **The spurious growing bulk flow is a major unresolved defect.** Both reviews independently treat the non-zero mean velocities as an important problem, not just a future-work curiosity.

9. **The chapter is commendably honest about its limitations.** Both reviews strongly value the candour with which the chapter discusses velocity problems, interference, and unsatisfactory runs.

10. **The chapter is too long in places and could use tighter organisation.** Both reviews think the overall structure works, but the chapter would benefit from clearer signalling of what is new and a more explicit separation between solid results and unresolved issues.

11. **There are many straightforward copy-editing issues.** Both reviews independently flag repeated articles, Widrow/Widow typos, awkward phrasing, and colloquial wording.

12. **Pronoun inconsistency remains a problem.** Both reviews note the continued shift between `I` and `we`.

## GPT-Only Additions Worth Keeping

These points appear in the GPT review but not the Claude review and are worth preserving.

1. **`p113-114`**: The basis discussion around the matrix structure appears backwards. In position space the potential should be diagonal, while the finite-difference kinetic operator is banded. If another basis is intended, that must be stated.

2. **`p118-120`**: The chapter needs to distinguish more clearly between a truly implicit scheme and the practical semi-implicit/quasi-implicit implementation that uses the current-time potential instead of the advanced one.

3. **`p149-150`**: The top-hat discussion appears internally inconsistent about the role of `ν`, and this should be reconciled more carefully.

4. **`p133`**: The switch from Hydra in Chapter 4 to GADGET-2 here should be motivated more explicitly.

5. **`p153-156`**: The discussion of poor lower-resolution runs should separate more clearly what looks like a resolution issue from what might be a deeper framework or implementation issue.

6. **`p155-164`**: Claims that certain density behaviours are "inevitable" or "fundamental" to wave-mechanics are too strong unless demonstrated more directly.

7. **`p158-166`**: The smoothing discussion should state more explicitly whether Gaussian smoothing is just a diagnostic aid or part of the meaningful physical comparison.

8. **`p174`**: The final no-vorticity claim is interesting, but the threshold/diagnostic logic is still qualitative and should be framed more cautiously.

9. **Metadata**: the PDF metadata fields are blank here as well.

## Claude-Only Additions Worth Keeping

These points appear in the Claude review but not the GPT review and are worth preserving.

1. **`p113`**: The attribution around the "Cayley method" is unclear because the text seems to attribute it to Widrow & Kaiser in one place and Goldberg elsewhere.

2. **`p118`**: The unusual convention `m = 1/2` in the derivation should be justified briefly.

3. **`p119`**: Re-using `δ` for the timestep is confusing given its use elsewhere for density contrast and delta functions.

4. **`p115`**: The claim that the Cayley factors commute should be justified by stating that both factors are functions of the same matrix.

5. **`p116-117`**: The move from the naive higher-dimensional Goldberg extension to Suzuki/Strang splitting should be highlighted even more strongly as a key methodological step.

6. **`p116`**: The ad hoc split of `V/3` across the three dimensions should either be motivated more clearly or more explicitly treated as a failed intermediate idea.

7. **`p119-129`**: The recursion relations are central but dense; a worked example or sweep-direction illustration would help implementation-minded readers.

8. **`p123`**: The practical observation that two boundary-recursion passes suffice deserves either a convergence comment or a brief theoretical explanation.

9. **`p126-127`**: The dimensionless parameter `\mathcal{L}` is crucial but under-explained physically.

10. **`p132-133`**: The reason for choosing GADGET-2 rather than Hydra here could be stated more directly.

11. **`p139-141`**: The interference/dispersion timescale discussion could be made more quantitative.

12. **`p168`**: The GADGET velocity scaling should be double-checked as a possible source of the mismatch.

13. **`p169`**: A further possible explanation for the initial velocity mismatch is that the wavefunction phase was built from the gravitational potential rather than directly from the velocity potential.

14. **`p174-177`**: A dedicated "Known Issues and Limitations" subsection would improve the summary materially.

## Reconciled Priorities

1. Fix the Poisson-equation presentation first:
   add the mean-subtraction clarification for periodic boxes,
   remove the misleading integral formula,
   and replace `ψ_k ψ_k^*` with the Fourier transform of `|ψ|^2`.

2. Correct the conservation discussion:
   especially Eq. 5.69,
   and avoid implying that mass plus momentum conservation establishes energy conservation.

3. Reframe the velocity results much more explicitly as a major unresolved limitation:
   the scaling issue,
   the mismatch with GADGET,
   and the spurious net bulk flow.

4. Tighten the presentation of the numerical method:
   implicit vs semi-implicit language,
   the Goldberg-to-Strang transition,
   and the physical meaning of the key dimensionless parameters.

5. Keep highlighting the real achievements:
   periodic boundaries,
   cosmological generalisation,
   and the careful mass-conservation tests.

6. Add a short explicit limitations subsection near the end so the chapter's strongest and weakest results are cleanly separated.

## Bottom Line

Chapter 5 is the most important and most original chapter in the thesis, and both reviews treat it with corresponding seriousness. The reconciled verdict is that the chapter contains substantial real contributions, but it also carries several non-trivial equation/presentation errors and one major scientific limitation: the cosmological velocity field is not yet reliable. The density-side results and numerical implementation look promising; the velocity-side defect needs to be treated as central whenever the chapter is summarised.
