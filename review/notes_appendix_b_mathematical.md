# Appendix B: Mathematical Appendix - Reconciled Review Notes

## Source Coverage

- GPT review: `gpt-review/notes_appendix_b_mathematical.md`
- Claude review: `claude-review/notes_appendix_b_mathematical.md`
- Overall result: both reviews agree the appendix is useful, but they differ on whether the FPA velocity derivation is fully correct as printed

## Shared Findings

The two reviews independently agree on the following main points.

1. **This is a useful appendix and should be kept.** Both reviews treat Appendix B as a worthwhile mathematical support appendix for the main thesis.

2. **Appendix B is substantially stronger than Appendix A.** Both reviews see it as directly useful for the thesis, rather than just historically interesting.

3. **The three-part structure works well.** Both reviews recognise that:
   - B.1 supports Chapter 3,
   - B.2 supports Chapter 4,
   - B.3 supports Chapter 5.

4. **The opening sentence is garbled and should be rewritten.** Both reviews independently flag the wording around the "equivalence of the Schrödinger equation and the continuity equation" as bad English, and the GPT review also correctly notes that it is conceptually incomplete.

5. **The appendix would benefit from better signposting.** Both reviews think it should be clearer what each section is proving and what notation is being inherited from earlier chapters.

6. **Section B.3 is useful but compressed.** Both reviews agree that the EdS derivation is valuable, but it assumes the reader remembers the Chapter 5 normalisations and would benefit from a brief reminder of what `\chi`, `U`, `\mathbf{y}`, and `\mathcal{L}` mean.

7. **The note about Widrow & Kaiser needs one more sentence.** Both reviews independently say that if the appendix mentions a mistake in Widrow & Kaiser, it should state explicitly what that mistake is.

8. **There are straightforward copy-editing issues.** Both reviews flag:
   the broken opening sentence,
   the "imaginary party" typo,
   and the abruptness of the final `QED`-style statement in B.2.

## Resolved Disagreements

There is one substantive disagreement between the source reviews, and the appendix itself resolves it.

1. **The initial-velocity derivation in B.2 is not internally sign-consistent as printed.** A direct check of the appendix pages supports the GPT review over the Claude review on this point:
   - Eq. B.19 defines the wavefunction with a negative phase, `ψ_i = (1+δ_i)^{1/2} e^{-iφ_i/ν}`;
   - with the standard probability current, that convention leads naturally to `v_i = -∇ φ_i`, not `+∇ φ_i`;
   - Eq. B.33 nevertheless states `v_i = ∇ φ_i`;
   - and Eqs. B.18 and B.34-B.35 do not remain sign-consistent with each other.

2. **The intermediate algebra in B.29-B.30 supports the concern.** The Claude review already noticed a sign flip there, but then accepted the final result. The reconciled verdict is that the derivation needs correction, not just clarification.

## GPT-Only Additions Worth Keeping

These points appear in the GPT review but not the Claude review and are worth preserving.

1. **Section B.1 is conceptually overstated as written.** The appendix does not show equivalence to the continuity equation alone; it is the full Madelung pair that matters.

2. **`p230` / Eqs. B.8-B.10**: Reusing `V` for a Bernoulli-type quantity inside a derivation that already uses `V` for the external potential is confusing and should be avoided.

3. **`p230-231`**: Since the derivation explicitly omits pressure and potential terms in order to isolate continuity, the appendix should say more clearly that it is only extracting part of the Madelung system at that stage.

4. **`p235`**: The notation around the dimensionless parameter `\mathcal{L}` is visually awkward because the same letter is already being used for the physical box size. Even if mathematically correct, it deserves clearer explanation.

5. **Metadata**: the PDF metadata fields are blank here as well.

## Claude-Only Additions Worth Keeping

These points appear in the Claude review but not the GPT review and are worth preserving.

1. **Section B.1 is mathematically valuable once the exposition is cleaned up.** The Claude review's spot-checks are useful in establishing that most of the Madelung derivation itself is structurally sound.

2. **The `α_1` ... `α_4` substitutions in B.2 are serviceable, even if clunky.** The GPT review is right that they are cumbersome, but the Claude review is right that they do provide a full bookkeeping path through the probability-current algebra.

3. **Section B.3 should connect itself back to Chapter 5 more explicitly.** The Claude review is right that the appendix currently ends too abruptly after Eq. B.46.

4. **A short appendix-level introduction would help.** The Claude review is right that one sentence stating what B.1, B.2, and B.3 each support would improve usability.

## Reconciled Priorities

1. Fix the real mathematical/presentation problem first:
   make the B.2 phase convention, probability-current calculation, and final velocity formula mutually consistent.

2. Rewrite the opening of B.1 so it correctly states what the Madelung transform yields and so the English is clean.

3. Separate notation more carefully:
   especially `V` as external potential versus Bernoulli-like quantity,
   and `\mathcal{L}` versus the physical box length `L`.

4. Add one or two orienting sentences:
   what each appendix section is for,
   and what the Widrow & Kaiser mistake actually was.

5. Keep the appendix concise, but make the prose more explicit about what is being proved and what is being imported from earlier chapters.

## Bottom Line

Appendix B is mathematically useful and is worth keeping. The reconciled verdict is that most of it is serviceable, but the FPA initial-velocity derivation is not internally sign-consistent as printed and should be corrected. Once that is fixed, the remaining improvements are mostly about clarity, notation, and better integration with the main thesis chapters.
