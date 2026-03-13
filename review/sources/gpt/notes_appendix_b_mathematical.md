# Appendix B: Mathematical Appendix - Review Notes

## Overall Assessment

This appendix is useful and worth keeping. It gives supporting derivations for three important parts of the thesis:

- the Madelung transform,
- the initial velocity used in the FPA,
- and the Schrödinger/Poisson equations in the Einstein-de Sitter model.

Of the two appendices, this one is the more directly useful mathematically. The main issue is not that the derivations are worthless, but that a few parts are presented too loosely. There is also one definite sign inconsistency in the FPA initial-velocity derivation that should be fixed.

---

## Higher-Priority Issues

1. **p229 / opening sentence**: The appendix says it shows the "mathematical equivalence of the Schrödinger equation and the Continuity equation of fluid dynamics". That is not the right statement. The Madelung transform gives a fluid-like system consisting of:
   - a continuity equation, and
   - a Bernoulli / Hamilton-Jacobi / Euler-type equation with the quantum-pressure term.
   The Schrödinger equation is not equivalent to the continuity equation alone.

2. **p230 / Eqs. B.8-B.10**: The notation becomes confusing because `V` is reused. At the top of the appendix, `V` is the potential in the Schrödinger equation, but in Eq. B.8 it is being used for a Bernoulli-type expression built from `\phi`. Those should not share the same symbol inside one derivation.

3. **p230-231 / Madelung section**: The derivation explicitly says "The pressure and potential terms have been omitted here" before arriving at the continuity equation. That is acceptable if the goal is only to extract continuity, but it means this section is not really showing the full Madelung equivalence as advertised in the opening sentence. The appendix should say that more clearly.

4. **p232-234 / Section B.2**: There is a definite sign inconsistency in the initial-velocity derivation.
   - Eq. B.19 defines the initial wavefunction with a **negative** phase:
     `\psi_i = (1+\delta_i)^{1/2} e^{-i\phi_i/\nu}`.
   - With the standard probability-current formula, that convention gives
     `v_i = -\nabla \phi_i`,
     not `v_i = \nabla \phi_i`.
   - However the appendix concludes in Eq. B.33 that `v_i = \nabla \phi_i`.
   This should be corrected or the phase convention changed.

5. **p232 vs p234 / Eqs. B.18 and B.34-B.35**: The sign inconsistency shows up a second time.
   - Eq. B.18 defines
     `\delta_i = -\delta_a \cos(2\pi x/p)` and `\phi_i = -(p/2\pi)^2 \delta_i`,
     so `\phi_i` is proportional to `+\cos(2\pi x/p)`.
   - But Eq. B.34 writes the velocity as the gradient of
     `-(p/2\pi)^2 \delta_a \cos(2\pi x/p)`,
     which corresponds to the opposite sign choice.
   So B.18 and B.34 are not internally consistent as printed.

6. **p235 / Eq. B.46**: The final EdS Schrödinger equation is useful, but the definition of the dimensionless parameter `\mathcal{L}` is very easy to misread because the same letter `L` is also used for the physical box length earlier in the derivation. Even if this is mathematically correct, the notation is visually awkward and should be explained more clearly.

---

## Clarity / Presentation Issues

7. **p229**: "The role of the pressure term is discussed elsewhere and is extracted out from this derivation in the appropriate place." This is awkwardly phrased. Better to say explicitly whether the appendix is:
   - deriving only the continuity equation,
   - or deriving the full Madelung pair and then isolating the quantum-pressure term.

8. **p232**: "The subscript i refers to the initial states and does not the imaginary party of a complex number." This is a clear typo and should be something like:
   "does not denote the imaginary part of a complex number."

9. **p232**: "Given both the density potential velocity field" is ungrammatical. This should be something like:
   "Given both the density field and the velocity potential..."

10. **p233**: The `\alpha_1, \alpha_2, \alpha_2^*, \alpha_3, \alpha_4` bookkeeping works, but it is more cumbersome than illuminating. The derivation is correct in spirit, but the notation obscures the main point.

11. **p233**: "Essentially, `\alpha_4 = i/\nu \nabla \phi_i`" is exactly where the sign issue enters. With the stated phase convention and Eq. B.18, this should be checked carefully.

12. **p234-235 / Section B.3**: This section is concise and useful, but it assumes the reader already remembers all the normalisations from Chapter 5. A one-line reminder of what `\chi`, `U`, `\mathbf{y}`, and `\mathcal{L}` represent would help.

13. **Widrow & Kaiser correction**: The appendix says there is a mistake in their version of the coupled equations, but it does not briefly restate what the mistake is. A short explicit note would improve readability.

---

## Minor Copy-Editing Issues

14. **p229**: "Continuity equation of fluid dynamics Schrödinger equation" is garbled and needs rewriting regardless of the conceptual issue noted above.

15. **p231**: "Re-call" should be `Recall`.

16. **p234**: "QED. This matches the result of Short in his PhD thesis." is fine, but slightly abrupt; a short clause noting that this is the expected FPA initial velocity would read more smoothly.

17. **Throughout**: The appendix would benefit from more explicit signposting between "what is being proved" and "what is being assumed from earlier chapters."

18. **Metadata**: As with the other PDFs, the metadata fields are blank (`Title`, `Author`, etc.). Minor, but still worth fixing.

---

## Bottom Line

Appendix B is mathematically useful and substantially stronger than Appendix A, but it needs one real mathematical fix and a small amount of cleanup.

The highest-priority corrections are:

- fix the opening claim about what the Madelung transform actually proves,
- stop reusing `V` for two different meanings in the same derivation,
- correct the sign inconsistency in the initial-velocity derivation,
- and make Eqs. B.18, B.19, B.33, and B.34-B.35 mutually consistent.

Once those are corrected, the appendix should work well as a supporting derivation appendix for the main thesis.
