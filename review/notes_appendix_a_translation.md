# Appendix A: Translation — Review Notes

## Overall Assessment
A short appendix (6 pages, pp223-228) containing an English translation of Madelung's 1927 paper "Quantentheorie in hydrodynamischer Form" (Zeitschrift für Physik, 40, 322), translated by T. Krämer and E. Thomson on 14th August 2008. This is a valuable scholarly contribution — Madelung's original paper is in German and not widely accessible in English. The translation is central to the thesis since the Madelung transform (Chapter 3) is the theoretical foundation of the entire wave-mechanical approach. The translation reads clearly and the mathematical content is well-preserved. There are a few minor issues.

---

## Translation Quality

1. **p223**: "abstract" is presented in lower case — should be capitalised to "**Abstract**" for consistency with standard English academic formatting. (Though this may be preserving the original German formatting where "Zusammenfassung" is sometimes not capitalised.)

2. **p224**: "ansätze" — correctly italicised as a German loanword. Good.

3. **p224**: Equation A.2: ∇²ψ - (8π²m/h²)Uψ - i(4πm/h)(∂ψ/∂t) = 0 — this is the time-dependent Schrödinger equation in Madelung's notation with h (not ℏ). In modern notation with ℏ = h/2π, this becomes ∇²ψ - (2m/ℏ²)Uψ - (2mi/ℏ)(∂ψ/∂t) = 0, which rearranges to iℏ(∂ψ/∂t) = -(ℏ²/2m)∇²ψ + Uψ. Correct.

4. **p224**: The substitution ψ = αe^{iβ} and the resulting equations A.3 and A.4 — these are the real and imaginary parts of the Schrödinger equation after the Madelung ansatz, directly corresponding to Eqs 3.4-3.8 in Chapter 3. The connection should be noted (and is, in the Chapter 3 text).

5. **p225**: Equation (3'): ∂φ/∂t + ½(∇φ)² + U/m - (∇²α/α)(h²/8π²m²) = 0 — this is the quantum Hamilton-Jacobi equation (Bernoulli equation with quantum pressure). The last term is the quantum potential Q = -(ℏ²/2m)(∇²α/α), written in Madelung's h-notation. Correct.

6. **p225**: Equation (3''): ∂u/∂t + ½∇u² = du/dt = -∇U/m + ∇(∇²α/α)(h²/8π²m²) — the Euler equation with quantum pressure force. The step from ∂u/∂t + ½∇u² to du/dt is correct only for irrotational flow (u = ∇φ implies u·∇u = ½∇u²). Good that Madelung notes ∇ × u = 0.

7. **p225**: The translator's footnotes (footnotes 1-5) are helpful, providing the original German terms (Besonderheit, Strömungszustände, statische Gebilde, Strömung). These are a nice touch for readers who might want to consult the original.

8. **p226**: Equation A.5 and (5'): The energy expression W = (m/2)∇φ² + U - (∇²α/α)(h²/8π²m) — this is the quantum Hamilton-Jacobi energy. The integral form (5') gives the total energy as kinetic + potential + quantum pressure energy. Correct.

9. **p227**: Madelung's discussion of the many-electron problem is remarkably prescient. His three options (electrons fuse, exclude each other, or penetrate each other) anticipate later developments in density functional theory and many-body quantum mechanics. The translation of this passage reads naturally.

10. **p227**: "Schrödinger doesn't give a specific form" — the contraction "doesn't" feels too informal for a translation of a 1927 academic paper. Consider "does not."

---

## Comments on the Appendix's Role in the Thesis

11. This translation directly supports Chapter 3's derivation of the Madelung/hydrodynamic equations. Having the original source available in English strengthens the thesis's scholarly foundation.

12. The translation preserves Madelung's original equation numbering (A.1, A.2, ..., with primed numbers 3', 4', 5' for derived equations). This is helpful for cross-referencing with the original paper.

13. It would be useful to include a brief translator's note at the beginning explaining any choices made in the translation (e.g., how mathematical notation was modernised, whether any passages were particularly difficult to translate, how the primed equation numbers relate to Madelung's originals).

14. The bibliography for this appendix contains only one reference (Schrödinger's papers). This is correct — Madelung's original paper cited very little.

---

## Equations Check

- **Eq A.1**: ∇²ψ₀ + (8π²m/h²)(W-U)ψ₀ = 0 — correct (time-independent Schrödinger equation in Madelung's notation)
- **Eq A.2**: ∇²ψ - (8π²m/h²)Uψ - i(4πm/h)(∂ψ/∂t) = 0 — correct (time-dependent Schrödinger equation)
- **Eq A.3**: ∇²α - α(∇β)² - (8π²m/h²)Uα - (4πm/h)α(∂β/∂t) = 0 — correct (real part)
- **Eq A.4**: α∇²β + 2(∇α·∇β) - (4πm/h)(∂α/∂t) = 0 — correct (imaginary part, continuity equation)
- **Eq (4')**: ∇·(α²∇φ) + ∂α²/∂t = 0 — correct (continuity equation with φ = -2βh/2πm)
- **Eq (3')**: ∂φ/∂t + ½(∇φ)² + U/m - (∇²α/α)(h²/8π²m²) = 0 — correct (quantum Bernoulli/Hamilton-Jacobi)
- **Eq (3'')**: Euler equation with quantum pressure — correct
- **Eq A.5**: W = (m/2)∇φ² + U - (∇²α/α)(h²/8π²m) — correct (energy expression)
- **Eq (5')**: Integral energy form — correct
