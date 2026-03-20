# Chapter 2: Review of Numerical Simulations — Checklist

## Copy-editing (Shared #7)
- [x] p38: "under go" → "undergo"
- [x] p41: "computational practise" → "computational practice"
- [x] p43: "Direction summation" → "Direct summation"

## Copy-editing (Claude-only)
- [x] p37: Millennium Simulation notation — already clear as `$N = 1024^3$`
- [x] p40: "(Where N is the number of bodies)" — lowercase "where"

## Section 2.2 — Major rework (Shared #1–6)
- [x] Symplectic discussion too compressed/abstract — rewritten paragraph by paragraph
- [x] Symplectic/Hamiltonian wording overstated or inaccurate — "synonymous" removed, relationship clarified
- [x] Simply-connected claim is wrong — removed, replaced with continuous Hamiltonian flow description
- [x] Poincaré/symplectic group paragraph not defensible — rewritten around Noether's theorem properly
- [x] Refocus on practical reasons leapfrog/symplectic schemes matter — shadow Hamiltonian, bounded energy error
- [x] p46: "for real quantities this is always the case" (non-degeneracy) — replaced with correct null-direction explanation
- [x] p46: `dω = 0` should be called "closed", not "divergenceless" — fixed

## Conceptual (GPT-only)
- [x] "The word symplectic is synonymous with the word Hamiltonian" — corrected in Section 2.2 rewrite
- [x] Phase-space volume conservation ≠ Hamiltonian conservation — corrected, shadow Hamiltonian explained
- [x] p40: Sundman's regularisation / zero angular momentum claim — corrected to "zero total angular momentum for triple-collision case"
- [x] p37: GADGET-2 description blurs tree, PM, P3M, SPH — rewritten as TreePM + SPH
- [x] p42: "softening parameter is non-derivable" — rephrased to "phenomenological"
- [x] p44: `ρ = ρ̄ δ` — corrected, now defines δ = (ρ−ρ̄)/ρ̄ explicitly as overdensity
- [x] p45: "dominant calculation is still the PM force routine" — qualified ("typically dominates... because it operates on all particles via the FFT")
- [x] p48: "most perturbation modes will evolve linearly" — restricted to "at early times and on sufficiently large scales"
- [x] p50: "evolving the Euler equation" — corrected to "continuity, Euler, and Poisson equations"

## Conceptual (Claude-only)
- [x] p43: leapfrog accuracy — clarified as "second-order accurate method (global error O(h²))", symplectic description updated
- [x] p42 / Eq 2.2: `|ε| ≪ 1` — changed to `ε ≪ 1` described as "small positive length"
- [x] p44 / Eq 2.3: stated explicitly that comoving Poisson is sourced by overdensity ρ̄δ
- [ ] p44 / Eq 2.4: cross-check discrete Poisson kernel with Chapter 5 — deferred to Ch5 review
- [x] p48: linearisation — explicitly stated that `v · ∇v` has been dropped
- [x] p50: closing comparison of N-body / phase-space / fluid — reviewed, acceptable as prose; table format optional

## Superparticle / mass-resolution (Shared #8)
- [x] Reviewed — some repetition is acceptable in a review chapter for clarity

## Metadata
- [ ] PDF metadata fields blank
