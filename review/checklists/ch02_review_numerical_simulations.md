# Chapter 2: Review of Numerical Simulations — Checklist

## Copy-editing (Shared #7)
- [ ] p38: "under go" → "undergo"
- [ ] p41: "computational practise" → "computational practice"
- [ ] p43: "Direction summation" → "Direct summation"

## Copy-editing (Claude-only)
- [ ] p37: Millennium Simulation notation — clarify `1024^3`
- [ ] p40: "(Where N is the number of bodies)" — lowercase "Where" mid-sentence

## Section 2.2 — Major rework (Shared #1–6)
- [ ] Symplectic discussion too compressed/abstract — rewrite or substantially cut back
- [ ] Symplectic/Hamiltonian wording overstated or inaccurate — correct
- [ ] Simply-connected claim is wrong — remove or fix
- [ ] Poincaré/symplectic group paragraph not defensible — major repair or remove
- [ ] Refocus on practical reasons leapfrog/symplectic schemes matter in N-body work
- [ ] p46: "for real quantities this is always the case" (non-degeneracy) — replace
- [ ] p46: `dω = 0` should be called "closed", not "divergenceless"

## Conceptual (GPT-only)
- [ ] "The word symplectic is synonymous with the word Hamiltonian" — correct
- [ ] Phase-space volume conservation ≠ Hamiltonian conservation — rewrite or remove
- [ ] p40: Sundman's regularisation / zero angular momentum claim is misleading
- [ ] p37: GADGET-2 description blurs tree, PM, P3M, SPH — clean up
- [ ] p42: "softening parameter is non-derivable" — rephrase ("ad hoc" or "phenomenological")
- [ ] p44: `ρ = ρ̄ δ` — check notation against thesis conventions
- [ ] p45: "dominant calculation is still the PM force routine" — qualify or source
- [ ] p48: "most perturbation modes will evolve linearly" — restrict to early times / large scales
- [ ] p50: "evolving the Euler equation" — also includes continuity and Poisson

## Conceptual (Claude-only)
- [ ] p43: leapfrog accuracy — clarify as second-order accurate globally
- [ ] p42 / Eq 2.2: `|ε| ≪ 1` dimensionally awkward; absolute value bars unnecessary
- [ ] p44 / Eq 2.3: state explicitly that comoving Poisson is sourced by overdensity
- [ ] p44 / Eq 2.4: cross-check discrete Poisson kernel with Chapter 5
- [ ] p48: linearisation — explicitly state that `v · ∇v` has been dropped
- [ ] p50: closing comparison of N-body / phase-space / fluid — consider table format

## Superparticle / mass-resolution (Shared #8)
- [ ] Tighten wording; superparticle explanation somewhat repetitive

## Metadata
- [ ] PDF metadata fields blank
