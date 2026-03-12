# Chapter 7: Epilogue — Vorticity and Spin - Review Notes

## Overall Assessment
A speculative and ambitious chapter (22 pages of content, pp191-213) that explores two novel ideas: (1) incorporating vorticity via gravitoelectromagnetism (GEM), and (2) modelling spinning extended objects via a Pauli-like equation derived from a multipole expansion. The chapter is explicitly labelled as "not part of the work of this thesis" and is more exploratory in nature. The GEM section (7.1) is well-grounded in existing physics, while the spin/Pauli section (7.2) is genuinely original and intellectually impressive — deriving a classical gravitational analogue of the Pauli equation from first principles using Clifford algebra. However, the speculative nature leads to some loose reasoning, and the mathematical presentation in the later sections becomes dense and would benefit from more careful exposition.

---

## Typos and Grammatical Errors

1. **p191**: "ψ = (0,0), ψ ∈ ℂ" — as noted in Chapter 4 review, writing ψ = (0,0) for ψ = 0 is unusual. Simply "ψ = 0" is clearer. This notation is repeated from section 3.5.

2. **p191**: "centrepoint" — not standard; use "centre point" or simply "centre."

3. **p192**: "focussed" — British English accepts both "focused" and "focussed," but the thesis uses "focused" elsewhere. Should be consistent.

4. **p193**: "We could we propose a phenomenological fix" — extra "we": should be "We could propose a phenomenological fix" or "We propose a phenomenological fix."

5. **p194**: "finding the centre of the vortex could be a possible method of identifying spiral galaxies" — "centre" is spelled correctly (British English), but note the earlier use of "centrepoint" (p191). Minor inconsistency.

6. **p194**: "the wavefunction would take a form ψ = (a,b) where a ≠ 0, b = 0; a,b ∈ ℝ; ψ ∈ ℂ" — this is confusing. If ψ ∈ ℂ then ψ = a + ib = a (a real number). Saying the velocity is "undefined" when ψ is real and non-zero is incorrect — v = (ℏ/m)ℑ(∇ψ/ψ), and if ψ is real, ∇ψ is real, so ℑ(∇ψ/ψ) = 0, giving v = 0 (not undefined). The velocity is undefined only when ψ = 0. This needs correction.

7. **p195**: "handing the vector potential" → "**handling** the vector potential" (typo)

8. **p196**: "straight forward" → "**straightforward**" (one word; recurring issue)

9. **p197**: "In the previous section we mentioned that simulations of Large Scale Structure already display regions under going rotational motion" → "**undergoing**" (one word)

10. **p198**: "I will re-iterate by we 'omitted' a large divisor" → "I will re-iterate **that** we 'omitted' a large divisor" (missing word). Also "re-iterate" → "reiterate" (one word).

11. **p199**: "we don't expect" — contraction in formal thesis. "we do not expect."

12. **p199**: "The result is likely to be independent of the simulation method used" — this is an assertion without justification. The GEM field calculation depends on the particle data; whether the result is "independent of method" depends on whether the methods agree on the velocity field.

13. **p200**: "description)." followed by "Hence, with this idea I propose a 'minimal' extension" — the paragraph break is awkward. These should be connected or the transition improved.

14. **p205**: "|r| << 1" — this should have units or be relative to some scale. "|r| ≪ R" would make more physical sense (the internal separation is much smaller than the distance to the centre of mass).

15. **p207**: "the force is no longer the simple form that is has from Newton's theory" → "**it** has" (typo: "is" → "it")

16. **p208**: "the the effect" — if this appears as a double "the", it's a typo.

17. **p209**: "the quantity d_i is a scalar" — d_i is a component of the vector **d**, so calling it "a scalar" is slightly misleading. Better: "d_i is the i-th component of **d**."

---

## Clarity and Content Issues

### Introduction (pp191-192)

18. **p191**: The opening statement "This chapter is not part of the work of this thesis" is honest but slightly undermining. Consider: "This chapter presents speculative extensions that go beyond the main results of this thesis" — this frames it more positively while still managing expectations.

19. **p191**: "∇ × F_g = ∇ × ∇φ_g ≡ 0" — correct. The curl of a gradient is identically zero. Good clear statement.

20. **p191-192**: The distinction between vorticity (∇ × v ≠ 0) and the circulation theorem (∇ × F = 0) is important. The text correctly notes that Kelvin's theorem holds before shell crossing but not after. However, the text could be clearer that in wave-mechanics, the velocity field v = (ℏ/m)ℑ(∇ψ/ψ) is irrotational everywhere *except* at ψ = 0 nodes — this is a topological constraint, not a dynamical one. The vorticity in wave-mechanics is quantised and confined to vortex lines, unlike classical fluid vorticity.

21. **p192**: "we expect the vorticity to be zero in the Universe now as it is related to the decaying mode of density perturbations" — this is correct for the linear regime, but tidal torque theory (Peebles 1969, White 1984) shows that angular momentum is generated at second order in perturbation theory. The initial vorticity is zero, but angular momentum develops naturally. This well-established result should be mentioned.

### Section 7.1 — Vorticity

22. **p193**: The proposal of a gravitational vector potential V_c to allow ∇ × F_g ≠ 0 is interesting but the transition from "phenomenological fix" to GEM is somewhat abrupt. The text should be clearer: GEM provides a *physical* (not merely phenomenological) vector potential that arises naturally from general relativity in the weak-field limit.

23. **p193-194**: The discussion of identifying vortex centres with galaxy locations is creative but highly speculative. The connection between wavefunction nodes (ψ = 0) and galaxy positions is not straightforward — in the wave-mechanics approach, galaxies form at density *peaks* (|ψ|² maxima), not at density zeros.

24. **p194**: The statement about the Halo Model being "only appropriate to particle simulations" is debatable. The halo model is fundamentally a statistical description of the density field and could in principle be applied to any density field, including |ψ|².

### Section 7.1.1 — Gravitoelectromagnetism

25. **p195**: The GEM equations table is well-presented and the analogy with Maxwell's equations is clear. One important caveat missing: the GEM equations are only valid in the weak-field, slow-motion limit of GR. For cosmological simulations where v ≪ c everywhere, this is well-satisfied, but it should be stated.

26. **p195**: "ρ_g is simply the mass density as it appeared earlier in this thesis" — but in the GEM context, ρ_g is the source for the gravitoelectric field (Newtonian gravity), while J_g = ρ_g v is the source for the gravitomagnetic field. The text should note that J_g is the mass current density, directly analogous to the probability current J = (ℏ/2mi)(ψ*∇ψ - ψ∇ψ*) already computed in the thesis.

27. **p196**: Equation 7.1: F_g = m(E_g + v × 2B_g) — the factor of 2 in front of B_g is specific to the GEM formulation (it differs from the EM Lorentz force by this factor due to the spin-2 nature of gravity vs. spin-1 of electromagnetism). This factor should be noted and explained, as it's easy to miss.

28. **p196**: Equation 7.2: iℏ(∂/∂t)ψ = -(ℏ²/2m)(∇ - ie/ℏ A)²ψ — correct (minimal coupling). For the GEM case, the replacement would be e → m (mass plays the role of charge) and A → A_g (gravitomagnetic vector potential). But there's a subtlety: in EM, the coupling is e/ℏ, while in GEM it would be m/ℏ = 1/ν. This connection to the existing ν parameter should be highlighted.

29. **p196**: "the magnetic field just changes the phase of the wavefunction, so it is very easy to compute" — this is only true if A is spatially uniform. If A varies in space (as it would in any realistic gravitomagnetic field), the (∇ - ieA/ℏ)² operator couples the vector potential to the kinetic energy non-trivially. This oversimplification should be corrected.

### Section 7.1.1.1 — Test of Gravitomagnetism

30. **p197**: Equation 7.3: B_g = (G/2c²)(L - 3(L·r/r)r/r)/r³ — this is the gravitomagnetic field of a spinning mass (dipole approximation). However, L here is the angular momentum of *the system*, not of individual particles. The text says "L is the angular momentum of a particle in the simulation" but then computes it as a global quantity from all particles. The distinction between the field *produced by* a spinning mass and the field *experienced by* test particles should be clearer.

31. **p197**: Equation 7.4: L = ∫ρ(v × r)dV = Σm(v × r) — the angular momentum should be L = ∫ρ(r × v)dV, not (v × r). The cross product is anti-commutative, so r × v = -(v × r). This sign error propagates to L having the wrong sign. Standard convention: L = r × p = m(r × v).

32. **p198**: The table of angular momentum and gravitomagnetic field values is interesting. The key observation — that max(F_{Bg}/F_{Eg}) decreases over time — is noted but not fully explained. In an expanding universe, peculiar velocities grow but distances also grow; the ratio v/c decreases in the linear regime, so B_g/E_g should indeed decrease. Worth stating this explicitly.

33. **p198**: "the ratios in the above table may seem higher than expected but I will re-iterate by we 'omitted' a large divisor of the F_{Bg} calculation: 1/c² = 1" — setting c = 1 in a non-relativistic simulation is confusing. If the simulation uses physical units (as it should for comparing with Hydra), then c has a definite value and the gravitomagnetic force is suppressed by 1/c². The statement that 1/c² = 1 suggests natural units, but the Hydra simulation presumably uses physical units. This needs clarification.

34. **p198-199**: The conclusion that "we don't expect that the field will be significant in a cosmological simulation but rather in the simulation of a system that has stronger gravity" is reasonable. But the more important point is that GEM effects are O(v²/c²) corrections to Newtonian gravity — for cosmological peculiar velocities of ~10³ km/s, this gives corrections of order (10³/3×10⁵)² ~ 10⁻⁵, which is indeed negligible.

### Section 7.2 — Spinning Objects

35. **p199**: "The Pauli-like equation that we eventually derive for a wave-mechanics implementation of spinning objects is entirely new, as far as we are aware." — This is a strong claim. The derivation is indeed creative, but the reader should be warned that "new" here means "speculative and untested."

36. **p199-200**: The key conceptual distinction — spin (intrinsic rotation of extended objects) vs. vorticity (circulation in the velocity field) — is well-articulated. The point that spin requires only a differential gravitational field (tidal forces), not a vector potential, is important and correct.

37. **p200**: "We believe the total spin of the Universe should be conserved" — this is correct by angular momentum conservation, but the claim that "naively we expect there to be a measurable difference in the results of a simulation of a universe with zero net-spin from that with non-zero net-spin" needs qualification. If total spin is zero, local spin can still develop through tidal interactions — it's the *distribution* of spin that changes, not the total.

38. **p200-201**: The approach of adapting Jackson's multipole expansion for electrostatics to gravitational systems is well-motivated and clearly attributed.

### Section 7.2 — Multipole Expansion (pp201-206)

39. **p201**: Equation 7.5: U(R+r) = φ₀∫ρ dV + ∂ᵢφ∫ρrᵢ dV + ∂ᵢⱼφ∫ρrᵢrⱼ dV + ... — correct (Taylor expansion of potential). The notation is clear.

40. **p202**: The argument that the gravitational dipole moment is always zero (because mass has only one sign) is correct and well-explained. The contrast with electromagnetism (where opposite charges create non-zero dipoles) is a nice pedagogical touch.

41. **p202-203**: Equations 7.6-7.9: The dipole and quadrupole calculations for the two-mass "dumbbell" system are correct and clearly presented.

42. **p203**: Equations 7.10-7.12: The dipole potential U^(d) = ∂ᵢφ · dᵢ is correct but since the dipole moment is zero for gravity, this term vanishes. The text acknowledges this but then uses the *difference vector* d = r^(1) - r^(2) as a proxy — this is a clever trick but should be more clearly distinguished from the actual dipole moment.

43. **p204**: Equations 7.13-7.16: The quadrupole potential U^(q) = ∂ᵢⱼφ · Qᵢⱼ with the full 3×3 expansion — correct. The notation is clear and the expansion is systematic.

44. **p204**: Equation 7.17: F(R+r) = ∇⊗U(R+r) = ∂ᵢ(φ₀)∫ρ dV + ∂ᵢ(∂ⱼφ)∫ρrⱼ dV + ... — correct (force as gradient of potential). The notation ∇⊗U using the tensor product is non-standard; ∇U or simply F = -∂U/∂R would be more conventional.

45. **p204-205**: Equations 7.18-7.20: The quadrupole contribution to the force F^(q) = ∂ᵢ(∂ⱼₖφ)∫ρrⱼrₖ dV involves the third derivative of the potential. This is correct but computationally expensive — it requires computing ∂³φ/∂xᵢ∂xⱼ∂xₖ, a rank-3 tensor with up to 10 independent components. This practical cost should be noted.

46. **p205**: Equations 7.21-7.24: The force on each sub-particle and the addition/subtraction to get centre-of-mass motion and relative (spin) motion — this is a clear and elegant decomposition. The key result: the centre-of-mass feels only the monopole force (Eq 7.23 reduces to dipole = 0), while the *difference* between the sub-particle forces (Eq 7.24) drives spin.

47. **p205**: Equation 7.25: 2m(r_j^(1) - r_j^(2)) = 2m[(x,0,0) - (-x,0,0)] = 4m(x,0,0) — correct for the specific dumbbell geometry.

### Section 7.2.2 — Pauli Equation (pp207-210)

48. **p208**: Equation 7.27: U = ∂ᵢ(φ)∫ρrᵢ dV ≃ B·S = B_xσ_x + B_yσ_y + B_zσ_z — the identification of the dipole potential with a "B·S" spin-coupling term is the central creative leap of the chapter. The analogy is: B ↔ ∇φ (gradient of gravitational potential), S ↔ ∫ρd dV (mass-weighted displacement, mapped to Pauli matrices via Clifford algebra). This is intellectually bold. However, calling the gradient of the gravitational potential "B" is confusing because in the GEM section (7.1), B was the gravitomagnetic field. The notation should be distinct.

49. **p208**: The use of the Clifford algebra isomorphism Cl₃ ≃ Mat(2,ℂ) to map the displacement vector d = dᵢeᵢ to dᵢσᵢ is mathematically sound. The key property used is that the Pauli matrices form a basis for the traceless Hermitian 2×2 matrices, which represent ℝ³ vectors in the spinor formalism.

50. **p209**: Equations 7.32-7.34: The Pauli-like B·S term expanded as a 2×2 matrix:
   U = Σmⱼ (B_zd_z, B_xd_x - iB_yd_y; B_xd_x + iB_yd_y, -B_zd_z)
   — let me check this. Using σ_x = (0,1;1,0), σ_y = (0,-i;i,0), σ_z = (1,0;0,-1):
   B_xd_xσ_x + B_yd_yσ_y + B_zd_zσ_z = (B_zd_z, B_xd_x - iB_yd_y; B_xd_x + iB_yd_y, -B_zd_z)
   This is correct. The matrix is Hermitian as required.

51. **p209**: Equation 7.35: The Pauli matrices are listed as σ₁ = σ_x = (0,1;1,0), σ₂ = σ_y = (0,-i;i,0), σ₃ = σ_z = (1,0;0,-1) — correct, standard definitions.

52. **p209-210**: Equation 7.36: iℏ(∂/∂t)ψ± = (-ℏ²/2m ∇² + mV)Iψ± - μB·Sψ± — this is the proposed Pauli-like equation. The first term is the standard Schrödinger equation (with V being the gravitational potential); the second term couples the "spin" degrees of freedom to the tidal field. The wavefunction is now a 2-spinor ψ± = (ψ₊, ψ₋)ᵀ.

53. **p210**: "In the actual Pauli equation μ = qℏ/2m, in this formulation I expect the multiplier of U to be the same as it is for V: hence, μ = m." — This reasoning is not convincing. In quantum mechanics, μ = qℏ/2m has a specific derivation from minimal coupling. Simply asserting μ = m by analogy doesn't follow. The coupling constant should be derived from the equations of motion, not guessed. The dimensional analysis also needs checking: if B = ∂ᵢφ has dimensions of acceleration (m/s²) and S has dimensions of mass × length (from ∫ρd dV), then B·S has dimensions of force × length = energy, which is correct for a potential energy term. But the coupling μ needs to make the units work in the Schrödinger equation.

54. **p210**: "If the new term, U = B·S, is zero then the above equation would reduce to two Schrödinger equations, which are essentially independent except for the mutual gravitational interaction." — Correct, and this is an important limiting case. It connects back to Johnston's two-fluid approach.

### Section 7.2.3 — Quadrupole Term and Higher (pp210-213)

55. **p211**: Equations 7.38-7.40: The quadrupole extension with B = ∂ᵢⱼ(φ) and S = ∫ρ(r eᵢ)(r eⱼ) dV — the approach of using the Clifford algebra isomorphism again to map the rank-2 quadrupole tensor to a 2×2 matrix via products of Pauli matrices is creative.

56. **p211**: Equation 7.41: The quadrupole matrix Q_ij with entries like r_xxσ_xx, r_xyσ_xy, etc. — the notation is becoming very compressed. The entries involve *products* of Pauli matrices (e.g., σ_xxσ = σ_x² = I, σ_xyσ = σ_xσ_y = iσ_z). The reader would benefit from seeing at least one explicit example worked out.

57. **p212**: Equations 7.42-7.45: The Pauli matrix multiplication rules — correct:
   σ₁² = I, σ₁σ₂ = iσ₃ = -σ₂σ₁, σ₂σ₃ = iσ₁ = -σ₃σ₂, σ₃σ₁ = iσ₂ = -σ₁σ₃.
   These are the standard su(2) / quaternion relations.

58. **p212**: Equation 7.46: Q_ij with entries like r_xxI (diagonal) and ir_xyσ_z (off-diagonal) — this follows from applying Eqs 7.42-7.45 to Eq 7.41. The key result: diagonal entries are proportional to the identity (scalar), off-diagonal entries are proportional to Pauli matrices (vector-like). This means the quadrupole contribution has both "scalar" and "spin" parts.

59. **p212**: Equations 7.47-7.51: The 2,2-spinor components X^{AÛ} — this notation from Misner, Thorne & Wheeler is heavy for this context. The result (Eq 7.52) is that the final Pauli equation has the same form as Eq 7.36 but with a more complicated B·S term. The presentation could be streamlined by simply stating the result and pointing to the derivation steps.

60. **p212**: Equation 7.52: iℏ(∂/∂t)ψ± = (-ℏ²/2m ∇² + mV)Iψ± - μB·Sψ± — identical in form to Eq 7.36, confirming that the quadrupole just modifies the B·S coupling. This is an elegant result.

61. **p213**: The chapter ends rather abruptly with "Achieving the same physical system in an N-body code could be harder; the force at quadrupole order (equation 7.18) requires a third derivative to be applied to the gravitational potential φ." — While true, a proper concluding paragraph summarising the chapter's contributions would be more satisfying.

---

## Equations Check

- **Eq 7.1**: F_g = m(E_g + v × 2B_g) — correct GEM Lorentz force (factor of 2 is standard in GEM)
- **Eq 7.2**: iℏ∂ψ/∂t = -(ℏ²/2m)(∇ - ie/ℏ A)²ψ — correct (minimal coupling Schrödinger equation)
- **Eq 7.3**: B_g = (G/2c²)(L - 3(L·r̂)r̂)/r³ — correct (gravitomagnetic dipole field)
- **Eq 7.4**: L = ∫ρ(v × r)dV = Σm(v × r) — **sign error**: should be r × v, not v × r. L = ∫ρ(r × v)dV
- **Eq 7.5**: U(R+r) = φ₀∫ρ dV + ∂ᵢφ∫ρrᵢ dV + ... — correct (multipole expansion)
- **Eq 7.6-7.7**: Dipole moment for two equal masses at ±x — correct, gives D = 0
- **Eq 7.8-7.9**: Quadrupole moment — correct, gives Q_xx = 2mx²
- **Eq 7.10-7.12**: Dipole potential U^(d) = ∂ᵢφ · dᵢ — correct
- **Eq 7.13-7.16**: Quadrupole potential U^(q) = ∂ᵢⱼφ · Q_ij — correct (full expansion)
- **Eq 7.17**: Force equation F = ∇⊗U — correct in form
- **Eq 7.18-7.20**: Quadrupole force F^(q) — correct, involves third derivative of φ
- **Eq 7.21-7.22**: Force on each sub-particle — correct
- **Eq 7.23-7.24**: Sum and difference of forces — correct (centre-of-mass + relative motion)
- **Eq 7.25**: Difference vector for dumbbell — correct
- **Eq 7.26**: d = r_j^(1) - r_j^(2) — correct (definition)
- **Eq 7.27**: U ≃ B·S — analogy, not a derivation; see note #48
- **Eq 7.28-7.29**: B = ∂ᵢφ, S = ∫ρd dV — definitions
- **Eq 7.30**: S ~ ∫ρdᵢσᵢ dV — Clifford algebra mapping
- **Eq 7.31**: B·S ~ ∫ρ∂ᵢ(φ)(d σ)ᵢ dV — correct
- **Eq 7.32-7.34**: Matrix form of B·S — **correct** (verified by explicit Pauli matrix computation)
- **Eq 7.35**: Pauli matrices — correct (standard definitions)
- **Eq 7.36**: Pauli-like Schrödinger equation — correct in form
- **Eq 7.37**: 2-spinor ψ± = (ψ₊, ψ₋)ᵀ — correct
- **Eq 7.38-7.39**: Quadrupole B and S — correct (rank-2 generalisations)
- **Eq 7.40**: Quadrupole potential U^q = B··S — correct
- **Eq 7.41**: Q_ij with Pauli matrix products — correct in structure
- **Eq 7.42-7.45**: Pauli algebra — correct
- **Eq 7.46**: Q_ij after Pauli algebra simplification — correct
- **Eq 7.47-7.51**: 2,2-spinor components — plausible, dense calculation not fully verified
- **Eq 7.52**: Final Pauli equation with quadrupole — correct in form (same as 7.36)

---

## Structural Comments

62. The chapter has a clear two-part structure: vorticity/GEM (7.1) and spin/Pauli (7.2). Both sections are self-contained and could be read independently.

63. The GEM section (7.1) is shorter and more grounded — it connects to existing physics (GR weak-field limit) and includes a numerical test using Hydra data. This is the stronger section.

64. The Pauli equation section (7.2) is more ambitious and original. The progression from multipole expansion → dipole model → Clifford algebra → Pauli equation is logical and creative. However, the derivation is quite compressed and would benefit from more intermediate steps, particularly in the Clifford algebra sections.

65. **Missing**: The chapter has no proper conclusion section. It ends mid-discussion after stating the quadrupole result. A concluding section summarising both ideas, their potential impact, and concrete next steps for implementation would greatly improve the chapter.

66. **Missing**: No discussion of how the Pauli equation (7.36/7.52) would actually be *solved* numerically. The Schrödinger equation for a scalar wavefunction was challenging enough (Chapter 5); a 2-spinor equation doubles the degrees of freedom and introduces coupling between components. What would the Cayley scheme look like for a 2×2 matrix equation?

67. **Connection to rest of thesis**: The chapter effectively bridges two separate research directions (GEM from GR, Pauli equation from quantum mechanics) and maps them onto the wave-mechanics framework developed in Chapters 3-5. This is intellectually ambitious and, despite the speculative nature, demonstrates the versatility of the wave-mechanical approach.

68. The notation B is overloaded: in section 7.1, B_g is the gravitomagnetic field; in section 7.2, B = ∂ᵢφ is the gradient of the scalar gravitational potential (which is actually E_g in GEM notation!). This is confusing and should use different symbols.
