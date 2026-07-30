# Future direction: rotor and spinor internal states for structured gravitating objects

**Status: speculative proposal, not a result.** This note develops a possible future
direction for the spin extension proposed in Chapter 7 (Epilogue: Vorticity and Spin)
of the thesis. The multipole setup, the vanishing gravitational dipole, and the tidal
torque are standard and reliable. The rotor formulation is a reframing of the chapter's
own idea on firmer footing. The multicomponent Schrodinger-Poisson equations at the end
are schematic: they are written to show what a self-consistent implementation might look
like, and still need to be derived carefully, given coupling constants, and tested
numerically. Nothing here has been implemented.

## Relationship to the thesis

Chapter 7 proposes that a wave-mechanical super-particle could carry internal spin,
built by analogy with the quantum-mechanical Pauli equation: a two-component spinor
wavefunction coupled to a field through a term of the form `mu * B . S`. That chapter
identifies the right physical ingredients (the gravitational dipole vanishes in the
centre-of-mass frame, so the first non-trivial coupling is at quadrupole order through
the tidal tensor) but forces the result into a Pauli form and fixes the coupling
constant `mu = m` by analogy rather than derivation.

The purpose of this note is to make one correction of principle and follow it through:

- Gravity couples an extended body to the **rank-2 tidal tensor** `T_ij`, not to a
  rank-1 field. A linear `B . sigma` Pauli term is therefore not the natural
  fundamental object.
- A 2-component spinor carries only **two** physical rotational degrees of freedom (a
  single axis). That is exactly enough for an axisymmetric body and not enough for a
  general triaxial one.
- The honest internal variable is a **rotor** (unit quaternion) that tracks the full
  body frame. The Pauli 2-spinor then reappears as the axisymmetric special case, not
  as the starting point.

## Multipole setup

Let `R` be the centre of mass of a super-particle and `s` the internal displacement
from it. With internal mass distribution `rho(s)` in an external potential
`phi_ext(R + s)`, the interaction energy is

$$
U(\mathbf{R}) = \int \rho(\mathbf{s})\, \phi_{\mathrm{ext}}(\mathbf{R}+\mathbf{s})\, d^3 s
= M\,\phi_{\mathrm{ext}}(\mathbf{R})
+ d_i\,\partial_i \phi_{\mathrm{ext}}
+ \tfrac{1}{2} Q_{ij}\,\partial_i \partial_j \phi_{\mathrm{ext}}
+ \cdots ,
$$

$$
M = \int \rho\, d^3 s, \qquad
d_i = \int \rho\, s_i\, d^3 s, \qquad
Q_{ij} = \int \rho\, s_i s_j\, d^3 s .
$$

For a gravitational mass distribution in its own centre-of-mass frame the dipole moment
vanishes, `d_i = 0`, because all mass "charges" have the same sign (Jackson's
electromagnetic treatment of an extended charged body, adapted to gravity). The first
non-trivial coupling to internal structure is therefore the quadrupole term, built from
the tidal tensor `T_ij = partial_i partial_j phi_ext`.

## Torque and the axisymmetric reduction

Writing the local force density as `-rho(s) grad phi_ext(R + s)` and expanding the
gradient to first order in `s`,

$$
\nabla \phi_{\mathrm{ext}}(\mathbf{R}+\mathbf{s}) \approx
\nabla \phi_{\mathrm{ext}}(\mathbf{R}) + T\cdot\mathbf{s},
\qquad T_{ij} = \partial_i \partial_j \phi_{\mathrm{ext}}(\mathbf{R}),
$$

the zeroth-order term produces no torque (the dipole vanishes), and the leading torque is

$$
\tau_i = -\,\epsilon_{ijk}\, T_{k\ell}\, Q_{j\ell}.
$$

This is the clean Newtonian statement of the mechanism: a non-spherical body spins up
because the tidal tensor is not aligned with its quadrupole tensor. (Note this is the
correct torque law; it does not carry the spurious factor of two that appears in the
thesis chapter's dumbbell equation of motion, which has since been corrected in the
published chapter.)

For an axisymmetric toy model with symmetry axis `n`,

$$
Q_{ij} = q\, n_i n_j, \qquad
U_{\mathrm{spin}} = \tfrac{1}{2} q\, n_i T_{ij} n_j, \qquad
\boldsymbol{\tau} = -\,q\, \mathbf{n}\times(T\cdot\mathbf{n}).
$$

This provides a minimal internal degree of freedom (one axis) that can be evolved
numerically, either by attaching an orientation vector to each super-particle in an
N-body code, or, in wave-mechanics, by carrying it in an auxiliary internal state.

## Rotor formulation (the general case)

For a fully oriented body the honest internal variable is a rotor rather than a
2-spinor. In three dimensions the even Clifford algebra is isomorphic to the
quaternions, `Cl_3^+ ~ H`, and unit rotors double-cover `SO(3)` (Lounesto). A rotor `R`
satisfies `R R~ = 1` (with `R~` the reversion), and a body-fixed orthonormal frame
`{e_1, e_2, e_3}` maps to lab-frame axes

$$
E_a = R\, e_a\, \tilde{R}, \qquad a = 1,2,3 .
$$

The lab-frame quadrupole of a body with principal-frame quadrupole
`Q_body = diag(q_1, q_2, q_3)` is then

$$
Q_{ij}(R) = \sum_{a=1}^{3} q_a\, (E_a)_i (E_a)_j
= \big(\mathcal{R}(R)\, Q_{\mathrm{body}}\, \mathcal{R}(R)^{T}\big)_{ij},
$$

with `R(R)` the `SO(3)` matrix of the rotor, and the interaction energy keeps the same
quadrupole form,

$$
U_{\mathrm{int}}[R] = \tfrac{1}{2} Q_{ij}(R)\, T_{ij}.
$$

This is the mathematically honest generalisation of the axisymmetric model: the internal
variable is a full orientation state and the force and torque laws still descend from
the tidal coupling rather than from a guessed linear `B . sigma` term. The Clifford /
Pauli-matrix representation (`Cl_3 ~ Mat(2, C)`, `e_i <-> sigma_i`; see Lounesto, and
Misner, Thorne & Wheeler for the spinor conventions) remains useful as a representation
tool, but it no longer carries the physics on its own.

### Why the 2-spinor only covers the axisymmetric case

A 2-component complex spinor has four real parameters; removing the overall norm and
phase leaves two physical degrees of freedom, a point on the Bloch sphere, i.e. a single
axis direction `n_i = (xi^dagger sigma_i xi)/(xi^dagger xi)`. That is exactly an
axisymmetric body with `Q_ij = q n_i n_j`. A general triaxial body needs three
rotational degrees of freedom plus two trace-free shape parameters, which a 2-spinor
cannot encode. Hence: 2-spinor for axisymmetric objects, rotor/quaternion (plus
principal quadrupole strengths) for the general case.

## Multicomponent Schrodinger-Poisson extension (schematic)

**The equations in this section are conjectural.** They indicate the shape a
self-consistent wave-mechanical implementation might take; the internal Hamiltonian
`H_int`, the measure `d mu(R)` over rotor space, and the coupling constants are not yet
derived, and the reductions should be checked against the Chapter 5 formulation before
any use.

In physical variables, the scalar Schrodinger-Poisson system of the thesis is

$$
i\hbar\,\partial_t \psi = \Big[-\tfrac{\hbar^2}{2m}\nabla^2 + m\Phi\Big]\psi,
\qquad \nabla^2 \Phi = 4\pi G \rho, \quad \rho = m|\psi|^2 .
$$

A natural multicomponent generalisation carries an internal rotor coordinate `R`:

$$
i\hbar\,\partial_t \Psi = \Big[-\tfrac{\hbar^2}{2m}\nabla^2 + m\Phi
+ \hat{H}_{\mathrm{int}}\big(T[\Phi], R\big)\Big]\Psi,
\qquad
\nabla^2 \Phi = 4\pi G\, m \int \Psi^\dagger(\mathbf{x},t,R)\,\Psi(\mathbf{x},t,R)\, d\mu(R),
$$

with the internal Hamiltonian built from the self-consistent tidal tensor
`T_ij[Phi] = partial_i partial_j Phi`, encoding the quadrupole coupling
`U_int = (1/2) Q_ij(R) T_ij[Phi]`. If the internal state factorises,
`Psi = psi(x,t) Xi(R,t)`, with unit norm `int Xi^dagger Xi d mu(R) = 1`, and `H_int` is
absent or frozen, the density reduces to `rho = m |psi|^2` and the scalar
Schrodinger-Poisson system is recovered. The proposal is thus an **extension** of the
thesis framework, not a replacement.

### Comoving form (to match Chapter 5)

In the rescaled comoving variables of Chapter 5 (comoving coordinate `y`, comoving
wavefunction `chi`, comoving potential `U`, dimensionless parameter `L`), the scalar
system is, schematically,

$$
2 i \mathcal{L}\big[\Omega_{m0} + (1-\Omega_{m0})a^3\big]^{1/2}
\frac{\partial \chi}{\partial \ln a}
= \Big[-\nabla_{\mathbf{y}}^2 + 3\Omega_{m0}\mathcal{L}^2 U\Big]\chi,
\qquad \nabla_{\mathbf{y}}^2 U = \chi\chi^\ast - 1 ,
$$

and the rotor-extended version would add `H_int(T[U], R)` inside the bracket and replace
`chi chi*` in the Poisson source by `int Psi^dagger Psi d mu(R)`. The Einstein-de Sitter
case is `Omega_m0 = 1`. These forms should be verified against the exact normalisation
in Chapter 5 before implementation.

## Open problems

- Derive `H_int` (and its coupling constant) from the rigid-body / tidal model rather
  than by analogy; the thesis chapter's `mu = m` is a placeholder.
- Define the internal state space and measure `d mu(R)` concretely (a sharp rotor state,
  a wavepacket on `SO(3)`, or a small basis).
- Reduce the rotor model to the Pauli 2-spinor explicitly for the axisymmetric case, by
  linearising the tidal interaction about a reference orientation, to confirm the
  chapter's Pauli equation as a controlled limit.
- Verify the multicomponent Schrodinger-Poisson forms against the Chapter 5
  normalisation, then test the torque law on simple controlled configurations (a single
  dumbbell in a fixed external tide) before attempting a cosmological run.
- Assess whether the effect is ever dynamically relevant for large-scale structure, or
  only for compact / rapidly tidally forced systems.

## References

These correspond to entries in the thesis bibliography (`biblio.bib`); keys are given so
they can be wired up directly if this note is ever promoted to LaTeX.

- J. D. Jackson, *Classical Electrodynamics* (multipole expansion of an extended body) — `bib_jack`
- P. Lounesto, *Clifford Algebras and Spinors* (the isomorphisms `Cl_3^+ ~ H`, `Cl_3 ~ Mat(2, C)`) — `bib_lou`
- C. W. Misner, K. S. Thorne, J. A. Wheeler, *Gravitation* (spinor conventions, chapter 41) — `bib_grav`
- Watanabe & Tsukada (2000), Schrodinger evolution with a vector potential via splitting operators — `bib_wat2`
- C. J. Short, PhD thesis (circulation theorem and vorticity in the wave-mechanical picture) — `bib_short`
- R. Johnston, two-component ("two-fluid") Schrodinger evolution — `bib_rj`

See Chapter 7 of the thesis (`thesis/latex/ch07_epilogue_vorticity_spin.tex`) for the
original proposal this note extends.
