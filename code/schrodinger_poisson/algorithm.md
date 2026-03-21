# Schrodinger-Poisson Solver: Algorithm Reference

Derived from Chapter 5 of Thomson (2011), sections 5.0-5.5.
This document captures the computational algorithm for solving the coupled
Schrodinger-Poisson system with:

1. 3D coordinates
2. Self-consistent gravity
3. Expanding coordinates (flat LCDM)
4. Periodic boundaries (double recursion method)
5. Mass conservation (unitary Cayley transform)

## Key Parameters

```
N         = number of grid points per dimension
epsilon   = 1/N                          (grid spacing, box units [0,1))
delta     = d(ln a)                      (timestep in ln(a))
lambda    = 4 * epsilon^2 / delta        (Courant-like parameter)
L_phys    = physical box size in Mpc
nu        = hbar / m                     (de Broglie parameter)
```

### Expansion-dependent quantities (flat LCDM)

From thesis Eq 5.55 (section 5.2.2):

```
H(a)      = H0 * sqrt(Omega_m0 / a^3 + Omega_Lambda)
L_param   = m * sqrt(a) * L_phys^2 * H0 / hbar
```

### Dimensionless variables

From thesis Eqs 5.53-5.54:

```
chi       = sqrt(8*pi*G / (3*H0^2*Omega_m0)) * psi
U         = (2*a / (3*L_phys^2*H0^2*Omega_m0)) * V
y         = x / L_phys
```

The dimensionless Schrodinger equation (thesis Eq 5.55):

```
2i * L_param * [Omega_m0 + (1 - Omega_m0)*a^3]^{1/2} * d(chi)/d(ln a)
    = [-nabla^2 + 3*Omega_m0*L_param^2 * U] * chi
```

The dimensionless Poisson equation:

```
nabla^2 U = chi*chi* - 1
```

## Time-stepping: Strang Splitting (KVK)

The evolution operator is split as (thesis Eq 5.35, section 5.3.1):

```
psi(t+dt) = exp(-i*K*dt/2) * exp(-i*V*dt) * exp(-i*K*dt/2) * psi(t)
```

Since [K_x, K_y] = [K_x, K_z] = [K_y, K_z] = 0, each kinetic half-step
decomposes exactly into three sequential 1D sweeps (thesis Eq 5.37):

```
exp(-i*K*dt/2) = exp(-i*K_x*dt/2) * exp(-i*K_y*dt/2) * exp(-i*K_z*dt/2)
```

## Computational Algorithm

From thesis section 5.4 (pp 130-131):

```
1. Construct initial wavefunction
2. Start time loop
3.   First half-step kinetic energy (Cayley, x then y then z):
       a. calculate auxiliary function e  (double recursion for periodic BCs)
       b. calculate Omega
       c. calculate auxiliary function f  (double recursion for periodic BCs)
       d. update psi                      (double backward recursion)
4.   Full-step potential energy:
       a. compute density rho = |psi|^2
       b. solve Poisson equation via FFT for V
       c. update psi via diagonal Cayley rotation
5.   Second half-step kinetic energy (same as step 3)
6. End time loop
```

## Kinetic Half-Step: Goldberg Method with Periodic BCs

Each kinetic half-step is applied as a 1D sweep along one dimension.
For a 3D grid, we loop over all "pencils" (1D lines) along that dimension.

The Cayley transform for the kinetic-only operator gives (thesis Eq 5.38):

```
e = 2 - i*lambda       (no V term, since V is handled separately)
```

### Omega (no recursion needed)

Thesis Eq 5.15 (adapted for kinetic-only, no V term):

```
for j = 0 to N-1:
    j_plus  = (j + 1) mod N
    j_minus = (j - 1 + N) mod N
    Omega[j] = -psi[j_plus] + (i*lambda + 2)*psi[j] - psi[j_minus]
```

Note: the indices wrap periodically since we use periodic BCs.

### Auxiliary function e (double recursion)

First pass (thesis Eqs 5.20-5.22):

```
e[0] = 2 - i*lambda
for j = 1 to N-1:
    e[j] = 2 - i*lambda - 1/e[j-1]
```

Second pass -- periodic wrapping (thesis Eqs 5.23-5.25).
Seed with last value from first pass:

```
e[0] = 2 - i*lambda - 1/e[N-1]
for j = 1 to N-1:
    e[j] = 2 - i*lambda - 1/e[j-1]
```

### Auxiliary function f (double recursion)

First pass (thesis Eqs 5.26-5.27):

```
f[0] = Omega[0]
for j = 1 to N-1:
    f[j] = Omega[j] + f[j-1] / e[j-1]
```

Second pass -- periodic wrapping (thesis Eqs 5.28-5.29):

```
f[0] = Omega[0] + f[N-1] / e[N-1]
for j = 1 to N-1:
    f[j] = Omega[j] + f[j-1] / e[j-1]
```

### Wavefunction update (double backward recursion)

First backward pass (thesis Eq 5.30):

```
psi_new[N-1] = (psi_new[0] - f[N-1]) / e[N-1]
for j = N-2 down to 0:
    psi_new[j] = (psi_new[j+1] - f[j]) / e[j]
```

Second backward pass (periodic convergence):

```
psi_new[N-1] = (psi_new[0] - f[N-1]) / e[N-1]
for j = N-2 down to 0:
    psi_new[j] = (psi_new[j+1] - f[j]) / e[j]
```

The double recursion is the thesis author's original contribution for
periodic boundaries. Two passes suffice; a third pass gives identical
results (demonstrated in thesis Figure 5.1).

## Potential Energy Step

### Poisson solve via FFT

Compute density and solve for the potential in Fourier space.
Uses the discrete cosine correction (thesis Eq 5.33, section 5.5):

```
rho[j] = |psi[j]|^2
rho_mean = mean(rho)

FFT: rho_hat = FFT(rho - rho_mean)

for each mode (n, m, o):
    kappa = cos(2*pi*n/N) + cos(2*pi*m/N) + cos(2*pi*o/N)
    V_hat[n,m,o] = rho_hat[n,m,o] / (2*kappa - 3)

V_hat[0,0,0] = 0    (mean potential is arbitrary)

IFFT: V = IFFT(V_hat)

Scale V by dimensionless factors as appropriate for the expansion model.
```

### Wavefunction update (diagonal Cayley)

The potential operator is diagonal in real space, so the Cayley transform
is a simple pointwise operation (thesis Eq 5.36):

```
for each grid point j:
    psi[j] *= (1 - i*delta/2 * V[j]) / (1 + i*delta/2 * V[j])
```

This is a pure phase rotation of magnitude 1, so it is exactly unitary.

## Notes

- The kinetic step uses finite-difference tridiagonal recursion (Goldberg
  method), NOT spectral/FFT methods. Only the Poisson solve uses FFTs.
- The potential V does not appear in the kinetic step auxiliary functions
  (e, f, Omega) because the splitting operator separates K and V.
- Consequently, the `e` recursion for the kinetic step is particularly
  simple: `e[j] = 2 - i*lambda - 1/e[j-1]` with no spatially varying term.
- Expansion enters through L_param(a), which changes each timestep.
- The timestep delta = d(ln a) is uniform in ln(a), not in a or t.
- Mass conservation: guaranteed by unitarity of the Cayley transform,
  provided the double recursion is used for periodic BCs.
