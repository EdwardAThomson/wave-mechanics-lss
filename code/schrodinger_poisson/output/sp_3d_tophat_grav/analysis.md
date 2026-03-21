# Tophat Collapse Test: Analysis

## Parameters

- Grid: N = 64^3, box = [0,1)^3, epsilon = 1/64
- Timestep: delta = 1e-4, lambda = 9.77
- G_eff = 5000
- Tophat: R = 0.2, delta = 0.5 (rho_in = 1.5, rho_out = 0.983)
- Edge smoothing: 3 grid cells (tanh profile)
- No expansion (static universe)
- Total time: 0.05 (500 steps)

## Results

- Mass conserved to ~5e-13 throughout
- Tophat contracts under self-gravity as expected
- Central density increases from ~1.5 to ~1200 by t = 0.05
- Collapse maintains spherical symmetry despite split-operator treatment of x, y, z
- Central spike has finite width (~few cells) due to quantum pressure

## Analytical comparison

### Freefall time

For a uniform overdense sphere, the gravitational potential inside satisfies:

    nabla^2 V = rho - rho_mean = delta_rho = 0.5
    V(r) = delta_rho * r^2 / 6     (spherical solution)

Acceleration at the surface:

    a = -G_eff * dV/dr|_{r=R} = -G_eff * delta_rho * R / 3
      = -5000 * 0.5 * 0.2 / 3
      = -167

Estimated freefall time (constant-acceleration, surface to centre):

    t_ff = sqrt(2R / |a|) = sqrt(2 * 0.2 / 167) = sqrt(0.0024) ~ 0.049

This matches the observed collapse timescale well: the strong central
concentration builds up between t = 0.03 and t = 0.05.

### Quantum Jeans length

The quantum Jeans length determines the scale below which quantum
pressure resists gravitational collapse:

    lambda_J ~ sqrt(4 / (G_eff * rho)) = sqrt(4 / (5000 * 1.5)) ~ 0.023

In grid cells: 0.023 * 64 ~ 1.5 cells.

Since the tophat radius (R = 0.2 ~ 13 cells) is much larger than the
Jeans length, gravitational collapse proceeds at nearly the classical
rate. Quantum pressure only becomes important at the centre, giving
the density spike a finite width rather than a true singularity.

## Conclusions

The gravitational solver produces physically correct collapse dynamics
in the non-expanding case. The collapse timescale and the role of
quantum pressure both match analytical expectations.

Next steps: add cosmological expansion and repeat the tophat test to
validate the combined gravity + expansion system.
