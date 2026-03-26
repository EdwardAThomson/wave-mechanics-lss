# GADGET-4 vs Schrödinger-Poisson: Statistical Comparison

## Setup

Both codes evolve the same initial conditions (extracted from GADGET-4's
NGENIC IC generator via TSC deposit + Poisson-derived velocity potential)
through identical LCDM cosmology.

| Parameter | Value |
|---|---|
| Box size | 32 Mpc/h |
| Grid / particles | 128³ |
| Omega_m | 0.27 |
| Omega_Lambda | 0.73 |
| h | 0.71 |
| sigma_8 | 0.81 |
| z_init | 50 |
| S-P nu | 3e-4 |
| S-P method | Spectral kinetic step + Cayley potential |
| GADGET-4 | TreePM, 2LPT ICs, adaptive timestepping |
| GADGET density | TSC deposit of particles onto 128³ grid |

## Available snapshots

| Snapshot | z (GADGET) | z (S-P matched) | r (correlation) |
|---|---|---|---|
| 001 | 19.0 | 19.0 (a=0.043) | 0.615 |
| 002 | 9.0 | 9.6 (a=0.095) | 0.421 |
| 003 | 5.7 | 6.1 (a=0.140) | 0.177 |
| 004 | 4.0 | 3.8 (a=0.208) | 0.138 |
| 005 | 2.3 | 2.3 (a=0.307) | 0.182 |
| 006+ | Pending | — | — |

GADGET-4 is still running (currently at z~1). Snapshots at a=0.4
through a=1.0 (z=1.5 to z=0) pending.

## Results

### Power spectrum P(k)

At **z=19** (most linear epoch):
- P(k) agrees well at low k (large scales, k < 0.2 h/Mpc)
- S-P shows excess power at high k (small scales)
- The cross-correlation r(k) is close to 1 at low k

At **z < 10** (nonlinear regime):
- S-P develops progressively more power at all scales
- The transfer function T(k) = sqrt(P_SP/P_G4) exceeds 1 everywhere
- r(k) drops below 0.5 at most scales by z=4

Interpretation: the S-P code produces more small-scale structure than
GADGET-4. This is consistent with two effects:
1. **Phase saturation** at nu=3e-4: the wavefunction cannot represent
   velocities beyond v_max = nu × pi × N = 0.12, causing kinetic energy
   to be redistributed as density structure rather than bulk flow.
2. **Quantum pressure** in underdense regions: the Bohm potential
   evacuates voids faster than classical gravity.

### Density PDF

At **z=19**: both codes produce nearly identical Gaussian distributions
centred on delta=0. This confirms the ICs are correctly matched.

At **z < 5**: the distributions diverge:
- S-P develops a broader high-density tail (sharper peaks, more
  concentrated structures)
- S-P has a sharper low-density cutoff (emptier voids)
- GADGET maintains a more symmetric, gradually broadening distribution

This is qualitatively similar to what the thesis (2011) reported in
the comparison with GADGET-2 (Chapter 5, Figures 5.10-5.11), but now
with matched ICs we can confirm it is a real physical difference rather
than a random-seed artefact.

### rms(delta) growth

- Both codes agree at z > 15 (rms ~ 0.3)
- S-P grows faster than GADGET from z=10 onwards
- By z=2.3: rms_SP = 6.2 vs rms_G4 = 4.0 (factor 1.5x)
- The excess growth is due to both phase saturation and the lack of
  sub-grid force resolution in the S-P code

### Point-by-point correlation

- Peaks at r=0.62 (z=19)
- Drops to r~0.15-0.2 by z=4
- Stabilises at late times (the large-scale structure is set)
- The remaining correlation comes from the fundamental mode (k=1)
  which both codes evolve correctly

## Key findings

1. **The codes agree in the linear regime.** At z=19, both the power
   spectrum and density PDF match, confirming the IC matching and the
   S-P growth calibration are correct.

2. **The S-P produces excess small-scale power.** This is primarily
   the nu×N phase saturation effect: kinetic energy that cannot be
   represented as velocity is redistributed as density structure.

3. **Voids evacuate faster in S-P.** The quantum pressure (Bohm
   potential) is physical at this nu but absent in GADGET. This is
   a prediction of the wave-mechanical framework.

4. **The density PDF comparison reproduces the thesis finding.** The
   wave-mechanical code gives "more lows and more highs" than the
   N-body code, exactly as reported in Thomson (2011) Chapter 5.

## Next steps

- Complete GADGET-4 run to z=0 and update comparison
- Run S-P with larger nu (1e-3) on matched ICs for comparison with
  less phase saturation
- Compute velocity statistics once GADGET velocity units are resolved
- Try N=256 to push nu smaller and reduce phase saturation
- Power spectrum ratio P_SP/P_G4 as a function of both k and z
  (2D transfer function plot)
