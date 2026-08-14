// Unit system for the disk bending-wave solver.
//
// Length   : kpc
// Velocity : km/s
// Mass     : Msun
// Time     : kpc / (km/s) = 977.792 Myr
//
// In these units hbar_eff = hbar/m has dimensions of kpc*km/s, which is the
// natural phase-space coarse-graining scale: a phase-space cell of area
// hbar_eff in the (z, v_z) plane.
#pragma once

namespace units {

constexpr double PI = 3.14159265358979323846;

// Newton's constant in kpc (km/s)^2 / Msun
constexpr double G = 4.300917270e-6;

// One time unit in Myr
constexpr double TIME_MYR = 977.79222;
constexpr double TIME_GYR = 0.97779222;

// Convenience: Msun/pc^2 -> Msun/kpc^2
constexpr double MSUN_PC2 = 1.0e6;
// Msun/pc^3 -> Msun/kpc^3
constexpr double MSUN_PC3 = 1.0e9;

}  // namespace units
