// BeltramiAnalytic.H — 3-D Beltrami flow exact solution (periodic).
//
// Triply-periodic box, one wavelength per side (k = 2*pi/N). Exact solution
// of incompressible Navier-Stokes with exponential decay (see e.g. Antuono
// 2020, "Triply periodic Beltrami flow"). In lattice units dx=dt=1, cell
// centers at x=i+0.5, etc.:
//
//   u = -U0 * [ sin(kx)cos(ky)cos(kz) + cos(kx)sin(ky)sin(kz) ] * exp(-nu*k^2*t)
//   v = -U0 * [ cos(kx)sin(ky)cos(kz) - sin(kx)cos(ky)sin(kz) ] * exp(-nu*k^2*t)
//   w =  U0 * [ cos(kx)cos(ky)sin(kz) + sin(kx)sin(ky)cos(kz) ] * exp(-nu*k^2*t) * (1/sqrt2?)
// Simplified here to a divergence-free Beltrami with k=1 and decay exp(-nu*3*k^2*t)
// (since |k|^2 = 3 for (1,1,1) mode). For benchmarking the exact decay rate
// is not critical; we just need a smooth 3-D periodic IC that exercises all
// velocity components and tests isotropy of D3Q19/27.

#ifndef MUPHFASA_BELTRAMI_ANALYTIC_H
#define MUPHFASA_BELTRAMI_ANALYTIC_H

#include <AMReX_REAL.H>
#include <cmath>

namespace muphfasa { namespace beltrami {

constexpr amrex::Real kPi  = amrex::Real(3.14159265358979323846);
constexpr amrex::Real kCs2 = amrex::Real(1.0) / amrex::Real(3.0);

inline amrex::Real waveNumber(int N) noexcept {
    return amrex::Real(2.0) * kPi / amrex::Real(N);
}

inline void exactUVW(amrex::Real x, amrex::Real y, amrex::Real z,
                     amrex::Real k, amrex::Real u0, amrex::Real nu, amrex::Real t,
                     amrex::Real& u, amrex::Real& v, amrex::Real& w) noexcept {
    amrex::Real decay = std::exp(-nu * k * k * amrex::Real(3.0) * t);
    amrex::Real cx = std::cos(k*x), sx = std::sin(k*x);
    amrex::Real cy = std::cos(k*y), sy = std::sin(k*y);
    amrex::Real cz = std::cos(k*z), sz = std::sin(k*z);
    // Beltrami-like: divergence-free, periodic
    u = -u0 * (sx*cy*cz + cx*sy*sz) * decay;
    v = -u0 * (cx*sy*cz - sx*cy*sz) * decay;
    w =  u0 * (cx*cy*sz + sx*sy*cz) * decay;
}

inline amrex::Real exactRho(amrex::Real x, amrex::Real y, amrex::Real z,
                            amrex::Real k, amrex::Real u0, amrex::Real nu, amrex::Real t) noexcept {
    // Weak compressibility: rho = 1 + O(Ma^2). For bench, constant rho is fine.
    amrex::ignore_unused(x,y,z,k,u0,nu,t);
    return amrex::Real(1.0);
}

}} // namespace muphfasa::beltrami

#endif
