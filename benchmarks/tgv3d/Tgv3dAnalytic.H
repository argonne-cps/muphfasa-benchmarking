// Tgv3dAnalytic.H — Shared analytic solution for the 3-D Taylor-Green vortex.
//
// Triply-periodic domain, one vortex wavelength per domain (k=2*pi/N, N
// cells per side). 3-D extension of benchmarks/tgv2d/TgvAnalytic.H (Kruger
// ch. 10.4, Brachet et al. 1983): classic 3-D TGV with w=0 at t=0.
//
//   u(x,y,z,t) =  U0 * sin(k*x) * cos(k*y) * cos(k*z) * exp(-2*nu*k^2*t)
//   v(x,y,z,t) = -U0 * cos(k*x) * sin(k*y) * cos(k*z) * exp(-2*nu*k^2*t)
//   w(x,y,z,t) =  0
//   rho(x,y,z,t) = 1 - (U0^2/(4*cs2)) * (cos(2*k*x)+cos(2*k*y)) * (cos(2*k*z)+2)/? * exp(-4*nu*k^2*t)
//                For bench we use rho=1 (weak compressibility, Ma<<1) like beltrami.
//
// Used by both InitData.cpp (t=0) and validate/test_tgv3d.cpp (t>0).
// Lattice units dx=dt=1, cell centers at x=i+0.5 etc.

#ifndef MUPHFASA_TGV3D_ANALYTIC_H
#define MUPHFASA_TGV3D_ANALYTIC_H

#include <AMReX_REAL.H>
#include <cmath>

namespace muphfasa { namespace tgv3d {

constexpr amrex::Real kPi  = amrex::Real(3.14159265358979323846);
constexpr amrex::Real kCs2 = amrex::Real(1.0) / amrex::Real(3.0);

inline amrex::Real waveNumber(int N) noexcept {
    return amrex::Real(2.0) * kPi / amrex::Real(N);
}

inline void exactUVW(amrex::Real x, amrex::Real y, amrex::Real z,
                     amrex::Real k, amrex::Real u0, amrex::Real nu, amrex::Real t,
                     amrex::Real& u, amrex::Real& v, amrex::Real& w) noexcept {
    const amrex::Real decay = std::exp(-amrex::Real(2.0) * nu * k * k * t);
    u =  u0 * std::sin(k*x) * std::cos(k*y) * std::cos(k*z) * decay;
    v = -u0 * std::cos(k*x) * std::sin(k*y) * std::cos(k*z) * decay;
    w = amrex::Real(0.0);
}

inline amrex::Real exactRho(amrex::Real x, amrex::Real y, amrex::Real z,
                            amrex::Real k, amrex::Real u0, amrex::Real nu, amrex::Real t) noexcept {
    amrex::ignore_unused(x,y,z,k,u0,nu,t);
    return amrex::Real(1.0);
}

}} // namespace muphfasa::tgv3d

#endif
