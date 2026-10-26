// Tgv3dAnalytic.H — classic 3-D Taylor-Green vortex initial condition, plus
//
// Triply-periodic domain, one vortex wavelength per domain (k=2*pi/N, N
// cells per side). The t=0 field below is the classic 3-D TGV IC (Taylor &
// Green 1937; Brachet et al. 1983), w=0:
//
//   u(x,y,z,0) =  U0 * sin(k*x) * cos(k*y) * cos(k*z)
//   v(x,y,z,0) = -U0 * cos(k*x) * sin(k*y) * cos(k*z)
//   w(x,y,z,0) =  0
//
// IMPORTANT: unlike the 2-D TGV (benchmarks/tgv2d/TgvAnalytic.H) and the
// Beltrami flow (benchmarks/beltrami3d/BeltramiAnalytic.H), this IC is NOT
// an exact solution of 3-D Navier-Stokes for t>0. Both of those cases have
// vorticity parallel/perpendicular in a way that makes the nonlinear
// advection term (u.grad)u a pure gradient (absorbed into pressure), so a
// single Fourier mode with simple exponential viscous decay solves NS
// exactly for all t. The classic 3-D TGV has no such property: vortex
// stretching is nonzero (curl(u) is not aligned with u, and is not
// perpendicular to the domain like the 2-D case), so energy cascades into
// other modes -- this is precisely *why* 3-D TGV is the canonical DNS test
// case for transition to turbulence (Brachet et al. 1983), and no closed
// form exists for t>0.
//
// initUVW() below is only the leading-order/linearized (Stokes) viscous
// decay of the t=0 mode, exp(-2*nu*k^2*t) applied to the unchanged spatial
// pattern -- it ignores the nonlinear mode coupling entirely, and is only
// ever evaluated at t=0 (by ../InitData.cpp), where it's exact.
// It is NOT used for t>0 validation: validate/test_tgv3d.cpp instead checks
// the global kinetic-energy balance KE(0)-KE(T) = nu*integral(|curl(u)|^2)dt,
// an identity that holds exactly for incompressible Navier-Stokes at any
// Re/amplitude/duration despite the vortex stretching described above (see
// validate/README.md and the comment in validate/test_tgv3d.cpp).
//
// rho is held at 1 (weak-compressibility assumption, Ma<<1), as in
// BeltramiAnalytic.H, rather than attempting the 2-D TGV's second-order
// density correction.
//
// Used by InitData.cpp (t=0)
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

inline void initUVW(amrex::Real x, amrex::Real y, amrex::Real z,
                     amrex::Real k, amrex::Real u0, amrex::Real nu, amrex::Real t,
                     amrex::Real& u, amrex::Real& v, amrex::Real& w) noexcept {
    const amrex::Real decay = std::exp(-amrex::Real(2.0) * nu * k * k * t);
    u =  u0 * std::sin(k*x) * std::cos(k*y) * std::cos(k*z) * decay;
    v = -u0 * std::cos(k*x) * std::sin(k*y) * std::cos(k*z) * decay;
    w = amrex::Real(0.0);
}

inline amrex::Real intRho(amrex::Real x, amrex::Real y, amrex::Real z,
                            amrex::Real k, amrex::Real u0, amrex::Real nu, amrex::Real t) noexcept {
    amrex::ignore_unused(x,y,z,k,u0,nu,t);
    return amrex::Real(1.0);
}

}} // namespace muphfasa::tgv3d

#endif
