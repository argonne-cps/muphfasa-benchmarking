// InitData.cpp — 3-D Taylor-Green vortex initial condition (t=0).
// Overrides default quiescent src/Simulation/InitData.cpp via VPATH.

#include "InitData.H"
#include "Tgv3dAnalytic.H"
#include "IO/InputParser.H"
#include "LBM/Units.H"

#include <AMReX_MultiFab.H>

namespace muphfasa {

void initData(amrex::MultiFab& rho, amrex::MultiFab& vel, const UnitConverter& units) {
    InputParser params;
    amrex::Real u0 = units.toLatticeVelocity(params.refVelocity());
    const amrex::Box domain = rho.boxArray().minimalBox();
    const int N = domain.length(0); // cubic
    const amrex::Real k = tgv3d::waveNumber(N);
    for (amrex::MFIter mfi(rho); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto r = rho.array(mfi);
        auto u = vel.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k_) {
            amrex::Real x = amrex::Real(i) + 0.5;
            amrex::Real y = amrex::Real(j) + 0.5;
            amrex::Real z = amrex::Real(k_) + 0.5;
            amrex::Real ux, uy, uz;
            tgv3d::exactUVW(x,y,z,k,u0,amrex::Real(0),amrex::Real(0),ux,uy,uz);
            r(i,j,k_,0) = tgv3d::exactRho(x,y,z,k,u0,amrex::Real(0),amrex::Real(0));
            u(i,j,k_,0) = ux;
            u(i,j,k_,1) = uy;
            u(i,j,k_,2) = uz;
        });
    }
}

} // namespace muphfasa
