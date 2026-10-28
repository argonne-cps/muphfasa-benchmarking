// test_tgv3d.cpp — Correctness test: 3-D TGV global kinetic-energy balance.
//
// Mirrors external/muphfasa/tests/regression/tgv2d/test_tgv2d.cpp's use of the
// normal InputParser/ISimulation/advanceFused() machinery with the non-quiescent
// initial condition supplied by ../InitData.cpp (see ../Tgv3dAnalytic.H for
// the analytic IC). This is the validation counterpart to
// benchmarks/tgv3d/inputs.bench (performance-only, no correctness check) --
// see benchmarks/tgv3d/../../README.md's "Validate vs Perf" section.
//
// Unlike tgv2d/beltrami3d, the classic 3-D TGV has no closed-form solution
// for t>0 (vortex stretching moves energy into other Fourier modes -- see
// ../Tgv3dAnalytic.H), so this test can't fcompare the field against an
// analytic plotfile. Instead it checks the one identity that IS exact for
// incompressible Navier-Stokes regardless of vortex stretching: the global
// kinetic-energy balance
//
//   d/dt [ (1/2) integral |u|^2 dV ] = -nu * integral |curl(u)|^2 dV
//
// Nonlinear advection and the pressure gradient do no net work on a
// periodic, divergence-free field -- they only move energy between Fourier
// modes, never create or destroy it -- so kinetic energy can only change via
// viscous dissipation, exactly, at any Re/amplitude/duration. Integrating
// over the run:
//
//   KE(0) - KE(T) = nu * integral_0^T ( integral |curl(u)|^2 dV ) dt
//
// At every amr.plot_int step this computes, in lattice units straight off
// sim->vel() (dx=dt=1 so no unit conversion needed for the balance itself):
//   KE(t) = (1/2) * sum_cells |u|^2
//   W(t)  =         sum_cells |curl(u)|^2   (central differences on a
//                                             periodic-ghost-filled copy)
// then compares KE(0)-KE(T) against nu * trapz(W, t) and reports the
// relative error. Remaining error here comes only from the trapezoidal time
// integration and the solver's own truncation error, not from an
// approximation to the physics -- unlike the old pointwise check against a
// linearized decay, which only held in a regime too weak to ever exercise
// the nonlinear term it was meant to be validating.

#include "IO/InputParser.H"
#include "Simulation/ISimulation.H"
#include "Simulation/SimulationFactory.H"
#include "Tgv3dAnalytic.H"

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cmath>
#include <vector>

using namespace muphfasa;

namespace {

// KE = (1/2) sum|u|^2 ; curl_sq = sum|curl(u)|^2, both in lattice units
// (dx=dt=1), via central differences on a periodic-ghost-filled copy of vel.
void computeEnergyAndCurlSq(const amrex::MultiFab& vel, const amrex::Geometry& geom,
                            amrex::Real& ke, amrex::Real& curl_sq)
{
    ke = amrex::Real(0.5) * amrex::MultiFab::Dot(vel, 0, vel.nComp(), 0);

    amrex::MultiFab vel_g(vel.boxArray(), vel.DistributionMap(), vel.nComp(), 1);
    amrex::MultiFab::Copy(vel_g, vel, 0, 0, vel.nComp(), 0);
    vel_g.FillBoundary(geom.periodicity());

    amrex::MultiFab curl(vel.boxArray(), vel.DistributionMap(), 3, 0);
    for (amrex::MFIter mfi(vel_g); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        amrex::Array4<amrex::Real const> const& u = vel_g.const_array(mfi);
        amrex::Array4<amrex::Real> const& w = curl.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const amrex::Real dudy = (u(i,j+1,k,0) - u(i,j-1,k,0)) * amrex::Real(0.5);
            const amrex::Real dudz = (u(i,j,k+1,0) - u(i,j,k-1,0)) * amrex::Real(0.5);
            const amrex::Real dvdx = (u(i+1,j,k,1) - u(i-1,j,k,1)) * amrex::Real(0.5);
            const amrex::Real dvdz = (u(i,j,k+1,1) - u(i,j,k-1,1)) * amrex::Real(0.5);
            const amrex::Real dwdx = (u(i+1,j,k,2) - u(i-1,j,k,2)) * amrex::Real(0.5);
            const amrex::Real dwdy = (u(i,j+1,k,2) - u(i,j-1,k,2)) * amrex::Real(0.5);
            w(i,j,k,0) = dwdy - dvdz;
            w(i,j,k,1) = dudz - dwdx;
            w(i,j,k,2) = dvdx - dudy;
        });
    }
    curl_sq = amrex::MultiFab::Dot(curl, 0, curl.nComp(), 0);
}

} // namespace

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        InputParser params;
        std::unique_ptr<ISimulation> sim = makeSimulation(params);
        sim->initialize();

        const amrex::Real tau = sim->units().relaxationTime();
        const amrex::Real nu  = (tau - amrex::Real(0.5)) * tgv3d::kCs2;

        const int max_step = params.maxStep();
        const amrex::Geometry& geom = sim->geom();

        std::vector<amrex::Real> times, ke_hist, curlsq_hist;

        for (int step = 0; step <= max_step; ++step) {
            sim->advanceFused(step);
            sim->writeOutput(step);

            if (step % params.plotInterval() == 0) {
                amrex::Real ke, curl_sq;
                computeEnergyAndCurlSq(sim->vel(), geom, ke, curl_sq);
                times.push_back(sim->time());
                ke_hist.push_back(ke);
                curlsq_hist.push_back(curl_sq);
                amrex::Print() << "At Time=" << sim->time()
                               << ": KE=" << ke << " |curl(u)|^2=" << curl_sq << "\n";
            }
        }

        // trapezoidal integral of nu*|curl(u)|^2 over the recorded times
        amrex::Real dissipation_integral = amrex::Real(0.0);
        for (std::size_t n = 1; n < times.size(); ++n) {
            const amrex::Real dt = times[n] - times[n - 1];
            dissipation_integral += amrex::Real(0.5) * dt * nu * (curlsq_hist[n] + curlsq_hist[n - 1]);
        }

        const amrex::Real ke_drop = ke_hist.front() - ke_hist.back();
        const amrex::Real denom = std::max(std::abs(ke_drop), amrex::Real(1e-30));
        const amrex::Real rel_err = std::abs(ke_drop - dissipation_integral) / denom;

        amrex::Print() << "KE(0)=" << ke_hist.front() << " KE(T)=" << ke_hist.back()
                       << " dissipated=" << ke_drop
                       << " nu*integral(|curl(u)|^2)dt=" << dissipation_integral
                       << " relative error=" << rel_err << "\n";
    }
    amrex::Finalize();
    return 0;
}
