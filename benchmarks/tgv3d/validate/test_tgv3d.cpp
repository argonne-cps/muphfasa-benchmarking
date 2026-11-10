// test_tgv3d.cpp — TGV Re=1600 validation vs. Dairay et al. 2017 DNS (TGV_Re1600.dat).
//
// Triply-periodic 3-D Taylor-Green vortex at Re=1600. The classic TGV has no
// closed-form t>0 solution (vortex stretching cascades energy), so two
// complementary checks are performed:
//
//  1. Global energy balance (in-situ, always):  KE(0)-KE(T) = nu * int |curl u|^2 dt
//     Exact for periodic incompressible NS; residual is trapz + truncation error.
//  2. DNS benchmark (post-processing): E_k(t), dissipation epsilon(t), and
//     enstrophy vs. TGV_Re1600.dat, evaluated by analyze_tgv.py reading the
//     AMReX plotfiles via yt. This driver writes the plotfiles and an in-situ
//     CSV (tgv_diagnostics.csv) with the same quantities for cross-check.
//
// Physical units: L_ref=1 (one wavelength = domain length), U0=1, Re=1600
// => nu=6.25e-4. Domain [0,1]^3 (equivalent to benchmark (2pi)^3 up to L scale).
// Time t=20 in L/U units (benchmark range). The inputs file pins lbm.u_lat
// so tau is derived consistently; see inputs for the dt/step budget.

#include "IO/InputParser.H"
#include "Simulation/ISimulation.H"
#include "Simulation/SimulationFactory.H"
#include "Tgv3dAnalytic.H"

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cmath>
#include <fstream>
#include <vector>

using namespace muphfasa;

namespace {

struct Diagnostics {
    amrex::Real time_phys = 0.0;
    int step = 0;
    amrex::Real ke = 0.0;          // 0.5 * mean(|u_phys|^2)  volume-averaged
    amrex::Real enstrophy = 0.0;   // 0.5 * mean(|vort_phys|^2)
    amrex::Real dissipation = 0.0; // nu * mean(|grad u_phys|^2) = 2*nu*enstrophy (periodic)
    amrex::Real dissipation_alt = 0.0; // -dEk/dt via finite diff (check)
};

void computeDiagnostics(const amrex::MultiFab& vel_lat, const amrex::Geometry& geom,
                        const UnitConverter& units, Diagnostics& d)
{
    const amrex::Real C_vel = units.dx() / units.dt(); // u_phys = u_lat * C_vel
    const amrex::Real C_vel2 = C_vel * C_vel;
    const amrex::Real dx_phys = units.dx();
    const amrex::Real nu_lat = (units.relaxationTime() - amrex::Real(0.5)) * tgv3d::kCs2;
    const amrex::Real nu_phys = units.toPhysicalViscosity(nu_lat);

    const long n_cells = vel_lat.boxArray().numPts();
    const amrex::Real inv_n = amrex::Real(1.0) / amrex::Real(n_cells);

    amrex::Real ke_lat_sum = amrex::Real(0.5) * amrex::MultiFab::Dot(vel_lat, 0, vel_lat.nComp(), 0);
    d.ke = ke_lat_sum * inv_n * C_vel2;

    amrex::MultiFab vel_g(vel_lat.boxArray(), vel_lat.DistributionMap(), vel_lat.nComp(), 1);
    amrex::MultiFab::Copy(vel_g, vel_lat, 0, 0, vel_lat.nComp(), 0);
    vel_g.FillBoundary(geom.periodicity());

    amrex::MultiFab grad_sq(vel_lat.boxArray(), vel_lat.DistributionMap(), 1, 0);
    amrex::MultiFab vort_sq(vel_lat.boxArray(), vel_lat.DistributionMap(), 1, 0);

    for (amrex::MFIter mfi(vel_g); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        amrex::Array4<amrex::Real const> const& u = vel_g.const_array(mfi);
        amrex::Array4<amrex::Real> const& g2 = grad_sq.array(mfi);
        amrex::Array4<amrex::Real> const& v2 = vort_sq.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const amrex::Real h = amrex::Real(0.5);
            const amrex::Real dudx = h * (u(i+1,j,k,0) - u(i-1,j,k,0));
            const amrex::Real dudy = h * (u(i,j+1,k,0) - u(i,j-1,k,0));
            const amrex::Real dudz = h * (u(i,j,k+1,0) - u(i,j,k-1,0));
            const amrex::Real dvdx = h * (u(i+1,j,k,1) - u(i-1,j,k,1));
            const amrex::Real dvdy = h * (u(i,j+1,k,1) - u(i,j-1,k,1));
            const amrex::Real dvdz = h * (u(i,j,k+1,1) - u(i,j,k-1,1));
            const amrex::Real dwdx = h * (u(i+1,j,k,2) - u(i-1,j,k,2));
            const amrex::Real dwdy = h * (u(i,j+1,k,2) - u(i,j-1,k,2));
            const amrex::Real dwdz = h * (u(i,j,k+1,2) - u(i,j,k-1,2));

            g2(i,j,k) = dudx*dudx + dudy*dudy + dudz*dudz
                      + dvdx*dvdx + dvdy*dvdy + dvdz*dvdz
                      + dwdx*dwdx + dwdy*dwdy + dwdz*dwdz;

            const amrex::Real wx = dwdy - dvdz;
            const amrex::Real wy = dudz - dwdx;
            const amrex::Real wz = dvdx - dudy;
            v2(i,j,k) = wx*wx + wy*wy + wz*wz;
        });
    }

    // grad_phys = grad_lat / dt  (since grad_phys = grad_lat*C_vel/dx_phys = grad_lat/dt)
    // vort_phys = vort_lat / dt
    const amrex::Real dt = units.dt();
    const amrex::Real inv_dt2 = amrex::Real(1.0) / (dt * dt);
    amrex::Real g2_sum = grad_sq.sum(0);
    amrex::Real v2_sum = vort_sq.sum(0);

    amrex::Real g2_phys_mean = g2_sum * inv_n * inv_dt2;
    amrex::Real v2_phys_mean = v2_sum * inv_n * inv_dt2;

    d.dissipation = nu_phys * g2_phys_mean;
    d.enstrophy = amrex::Real(0.5) * v2_phys_mean;
}

} // namespace

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        InputParser params;
        std::unique_ptr<ISimulation> sim = makeSimulation(params);
        sim->initialize();

        const UnitConverter& units = sim->units();
        const amrex::Real tau = units.relaxationTime();
        const amrex::Real nu_lat = (tau - amrex::Real(0.5)) * tgv3d::kCs2;
        const amrex::Real nu_phys = units.toPhysicalViscosity(nu_lat);
        const amrex::Geometry& geom = sim->geom();

        int max_step = params.maxStep();
        const amrex::Real stop_time_phys = params.stopTime();
        if (stop_time_phys > 0.0) {
            int steps_from_time = static_cast<int>(std::ceil(units.toLatticeTime(stop_time_phys)));
            if (max_step <= 0 || steps_from_time < max_step) {
                max_step = steps_from_time;
                amrex::Print() << "Using stop_time=" << stop_time_phys
                               << " => max_step=" << max_step
                               << " (dt=" << units.dt() << " dx=" << units.dx() << ")\n";
            }
        }

        int plot_int = params.plotInterval();
        int sample_int=-1; 
        amrex::ParmParse ppbench("bench");
        ppbench.query("sample_int", sample_int);

        if (plot_int <= 0) plot_int = max_step;

        amrex::Print() << "TGV Re=1600 validation: N=" << params.ncells()[0]
                       << " Re=" << params.reynoldsNumber()
                       << " U0=" << params.refVelocity()
                       << " L=" << params.refLength()
                       << " nu_phys=" << nu_phys
                       << " tau=" << tau
                       << " u_lat=" << units.toLatticeVelocity(params.refVelocity())
                       << " Ma=" << units.machNumber()
                       << " dx=" << units.dx() << " dt=" << units.dt()
                       << " max_step=" << max_step
                       << " plot_int=" << plot_int << "\n";

        std::vector<Diagnostics> hist;
        hist.reserve(max_step / std::max(1, plot_int) + 2);

        const amrex::Real L_ref = params.refLength();
        const amrex::Real U0 = params.refVelocity();

        std::ofstream csv("tgv_diagnostics.csv", std::ios::out);
        csv << "step,time_phys,time,Ek,dissipation,enstrophy\n";

        auto doSample = [&](int step) {
            sim->updateMacroscopic();
            Diagnostics d;
            d.step = step;
            d.time_phys = units.toPhysicalTime(sim->time());
            computeDiagnostics(sim->vel(), geom, units, d);
            hist.push_back(d);
            const amrex::Real t_bench = d.time_phys * U0 / L_ref;
            const amrex::Real ek_bench = d.ke / (U0 * U0);
            const amrex::Real eps_bench = d.dissipation * L_ref / (U0 * U0 * U0);
            const amrex::Real enst_bench = d.enstrophy * L_ref * L_ref / (U0 * U0);
            amrex::Print() << "sample step=" << step << " t_phys=" << d.time_phys
                           << " t_bench=" << t_bench
                           << " Ek=" << d.ke << " eps=" << d.dissipation
                           << " enstrophy=" << d.enstrophy << "\n";
            csv << d.step << "," << d.time_phys << "," << t_bench << ","
                << ek_bench << "," << eps_bench << "," << enst_bench << "\n";
            csv.flush();
        };

        // Initial sample at t=0 before any advance
        doSample(0);

        for (int step = 0; step <= max_step; ++step) {
            sim->advanceFused(step);
            sim->writeOutput(step);

            // Sample at plot interval (avoid duplicating the t=0 sample)
            if (step > 0 && (step % sample_int == 0)) {
                doSample(step);
            }
            if (stop_time_phys > 0.0 && units.toPhysicalTime(sim->time()) >= stop_time_phys) {
                if (step % plot_int != 0) doSample(step);
                break;
            }
        }
        csv.close();

        if (hist.size() < 2) {
            amrex::Print() << "WARNING: fewer than 2 diagnostic samples; "
                              "reduce amr.plot_int for a finer time history.\n";
        } else {
            amrex::Real diss_int = 0.0;
            for (std::size_t n = 1; n < hist.size(); ++n) {
                amrex::Real dt = hist[n].time_phys - hist[n-1].time_phys;
                diss_int += 0.5 * dt * (hist[n].dissipation + hist[n-1].dissipation);
            }
            amrex::Real ke_drop = hist.front().ke - hist.back().ke;
            amrex::Real denom = std::max(std::abs(ke_drop), amrex::Real(1e-30));
            amrex::Real rel_err = std::abs(ke_drop - diss_int) / denom;

            amrex::Print() << "\n--- Energy balance (physical units) ---\n"
                           << "KE(0)=" << hist.front().ke
                           << " KE(T)=" << hist.back().ke
                           << " KE(0)-KE(T)=" << ke_drop << "\n"
                           << "int eps dt=" << diss_int
                           << " rel_err=" << rel_err << "\n";

            const amrex::Real tol = 0.10;
            if (rel_err < tol) {
                amrex::Print() << "PASS energy balance (tol " << tol << ")\n";
            } else {
                amrex::Print() << "FAIL energy balance: rel_err " << rel_err
                               << " exceeds tol " << tol
                               << " (try smaller amr.plot_int or larger N)\n";
            }

            amrex::Print() << "\nDiagnostics written to tgv_diagnostics.csv ("
                           << hist.size() << " samples).\n"
                           << "For DNS comparison run:\n"
                           << "  python3 analyze_tgv.py --run-dir . --benchmark TGV_Re1600.dat\n";
        }
    }
    amrex::Finalize();
    return 0;
}
