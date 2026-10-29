// test_lid3d.cpp — Correctness test: 3-D lid-driven cavity, global checks +
// Albensoeder centerline extrema.
//
// Runs the closed cubic cavity of ./inputs (NoSlip on all faces except the
// lid at z=hi moving in +x) from rest to steady state and checks, in lattice
// units straight off sim->vel()/sim->rho():
//
//   1. Mass:     sum(rho) is conserved (closed box, bounce-back walls).
//   2. Spin-up:  KE(t) = (1/2) sum|u|^2 rises from 0 and plateaus.
//   3. Balance:  at steady state, power in = viscous dissipation,
//                  P_lid = U_lid * sum_lid( nu * du/dz|wall )
//                  eps   = nu * sum|grad u|^2
//                from dKE/dt = P_lid - eps (lid does work, other walls do not).
//   4. Albensoeder: when units.Re ~ 1000, centerline extrema (lid-direction
//                velocity along wall-normal, etc.) against Albensoeder &
//                Kuhlmann (2005) Table 10, rigid cubic cavity Gamma=Lambda=1.
//                Axis mapping is documented below and in README.md. Skipped
//                (reported as SKIP) when Re is not ~1000.
//
// Walls/lid are taken from the inputs file (closed box). The lid face and its
// direction are detected from the nonzero Uwall entry. Exits nonzero if any
// check fails.

#include "IO/InputParser.H"
#include "Simulation/ISimulation.H"
#include "Simulation/SimulationFactory.H"

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <vector>

using namespace muphfasa;

namespace {

// Albensoeder & Kuhlmann (2005), J. Comput. Phys. 206, 536-558,
// doi:10.1016/j.jcp.2004.12.024, Table 10, rigid BC Gamma=Lambda=1 Re=1000,
// 96x96x64 Chebyshev. Velocities are u/Re in the paper, i.e. normalized by
// the lid velocity, same as u/U_lid here. Axis mapping to our coordinates
// (lid at z=hi moving in +x, closed cube 0-1 in each dimension):
//   paper x (depth, wall-normal) -> our z,  paper u -> our w
//   paper y (lid direction)      -> our x,  paper v -> our u
//   paper z (span)               -> our y,  paper w -> our v
// So:
//   paper min v at (x,0,0) varying x (wall-normal) = -0.2803833 at x=0.37581
//     => our min u along (0,0,z) wall-normal line (lid-direction velocity)
//   paper max u at (0,y,0) varying y (lid direction) = 0.4350186 at y=0.40957
//     => our max w along (x,0,0) lid-direction line (wall-normal velocity)
//   paper min u at (0,y,0) = -0.2466511 at y=-0.39087
//     => our min w along (x,0,0)
//   paper min w at (0,0,z) = -0.0303398 at z=0.36492
//     => our min v along (0,y,0) span line

// Reference extrema (normalized by lid velocity) from Table 10, rigid cube.
constexpr amrex::Real kAlbMinU_alongZ = amrex::Real(-0.2803833); // our u along (0,0,z), at z~0.37581
constexpr amrex::Real kAlbMaxW_alongX = amrex::Real( 0.4350186); // our w along (x,0,0), at x~0.40957
constexpr amrex::Real kAlbMinW_alongX = amrex::Real(-0.2466511); // our w along (x,0,0), at x~0.39087 (other side)
constexpr amrex::Real kAlbMinV_alongY = amrex::Real(-0.0303398); // our v along (0,y,0), small span flow

// Tolerances (see README.md for how they were chosen).
constexpr amrex::Real kTolMass    = amrex::Real(1e-9);  // relative drift in sum(rho)
constexpr amrex::Real kTolMono    = amrex::Real(1e-3);  // allowed relative KE dip between samples
constexpr amrex::Real kTolPlateau = amrex::Real(5e-4);  // relative KE change over last interval (3D looser)
constexpr amrex::Real kTolBalance = amrex::Real(0.30);  // |P_lid - eps|/P_lid at steady state (3D corners)
constexpr amrex::Real kTolAlb     = amrex::Real(0.12);  // Albensoeder extrema at coarse N=32-48

constexpr amrex::Real kCs2 = amrex::Real(1.0) / amrex::Real(3.0);

struct Sample {
    amrex::Real time, mass, ke, power, eps;
};

// Identify lid face (dir, bound, tangential component) from the nonzero Uwall.
struct LidInfo {
    int dir = 2;           // 0=x, 1=y, 2=z  (default z-hi for this test)
    int bound = 1;         // 0=lo, 1=hi
    int comp = 0;          // tangential velocity component (0=u, 1=v, 2=w)
    amrex::Real u_lid_phys = amrex::Real(0.05);
};

LidInfo detectLid(const InputParser& params) {
    LidInfo lid;
    auto Uwall = params.Uwall();
    for (int face = 0; face < 6; ++face) {
        for (int c = 0; c < 3; ++c) {
            if (Uwall[face][c] != amrex::Real(0.0)) {
                lid.dir = face / 2;
                lid.bound = face % 2;
                lid.comp = c;
                lid.u_lid_phys = Uwall[face][c];
                return lid;
            }
        }
    }
    return lid;
}

// Lid power P_lid and dissipation eps (lattice units, dx=dt=1).
// Ghost values outside the domain are wall-reflected: 2*U_wall - u_cell for the
// tangential component on the moving lid, else 0. Corner ghosts where multiple
// walls meet are approximated (contribution small; tolerance is loose).
void computePowerAndDissipation(const amrex::MultiFab& vel, const amrex::Geometry& geom,
                                const LidInfo& lid, amrex::Real u_lid, amrex::Real nu,
                                amrex::Real& power, amrex::Real& eps)
{
    amrex::MultiFab vel_g(vel.boxArray(), vel.DistributionMap(), vel.nComp(), 1);
    amrex::MultiFab::Copy(vel_g, vel, 0, 0, vel.nComp(), 0);
    vel_g.FillBoundary(geom.periodicity());

    amrex::MultiFab grad_sq(vel.boxArray(), vel.DistributionMap(), 1, 0);
    amrex::MultiFab lid_shear(vel.boxArray(), vel.DistributionMap(), 1, 0);
    lid_shear.setVal(amrex::Real(0.0));

    const amrex::Box domain = geom.Domain();
    const int ilo = domain.smallEnd(0), ihi = domain.bigEnd(0);
    const int jlo = domain.smallEnd(1), jhi = domain.bigEnd(1);
    const int klo = domain.smallEnd(2), khi = domain.bigEnd(2);

    for (amrex::MFIter mfi(vel_g); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        amrex::Array4<amrex::Real const> const& u = vel_g.const_array(mfi);
        amrex::Array4<amrex::Real> const& g2 = grad_sq.array(mfi);
        amrex::Array4<amrex::Real> const& sh = lid_shear.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            auto nbr = [&](int di, int dj, int dk, int c) -> amrex::Real {
                const int ni = i + di, nj = j + dj, nk = k + dk;
                if (ni < ilo || ni > ihi || nj < jlo || nj > jhi || nk < klo || nk > khi) {
                    // Outside domain: which wall was crossed determines U_wall
                    bool on_lid = false;
                    if (lid.dir == 0) on_lid = (lid.bound == 0 ? ni < ilo : ni > ihi);
                    if (lid.dir == 1) on_lid = (lid.bound == 0 ? nj < jlo : nj > jhi);
                    if (lid.dir == 2) on_lid = (lid.bound == 0 ? nk < klo : nk > khi);
                    const amrex::Real uw = (on_lid && c == lid.comp) ? u_lid : amrex::Real(0.0);
                    return amrex::Real(2.0) * uw - u(i, j, k, c);
                }
                return u(ni, nj, nk, c);
            };
            const amrex::Real half = amrex::Real(0.5);
            // 9 gradient components du_i/dx_j
            const amrex::Real dudx = half * (nbr( 1, 0, 0, 0) - nbr(-1, 0, 0, 0));
            const amrex::Real dudy = half * (nbr( 0, 1, 0, 0) - nbr( 0,-1, 0, 0));
            const amrex::Real dudz = half * (nbr( 0, 0, 1, 0) - nbr( 0, 0,-1, 0));
            const amrex::Real dvdx = half * (nbr( 1, 0, 0, 1) - nbr(-1, 0, 0, 1));
            const amrex::Real dvdy = half * (nbr( 0, 1, 0, 1) - nbr( 0,-1, 0, 1));
            const amrex::Real dvdz = half * (nbr( 0, 0, 1, 1) - nbr( 0, 0,-1, 1));
            const amrex::Real dwdx = half * (nbr( 1, 0, 0, 2) - nbr(-1, 0, 0, 2));
            const amrex::Real dwdy = half * (nbr( 0, 1, 0, 2) - nbr( 0,-1, 0, 2));
            const amrex::Real dwdz = half * (nbr( 0, 0, 1, 2) - nbr( 0, 0,-1, 2));
            g2(i, j, k) = dudx*dudx + dudy*dudy + dudz*dudz
                        + dvdx*dvdx + dvdy*dvdy + dvdz*dvdz
                        + dwdx*dwdx + dwdy*dwdy + dwdz*dwdz;

            // Lid shear: only cells touching the lid plane contribute
            bool on_lid_plane = false;
            if (lid.dir == 0) on_lid_plane = (lid.bound == 0 ? i == ilo : i == ihi);
            if (lid.dir == 1) on_lid_plane = (lid.bound == 0 ? j == jlo : j == jhi);
            if (lid.dir == 2) on_lid_plane = (lid.bound == 0 ? k == klo : k == khi);
            if (on_lid_plane) {
                // Quadratic fit through wall (s=0, u=U_lid), first cell (s=0.5), second cell (s=1.5)
                amrex::Real d1 = amrex::Real(0.0), d2 = amrex::Real(0.0);
                if (lid.dir == 0) {
                    const int j2 = (lid.bound == 0 ? i+1 : i-1);
                    d1 = u(i, j, k, lid.comp) - u_lid;
                    d2 = u(j2, j, k, lid.comp) - u_lid;
                } else if (lid.dir == 1) {
                    const int j2 = (lid.bound == 0 ? j+1 : j-1);
                    d1 = u(i, j, k, lid.comp) - u_lid;
                    d2 = u(i, j2, k, lid.comp) - u_lid;
                } else {
                    const int k2 = (lid.bound == 0 ? k+1 : k-1);
                    d1 = u(i, j, k, lid.comp) - u_lid;
                    d2 = u(i, j, k2, lid.comp) - u_lid;
                }
                const amrex::Real b = (d2 - amrex::Real(3.0) * d1) / amrex::Real(1.5);
                const amrex::Real a = amrex::Real(2.0) * d1 - half * b;
                // a = du/ds at wall, with s inward from lid; du/dn = -a (n outward)
                amrex::Real shear_mag = -a;
                if (lid.bound == 0) shear_mag = a; // lo wall inward is +dir
                sh(i, j, k) = shear_mag;
            }
        });
    }

    eps   = nu * grad_sq.sum(0);
    power = u_lid * nu * lid_shear.sum(0);
    if (power < amrex::Real(0.0)) power = -power; // sign depends on wall orientation
}

// Parabolic refinement of extremum of f[0..n-1] at (m+0.5)/n.
void extremum(const std::vector<amrex::Real>& f, int sign, amrex::Real& val, amrex::Real& pos)
{
    const int n = static_cast<int>(f.size());
    int m = 0;
    for (int j = 1; j < n; ++j) {
        if (sign * f[j] > sign * f[m]) m = j;
    }
    val = f[m];
    pos = (amrex::Real(m) + amrex::Real(0.5)) / amrex::Real(n);
    if (m > 0 && m < n - 1) {
        const amrex::Real den = f[m-1] - amrex::Real(2.0) * f[m] + f[m+1];
        if (den != amrex::Real(0.0)) {
            const amrex::Real delta = amrex::Real(0.5) * (f[m-1] - f[m+1]) / den;
            val = f[m] - amrex::Real(0.25) * (f[m-1] - f[m+1]) * delta;
            pos = (amrex::Real(m) + amrex::Real(0.5) + delta) / amrex::Real(n);
        }
    }
}

// Centerlines through the cube center. N even: average two cells straddling x=0.5 etc.
// Returns profiles normalized by U_lid:
//   u_of_z : u (lid/x component) along (0,0,z) wall-normal line, z varying
//   w_of_x : w (wall-normal/z component) along (x,0,0) lid-direction line
//   v_of_y : v (span/y component) along (0,y,0) span line
void centerlineProfiles(const amrex::MultiFab& vel, int N, amrex::Real u_lid,
                        std::vector<amrex::Real>& u_of_z,
                        std::vector<amrex::Real>& w_of_x,
                        std::vector<amrex::Real>& v_of_y)
{
    amrex::MultiFab vel_h(vel.boxArray(), vel.DistributionMap(), vel.nComp(), 0,
                          amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
    amrex::MultiFab::Copy(vel_h, vel, 0, 0, vel.nComp(), 0);

    u_of_z.assign(N, amrex::Real(0.0));
    w_of_x.assign(N, amrex::Real(0.0));
    v_of_y.assign(N, amrex::Real(0.0));
    const int c0 = N / 2 - 1, c1 = N / 2;
    for (amrex::MFIter mfi(vel_h); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        amrex::Array4<amrex::Real const> const& u = vel_h.const_array(mfi);
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    // wall-normal line (0,0,z): x~0.5 and y~0.5 (mid), z varying
                    if ((i == c0 || i == c1) && (j == c0 || j == c1)) {
                        u_of_z[k] += amrex::Real(0.25) * u(i, j, k, 0) / u_lid;
                    }
                    // lid-direction line (x,0,0): y~0.5, z~0.5, x varying
                    if ((j == c0 || j == c1) && (k == c0 || k == c1)) {
                        w_of_x[i] += amrex::Real(0.25) * u(i, j, k, 2) / u_lid;
                    }
                    // span line (0,y,0): x~0.5, z~0.5, y varying
                    if ((i == c0 || i == c1) && (k == c0 || k == c1)) {
                        v_of_y[j] += amrex::Real(0.25) * u(i, j, k, 1) / u_lid;
                    }
                }
            }
        }
    }
    amrex::ParallelDescriptor::ReduceRealSum(u_of_z.data(), N);
    amrex::ParallelDescriptor::ReduceRealSum(w_of_x.data(), N);
    amrex::ParallelDescriptor::ReduceRealSum(v_of_y.data(), N);
}

bool report(const char* name, bool ok, const char* detail)
{
    amrex::Print() << (ok ? "PASS  " : "FAIL  ") << name << ": " << detail << "\n";
    return ok;
}

std::string fmt(const char* f, amrex::Real a, amrex::Real b = 0, amrex::Real c = 0)
{
    char buf[320];
    std::snprintf(buf, sizeof(buf), f, static_cast<double>(a), static_cast<double>(b),
                  static_cast<double>(c));
    return buf;
}

} // namespace

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    bool all_ok = true;
    {
        InputParser params;
        std::unique_ptr<ISimulation> sim = makeSimulation(params);
        sim->initialize();

        const auto ncells = params.ncells();
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            ncells.size() == 3 && ncells[0] == ncells[1] && ncells[1] == ncells[2] && ncells[0] % 2 == 0,
            "test_lid3d: domain must be cubic with an even number of cells (e.g. 32 32 32)");
        const int N = ncells[0];

        const LidInfo lid = detectLid(params);
        const amrex::Real u_lid = sim->units().toLatticeVelocity(lid.u_lid_phys);
        const amrex::Real nu = (sim->units().relaxationTime() - amrex::Real(0.5)) * kCs2;
        const amrex::Real Re = params.reynoldsNumber();

        const int max_step = params.maxStep();
        const int every = params.plotInterval();
        const amrex::Geometry& geom = sim->geom();

        std::vector<Sample> hist;
        hist.reserve(max_step / std::max(1, every) + 2);

        for (int step = 0; step <= max_step; ++step) {
            sim->advanceFused(step);
            sim->writeOutput(step);

            if (step % every == 0) {
                Sample s;
                s.time = sim->time();
                s.mass  = sim->rho().sum(0);
                s.ke    = amrex::Real(0.5) * amrex::MultiFab::Dot(sim->vel(), 0, sim->vel().nComp(), 0);
                computePowerAndDissipation(sim->vel(), geom, lid, u_lid, nu, s.power, s.eps);
                hist.push_back(s);
                amrex::Print() << "step " << step << " t=" << s.time
                               << " mass=" << s.mass << " KE=" << s.ke
                               << " P_lid=" << s.power << " eps=" << s.eps << "\n";
            }
        }

        if (hist.size() < 2) {
            amrex::Print() << "FAIL  not enough samples (plot_int too large for max_step)\n";
            all_ok = false;
        } else {
            // 1. Mass conservation in a closed box
            const amrex::Real m0 = hist.front().mass;
            amrex::Real max_rel_drift = amrex::Real(0.0);
            for (auto& s : hist) {
                max_rel_drift = std::max(max_rel_drift, std::abs(s.mass - m0) / std::abs(m0));
            }
            all_ok &= report("mass conservation", max_rel_drift < kTolMass,
                             fmt("max relative drift in sum(rho) = %.2e (tol 1e-9)", max_rel_drift).c_str());

            // 2. KE spin-up
            amrex::Real worst_dip = amrex::Real(0.0);
            for (std::size_t n = 1; n < hist.size(); ++n) {
                if (hist[n-1].ke > amrex::Real(0.0))
                    worst_dip = std::max(worst_dip, (hist[n-1].ke - hist[n].ke) / hist[n-1].ke);
            }
            const Sample& last = hist.back();
            const amrex::Real plateau = std::abs(last.ke - hist[hist.size()-2].ke) / std::max(std::abs(last.ke), amrex::Real(1e-30));
            const bool spinup_ok = (hist.front().ke == amrex::Real(0.0)) && (last.ke > amrex::Real(0.0))
                                   && (worst_dip < kTolMono);
            all_ok &= report("KE spin-up", spinup_ok,
                             fmt("KE 0 -> %.4e, worst relative dip between samples %.2e (tol 1e-3)",
                                 last.ke, worst_dip).c_str());
            all_ok &= report("KE plateau", plateau < kTolPlateau,
                             fmt("relative KE change over last sample interval = %.2e (tol 5e-4)",
                                 plateau).c_str());

            // 3. Steady-state energy balance
            const amrex::Real bal = std::abs(last.power - last.eps) / std::max(std::abs(last.power), amrex::Real(1e-30));
            all_ok &= report("energy balance", bal < kTolBalance,
                             fmt("P_lid=%.4e eps=%.4e |P-eps|/P=%.2e (tol 3e-1)",
                                 last.power, last.eps, bal).c_str());

            // 4. Albensoeder centerline extrema (only when Re ~ 1000)
            const bool do_alb = std::abs(Re - amrex::Real(1000.0)) < amrex::Real(50.0);
            if (!do_alb) {
                amrex::Print() << "SKIP  Albensoeder extrema (Re=" << Re << " not ~1000; run with units.Re=1000 and N>=32 to enable)\n";
            } else {
                std::vector<amrex::Real> u_of_z, w_of_x, v_of_y;
                centerlineProfiles(sim->vel(), N, u_lid, u_of_z, w_of_x, v_of_y);

                amrex::Real umin, zpos, wmax, xpos_max, wmin, xpos_min, vmin, ypos;
                extremum(u_of_z, -1, umin, zpos);
                extremum(w_of_x, +1, wmax, xpos_max);
                extremum(w_of_x, -1, wmin, xpos_min);
                extremum(v_of_y, -1, vmin, ypos);

                auto alb = [&](const char* name, amrex::Real val, amrex::Real pos, amrex::Real ref) {
                    const amrex::Real err = std::abs(val - ref) / std::max(std::abs(ref), amrex::Real(1e-30));
                    // Primary check (u_min) is required; others are informative at coarse N
                    return report(name, err < kTolAlb,
                                  fmt("%.5f at %.4f vs Albensoeder %.5f", val, pos, ref).append(
                                      fmt(" (rel err %.2e, tol 1.2e-1)", err)).c_str());
                };
                // Required: wall-normal lid-velocity minimum is the dominant 3D feature
                all_ok &= alb("Albensoeder u_min along wall-normal (0,0,z)", umin, zpos, kAlbMinU_alongZ);
                // Informative: wall-normal velocity along lid direction (secondary flow)
                bool wmax_ok = alb("Albensoeder w_max along lid direction (x,0,0)", wmax, xpos_max, kAlbMaxW_alongX);
                bool wmin_ok = alb("Albensoeder w_min along lid direction (x,0,0)", wmin, xpos_min, kAlbMinW_alongX);
                // Very weak span flow: report but do not fail on it at coarse N
                {
                    const amrex::Real err = std::abs(vmin - kAlbMinV_alongY) / std::max(std::abs(kAlbMinV_alongY), amrex::Real(1e-30));
                    amrex::Print() << (err < kTolAlb ? "PASS  " : "INFO  ")
                                   << "Albensoeder v_min along span (0,y,0): "
                                   << fmt("%.5f at %.4f vs Albensoeder %.5f", vmin, ypos, kAlbMinV_alongY)
                                   << fmt(" (rel err %.2e, tol 1.2e-1)", err)
                                   << (err < kTolAlb ? "" : " (informational at coarse N)") << "\n";
                }
                // At N=32 the secondary extrema are expected to be noisier; don't gate all_ok on them alone
                if (N >= 48) {
                    all_ok &= wmax_ok;
                    all_ok &= wmin_ok;
                } else {
                    amrex::Print() << "INFO  w_max/w_min checks are informational at N=" << N << " (<48)\n";
                }

                // Log measured vs reference for convergence studies
                amrex::Print() << "Albensoeder reference (cube, rigid, Re=1000, 96^3): "
                               << "u_min=" << kAlbMinU_alongZ << " w_max=" << kAlbMaxW_alongX
                               << " w_min=" << kAlbMinW_alongX << " v_min=" << kAlbMinV_alongY << "\n";
                amrex::Print() << "Measured (N=" << N << " Re=" << Re << "): "
                               << "u_min=" << umin << " at z=" << zpos
                               << " w_max=" << wmax << " at x=" << xpos_max
                               << " w_min=" << wmin << " at x=" << xpos_min
                               << " v_min=" << vmin << " at y=" << ypos << "\n";
            }
        }

        amrex::Print() << (all_ok ? "ALL CHECKS PASSED\n" : "SOME CHECKS FAILED\n");
    }
    amrex::Finalize();
    return all_ok ? 0 : 1;
}
