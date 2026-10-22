// bench.cpp — Generic LBM benchmark driver (no I/O, timed).
//
// Thin clone of external/muphfasa/tests/regression/couette/test_couette.cpp
// and tgv2d/test_tgv2d.cpp: uses InputParser/ISimulation/advance() with the
// benchmark's own InitData.cpp (if any) overriding the default quiescent
// start via VPATH precedence. Measures MLUPS without plotfile overhead.

#include "IO/InputParser.H"
#include "Simulation/ISimulation.H"
#include "Simulation/SimulationFactory.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <cmath>

using namespace muphfasa;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        InputParser params;
        auto sim = makeSimulation(params);
        sim->initialize();

        const auto ncells = params.ncells();
        long long total_cells = 1;
        int Q = 0;
        // Lattice Q depends on DIM and inputs; get from sim's factory via params
        // For MLUPS we need Q: D2Q9=9, D3Q15/19/27. Use lattice_type string.
        std::string lattice = "unknown";
        {
            amrex::ParmParse pp("lbm");
            pp.query("lattice_type", lattice);
            if (lattice == "D2Q9") Q = 9;
            else if (lattice == "D3Q15") Q = 15;
            else if (lattice == "D3Q19") Q = 19;
            else if (lattice == "D3Q27") Q = 27;
            else Q = 9; // fallback
        }
        for (auto n : ncells) total_cells *= n;
        // Account for MPI ranks: total_cells is global (from inputs), not per-rank
        // For MLUPS we report global lattice updates per second.

        const int max_step = params.maxStep();
        const int warmup = std::min(10, max_step/10);
        const int bench_steps = max_step - warmup;

        amrex::Print() << "\n=== Benchmark: " << lattice << " "
                       << "cells=" << total_cells << " Q=" << Q
                       << " max_step=" << max_step
                       << " warmup=" << warmup << " bench_steps=" << bench_steps << " ===\n";

        // Warmup (not timed, also fills caches)
        for (int step = 0; step < warmup; ++step) {
            sim->advance(step);
        }
        amrex::ParallelDescriptor::Barrier();
        double t0 = amrex::second();

        for (int step = warmup; step < max_step; ++step) {
            sim->advance(step);
            // No I/O in bench mode: inputs.bench should have amr.plot_int=-1
        }
        amrex::ParallelDescriptor::Barrier();
        double t1 = amrex::second();
        double elapsed = t1 - t0;
        // Reduce across ranks (max, since all ranks should be similar)
        amrex::ParallelDescriptor::ReduceRealMax(elapsed);

        // MLUPS = Million Lattice (cell) Updates Per Second, the standard LBM
        // community definition (Palabos/waLBerla/HemeLB/OpenLB): cells*steps
        // over wall time, WITHOUT multiplying by the number of discrete
        // velocities Q. Earlier versions of this driver multiplied by Q,
        // which inflated the number by 9x/19x/27x and made it incomparable
        // to MLUPS figures reported elsewhere. MLUDPS (lattice-Direction
        // updates/s) is kept alongside for anyone who wants the Q-scaled
        // (roughly memory-bandwidth-proportional) figure.
        double mlups = (double)total_cells * bench_steps / (elapsed * 1e6);
        double mludps = (double)total_cells * Q * bench_steps / (elapsed * 1e6);
        double ms_per_step = (elapsed / bench_steps) * 1e3;

        amrex::Print() << "Benchmark result: elapsed=" << elapsed << " s for "
                       << bench_steps << " steps\n";
        amrex::Print() << "  MLUPS=" << mlups << "  MLUDPS=" << mludps
                       << "  ms/step=" << ms_per_step << "\n";
        amrex::Print() << "  Config: lattice=" << lattice
                       << " ncells=";
        for (size_t i=0;i<ncells.size();++i) amrex::Print() << ncells[i] << (i+1<ncells.size()?"x":"");
        amrex::Print() << " MPI ranks=" << amrex::ParallelDescriptor::NProcs() << "\n";

        // Also print in parseable form for run_scaling.py / perf.py
        amrex::Print() << "BENCH_RESULT mlups=" << mlups
                       << " ms_per_step=" << ms_per_step
                       << " elapsed=" << elapsed
                       << " mludps=" << mludps << "\n";
    }
    amrex::Finalize();
    return 0;
}
