// bench.cpp — Generic LBM benchmark driver (no I/O, timed).
//
// Thin clone of external/muphfasa/tests/regression/couette/test_couette.cpp
// and tgv2d/test_tgv2d.cpp: uses InputParser/ISimulation/advance() with the
// benchmark's own InitData.cpp (if any) overriding the default quiescent
// start via VPATH precedence. Measures MLUPS without plotfile overhead.
//
// Step counting (fresh start or restart from a checkpoint via amr.restart):
//   bench.warmup untimed steps, then bench.steps timed steps, both counted from
//   the step the run starts at (0, or just after the checkpoint). If bench.steps
//   is not given, the legacy behaviour applies: sim.max_step is the absolute
//   step count to reach, so bench.steps = sim.max_step - steps already done
//   - bench.warmup.

#include "bench.H"
#include "IO/Checkpoint.H"
#include "IO/InputParser.H"
#include "Simulation/ISimulation.H"
#include "Simulation/SimulationFactory.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>

using namespace muphfasa;

namespace bench {

void do_bench(InputParser& params)
{
    auto sim = makeSimulation(params);
    sim->initialize();

    const auto ncells = params.ncells();
    long long total_cells = 1;

    std::string lattice = params.latticeType();

    for (auto n : ncells) total_cells *= n;
    // For MLUPS we report global lattice updates per second.

    // Lattice steps already advanced: 0 for a fresh start, the checkpoint's
    // time after a restart (the simulation time advances by 1 per step).
    const std::string restart = params.restartFile();
    const long long steps_done = std::llround(sim->time());

    const int max_step = params.maxStep();
    int bench_steps = -1;
    int warmup = -1;
    int save_checkpoint = 0;
    std::string check_file = params.checkpointPrefix();
    {
        amrex::ParmParse ppbench("bench");
        ppbench.query("steps", bench_steps);
        ppbench.query("warmup", warmup);
        ppbench.query("save_checkpoint", save_checkpoint);
        ppbench.query("checkpoint", save_checkpoint);
        ppbench.query("check_file", check_file);
    }
    {
        amrex::ParmParse ppamr("amr");
        ppamr.query("check_file", check_file);
    }

    if (bench_steps <= 0) {
        // legacy: sim.max_step is the absolute step to run up to
        const long long remaining = static_cast<long long>(max_step) - steps_done;
        if (warmup < 0) warmup = static_cast<int>(std::min<long long>(10, remaining / 10));
        bench_steps = static_cast<int>(remaining) - warmup;
    } else if (warmup < 0) {
        warmup = std::min(10, bench_steps / 10);
    }
    if (bench_steps <= 0) {
        amrex::Abort("bench: no timed steps left (set bench.steps, or raise sim.max_step above "
                     "the steps already done + bench.warmup)");
    }
    const long long end_steps = steps_done + warmup + bench_steps;

    amrex::Print() << "\n=== Benchmark: " << lattice << " "
                   << "cells=" << total_cells
                   << " warmup=" << warmup << " bench_steps=" << bench_steps << " ===\n";
    if (!restart.empty()) {
        amrex::Print() << "Restarted from " << restart << ": " << steps_done
                       << " steps already done, run continues to step " << end_steps << "\n";
    }

    // Warmup (not timed, also fills caches)
    for (int step = 0; step < warmup; ++step) {
        sim->advance();
    }
    amrex::ParallelDescriptor::Barrier();
    double t0 = amrex::second();

    for (int step = 0; step < bench_steps; ++step) {
        sim->advance();
        // No I/O in bench mode: inputs.bench should have amr.plot_int=-1
    }
    amrex::ParallelDescriptor::Barrier();
    double t1 = amrex::second();
    double elapsed = t1 - t0;
    // Reduce across ranks (max, since all ranks should be similar)
    amrex::ParallelDescriptor::ReduceRealMax(elapsed);

    // MLUPS = Million Lattice (cell) Updates Per Second, the standard LBM
    // community definition (Palabos/waLBerla/HemeLB/OpenLB): cells*steps
    // over wall time.
    double mlups = (double)total_cells * bench_steps / (elapsed * 1e6);
    double ms_per_step = (elapsed / bench_steps) * 1e3;

    amrex::Print() << "Benchmark result: elapsed=" << elapsed << " s for "
                   << bench_steps << " steps\n";
    amrex::Print() << "  MLUPS=" << mlups
                   << "  ms/step=" << ms_per_step << "\n";
    amrex::Print() << "  Config: lattice=" << lattice
                   << " ncells=";
    for (size_t i=0;i<ncells.size();++i) amrex::Print() << ncells[i] << (i+1<ncells.size()?"x":"");
    amrex::Print() << " MPI ranks=" << amrex::ParallelDescriptor::NProcs() << "\n";

    // Also print in parseable form for run_scaling.py / perf.py
    amrex::Print() << "BENCH_RESULT mlups=" << mlups
                   << " ms_per_step=" << ms_per_step
                   << " elapsed=" << elapsed
                   << "\n";

    if (save_checkpoint) {
        // Label by steps advanced so chained restarts keep counting correctly
        // (equals sim.max_step for a fresh run with the legacy step counting).
        const int chk_step = static_cast<int>(std::llround(sim->time()));
        std::string cname = checkpointName(check_file, chk_step);
        amrex::Print() << "Saving final checkpoint: " << cname
                       << " at step " << chk_step << "\n";
        writeCheckpoint(cname, sim->f(), sim->geom(), chk_step, sim->time());
        amrex::ParallelDescriptor::Barrier();
        amrex::Print() << "Checkpoint done: " << cname << "\n";
    }
}

} // namespace bench
