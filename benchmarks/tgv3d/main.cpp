// bench.cpp — Generic LBM benchmark driver (no I/O, timed).
//
// Thin clone of external/muphfasa/tests/regression/couette/test_couette.cpp
// and tgv2d/test_tgv2d.cpp: uses InputParser/ISimulation/advanceFused() with the
// benchmark's own InitData.cpp (if any) overriding the default quiescent
// start via VPATH precedence. Measures MLUPS without plotfile overhead.

#include "IO/Checkpoint.H"
#include "IO/InputParser.H"
#include "Simulation/ISimulation.H"
#include "Simulation/SimulationFactory.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include "bench.H"

#include <cmath>

using namespace muphfasa;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
		amrex::Print() << "Benchmarking TGV 3D" << std::endl;
        InputParser params;
		bench::do_bench(params);
    }
    amrex::Finalize();
    return 0;
}
