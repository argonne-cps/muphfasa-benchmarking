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
		amrex::Print() << "Benchmarking Lid-Driven Cavity 2D" << std::endl;
        InputParser params;
		bench::do_bench(params);
    }
    amrex::Finalize();
    return 0;
}
