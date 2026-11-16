// main.cpp — Entry point for the Beltrami 3D benchmark.
// Overrides the solver's default src/main.cpp (a same-named file in the
// build directory wins over VPATH_LOCATIONS) so the timed driver in
// ../bench.cpp runs instead of the plain time-advance loop.

#include "bench.H"
#include "IO/InputParser.H"

#include <AMReX.H>
#include <AMReX_Print.H>

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        amrex::Print() << "Benchmarking Beltrami 3D" << std::endl;
        muphfasa::InputParser params;
        bench::do_bench(params);
    }
    amrex::Finalize();
    return 0;
}
