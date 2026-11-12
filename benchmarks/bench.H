#ifndef MUPHFASA_BENCH_H
#define MUPHFASA_BENCH_H

#include <AMReX.H>
#include <AMReX_ParmParse.H>

#include "IO/InputParser.H"

namespace bench {

void do_bench(muphfasa::InputParser& params);

} // namespace bench

#endif
