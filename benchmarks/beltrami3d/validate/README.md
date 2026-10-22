# Beltrami 3-D correctness check

Minimal validation counterpart to `../inputs.bench` (perf-only, no
correctness check): a triply-periodic cubic domain seeded with the exact
Beltrami solution (`../BeltramiAnalytic.H`). `../InitData.cpp` sets this as
the t=0 initial condition; `test_beltrami3d.cpp` then compares the solver's
evolution against the same analytic solution at later times, mirroring
`external/muphfasa/tests/regression/tgv2d/test_tgv2d.cpp`.

## Build and run

```sh
make                             # builds test_beltrami3d3d.gnu.ex (DIM=3, no MPI)
mkdir -p run && cd run
../test_beltrami3d3d.gnu.ex ../inputs
```

Prints the max velocity/density error against the analytic solution at each
`amr.plot_int` step (32³, 200 steps, tau=0.8 by default), and writes
`plt<step>`/`analytic_plt<step>` pairs that can be diffed with AMReX's
`fcompare`:

```sh
python3 -c "
from helper import ensure_fcompare
print(ensure_fcompare(dim=3))" # from repo root, builds+caches fcompare3d
<fcompare_path> -n 0 run/plt0000200 run/analytic_plt0000200
```

## Expected result

Absolute velocity/density errors stay ~1e-4-1e-3 over the 200-step run
(Ma≈0.005, well into the low-Mach-error regime for `D3Q19`/`tau=0.8`/N=32).
Both `plt<step>` and `analytic_plt<step>` are written in **physical (SI)
units** (muphfasa converts lattice → physical in `SimulationSetup::writeOutput`
and this test mirrors that conversion for its own analytic/diff fields before
writing and diffing, matching upstream `tests/regression/tgv2d/test_tgv2d.cpp`
since the "output in physical units" change) — with `units.ref_length=32`,
`dx=1` here, so the physical and lattice magnitudes happen to coincide, but
don't assume that holds for other `inputs` configs.
This is a sanity check, not a resolution/decay-rate convergence study (see
`tgv2d`'s `convergence.py` for that pattern if one is wanted for Beltrami
later) -- its purpose is only to catch a grossly wrong IC/BC in
`../InitData.cpp`/`../BeltramiAnalytic.H` before trusting `../inputs.bench`
performance numbers.
