# TGV 3-D correctness check

Minimal validation counterpart to `../inputs.bench` (perf-only): triply-periodic cubic domain seeded with exact 3-D TGV solution (`../Tgv3dAnalytic.H`). `../InitData.cpp` sets t=0 IC; `test_tgv3d.cpp` compares evolution against same analytic, mirroring `benchmarks/beltrami3d/validate/` and `external/muphfasa/tests/regression/tgv2d/`.

## Build and run
```sh
make                             # builds test_tgv3d3d.gnu.ex (DIM=3)
mkdir -p run && cd run
../test_tgv3d3d.gnu.ex ../inputs
```
Prints max velocity/density error at each `amr.plot_int` step (32³, 200 steps, tau=0.8) and writes `plt`/`analytic_plt` pairs for `fcompare` (cached `fcompare3d` via `helper.py`).

## Expected result
Errors ~1e-4-1e-3 over 200 steps (Ma≈0.005, D3Q19/tau=0.8/N=32). Both plotfiles in physical units (UnitConverter::toPhysical), so don't compare lattice analytic directly.
