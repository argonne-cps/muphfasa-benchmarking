# muphfasa-benchmarking

Scale benchmarks for [muphfasa](https://github.com/biboyd/muphfasa) on Aurora (Intel Data Center GPU Max, SYCL + MPI). This repo pulls `muphfasa` as a submodule (`external/muphfasa`) and adds large-scale benchmarks that don't belong in the core repo's regression tests.

## Benchmarks

| Case | Dim | Lattice | Sizes (weak) | Physics |
|------|-----|---------|--------------|---------|
| `tgv2d` | 2D | D2Q9 | 1024² → 8192² (1024²/GPU) | Taylor-Green vortex, periodic, analytical decay `exp(-2νk²t)` |
| `lid3d` | 3D | D3Q19/27 | 256³ → 1024³ (256³/GPU) | Lid-driven cavity, Re=1000/3200, `MovingWall` |
| `beltrami3d` | 3D | D3Q19/27 | 128³ → 512³ (128³/GPU) | 3D Beltrami flow, fully periodic, exact N-S |

All share a single `bench.cpp` driver (thin clone of `external/muphfasa/tests/regression/couette/test_couette.cpp`) that times `advance()` without I/O and reports **MLUPS**.

## Quick Start (local, no Aurora)

```sh
git clone --recurse-submodules git@github.com:biboyd/muphfasa-benchmarking.git
cd muphfasa-benchmarking

# 2D TGV, 1024², 100 steps, no plotfiles
make -C benchmarks/tgv2d -j6 DIM=2
./benchmarks/tgv2d/bench2d.gnu.ex benchmarks/tgv2d/inputs.bench

# 3D Lid, 256³
make -C benchmarks/lid3d -j6 DIM=3
./benchmarks/lid3d/bench3d.gnu.ex benchmarks/lid3d/inputs.bench

# Weak scaling sweep (uses helper.py like couette)
python3 run_scaling.py --case tgv2d --sizes 1024 2048 4096 --nodes 1 2 4
```

## Aurora

```sh
# Build for Aurora (SYCL, MPI, Intel)
make -C benchmarks/tgv2d DIM=2 USE_MPI=TRUE USE_SYCL=TRUE COMP=intel
make -C benchmarks/lid3d DIM=3 USE_MPI=TRUE USE_SYCL=TRUE COMP=intel

# Batch (PBS, 12 ranks/node, 6 GPUs/node)
qsub -l select=4:system=aurora -l walltime=00:30:00 batch/aurora_bench.sh
# aurora_bench.sh: mpiexec -n $NNODES --ppn 12 --gpus-per-node 6 ./bench*.ex inputs.bench
```

`run_scaling.py` parses `run.log` for `MLUPS` and `ms/step` and plots weak/strong efficiency vs. 1-node baseline (see `benchmarks/README.md`).

## Layout

```
muphfasa-benchmarking/
  external/muphfasa/          # submodule, pinned SHA (update: git submodule update --remote)
  benchmarks/
    bench.cpp                 # generic timed driver (no analytic plotfiles)
    tgv2d/   GNUmakefile, inputs.bench, InitData.cpp (TGV, overrides src/Simulation/InitData.cpp)
    lid3d/   GNUmakefile, inputs.bench (lid-driven, z-wall like couette)
    beltrami3d/ GNUmakefile, inputs.bench, InitData.cpp, BeltramiAnalytic.H
  batch/
    aurora_bench.sh
  helper.py                   # shared build/run/parse helpers (from couette)
  run_scaling.py
```

`Make.muphfasa` (`external/muphfasa/Make.muphfasa`) holds all solver sources and AMReX modules; each benchmark's `GNUmakefile` only sets `TOP`, `DIM`, `EBASE`, `CEXE_sources` and `include $(TOP)/Make.muphfasa`. Local `InitData.cpp` wins via `VPATH` precedence.

## Validate vs Perf

- **Validate** (small, correctness): `inputs` with `amr.plot_int=200`, fcompare vs analytic (TGV/Beltrami) – same as regression.
- **Perf** (large, no I/O): `inputs.bench` with `amr.plot_int=-1`, `sim.max_step=100`, `lbm.tau=0.8`, measure `MLUPS`.

See `benchmarks/README.md` for problem sizes, Re/tau, and expected MLUPS on Aurora.
