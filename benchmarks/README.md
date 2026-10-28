# Benchmarks

All benchmarks use `bench.cpp` (no I/O) and report `MLUPS`.

## Cases

| Case | Grid (bench default) | Aurora weak (per GPU) | Lattice | BC |
|------|---------------------|----------------------|---------|----|
| `tgv2d` | 64² (local) → 1024² | 1024²/GPU → 8192² | D2Q9 | Periodic |
| `tgv3d` | 64³ → 128³ | 128³/GPU → 512³ | D3Q19/27 | Periodic xyz (TGV) |
| `lid2d` | 64² (local) → 1024² | 1024²/GPU → 4096² | D2Q9 | Closed cavity: NoSlip x-lo/x-hi/y-lo, MovingWall lid at y-hi (Re=100) |
| `lid3d` | 64³ → 256³ | 256³/GPU → 1024³ | D3Q19/27 | Periodic x,y, z walls (NoSlip/MovingWall) |
| `beltrami3d` | 64³ → 128³ | 128³/GPU → 512³ | D3Q19/27 | Periodic xyz (Beltrami) |

## Build (uses `external/muphfasa/Make.muphfasa`)

```sh
make -C benchmarks/tgv2d -j6 DIM=2
make -C benchmarks/tgv3d -j6 DIM=3
make -C benchmarks/lid2d -j6 DIM=2
make -C benchmarks/lid3d -j6 DIM=3
make -C benchmarks/beltrami3d -j6 DIM=3
# Aurora SYCL
make -C benchmarks/lid2d DIM=2 USE_MPI=TRUE USE_SYCL=TRUE COMP=intel -j6
make -C benchmarks/lid3d DIM=3 USE_MPI=TRUE USE_SYCL=TRUE COMP=intel -j6
make -C benchmarks/tgv3d DIM=3 USE_MPI=TRUE USE_SYCL=TRUE COMP=intel -j6
```

`VPATH` ensures `benchmarks/tgv2d/InitData.cpp` (TGV) and `benchmarks/beltrami3d/InitData.cpp` override `external/muphfasa/src/Simulation/InitData.cpp`.

## Run

```sh
# Local, no plotfiles: 10 warmup + 100 timed steps
./benchmarks/tgv2d/bench2d.gnu.ex benchmarks/tgv2d/inputs.bench
# Override size/lattice inline
./benchmarks/lid3d/bench3d.gnu.ex benchmarks/lid3d/inputs.bench grid.ncells="128 128 128" lbm.lattice_type=D3Q27 amr.plot_int=-1

# Scaling sweep (helper.py)
python3 run_scaling.py --case tgv2d --sizes 1024 2048 4096 --nodes 1 2 4 --plot tgv2d_weak.png
python3 run_scaling.py --case lid2d --sizes 512 1024 2048 --nodes 1 2 4 --mode weak
python3 run_scaling.py --case lid3d --sizes 128 256 512 --lattice D3Q19 --mode weak
```

`lid2d`/`lid3d` pin `lbm.tau` with `units.Re=100`, so the lattice lid velocity is derived as `u_lat = nu_lat*Re/N` and shrinks as the grid grows (0.1 at N=64 for `lid2d`). Sweeping `lid2d` *below* N≈43 trips the `u_lat>0.15` warning, and below N≈22 the `u_lat<0.3` assert aborts the run.

Correctness check for `lid2d` (mass, KE, power balance, Ghia Re=100 extrema): `benchmarks/lid2d/validate/` (see its README; currently fails on the pinned `muphfasa` due to the bounce-back diagonal-link bug in `muphfasa_issues.md` item 8).

## Aurora

See `batch/aurora_bench.sh` for PBS + `mpiexec --gpus-per-node 6`. Each rank = 1 GPU tile (PVC has 2 tiles/GPU, 6 GPUs/node = 12 tiles). Use `MPICH_GPU_SUPPORT_ENABLED=1`.

## Metrics

`bench.cpp` prints `BENCH_RESULT mlups=... ms_per_step=... elapsed=...` (global `Nx*Ny*Nz*Q*steps / 1e6 / sec`). `run_scaling.py` plots weak (fixed per-GPU) vs strong (fixed global) efficiency.

## Validate vs Perf

- Validate: `inputs` (small, `amr.plot_int=200`, fcompare vs analytic) – same as `external/muphfasa/tests/regression/tgv2d`.
- Perf: `inputs.bench` (`amr.plot_int=-1`, `sim.max_step=110`, `lbm.tau=0.8`).
