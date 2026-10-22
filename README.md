# muphfasa-benchmarking

Scale benchmarks for [muphfasa](https://github.com/biboyd/muphfasa) on Aurora (Intel Data Center GPU Max, SYCL + MPI). This repo pulls `muphfasa` as a submodule (`external/muphfasa`) and adds large-scale benchmarks that don't belong in the core repo's regression tests.

## Benchmarks

| Case | Dim | Lattice | Sizes (weak) | Physics |
|------|-----|---------|--------------|---------|
| `tgv2d` | 2D | D2Q9 | 1024² → 8192² (1024²/GPU) | Taylor-Green vortex, periodic, analytical decay `exp(-2νk²t)` |
| `lid3d` | 3D | D3Q19/27 | 256³ → 1024³ (256³/GPU) | Lid-driven cavity, Re=100 (`inputs.bench`; override `units.Re`/`lbm.tau` for other Re), `MovingWall` |
| `beltrami3d` | 3D | D3Q19/27 | 128³ → 512³ (128³/GPU) | 3D Beltrami flow, fully periodic, exact N-S |

All share a single `bench.cpp` driver (thin clone of `external/muphfasa/tests/regression/couette/test_couette.cpp`) that times `advance()` without I/O and reports **MLUPS** (Million Lattice/cell Updates Per Second — the standard LBM community definition, `cells*steps/elapsed`, with no multiplication by the lattice's number of discrete velocities Q) alongside **MLUDPS** (the Q-scaled, roughly memory-bandwidth-proportional variant) for anyone who wants that figure instead.

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

# Weak scaling sweep (uses helper.py like couette). Each sweep point actually
# launches nprocs=nodes*ppn ranks via mpiexec; --max-grid-size is a fixed AMReX
# box tile size kept independent of --sizes so the domain is decomposed into
# enough boxes for those ranks to do real work (don't set it equal to the
# domain size, or every rank past the first sits idle on an empty box).
python3 run_scaling.py --case tgv2d --sizes 1024 2048 4096 --nodes 1 2 4 --ppn 1 --max-grid-size 64
```

## Aurora

```sh
# Build for Aurora (SYCL, MPI, Intel)
make -C benchmarks/tgv2d DIM=2 USE_MPI=TRUE USE_SYCL=TRUE COMP=intel
make -C benchmarks/lid3d DIM=3 USE_MPI=TRUE USE_SYCL=TRUE COMP=intel

# Batch (PBS, 12 ranks/node, 6 GPUs/node) — requires -l filesystems=flare, see skills/aurora-overview.md
qsub -A <ALCF_PROJECT> -q debug -l select=2:system=aurora -l walltime=00:30:00 -l filesystems=flare batch/aurora_bench.sh
# Check queue first: qstat -u $USER; debug max 1 job/user, prod routes small/medium/large
# aurora_bench.sh: mpiexec -n $NNODES --ppn 12 --gpus-per-node 6 ./bench*.ex inputs.bench
# Do not submit from /soft/modulefiles; prefer /opt/aurora APE; proxy http://proxy.alcf.anl.gov:3128 set in job
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

`Make.muphfasa` (`external/muphfasa/Make.muphfasa`) holds all solver sources and AMReX modules; each benchmark's `GNUmakefile` only sets `TOP`, `DIM`, `EBASE`, `CEXE_sources` and `include $(TOP)/Make.muphfasa`. Local `InitData.cpp` wins via `VPATH` precedence (`VPATH_LOCATIONS+=..` before include).

## Aurora skills

This repo vendors the relevant `alcf-agent-skillset` docs in `skills/`:
- `aurora-overview.md` — 10,624 nodes, 104 cores + 6 PVC (12 tiles)/node, flare/DAOS mounts, queues (`debug` 1-2 nodes, `debug-scaling` 2-256, `prod`/`capacity`), `-l filesystems=flare`, proxy `proxy.alcf.anl.gov:3128`, cores 0/52 reserved.
- `etiquette.md` / `safety.md` — do not build/run heavy work on login nodes, `qstat -u $USER` before submit, `BENCH_J=4` for local makes, no credentials in job scripts, never `rm -rf` shared dirs.

`batch/aurora_bench.sh` follows these: `#PBS -l filesystems=flare`, `module use /soft/modulefiles` with note to prefer `/opt/aurora`, `MPICH_GPU_SUPPORT_ENABLED=1`, `FI_CXI_ATS=0`, lazy build on compute, and checkpoint note for node failures.

## Validate vs Perf

Plotfiles (`plt<step>`) are written in **physical (SI) units**, not lattice units —
`SimulationSetup::writeOutput` converts density/velocity via `UnitConverter::toPhysical`
before writing (upstream `muphfasa` commit `c736661`, "output in physical units"). The
solver's internal state stays in lattice units; only I/O is converted. Any `InitData.cpp`
override (`tgv2d`, `beltrami3d`) takes `const UnitConverter&` and must convert its IC
amplitude from the physical `units.ref_velocity` via `units.toLatticeVelocity()` — see
either file for the pattern. A correctness-check driver that writes its own analytic
plotfile (like `beltrami3d/validate/test_beltrami3d.cpp`) must likewise convert its
analytic solution to physical units before `fcompare`, or the comparison will be
comparing physical sim output against a lattice-unit analytic field and report bogus
errors.

- **Validate** (small, correctness): `inputs` with `amr.plot_int=200`, fcompare vs analytic – same as regression.
  - `tgv2d`: validated upstream, `external/muphfasa/tests/regression/tgv2d` (exact 2D TGV decay).
  - `beltrami3d`: validated in this repo, `benchmarks/beltrami3d/validate/` (exact 3D Beltrami decay,
    mirrors the tgv2d test — see its README for build/run/fcompare steps).
  - `lid3d`: no closed-form solution exists for lid-driven cavity (Ghia et al. is a numerical
    reference, not analytic), so there's no `fcompare`-vs-analytic path here; at minimum, compare
    against the upstream 2D `tests/regression/lid_driven_cavity` steady-state profile or a
    literature Re=100 reference before trusting `lid3d`'s perf numbers.
- **Perf** (large, no I/O): `inputs.bench` with `amr.plot_int=-1`, `sim.max_step=100`, `lbm.tau=0.8`, measure `MLUPS`.

## Performance History (`results/history.jsonl`)

`perf.py log` (and `run_scaling.py`, which calls it) appends a row to the git-tracked
`results/history.jsonl` + `.cache/perf.db` whenever it runs on the `main` branch — there is
**no CI gate yet**. In practice this means `results/history.jsonl` today is "whatever anyone
happened to run locally on `main`," not a curated/vetted performance record: a noisy laptop run
or a one-off wrong build flag becomes part of the permanent history just like a clean Aurora
run. Each row's `ci` field (0 for all current rows) and `hostname` (part of the dedup key, so
different machines don't collide) let you filter by provenance when reading the history. Wiring
up CI-only ingestion so `main` pushes are the only trusted source is tracked as future work in
`/home/biboyd/.opencode/plan/benchmark-perf-db.md` (§6) — until then, treat `history.jsonl`
rows as indicative, not authoritative, unless you know who ran them and on what.

See `benchmarks/README.md` for problem sizes, Re/tau, and expected MLUPS on Aurora.
