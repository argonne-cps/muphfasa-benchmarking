# TGV 3-D validation — Re=1600 vs. Dairay et al. 2017 DNS

Triply-periodic 3-D Taylor-Green vortex (Taylor & Green 1937; Brachet et al. 1983)
at **Re = 1600**, matching the DNS reference in `TGV_Re1600.dat` (data from
Incompact3d 512³, domain (2π)³, doi:10.1016/j.jcp.2017.02.035 — see file header and the
citation notice there). `validate` is the mirror of `../inputs.bench` (perf-only):
same IC (`../Tgv3dAnalytic.H` via `../InitData.cpp`, one wavelength per domain) but
run as a correctness/validation case rather than a timed benchmark.

## What is checked

Two complementary layers:

1. **In-situ energy balance** (`test_tgv3d.cpp`): the identity
   `KE(0)-KE(T) = ∫ ε dt` (`ε = ν⟨|∇u|²⟩ = 2νℰ`) is exact for periodic
   incompressible NS at any Re. The driver samples `E_k = ½⟨|u|²⟩`, `ε`, and
   enstrophy `ℰ = ½⟨|ω|²⟩` at every `amr.plot_int` step (central differences,
   periodic ghost fill, converted to physical units) and reports the balance
   error plus a CSV. This holds regardless of vortex stretching.

2. **DNS benchmark** (`analyze_tgv.py`): reads the AMReX plotfiles via `yt`
   (or the in-situ `tgv_diagnostics.csv`) and compares `E_k(t)`, `ε(t)`,
   `ℰ(t)` against `TGV_Re1600.dat` — L2 errors, peak dissipation magnitude
   and timing, and overlay plots. This is the Re=1600 validation proper.

## Inputs

- `inputs` — default **N=512, Re=1600, t=20, u_lat=0.1** (matches the
  512³ DNS). Runs in the benchmark's own units: domain `[0,2π]³`,
  `units.ref_length=2π`, `U₀=1`, `ν=1/1600`, so results compare 1:1 with
  `TGV_Re1600.dat` (t=0: `E_k=0.125`, `⟨(∂u/∂x)²⟩=0.125`). Note
  `units.Re = U₀·ref_length/ν = 2π·1600` (muphfasa's Re is built on
  `ref_length`). At N=512 `dx=2π/512`, `dt=u_lat·dx/U₀≈1.23e-3`, `τ≈0.515`,
  so `t=20` takes ~16.3k steps; `amr.plot_int=163` gives ~100 outputs
  (Δt≈0.2, ~300 GB — increase if disk-limited). `bench.sample_int` sets the
  in-situ CSV cadence (defaults to `plot_int`).
- `inputs.small` — N=64, Re_bench=100 smoke test (t=2, ~200 steps) for local
  sanity without HPC (`units.Re = 2π·100`). Energy-balance error is numerical
  dissipation of the coarse LBM grid: ~19% (N=32), ~7% (N=64), ~3% (N=128).
  Re=1600 at small N is unstable (`τ→0.5`); lower Re for small N.

`sim.stop_time` (physical, `L/U` units) is honored when >0: the driver
converts it to lattice steps via `UnitConverter::toLatticeTime`, so
`sim.max_step` can be left as a cap. `lbm.u_lat` pins the lattice velocity
(`Ma = u_lat/c_s`); `τ` is derived so `Re`/`N`/`u_lat` stay consistent.

## Build and run

```sh
make                               # builds test_tgv3d3d.gnu.ex (DIM=3)
mkdir -p run && cd run
../test_tgv3d3d.gnu.ex ../inputs
# or smoke test:
../test_tgv3d3d.gnu.ex ../inputs.small
```

Each `amr.plot_int` step prints `t  E_k  ε  ℰ`; at the end the energy-balance
`rel_err` and the hint to run `analyze_tgv.py`. Plotfiles `plt*/` are in
physical units; `tgv_diagnostics.csv` is written in `run/`.

## Post-processing

```sh
# Activate the venv that has yt (per project instructions)
source /home/biboyd/.venv/bin/activate

# From plotfiles (yt path — reads plt*/ in run/):
python3 ../analyze_tgv.py --run-dir run --benchmark ../TGV_Re1600.dat

# From the in-situ CSV (no yt needed, faster):
python3 ../analyze_tgv.py --csv run/tgv_diagnostics.csv --benchmark ../TGV_Re1600.dat

# Save plots elsewhere / tune tolerances:
python3 ../analyze_tgv.py --run-dir run --out-prefix ./tgv_compare --peak-eps-tol 0.20
```

Outputs `tgv_compare_Ek.png`, `*_dissipation.png`, `*_enstrophy.png`,
`*_combined.png`, plus a console PASS/FAIL on peak-ε and KE-L2. Full
comparison needs the run to reach t≈9 where dissipation peaks; shorter
smoke runs will show truncated curves.

## Expected result

- **Energy balance:** <10% `rel_err` (trapezoidal time integration; shrinks
  with smaller `plot_int`).
- **DNS comparison (N=512):** peak `ε≈0.0128` at `t≈9`; KE decay and
  enstrophy within ~15–20% depending on resolution/compressibility. Coarser
  N will deviate more — the diagnostic is for convergence, not a hard gate.
