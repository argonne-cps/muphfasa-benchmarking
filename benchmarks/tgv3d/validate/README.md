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

- `inputs` — default **N=512, Re=1600, t=20, u_lat=0.04** (matches the
  512³ DNS). At N=512, `dt = u_lat·dx/U₀ = 0.04/512 ≈ 7.8e-5`, so `t=20`
  needs ~256k steps; `amr.plot_int=2560` gives 100 outputs (Δt≈0.2). Full
  run stores ~300 GB of plotfiles — increase `plot_int` if disk-limited;
  the CSV is always written.
- `inputs.small` — N=32, Re=100 smoke test (t=2, 800 steps) for local
  sanity without HPC. Re=1600 at N=32 is unstable (`τ→0.5`) and
  under-resolved; lower Re for small N or override on the CLI:
  ```sh
  ./test_tgv3d3d.gnu.ex inputs grid.ncells="64 64 64" grid.max_grid_size=64 \
      units.Re=400 sim.max_step=5000 amr.plot_int=500
  ```

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
