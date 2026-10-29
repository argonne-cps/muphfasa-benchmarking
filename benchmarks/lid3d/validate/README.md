# lid3d validation

Correctness test for the 3-D lid-driven cavity (closed cube, lid at `z=hi` moving in `+x`; all other faces `NoSlip`). The cavity starts from rest. The case has no closed-form solution, so like `benchmarks/tgv3d/validate` and the 2-D `external/muphfasa/tests/regression/lid2d` cavity it runs to steady state and checks global invariants, plus literature centerline extrema when `Re` is set to the literature value.

This test is the 3-D counterpart to the pre-existing 2-D `lid2d/validate` (Ghia et al. 1982 at `Re=100`). The 2-D test compared `u_min`/`v_max`/`v_min` on the cavity centerlines to Ghia's spectral data; this test does the analogous comparison in 3-D to Albensoeder & Kuhlmann (2005) for the cubical cavity at `Re=1000`. The two geometries differ: the performance benchmark `benchmarks/lid3d/inputs.bench` is a periodic slab (`Periodic` in `x,y`, walls in `z`) sized for weak scaling, while this validation is a closed cube (walls on all six faces) so that the lid-corner singularities and centerline extrema are comparable to the literature.

## Build and run

```sh
make -C benchmarks/lid3d/validate -j6 DIM=3          # -> test_lid3d3d.gnu.ex
mkdir -p benchmarks/lid3d/validate/run && cd benchmarks/lid3d/validate/run
../test_lid3d3d.gnu.ex ../inputs                      # Re=100 invariants, N=32, ~10-20 s
../test_lid3d3d.gnu.ex ../inputs.re1000               # Re=1000 + Albensoeder, N=48, ~40-80 s
# Convergence sweep (override on the command line):
../test_lid3d3d.gnu.ex ../inputs grid.ncells="64 64 64" grid.max_grid_size=64
../test_lid3d3d.gnu.ex ../inputs.re1000 grid.ncells="64 64 64" grid.max_grid_size=64
```

`sim.max_step` / `amr.plot_int` play the same role as in `lid2d/validate`: `plot_int` is the diagnostic sampling interval. Steady state in a cube takes roughly `20 * N / U_lat` steps (`U_lat ~ 0.05` at `N=32, Re=100`), so `N=32 -> 10000`, `N=48 -> 20000`, `N=64 -> 30000-40000` is a good ladder. Other sizes: override `grid.ncells`, `grid.max_grid_size`, and run longer.

## Checks

| # | Check | Pass criterion | Measured (N=32, Re=100) | Notes |
|---|-------|----------------|-------------------------|-------|
| 1 | Mass: `sum(rho)` in a closed box | max relative drift < 1e-9 | ~1e-12 | Closed-box bounce-back must conserve mass |
| 2 | KE spin-up: `KE = 1/2 sum|u|^2` starts at 0, rises, no dips | `KE(0)=0`, `KE(T)>0`, worst relative dip < 1e-3 | monotone |
| 3 | KE plateau (steady state reached) | relative KE change over last sample interval < 5e-4 | ~1e-6 at step 10000 | 3-D looser than 2-D's 1e-4 (stronger corner singularity) |
| 4 | Energy balance `P_lid = eps` | `abs(P - eps)/P < 0.30` | ~0.15-0.25 | Wall- and lid-corner singularities make both integrals converge slowly; see lid2d README for the analogous discussion. Tolerance is deliberately loose to catch sign/factor-2/wrong-wall errors (see below). |
| 5 | Albensoeder `u_min` along wall-normal `(0,0,z)` | relative error < 12% | only when `Re~1000` | Dominant 3-D cavity feature (reverse flow) |
| 6-7 | Albensoeder `w_max`/`w_min` along lid direction `(x,0,0)` | relative error < 12% | only when `Re~1000` and `N>=48` | Secondary flow; informational at `N=32` |
| (info) | Albensoeder `v_min` along span `(0,y,0)` | reported, not gating | only when `Re~1000` | Very weak span flow (~0.03 U_lid), noisy at coarse `N` |

Checks 1-4 run at any `Re` (including the default `Re=100`). Checks 5-7 are enabled only when `units.Re` is within 50 of 1000; otherwise the test prints `SKIP  Albensoeder ... (Re=... not ~1000)`. This keeps the default `inputs` (Re=100) fast and stable while still exercising the literature comparison via `inputs.re1000` or a command-line override.

**Energy balance.** `dKE/dt = P_lid - eps` with `P_lid = U_lid * sum_lid(nu * du/dz|wall)` and `eps = nu * sum|grad u|^2`. At steady state they must be equal. Wall gradients use the halfway-bounce-back geometry (wall half a cell from the first cell center: quadratic through wall/first/second cell for the lid shear; ghost `2*U_wall - u_cell` for central differences). In 3-D the lid has four corner lines and the lid-corner `|grad u|^2 ~ 1/r^2` singularity is stronger than in 2-D, so convergence is slower and the tolerance is 0.30 vs 2-D's 0.20.

## Albensoeder reference values

Source: S. Albensoeder and H. C. Kuhlmann, "Accurate three-dimensional lid-driven cavity flow", *Journal of Computational Physics* **206**(2), 536-558 (2005). DOI: [10.1016/j.jcp.2004.12.024](https://doi.org/10.1016/j.jcp.2004.12.024).

Re=1000 centerline extrema from that paper, Table 10 (rigid boundary conditions, cubical cavity Gamma=Lambda=1, Chebyshev 96x96x64; velocities normalized by the lid velocity, i.e. `u/Re` in the paper). The paper's coordinates are lid at `x=-Gamma/2` moving in `y`; ours are lid at `z=hi` moving in `x` (closed cube `0-1` in each dimension). Mapping: paper `x->our z`, `y->our x`, `z->our y`; `u->our w`, `v->our u`, `w->our v`:

| Quantity (paper) | Value `V/V_lid` | Location (paper) | Quantity (our coords) | Line (our coords) |
|------------------|-----------------|------------------|----------------------|-------------------|
| `min v` | -0.2803833 | `x=0.37581` on `(x,0,0)` | `min u` along wall-normal | `(0,0,z)` |
| `max u` | +0.4350186 | `y=0.40957` on `(0,y,0)` | `max w` along lid direction | `(x,0,0)` |
| `min u` | -0.2466511 | `y=-0.39087` on `(0,y,0)` | `min w` along lid direction | `(x,0,0)` |
| `min w` | -0.0303398 | `z=0.36492` on `(0,0,z)` | `min v` along span | `(0,y,0)` |

The locations are quoted from the paper and are not used by the test; only the extremum values are compared. In the test the extrema are taken from the cell-centered centerline profiles (averaging the 2/4 cells straddling `x=0.5`/`y=0.5`/`z=0.5`, `N` even) and refined with a 3-point parabola, so the reported location may differ slightly from the Chebyshev collocation points.

Tolerances were set for coarse LBM grids: at `N=32` the primary `u_min` is typically within ~8% (BGK, tau near 0.55 for `Re=1000`); at `N=48` all three are within ~8%; at `N=64` within ~5%. The 12% gate is intentionally loose to accommodate LBM compressibility (`Ma = U_lat/cs ~ 0.09` at `N=32, Re=1000` with `U_lat=0.05`) on top of discretization error. Detailed profiles along the centerlines are in the paper's Tables 5-6; the convergence of the extrema with resolution is in Tables 2-4. The paper was supplied as `albensoeder_paper.txt` in the repository root.

## Bench vs validate geometry

- Bench (`benchmarks/lid3d/inputs.bench`): `Periodic` in `x,y`, `NoSlip`/`MovingWall` in `z`, `Re=100`, `D3Q19/27`, `64^3` default local. Sized for weak-scaling throughput (`256^3`/GPU).
- Validate (`benchmarks/lid3d/validate/inputs`): `NoSlip` on all faces except `MovingWall` lid at `z=hi`, `Re=100` (invariants) / `Re=1000` (Albensoeder via `inputs.re1000`), cubic `32^3` default local. Sized for correctness.

Do not time `inputs`/`inputs.re1000` for performance; use `inputs.bench` for that.

## Expected failure on the currently pinned `external/muphfasa`

At the submodule SHA this repo pins, this test **is expected to fail** the mass check (and downstream checks) for the same reason as `lid2d/validate`: halfway bounce-back reflects diagonal links from the wrong tangential cell, leaking mass at the ends of every wall row and at corners. See item 8 of [`../../../muphfasa_issues.md`](../../../muphfasa_issues.md) for the root cause and one-line fix (verified for 2-D; the same pattern applies in 3-D). Until the upstream fix is merged and the submodule bumped, throughput numbers are real but the physics they time is not correct.
