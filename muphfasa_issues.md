# muphfasa issues found during benchmarking validation (do not fix here)

These require changes in `muphfasa` (not `muphfasa-benchmarking`). Listed for later.

1. **InputParser ignores `units.ref_viscosity` but benchmark templated it**
   - `benchmarks/lid3d/inputs.bench` originally used `units.ref_viscosity = 0.01` and `units.Ma = 0.05`. `InputParser` has no `ref_viscosity()` accessor; `Units` derives `nu_ref = U_ref * L_ref / Re`. Unused ParmParse keys emit `ParmParse: unused key` warnings and confuse users. Suggest either add `ref_viscosity` accessor or document that viscosity is Re-derived and warn on unknown keys.

2. **Duplicate `sim.max_step` in TGV bench inputs**
   - Imported `tests/regression/tgv2d/inputs` had `sim.max_step=2000` then later overwritten to `110` for bench. AMReX overwrites silently; the duplicate is confusing. Clean source comment already fixed in benchmarking copy. Upstream could provide `inputs.bench` variant.

3. **Make.muphfasa hardcodes EB module**
   - `USE_EB=TRUE` always pulls `EB` (requires `EB2` etc.) even for periodic benchmarks that don't need EB. Adds compile time. Suggest make EB opt-in per-executable or provide fast-path for bench.

4. **UnitConverter requires Re even when tau pinned**
   - When `lbm.tau` is pinned (TGV, bench), `UnitConverter` still computes `nu_ref = U*L/Re` and warns about `u_lat` derived from tau. For bench we set dummy `units.Re=1.0`. Could allow `Re=0` sentinel when tau is explicit, or document requirement.

5. **BC parsing assumes exactly AMREX_SPACEDIM entries**
   - `bc.hiz_Uwall` style keys are per-face, but error messages about “wall normal must be zero” are cryptic for new users. Benchmark lid3d hit this when specifying MovingWall. Better docs/example for lid3d.

6. **bench.cpp ReduceRealMax semantics**
   - Using `ReduceRealMax(elapsed)` after `ParallelDescriptor::Barrier()` is correct for max elapsed, but if `USE_MPI=FALSE` the reduction is no-op. Document or add `USE_MPI` guard/conditional.

7. **Scaling helper: `grid.max_grid_size` tiling**
   - Fixed to `N` defeats AMReX load balance for multi-rank runs. Aurora runs need smaller `max_grid_size` (e.g., 32) for distribution. Suggest benchmark inputs use `max_grid_size=32` default and override for single-rank.

8. **Halfway bounce-back shifts diagonal links tangentially (mass leak, wrong walls) — found by `benchmarks/lid2d/validate`**
   - `src/LBM/BCTypes/BounceBack.H` fills a ghost cell `g` with `f(g - bound*e_dir, opp(n))`, i.e. the fluid cell at the *same tangential position*. Streaming is pull (`f_new(x,n) = f(x - c_n, n)`), so the fluid cell `x` that pulls population `n` from ghost `g = x - c_n` needs its *own* outgoing `opp(n)` population, which sits at `g + c_n`. For axis-aligned links `c_n` has no tangential part so the two agree; for diagonal links (all of D2Q9's corners, D3Q15/19/27's edge/corner links) the reflected population is shifted one cell along the wall. Interior rows still conserve mass, but at the ends of every wall row (and at domain corners, where the neighbouring ghost's matching population is never filled) populations are lost or taken from stale memory.
   - Symptom on a closed 2-D lid-driven cavity (N=64, Re=100): total mass drops by exactly 4/36 on step 1 (four corner links) and keeps leaking ~2e-5 of the mass per step; KE grows without bound and the lid power goes negative (flow diverges by ~step 10000).
   - One-line fix (verified in a scratch copy: mass conserved to ~3e-12 over 40000 steps, KE plateaus to 1e-13, Ghia Re=100 extrema within ~1.4% at N=64): in `BounceBack::apply`, read the source from `iv + c_n` instead of `iv - bound*e_dir`, for every dir (2-D and 3-D branches), e.g. 2-D `f(iv, n) = f(iv[0]+LatticeTag::cx(n), iv[1]+LatticeTag::cy(n), opp_n) + contrib(opp_n, rho, U);`. Same pattern in the 3-D branches (add `cz`). Also affects `lid3d`, `couette`, and every NoSlip/MovingWall/Inflow face.
   - Until fixed upstream and the submodule is bumped, `benchmarks/lid2d/validate` is expected to FAIL against the pinned `external/muphfasa`.
