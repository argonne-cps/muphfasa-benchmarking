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
