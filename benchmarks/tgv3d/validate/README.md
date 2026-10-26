# TGV 3-D correctness check

Minimal validation counterpart to `../inputs.bench` (perf-only): triply-periodic domain
seeded with the classic 3-D TGV IC (`../Tgv3dAnalytic.H`, via `../InitData.cpp`).

**Unlike `tgv2d`/`beltrami3d`, this isn't an `fcompare`-vs-analytic test.** The classic
3-D TGV field has no closed-form solution for t>0 — vortex stretching moves energy into
other Fourier modes (see the comment in `../Tgv3dAnalytic.H`) — so there's no analytic
field to diff against. Instead `test_tgv3d.cpp` checks the **global kinetic-energy
balance**, an identity that *is* exact for incompressible Navier-Stokes regardless of
vortex stretching: nonlinear advection and pressure do no net work on a periodic,
divergence-free field, so kinetic energy can only drop via viscous dissipation —

```
KE(0) - KE(T) = nu * integral_0^T ( integral |curl(u)|^2 dV ) dt
```

At every `amr.plot_int` step it computes `KE(t) = (1/2) sum|u|^2` and
`W(t) = sum|curl(u)|^2` (lattice units, central differences with periodic ghost fill),
then at the end compares the measured energy drop `KE(0)-KE(T)` against
`nu * trapz(W, t)`. Unlike the old pointwise/linearized-decay check, this holds at any
Re/amplitude/duration — it's a genuine correctness check on the solver, not just on the
IC wiring, and residual error comes only from the trapezoidal time integration and the
solver's own truncation error.

## Build and run
```sh
make                             # builds test_tgv3d3d.gnu.ex (DIM=3)
mkdir -p run && cd run
../test_tgv3d3d.gnu.ex ../inputs
```
Prints `KE`/`|curl(u)|^2` at each `amr.plot_int` step, then the final energy-balance
relative error.

## Expected result
~1% relative error between measured KE drop and `nu * integral(|curl(u)|^2) dt`
(D3Q19/tau=0.8/N=32, `amr.plot_int=10` → 21 samples over 200 steps). This is purely
trapezoidal-rule error from under-sampling the (rapidly decaying) `|curl(u)|^2` curve —
not an approximation to the physics — so it shrinks with finer `amr.plot_int` (e.g.
~9% at `plot_int=50`, ~1.5% at `plot_int=5`) and should hold at any Re/amplitude/run
length, unlike the old linearized-decay check.
