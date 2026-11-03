#!/usr/bin/env python3
"""analyze_tgv.py — Compare TGV Re=1600 simulation to Dairay et al. 2017 DNS.

Reads AMReX plotfiles via yt, computes volume-averaged kinetic energy,
dissipation, and enstrophy at each output time, and compares to the
reference time series in TGV_Re1600.dat.

Outputs:
  - Console summary with L2 / peak errors and PASS/FAIL verdict
  - Comparison plots (Ek, dissipation, enstrophy vs time)
  - Optional CSV of the simulation-derived time series

Usage:
  # After running test_tgv3d (produces plt*/ in run_dir):
  python3 analyze_tgv.py --run-dir run --benchmark TGV_Re1600.dat

  # Use in-situ CSV instead of yt (faster, no yt needed):
  python3 analyze_tgv.py --csv run/tgv_diagnostics.csv --benchmark TGV_Re1600.dat

  # Custom thresholds:
  python3 analyze_tgv.py --run-dir run --peak-eps-tol 0.15 --ke-l2-tol 0.10

Requires: numpy, matplotlib, yt (for plotfile mode). yt is optional if
--csv is used.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

DEFAULT_TOL_PEAK_EPS = 0.20   # 20% error on peak dissipation admitted at coarse N
DEFAULT_TOL_KE_L2 = 0.10
DEFAULT_TOL_ENST_PEAK = 0.20

_SCRIPT_DIR = Path(__file__).resolve().parent

# ---------------------------------------------------------------------------
# Benchmark loader
# ---------------------------------------------------------------------------

def load_benchmark(path: Path):
    """Load TGV_Re1600.dat: skip comment lines, return dict of columns."""
    data = np.loadtxt(str(path))
    # Columns per header: 1:t 2:Ek 3:eps_t 4:eps 5:Dzeta 6-8:u2/v2/w2 9-17:grad2
    return {
        "t": data[:, 0],
        "Ek": data[:, 1],
        "eps_t": data[:, 2],
        "eps": data[:, 3],
        "enstrophy": data[:, 4],
        "data": data,
    }


def load_sim_csv(path: Path):
    header = Path(path).read_text().splitlines()[0].strip() if Path(path).exists() else ""
    arr = np.loadtxt(str(path), delimiter=",", skiprows=1)
    if arr.ndim == 1:
        arr = arr[None, :]
    # New CSV: step,time_phys,time,Ek,eps,enstrophy  (6 cols, time is bench-normalized)
    # Old CSV: step,time,Ek,eps,enstrophy (5 cols)
    if arr.shape[1] >= 6:
        # prefer bench time if present
        t = arr[:, 2] if "time_phys" in header else arr[:, 1]
        return {"t": t, "Ek": arr[:, 3], "eps": arr[:, 4], "enstrophy": arr[:, 5],
                "step": arr[:, 0].astype(int), "t_phys": arr[:, 1]}
    return {
        "t": arr[:, 1],
        "Ek": arr[:, 2],
        "eps": arr[:, 3],
        "enstrophy": arr[:, 4],
        "step": arr[:, 0].astype(int),
    }


# ---------------------------------------------------------------------------
# yt-based diagnostics (from plotfiles)
# ---------------------------------------------------------------------------

def _yt_available() -> bool:
    try:
        import yt  # noqa: F401
        return True
    except ImportError:
        return False


def compute_yt_diagnostics(pltdirs: list[Path], nu_phys: float = 1.0 / 1600.0,
                           L_ref: float = 6.283185307179586, U0: float = 1.0):
    """For each plotfile, load via yt and compute Ek, eps, enstrophy.

    yt fields are in physical units (SimulationImpl::writeOutput converts).
    Gradients are central differences with periodic BC; dx from dataset.
    If L_ref != 1 the time/eps/enstrophy are converted to benchmark
    nondim (t_bench=t_phys*U0/L_ref, etc.) to match TGV_Re1600.dat.
    """
    import yt
    yt.set_log_level(40)

    results = {"t": [], "Ek": [], "eps": [], "enstrophy": []}

    for pdir in pltdirs:
        ds = yt.load(str(pdir))
        # physical domain and cell size
        dims = np.array(ds.domain_dimensions, dtype=float)
        left = np.array(ds.domain_left_edge, dtype=float)
        right = np.array(ds.domain_right_edge, dtype=float)
        L = right - left
        dx = float(L[0] / dims[0])
        t = float(ds.current_time)

        cg = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
        # AMReX plotfile vars: density, vx, vy, vz in physical units
        # yt frontend exposes them as ("boxlib", var); covering_grid returns with ghost handling
        try:
            u = np.asarray(cg["boxlib", "vx"])
            v = np.asarray(cg["boxlib", "vy"])
            w = np.asarray(cg["boxlib", "vz"])
        except Exception:
            # fallback: try "stream" or bare name
            u = np.asarray(cg["vx"])
            v = np.asarray(cg["vy"])
            w = np.asarray(cg["vz"])

        # yt covering_grid may return shape (nx, ny, nz) with Fortran order;
        # squeeze any extra dims and ensure 3-D
        # AMReX ordering is (z,y,x) or (x,y,z) depending; yt doc says (x,y,z) with cg.
        # We treat axis 0 as x for gradient purposes consistently using np.roll.
        # Verify shape matches dims
        # If 4-D (with component), squeeze; if 2-D slice, skip
        # Normal case: shape == (nx, ny, nz) or (nx, ny, nz, 1)
        # For safety, just work with numpy arrays regardless of axis labels —
        # periodic central diff is symmetric, and isotropic averaging is invariant.

        # Ensure 3-D
        if u.ndim == 4:
            u = u[:, :, :, 0]
            v = v[:, :, :, 0]
            w = w[:, :, :, 0]
        # covering_grid pads one level: dims may be off by 1; trim if needed
        # Use actual array shape for N
        n_cells = u.size

        # KE: 0.5 * mean(u^2+v^2+w^2)
        ke = 0.5 * float(np.mean(u**2 + v**2 + w**2))

        # Gradients: 9 components, central diff, periodic, physical dx
        # Determine axis order: yt covering_grid returns array indexed [x, y, z]
        # (consistent with AMReX Box ordering). We assume axis 0=x, 1=y, 2=z.
        inv_2dx = 1.0 / (2.0 * dx)

        def grad(field, axis):
            return (np.roll(field, -1, axis=axis) - np.roll(field, 1, axis=axis)) * inv_2dx

        dudx = grad(u, 0); dudy = grad(u, 1); dudz = grad(u, 2)
        dvdx = grad(v, 0); dvdy = grad(v, 1); dvdz = grad(v, 2)
        dwdx = grad(w, 0); dwdy = grad(w, 1); dwdz = grad(w, 2)

        g2_mean = float(np.mean(
            dudx**2 + dudy**2 + dudz**2 +
            dvdx**2 + dvdy**2 + dvdz**2 +
            dwdx**2 + dwdy**2 + dwdz**2
        ))
        eps = nu_phys * g2_mean

        # vorticity and enstrophy
        wx = dwdy - dvdz
        wy = dudz - dwdx
        wz = dvdx - dudy
        enst = 0.5 * float(np.mean(wx**2 + wy**2 + wz**2))

        # Convert physical -> benchmark nondim if L_ref != 1
        t_bench = t * U0 / L_ref if L_ref != 1.0 else t
        ke_bench = ke / (U0 * U0) if U0 != 1.0 else ke
        eps_bench = eps * L_ref / (U0 ** 3) if L_ref != 1.0 else eps
        enst_bench = enst * L_ref * L_ref / (U0 * U0) if L_ref != 1.0 else enst

        results["t"].append(t_bench)
        results["Ek"].append(ke_bench)
        results["eps"].append(eps_bench)
        results["enstrophy"].append(enst_bench)

    for k in results:
        results[k] = np.array(results[k])
    # sort by time
    order = np.argsort(results["t"])
    for k in results:
        results[k] = results[k][order]
    return results


def find_plotfiles(run_dir: Path):
    return sorted(
        (p for p in run_dir.glob("plt*") if re.fullmatch(r"plt\d+", p.name)),
        key=lambda p: int(p.name[3:]),
    )


# ---------------------------------------------------------------------------
# Comparison / plotting
# ---------------------------------------------------------------------------

def interpolate_benchmark(bench, t_sim):
    """Linear interpolation of benchmark onto simulation time grid."""
    return {
        "Ek": np.interp(t_sim, bench["t"], bench["Ek"]),
        "eps": np.interp(t_sim, bench["t"], bench["eps"]),
        "enstrophy": np.interp(t_sim, bench["t"], bench["enstrophy"]),
    }


def compute_errors(sim, bench_interp, bench):
    """Return dict of error metrics."""
    def l2(a, b):
        denom = np.sqrt(np.mean(b**2))
        if denom == 0:
            return float("nan")
        return float(np.sqrt(np.mean((a - b)**2)) / denom)

    def rel(a, b):
        return float(np.abs(a - b) / max(np.abs(b), 1e-30))

    # L2 over overlapping time window
    out = {}
    for key in ("Ek", "eps", "enstrophy"):
        out[f"{key}_l2"] = l2(sim[key], bench_interp[key])

    # Peak dissipation (key DNS target: timing + magnitude of eps max)
    i_sim_peak = int(np.argmax(sim["eps"]))
    i_bench_peak = int(np.argmax(bench["eps"]))
    out["peak_eps_sim"] = float(sim["eps"][i_sim_peak])
    out["peak_eps_bench"] = float(bench["eps"][i_bench_peak])
    out["peak_eps_t_sim"] = float(sim["t"][i_sim_peak])
    out["peak_eps_t_bench"] = float(bench["t"][i_bench_peak])
    out["peak_eps_rel_err"] = rel(out["peak_eps_sim"], out["peak_eps_bench"])
    out["peak_t_abs_err"] = float(abs(out["peak_eps_t_sim"] - out["peak_eps_t_bench"]))

    # Enstrophy peak similarly
    i_sim_e = int(np.argmax(sim["enstrophy"]))
    i_bench_e = int(np.argmax(bench["enstrophy"]))
    out["peak_enst_sim"] = float(sim["enstrophy"][i_sim_e])
    out["peak_enst_bench"] = float(bench["enstrophy"][i_bench_e])
    out["peak_enst_rel_err"] = rel(out["peak_enst_sim"], out["peak_enst_bench"])

    return out


def make_plots(sim, bench, bench_interp, out_prefix: Path, show: bool = False):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # --- Ek(t) ---
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(bench["t"], bench["Ek"], label="DNS (Dairay 2017)", color="k", lw=1.2)
    ax.plot(sim["t"], sim["Ek"], "o-", ms=3, label="muphfasa", color="C0")
    ax.set_xlabel("t")
    ax.set_ylabel(r"$E_k = \frac12\langle|u|^2\rangle$")
    ax.set_title("TGV Re=1600: kinetic energy")
    ax.legend()
    fig.tight_layout()
    fig.savefig(str(out_prefix) + "_Ek.png", dpi=150)
    plt.close(fig)
    print(f"Wrote {out_prefix}_Ek.png")

    # --- dissipation ---
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(bench["t"], bench["eps"], label="DNS", color="k", lw=1.2)
    ax.plot(sim["t"], sim["eps"], "o-", ms=3, label="muphfasa", color="C1")
    ax.set_xlabel("t")
    ax.set_ylabel(r"$\epsilon = \nu\langle|\nabla u|^2\rangle$")
    ax.set_title("TGV Re=1600: dissipation")
    ax.legend()
    fig.tight_layout()
    fig.savefig(str(out_prefix) + "_dissipation.png", dpi=150)
    plt.close(fig)
    print(f"Wrote {out_prefix}_dissipation.png")

    # --- enstrophy ---
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(bench["t"], bench["enstrophy"], label="DNS", color="k", lw=1.2)
    ax.plot(sim["t"], sim["enstrophy"], "o-", ms=3, label="muphfasa", color="C2")
    ax.set_xlabel("t")
    ax.set_ylabel(r"$\mathcal{E} = \frac12\langle|\omega|^2\rangle$")
    ax.set_title("TGV Re=1600: enstrophy")
    ax.legend()
    fig.tight_layout()
    fig.savefig(str(out_prefix) + "_enstrophy.png", dpi=150)
    plt.close(fig)
    print(f"Wrote {out_prefix}_enstrophy.png")

    # --- combined ---
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), sharex=True)
    for a, key, title, col in [
        (axes[0], "Ek", r"$E_k$", "C0"),
        (axes[1], "eps", r"$\epsilon$", "C1"),
        (axes[2], "enstrophy", r"$\mathcal{E}$", "C2"),
    ]:
        a.plot(bench["t"], bench[key], color="k", lw=1.2, label="DNS" if a == axes[0] else None)
        a.plot(sim["t"], sim[key], "o-", ms=2, color=col, label="muphfasa" if a == axes[0] else None)
        a.set_title(title)
        a.set_xlabel("t")
    axes[0].set_ylabel("value")
    axes[0].legend(fontsize=8)
    fig.suptitle("TGV Re=1600", fontsize=11)
    fig.tight_layout()
    fig.savefig(str(out_prefix) + "_combined.png", dpi=150)
    plt.close(fig)
    print(f"Wrote {out_prefix}_combined.png")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-dir", type=Path, default=None,
                   help="directory containing plt*/ plotfiles")
    p.add_argument("--csv", type=Path, default=None,
                   help="path to tgv_diagnostics.csv (in-situ, no yt needed)")
    p.add_argument("--benchmark", type=Path, default=_SCRIPT_DIR / "TGV_Re1600.dat",
                   help="path to TGV_Re1600.dat")
    p.add_argument("--nu", type=float, default=1.0 / 1600.0,
                   help="physical viscosity (default 1/1600, benchmark nondim)")
    p.add_argument("--Lref", type=float, default=1.0,
                   help="reference length for bench scaling (default 1.0)")
    p.add_argument("--U0", type=float, default=1.0,
                   help="reference velocity for bench scaling")
    p.add_argument("--out-prefix", type=Path, default=None,
                   help="prefix for output plots/CSVs (default: <run-dir>/tgv_compare)")
    p.add_argument("--peak-eps-tol", type=float, default=DEFAULT_TOL_PEAK_EPS,
                   help="tolerance on |eps_peak_sim - eps_peak_dns|/eps_peak_dns")
    p.add_argument("--ke-l2-tol", type=float, default=DEFAULT_TOL_KE_L2,
                   help="tolerance on L2(Ek) relative error")
    p.add_argument("--plot", action="store_true", default=True,
                   help="generate comparison plots (default: on)")
    p.add_argument("--no-plot", dest="plot", action="store_false",
                   help="skip plot generation")
    p.add_argument("--save-csv", type=Path, default=None,
                   help="write simulation time series to CSV")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    # Resolve benchmark path (allow relative to script dir)
    bench_path = args.benchmark
    if not bench_path.is_file():
        alt = _SCRIPT_DIR / bench_path.name
        if alt.is_file():
            bench_path = alt
    if not bench_path.is_file():
        print(f"Benchmark file not found: {args.benchmark}", file=sys.stderr)
        return 2

    bench = load_benchmark(bench_path)
    print(f"Loaded benchmark {bench_path}: {len(bench['t'])} steps t=0..{bench['t'][-1]:.1f} "
          f"peak eps={bench['eps'].max():.5f} at t={bench['t'][np.argmax(bench['eps'])]:.2f}")

    # Load simulation data
    sim = None
    if args.csv is not None:
        if not args.csv.is_file():
            # try relative to run_dir
            if args.run_dir and (args.run_dir / args.csv.name).is_file():
                args.csv = args.run_dir / args.csv.name
            elif (Path.cwd() / args.csv).is_file():
                args.csv = Path.cwd() / args.csv
            else:
                print(f"CSV not found: {args.csv}", file=sys.stderr)
                return 2
        sim = load_sim_csv(args.csv)
        print(f"Loaded simulation CSV {args.csv}: {len(sim['t'])} samples "
              f"t={sim['t'][0]:.3f}..{sim['t'][-1]:.3f}")
    elif args.run_dir is not None:
        if not args.run_dir.is_dir():
            print(f"run-dir not found: {args.run_dir}", file=sys.stderr)
            return 2
        # prefer CSV inside run_dir if present (faster, no yt)
        csv_in_run = args.run_dir / "tgv_diagnostics.csv"
        if csv_in_run.is_file():
            print(f"Found in-situ CSV {csv_in_run}; using it (pass --csv to force yt mode)")
            sim = load_sim_csv(csv_in_run)
            print(f"Loaded simulation CSV: {len(sim['t'])} samples")
        else:
            if not _yt_available():
                print("yt not installed and no tgv_diagnostics.csv found. "
                      "Install yt or run with --csv.", file=sys.stderr)
                return 2
            plts = find_plotfiles(args.run_dir)
            if not plts:
                print(f"No plt*/ found in {args.run_dir}", file=sys.stderr)
                return 2
            print(f"Found {len(plts)} plotfiles; computing diagnostics via yt "
                  f"(nu={args.nu:.6f}, Lref={args.Lref:.4f})...")
            # yt physical nu = U*L/Re; for bench scaling pass physical nu then convert inside
            nu_phys_plot = args.U0 * args.Lref / 1600.0 if args.Lref != 1.0 else args.nu
            sim = compute_yt_diagnostics(plts, nu_phys=nu_phys_plot, L_ref=args.Lref, U0=args.U0)
            print(f"Diagnostics: t={sim['t'][0]:.3f}..{sim['t'][-1]:.3f}")
    else:
        # auto-detect: look for CSV or plt* in CWD or script dir
        for cand in [Path.cwd(), _SCRIPT_DIR]:
            if (cand / "tgv_diagnostics.csv").is_file():
                sim = load_sim_csv(cand / "tgv_diagnostics.csv")
                print(f"Auto-found {cand / 'tgv_diagnostics.csv'}")
                break
            if find_plotfiles(cand):
                if _yt_available():
                    plts = find_plotfiles(cand)
                    print(f"Auto-found {len(plts)} plotfiles in {cand}")
                    nu_phys_plot = args.U0 * args.Lref / 1600.0 if args.Lref != 1.0 else args.nu
                    sim = compute_yt_diagnostics(plts, nu_phys=nu_phys_plot, L_ref=args.Lref, U0=args.U0)
                    break
        if sim is None:
            print("No simulation data found. Pass --run-dir or --csv.", file=sys.stderr)
            return 2

    # Restrict benchmark interp to simulation time window
    bench_interp = interpolate_benchmark(bench, sim["t"])
    errors = compute_errors(sim, bench_interp, bench)

    print("\n--- Comparison ---")
    print(f"  KE L2 rel err:       {errors['Ek_l2']:.3e}")
    print(f"  eps L2 rel err:      {errors['eps_l2']:.3e}")
    print(f"  enstrophy L2 rel err:{errors['enstrophy_l2']:.3e}")
    print(f"  peak eps: sim={errors['peak_eps_sim']:.5f} at t={errors['peak_eps_t_sim']:.2f}  "
          f"dns={errors['peak_eps_bench']:.5f} at t={errors['peak_eps_t_bench']:.2f}  "
          f"rel_err={errors['peak_eps_rel_err']:.3e}  dt_err={errors['peak_t_abs_err']:.2f}")
    print(f"  peak enstrophy: sim={errors['peak_enst_sim']:.5f} dns={errors['peak_enst_bench']:.5f} "
          f"rel_err={errors['peak_enst_rel_err']:.3e}")

    # Verdict
    ok = True
    msgs = []
    if errors["peak_eps_rel_err"] > args.peak_eps_tol:
        msgs.append(f"FAIL peak eps rel err {errors['peak_eps_rel_err']:.3e} > {args.peak_eps_tol}")
        ok = False
    else:
        msgs.append(f"PASS peak eps (tol {args.peak_eps_tol})")
    if errors["Ek_l2"] > args.ke_l2_tol:
        msgs.append(f"FAIL KE L2 {errors['Ek_l2']:.3e} > {args.ke_l2_tol}")
        ok = False
    else:
        msgs.append(f"PASS KE L2 (tol {args.ke_l2_tol})")
    for m in msgs:
        print(m)
    print("OVERALL:", "PASS" if ok else "FAIL")

    # Save CSV if requested
    if args.save_csv:
        arr = np.column_stack([sim["t"], sim["Ek"], sim["eps"], sim["enstrophy"]])
        header = "t,Ek,eps,enstrophy"
        np.savetxt(str(args.save_csv), arr, delimiter=",", header=header, comments="")
        print(f"Wrote {args.save_csv}")

    # Plots
    if args.plot:
        out_prefix = args.out_prefix
        if out_prefix is None:
            base = args.run_dir if args.run_dir else Path.cwd()
            # if base is a file path (csv), use its parent
            if base.is_file():
                base = base.parent
            out_prefix = base / "tgv_compare"
            # also allow script dir default
            if not base.is_dir():
                out_prefix = Path("tgv_compare")
        try:
            make_plots(sim, bench, bench_interp, out_prefix)
        except Exception as e:
            print(f"Plotting failed: {e}", file=sys.stderr)

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
