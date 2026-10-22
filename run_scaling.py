#!/usr/bin/env python3
"""run_scaling.py — Weak/strong scaling sweeps for benchmarks.

Parses BENCH_RESULT from bench.cpp (MLUPS, ms/step) and plots efficiency.

Usage:
  python3 run_scaling.py --case tgv2d --mode weak --sizes 1024 2048 4096 --nodes 1 2 4
  python3 run_scaling.py --case lid3d --mode weak --sizes 64 128 256 --lattice D3Q19 --plot scaling.png
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

from helper import build_bench, parse_bench_result, run_bench

CASES = {
    "tgv2d": {"dim": 2, "default_lattice": "D2Q9", "bench_dir": "benchmarks/tgv2d"},
    "lid3d": {"dim": 3, "default_lattice": "D3Q19", "bench_dir": "benchmarks/lid3d"},
    "beltrami3d": {"dim": 3, "default_lattice": "D3Q19", "bench_dir": "benchmarks/beltrami3d"},
}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=CASES.keys(), required=True)
    parser.add_argument("--mode", choices=["weak","strong"], default="weak")
    parser.add_argument("--sizes", type=int, nargs="+", required=True, help="N per side (e.g. 1024 2048)")
    parser.add_argument("--nodes", type=int, nargs="+", default=[1], help="MPI ranks or nodes")
    parser.add_argument("--lattice", type=str, default=None)
    parser.add_argument("--plot", type=Path, default=None)
    parser.add_argument("--no-build", action="store_true")
    args = parser.parse_args()

    cfg = CASES[args.case]
    dim = cfg["dim"]
    lattice = args.lattice or cfg["default_lattice"]
    from helper import BENCH_DIR
    bench_dir = BENCH_DIR / cfg["bench_dir"]
    inputs = bench_dir / "inputs.bench"

    if not args.no_build:
        build_bench(args.case, dim, use_mpi=True)

    # For now, just run locally with 1 rank per size (Aurora batch handles MPI)
    # Weak: size scales with nodes, Strong: fixed size, vary nodes
    results = []
    import itertools
    nodes_cycle = itertools.islice(itertools.cycle(args.nodes), len(args.sizes))
    for size, nodes in zip(args.sizes, nodes_cycle):
        # Build overrides: grid.ncells and lattice
        if dim == 2:
            overrides = [f"grid.ncells={size} {size}", f"grid.max_grid_size={size}", f"lbm.lattice_type={lattice}"]
        else:
            overrides = [f"grid.ncells={size} {size} {size}", f"grid.max_grid_size={size}", f"lbm.lattice_type={lattice}"]
        # For MPI, we would do mpiexec -n nodes*ppn, but locally just 1 rank
        solver_exe = bench_dir / f"bench{dim}d.gnu.ex" if (bench_dir / f"bench{dim}d.gnu.ex").exists() else find_bench_exe(bench_dir)
        run_dir = Path(f"scaling_{args.case}_{size}_{nodes}")
        log = run_bench(solver_exe, inputs, run_dir, overrides)
        res = parse_bench_result(log)
        res["size"] = size
        res["nodes"] = nodes
        results.append(res)
        print(f"size={size} nodes={nodes} MLUPS={res['mlups']:.1f} ms/step={res['ms_per_step']:.2f}")

    if args.plot and results:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        sizes = [r["size"] for r in results]
        mlups = [r["mlups"] for r in results]
        plt.figure()
        plt.loglog(sizes, mlups, "o-")
        plt.xlabel("N")
        plt.ylabel("MLUPS")
        plt.title(f"{args.case} {args.mode} scaling ({lattice})")
        plt.savefig(args.plot, dpi=150)
        print(f"Wrote {args.plot}")

def find_bench_exe(bench_dir: Path) -> Path:
    exes = sorted(bench_dir.glob("bench*.ex")) + sorted(bench_dir.glob("*.ex"))
    if not exes:
        raise FileNotFoundError(f"no bench exe in {bench_dir}")
    return exes[0]

if __name__ == "__main__":
    sys.exit(main())
