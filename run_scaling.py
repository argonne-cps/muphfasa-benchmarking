#!/usr/bin/env python3
"""run_scaling.py — Weak/strong scaling sweeps for benchmarks.

Parses BENCH_RESULT from bench.cpp (MLUPS, ms/step) and plots efficiency.

Each sweep point actually launches nprocs = nodes*ppn ranks via mpiexec (see
helper.run_bench); --max-grid-size is fixed independent of --sizes so the
domain is decomposed into enough boxes for those ranks to do real work.

Usage:
  python3 run_scaling.py --case tgv2d --mode weak --sizes 1024 2048 4096 --nodes 1 2 4 --ppn 1
  python3 run_scaling.py --case lid3d --mode weak --sizes 64 128 256 --lattice D3Q19 --plot scaling.png
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

from helper import BENCH_DIR, build_bench, find_bench_exe, parse_bench_result, run_bench

CASES = {
    "tgv2d": {"dim": 2, "default_lattice": "D2Q9", "bench_dir": "benchmarks/tgv2d"},
    "tgv3d": {"dim": 3, "default_lattice": "D3Q19", "bench_dir": "benchmarks/tgv3d"},
    "lid3d": {"dim": 3, "default_lattice": "D3Q19", "bench_dir": "benchmarks/lid3d"},
    "beltrami3d": {"dim": 3, "default_lattice": "D3Q19", "bench_dir": "benchmarks/beltrami3d"},
}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=CASES.keys(), required=True)
    parser.add_argument("--mode", choices=["weak","strong"], default="weak")
    parser.add_argument("--sizes", type=int, nargs="+", required=True, help="N per side (e.g. 1024 2048)")
    parser.add_argument("--nodes", type=int, nargs="+", default=[1], help="MPI ranks or nodes")
    parser.add_argument("--ppn", type=int, default=1,
                         help="MPI ranks per node (default 1; Aurora is 12). "
                              "nprocs for each sweep point = nodes*ppn, actually launched via mpiexec.")
    parser.add_argument("--max-grid-size", type=int, default=32,
                         help="AMReX box tile size (default 32). Kept independent of --sizes so the "
                              "domain is actually decomposed into multiple boxes for multi-rank runs "
                              "-- setting this equal to the domain size (the old default) produces a "
                              "single box that can only run on one rank regardless of --nodes/--ppn.")
    parser.add_argument("--lattice", type=str, default=None)
    parser.add_argument("--plot", type=Path, default=None)
    parser.add_argument("--no-build", action="store_true")
    parser.add_argument("--no-log", action="store_true", help="skip perf DB logging")
    parser.add_argument("--force-log", action="store_true", help="bypass dedup for perf log")
    args = parser.parse_args()

    cfg = CASES[args.case]
    dim = cfg["dim"]
    lattice = args.lattice or cfg["default_lattice"]
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
        nprocs = max(1, nodes) * max(1, args.ppn)
        # Build overrides: grid.ncells and lattice. max_grid_size is fixed and
        # independent of the domain size so multi-rank runs actually get more
        # than one box to distribute (see --max-grid-size help).
        if dim == 2:
            overrides = [f"grid.ncells={size} {size}", f"grid.max_grid_size={args.max_grid_size}", f"lbm.lattice_type={lattice}"]
        else:
            overrides = [f"grid.ncells={size} {size} {size}", f"grid.max_grid_size={args.max_grid_size}", f"lbm.lattice_type={lattice}"]
        solver_exe = find_bench_exe(bench_dir, dim, use_mpi=(nprocs > 1))
        run_dir = Path(f"scaling_{args.case}_{size}_{nodes}")
        log = run_bench(solver_exe, inputs, run_dir, overrides, nprocs=nprocs, ppn=args.ppn)
        res = parse_bench_result(log)
        res["size"] = size
        res["nodes"] = nodes
        res["nprocs"] = nprocs
        results.append(res)
        print(f"size={size} nodes={nodes} nprocs={nprocs} MLUPS={res['mlups']:.1f} ms/step={res['ms_per_step']:.2f}")
        if not args.no_log:
            # Log to perf DB: test_problem, lattice, ncells from overrides, real nprocs
            ncells_override = next((o.split("=",1)[1] for o in overrides if o.startswith("grid.ncells=")), "")
            log_path = run_dir / "run.log"
            try:
                cmd = [sys.executable, str(BENCH_DIR / "perf.py"), "log", "--test", args.case,
                       "--lattice", lattice, "--ncells", ncells_override, "--nprocs", str(nprocs),
                       "--ppn", str(args.ppn), "--nnodes", str(nodes),
                       "--inputs", str(inputs), "--log", str(log_path)]
                if args.force_log:
                    cmd.append("--force")
                subprocess.run(cmd, check=False)
            except Exception as e:
                print(f"perf log failed: {e}", file=sys.stderr)

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

if __name__ == "__main__":
    sys.exit(main())
