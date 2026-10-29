#!/usr/bin/env python3
"""run_scaling.py — Weak/strong scaling sweeps for benchmarks.

Batch mode (default): delegates to batch/generate.py to emit Aurora PBS
scripts (one per sweep point). Submission is then
`qsub batch/generated/<case>_<size>_<nodes>nodes.sh`. This is the preferred
path for Aurora.

Local mode (--local): builds and runs via helper.run_bench (mpiexec), parses
BENCH_RESULT from bench.cpp (MLUPS, ms/step) and plots efficiency. Useful for
quick dev on a workstation or single-node test.

Usage (batch/Aurora, default):
  python3 run_scaling.py --case lid2d --sizes 1024 2048 4096 --nodes 1 2 4 --dry-run
  python3 run_scaling.py --case lid3d --sizes 128 256 512 --lattice D3Q19 --queue debug-scaling
  python3 run_scaling.py --case tgv2d --sizes 2048 4096 --nodes 2 4 --stages performance --walltime 00:30:00
  python3 run_scaling.py --case lid2d --sizes 2048 --nodes 2  # full pipeline, writes batch/generated/

Usage (local):
  python3 run_scaling.py --case tgv2d --mode weak --sizes 1024 2048 4096 --nodes 1 2 4 --ppn 1 --local
  python3 run_scaling.py --case lid2d --mode weak --sizes 512 1024 2048 --nodes 1 2 4 --ppn 1 --local
  python3 run_scaling.py --case lid3d --mode weak --sizes 64 128 256 --lattice D3Q19 --plot scaling.png --local
"""

import argparse
import itertools
import subprocess
import sys
from pathlib import Path

from helper import BENCH_DIR, build_bench, find_bench_exe, parse_bench_result, run_bench

CASES = {
    "tgv2d": {"dim": 2, "default_lattice": "D2Q9", "bench_dir": "benchmarks/tgv2d"},
    "tgv3d": {"dim": 3, "default_lattice": "D3Q19", "bench_dir": "benchmarks/tgv3d"},
    "lid2d": {"dim": 2, "default_lattice": "D2Q9", "bench_dir": "benchmarks/lid2d"},
    "lid3d": {"dim": 3, "default_lattice": "D3Q19", "bench_dir": "benchmarks/lid3d"},
    "beltrami3d": {"dim": 3, "default_lattice": "D3Q19", "bench_dir": "benchmarks/beltrami3d"},
    # Validate-only (no bench.cpp) — batch/generate.py handles them for
    # the validating stage; use default batch mode to emit their scripts.
    "couette2d": {"dim": 2, "default_lattice": "D2Q9", "bench_dir": None},
    "couette3d": {"dim": 3, "default_lattice": "D3Q19", "bench_dir": None},
    "poiseuille2d": {"dim": 2, "default_lattice": "D2Q9", "bench_dir": None},
    "poiseuille3d": {"dim": 3, "default_lattice": "D3Q19", "bench_dir": None},
    "uniform2d": {"dim": 2, "default_lattice": "D2Q9", "bench_dir": None},
    "uniform3d": {"dim": 3, "default_lattice": "D3Q19", "bench_dir": None},
}


def _run_batch_generate(args, dim: int, lattice: str):
    """Delegate sweep to batch/generate.py: one PBS script per (size, nodes) pair.

    Mirrors the local zip-cycle pairing (sizes[i] ↔ nodes[i%len(nodes)]) so
    `run_scaling.py --sizes A B C --nodes 1 2 4` and
    `run_scaling.py --sizes A B C --nodes 1 2 4 --local` (local mode) produce
    the same three sweep points (2048 on 1 node, 4096 on 2, …). Each point is
    emitted via a separate generate.py invocation ( --case X --nodes N --size S )
    which is the most explicit and avoids product/zip ambiguity in the generator.
    """
    gen = BENCH_DIR / "batch" / "generate.py"
    if not gen.is_file():
        print(f"batch generator not found: {gen}", file=sys.stderr)
        sys.exit(1)

    pairs = list(zip(args.sizes, itertools.islice(itertools.cycle(args.nodes), len(args.sizes))))
    print(f"Batch mode: {args.case} dim={dim} lattice={lattice} ppn={args.ppn} queue={args.queue or 'auto'}")
    print(f"  sweep pairs (size ↔ nodes): {pairs}")
    if args.max_grid_size != 32:
        print("  note: --max-grid-size is a local-run AMReX box tiling knob; batch scripts use the bench default", file=sys.stderr)

    written = []
    for size, nodes in pairs:
        cmd = [
            sys.executable, str(gen),
            "--case", args.case,
            "--nodes", str(nodes),
            "--size", str(size),
            "--ppn", str(args.ppn),
        ]
        if lattice:
            cmd += ["--lattice", lattice]
        if args.queue:
            cmd += ["--queue", args.queue]
        if args.walltime:
            cmd += ["--walltime", args.walltime]
        if args.project:
            cmd += ["--project", args.project]
        if args.stages:
            cmd += ["--stages"] + args.stages
        if args.warmup is not None:
            cmd += ["--warmup", str(args.warmup)]
        if args.steps is not None:
            cmd += ["--steps", str(args.steps)]
        if args.flare_root:
            cmd += ["--flare-root", args.flare_root]
        if args.outdir:
            cmd += ["--outdir", str(args.outdir)]
        if args.dry_run:
            cmd.append("--dry-run")

        print(f"  $ {' '.join(cmd)}")
        result = subprocess.run(cmd)
        if result.returncode != 0:
            print(f"generate failed for size={size} nodes={nodes} (exit {result.returncode})", file=sys.stderr)
            sys.exit(result.returncode)
        written.append((size, nodes))

    if not args.dry_run:
        print(f"\nBatch scripts ready for {len(written)} sweep point(s). Submit with e.g.:")
        for size, nodes in written:
            tag = ""
            if args.stages and set(args.stages) != {"compiling", "validating", "performance"}:
                tag = "_" + "_".join(sorted(args.stages, key=lambda s: ["compiling", "validating", "performance"].index(s)))
            print(f"  qsub batch/generated/{args.case}_{size}_{nodes}nodes{tag}.sh")
        print(f"  or: qsub batch/generated/{args.case}_*.sh  (all)")
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--case", choices=CASES.keys(), required=True)
    parser.add_argument("--mode", choices=["weak", "strong"], default="weak")
    parser.add_argument("--sizes", type=int, nargs="+", required=True, help="N per side (e.g. 1024 2048)")
    parser.add_argument("--nodes", type=int, nargs="+", default=[1], help="MPI ranks or nodes")
    parser.add_argument("--ppn", type=int, default=1,
                        help="MPI ranks per node (default 1; Aurora is 12). "
                             "nprocs for each sweep point = nodes*ppn, actually launched via mpiexec.")
    parser.add_argument("--max-grid-size", type=int, default=32,
                        help="AMReX box tile size (default 32). Kept independent of --sizes so the "
                             "domain is decomposed into multiple boxes for multi-rank runs. Local mode only; "
                             "batch mode uses the bench default.")
    parser.add_argument("--lattice", type=str, default=None)
    parser.add_argument("--plot", type=Path, default=None)
    parser.add_argument("--no-build", action="store_true", help="local only: skip build step")
    parser.add_argument("--no-log", action="store_true", help="local only: skip perf DB logging")
    parser.add_argument("--force-log", action="store_true", help="local only: bypass dedup for perf log")
    # Execution mode: batch (default) vs local
    parser.add_argument("--local", action="store_true",
                        help="run locally via helper.run_bench (mpiexec) instead of generating Aurora PBS scripts "
                             "(default is batch mode via batch/generate.py)")
    parser.add_argument("--queue", type=str, default=None, help="PBS queue (default: auto from nodes)")
    parser.add_argument("--walltime", type=str, default=None, help="PBS walltime (default 00:30:00)")
    parser.add_argument("--project", type=str, default=None, help="ALCF project (default muphfasa)")
    parser.add_argument("--stages", nargs="+", choices=["compiling", "validating", "performance"], default=None,
                        help="stages to include (default: all three; e.g. --stages performance)")
    parser.add_argument("--warmup", type=int, default=None, help="warmup steps (default per dim: 2D=100, 3D=50)")
    parser.add_argument("--steps", type=int, default=None, help="timed steps (default per dim: 2D=2000, 3D=500)")
    parser.add_argument("--flare-root", type=str, default=None, help="flare run root")
    parser.add_argument("--outdir", type=Path, default=None, help="output dir (default batch/generated/)")
    parser.add_argument("--dry-run", action="store_true", help="print what would be generated without writing (batch mode)")
    args = parser.parse_args()

    cfg = CASES[args.case]
    dim = cfg["dim"]
    lattice = args.lattice or cfg["default_lattice"]

    # Batch path (default): emit PBS scripts and exit (no local build/run)
    if not args.local:
        return _run_batch_generate(args, dim, lattice)

    bench_dir = BENCH_DIR / cfg["bench_dir"] if cfg["bench_dir"] else None
    if bench_dir is None or not bench_dir.is_dir():
        print(f"case '{args.case}' has no local bench (validate-only); remove --local to generate its validate scripts in batch mode", file=sys.stderr)
        sys.exit(1)
    inputs = bench_dir / "inputs.bench"

    if not args.no_build:
        build_bench(args.case, dim, use_mpi=True)

    results = []
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
            ncells_override = next((o.split("=", 1)[1] for o in overrides if o.startswith("grid.ncells=")), "")
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
