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

Restart the timed region from a checkpoint (e.g. one written by --save-checkpoint):
  python3 run_scaling.py --case tgv3d --sizes 1024 --nodes 4 --stages performance \\
      --restart /flare/proj/$USER/muphfasa-benchmarking-runs/tgv3d_1024_4nodes_123/chk0000550
  # --warmup/--steps then count from the checkpoint's step; the size must match the checkpoint.

Queue limits (batch mode): every sweep point's node count and --walltime are checked against the Aurora
queue it will go to (--queue, else chosen from the node count: debug 1-2, debug-scaling 2-256, prod 256+;
prod's walltime cap depends on the node count). All violations are listed and nothing is generated:
  python3 run_scaling.py --case tgv3d --sizes 1024 --nodes 4 --walltime 02:00:00   # debug-scaling max is 1 h
  python3 run_scaling.py --case tgv3d --sizes 1024 --nodes 4 --walltime 02:00:00 --queue capacity   # ok
Pass --no-queue-check to skip (e.g. after ALCF changes the limits; edit QUEUE_LIMITS in batch/generate.py).

Check the settings with muphfasa's input_info (util/input_info, must be built there):
  python3 run_scaling.py --case tgv3d --mode strong --sizes 1024 --nodes 1 2 4 8 --stages performance \\
      --dry-run --input-info
  # Prints input_info's full report for the first sweep point (grid, dt/tau/u_lat, steps, disk) with
  # the same overrides the run would use, then a one-line memory check per point against the 768 GiB
  # of GPU memory per Aurora node. Works with or without --dry-run, in batch and --local mode.
  # Build the tool first:  make -C external/muphfasa/util/input_info DIM=3   (DIM=2 for 2-D cases)

Usage (local):
  python3 run_scaling.py --case tgv2d --mode weak --sizes 1024 2048 4096 --nodes 1 2 4 --ppn 1 --local
  python3 run_scaling.py --case lid2d --mode weak --sizes 512 1024 2048 --nodes 1 2 4 --ppn 1 --local
  python3 run_scaling.py --case lid3d --mode weak --sizes 64 128 256 --lattice D3Q19 --plot scaling.png --local
"""

import argparse
import itertools
import os
import re
import subprocess
import sys
from pathlib import Path

from batch import generate as batch_generate
from helper import (BENCH_DIR, INPUT_INFO_DIR, build_bench, drop_bench_keys, find_bench_exe, find_input_info,
                    bench_keys_read, node_memory_check, parse_bench_result, run_bench, run_input_info,
                    unknown_bench_keys)

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


def _bench_overrides(args, dim: int, size: int, lattice: str, restart: str, local: bool) -> list:
    """key=value settings applied on top of benchmarks/<case>/inputs.bench for one sweep point.

    local: what run_scaling.py --local passes to the exe. Otherwise: what the batch script written by
    batch/generate.py puts in the run directory's inputs file (keep in step with its performance block;
    the script additionally sets sim.max_step, which bench.cpp ignores when bench.steps is given)."""
    ov = [f"grid.ncells={' '.join([str(size)] * dim)}"]
    if local:
        ov.append(f"grid.max_grid_size={args.max_grid_size}")
    ov.append(f"lbm.lattice_type={lattice}")
    if not local:
        ov.append("amr.plot_int=-1")
    if args.save_checkpoint:
        ov.append("bench.save_checkpoint=1")
    elif not local:
        ov.append("bench.save_checkpoint=0")
    if restart:
        ov.append(f"amr.restart={restart}")
    # bench.steps timed steps after bench.warmup warmup steps, counted from the start (or checkpoint)
    # step; locally, without these the inputs file's sim.max_step/bench.warmup apply.
    if not local or restart or args.warmup is not None or args.steps is not None:
        ov.append(f"bench.warmup={args.warmup if args.warmup is not None else batch_generate.DEFAULT_WARMUP[dim]}")
        ov.append(f"bench.steps={args.steps if args.steps is not None else batch_generate.DEFAULT_STEPS[dim]}")
    return ov


def _known_bench_keys(args) -> set:
    """bench.* keys the performance driver reads: benchmarks/bench.cpp and the case's main.cpp."""
    bench_rel = CASES[args.case]["bench_dir"]
    return bench_keys_read(BENCH_DIR / "benchmarks" / "bench.cpp", BENCH_DIR / bench_rel / "main.cpp")


def _check_bench_keys(args, dim, lattice, size, restart, local, inputs) -> int:
    """Warn about bench.* keys in the inputs file / overrides that the driver never reads.

    input_info cannot judge these (the driver, not the solver, reads them) and they are hidden from
    its "never reads" note, so a typo such as bench.warmpu would otherwise pass silently. The known
    keys come from the driver sources. Returns the number of problems found."""
    known = _known_bench_keys(args)
    if not known:
        print("warning: --input-info: could not read the bench.* keys from benchmarks/bench.cpp; "
              "bench.* typos are not checked", file=sys.stderr)
        return 0
    bad = unknown_bench_keys(inputs, _bench_overrides(args, dim, size, lattice, restart, local), known)
    for key, suggestion in bad:
        hint = f"; did you mean {suggestion}?" if suggestion else ""
        print(f"WARNING: {key} is not read by the benchmark driver (bench.cpp) and has no effect{hint}")
    return len(bad)


def _input_info_checks(args, dim: int, lattice: str, points: list, restart: str, local: bool):
    """Run muphfasa's input_info on the benchmark inputs with each sweep point's overrides.

    The full report is printed for the first point only (a sweep varies size/nodes, not the physics);
    every further size gets its tool warnings, and every (size, nodes) point a memory line. In batch
    mode that line compares the estimate with the GPU memory of the requested nodes. Problems are only
    reported, never fatal. points: [(size, nodes), ...]."""
    bench_rel = CASES[args.case]["bench_dir"]
    if not bench_rel:
        print(f"input_info: '{args.case}' has no inputs.bench (validate-only case); skipping", file=sys.stderr)
        return
    inputs = BENCH_DIR / bench_rel / "inputs.bench"
    issues = _check_bench_keys(args, dim, lattice, points[0][0], restart, local, inputs)

    exe = find_input_info(dim, args.input_info_exe)
    if exe is None:
        looked = args.input_info_exe or os.path.relpath(INPUT_INFO_DIR / f"input_info{dim}d.*.ex", BENCH_DIR)
        print(f"warning: --input-info: input_info not found ({looked}); skipping the input check.\n"
              f"         Build it with:  make -C {os.path.relpath(INPUT_INFO_DIR, BENCH_DIR)} DIM={dim}",
              file=sys.stderr)
        return
    restart_here = bool(restart) and os.path.isfile(os.path.join(restart, "Header"))
    print(f"\n=== input_info: {Path(exe).name} on {os.path.relpath(inputs, BENCH_DIR)} ===")
    if restart and not restart_here:
        print("  note: restart checkpoint not found on this machine; it is not checked here "
              "(the compute node checks it at run time)")

    results = {}
    known_bench = _known_bench_keys(args)
    for size, nodes in points:
        first_of_size = size not in results
        if first_of_size:
            ov = _bench_overrides(args, dim, size, lattice, restart, local)
            bench = dict(o.split("=", 1) for o in ov)
            if restart and not restart_here:
                ov = [o for o in ov if not o.startswith("amr.restart=")]
            if "bench.steps" in bench:
                # input_info knows only the absolute sim.max_step; bench.cpp runs warmup + steps from
                # the start step (a bench checkpoint is named chk<steps done>)
                m = re.search(r"(\d+)/*$", restart) if restart else None
                if not restart or m:
                    start = int(m.group(1)) if m else 0
                    ov.append(f"sim.max_step={start + int(bench['bench.warmup']) + int(bench['bench.steps'])}")
            results[size] = run_input_info(exe, inputs, ov)
            res = results[size]
            if len(results) == 1:
                print(f"  first point: size={size}  overrides: {' '.join(repr(o) if ' ' in o else o for o in ov)}")
                print("  " + "-" * 72)
                print(drop_bench_keys(res["text"], known_bench).rstrip("\n"))
                print("  " + "-" * 72)
            elif res["warnings"]:
                print(f"  size={size}: input_info warnings:")
                for w in res["warnings"]:
                    print(f"      WARNING: {w}")
            if res["rc"] != 0:
                issues += 1
                if len(results) > 1:
                    print(f"  size={size}: input_info failed (exit {res['rc']}):\n{res['text'].rstrip()}")
            issues += len(res["warnings"])
        res = results[size]
        if res["memory_bytes"] is None:
            continue
        if local:
            if first_of_size:
                print(f"  size={size}: est. memory {res['memory_bytes'] / 2**30:.4g} GiB")
        else:
            level, msg = node_memory_check(res["memory_bytes"], nodes)
            print(f"  size={size} nodes={nodes}: {'WARNING: ' if level != 'ok' else ''}{msg}")
            issues += level != "ok"
    print(f"=== input_info: {issues} warning/error(s) above (reported only, nothing was blocked) ===\n"
          if issues else "=== input_info: no problems found ===\n")


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

    # if strong scaling then match len size to nodes
    if args.mode == 'strong':
        sizes = [args.sizes[0]]*len(args.nodes)
    else:
        sizes = args.sizes
    pairs = list(zip(sizes, itertools.islice(itertools.cycle(args.nodes), len(sizes))))
    print(f"Batch mode: {args.case} dim={dim} lattice={lattice} ppn={args.ppn} queue={args.queue or 'auto'}")
    print(f"  sweep pairs (size ↔ nodes): {pairs}")
    if args.max_grid_size != 32:
        print("  note: --max-grid-size is a local-run AMReX box tiling knob; batch scripts use the bench default", file=sys.stderr)

    walltime = args.walltime or batch_generate.DEFAULT_WALLTIME
    queues = [args.queue or batch_generate.auto_queue(nodes) for _, nodes in pairs]
    if not args.sunspot and not args.no_queue_check:
        # before anything is generated: PBS rejects a request outside the queue's limits
        batch_generate.validate_queue_requests(
            [(q, nodes, f"size={size} nodes={nodes}") for q, (size, nodes) in zip(queues, pairs)], walltime)

    if args.input_info:
        restart = os.path.abspath(os.path.expanduser(args.restart)) if args.restart else ""
        if "performance" not in (args.stages or batch_generate.ALL_STAGES):
            restart = ""
        _input_info_checks(args, dim, lattice, pairs, restart, local=False)

    written = []
    for size, nodes in pairs:
        cmd = [
            sys.executable, str(gen),
            "--case", args.case,
            "--nodes", str(nodes),
            "--size", str(size),
            "--ppn", str(args.ppn),
            "--no-submit-hint",   # the sweep-wide loop is printed below
        ]
        if lattice:
            cmd += ["--lattice", lattice]
        if args.queue:
            cmd += ["--queue", args.queue]
        if args.walltime:
            cmd += ["--walltime", args.walltime]
        if args.no_queue_check:
            cmd += ["--no-queue-check"]
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
        if args.sunspot:
            cmd += ["--sunspot"]
        if args.save_checkpoint:
            cmd += ["--save-checkpoint"]
        if args.restart:
            cmd += ["--restart", args.restart]
        if args.outdir:
            cmd += ["--outdir", str(args.outdir)]
        if args.dry_run:
            cmd.append("--dry-run")

        print(f"  $ {' '.join(cmd)}", flush=True)
        result = subprocess.run(cmd)
        if result.returncode != 0:
            print(f"generate failed for size={size} nodes={nodes} (exit {result.returncode})", file=sys.stderr)
            sys.exit(result.returncode)
        written.append((size, nodes))

    stages = sorted(set(args.stages or batch_generate.ALL_STAGES), key=batch_generate.ALL_STAGES.index)
    restart = os.path.abspath(os.path.expanduser(args.restart)) if args.restart else ""
    if "performance" not in stages:
        restart = ""
    tail = batch_generate.name_tail(stages, args.save_checkpoint, restart)
    outdir = Path(args.outdir) if args.outdir else batch_generate.REPO_ROOT / "batch" / "generated"
    outdir_rel = os.path.relpath(outdir, batch_generate.REPO_ROOT)
    records = [(args.case, size, lattice, nodes) for size, nodes in written]
    verb = "would be written" if args.dry_run else "ready"
    print()
    if not args.sunspot:
        for note in batch_generate.queue_notes(queues, len(written)):
            print(f"note: {note}")
    print(f"Batch scripts {verb} for {len(written)} sweep point(s). Submit with (one at a time; "
          f"each qsub blocks until its job finishes):")
    print(batch_generate.submit_hint(records, tail, outdir_rel))
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
    parser.add_argument("--queue", type=str, default=None,
                        help="PBS queue (default: auto from nodes: debug 1-2, debug-scaling 2-256, prod 256+). "
                             "Also checked: nodes and walltime must fit the queue (see --no-queue-check)")
    parser.add_argument("--walltime", type=str, default=None,
                        help="PBS walltime, [[HH:]MM:]SS (default 00:30:00); checked against the queue's limits")
    parser.add_argument("--no-queue-check", action="store_true",
                        help="batch mode: skip the check of nodes / walltime against the Aurora queue limits "
                             "(QUEUE_LIMITS in batch/generate.py)")
    parser.add_argument("--project", type=str, default=None, help="ALCF project (default muphfasa)")
    parser.add_argument("--stages", nargs="+", choices=["compiling", "validating", "performance"], default=None,
                        help="stages to include (default: all three; e.g. --stages performance)")
    parser.add_argument("--warmup", type=int, default=None, help="warmup steps (default per dim: 2D=100, 3D=50)")
    parser.add_argument("--steps", type=int, default=None, help="timed steps (default per dim: 2D=2000, 3D=500)")
    parser.add_argument("--flare-root", type=str, default=None, help="flare run root")
    parser.add_argument("--sunspot", action="store_true", default=False, help="target Sunspot (queue=workq, project=catalyst, filesystems=home:tegu)")
    parser.add_argument("--save-checkpoint", action="store_true", default=False, help="save final checkpoint (chk<max_step>) — batch mode passes bench.save_checkpoint=1; local mode passes bench.save_checkpoint=1 override")
    parser.add_argument("--restart", type=str, default=None, metavar="CHK_DIR",
                        help="restart the timed run from this checkpoint directory (amr.restart). --warmup/--steps "
                             "count from the checkpoint's step (default per dim: 2D=100/2000, 3D=50/500). The "
                             "checkpoint must have the same domain as the run (--sizes); batch mode checks it on "
                             "the compute node, local mode checks it up front.")
    parser.add_argument("--input-info", action="store_true",
                        help="check the run settings with muphfasa's input_info (external/muphfasa/util/input_info, "
                             "must be built there; a warning is printed if not): full report for the first sweep "
                             "point, a memory line per point (batch mode: vs the 768 GiB of GPU memory per node). "
                             "Works with or without --dry-run, in batch and --local mode; warnings are never fatal.")
    parser.add_argument("--input-info-exe", type=str, default=None, metavar="PATH",
                        help="input_info executable to use instead of the one built in external/muphfasa/util/input_info "
                             "(must match the case's dimensionality)")
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

    restart = ""
    if args.restart:
        restart = os.path.abspath(os.path.expanduser(args.restart))
        if not os.path.isfile(os.path.join(restart, "Header")):
            print(f"restart checkpoint not found (no Header file): {restart}", file=sys.stderr)
            sys.exit(1)
        if len(set(args.sizes)) > 1:
            print(f"warning: --restart with several sizes {args.sizes}: a checkpoint has one fixed domain size, "
                  "so runs whose size differs from the checkpoint will abort", file=sys.stderr)

    # if strong scaling then match len size to nodes
    if args.mode == 'strong':
        sizes = [args.sizes[0]]*len(args.nodes)
    else:
        sizes = args.sizes
    points = list(zip(sizes, itertools.islice(itertools.cycle(args.nodes), len(sizes))))

    if args.input_info:
        _input_info_checks(args, dim, lattice, points, restart, local=True)

    if not args.no_build:
        build_bench(args.case, dim, use_mpi=True)

    results = []

    for size, nodes in points:
        nprocs = max(1, nodes) * max(1, args.ppn)
        # grid.ncells, lattice and max_grid_size (fixed and independent of the domain size so
        # multi-rank runs get more than one box to distribute, see --max-grid-size help)
        overrides = _bench_overrides(args, dim, size, lattice, restart, local=True)
        solver_exe = find_bench_exe(bench_dir, dim, use_mpi=(nprocs > 1))
        run_dir = Path(f"scaling_{args.case}_{size}_{nodes}" + ("_restart" if restart else ""))
        # run_bench() wipes run_dir first: never let that delete the checkpoint being restarted from
        # (e.g. restarting a restart, whose checkpoint lives in the *_restart directory)
        n_try = 1
        while restart and Path(restart).is_relative_to(run_dir.resolve()):
            n_try += 1
            run_dir = Path(f"scaling_{args.case}_{size}_{nodes}_restart{n_try}")
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
                       "--inputs", str(run_dir / inputs.name), "--log", str(log_path)]
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
