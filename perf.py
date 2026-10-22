#!/usr/bin/env python3
"""perf.py — Performance history DB for muphfasa-benchmarking.

SQLite primary (.cache/perf.db) + git-tracked JSONL (results/history.jsonl) source of truth on main.
See /home/biboyd/.opencode/plan/benchmark-perf-db.md
"""
import argparse
import csv
import datetime
import hashlib
import json
import os
import re
import shutil
import socket
import sqlite3
import subprocess
import sys
from pathlib import Path

BENCH_DIR = Path(__file__).resolve().parent
REPO_ROOT = BENCH_DIR
DEFAULT_CACHE_DB = REPO_ROOT / ".cache" / "perf.db"
DEFAULT_JSONL = REPO_ROOT / "results" / "history.jsonl"
DEFAULT_RESULTS_DB = REPO_ROOT / "results" / "perf.db"

SCHEMA = """
CREATE TABLE IF NOT EXISTS runs(
  id INTEGER PRIMARY KEY,
  timestamp TEXT NOT NULL,
  git_sha TEXT NOT NULL, git_dirty INT,
  test_problem TEXT NOT NULL,
  dim INT, lattice TEXT, Q INT,
  ncells TEXT, total_cells INT,
  architecture TEXT, hostname TEXT, queue TEXT,
  mpi INT, nprocs INT, ppn INT, nnodes INT,
  compiler TEXT, comp_flags TEXT, amrex_version TEXT,
  tau REAL, re REAL, u_lat REAL,
  steps INT, warmup INT, elapsed REAL,
  mlups REAL, ms_per_step REAL, mlups_per_rank REAL,
  inputs_hash TEXT, bench_exe_mtime TEXT,
  ci INT
);
"""

# hostname is included so two different real machines both reporting
# architecture="local" don't collide and silently drop a legitimate second
# run (see /home/biboyd/.claude/plans/this-repository-has-been-precious-lamport.md).
DEDUP_COLS = ["git_sha","test_problem","lattice","ncells","nprocs","architecture","hostname","mpi","compiler"]
# Unique index on those + date(timestamp) created in init_db

LATTICE_Q = {"D2Q9":9,"D3Q15":15,"D3Q19":19,"D3Q27":27}

BENCH_RESULT_RE = re.compile(r"BENCH_RESULT\s+mlups=([0-9.eE+-]+)\s+ms_per_step=([0-9.eE+-]+)\s+elapsed=([0-9.eE+-]+)")

def get_git_info():
    try:
        sha = subprocess.check_output(["git","rev-parse","HEAD"], cwd=REPO_ROOT, text=True).strip()
    except Exception:
        sha = "unknown"
    try:
        dirty = 1 if subprocess.check_output(["git","status","--porcelain"], cwd=REPO_ROOT, text=True).strip() else 0
    except Exception:
        dirty = 0
    try:
        amrex_sha = subprocess.check_output(["git","rev-parse","HEAD"], cwd=REPO_ROOT/"external"/"amrex", text=True).strip()[:12]
    except Exception:
        amrex_sha = ""
    return sha, dirty, amrex_sha

def get_architecture():
    hostname = socket.gethostname()
    arch = "local"
    if "aurora" in hostname.lower() or os.getenv("PBS_O_HOST"):
        # also check queue env
        arch = "aurora"
    elif "polaris" in hostname.lower():
        arch = "polaris"
    return arch, hostname

def is_ci():
    """True when running under a recognized CI system (currently only GitHub
    Actions, since that's what benchmark-perf-db.md's v1 targets). There is
    no CI workflow wired up yet (see plan doc "Future Work"), so today this
    is always False in practice and every write to results/history.jsonl on
    `main` comes from a local run -- the `ci` column on each row makes that
    provenance explicit/queryable instead of silently assumed."""
    return os.getenv("GITHUB_ACTIONS") == "true"

def is_main_branch():
    if os.getenv("MUPHFASA_BENCH_TRACK_MAIN") == "1":
        return True
    if os.getenv("MUPHFASA_BENCH_TRACK_MAIN") == "0":
        return False
    if os.getenv("GITHUB_REF") == "refs/heads/main":
        return True
    try:
        branch = subprocess.check_output(["git","rev-parse","--abbrev-ref","HEAD"], cwd=REPO_ROOT, text=True).strip()
        if branch == "main":
            return True
    except Exception:
        pass
    # fallback: check if HEAD == refs/heads/main
    try:
        head = (REPO_ROOT/".git"/"refs"/"heads"/"main").read_text().strip()
        sha,_ ,_ = get_git_info()
        if head == sha:
            return True
    except Exception:
        pass
    return False

def init_db(db_path: Path):
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(db_path)
    conn.executescript(SCHEMA)
    conn.execute("CREATE INDEX IF NOT EXISTS idx_query ON runs(test_problem, architecture, mpi, nprocs, lattice);")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_dedup ON runs(git_sha, test_problem, lattice, ncells, nprocs, architecture, hostname, mpi, compiler, date(timestamp));")
    conn.commit()
    return conn

def parse_inputs_meta(inputs_path: Path, overrides=None):
    text = inputs_path.read_text() if inputs_path and inputs_path.exists() else ""
    # apply overrides to text for parsing
    if overrides:
        for ov in overrides:
            # ov like grid.ncells=64 64 64
            k,v = ov.split("=",1) if "=" in ov else (ov,"")
            # simple replace: if key exists in text, replace line else append
            pat = re.compile(rf"^\s*{re.escape(k)}\s*=.*$", re.M)
            if pat.search(text):
                text = pat.sub(f"{k} = {v}", text)
            else:
                text += f"\n{k} = {v}\n"
    m = re.search(r"^\s*grid\.ncells\s*=\s*([^\n#]+)", text, re.M)
    ncells = m.group(1).strip() if m else ""
    try:
        total = 1
        for n in ncells.split():
            total *= int(n)
    except Exception:
        total = 0
    dim = len(ncells.split()) if ncells else 0
    m = re.search(r"^\s*lbm\.lattice_type\s*=\s*(\S+)", text, re.M)
    lattice = m.group(1).strip() if m else ""
    q = LATTICE_Q.get(lattice,0)
    m = re.search(r"^\s*lbm\.tau\s*=\s*([0-9.eE+-]+)", text, re.M)
    tau = float(m.group(1)) if m else None
    m = re.search(r"^\s*units\.Re\s*=\s*([0-9.eE+-]+)", text, re.M)
    re_val = float(m.group(1)) if m else None
    inputs_hash = hashlib.sha256(text.encode()).hexdigest()[:12] if text else ""
    return {"ncells":ncells,"total_cells":total,"dim":dim,"lattice":lattice,"Q":q,"tau":tau,"re":re_val,"inputs_hash":inputs_hash}

def dedup_exists(conn, rec):
    cur = conn.execute(
        "SELECT id FROM runs WHERE git_sha=? AND test_problem=? AND lattice=? AND ncells=? AND nprocs=? "
        "AND architecture=? AND hostname=? AND mpi=? AND compiler=? AND date(timestamp)=date('now')",
        (rec["git_sha"], rec["test_problem"], rec["lattice"], rec["ncells"], rec["nprocs"],
         rec["architecture"], rec["hostname"], rec["mpi"], rec["compiler"])
    )
    return cur.fetchone() is not None

def log_run(args):
    # Parse BENCH_RESULT from log file or direct values
    log_text = ""
    if args.log:
        log_text = Path(args.log).read_text()
    elif args.mlups is None:
        print("Need --log or --mlups", file=sys.stderr)
        sys.exit(2)
    if args.mlups is not None:
        mlups = args.mlups; ms = args.ms_per_step; elapsed = args.elapsed
    else:
        m = BENCH_RESULT_RE.search(log_text)
        if not m:
            print(f"Could not find BENCH_RESULT in {args.log}", file=sys.stderr)
            sys.exit(2)
        mlups = float(m.group(1)); ms = float(m.group(2)); elapsed = float(m.group(3))

    git_sha, git_dirty, amrex_version = get_git_info()
    if args.git_sha:
        git_sha = args.git_sha
    arch, hostname = get_architecture()
    if args.arch:
        arch = args.arch
    # inputs meta
    inputs_path = Path(args.inputs) if args.inputs else None
    meta = parse_inputs_meta(inputs_path, args.override)
    ncells = args.ncells or meta["ncells"]
    lattice = args.lattice or meta["lattice"]
    dim = meta["dim"] if not args.test else {"tgv2d":2,"lid3d":3,"beltrami3d":3}.get(args.test, meta["dim"])
    total_cells = meta["total_cells"]
    if ncells:
        try:
            total_cells = 1
            for n in ncells.split():
                total_cells *= int(n)
            dim = len(ncells.split())
        except: pass
    q = LATTICE_Q.get(lattice, meta["Q"])
    nprocs = args.nprocs or int(os.getenv("PBS_NNODES","0") or 0)*12 or 1  # fallback
    # PBS/queue
    queue = os.getenv("PBS_QUEUE") or os.getenv("SLURM_PARTITION") or ""
    # compiler
    compiler = args.compiler or os.getenv("COMP","gnu")
    timestamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
    rec = {
        "timestamp": timestamp,
        "git_sha": git_sha, "git_dirty": git_dirty,
        "test_problem": args.test,
        "dim": dim, "lattice": lattice, "Q": q,
        "ncells": ncells, "total_cells": total_cells,
        "architecture": arch, "hostname": hostname, "queue": queue,
        "mpi": int(args.mpi) if args.mpi is not None else (1 if nprocs>1 else 0),
        "nprocs": nprocs, "ppn": args.ppn or 0, "nnodes": args.nnodes or 0,
        "compiler": compiler, "comp_flags": args.comp_flags or "",
        "amrex_version": amrex_version,
        "tau": meta["tau"], "re": meta["re"], "u_lat": None,
        "steps": args.steps or 0, "warmup": args.warmup or 0, "elapsed": elapsed,
        "mlups": mlups, "ms_per_step": ms, "mlups_per_rank": mlups/max(1,nprocs),
        "inputs_hash": meta["inputs_hash"], "bench_exe_mtime": "",
        "ci": 1 if is_ci() else 0,
    }
    # DB handling
    cache_db = Path(args.db) if args.db else DEFAULT_CACHE_DB
    jsonl = Path(args.jsonl) if args.jsonl else DEFAULT_JSONL
    conn = init_db(cache_db)
    if not args.force and dedup_exists(conn, rec):
        print(f"WARN duplicate for {rec['test_problem']} {rec['lattice']} {rec['ncells']} nprocs={rec['nprocs']} sha={rec['git_sha'][:7]} today; use --force to override. Skipping.")
        return 0
    cols = [k for k in rec.keys()]
    vals = [rec[k] for k in cols]
    placeholders = ",".join(["?"]*len(cols))
    conn.execute(f"INSERT INTO runs ({','.join(cols)}) VALUES ({placeholders})", vals)
    conn.commit()
    print(f"Logged to {cache_db}: {rec['test_problem']} {rec['lattice']} {rec['ncells']} mlups={mlups:.1f} sha={git_sha[:7]}")
    if is_main_branch():
        jsonl.parent.mkdir(parents=True, exist_ok=True)
        with open(jsonl,"a") as f:
            f.write(json.dumps(rec, sort_keys=True)+"\n")
        origin = "CI" if rec["ci"] else "local run on main -- NOT a vetted/curated measurement"
        print(f"Appended to {jsonl} ({origin})")
        # also ensure results/perf.db is synced
        sync_db(jsonl, DEFAULT_RESULTS_DB)
    else:
        print(f"Not on main branch (use MUPHFASA_BENCH_TRACK_MAIN=1 or --force to write history); only .cache updated.")
    return 0

def sync_db(jsonl_path: Path, db_path: Path):
    jsonl_path = Path(jsonl_path); db_path = Path(db_path)
    if not jsonl_path.exists():
        print(f"No {jsonl_path} to sync")
        return
    db_path.parent.mkdir(parents=True, exist_ok=True)
    if db_path.exists():
        db_path.unlink()
    conn = init_db(db_path)
    count=0
    with open(jsonl_path) as f:
        for line in f:
            line=line.strip()
            if not line: continue
            rec=json.loads(line)
            # remove id if present to let autoincrement
            rec.pop("id", None)
            cols=list(rec.keys()); vals=[rec[k] for k in cols]
            try:
                conn.execute(f"INSERT INTO runs ({','.join(cols)}) VALUES ({','.join(['?']*len(cols))})", vals)
                count+=1
            except sqlite3.IntegrityError:
                # dedup violation in jsonl (should not happen) skip
                pass
    conn.commit()
    print(f"Synced {count} rows from {jsonl_path} to {db_path}")

def export_db(args):
    db_path = Path(args.db) if args.db else DEFAULT_CACHE_DB
    jsonl_path = Path(args.jsonl) if args.jsonl else DEFAULT_JSONL
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    rows = conn.execute("SELECT * FROM runs ORDER BY timestamp").fetchall()
    jsonl_path.parent.mkdir(parents=True, exist_ok=True)
    with open(jsonl_path,"w") as f:
        for r in rows:
            d=dict(r)
            f.write(json.dumps(d, sort_keys=True)+"\n")
    print(f"Exported {len(rows)} rows to {jsonl_path}")

def query_db(args):
    db_path = Path(args.db) if args.db else (DEFAULT_CACHE_DB if not Path(DEFAULT_RESULTS_DB).exists() else DEFAULT_RESULTS_DB)
    if not db_path.exists():
        # fallback to cache
        db_path = DEFAULT_CACHE_DB
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    where=[]; params=[]
    if args.test:
        where.append("test_problem=?"); params.append(args.test)
    if args.arch:
        where.append("architecture=?"); params.append(args.arch)
    if args.lattice:
        where.append("lattice=?"); params.append(args.lattice)
    if args.mpi is not None:
        where.append("mpi=?"); params.append(1 if args.mpi else 0)
    if args.nprocs:
        where.append("nprocs=?"); params.append(args.nprocs)
    if args.since:
        where.append("timestamp >= ?"); params.append(args.since)
    sql="SELECT * FROM runs"
    if where:
        sql+=" WHERE "+" AND ".join(where)
    sql+=" ORDER BY timestamp"
    rows=conn.execute(sql, params).fetchall()
    fmt=args.format or "table"
    if fmt=="json":
        for r in rows:
            print(json.dumps(dict(r)))
    elif fmt=="csv":
        if rows:
            w=csv.DictWriter(sys.stdout, fieldnames=rows[0].keys())
            w.writeheader()
            for r in rows: w.writerow(dict(r))
    else:
        if not rows:
            print("No rows")
            return
        # simple table
        keys=["timestamp","git_sha","test_problem","lattice","ncells","nprocs","architecture","mlups","ms_per_step"]
        print("\t".join(keys))
        for r in rows:
            print("\t".join(str(r[k])[:19] if k=="timestamp" else (r[k][:7] if k=="git_sha" and r[k] else "") if k=="git_sha" else str(r[k]) for k in keys))
        print(f"\n{len(rows)} rows from {db_path}")

def main():
    ap=argparse.ArgumentParser(description="perf history DB")
    sub=ap.add_subparsers(dest="cmd", required=True)
    # log
    a=sub.add_parser("log", help="log a bench result")
    a.add_argument("--test", required=True, choices=["tgv2d","lid3d","beltrami3d"])
    a.add_argument("--log", type=Path, help="run.log with BENCH_RESULT")
    a.add_argument("--mlups", type=float); a.add_argument("--ms-per-step", dest="ms_per_step", type=float); a.add_argument("--elapsed", type=float)
    a.add_argument("--lattice", type=str); a.add_argument("--ncells", type=str); a.add_argument("--nprocs", type=int)
    a.add_argument("--ppn", type=int); a.add_argument("--nnodes", type=int)
    a.add_argument("--arch", type=str); a.add_argument("--mpi", action=argparse.BooleanOptionalAction, default=None)
    a.add_argument("--inputs", type=Path); a.add_argument("--override", action="append", default=[])
    a.add_argument("--compiler", type=str); a.add_argument("--comp-flags", type=str)
    a.add_argument("--steps", type=int); a.add_argument("--warmup", type=int)
    a.add_argument("--git-sha", type=str)
    a.add_argument("--db", type=Path); a.add_argument("--jsonl", type=Path)
    a.add_argument("--force", action="store_true", help="bypass dedup")
    a.set_defaults(func=lambda a: log_run(a))
    # sync
    b=sub.add_parser("sync", help="rebuild results/perf.db from history.jsonl")
    b.add_argument("--jsonl", type=Path, default=DEFAULT_JSONL); b.add_argument("--db", type=Path, default=DEFAULT_RESULTS_DB)
    b.set_defaults(func=lambda a: (sync_db(a.jsonl, a.db),0)[1])
    # export
    c=sub.add_parser("export", help="export DB to JSONL")
    c.add_argument("--db", type=Path, default=DEFAULT_CACHE_DB); c.add_argument("--jsonl", type=Path, default=DEFAULT_JSONL)
    c.set_defaults(func=lambda a: (export_db(a),0)[1])
    # query
    d=sub.add_parser("query", help="query runs")
    d.add_argument("--test", required=False); d.add_argument("--arch", required=False); d.add_argument("--lattice", required=False)
    d.add_argument("--mpi", action=argparse.BooleanOptionalAction, default=None); d.add_argument("--nprocs", type=int)
    d.add_argument("--since", type=str); d.add_argument("--format", choices=["table","csv","json"], default="table")
    d.add_argument("--db", type=Path)
    d.set_defaults(func=lambda a: (query_db(a),0)[1])
    # plot
    e=sub.add_parser("plot", help="plot mlups vs time")
    e.add_argument("--test", required=True); e.add_argument("--metric", default="mlups", choices=["mlups","ms_per_step"])
    e.add_argument("--group-by", default="nprocs"); e.add_argument("--out", type=Path, required=True)
    e.add_argument("--db", type=Path)
    def do_plot(a):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        db_path = Path(a.db) if a.db else (DEFAULT_CACHE_DB if not DEFAULT_RESULTS_DB.exists() else DEFAULT_RESULTS_DB)
        conn=sqlite3.connect(db_path); conn.row_factory=sqlite3.Row
        rows=conn.execute(f"SELECT * FROM runs WHERE test_problem=? ORDER BY timestamp", (a.test,)).fetchall()
        if not rows:
            print("No rows for plot"); return 1
        # group
        groups={}
        for r in rows:
            g=r[a.group_by] if a.group_by in r.keys() else "all"
            groups.setdefault(g, []).append(r)
        fig,ax=plt.subplots()
        for g, rs in groups.items():
            xs=[r["timestamp"][:10] for r in rs]; ys=[r[a.metric] for r in rs]
            ax.plot(xs, ys, "o-", label=str(g))
        ax.set_xlabel("date"); ax.set_ylabel(a.metric); ax.set_title(f"{a.test} {a.metric} by {a.group_by}")
        ax.legend(); ax.tick_params(axis='x', rotation=30)
        fig.tight_layout(); fig.savefig(a.out, dpi=150); print(f"Wrote {a.out}")
        return 0
    e.set_defaults(func=do_plot)

    args=ap.parse_args()
    sys.exit(args.func(args))

if __name__=="__main__":
    main()
