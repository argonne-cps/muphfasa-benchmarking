"""helper.py — Shared build/run helpers for benchmarks (from couette)."""

import os
import re
import shutil
import subprocess
from pathlib import Path

BENCH_DIR = Path(__file__).resolve().parent
REPO_ROOT = BENCH_DIR
MUPHFASA_DIR = REPO_ROOT / "external" / "muphfasa"
AMREX_HOME = Path(os.environ.get("AMREX_HOME", MUPHFASA_DIR / "external" / "amrex"))
PLOTFILE_TOOLS_DIR = AMREX_HOME / "Tools" / "Plotfile"

BENCH_RESULT_RE = re.compile(r"BENCH_RESULT\s+mlups=([0-9.eE+-]+)\s+ms_per_step=([0-9.eE+-]+)\s+elapsed=([0-9.eE+-]+)")

def find_one(pattern: str, directory: Path) -> Path:
    m = sorted(Path(directory).glob(pattern))
    if not m:
        raise FileNotFoundError(f"no file matching {pattern!r} in {directory}")
    return m[0]

def build_bench(case: str, dim: int, comp: str = "gnu", use_mpi: bool = True, use_sycl: bool = False) -> Path:
    bench_dir = BENCH_DIR / "benchmarks" / case
    flags = [f"DIM={dim}", f"COMP={comp}", f"USE_MPI={'TRUE' if use_mpi else 'FALSE'}"]
    if use_sycl:
        flags += ["USE_SYCL=TRUE"]
    print(f"Building {case} DIM={dim} {' '.join(flags)}...")
    subprocess.run(["make", "-j6"] + flags, cwd=bench_dir, check=True)
    # Find the built exe: bench2d.gnu.ex / bench3d.gnu.ex etc.
    exes = sorted(bench_dir.glob("bench*.ex")) + sorted(bench_dir.glob("bench*.gnu.ex"))
    if not exes:
        exes = sorted(bench_dir.glob("*.ex"))
    if not exes:
        raise FileNotFoundError(f"no bench exe in {bench_dir}")
    return exes[0]

def run_bench(solver_exe: Path, inputs_path: Path, run_dir: Path, overrides: list = None) -> str:
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True)
    log_path = run_dir / "run.log"
    cmd = [str(solver_exe), str(inputs_path)] + (overrides or [])
    print(f"Running {' '.join(cmd)} in {run_dir}")
    with open(log_path, "w") as log_file:
        subprocess.run(cmd, cwd=run_dir, stdout=log_file, stderr=subprocess.STDOUT, check=True)
    return log_path.read_text()

def parse_bench_result(log_text: str) -> dict:
    m = BENCH_RESULT_RE.search(log_text)
    if not m:
        raise ValueError("could not find BENCH_RESULT in log")
    return {"mlups": float(m.group(1)), "ms_per_step": float(m.group(2)), "elapsed": float(m.group(3))}
