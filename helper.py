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
# Shared cache: build once per DIM, copy to .cache/amrex_plotfile_tools/fcompare2d.gnu.ex etc.
SHARED_TOOLS_DIR = REPO_ROOT / ".cache" / "amrex_plotfile_tools"

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
    # Per skills/etiquette.md: do not run heavy -j6 builds on login nodes; use -j4 modest, or qsub -I for full builds
    jflag = os.environ.get("BENCH_J", "4")
    print(f"Building {case} DIM={dim} {' '.join(flags)} (BENCH_J={jflag})...")
    subprocess.run(["make", f"-j{jflag}"] + flags, cwd=bench_dir, check=True)
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
    # Use absolute paths so cwd=run_dir still finds exe/inputs (fixes scaling runs)
    solver_exe = Path(solver_exe).resolve()
    inputs_path = Path(inputs_path).resolve()
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


# ── AMReX plotfile tools with shared cache (same logic as muphfasa/tests/regression/couette/helper.py) ──

def _cached_tool_path(tool: str, dim: int) -> Path:
    return SHARED_TOOLS_DIR / f"{tool}{dim}d.gnu.ex"


def _verify_tool_dim(exe: Path, expected_dim: int, tool: str):
    cfg = PLOTFILE_TOOLS_DIR / "tmp_build_dir" / "s" / f"{expected_dim}d.gnu.EXE" / "AMReX_Config.H"
    if cfg.is_file():
        txt = cfg.read_text()
        if f"#define AMREX_SPACEDIM {expected_dim}" not in txt:
            raise RuntimeError(f"{tool} DIM check failed: {cfg} does not contain AMREX_SPACEDIM {expected_dim}")
        return
    try:
        out = subprocess.run(["strings", str(exe)], capture_output=True, text=True, timeout=5).stdout
        if f"AMREX_SPACEDIM {expected_dim}" not in out and f"BL_SPACEDIM {expected_dim}" not in out:
            other = 3 if expected_dim == 2 else 2
            if f"AMREX_SPACEDIM {other}" in out or f"BL_SPACEDIM {other}" in out:
                raise RuntimeError(f"{tool} {exe} appears to be DIM={other}, expected DIM={expected_dim}")
            print(f"Warning: could not verify {tool} DIM={expected_dim} via strings; assuming build is correct")
    except RuntimeError:
        raise
    except Exception:
        pass


def ensure_fcompare(dim: int = 2) -> Path:
    cached = _cached_tool_path("fcompare", dim)
    if cached.is_file():
        try:
            _verify_tool_dim(cached, dim, "fcompare")
            return cached
        except Exception as e:
            print(f"Cached fcompare {cached} failed verification ({e}), rebuilding...")
            try:
                cached.unlink()
            except FileNotFoundError:
                pass
    print(f"Building fcompare (DIM={dim})...")
    subprocess.run(["make", f"DIM={dim}", "programs=fcompare"], cwd=PLOTFILE_TOOLS_DIR, check=True)
    built = find_one("fcompare*.ex", PLOTFILE_TOOLS_DIR)
    _verify_tool_dim(built, dim, "fcompare")
    SHARED_TOOLS_DIR.mkdir(parents=True, exist_ok=True)
    shutil.copy2(built, cached)
    try:
        cached.chmod(0o755)
    except Exception:
        pass
    return cached


def ensure_faverage(dim: int = 2) -> Path:
    cached = _cached_tool_path("faverage", dim)
    if cached.is_file():
        try:
            _verify_tool_dim(cached, dim, "faverage")
            return cached
        except Exception as e:
            print(f"Cached faverage {cached} failed verification ({e}), rebuilding...")
            try:
                cached.unlink()
            except FileNotFoundError:
                pass
    print(f"Building faverage (DIM={dim})...")
    subprocess.run(["make", f"DIM={dim}", "programs=faverage"], cwd=PLOTFILE_TOOLS_DIR, check=True)
    built = find_one("faverage*.ex", PLOTFILE_TOOLS_DIR)
    _verify_tool_dim(built, dim, "faverage")
    SHARED_TOOLS_DIR.mkdir(parents=True, exist_ok=True)
    shutil.copy2(built, cached)
    try:
        cached.chmod(0o755)
    except Exception:
        pass
    return cached
