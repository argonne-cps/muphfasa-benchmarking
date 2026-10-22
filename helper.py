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

BENCH_RESULT_RE = re.compile(
    r"BENCH_RESULT\s+mlups=([0-9.eE+-]+)\s+ms_per_step=([0-9.eE+-]+)\s+elapsed=([0-9.eE+-]+)"
    r"(?:\s+mludps=([0-9.eE+-]+))?"
)

# MPI launcher used for multi-rank local/cluster runs; Aurora's own batch
# script (batch/aurora_bench.sh) calls mpiexec directly and does not go
# through run_bench(), so this only matters for run_scaling.py.
MPI_LAUNCHER = os.environ.get("MUPHFASA_BENCH_LAUNCHER", "mpiexec")

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
    return find_bench_exe(bench_dir, dim, use_mpi=use_mpi)


def find_bench_exe(bench_dir: Path, dim: int, use_mpi: bool = True) -> Path:
    """Pick the right bench executable out of bench_dir.

    AMReX's GNUmake appends ".MPI" to the exe name when USE_MPI=TRUE (e.g.
    bench2d.gnu.MPI.ex vs bench2d.gnu.ex). A non-MPI-linked exe launched under
    mpiexec -n N silently runs as N uncoordinated single-rank copies instead
    of a real N-rank job (each prints "MPI ranks=1"), so this must prefer the
    MPI-linked exe whenever MPI is wanted, not just glob-sort and hope.
    """
    mpi_exes = sorted(bench_dir.glob(f"bench{dim}d.*.MPI.ex"))
    plain_exes = sorted(bench_dir.glob(f"bench{dim}d.*.ex"))
    plain_exes = [e for e in plain_exes if e not in mpi_exes]
    if use_mpi:
        if mpi_exes:
            return mpi_exes[0]
        if plain_exes:
            raise FileNotFoundError(
                f"no MPI-linked bench exe in {bench_dir} (found only {plain_exes[0].name}); "
                f"rebuild with USE_MPI=TRUE"
            )
    else:
        if plain_exes:
            return plain_exes[0]
        if mpi_exes:
            return mpi_exes[0]
    exes = sorted(bench_dir.glob("bench*.ex")) or sorted(bench_dir.glob("*.ex"))
    if not exes:
        raise FileNotFoundError(f"no bench exe in {bench_dir}")
    return exes[0]


_launcher_supports_ppn_cache = {}


def _launcher_supports_ppn(launcher: str) -> bool:
    """--ppn is an MPICH/Cray mpiexec flag (used on Aurora); Open MPI's
    mpiexec does not recognize it and exits nonzero. Detect by launcher
    identity so local Open MPI dev runs don't crash, while Aurora's MPICH
    keeps getting --ppn for real per-node placement control."""
    if launcher not in _launcher_supports_ppn_cache:
        try:
            out = subprocess.run([launcher, "--version"], capture_output=True, text=True, timeout=5)
            supports = "Open MPI" not in (out.stdout + out.stderr)
        except Exception:
            supports = False
        _launcher_supports_ppn_cache[launcher] = supports
    return _launcher_supports_ppn_cache[launcher]


def run_bench(solver_exe: Path, inputs_path: Path, run_dir: Path, overrides: list = None,
              nprocs: int = 1, ppn: int = None) -> str:
    """Run bench exe. nprocs>1 launches via MPI_LAUNCHER (mpiexec by default);
    nprocs==1 runs the exe directly (no MPI wrapper needed/assumed)."""
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True)
    # Use absolute paths so cwd=run_dir still finds exe/inputs (fixes scaling runs)
    solver_exe = Path(solver_exe).resolve()
    inputs_path = Path(inputs_path).resolve()
    log_path = run_dir / "run.log"
    if nprocs and nprocs > 1:
        launch = [MPI_LAUNCHER, "-n", str(nprocs)]
        if ppn and _launcher_supports_ppn(MPI_LAUNCHER):
            launch += ["--ppn", str(ppn)]
        cmd = launch + [str(solver_exe), str(inputs_path)] + (overrides or [])
    else:
        cmd = [str(solver_exe), str(inputs_path)] + (overrides or [])
    print(f"Running {' '.join(cmd)} in {run_dir}")
    with open(log_path, "w") as log_file:
        subprocess.run(cmd, cwd=run_dir, stdout=log_file, stderr=subprocess.STDOUT, check=True)
    return log_path.read_text()

def parse_bench_result(log_text: str) -> dict:
    m = BENCH_RESULT_RE.search(log_text)
    if not m:
        raise ValueError("could not find BENCH_RESULT in log")
    result = {"mlups": float(m.group(1)), "ms_per_step": float(m.group(2)), "elapsed": float(m.group(3))}
    if m.group(4) is not None:
        result["mludps"] = float(m.group(4))
    return result


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
