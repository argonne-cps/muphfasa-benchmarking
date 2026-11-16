"""helper.py — Shared build/run helpers for benchmarks (from couette)."""

import difflib
import math
import os
import re
import shlex
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

# Standalone input-file inspector shipped with muphfasa; built per dimension
# (make -C <INPUT_INFO_DIR> DIM=3 -> input_info3d.gnu.ex).
INPUT_INFO_DIR = MUPHFASA_DIR / "util" / "input_info"

# Aurora / Sunspot node: 6 PVC GPUs x 128 GiB HBM.
GPU_MEM_PER_NODE_GIB = 768
# input_info's memory estimate covers the distributions, rho and u only. MPI buffers,
# the AMReX arena and temporaries need headroom, so warn above this fraction of a node.
GPU_MEM_USABLE_FRACTION = 0.8

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


def apply_overrides(text: str, overrides) -> str:
    """Apply "key=value" overrides to the text of an AMReX inputs file.

    Any existing line defining a key is dropped and the new definitions are
    appended in a labelled block, so the result is a self-contained inputs file
    (no duplicate keys, nothing left for the command line to override).
    """
    if not overrides:
        return text
    new_lines = []
    for ov in overrides:
        key, sep, value = ov.partition("=")
        key, value = key.strip(), value.strip()
        if not sep or not key:
            raise ValueError(f"override must look like key=value, got {ov!r}")
        text = re.sub(rf"^[ \t]*{re.escape(key)}[ \t]*=.*\n?", "", text, flags=re.MULTILINE)
        new_lines.append(f"{key} = {value}")
    if not text.endswith("\n"):
        text += "\n"
    return text + "\n# ── Run settings written by run_scaling.py ──\n" + "\n".join(new_lines) + "\n"


def run_bench(solver_exe: Path, inputs_path: Path, run_dir: Path, overrides: list = None,
              nprocs: int = 1, ppn: int = None) -> str:
    """Run bench exe. nprocs>1 launches via MPI_LAUNCHER (mpiexec by default);
    nprocs==1 runs the exe directly (no MPI wrapper needed/assumed).

    The overrides are written into run_dir's copy of the inputs file (same file
    name), so `<exe> <that file>` reproduces the run with no extra arguments;
    the exact command is also the first line of run.log.
    """
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True)
    # Use absolute paths so cwd=run_dir still finds exe/inputs (fixes scaling runs)
    solver_exe = Path(solver_exe).resolve()
    inputs_path = Path(inputs_path).resolve()
    log_path = run_dir / "run.log"
    run_inputs = run_dir / inputs_path.name
    run_inputs.write_text(apply_overrides(inputs_path.read_text(), overrides))
    if nprocs and nprocs > 1:
        launch = [MPI_LAUNCHER, "-n", str(nprocs)]
        if ppn and _launcher_supports_ppn(MPI_LAUNCHER):
            launch += ["--ppn", str(ppn)]
        cmd = launch + [str(solver_exe), run_inputs.name]
    else:
        cmd = [str(solver_exe), run_inputs.name]
    print(f"Running {shlex.join(cmd)} in {run_dir}")
    with open(log_path, "w") as log_file:
        log_file.write(f"# Command: {shlex.join(cmd)}\n# cwd: {run_dir.resolve()}\n")
        log_file.flush()
        subprocess.run(cmd, cwd=run_dir, stdout=log_file, stderr=subprocess.STDOUT, check=True)
    return log_path.read_text()

# ── input_info (muphfasa/util/input_info) ────────────────────────────────────

_SIZE_UNITS = {"B": 1, "KiB": 2**10, "MiB": 2**20, "GiB": 2**30, "TiB": 2**40}
_MEMORY_RE = re.compile(r"^\s*memory \(est\.\)\s+~([0-9.eE+-]+)\s*(B|KiB|MiB|GiB|TiB)\b", re.MULTILINE)
_WARNING_RE = re.compile(r"^\s*WARNING: (.+)$", re.MULTILINE)


def find_input_info(dim: int, exe=None):
    """Path of the built input_info for this dimensionality, or None if not built.

    exe overrides the search. Otherwise look for input_info<dim>d.*.ex in
    external/muphfasa/util/input_info, preferring the plain serial build
    (shortest name, i.e. no .MPI/.DEBUG suffix)."""
    if exe:
        exe = Path(exe).expanduser()
        return exe if exe.is_file() else None
    found = sorted(INPUT_INFO_DIR.glob(f"input_info{dim}d.*.ex"), key=lambda p: (len(p.name), p.name))
    return found[0] if found else None


def run_input_info(exe, inputs_path, overrides=()) -> dict:
    """Run input_info on inputs_path with key=value overrides (AMReX command-line style).

    Returns {"cmd", "rc", "text", "memory_bytes" (None if not reported), "warnings" (headlines)}."""
    cmd = [str(Path(exe).resolve()), str(inputs_path)] + list(overrides)
    r = subprocess.run(cmd, capture_output=True, text=True)
    text = r.stdout + r.stderr
    m = _MEMORY_RE.search(text)
    return {
        "cmd": cmd,
        "rc": r.returncode,
        "text": text,
        "memory_bytes": float(m.group(1)) * _SIZE_UNITS[m.group(2)] if m else None,
        "warnings": _WARNING_RE.findall(text),
    }


def bench_keys_read(*sources) -> set:
    """Names of the bench.* keys that the given driver sources read.

    Finds each `ParmParse <var>("bench")` in the C++ files and collects the key of every
    `<var>.<method>("key", ...)` call, so the set follows the code (nothing to keep in step).
    Missing files are skipped; an empty result means the sources could not be analysed."""
    keys = set()
    for src in sources:
        src = Path(src)
        if not src.is_file():
            continue
        text = src.read_text()
        for var in re.findall(r'ParmParse\s+(\w+)\s*\(\s*"bench"\s*\)', text):
            keys |= set(re.findall(rf'\b{re.escape(var)}\s*\.\s*\w+\s*\(\s*"([^"]+)"', text))
    return keys


def unknown_bench_keys(inputs_path, overrides, known: set) -> list:
    """bench.* keys set in the inputs file or overrides that the driver never reads.

    Returns [(full_key, suggestion or None)] (suggestion: closest key that is read)."""
    seen = []
    for line in Path(inputs_path).read_text().splitlines():
        m = re.match(r"\s*bench\.([^\s=#]+)\s*=", line)
        if m:
            seen.append(m.group(1))
    for ov in overrides:
        m = re.match(r"\s*bench\.([^\s=]+)\s*=", ov)
        if m:
            seen.append(m.group(1))
    bad = []
    for name in dict.fromkeys(seen):
        if name not in known:
            close = difflib.get_close_matches(name, sorted(known), n=1, cutoff=0.6)
            bad.append((f"bench.{name}", f"bench.{close[0]}" if close else None))
    return bad


def drop_bench_keys(text: str, known: set = None) -> str:
    """Remove bench.* entries from input_info's "keys the solver never reads" note.

    The driver (bench.cpp), not the solver, reads those keys, so the note flags them wrongly.
    Only keys in `known` (see bench_keys_read) are dropped; anything else stays in the note so a
    typo is still listed. The note is dropped if nothing is left in it."""
    known = known or set()

    def is_known(line):
        key = line.strip().split("=")[0].strip()
        return key.startswith("bench.") and key[len("bench."):] in known

    out, block = [], None   # block: index in out of the note's headline while inside it
    for line in text.splitlines():
        if line.lstrip().startswith("NOTE: Input keys the solver never reads"):
            block = len(out)
            out.append(line)
        elif block is not None and line.startswith("      "):
            if not is_known(line):
                out.append(line)
        else:
            if block is not None and block == len(out) - 1:
                out.pop()            # note had only bench.* entries
            block = None
            out.append(line)
    if block is not None and block == len(out) - 1:
        out.pop()
    return "\n".join(out) + ("\n" if text.endswith("\n") else "")


def node_memory_check(memory_bytes: float, nodes: int):
    """Compare an input_info memory estimate with the GPU memory of `nodes` Aurora nodes.

    Returns (level, message); level is "ok", "warn" (above GPU_MEM_USABLE_FRACTION of the
    nodes' memory) or "over" (above all of it). Assumes the domain is spread evenly."""
    cap = GPU_MEM_PER_NODE_GIB * 2**30
    gib = 2**30
    per_node = memory_bytes / nodes
    frac = per_node / cap
    min_nodes = max(1, math.ceil(memory_bytes / (cap * GPU_MEM_USABLE_FRACTION)))
    base = (f"est. {memory_bytes / gib:.4g} GiB total = {per_node / gib:.4g} GiB/node "
            f"= {100 * frac:.0f}% of {GPU_MEM_PER_NODE_GIB} GiB")
    if frac > 1.0:
        return "over", f"{base}: does not fit, needs at least {min_nodes} node(s) (to stay under {GPU_MEM_USABLE_FRACTION:.0%})"
    if frac > GPU_MEM_USABLE_FRACTION:
        return "warn", f"{base}: little headroom for MPI buffers/arena, consider >= {min_nodes} node(s)"
    return "ok", base


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
