#!/bin/bash -l
#PBS -l select=4:system=aurora
#PBS -l walltime=00:30:00
#PBS -l filesystems=flare
#PBS -q prod
#PBS -A <ALCF_PROJECT>
# See skills/aurora-overview.md and skills/etiquette.md: use flare (or daos_user),
# do not submit from /soft/modulefiles, queue limits: prod routes small/medium/large,
# debug=1-2 nodes max 1 job/user, capacity=1-16 nodes. Check qstat before submit.

# Aurora benchmark batch: runs bench.cpp for tgv2d/tgv3d/lid3d/beltrami3d with MPI+SYCL
# Usage:
#   qsub -l select=4:system=aurora -l filesystems=flare -A <proj> -q prod \
#        batch/aurora_bench.sh -- --nodes=4 --case=lid3d --size=256 --lattice=D3Q19
#   (args after `--` are forwarded to the script; `NNODES=4 CASE=lid3d qsub ...`
#    env-var form still works, CLI args take precedence)

set -e
cd ${PBS_O_WORKDIR:-$(dirname $0)/..}

# Defaults (env-vars kept for backwards compat with `NNODES=4 qsub ...`)
NNODES=${NNODES:-4}
CASE=${CASE:-lid3d}
SIZE=${SIZE:-256}
LATTICE=${LATTICE:-D3Q19}
PPN=${PPN:-12}  # ranks per node (12 on Aurora, 6 GPUs x 2 tiles)

# CLI args: --nodes=, --nnodes=, --case=, --size=, --lattice=, --ppn=
# Also accepts --nodes 4 / -n 4 style for convenience. Unknown args warn.
while [[ $# -gt 0 ]]; do
  case "$1" in
    --nodes=*|--nnodes=*) NNODES="${1#*=}"; shift ;;
    --nodes|--nnodes|-n)   NNODES="$2"; shift 2 ;;
    --case=*)             CASE="${1#*=}"; shift ;;
    --case|-c)            CASE="$2"; shift 2 ;;
    --size=*)             SIZE="${1#*=}"; shift ;;
    --size|-s)            SIZE="$2"; shift 2 ;;
    --lattice=*)          LATTICE="${1#*=}"; shift ;;
    --lattice)            LATTICE="$2"; shift 2 ;;
    --ppn=*)              PPN="${1#*=}"; shift ;;
    --ppn)                PPN="$2"; shift 2 ;;
    --help|-h)
      echo "Usage: qsub {qsub args} batch/aurora_bench.sh -- --nodes=N --case={tgv2d|tgv3d|lid3d|beltrami3d} --size=N [--lattice=D3Q19] [--ppn=12]"
      echo "  Env-var form still works: NNODES=4 CASE=lid3d SIZE=256 LATTICE=D3Q19 PPN=12 qsub ..."
      exit 0 ;;
    --) shift; ;;  # qsub delimiter: `qsub ... script.sh -- --nodes=4`
    *) echo "WARNING: unknown arg $1" >&2; shift ;;
  esac
done

# Aurora modules: prefer /opt/aurora APE, avoid NFS /soft at scale (>1000 nodes)
# per skills/aurora-overview.md: /soft causes scaling issues, do not submit from /soft/modulefiles
if command -v module &>/dev/null; then
  module use /soft/modulefiles || true
  module load cmake || true
  # Intel oneAPI DPC++ for SYCL on PVC (e.g. oneapi/eng-compiler/2024.2.1 or aurora-PE)
  # module load oneapi
fi

export MPICH_GPU_SUPPORT_ENABLED=1
export FI_CXI_ATS=0
# Proxy required on compute nodes (see aurora-overview.md)
export HTTP_PROXY=http://proxy.alcf.anl.gov:3128
export HTTPS_PROXY=http://proxy.alcf.anl.gov:3128
export http_proxy=http://proxy.alcf.anl.gov:3128
export https_proxy=http://proxy.alcf.anl.gov:3128
# Aurora: cores 0 and 52 reserved for OS since Mar 2025; do not bind to them
# mpiexec will handle via --cpu-bind; avoid explicit binding to 0,52

echo "Aurora bench: case=$CASE size=$SIZE lattice=$LATTICE nodes=$NNODES ppn=$PPN"

# Build: do NOT compile on login nodes for large jobs; pre-build or build once per allocation.
# Per skills/etiquette.md: login nodes are shared, heavy compiles belong on compute (qsub -I) or via batch.
# If binary missing, build with modest -j (not -j6 on login) to avoid disrupting others.
if [ ! -f "benchmarks/$CASE/bench3d.gnu.ex" ] && [ ! -f "benchmarks/$CASE/bench2d.gnu.ex" ]; then
  echo "WARNING: bench exe not found — building on compute node (use login only for -j2 small tests)"
  make -C benchmarks/$CASE -j4 DIM=$( [ "$CASE" = "tgv2d" ] && echo 2 || echo 3 ) USE_MPI=TRUE USE_SYCL=TRUE COMP=intel
fi
# Checkpointing: Aurora node failures common (aurora-overview.md); bench is short (110 steps) so no chk needed,
# but long production runs should set amr.chk_int and write to /lus/flare (not DAOS) every 15-60 min.

# Run
BENCH_EXE=$(ls benchmarks/$CASE/bench*.ex benchmarks/$CASE/*.ex 2>/dev/null | head -1)
INPUTS=benchmarks/$CASE/inputs.bench
RUN_DIR=run_${CASE}_${SIZE}_${NNODES}

# mpiexec with GPU awareness: 12 ranks/node, 1 GPU per rank (6 GPUs * 2 tiles)
mpiexec -n $((NNODES * PPN)) --ppn $PPN --gpus-per-node 6 --gpu-bind closest \
  $BENCH_EXE $INPUTS grid.ncells="$SIZE $SIZE $([ "$CASE" = "tgv2d" ] && echo "" || echo $SIZE)" lbm.lattice_type=$LATTICE amr.plot_int=-1 2>&1 | tee $RUN_DIR.log

# Parse MLUPS
grep BENCH_RESULT $RUN_DIR.log || echo "No BENCH_RESULT found"

# Log to perf DB (only appends to results/history.jsonl on main branch; use --force to override dedup)
python3 "$(dirname "$0")/../perf.py" log --test "$CASE" --lattice "$LATTICE" --ncells "$SIZE $SIZE $([ "$CASE" = "tgv2d" ] && echo "" || echo $SIZE)" --nprocs $((NNODES * PPN)) --inputs "$INPUTS" --log "$RUN_DIR.log" || echo "perf log failed (non-main branch skips history.jsonl by design)"
