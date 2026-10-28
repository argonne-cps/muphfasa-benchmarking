#!/bin/bash -l
#PBS -l select=4
#PBS -l walltime=00:30:00
#PBS -l filesystems=flare
#PBS -q debug-scaling
#PBS -A muphfasa
# See skills/aurora-overview.md and skills/etiquette.md: use flare (or daos_user),
# do not submit from /soft/modulefiles, queue limits: prod routes small/medium/large,
# debug=1-2 nodes max 1 job/user, capacity=1-16 nodes. Check qstat before submit.

# Aurora benchmark batch: runs bench.cpp for tgv2d/tgv3d/lid2d/lid3d/beltrami3d with MPI+SYCL
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
export ZE_FLAT_DEVICE_HIERARCHY=COMPOSITE
export OMP_NUM_THREADS=1
# Proxy required on compute nodes (see aurora-overview.md)
export HTTP_PROXY=http://proxy.alcf.anl.gov:3128
export HTTPS_PROXY=http://proxy.alcf.anl.gov:3128
export http_proxy=http://proxy.alcf.anl.gov:3128
export https_proxy=http://proxy.alcf.anl.gov:3128
# Aurora: cores 0 and 52 reserved for OS since Mar 2025; do not bind to them
# mpiexec will handle via --cpu-bind; avoid explicit binding to 0,52

# Spatial dimension and grid.ncells string per case (2-D cases take "N N", 3-D "N N N")
case "$CASE" in
  tgv2d|lid2d) DIM=2; NCELLS="$SIZE $SIZE" ;;
  *)           DIM=3; NCELLS="$SIZE $SIZE $SIZE" ;;
esac

echo "Aurora bench: case=$CASE size=$SIZE lattice=$LATTICE nodes=$NNODES ppn=$PPN"

# Build: do NOT compile on login nodes for large jobs; pre-build or build once per allocation.
# Per skills/etiquette.md: login nodes are shared, heavy compiles belong on compute (qsub -I) or via batch.
# If binary missing, build with modest -j (not -j6 on login) to avoid disrupting others.
if [ ! -f "benchmarks/$CASE/bench3d.gnu.ex" ] && [ ! -f "benchmarks/$CASE/bench2d.gnu.ex" ]; then
  echo "WARNING: bench exe not found — building on compute node (use login only for -j2 small tests)"
  make -C benchmarks/$CASE -j4 DIM=$DIM USE_MPI=TRUE USE_SYCL=TRUE COMP=intel
fi
# Checkpointing: Aurora node failures common (aurora-overview.md); bench is short (110 steps) so no chk needed,
# but long production runs should set amr.chk_int and write to /lus/flare (not DAOS) every 15-60 min.

# Run
BENCH_EXE=$(ls benchmarks/$CASE/bench*.ex benchmarks/$CASE/*.ex 2>/dev/null | head -1)
INPUTS=benchmarks/$CASE/inputs.bench
RUN_DIR=run_${CASE}_${SIZE}_${NNODES}

# GPU binding per https://docs.alcf.anl.gov/aurora/running-jobs-aurora/#binding-mpi-ranks-to-gpus
# 6 PVC GPUs x 2 tiles = 12 tiles/node. PPN=12 uses gpu_tile_compact.sh (1 rank/tile),
# PPN=6 uses gpu_dev_compact.sh (1 rank/GPU). Wrapper sets ZE_AFFINITY_MASK.
if [ "$PPN" -eq 12 ]; then
  CPU_BIND="--cpu-bind=list:1-8:9-16:17-24:25-32:33-40:41-48:53-60:61-68:69-76:77-84:85-92:93-100"
  GPU_WRAPPER="gpu_tile_compact.sh"
elif [ "$PPN" -eq 6 ]; then
  CPU_BIND="--cpu-bind=list:1-16:17-32:33-48:53-68:69-84:85-100"
  GPU_WRAPPER="gpu_dev_compact.sh"
else
  CPU_BIND="--cpu-bind=depth --depth=8"
  GPU_WRAPPER="gpu_tile_compact.sh"
fi
mpiexec -n $((NNODES * PPN)) --ppn $PPN $CPU_BIND $GPU_WRAPPER \
  $BENCH_EXE $INPUTS grid.ncells="$NCELLS" lbm.lattice_type=$LATTICE amr.plot_int=-1 2>&1 | tee $RUN_DIR.log

# Parse MLUPS
grep BENCH_RESULT $RUN_DIR.log || echo "No BENCH_RESULT found"

# Log to perf DB (only appends to results/history.jsonl on main branch; use --force to override dedup)
python3 "$(dirname "$0")/../perf.py" log --test "$CASE" --lattice "$LATTICE" --ncells "$NCELLS" --nprocs $((NNODES * PPN)) --inputs "$INPUTS" --log "$RUN_DIR.log" || echo "perf log failed (non-main branch skips history.jsonl by design)"
