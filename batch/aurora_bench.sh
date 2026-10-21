#!/bin/bash -l
#PBS -l select=4:system=aurora
#PBS -l walltime=00:30:00
#PBS -q prod
#PBS -A <ALCF_PROJECT>

# Aurora benchmark batch: runs bench.cpp for tgv2d/lid3d/beltrami3d with MPI+SYCL
# Usage: qsub -l select=4:system=aurora batch/aurora_bench.sh
# Override: NNODES=8 CASE=lid3d SIZE=256 LATTICE=D3Q19 qsub ...

set -e
cd ${PBS_O_WORKDIR:-$(dirname $0)/..}

NNODES=${NNODES:-4}
CASE=${CASE:-lid3d}
SIZE=${SIZE:-256}
LATTICE=${LATTICE:-D3Q19}
PPN=${PPN:-12}  # ranks per node (12 on Aurora, 6 GPUs x 2 tiles)

module use /soft/modulefiles
module load cmake
# Intel OneAPI for SYCL
# module load oneapi/eng-compiler/2024.1

export MPICH_GPU_SUPPORT_ENABLED=1
export FI_CXI_ATS=0

echo "Aurora bench: case=$CASE size=$SIZE lattice=$LATTICE nodes=$NNODES ppn=$PPN"

# Build (once per node, cached)
make -C benchmarks/$CASE -j6 DIM=$( [ "$CASE" = "tgv2d" ] && echo 2 || echo 3 ) USE_MPI=TRUE USE_SYCL=TRUE COMP=intel

# Run
BENCH_EXE=$(ls benchmarks/$CASE/bench*.ex benchmarks/$CASE/*.ex 2>/dev/null | head -1)
INPUTS=benchmarks/$CASE/inputs.bench
RUN_DIR=run_${CASE}_${SIZE}_${NNODES}

# mpiexec with GPU awareness: 12 ranks/node, 1 GPU per rank (6 GPUs * 2 tiles)
mpiexec -n $((NNODES * PPN)) --ppn $PPN --gpus-per-node 6 --gpu-bind closest \
  $BENCH_EXE $INPUTS grid.ncells="$SIZE $SIZE $([ "$CASE" = "tgv2d" ] && echo "" || echo $SIZE)" lbm.lattice_type=$LATTICE amr.plot_int=-1 2>&1 | tee $RUN_DIR.log

# Parse MLUPS
grep BENCH_RESULT $RUN_DIR.log || echo "No BENCH_RESULT found"
