#!/bin/bash
#PBS -P jk72
#PBS -q normal
#PBS -l walltime=06:00:00
#PBS -l ncpus=1
#PBS -l mem=8GB
#PBS -l storage=scratch/jk72
#PBS -N radi_sweep
#PBS -j oe
# NOTE: ncpus, mem, walltime, and job name are overridden by qsub command-line
# flags in submit_benchmark_sweep.sh. Values above are fallback placeholders only.

# -------------------------------------------------------------------
# Parametric PBS script for RADI sweep benchmarking.
# Submit via submit_benchmark_sweep.sh (or manually with qsub -v).
#
# Required variables (pass with qsub -v):
#   NCPUS     — number of CPU cores  (e.g. 8)
#   MEM       — memory in GB         (e.g. 32)
#   IC_FILE   — path to IC file      (e.g. setup/IC_HF2_shallow_fact_10yr.jl)
#   ABSTOL    — ODE absolute tol     (e.g. 1e-7)
#   RELTOL    — ODE relative tol     (e.g. 1e-5)
#   RUN_TAG   — label for job name   (e.g. HF2_10yr_8c_tol1)
# -------------------------------------------------------------------

module load julia/1.11.0

cd $PBS_O_WORKDIR

export JULIA_NUM_THREADS=${NCPUS}
export OPENBLAS_NUM_THREADS=1

echo "========================================"
echo "RUN_TAG  : ${RUN_TAG}"
echo "IC_FILE  : ${IC_FILE}"
echo "NCPUS    : ${NCPUS}"
echo "ABSTOL   : ${ABSTOL}"
echo "RELTOL   : ${RELTOL}"
echo "Start    : $(date)"
echo "========================================"

time julia --project=. runfile_ensemble_sweep.jl "${IC_FILE}" "${ABSTOL}" "${RELTOL}"

echo "========================================"
echo "End      : $(date)"
echo "========================================"
