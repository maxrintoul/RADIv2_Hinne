#!/bin/bash
#PBS -P jk72
#PBS -q normal
#PBS -l walltime=06:00:00
#PBS -l ncpus=8
#PBS -l mem=32GB
#PBS -l storage=scratch/jk72
#PBS -N radi_run_8cores
#PBS -j oe

# Load Julia module
module load julia/1.11.0

# Change to your working directory
cd $PBS_O_WORKDIR

export JULIA_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=1

julia --project=. runfile_ensemble.jl