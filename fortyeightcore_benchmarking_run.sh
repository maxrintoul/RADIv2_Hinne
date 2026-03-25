#!/bin/bash
#PBS -P jk72
#PBS -q normal
#PBS -l walltime=02:00:00
#PBS -l ncpus=48
#PBS -l mem=192GB
#PBS -l storage=scratch/jk72
#PBS -N radi_run_48cores
#PBS -j oe

# Load Julia module
module load julia/1.11.0

# Change to your working directory
cd $PBS_O_WORKDIR

export JULIA_NUM_THREADS=48
export OPENBLAS_NUM_THREADS=1

julia --project=. 'include("runfile_ensemble.jl")'