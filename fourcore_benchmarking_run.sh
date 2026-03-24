#!/bin/bash
#PBS -P jk72
#PBS -q normal
#PBS -l walltime=02:00:00
#PBS -l ncpus=4
#PBS -l mem=16GB
#PBS -l storage=scratch/jk72
#PBS -N radi_run_4cores
#PBS -j oe

# Load Julia module
module load julia/1.11.0

# Change to your working directory
cd $PBS_O_WORKDIR

export JULIA_NUM_THREADS=4
julia -e 'include("runfile_ensemble.jl")'