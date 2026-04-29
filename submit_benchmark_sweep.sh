#!/bin/bash
# submit_benchmark_sweep.sh
# Submits the full (IC file) x (core count) x (tolerance) benchmark factorial.
#
# Edit the arrays below to select what to sweep, then run:
#   chmod +x submit_benchmark_sweep.sh && ./submit_benchmark_sweep.sh
#
# Dry-run (print qsub commands without submitting):
#   DRY_RUN=1 ./submit_benchmark_sweep.sh

DRY_RUN=${DRY_RUN:-0}

# Absolute path to this script's directory — so qsub finds the PBS script
# regardless of where the user runs submit_benchmark_sweep.sh from.
# SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# echo $SCRIPT_DIR

# ---- IC files (timespan controls steady-state approach) ---------------
IC_FILES=(
    # "/home/581/mr9897/RADIv2_Hinne/setup/IC_HF2_shallow_fact_1yr.jl"
    "/home/581/mr9897/RADIv2_Hinne/setup/IC_HF2_shallow_fact_10yr.jl"
    "/home/581/mr9897/RADIv2_Hinne/setup/IC_HF2_shallow_fact_100yr.jl"
    "/home/581/mr9897/RADIv2_Hinne/setup/IC_HF2_shallow_fact_500yr.jl"
)

# ---- Core/memory pairs (ncpus mem_GB) ---------------------------------
# Format: "ncpus mem_GB"
CORE_CONFIGS=(
    # "4 16"
    # "8 32"
    # "16 64"
    # "32 128"
    "48 192"
)

# ---- Tolerance pairs (abstol reltol) ----------------------------------
# Format: "abstol reltol"  label
TOL_CONFIGS=(
    # "1e-6 1e-4 loose"
    "2e-8 2e-6 medium"
    # "1e-8 1e-6 tight"
)

# ---- Walltime (adjust if needed) -------------------------------------
WALLTIME="24:00:00"

# ---- Submit ----------------------------------------------------------
total=0
for ic in "${IC_FILES[@]}"; do
    ic_tag=$(basename "$ic" .jl | sed 's/IC_//')   # e.g. HF2_shallow_fact_10yr
    for cfg in "${CORE_CONFIGS[@]}"; do
        ncpus=$(echo "$cfg" | awk '{print $1}')
        mem=$(echo   "$cfg" | awk '{print $2}')
        for tol in "${TOL_CONFIGS[@]}"; do
            abstol=$(echo "$tol" | awk '{print $1}')
            reltol=$(echo "$tol" | awk '{print $2}')
            tol_label=$(echo "$tol" | awk '{print $3}')

            tag="${ic_tag}_${ncpus}c_${tol_label}"

            # if [[ "$DRY_RUN" == "1" ]]; then
            #     echo "[DRY_RUN] qsub -N radi_${tag} -l ncpus=${ncpus} -l mem=${mem}GB -l walltime=${WALLTIME} -v NCPUS=${ncpus},MEM=${mem},IC_FILE=${ic},ABSTOL=${abstol},RELTOL=${reltol},RUN_TAG=${tag} benchmarking_sweep_run.sh"
            # else
            #     echo "Submitting: $tag"
            #     qsub \
            #       -N "radi_${tag}" \
            #       -l "ncpus=${ncpus}" \
            #       -l "mem=${mem}GB" \
            #       -l "walltime=${WALLTIME}" \
            #       -v "NCPUS=${ncpus},MEM=${mem},IC_FILE=${ic},ABSTOL=${abstol},RELTOL=${reltol},RUN_TAG=${tag}" \
            #       "benchmarking_sweep_run.sh"
            # fi

            if [[ "$DRY_RUN" == "1" ]]; then
                echo "[DRY_RUN] qsub -N radi_${tag} -l ncpus=${ncpus} -l mem=${mem}GB -l walltime=${WALLTIME} -v NCPUS=${ncpus},MEM=${mem},IC_FILE=${ic},ABSTOL=${abstol},RELTOL=${reltol},RUN_TAG=${tag} /home/581/mr9897/RADIv2_Hinne/benchmarking_sweep_run.sh"
            else
                echo "Submitting: $tag"
                echo "qsub -P jk72 -q normal -l walltime=${WALLTIME} -l ncpus=${ncpus} -l mem=${mem}GB -l storage=scratch/jk72 -N radi_${tag} -v "IC_FILE=${ic},ABSTOL=${abstol},RELTOL=${reltol},RUN_TAG=${tag}" "/home/581/mr9897/RADIv2_Hinne/benchmarking_sweep_run.sh""
                echo "qsub \
                    -P jk72 \
                    -q normal \
                    -l walltime=${WALLTIME} \
                    -l ncpus=${ncpus} \
                    -l mem=${mem}GB \
                    -l storage=scratch/jk72 \
                    -N radi_${tag} \
                    -v "IC_FILE=${ic},ABSTOL=${abstol},RELTOL=${reltol},RUN_TAG=${tag}" \
                    "/home/581/mr9897/RADIv2_Hinne/benchmarking_sweep_run.sh" "
                qsub \
                    -P jk72 \
                    -q normal \
                    -l walltime=${WALLTIME} \
                    -l ncpus=${ncpus} \
                    -l mem=${mem}GB \
                    -l storage=scratch/jk72 \
                    -N radi_${tag} \
                    -v "IC_FILE=${ic},ABSTOL=${abstol},RELTOL=${reltol},RUN_TAG=${tag}" \
                    "/home/581/mr9897/RADIv2_Hinne/benchmarking_sweep_run.sh"
            fi
            ((total++))
        done
    done
done

echo ""
echo "Total jobs: $total"
