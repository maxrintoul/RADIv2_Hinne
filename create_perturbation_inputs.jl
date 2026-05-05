## Specify that this is a perturbation run
PERTURBATION_RUN = true

## Load model output
using MAT

model_output = matread("/Users/maxrintoul/repos/RADIv2_09_03_2026/gadi_output/time_cost_tests/sols_all_20260430_055053_725_job167413397-gadi-pbs_pid374883_7a87cf90-2ab6-4df3-a616-e86a2cd4591e.mat");

## Get final model output for each run
fin_sols = model_output["u"][:, :, end, :];

## Get model inputs for each run 

# Input file name

if isdir("/home/581/mr9897/RADIv2_Hinne/setup/")
    input_file = "/home/581/mr9897/RADIv2_Hinne/setup/IC_HF2_shallow_fact_500yr.jl"
elseif isdir("/Users/maxrintoul/repos/RADIv2_09_03_2026/")
    input_file = "/Users/maxrintoul/repos/RADIv2_09_03_2026/setup/IC_HF2_shallow_fact_500yr.jl"
end

Base.include(Main, joinpath(@__DIR__, "modules", "gsw_rho.jl"))

# Load the input file to get 
include(input_file)






