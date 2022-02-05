using DelimitedFiles
using Plots
pyplot()
include("../analysis_package.jl")

fileIn = "results/dtbase_paper/lfp_15.1.out"
data = readdlm(fileIn);
dt = 0.05

v_local_a = data[:,2]
v_lg_a = data[:,4]
# v_tot_a = v_local_a + v_lg_a
v_tot_a = v_local_a
v_local_g = data[:,5]

v_lfp = abs.(v_local_g) + abs.(v_tot_a)
# plot(v_lfp) 

# analyseLFP(v_lfp, dt)
v_fr_i = data[:,11]
v_fr_e = data[:,12]

analyseLFP(v_fr_e, dt)