using DelimitedFiles
include("../plot_package.jl")
include("../network_package.jl")
include("../analysis_package.jl")

# Ni = 1; Ne = 2;
Ni = 1000; Ne = 4000;
conn_prob = 0.2
seed = 0
cond_version = "1"
code_version = "2.1"
const m3_adj, v_inh, v_exc = read_connection_parameters(Ni, Ne, conn_prob, seed, cond_version);
t_exec = 3000.
t_trans = 50.
ν0 = 2.4

# fileOut = "results/spikeTimes/spike_times_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_conn_prob_" * string(conn_prob) * "_seed_" * string(seed) *"_conductances_2.dat"
# fileOut = "results/spikeTimes/spike_times_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_conn_prob_" * string(conn_prob) * "_seed_" * string(seed) *"_conductances_2_30.dat"
# fileOut = "results/spikeTimes/spike_times_LIF_v2.0_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_conn_prob_" * string(conn_prob) * "_seed_" * string(seed) *"_conductances_2.dat"
fileOut = "results/spikeTimes/spike_times_LIF_v" * string(code_version) * "_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_conn_prob_" * string(conn_prob) * "_seed_" * string(seed) *"_conductances_2_rate_" * string(ν0) * ".dat"
tt = readdlm(fileOut)
N = length(tt[:,1])
tt = [filter!(x->(x≠""), tt[i,:]) for i=1:N]


# plotSpikeRate(m_tS_ext, 1, 0, t_exec, "results/spikeTimes/spikeRateNoise_LIF_v" * code_version * "_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_seed_" * string(seed) * "_conductances_1.png")
# plotRaster(tt, v_inh, "results/RP/RP_LIF_v" * code_version * "_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_seed_" * string(seed) * "_conductances_1.png")
# plotSpikeRate_IE(tt, 1, 0, t_exec, v_inh, "results/spikeTimes/spikeRateNetwork_LIF_v" * code_version * "_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_seed_" * string(seed) * "_conductances_1_rate_" * string(ν0) * "_10.png")
# plotRaster(tt, v_inh, t_trans, t_exec, "results/RP/RP_LIF_v" * code_version * "_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_seed_" * string(seed) * "_conductances_1_rate_" * string(ν0) * "_10.png")

# plotSpikeRate_IE(tt, 1, t_trans, t_exec, v_inh, "results/spikeTimes/spikeRateNetwork_LIF_v" * code_version * "_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_seed_" * string(seed) * "_conductances_1_rate_" * string(ν0) * "_zoom.png")
plotSpikeRate_IE(tt, 1, t_trans, t_exec, v_inh, "results/spikeTimes/spikeRateNetwork_LIF_v" * code_version * "_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_seed_" * string(seed) * "_conductances_1_rate_" * string(ν0) * ".png", 0, 25, 0, 6)
# plotRaster(tt, v_inh, t_trans, t_exec, "results/RP/RP_LIF_v" * code_version * "_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_seed_" * string(seed) * "_conductances_1_rate_" * string(ν0) * "_10.png")


fileIn = "results/lfp/LFP_LIF_v" * string(code_version) * "_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_conn_prob_" * string(conn_prob) * "_seed_" * string(seed) *"_conductances_2_rate_" * string(ν0) * ".dat"
dt = 0.05
data = readdlm(fileIn);
v_t = data[1,:];
lfp = data[2,:];

plot(v_t, lfp)

analyseLFP(lfp, dt)


# #_---firing rate 
# m_sR, m_t = calcSpikeRate(tt,1, t_trans, t_exec)

# v_average_sR = mean(m_sR)


# analyseLFP(v_average_sR, dt)
