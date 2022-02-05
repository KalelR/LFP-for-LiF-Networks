include("LIF_package.jl")
include("network_package.jl")
include("noise_package.jl")
# include("plot_package.jl")
include("integra_RK4.jl")
include("analysis_package.jl")

#-----------------------------------network params
# Ni = 100; Ne = 400;
Ni = 1000; Ne = 4000;
conn_prob = 0.2;
Ni = 1; Ne = 2;
N = Ni + Ne
conn_prob = 0

#----------------------------------simulation params
const t_exec = 3000.; #ms
const t_trans = 0.; #ms
seed = 0
code_version = "2.1" #2.05 = 20; 2.1 = 30
cond_version = "1"

#-----------model params
v_τ_m = [20., 10.] #EXC = [1]; INH = [2]
const v_reset = 11.0 #mV, reset potential
const vth = 18.0 #mV
const I = 0.;
v_t_refrac = [2.0, 1.0] #vector with the refractory periods for INH [2] and EXC[1] neurons
v_I = collect(range(I, I, length=N))
u0 = collect(range(v_reset+1, vth-1, length=N))

#-----------connection_params
const m3_adj, v_inh, v_exc = read_connection_parameters(Ni, Ne, conn_prob, seed, cond_version);


#--------------------------------------external input
ν0 = 1.2 #spikes/ms
t_end = t_exec + t_trans
dt = 0.05 #ms ----------- careful
num_digits = 1;
seed = 0;
m_tS_ext = generateExternalInput(N, t_end, dt, ν0, num_digits, seed) #-------------------- CARE: MULTIPLIQUEI RUIDO POR 10 PARA O SINAL BATER MELHOR COM O DO PAPER!
# m3_adj[1][3] = [1,2]
# m3_adj[1][2] = [1]
# m_tS_ext[2] = [-1]; m_tS_ext[3] = [-1]
# m_tS_ext2 = deepcopy(m_tS_ext)

#--auxiliary
max_num_spikes = Int64((t_exec-t_trans)/minimum(v_t_refrac))*10 #maximum number of spikes that a neuron may fire
m_tS =  [ [0.0 for i=1:max_num_spikes] for j=1:N]
v_numS = zeros(Int64, N) #counter of number of spikes at current time t
v_numS_ext = zeros(Int64, N) #counter of number of spikes at current time t
v_tS_ext = [m_tS_ext[i][1] for i=1:N] #vector with current external spike times
m_Ia = zeros(Float64, (N, Int64(floor((t_exec)/dt))))
m_Ig = zeros(Float64, (N, Int64(floor((t_exec)/dt))))
# m_Ia = [ [] for i=1:N]
# m_Ig = [ [] for i=1:N]

#---------integrator params
tspan = [0.0, t_trans+t_exec]
p = [v_τ_m, N, v_I, m3_adj, v_inh, v_t_refrac, m_tS, v_numS, m_tS_ext, v_numS_ext, m_Ia, m_Ig]
# p = [v_τ_m, N, v_I, m3_adj, v_inh, v_t_refrac, m_tS, v_numS, m_tS_ext2]
dt = 0.05

#------------------ CALL SOLVER
println("Now I'll start integrating version " * string(code_version) * ". I have Ni = " * string(Ni) *", Ne = " * string(Ne) *". Running for t_exec = " * string(t_exec) * ", with t_trans = " * string(t_trans) * ".")
v_t, m_y =  @time solver(u0, tspan, dt, p, LIF!, v_reset, vth, N, m_tS, v_numS, m_tS_ext, v_numS_ext) #passar a matriz de novo eh por preguiça, nao precisaria..
println("Finished integration.")


m_tS = [filter!(x->(x≠0.0 && x>t_trans), m_tS[i]) for i=1:N]
fileOut = "results/spikeTimes/spike_times_LIF_v" * code_version * "_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_conn_prob_" * string(conn_prob) * "_seed_" * string(seed) *"_conductances_2_rate_" * string(ν0) * ".dat"
saveSpikeTimes(fileOut, m_tS)

t_trans = 0.
# plotRaster(m_tS, v_inh, t_trans, t_exec, "results/RP/RP_LIF_v" * code_version * "_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_seed_" * string(seed) * "_conductances_1_rate_" * string(ν0) * ".png")

# plotSpikeRate(m_tS_ext, 1, t_trans, t_exec, "results/spikeTimes/spikeRateNoise_LIF_v" * code_version * "_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_seed_" * string(seed) * "_conductances_1_rate_" * string(ν0) * ".png")
# plotSpikeRate_IE(m_tS, 1, t_trans, t_exec, v_inh, "results/spikeTimes/spikeRateNetwork_LIF_v" * code_version * "_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_seed_" * string(seed) * "_conductances_1_rate_" * string(ν0) * ".png")

# fileOut = "results/potentials/potentials_LIF_v" * code_version * "_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_seed_" * string(seed) * "_conductances_1.png"
# plotPotentials(sol, [], fileOut)
# plotPotentials(v_t, m_y, [], "")
# plotPotentials(sol)
# plotPotentialNeuron(sol, 1)


# plotRasterMoreActive(tt, v_inh, 200, t_exec, "")

# LFP = calcLFP(m_Ia[:,10000:end], m_Ig[:,10000:end])
LFP = calcLFP(m_Ia, m_Ig)
fileOut = "results/lfp/LFP_LIF_v" * code_version * "_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_conn_prob_" * string(conn_prob) * "_seed_" * string(seed) *"_conductances_2_rate_" * string(ν0) * ".dat"
# saveLFP(fileOut, LFP, v_t)
writedlm(fileOut, [v_t, LFP])

# analyseLFP(LFP, dt)
# plot(v_t, LFP, xlim=(t_trans, t_exec))



# plotRaster(m_tS_ext, [], t_trans, t_exec)
