using DifferentialEquations
# using Distributions
# using Random
include("LIF_package.jl")
include("network_package.jl")
include("noise_package.jl")
include("plot_package.jl")


#------network params
# Ni = 100; Ne = 400;
# conn_prob = 0.2;
Ni = 1; Ne = 2;
conn_prob = 0
N = Ni + Ne


#----------simulation params
const t_exec = 45.0; #ms
const t_trans = 0.; #ms
seed = 0
code_version = "0.2"
cond_version = "1"

#--model params
v_τ_m = [20., 10.] #EXC = [1]; INH = [2]
const vr = 11.0 #mV, reset potential
const vth = 18.0 #mV
const I = 0.;
v_t_refrac = [2.0, 1.0] #vector with the refractory periods for INH [2] and EXC[1] neurons
v_I = collect(range(I, I, length=N))
u0 = collect(range(vr+1, vth-1, length=N))

#---connection_params
const m3_adj, v_inh, v_exc = read_connection_parameters(Ni, Ne, conn_prob, seed, cond_version);

#--- external input
ν0 = 1.2 #spikes/ms
t_end = t_exec + t_trans
dt = 0.1 #ms ----------- careful
num_digits = 1;
seed = 0;
# m_tS_ext = generateExternalInput(N, t_end, dt, ν0, num_digits, seed)
m3_adj[1][2] = [1]
# m3_adj[2][3] = [1]
m3_adj[1][1] = [2]
m3_adj[1][3] = [1,2]
m_tS_ext = [[1.0,2.,3.,4.,5., 10.,11.,12.,13.,14.,30.,31.,32.,33.,34.], [20.,21.,22.,23.,24.,35.,36.,37.,38.,39.], [-10.]]
m_tS_ext2 = deepcopy(m_tS_ext)
# println(m_tS_ext)
# v_hasSpiked = [false for i=1:N]


#--auxiliary
max_num_spikes = Int64((t_exec-t_trans)/minimum(v_t_refrac))*10 #maximum number of spikes that a neuron may fire#could be more memory conversative by using refrac period of each neuron
tt =  [ [0.0 for i=1:max_num_spikes] for j=1:N]

v_numSpikes = zeros(Int64, N) #counter of number of spikes at current time t
# v_numSpikes_ext = ones(Int64, N) #counter of number of spikes at current time t
v_tS_ext = [m_tS_ext[i][1] for i=1:N] #vector with current external spike times


# m_I = [[ [] for j = 1:N ] for i=1:2]
# v_t = []


#---integrator params
tspan = [0.0, t_trans+t_exec]
p = [v_τ_m, N, v_I, m3_adj, v_inh, v_t_refrac, tt, v_numSpikes, m_tS_ext2]



#------define solver
println("Now I'll start integrating version " * string(code_version) * ". I have Ni = " * string(Ni) *", Ne = " * string(Ne) *". Running for t_exec = " * string(t_exec) * ", with t_trans = " * string(t_trans) * ".")
cb = VectorContinuousCallback(reset_cond, reset_affect!, N)
prob = ODEProblem(LIF!, u0, tspan, p)
# sol = @time solve(prob, Tsit5(), callback=cb)
# sol = @time solve(prob, BS3(), callback=cb)
sol = @time solve(prob, Euler(), callback=cb, adaptative=false, dt=0.01, dense=false)

println("Finished integration")





tt = [filter!(x->(x≠0.0 && x>t_trans), tt[i]) for i=1:N]
fileOut = "results/spikeTimes/spike_times_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_conn_prob_" * string(conn_prob) * "_seed_" * string(seed) *"_conductances_2.dat"
saveSpikeTimes(fileOut, tt)


# m_SR, m_t = calcSpikeRate(m_tS_ext, 1, 0, 50)
# plotSpikeRate(m_tS_ext, 1, 0, t_exec) #ruido
# plotSpikeRate(tt, 1, 0, t_exec) #media da rede
# plotSpikeRate_IE(tt, 1, t_trans, t_exec, v_inh)

# fileOut = "results/potentials/test2_potentials_LIF_v" * code_version * "_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_seed_" * string(seed) * "_conductances_1.png"
# plotPotentials(sol, [], fileOut)
plotPotentials(sol)
# plotPotentialNeuron(sol, 1)
# plotRaster(tt, v_inh)
