using DifferentialEquations
using Plots
using DifferentialEquations
using Plots
pyplot()
# using Distributions
# using Random
include("LIF_package.jl")
include("network_package.jl")
include("noise_package.jl")



#------network params
Ni = 10; Ne = 40;
N = Ni + Ne
conn_prob = 0.2;


#----------simulation params
const t_exec = 30.0; #ms
const t_trans = 0.; #ms
seed = 0
code_version = "0.1"
cond_version = "1"

#--model params
v_τ_m = [20., 10.] #EXC = [1]; INH = [2]
const vr = 11.0 #mV, reset potential
const vth = 18.0 #mV
const I = 0.;
v_t_refrac = [2.0, 1.0] #vector with the refractory periods for INH [2] and EXC[1] neurons
v_I = collect(range(I, I, length=N))
u0 = collect(range(vr+1, vth-1, length=N))


#--auxiliary
max_num_spikes = Int64((t_exec-t_trans)/minimum(v_t_refrac))*10 #maximum number of spikes that a neuron may fire#could be more memory conversative by using refrac period of each neuron
tt =  [ [0.0 for i=1:max_num_spikes] for j=1:N]
v_numSpikes = zeros(Int64, N) #counter of number of spikes at current time t
v_numSpikes_ext = ones(Int64, N) #counter of number of spikes at current time t


#---connection_params
const m3_adj, v_inh, v_exc = read_connection_parameters(Ni, Ne, conn_prob, seed, cond_version);

#--- external input
ν0 = 0.05 #spikes/ms
t_end = t_exec + t_trans
dt = 0.1 #ms ----------- careful
num_digits = 1;
m_tS_ext = generateExternalInput(N, t_end, dt, ν0, num_digits)
v_hasSpiked = [false for i=1:N]

#---integrator params
tspan = [0.0, t_trans+t_exec]
p = [v_τ_m, N, v_I, m3_adj, v_inh, v_t_refrac, tt, v_numSpikes, m_tS_ext, v_numSpikes_ext, v_hasSpiked]



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


# plot(sol, marker=:circle, markerstrokewidth=0, alpha=0.8, markersize=2)
# v_aux = [zeros(length(v_t))*vth for v_t in tt ]
# for i=1:N
#     plot!(tt[i], v_aux[i], marker=:circle, markersize=10, linewidth=0)
# end
# plot!(size=(1200, 600), legend=false, xguide="t", yguide="V_i(t)", ylim = (vr-2, vth+2))

# savefig("results/potentials/LIF_v" * code_version * "_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_seed_" * string(seed) * "_conductances_1.png")


# i = 40;
# v_j = connectedToWhom(i, m_m_adj, N)

# plot(sol.t, sol[i,:])
# for j in v_j
#     println(j)
#  plot!(sol.t, sol[j,:])
# end
# gui()


# scatter()
# for i=1:N
#     t = tt[i]
#     v_aux = collect(range(i,i,length=length(t)))
#     if i ∈ v_inh mk_color = "blue" else mk_color = "red" end
#     scatter!(t, v_aux, markersize=2, markerstrokewidth=0, markercolor=mk_color, legend=false)
# end

# scatter!(xlims=(0, t_exec+t_trans), ylims=(0, N), xguide="t", yguide="neuron #")


# savefig("results/RP/LIF_v" * code_version * "_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_seed_" * string(seed) * "_conductances_2.png")
