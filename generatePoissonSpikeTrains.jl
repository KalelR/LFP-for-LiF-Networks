using DifferentialEquations
using Plots
pyplot()
using Random
using Statistics
# using DSP

function generateNoise(dt, t_end)
    τ_n = 16.0 #ms
    σ = 0.4 #spikes/ms

    θ = (1/τ_n)
    μ = 0.
    σ = (σ/τ_n)*sqrt(2/τ_n) 
    t0 = t_start = 0.
    W0 = 0.

    OHP =  OrnsteinUhlenbeckProcess(θ,μ,σ,t0,W0)
    prob = NoiseProblem(OHP, (t_start, t_end))
    sol = solve(prob, dt = dt)

    mean(sol.u)
    return sol.t, sol.u
end

#functio nto generatePoissonSpikeTrain, r is the vector containing the instantaneous rates and dt is the duration of each time bin, t_end is the total durationm 
#the spike train is generated for each neuron in the network
#returns the matrix in which each row contains the spike train for the i-th neuron
function generatePoissonSpikeTrain(r, dt, t_end, N)
    num_bins = floor(Int64(t_end/dt))
    m_t_spikes = [ Float64[] for i=1:N ]
    for i = 1:N #for each neuron 
        x = rand(Float64, (1, num_bins))
        for j=1:num_bins
            t = (j-1)*dt; #the time of the bin      
            if (x[j] ≤ r[j]*dt) 
                push!(m_t_spikes[i], t)
            end
        end
    end
    return m_t_spikes
end

# ---paramters
dt = 0.05; #ms 
t_end = 10^4; #ms
N = 5000;
num_digits = 2 #number of digits in dt

#--noise info
v_t, n = generateNoise(dt, t_end)
v_t = [round(t; digits=num_digits) for t in v_t] #
if (v_t[end] != t_end) pop!(v_t); pop!(n) end #now length is t_end/dt + 1

#--signal
ν0 = 1.2 #spikes/ms
ν_s = [ν0 for i=0:dt:t_end]

#--rate 
ν = ν_s
# ν = ν_s + n
ν = [x > 0 ? x : 0  for x in ν]

#--calculate spikes 
m_t_spikes = generatePoissonSpikeTrain(ν, dt, t_end, N)


#-----to validade method 


# plot spike times 
# m_idx_neurons = [ [i for j=1:length(m_t_spikes[i])] for i=1:N]
# plot(m_idx_neurons, m_t_spikes, marker=:hline)
# plot(m_t_spikes, m_idx_neurons, marker=:vline)
# savefig("PoissonSpikeTrains_rasterPlot_1.png")


# spike count histogram 
num_spikes = [length(m_t_spikes[i]) for i=1:N]
histogram(num_spikes, bins=100)
savefig("PoissonSpikeTrains_spikeCount_histogram_2.png")



#--ISI histogram
v_ISI = [v[i+1] - v[i] for v in m_t_spikes for i=1:(length(v)-1)]
Cv = std(v_ISI)/mean(v_ISI) #should be 1
histogram(v_ISI, title="CV = " * string(Cv), bins=200)
# savefig("results/tests/noise/PoissonSpikeTrains_histogram.png")
savefig("PoissonSpikeTrains_histogram_2.png")