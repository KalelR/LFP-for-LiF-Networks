# using DSP
using DifferentialEquations
using Random
# include("plot_package.jl")

function generateNoise(dt, t_end, num_digits)
    τ_n = 16.0 #ms
    σ = 0.4 #spikes/ms
    # σ = 4.0 #spikes/ms

    θ = (1/τ_n)
    μ = 0.
    σ = (σ/τ_n)*sqrt(2/τ_n)
    t0 = t_start = 0.
    W0 = 0.

    OHP =  OrnsteinUhlenbeckProcess(θ,μ,σ,t0,W0)
    prob = NoiseProblem(OHP, (t_start, t_end))
    sol = solve(prob, dt = dt)

    v_t, n = sol.t, sol.u
    v_t = [round(t; digits=num_digits) for t in v_t] #
    if (v_t[end] != t_end) pop!(v_t); pop!(n) end #now length is t_end/dt + 1


    # mean(sol.u)
    return v_t, n
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

function generateExternalInput(N, t_end, dt, ν0, num_digits, seed)
    Random.seed!(seed)
    #--noise info
    v_t, n = generateNoise(dt, t_end, num_digits)
    # n = n[1:(Int64(floor(t_end/dt)))]
    #--signal
    ν_s = [ν0 for i=0:dt:t_end]

    min_len = min(length(n), length(ν_s))
    n = n[1:min_len]; ν_s = ν_s[1:min_len]     
    # println(length(ν_s), " ", length(n))
    #--rate
    ν = ν_s + 10*n
    ν = [x > 0 ? x : 0  for x in ν]

    #--calculate spikes
    m_t_spikes = generatePoissonSpikeTrain(ν, dt, t_end, N)
    return m_t_spikes
end


function OUP(du,u,p,t)
    τ, σ = p
    η = randn()
    # η = 1.0
    du[1] = (-u[1] + σ*√(2/τ)*η)/τ
    return du
end

function generateNoise2(dt, t_end, num_digits)
    τ_n = 16.0 #ms
    σ = 0.4 #spikes/ms

    t_start = 0.
    p = [τ_n, σ]
    prob = ODEProblem(OUP, [0.0],  (t_start, t_end), p)
    sol = solve(prob, Euler(), adaptative=false, dt=dt)

    # sol = solve(prob)
    # plot(sol.t,sol.u)


    v_t, n = sol.t, [sol.u[i][1] for i=1:length(sol.u)]
    v_t = [round(t; digits=num_digits) for t in v_t] #
    if (v_t[end] != t_end) pop!(v_t); pop!(n) end #now length is t_end/dt + 1


    # mean(sol.u)
    return v_t, n
end
function generateExternalInput2(N, t_end, dt, ν0, num_digits, seed)
    Random.seed!(seed)
    #--noise info
    v_t, n = generateNoise2(dt, t_end, num_digits)
    n = n[1:end-1]
    #--signal
    ν_s = [ν0 for i=0:dt:t_end]

    #--rate
    # println(length(ν_s), " ", length(n))
    ν = ν_s + n
    ν = [x > 0 ? x : 0  for x in ν]

    #--calculate spikes
    m_t_spikes = generatePoissonSpikeTrain(ν, dt, t_end, N)
    return m_t_spikes
end
