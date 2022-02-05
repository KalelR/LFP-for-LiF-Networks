using DSP
using DifferentialEquations
using Random
# using FFTW
using Plots
pyplot()
using Statistics

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
    # σ = 1.5 #spikes/ms

    θ = (1/τ_n)
    μ = 0.
    # σ = (σ/τ_n)*sqrt(2/τ_n)
    t0 = t_start = 0.
    W0 = 0.


    p = [τ_n, σ]
    prob = ODEProblem(OUP, [0.0],  (t_start, t_end), p)
    sol = solve(prob)
    plot(sol.t,sol.u)


    v_t, n = sol.t, sol.u
    v_t = [round(t; digits=num_digits) for t in v_t] #
    if (v_t[end] != t_end) pop!(v_t); pop!(n) end #now length is t_end/dt + 1


    # mean(sol.u)
    return v_t, n
end


τ_n = 16.0 #ms
σ = 0.4 #spikes/ms
p = [τ_n, σ]
t_start = 0.
t_end = 10^5
prob = ODEProblem(OUP, [0.0],  (t_start, t_end), p)
dt = 0.01
sol = solve(prob, Euler(), adaptative=false, dt=dt)
# sol = solve(prob)

# plot(sol)
# savefig("exemplo_ruido_OUP_gera3.png")

# signal = [sol.u[i][1] for i=1:length(sol.u)]

signal = [(sol.u[i][1])^2 for i=1:length(sol.u)]
signal = signal .- mean(signal)

# 
fs = 1/dt;
pdg = DSP.Periodograms.periodogram(signal, fs=fs)

freqs = freq(pdg) .* 1000 #Hz
psd = power(pdg)

plot(freqs, psd, title = "Power Spectrum")
plot!(xlim=(0, 50)) #Hz



savefig("power_spectrum_OUP_3.png")
