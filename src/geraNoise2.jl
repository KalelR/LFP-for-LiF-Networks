using DifferentialEquations
using Plots
pyplot()
using Statistics
using FFTW

function OUP(u,p,t)
    # θ = (1/16) #1/16ms  
    # μ = 0
    # σ = 0.4*sqrt(2*16) 
    τ_n = 16 #ms
    du = (-u)/τ_n; 
end

function OUP2(u,p,t)
     τ_n = 16 #ms
    σ_n = 0.4
    du = (σ_n/τ_n)*sqrt(2/τ_n)
end
t0 = t_start = 0.
t_end = 100000.

prob = SDEProblem(OUP, OUP2, 0.0, (t_start, t_end))
sol = solve(prob, dt = 0.1)

# plot(sol.t, sol.u)

mean(sol.u)

Ts = 0.1 #sampling period
signal = sol.u .- mean(sol.u)
F = fft(signal) |> fftshift
freqs = fftfreq(Int64(t_end/0.1)+1, 1.0/Ts) |> fftshift

freq_domain = plot(freqs, abs.(F), title = "Spectrum") 
