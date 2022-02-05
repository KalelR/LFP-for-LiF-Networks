using DifferentialEquations
using Plots
pyplot()
using Statistics
# using FFTW
using DSP

# τ_n = 16.0/1000.0 #s
# σ_2 = 0.4*1000.0 #spikes/s
# dt = 0.001;
# t_end = 100000.

τ_n = 16.0 #ms
σ_2 = 0.4 #spikes/ms
dt = 0.1;
t_end = 10^6


θ = (1/τ_n)
μ = 0.
σ = (σ_2/τ_n)*sqrt(2/τ_n) 
W0 = nothing
t0 = t_start = 0.
W0 = 0.

OHP =  OrnsteinUhlenbeckProcess(θ,μ,σ,t0,W0)
prob = NoiseProblem(OHP, (t_start, t_end))
sol = solve(prob, dt = dt)

plot(sol.t, sol.u)
savefig("exemplo_ruido_OUP.png")
mean(sol.u)


# signal = sol.u .- mean(sol.u);
signal = sol.u
#--- METODO 1 (FFT)
# Ts = dt; #sampling period
# F = fft(signal) |> fftshift;
# freqs = fftfreq(Int64(t_end/0.1)+1, 1.0/Ts) |> fftshift;

# freq_domain = plot(freqs, abs.(F), title = "Spectrum") 



#--- METODO 2 (PSD)
fs = 1/dt;
p = DSP.Periodograms.periodogram(signal, fs=fs)

freqs = freq(p) .* 1000 #Hz
psd = power(p)

plot(freqs, psd, title = "Power Spectrum")
plot!(xlim=(0, 50)) #Hz

savefig("psd_OUP_second.png")

