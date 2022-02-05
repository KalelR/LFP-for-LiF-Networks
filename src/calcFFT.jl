using Plots
using FFTW

# Number of points 
N = 2^14 - 1 
# Sample period
Ts = 1 / (1.1 * N) 
# Start time 
t0 = 0 
tmax = t0 + N * Ts
# time coordinate
t = t0:Ts:tmax

# signal 
signal = sin.(2π * 60 .* t) # sin (2π f t) 


F = fft(signal) |> fftshift
freqs = fftfreq(length(t), 1.0/Ts) |> fftshift

freq_domain = plot(freqs, abs.(F), title = "Spectrum", xlim=(-1000, +1000))
