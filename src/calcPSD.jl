# calcula power spectrum density

using Plots
pyplot()
using Statistics
using DSP


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

fs = 1/Ts
p = DSP.Periodograms.periodogram(signal, fs=fs)

freqs = freq(p)
psd = power(p)

plot(freqs, psd, title = "Power Spectrum") #-- pico em f = 60 Hz!