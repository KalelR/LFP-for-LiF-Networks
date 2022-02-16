import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
import scipy.stats.distributions as dist

import nitime.algorithms as tsa
import nitime.utils as utils
from nitime.viz import winspect
from nitime.viz import plot_spectral_estimate


def dB(x, out=None):
    if out is None:
        return 10 * np.log10(x)
    else:
        np.log10(x, out)
        np.multiply(out, 10, out)

def spectralAnalysis(signal, dt):
    f, adaptive_psd_mt, nu = tsa.multi_taper_psd(signal, Fs=fs,  adaptive=True, jackknife=False)
    freqs, psd, adaptive_jk_var = tsa.multi_taper_psd(signal, Fs = fs, adaptive=True, jackknife=True)
    Kmax = max(nu)
    psd /= max(psd)
    freqs *= 1000; 
    return freqs, psd

dt = 0.05
fs = 1/dt 
v_rates = ["1.2", "2.4"]
code_version = "2.1"
Ni = 1000
Ne = 4000
conn_prob = 0.2
seed = 0
f0 = 1.2

for rate in v_rates:
    fileIn = "results/lfp/LFP_LIF_v" + str(code_version) + "_Ni_" + str(Ni) + "_Ne_" + str(Ne) + "_conn_prob_" + str(conn_prob) + "_seed_" + str(seed) +"_conductances_2_rate_" + str(rate) + ".dat"
    fileOut = "results/lfp/LFP_spectral_analysis_LIF_v" + str(code_version) + "_Ni_" + str(Ni) + "_Ne_" + str(Ne) + "_conn_prob_" + str(conn_prob) + "_seed_" + str(seed) +"_conductances_2_rate_" + str(rate) + ".dat"
    data = np.loadtxt(fileIn)
    v_t = data[0,:]
    signal = data[1,20000:]
    freqs, psd = spectralAnalysis(signal, dt)
    plt.semilogy(freqs, psd,label=rate)
    np.savetxt(fileOut, [freqs, psd])




plt.xlim(0, 150)
plt.ylim(10**(-3), 1)
plt.ylabel('LFP Power')
plt.xlabel('Frequency (Hz)')
plt.legend()
plt.savefig(fileIn.replace(".dat", ".png"))

