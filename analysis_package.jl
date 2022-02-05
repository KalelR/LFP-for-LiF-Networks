using DSP 

function calcSpikeRate(m_tS, binSize, t_start, t_end)
    m_SR, m_t = [], []
    for v_tS in m_tS
        # bins = [v_tS[1]+i*binSize for i=0:Int64(floor((v_tS[end]-v_tS[1])/binSize))]
        bins = [t_start+i*binSize for i=0:Int64(floor((t_end-t_start)/binSize))]
        hist = fit(Histogram, v_tS, bins)
        push!(m_SR, hist.weights)
        push!(m_t, hist.edges)
    end
    return m_SR, m_t
end


function calcLFP(m_Ia, m_Ig)
    println(length(m_Ia[1,:]))
    v_sum_Ia = abs.([sum(m_Ia[:,i]) for i=1:length(m_Ia[1,:])])
    v_sum_Ig = abs.([sum(m_Ig[:,i]) for i=1:length(m_Ig[1,:])])

    v_LFP = v_sum_Ia + v_sum_Ig
end


function analyseLFP(LFP, dt)
    fs = 1/dt;
    signal = deepcopy(LFP);
    # signal .-= mean(signal)
    pdg = DSP.Periodograms.periodogram(signal, fs=fs)

    v_freqs = freq(pdg) .* 1000 #Hz
    psd = power(pdg)
    psd = psd ./ maximum(psd)

    plot(v_freqs, psd, title = "Power Spectrum", yscale=:log10, yguide="LFP power", xguide="frequency (Hz)")
    # scatter(v_freqs, psd, title = "Power Spectrum", yscale=:log10, yguide="LFP power", xguide="frequency (Hz)")
    # plot(freqs[2:end], psd[2:end], title = "Power Spectrum", yscale=:log10, yguide="LFP power", xguide="frequency (Hz)", xscale=:log10)
    plot!(xlim=(0, 150)) #Hz
end