using StatsBase
using Plots
pyplot()

#function to calc Spike Rate com m_tS, matrix containging spike tiemes; also receives the binSize to calc spike rate
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


function plotSpikeRate(m_tS, binSize, t_start, t_end, fileOut="")
    m_SR, m_t = calcSpikeRate(m_tS, binSize, t_start, t_end);
    #average
    v_average_SR = mean(m_SR)
    println("media = ", mean(v_average_SR))
    plot(m_t[1][1][1:end-1], v_average_SR, label="input", color=:black)
    plot!(xguide="t", yguide="Spike rate")

    # for i=1:length(m_SR)
    #for each neuron

    #     plot(m_t[i], m_SR[i])
    # end
    if(fileOut != "")
        savefig(fileOut)
    end
end


function plotSpikeRate_IE(m_tS, binSize, t_start, t_end, v_inh, fileOut, ymin1=0, ymax1=0, ymin2=0, ymax2=0)
    m_tS_aux = [m_tS[i]  for i=1:length(m_tS) if i in v_inh ]
    m_SR, m_t = calcSpikeRate(m_tS_aux, binSize, t_start, t_end);
    #average
    v_average_SR = mean(m_SR)
    println("media = ", mean(v_average_SR))
    a=plot(m_t[1][1][1:end-1], v_average_SR*1000, label="interneuron", color=:blue)
    plot!(xguide="t", yguide="Spike rate(spikes/sec)")
    if(ymin1!=0 || ymax1!=0)    plot!(ylims=(ymin1, ymax1)) end


    m_tS_aux = [m_tS[i] for i=1:length(m_tS) if !(i in v_inh)  ]
    m_SR, m_t = calcSpikeRate(m_tS_aux, binSize, t_start, t_end);
    #average
    v_average_SR = mean(m_SR)
    println("media = ", mean(v_average_SR))
    b=plot(m_t[1][1][1:end-1], v_average_SR*1000, label="pyramidal", color=:red)
    plot!(xguide="t", yguide="Spike rate (spikes/sec)")
    if(ymin2!=0 || ymax2!=0)    plot!(ylims=(ymin2, ymax2)) end
    

    plot(a,b, layout=(2,1))
    if(fileOut != "")
        savefig(fileOut)
    end

end


function plotPotentials(v_t, m_V, tt=[], fileOut="")
    plot(v_t, m_V, marker=:circle, markerstrokewidth=0, alpha=0.8, markersize=2)
    plot!(legend=false, xguide="t", yguide="V_i(t)", ylim = (0, 20))
    # if(tt != [])
    #     v_aux = [zeros(length(v_t))*vth for v_t in tt ]
    #     for i=1:N
    #         plot!(tt[i], v_aux[i], marker=:circle, markersize=10, linewidth=0)
    #     end
    # end
    # if(fileOut!="")
    #     savefig(fileOut)
    # end
    # savefig("results/potentials/LIF_v" * code_version * "_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_seed_" * string(seed) * "_conductances_1.png")
end


function plotPotentials(sol, tt=[], fileOut="")
    plot(sol, marker=:circle, markerstrokewidth=0, alpha=0.8, markersize=2)
    plot!(legend=false, xguide="t", yguide="V_i(t)", ylim = (0, 20))
    # if(tt != [])
    #     v_aux = [zeros(length(v_t))*vth for v_t in tt ]
    #     for i=1:N
    #         plot!(tt[i], v_aux[i], marker=:circle, markersize=10, linewidth=0)
    #     end
    # end
    # if(fileOut!="")
    #     savefig(fileOut)
    # end
    # savefig("results/potentials/LIF_v" * code_version * "_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_seed_" * string(seed) * "_conductances_1.png")
end

function plotPotentialNeuron(sol, i)
    V = [sol[j][i] for j =1:length(sol)]
    # plot(sol.t, V)
    plot(sol.t, V/maximum(V))
    if(m_Iext!=[])
        plot!(sol.t, m_Iext[i]/maximum(m_Iext[i]), alpha=0.6)
    end
end



function plotRaster(tt, v_inh, t_trans, t_exec, fileOut)
    scatter()
    N = length(tt)
    for i=1:N
        t = tt[i]
        v_aux = collect(range(i,i,length=length(t)))
        if i ∈ v_inh mk_color = "blue" else mk_color = "red" end
        scatter!(t, v_aux, markersize=2, markerstrokewidth=0, markercolor=mk_color, legend=false)
    end

    scatter!(xlims=(t_trans, t_exec+t_trans), ylims=(0, N), xguide="t", yguide="neuron #")
    if(fileOut != "")
        savefig(fileOut)
    end
end



function plotRasterMoreActive(tt, v_inh, numActive, t_exec, fileOut)
    N = length(tt)
    m_SR, m_t = calcSpikeRate(tt, 1, 0, t_exec)
    v_average_SR = [mean(v_SR) for v_SR in m_SR] #average SR of each neuron
    # println(N, m_SR, v_average_SR)
    scatter()
    cont = 1
    for i=1:N
        if (v_average_SR[i] == 0) continue end
        t = tt[i]
        v_aux = collect(range(cont,cont,length=length(t)))
        if i ∈ v_inh mk_color = "blue" else mk_color = "red" end
        scatter!(t, v_aux, markersize=2, markerstrokewidth=0, markercolor=mk_color, legend=false)
        cont += 1;
    end

    # scatter!(xlims=(0, t_exec+t_trans), ylims=(0, N), xguide="t", yguide="neuron #")
    # if(fileOut != "")
        # savefig(fileOut)
    # end
end
