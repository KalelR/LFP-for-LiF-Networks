
#--funçao para calcular condutancia dupla exponencial
g(t, t0, g0, tau_d, tau_r, D) = g0*( exp(-(t-t0-D)/tau_d) - exp(-(t-t0-D)/tau_r) )

function calc_gij(t::Float64, tt::Array{Array{Float64, 1}, 1}, v_numSpikes::Array{Int64, 1}, j::Int64, D::Float64, g0::Float64, tau_d::Float64, tau_r::Float64)
    g_ij = 0.0::Float64
    if v_numSpikes[j]::Int64 > 0
        t0 = tt[j][v_numSpikes[j]]
        if (t-t0-D) >= 0 #--considers delay
            g_ij = g(t, t0, g0, tau_d, tau_r, D)
        end
    end
    return g_ij
end
function calc_gij!(t::Float64, t0::Float64, D::Float64, g0::Float64, tau_d::Float64, tau_r::Float64)
    g_ij = 0.0::Float64
        if (t-t0-D) >= 0 #--considers delay
            # if (t-t0-D > 1e-4 && j == 1) println(t-t0-D, " ", t, " ", t0) end
            g_ij = g(t, t0, g0, tau_d, tau_r, D)
            # println("t = ", t, "gij =  ",  g_ij, ", g0 = ", g0, ", (t-t0-D) =  ", (t-t0-D), " ", tau_d, " ", tau_r, " ", t0)
        end
    return g_ij
end

# m_Iext = [ [] for i=1:N]
# m_t = [ [] for i=1:N]

#funçao para calcular a corrente sinaptica chengaod no neuronio i; recebe o indice i, o potencial V(i), a matriz das matrizes de adj, o tempo e a matriz com os tempso dos spikes
function calcIsyn(i::Int64, N::Int64, V::Float64, m3_adj::Array{ Array{Array{Int64,1},1}, 1},  t::Float64, tt::Array{Array{Float64, 1}, 1}, v_numSpikes::Array{Int64, 1}, m_tS_ext::Array{Array{Float64,1},1})
    Isyn = 0.0

    #--for recurrent connections (inside)
    for idx_syn = 1:2    #loop sobre os tipos de receptors (1 = ampa, 2 = gaba)
        v_adj = m3_adj[idx_syn][i] #seleciona a matriz de adjacencia (com os k e o g_0 associado a cada sinapse)

        isExc = i in v_inh ? false : true
        #--choose synaptic times and synaptic efficacies
        if (idx_syn) == 1
            τ_r, τ_d = isExc ?  (0.4, 2.) : (0.2 , 1.)
            J = isExc ? 0.42 : 0.7
        else #GABA
            τ_r, τ_d = 0.25, 0.5
            J = isExc ? -1.7*0.5 : -2.7*1.5
        end

        τ_m = isExc ? 20 : 10
        g0 = (J*τ_m)/(τ_d - τ_r)
        Er = 0.;
        τ_L = 1. #ms

        for j ∈ v_adj
            if(j == -1) continue end #continue e break dao na mesma aqui
            g_ij = @fastmath calc_gij(t, tt, v_numSpikes, j, τ_L, g0, τ_d, τ_r)::Float64
            I_ij = g_ij*(V - Er)
            Isyn += I_ij;
        end
        # push!(v_I[1], Isyn)

    end
    #---now for external input
    isExc = i in v_inh ? false : true
    τ_r, τ_d = isExc ?  (0.4, 2.) : (0.2 , 1.)
    τ_L = 1.; Er = 0.;
    J = isExc ? 0.55 : 0.95
    τ_m = isExc ? 20 : 10
    g0 = (J*τ_m)/(τ_d - τ_r)
    g_ij = @fastmath calc_gij!(t, m_tS_ext[i][1], τ_L, g0, τ_d, τ_r) #calc g_ij, but make sure matrix has values
    I_ij = g_ij*(V - Er)
    # push!(m_Iext[i], I_ij)
    # push!(m_t[i], I_ij)
    Isyn += I_ij;
    # Isyn = I_ij; #only external input
    # if (i == 1) println(t, " "z, m_tS_ext[i][1], " ", t-m_tS_ext[i][1]-τ_L, " ", i," ", Isyn, " ", V, " ", g_ij, " ", (V - Er) ) end
    # if (i == 1) println(t, " ", Isyn, " ", V, " ", g_ij, " ", (V - Er) ) end
    # Isyn = 20.
    return Isyn
end



function LIF(du, u, p::Array{Any,1}, t::Float64)
    v_τ_m::Array{Float64, 1}, N::Int64, v_I::Array{Float64, 1}, m3_adj::Array{ Array{Array{Int64,1},1}, 1},  v_inh::Array{Int64, 1}, v_t_refrac::Array{Float64, 1}, tt::Array{Array{Float64, 1}, 1}, v_numSpikes::Array{Int64, 1},  m_tS_ext::Array{Array{Float64, 1}, 1}= p #alocar essa linha dá 3 alocaçoes e 48 bytes
    τ_L = 1.0 #ms
    @simd for i = 1:N
        if i ∈ v_inh idx = 2 else idx = 1 end #determine if neron is I (= 2) or E (= 1)
        t_refrac, τ_m = v_t_refrac[idx], v_τ_m[idx]
        @inbounds if v_numSpikes[i] > 0 t0 = tt[i][v_numSpikes[i]]::Float64 else t0 = -(t_refrac+1.0) end #find time of last spike of ith neuron

        if(length(m_tS_ext[i]) >= 2) if(t >= m_tS_ext[i][2] + τ_L) popfirst!(m_tS_ext[i]) end end #update matrix so that first column contains current external spike times; stops at last spike

        if (t >= (t0 + t_refrac))
            @inbounds V, I = u[1*i], v_I[1*i]
           # @inbounds @fastmath Isyn = calcIsyn(i, N, V, m3_adj, t, tt, v_numSpikes, m_tS_ext)
            # @inbounds  @fastmath du[i] = (-gl*(V - vr) + I - Isyn)/C
            Isyn = 20.
            @inbounds  @fastmath du[i] = ( -V + Isyn  )/τ_m
            # if (i == 1) println(t, " ", i, " ", du[i], " ", Isyn, " ", V) end
        else
            du[i] = 0
        end
    end
end

#condition for reset acontece qnd v > v_th
function reset_cond(out, u, t, integrator)
    vth = 18.
    @simd for i = 1:N
        # println(t, " ", u[i])
        @inbounds out[i] = vth - u[1*i]
    end
end

#effect of reset
function reset_affect!(integrator, event_index)
    # println(integrator.t, " ", event_index)
    # push!(tt[1*event_index], integrator.t) #save the spike times
    @inbounds v_numSpikes[event_index] += 1;
    @inbounds tt[1*event_index][v_numSpikes[event_index]] = integrator.t
    @inbounds integrator.u[1*event_index] = 11.0 #reset the potential
end



#============ for tests?
du = u = similar(u0)
t = 0.
@time LIF!(du,u,p,t)
 @time calcIsyn(1, N, 15., m3_adj, t, tt, v_numSpikes, m_tS_ext)
# =#
# function calcIsyn(i::Int64, N::Int64, V::Float64, m3_adj::Array{ Array{Array{Int64,1},1}, 1},  t::Float64, tt::Array{Array{Float64, 1}, 1}, v_numSpikes::Array{Int64, 1}, m_tS_ext::Array{Array{Float64,1},1})
#     Isyn = 0.0

#     #--for recurrent connections (inside)
#     for idx_syn = 1:2    #loop sobre os tipos de receptors (1 = ampa, 2 = gaba)
#         v_adj = m3_adj[idx_syn][i] #seleciona a matriz de adjacencia (com os k e o g_0 associado a cada sinapse)

#         isExc = i in v_inh ? false : true
#         #--choose synaptic times and synaptic efficacies
#         if (idx_syn) == 1
#             τ_r, τ_d = isExc ?  (0.4, 2.) : (0.2 , 1.)
#             g0 = isExc ? 0.42 : 0.7
#         else #GABA
#             τ_r, τ_d = 0.25, 0.5
#             g0 = isExc ? -1.7 : -2.7
#         end

#         Er = 0.;
#         τ_L = 1. #ms

#         for j ∈ v_adj
#             if(j == -1) continue end #continue e break dao na mesma aqui
#             g_ij = @fastmath calc_gij(t, tt, v_numSpikes, j, τ_L, g0, τ_d, τ_r)::Float64
#             I_ij = g_ij*(V - Er)
#             Isyn += I_ij;
#         end

#         if(idx_syn == 1) Ia = Isyn else Ig = Isyn - Ia end

#     end
    # push!(m_I[1][i], Ia)
    # push!(m_I[2][i], Ig)
    # push!(v_t, t)

#     #---now for external input
#     isExc = i in v_inh ? false : true
#     τ_r, τ_d = isExc ?  (0.4, 2.) : (0.2 , 1.)
#     τ_L = 1.; Er = 0.;
#     g0 = isExc ? 10*0.55 : 0.95
#     g_ij = @fastmath calc_gij!(t, m_tS_ext[i][1], τ_L, g0, τ_d, τ_r) #calc g_ij, but make sure matrix has values
#     I_ij = g_ij*(V - Er)
#     Isyn += I_ij;
#     # if (i == 1) println(t, " ", m_tS_ext[i][1], " ", t-m_tS_ext[i][1]-τ_L, " ", i," ", Isyn, " ", V, " ", g_ij, " ", (V - Er) ) end
#     return Isyn
# end
