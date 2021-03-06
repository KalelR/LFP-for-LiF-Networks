#--funçao para calcular condutancia dupla exponencial
g(t, t0, g0, tau_d, tau_r, D) = g0*( exp(-(t-t0-D)/tau_d) - exp(-(t-t0-D)/tau_r) )

#--second way, considers all previous spikes; for internal connections
# calculates conductance for all 30 previous spikes of neuron j
function calc_gij(t::Float64, tt::Array{Array{Float64, 1}, 1}, v_numSpikes::Array{Int64, 1}, j::Int64, D::Float64, g0::Float64, tau_d::Float64, tau_r::Float64, idx::Int64)
    g_ij = g_aux =  0.0::Float64
    if v_numSpikes[j]::Int64 > 0 #v_numSpikes[i] = number of spikes of unit i
        for idx_sp = v_numSpikes[j]:-1:v_numSpikes[j]-30
            if(idx_sp <= 0) break end #in case number of spikes is smaller than 30
            t0 = tt[j][idx_sp] #time of spike idx_sp of unit j
            if (t >= t0+D) #--considers delay
                g_aux = g(t, t0, g0, tau_d, tau_r, D)
                g_ij += g_aux;
            end
        end
    end
    return g_ij
end

#--second way, considers all previous spikes; for external input
function calc_gij2!(t::Float64, tt::Array{Array{Float64,1}, 1}, v_numSpikes::Array{Int64,1}, i::Int64, D::Float64, g0::Float64, tau_d::Float64, tau_r::Float64)
    g_ij = g_aux = 0.0::Float64
    if v_numSpikes[i]::Int64 > 0
        for idx_sp = v_numSpikes[i]:-1:v_numSpikes[i]-30
            if(idx_sp <= 0) break end
            t0 = tt[i][idx_sp];
            if (t-t0-D) >= 0 #--considers delay
                g_aux = g(t, t0, g0, tau_d, tau_r, D)
                g_ij += g_aux;
            end
        end
    end
    return g_ij
end

#funçao para calcular a corrente sinaptica chengaod no neuronio i; recebe o indice i, o potencial V(i), a matriz das matrizes de adj, o tempo e a matriz com os tempso dos spikes, idxf: how many times this function has been called for the same t (in 4th order RK4 it is called 4 times per t)
function calcIsyn!(i::Int64, N::Int64, V::Float64, m3_adj::Array{ Array{Array{Int64,1},1}, 1},  t::Float64, tt::Array{Array{Float64, 1}, 1}, v_numSpikes::Array{Int64, 1}, m_tS_ext::Array{Array{Float64,1},1}, v_numS_ext::Array{Int64,1}, m_Ia::Array{Float64,2}, m_Ig::Array{Float64,2}, aux::Int64, idxf::Int64)
# function calcIsyn!(i::Int64, N::Int64, V::Float64, m3_adj::Array{ Array{Array{Int64,1},1}, 1},  t::Float64, tt::Array{Array{Float64, 1}, 1}, v_numSpikes::Array{Int64, 1}, m_tS_ext::Array{Array{Float64,1},1})
    Isyn = 0.0
    Ia = 0.; Ig =  0.
    #--for recurrent connections (inside)
    for idx_syn = 1:2    #loop sobre os tipos de receptors (1 = ampa, 2 = gaba)
        v_adj = m3_adj[idx_syn][i] #seleciona a matriz de adjacencia do tipo de receptor (com os k e o g_0 associado a cada sinapse) 

        isExc = i in v_inh ? false : true
        #--choose synaptic times and synaptic efficacies
        if(idx_syn == 1) #AMPA
            τ_r, τ_d = isExc ?  (0.4, 2.) : (0.2 , 1.) #excitatory: talr = 0.4, tald = 2; inhibitory: talr = 0.2, tald = 1.
            J = isExc ? 0.42 : 0.7 #extitory: J = 0.42; inhibotiroy: J = 0.7
        else #GABA
            τ_r, τ_d = 0.25, 5.
            J = isExc ? -1.7 : -2.7
        end

        τ_m = isExc ? 20 : 10
        g0 = (J*τ_m)/(τ_d - τ_r)
        Er = 0.;
        τ_L = 1. #ms

        for j ∈ v_adj #coupling loop
            if(j == -1) continue end #continue e break dao na mesma aqui
            g_ij = @fastmath calc_gij(t, tt, v_numSpikes, j, τ_L, g0, τ_d, τ_r, idxf)::Float64
            # I_ij = g_ij*(V - Er)
            I_ij = g_ij
            Isyn += I_ij;
        end
        #save Ia and Ig for LFP
        if(idx_syn == 1) 
             Ia = Isyn 
        else 
                Ig = Isyn - Ia
        end


    end
    
    #---now for external input
    isExc = i in v_inh ? false : true
    τ_r, τ_d = isExc ?  (0.4, 2.) : (0.2 , 1.)
    τ_L = 1.; Er = 0.;
    J = isExc ? 0.55 : 0.95
    # J *= 0.8;
    τ_m = isExc ? 20 : 10
    g0 = (J*τ_m)/(τ_d - τ_r)
    g_ij = @fastmath calc_gij2!(t, m_tS_ext, v_numS_ext, i, τ_L, g0, τ_d, τ_r) 
    I_ij = g_ij
    # I_ij = g_ij*(V - Er)
    Isyn += I_ij;
    if (idxf == 1)
        m_Ia[i,aux] = Ia;
        m_Ig[i,aux] = Ig;
    end
    return Isyn
end

#--- idxf: how many times this function has been called for the same t (in 4th order RK4 it is called 4 times per t)
function LIF!(t::Float64, u::Array{Float64,1}, du::Array{Float64,1}, p::Array{Any,1}, aux::Int64, idxf::Int64)
    v_τ_m::Array{Float64, 1}, N::Int64, v_I::Array{Float64, 1}, m3_adj::Array{ Array{Array{Int64,1},1}, 1},  v_inh::Array{Int64, 1}, v_t_refrac::Array{Float64, 1}, tt::Array{Array{Float64, 1}, 1}, v_numSpikes::Array{Int64, 1},  m_tS_ext::Array{Array{Float64, 1}, 1}, v_numS_ext::Array{Int64, 1}, m_Ia::Array{Float64,2}, m_Ig::Array{Float64,2} = p 
    τ_L = 1.0 #ms
    Threads.@threads for i = 1:N
        if i ∈ v_inh idx = 2 else idx = 1 end #determine if neron is I (= 2) or E (= 1)
        t_refrac = v_t_refrac[idx]; τ_m = v_τ_m[idx]
        @inbounds if v_numSpikes[i] > 0 t0 = tt[i][v_numSpikes[i]]::Float64 else t0 = -(t_refrac+1.0) end #find time of last spike of ith neuron, used for refractory period

        if (t >= (t0 + t_refrac))
            @inbounds V, I = u[1*i], v_I[1*i]
            @inbounds @fastmath Isyn = calcIsyn!(i, N, V, m3_adj, t, tt, v_numSpikes, m_tS_ext, v_numS_ext, m_Ia, m_Ig, aux, idxf) #calculate synaptic current that arrives at neuron i
            @inbounds @fastmath du[i] = ( -V + Isyn  )/τ_m
        else
            du[i] = 0
        end
    end
    return du
end
