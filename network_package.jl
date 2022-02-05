using Printf
using DelimitedFiles


function read_connection_parameters(Ni, Ne, conn_prob, seed=0, cond_version="2")
    # m_params_syn, header = readdlm("data/synaptic_parameters.dat", ' ', header=true) #matrix with receptor parameters
    N = Ni + Ne

    #read adjancecy matrices; in old convention (starting at 0)
    v_nomes_receptores = ["ampa", "gaba"]
    m3_adj = convert(Array{Array{Array{Int64, 1}, 1}, 1}, []) #[1] = ampa, [2] = nmda, [3] = gaba
    # m3_w = convert(Array{Array{Array{Float64, 1}, 1}, 1}, []) #[1] = ampa, [2] = nmda, [3] = gaba
    for idx_syn =1:2
        receptor = v_nomes_receptores[idx_syn]        
        fileIn = "data/connection_data/matriz_conexoes_" * receptor * "_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_conn_prob_" * string(conn_prob) * "_seed_" * string(seed) *"_conductances_" * string(cond_version) * ".dat"
        aux = readdlm(fileIn)
        aux = [filter!(x->x≠"", aux[i,:]) for i=1:N]
        aux = convert(Array{Array{Int64, 1},1}, aux)
        push!(m3_adj, aux)
 
        #  fileIn = "data/connection_data/matriz_pesos_" * receptor * "_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_conn_prob_" * string(conn_prob) * "_seed_" * string(seed) *"_conductances_" * string(cond_version) * ".dat"
        # aux = readdlm(fileIn)
        # aux = [filter!(x->x≠"", aux[i,:]) for i=1:N]
        # aux = convert(Array{Array{Float64, 1},1}, aux)
        # push!(m3_w, aux)
    end

    #read inh and exc neurons (also in C convention) 
    m = readdlm("data/connection_data/vetores_inh_exc_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_conn_prob_" * string(conn_prob) * "_seed_" * string(seed) *"_conductances_" * string(cond_version) * ".dat")
    # println(m)
    v_inh = filter!(x->x≠"", m[1,:])
    v_exc = filter!(x->x≠"", m[2,:])

    v_inh = convert(Array{Int64, 1},v_inh)
    v_exc = convert(Array{Int64, 1},v_exc)

    # v_inh .+= 1; v_exc .+= 1 #convert to Julia 


    return  m3_adj, v_inh, v_exc
end
# function read_connection_parameters(Ni, Ne, conn_prob, seed=0, cond_version="2", fileIn="")
#     # m_params_syn, header = readdlm("data/synaptic_parameters.dat", ' ', header=true) #matrix with receptor parameters
#     N = Ni + Ne

#     #read adjancecy matrices; in old convention (starting at 0)
#     v_nomes_receptores = ["ampa", "gaba"]
#     m3_adj = convert(Array{Array{Array{Int64, 1}, 1}, 1}, []) #[1] = ampa, [2] = nmda, [3] = gaba
#     # m3_w = convert(Array{Array{Array{Float64, 1}, 1}, 1}, []) #[1] = ampa, [2] = nmda, [3] = gaba
#     for idx_syn =1:2
#         receptor = v_nomes_receptores[idx_syn]        
#         iffileIn == "" fileIn = "data/connection_data/matriz_conexoes_" * receptor * "_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_conn_prob_" * string(conn_prob) * "_seed_" * string(seed) *"_conductances_" * string(cond_version) * ".dat" end 
#         aux = readdlm(fileIn)
#         aux = [filter!(x->x≠"", aux[i,:]) for i=1:N]
#         aux = convert(Array{Array{Int64, 1},1}, aux)
#         push!(m3_adj, aux)
 
#         #  fileIn = "data/connection_data/matriz_pesos_" * receptor * "_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_conn_prob_" * string(conn_prob) * "_seed_" * string(seed) *"_conductances_" * string(cond_version) * ".dat"
#         # aux = readdlm(fileIn)
#         # aux = [filter!(x->x≠"", aux[i,:]) for i=1:N]
#         # aux = convert(Array{Array{Float64, 1},1}, aux)
#         # push!(m3_w, aux)
#     end

#     #read inh and exc neurons (also in C convention) 
#     m = readdlm("data/connection_data/vetores_inh_exc_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_conn_prob_" * string(conn_prob) * "_seed_" * string(seed) *"_conductances_" * string(cond_version) * ".dat")
#     # println(m)
#     v_inh = filter!(x->x≠"", m[1,:])
#     v_exc = filter!(x->x≠"", m[2,:])

#     v_inh = convert(Array{Int64, 1},v_inh)
#     v_exc = convert(Array{Int64, 1},v_exc)

#     # v_inh .+= 1; v_exc .+= 1 #convert to Julia 


#     return  m3_adj, v_inh, v_exc
# end

function saveSpikeTimes(fileOut, tt)
    open(fileOut, "w") do io 
        for line in tt 
            for item in line 
                @printf(io, "%f ", item)
            end
            if (length(line) == 0)
                @printf(io, "0")
            end
            @printf(io, "\n")
        end
    end
end

# function saveLFP(fileOut, LFP, v_t)
#     open(fileOut, "w") do io 
#         for i=1:length(LFP)
#             @printf(io, "%f %f\n", LFP[i], v_t[i])
#         end
#     end
# end


function connectedToWhom(i, m_m_adj, N)
    v_j = []
    for m_adj in m_m_adj 
        for (k, g0) in m_adj
            j = Int64(k - (i-1)*N + 1)::Int64 #select neuron j, which is connected to i; subtrai e soma o 1 para transformar da convençao de indice começando em 0 (que é o que uso no v_adj) para começando em 1
            if (j > N || j < 1) #checa se j conectado a esse i; 2a condiçao vale pra qnd acabar a matriz de conexoes, mas nao acabarem os i 
                continue end 
            println(j)
            push!(v_j, j)
        end
    end
    return v_j
end

