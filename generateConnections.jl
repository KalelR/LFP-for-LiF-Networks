using Random
using Distributions
using DelimitedFiles


#auxiliar function to choose a number from v_aux, used in generateConnections function
function chooseConn(v_aux)
    while(true)
        conn = rand(v_aux)
        if (conn ≠ -1)
            return conn
        end
    end
end

#function to calculate everything connection-related: picks inhibitory/excitatory neurons; picks the connections that will be made; picks which to which receptor each connection will be
function generateConnections(Ni, Ne, conn_p, m_params_recep, seed = 0)
    N = Ni + Ne

    #determine E and I
    Random.seed!(seed)
    v_exc = collect(0:(N-1))
    v_inh = sort([splice!(v_exc, rand(MersenneTwister(seed), eachindex(v_exc))) for i=1:Ni ]) #chooses random numbers from v_exc (which initially contains all neurons) and puts that number into v_inh; does that Ni times; the eachindex(v_exc) is the same as v_exc in this case

    println("----concluded step 1")

    #determine connections to be made (NOTE: indices start at 0 in this convention)
    Random.seed!(seed)
    Nconn = Int64(floor(conn_p*N*(N-1)))
    v_aux = collect( 0:(N*N-1) )
    for i = 0:(N-1)
        v_aux[i*(N+1)+1] = -1;
        # if (i % 100 == 0)  println(i) end
    end

    v_conn = []
    for i=1:Nconn
        conn = chooseConn(v_aux)
        push!(v_conn, conn);
        conn_idx = conn + 1
        v_aux[conn_idx] = -1;
        # println(conn, ' ', length(v_conn), ' ', length(v_aux), ' ' , (length(v_conn) + length(v_aux)) )
    end
    # filter!(x->x≠(-1), v_aux) #filter will be useful if p_conn increases
    v_conn = sort(v_conn);

    println("----concluded step 2")

    #assign receptors to each connections and create adjacency matrices for each receptor
    Random.seed!(seed)
    m3_adj = [ [ [-1] for j=1:N ]  for i=1:2]
    # m3_w =  [ [ [-1.] for j=1:N ]  for i=1:3]
    for conn in v_conn
        i, j = floor(Int, conn/N), floor(Int, conn % N) #i = postynatpc,  j = presynaptic

        #now to find g0 of the connection
        if j ∈ v_inh #if presynaptic is I
            idx_syn = 2 #GABA
            # v_g0_g = m_params_recep[idx_syn]
            # if i ∈ v_inh g0 = v_g0_g[1] else g0 = v_g0_g[2] end #pick which g0 (depends whether presynaptic is E([2]) or I([1]))

        else #if presynaptic is E
            # idx_syn = rand(1:2) #pick if AMPA or NDMA (ampa = 1, ndma = 2)
            idx_syn = 1
            # if i ∈ v_inh
                # g0 = 1.0 #nS
            # else
                # μ, σ = m_params_recep[idx_syn]
                # g0 = rand(Distributions.LogNormal(μ, σ))
            # end
        end

        #now save it to the appropriate matrix
        # println(i, conn)
        if (m3_adj[idx_syn][i+1][1] == -1)
             m3_adj[idx_syn][i+1][1] = j+1
        else
             push!(m3_adj[idx_syn][i+1], j+1)
        end

        # if (m3_w[idx_syn][i+1][1] == -1)
        #      m3_w[idx_syn][i+1][1] = g0;
        # else
        #      push!(m3_w[idx_syn][i+1],  g0)
        # end

    end

    v_inh .+= 1
    v_exc .+= 1

    # return v_inh, v_exc, m3_adj, m3_w
    return v_inh, v_exc, m3_adj
end



Ni = 1000
Ne = 4000
conn_prob = 0.2 #if this gets much higher, use the filter! in the generateConnections function
seed = 0
N = Ni + Ne
μ_a, σ_a = 10^(-0.31), 10^(-0.3)
μ_n, σ_n = μ_a, σ_a
m_params_recep = [[μ_a, σ_a], [μ_n, σ_n], [0, 0]]  # [3] contém g0 = 5 para inibitorio ([1]) g0 = 4 para excitatorio ([2]) no GABA; [1] e [2] contem μ e σ para o lognormal dos dois receptores
# m_params_recep = [[μ_a, σ_a], [μ_n, σ_n], [5, 4]]  # [3] contém g0 = 5 para inibitorio ([1]) g0 = 4 para excitatorio ([2]) no GABA; [1] e [2] contem μ e σ para o lognormal dos dois receptores
conductanes_version = "1"
v_nomes_receptores = ["ampa", "gaba"]


v_inh, v_exc, m3_adj = generateConnections(Ni, Ne, conn_prob, m_params_recep, seed)

fileOut = "data/connection_data/vetores_inh_exc_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_conn_prob_" * string(conn_prob) * "_seed_" * string(seed) *"_conductances_" * string(conductanes_version) * ".dat"
open(fileOut, "w") do io
    writedlm(io, [v_inh, v_exc], ' ')
end


# writedlm(io, m_m_adj, ' ')
for idx_syn =1:length(m3_adj)
    receptor = v_nomes_receptores[idx_syn]
    fileOut = "data/connection_data/matriz_conexoes_" * receptor * "_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_conn_prob_" * string(conn_prob) * "_seed_" * string(seed) *"_conductances_" * string(conductanes_version) *".dat"
    open(fileOut, "w") do io
        writedlm(io, m3_adj[idx_syn], " ")
    end
    # fileOut = "data/connection_data/matriz_pesos_" * receptor * "_Ni_" * string(Ni) * "_Ne_" * string(Ne) * "_conn_prob_" * string(conn_prob) * "_seed_" * string(seed) *"_conductances_" * string(conductanes_version) *".dat"
    # open(fileOut, "w") do io
    #     writedlm(io, m3_w[idx_syn], " ")
    # end

end
