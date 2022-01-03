function add_edges!(graph, edges::Vector, failure_default)
    for x in edges        
        failure = if length(x)==3 
            x[3]
        else
            failure_default
        end

        add_edge!(graph, x[1], x[2], failure)
    end
    return graph
end

function get_initial_state(graph)
    N_batch = 1

    N_edge = Graphs.ne(graph)
    N_node = Graphs.nv(graph)
    
    phi_label = zero_state(1; nbatch=N_batch)
    phi_edge = zero_state(N_edge; nbatch=N_batch)
    phi_node = zero_state(N_node; nbatch=N_batch)
    phi_aux = zero_state(1; nbatch=N_batch)
    
    return join(phi_aux, phi_node, phi_edge, phi_label)
end

function required_qbits(graph)
    N_edge = Graphs.ne(graph)
    N_node = Graphs.nv(graph)
    return 1 + N_edge + N_node + 1
end


function get_quantum_circuit(graph::AbstractSimpleWeightedGraph; steps=1:5)

    circuit = []

    for s in 1:5
        s ∉ steps && continue

        group = if s==1 # step 1: build circuit
            buildC_prepare_edge(graph)
        elseif s==2 # step 2: turn on one node
            buildC_turn_on_node(graph)
        elseif s==3 # step 3: propagate turned-on node across graph
            buildC_propagate_node(graph)
        elseif s==4 # step 4: assign label qubit
            buildC_connected(graph)
        elseif s==5 # step 5: measure label qubit
            N_qbits = required_qbits(graph)
            Yao.Measure(N_qbit; locs=1)       
        else
            continue
        end

        push!(circuit, group)
    end

    # join steps
    circuit = chain(circuit)

    return circuit
end


function R_init(failure::Real) 
    θ = 2*acos(sqrt(failure))
    return Ry(θ)
end


function buildC_prepare_edge(graph)
    n = required_qbits(graph)
    qbit_offset = 1

    return chain(n,
        put(i+qbit_offset => R_init(e.weight)) for (i,e) in enumerate(edges(graph))
    )
end


function buildC_turn_on_node(graph; node=1)
    n = required_qbits(graph)
    N_edge = Graphs.ne(graph)
    qbit_offset = 1 + N_edge

    return put(n, qbit_offset + node => X)
end


function buildC_propagate_node(graph)
    n = required_qbits(graph)
    N_node = Graphs.nv(graph)
    N_edge = Graphs.ne(graph)

    function buildC_turn_on_next_node(node_source, node_target, edge)
        offset_edge = 1
        offset_node = offset_edge + N_edge
        node_aux = offset_node + N_node + 1

        c1 = cnot(offset_node + node_target, node_aux)
        c2 = put(node_aux => X)

        c3 = cnot(
            (offset_node + node_source, node_aux, offset_edge + edge),
            offset_node + node_target
        )

        c4 = put(node_aux => H)
        m = Yao.Measure(; locs=node_aux, resetto=bit"0")

        return chain(c1, c2, c3, c4, m)
    end

    # circuit for propaget node over all edges
    elements = []
    for (i, e) in enumerate(edges(graph))
        circuit1 = buildC_turn_on_next_node(e.src, e.dst, i)
        circuit2 = buildC_turn_on_next_node(e.dst, e.src, i)
        circuit12 = chain(n, circuit1, circuit2)
        push!(elements, circuit12)
    end
    circuit = chain(elements)

    # repeat this propagtion N_node-1 times
    circuit = repeat([circuit], N_node-1)
    circuit = chain(circuit)

    return circuit
end


function buildC_connected(graph)
    n = required_qbits(graph)
    N_node = Graphs.nv(graph)
    N_edge = Graphs.ne(graph)
    qbit_offset = 1 + N_edge

    return cnot(n,
        qbit_offset .+ (1:N_node),
        1
    )
end


# Analyzing state
# ---------------

analzye(s::ArrayReg) = analzye(state(s)[:, 1])

function analzye(s::Vector)
    
    ψ = nothing

    # visualize only the label and edge qubits
    for (i, q) in enumerate(s)
        abs(q) == 0 && continue
        
        if imag(q) != 0
            @warn "state $i has an imaginary coefficient"
        else
            q = real(q)
        end

        # find bit representation
        # bs = bitstring(i-1)#[end-N_qbit+1:end]
        qbits = digits(i-1, base=2, pad=N_qbit)
        bs = join(qbits |> reverse)

        bs_label = bs[end:end]
        bs_edge = bs[end-N_edge:end-1]
        bs_rest = bs[1:end-N_edge-1]
        bs_interest = bs_edge * bs_label

        # # info 
        bs = [bs_rest, bs_edge, bs_label]
        println("""$i: $(join(bs, "_"))""")

        if isnothing(ψ)
            ψ = q * YaoSym.ket_m(bs_interest)
        else
            ψ += q * YaoSym.ket_m(bs_interest)
        end
    end

    return ψ
end



"""
Implement naive prute-force approach
"""
function network_reliability(g::Graphs.Graph, failure; level=0)::Real
    
    # if level > 2
    #     return 0.0
    # end

    ncomp = Graphs.connected_components(g)
    if length(ncomp)>1
        return 0.0
    end
    
    # probability of current subgraph

    r = (1.0-failure) ^ Graphs.ne(g)   
    ne = Graphs.ne(g)
    # @show r, ne

    # remove an edge and check if it is still fully connected
    for e in edges(g)
        g1 = copy(g)
        rem_edge!(g1, e)
        

        ncomp = Graphs.connected_components(g1)
        if length(ncomp)>1 # not connected
            continue

        else # subgraph is still fully connected

            r_sub = network_reliability(g1, failure; level=level+1) 
            # add reliability of subgraph (weighted by failure of current removed edge)            
            r += r_sub * failure / (level+1)  # divide by combinatorial factor

            # if r_sub>0
            #     @show ne, e, r, r_sub
            # end
        end
    end

    return r
end 