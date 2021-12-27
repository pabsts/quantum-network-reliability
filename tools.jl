
function R_init(f::Real) 
    θ = 2*acos(sqrt(f))
    return Ry(θ)
end


function buildC_prepare_edge(n)    
     offset = 1
    
    return chain(n, 
        put(i+offset => R_init(failure[i])) for i=1:N_edge
    )
end


function buildC_turn_on_node(n, node_id=1)
    qbit_offset = 1 + N_edge

    return put(n, qbit_offset + node_id => X)
end


function buildC_propagate_node(n)

    function buildC_turn_on_next_node(node_source, node_target, edge)
        offset_edge = 1
        offset_node = offset_edge + N_edge
        offset_node_aux = offset_node + N_node
        c1 = cnot(offset_node + node_target, offset_node_aux + node_target)

        c1_1 = put(offset_node_aux + node_target => X)


        c2 = cnot(
            (offset_node + node_source, offset_node_aux + node_target, offset_edge + edge),
            offset_node + node_target
        )

        c3 = put(offset_node_aux + node_target => H)
        m = Yao.Measure(; locs=offset_node_aux + node_target, resetto=bit"0")

        return chain(c1, c1_1, c2, c3, m)
    end

    # circuit for propaget node over all edges
    circuit = chain(buildC_turn_on_next_node(n1, n2, e) for (e, (n1,n2)) in enumerate(edge))
    circuit = circuit(n)

    # repeat this propagtion N_node-1 times
    circuit = repeat([circuit], N_node-1)
    circuit = chain(circuit)

    return circuit
end

function buildC_connected(n)
    offset_node = 1 + N_edge

    return cnot(n,
        offset_node .+ (1:N_node),
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
        # bs = [bs_rest, bs_edge, bs_label]
        # println("""$i: $(join(bs, "_"))""")

        if isnothing(ψ)
            ψ = q * YaoSym.ket_m(bs_interest)
        else
            ψ += q * YaoSym.ket_m(bs_interest)
        end
    end

    return ψ
end

