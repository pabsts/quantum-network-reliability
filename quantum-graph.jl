using Yao, YaoSym, SymEngine

using YaoPlots, YaoExtensions
import Cairo
using Compose

using Graphs, SimpleWeightedGraphs
import SparseArrays
# start with a graph


N_node  = 3
failure = 0.2

g = SimpleWeightedGraph(N_node)  # or use `SimpleWeightedDiGraph` for directed graphs
add_edge!(g, 1, 2, failure)
add_edge!(g, 2, 3, failure)
# add_edge!(g, 1, 3, 2.0)

N_edge = g.weights |> SparseArrays.nnz


# initial quantum state
phi_edge = szero_state(N_edge)
phi_node = szero_state(N_node)
phi_node_aux = szero_state(N_node)
phi_label = szero_state(1)


# build circuit
function R_init(f::Real) 
    θ = 2*asin(sqrt(f))
    return Ry(θ)
end

circuit = repeat(N_edge, R_init(failure))
plot(circuit)


s = apply(ket"0", R_init(failure))

(bra"0" * apply(ket"0", R_init(failure))) .^2


plot(H)


apply(phi_edge, circuit)
plot(qft_circuit(5))