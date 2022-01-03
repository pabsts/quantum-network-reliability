using Yao #, YaoBlocks
using YaoSym, SymEngine
using YaoPlots, YaoExtensions

import Cairo  # needed for plotting
using Graphs, SimpleWeightedGraphs, GraphPlot
#using Compose

# using Graphs, SimpleWeightedGraphs
# import SparseArrays

include("tools.jl")


failure_default = 1-0.999

# graphs
# -------
graph1 = SimpleWeightedGraph(3)
let g = graph1
    edges = [(1,2), (2,3), (1,3)]
    add_edges!(g, edges, failure_default)
end


graph2 = SimpleWeightedGraph(10)
let g = graph2
    edges = [
        (1,2), (1,3), (2,4), (3, 4), (1,4),
        (2, 5), (3, 6),
        (7, 8), (7, 9), (8, 10), (9, 10), (7, 10),
        (5, 8), (6, 9)
    ]

    add_edges!(g, edges, failure_default)
end
gplot(graph2)

graph21 = SimpleWeightedGraph(4)
let g = graph21
    edges = [(1,2), (1,3), (2,3), (2,4), (3, 4)]
    add_edges!(g, edges, failure_default)
end
gplot(graph21)

graph22 = SimpleWeightedGraph(4)
let g = graph22
    edges = [(1,2), (1,3), (2,3), (2,4), (3, 4)]
    add_edges!(g, edges, failure_default)
end
gplot(graph22)


graph = graph22

# initial quantum state
# ---------------------
phi = get_initial_state(graph)
# N_qbit = length(phi.state) |> log2 |> Int

circuit = get_quantum_circuit(graph; steps=1:4);
# plot(circuit)
@time phi_final = apply(phi, circuit)
# analzye(phi_final)

# last step calculate reduced density matrix of label qubit
@time rdm = (focus!(phi_final, [1]) |> density_matrix |> state)[:,:, 1]
display(rdm)
@info "network reliability is $(rdm[2,2])"
