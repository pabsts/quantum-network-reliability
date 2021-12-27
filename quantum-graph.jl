using Yao #, YaoBlocks
using YaoSym, SymEngine
using YaoPlots, YaoExtensions

import Cairo
#using Compose

# using Graphs, SimpleWeightedGraphs
# import SparseArrays

include("tools.jl")

# graph
# -------
N_node  = 3
node = 1:N_node
edge = [(1,2), (2,3)]
N_edge = length(edge)

failure_default = 0.2
failure = repeat([failure_default], N_edge)


# initial quantum state
# ---------------------
N_batch = 1
phi_edge = zero_state(N_edge; nbatch=N_batch)
phi_node = zero_state(N_node; nbatch=N_batch)
phi_node_aux = zero_state(N_node; nbatch=N_batch)
phi_label = zero_state(1; nbatch=N_batch)

phi = join(phi_node_aux, phi_node, phi_edge, phi_label)
N_qbit = length(phi.state) |> log2 |> Int


# step 1: build circuit
circuit = []
c = buildC_prepare_edge(N_qbit)
push!(circuit, c)

# step 2: turn on one node
c = buildC_turn_on_node(N_qbit)
push!(circuit, c)

# step 3: propagate turned-on node across graph
c = buildC_propagate_node(N_qbit)
push!(circuit, c)

# step 4: assign label qubit
c = buildC_connected(N_qbit)
push!(circuit, c)

# # step 5: measure label qubit
# c = Yao.Measure(N_qbit; locs=1)
# push!(circuit, c)


# analyze each step
# -----------------
# step 0: visualize state
analzye(phi)

# apply circuit up to step `istep`
istep = 4
c = circuit[1:istep] |> chain
plot(c)
phi_final = apply(phi, c)
analzye(phi_final)

# check `phi_final` is a pure-state
let ρ = density_matrix(phi_final).state[:,:,1]
    @show trnorm(ρ * ρ)
end


# last step calculate reduced density matrix of label qubit
rdm = (focus!(phi_final, [1]) |> density_matrix |> state)[:,:, 1]
if sum(abs.(imag.(rdm)))>0
    @warn "reduced density matrix has has imaginary coefficients"   
else
    rdm = real.(rdm)
end

@info "network reliability is $(rdm[2,2])"
display(rdm)