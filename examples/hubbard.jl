
using PauliOperators
using UnitaryPruning
using ACSE
using Printf
using LinearAlgebra

bfs_thresh  = 1e-5
grad_thresh = 1e-8

graph = zeros(Int,4,4)
graph[1,2] = 1
graph[2,3] = 1
graph[3,4] = 1
graph[1,4] = 1
H = ACSE.hubbard_hamiltonian(graph, 1, 1)
A = ACSE.hubbard_pool(graph)