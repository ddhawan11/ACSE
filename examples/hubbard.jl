
using PauliOperators
using UnitaryPruning
using ACSE
using Printf
using LinearAlgebra
using BlockDavidson

bfs_thresh  = 1e-8
grad_thresh = 1e-8

n_orbs = 8 
n_elec = 8 
N = 2*n_orbs

graph = zeros(Int, n_orbs, n_elec)
for i in 1:n_orbs
    graph[i, (i+1)%n_orbs+1] = 1
end

H = ACSE.hubbard_hamiltonian(graph, 1, 1)
A = ACSE.hubbard_pool(graph)
# A = ACSE.acse_residual_pool_test(n_orbs, n_orbs)
Na,Nb = ACSE.hubbard_number(graph)



ket = PauliOperators.KetBitString(N, 16890)
# ket = PauliOperators.KetBitString([1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1])
# ket = PauliOperators.KetBitString([1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0])
@show ket

@show ACSE.calc_energy(Na+Nb,ket)
@show ACSE.calc_energy(Na-Nb,ket)
## Reference Energy
reference_energy = ACSE.calc_energy(H, ket)
@printf("Reference energy: %10.8f\n", real(reference_energy))


## Find transformed Hamiltonian
H_transformed, energies, gradients, g_seq = ACSE.evolve_Hamiltonian(A, H, ket, bfs_thresh, grad_thresh, max_iter=500, verbose=1, alpha=.1)

M = 1
scr = zeros(ComplexF64, 2^N, M)
function mymatvec(v)
    fill!(scr,0.0) 
    mul!(scr, H_transformed, v)
    return scr 
end
lmat = LinOpMat{ComplexF64}(mymatvec, 2^N, true)
v0 = zeros(ComplexF64, 2^N,1)
v0[ket.v+1] += 1
dav = Davidson(lmat, T=ComplexF64, nroots=M)
@time e2, v2 = eigs(dav)