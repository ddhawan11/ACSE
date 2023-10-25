
using PauliOperators
using UnitaryPruning
using ACSE
using Printf
using LinearAlgebra
using BlockDavidson

bfs_thresh  = 1e-4
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
A2 = ACSE.acse_residual_pool_test(n_orbs, n_orbs)
println(" Length of Pool:", length(A))
# A = ACSE.zy_pool(N)
# t2 = jordan_wigner(7,N) * jordan_wigner(8,N) * jordan_wigner(9,N)' * jordan_wigner(10,N)'
# push!(A, t2-adjoint(t2))
# A = ACSE.acse_residual_pool_test(n_orbs, n_orbs)
#A,Amix = ACSE.hubbard_pool2(graph) 
Na,Nb = ACSE.hubbard_number(graph)


# ket = PauliOperators.KetBitString([1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1])
# ket = PauliOperators.KetBitString([1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0])
ket = KetBitString([i<=n_elec for i in 1:N])

H   = ACSE.particle_hole_transformation(H, ket)
A   = ACSE.particle_hole_transformation(A, ket)
Na  = ACSE.particle_hole_transformation(Na, ket)
Nb  = ACSE.particle_hole_transformation(Nb, ket)
A2  = ACSE.particle_hole_transformation(A2, ket)
ket = ACSE.particle_hole_transformation(ket, ket)

@show ket

@show ACSE.calc_energy(Na+Nb,ket)
@show ACSE.calc_energy(Na-Nb,ket)
## Reference Energy
reference_energy = ACSE.calc_energy(H, ket)
@printf("Reference energy: %10.8f\n", real(reference_energy))


H_transformed, energies, gradients, g_seq = ACSE.evolve_Hamiltonian(A, H, ket, bfs_thresh, grad_thresh, max_iter=100, verbose=1, alpha=.1)
Htmp = deepcopy(H_transformed)
H_transformed, energies, gradients, g_seq = ACSE.evolve_Hamiltonian(A2, H_transformed, ket, bfs_thresh, grad_thresh, max_iter=20, verbose=1, alpha=.1)

Htmp, energies, = ACSE.evolve_Hamiltonian(A2, g_seq, Htmp, ket, bfs_thresh, verbose=1)


# term_sizes = []
# for (pauli, coeff) in Htmp.ops
#     tmp, _, = ACSE.evolve_Hamiltonian(A2, g_seq, PauliSum(Pauli(pauli))*coeff, ket, bfs_thresh, verbose=0)
#     #println(length(tmp))
#     push!(term_sizes, length(tmp))
# end

# println(maximum(term_sizes))

# for i in 1:10
#     ## Find transformed Hamiltonian
#     global H_transformed, energies, gradients, g_seq = ACSE.evolve_Hamiltonian(A, H_transformed, ket, bfs_thresh, grad_thresh, max_iter=200, verbose=1, alpha=.01)
        
#     ## Find transformed Hamiltonian
#     global H_transformed, energies, gradients, g_seq = ACSE.evolve_Hamiltonian(A2, H_transformed, ket, bfs_thresh, grad_thresh, max_iter=6, verbose=1, alpha=.1)
# end

# M = 6
# scr = zeros(ComplexF64, 2^N, M)
# function mymatvec(v)
#    fill!(scr,0.0) 
#    mul!(scr, H, v)
#    return scr 
# end
# lmat = LinOpMat{ComplexF64}(mymatvec, 2^N, true)
# v0 = zeros(ComplexF64, 2^N,1)
# v0[ket.v+1] += 1
# dav = Davidson(lmat, T=ComplexF64, nroots=M)
# @time e2, v2 = eigs(dav)
