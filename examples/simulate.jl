using PauliOperators
using UnitaryPruning
using ACSE
using Printf
using LinearAlgebra
using StatProfilerHTML

bfs_thresh  = 1e-8
grad_thresh = 1e-8

## Number of orbitals and electrons 
norb = 4 ##Need to set this automatically
N_el = 4
N = 2*norb
## Transform Molecular Hamiltonian through JW Mapping

H = ACSE.transform_molecular_Hamiltonian("H4_0.75")
#H = ACSE.transform_contracted_Hamiltonian("H2_0.75", N_el)
#println("Hamiltonian", H)

Hmat = Matrix(H)
println("Exact Diagonalization for untransformed Hamiltonian: ") 
e_exact = eigvals(Hmat)
for i in 1:10
    @printf(" Exact Energies %2i %12.8f\n",i, e_exact[i])
end

# F = eigen(Hmat)
#println(F.values[1])
#println.(F.vectors[:,1])
#exit()
## Define HF reference state
ket = PauliOperators.KetBitString(N,0)
for i in 1:N_el
    global coeff, ket = Pauli(N,X=[i]) * ket
end
@show ket

## Reference Energy
reference_energy = ACSE.calc_energy(H, ket)
@printf("Reference energy: %10.8f\n", real(reference_energy))

## Generate Operator Pool
A = ACSE.acse_residual_pool_test(norb, norb)

## Find transformed Heisenberg Hamiltonian
H_transformed, energies, gradients, g_seq = ACSE.evolve_Hamiltonian(A, H, ket, bfs_thresh, grad_thresh, max_iter=50, verbose=1, alpha=.1)





# new_gen = []
# push!(new_gen, g_seq[1])
# for gi in 2:length(g_seq)
#     if new_gen[end][1] == g_seq[gi][1]
#         new_gen[end] = (new_gen[end][1], new_gen[end][2]+g_seq[gi][2])
#     else
#         push!(new_gen, g_seq[gi])
#     end
# end
# @show new_gen
# generators = [A[i]*j for (i,j) in new_gen]
# H_transformed = deepcopy(H)
# for gi in 1:length(generators)
#     g = generators[gi]
#     global H_transformed = ACSE.BFS(g, H_transformed, ket, thresh=bfs_thresh/10)
#     energyi = ACSE.calc_energy(H_transformed, ket)
#     @printf("Iter: %4i Energy: %10.8f + %10.8fi Gradient: %6.1e #Ops in H: %6i Op IDX: %4i\n", gi, real(energyi),imag(energyi), 0, length(H_transformed), 0)
# end
# @show new_energy = ACSE.calc_energy(H_transformed, ket)

# # #println(" Energy/gradient by iteration")
# # #for i in 1:length(energies)
# # #    @printf(" %4i %12.8f %12.8f\n", i, energies[i], gradients[i])
# # #end

# # ## Final ACSE Energy
# # energy = ACSE.calc_energy(H_transformed, ket)
# # @printf("Energy: %10.8f\n", real(energy))


