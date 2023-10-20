using PauliOperators
using UnitaryPruning
using ACSE
using Printf
using LinearAlgebra

bfs_thresh  = 1e-6
grad_thresh = 1e-8

## Number of orbitals and electrons 
norb = 4 ##Need to set this automatically
N_el = 4
N = 2*N_el
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

F = eigen(Hmat)
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
H_transformed, energies, gradients = ACSE.evolve_Hamiltonian(A, H, ket, bfs_thresh, grad_thresh, max_iter=40, verbose=1, alpha=.1)

#println(" Energy/gradient by iteration")
#for i in 1:length(energies)
#    @printf(" %4i %12.8f %12.8f\n", i, energies[i], gradients[i])
#end

## Final ACSE Energy
energy = ACSE.calc_energy(H_transformed, ket)
@printf("Energy: %10.8f\n", real(energy))


