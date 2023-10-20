using PauliOperators
using UnitaryPruning
using ACSE
using Printf
using LinearAlgebra

bfs_thresh  = 1e-5
grad_thresh = 1e-8

## Number of orbitals and electrons 
norb = 2 ##Need to set this automatically
N_el = 2
## Transform Molecular Hamiltonian through JW Mapping

#H = ACSE.transform_molecular_Hamiltonian("H4_0.75")
H = ACSE.transform_contracted_Hamiltonian("H2_0.75", N_el)
#println("Hamiltonian", H)

#println("Exact Diagonalization for untransformed Hamiltonian: ", eigvals(Matrix(H)))
#exit()
## Define HF reference state
ket = ACSE.reference_state(norb)

## Reference Energy
reference_energy = ACSE.calc_energy(H, ket)
@printf("Reference energy: %10.8f\n", real(reference_energy))

## Generate Operator Pool
A = ACSE.acse_residual_pool_test(norb, norb)

## Find transformed Heisenberg Hamiltonian
H_transformed = ACSE.evolve_Hamiltonian(A, H, ket, bfs_thresh, grad_thresh)

## Final ACSE Energy
energy = ACSE.calc_energy(H_transformed, ket)
@printf("Energy: %10.8f\n", real(energy))


