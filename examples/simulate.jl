using PauliOperators
using UnitaryPruning
using ACSE
using Printf

bfs_thresh  = 1e-6
grad_thresh = 1e-8

norb = 2

## Transform Molecular Hamiltonian through JW Mapping

H = ACSE.transform_molecular_Hamiltonian("H2_0.75")
#println("Hamiltonian", H)

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


