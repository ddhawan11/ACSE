using PauliOperators
using UnitaryPruning
using ACSE
using Printf

bfs_thresh  = 1e-6
grad_thresh = 1e-8

norb = 2
H = ACSE.transform_molecular_Hamiltonian("H2_0.75")
println("Hamiltonian", H)
ket = ACSE.reference_state(norb)

reference_energy = ACSE.calc_energy(H, ket)
@printf("Reference energy: %10.8f\n", real(reference_energy))

A = ACSE.acse_residual_pool(norb, norb)

H_transformed = ACSE.evolve_Hamiltonian(A, H, ket, bfs_thresh, grad_thresh)

energy = ACSE.calc_energy(H_transformed, ket)
@printf("Energy: %10.8f\n", real(energy))


