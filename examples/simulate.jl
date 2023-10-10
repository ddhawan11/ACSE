using PauliOperators
using UnitaryPruning
using ACSE
using Printf


norb = 6
H = ACSE.transform_molecular_Hamiltonian("H6_0.75")

ket = ACSE.reference_state(norb)

reference_energy = ACSE.calc_energy(H, ket)
@printf("Reference energy: %10.8f\n", real(reference_energy))

A = ACSE.acse_residual_pool(norb, norb)

generator, curr_grad = ACSE.find_generator(A, H, ket)
@printf("Number of Paulis in Hamiltonian operator: %i\n", length(H))

H_transformed = ACSE.BFS(generator, H, ket, curr_grad, thresh=1e-4)

