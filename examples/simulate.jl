using PauliOperators
using UnitaryPruning
using ACSE

norb = 6
H = ACSE.transform_molecular_Hamiltonian("H6_0.75")

ket = ACSE.reference_state(norb)

reference_energy = ACSE.calc_energy(H, ket)
println("Reference energy:", real(reference_energy))

A = ACSE.acse_residual_pool(norb, norb)

generator = ACSE.find_generator(A, H, ket)
println(generator)

