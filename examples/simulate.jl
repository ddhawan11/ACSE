using PauliOperators
using UnitaryPruning
using ACSE
using Printf
using LinearAlgebra
using StatProfilerHTML

ED = false

bfs_thresh     = 2e-5
grad_thresh    = 1e-8
energy_thresh  = 1e-8

##Need to set this automatically
norb = 6    # No. of spatial orbitals
N_el = 6    # No. of electrons
N = 2*norb  # No. of spin orbitals

## Transform Molecular Hamiltonian through JW Mapping
H = ACSE.transform_molecular_Hamiltonian("H6_0.75")
#H = ACSE.transform_contracted_Hamiltonian("H6_0.75", N_el)


if ED
    Hmat = Matrix(H)
    println("Exact Diagonalization for untransformed Hamiltonian: ") 
    e_exact = eigvals(Hmat)
    for i in 1:10
        @printf(" Exact Energies %2i %12.8f\n",i, e_exact[i])
    end
end

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
A = ACSE.qubit_pool(N)
#A = ACSE.acse_residual_pool(N)


## Find transformed Heisenberg Hamiltonian
H_transformed, energies, gradients, g_seq = ACSE.evolve_Hamiltonian(A, H, ket, bfs_thresh, grad_thresh, energy_thresh, max_iter=100, verbose=1, alpha=0.1)

ACSE.extrapolate_normerror(H, H_transformed)


