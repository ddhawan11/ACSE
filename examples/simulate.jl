using PauliOperators
using UnitaryPruning
using ACSE



H = ACSE.transform_molecular_Hamiltonian("data_h6")
A = ACSE.acse_residual_pool(6,6)

println(length(A))
exit()
ket = KetBitString(12, 31)
for op in A
    ACSE.commutator(H,op,ket)
end
