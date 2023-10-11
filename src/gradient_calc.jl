using PauliOperators
using ACSE
using UnitaryPruning

## \langle 0| [H, A]|0 \rangle
function compute_grad(A::PauliSum{N}, B::PauliSum{N}, ket) where N
    exp_val = zero(ComplexF64) 
    commutator = A*B - B*A

    for (key,value) in commutator.ops
        exp_val += value * PauliOperators.expectation_value(key, ket)
    end

    return exp_val
end


