using PauliOperators
using ACSE
using UnitaryPruning

function commutator(A::PauliSum{N}, B::PauliSum{N}, ket) where N
    exp_val = zero(ComplexF64) 
    commutator = A*B - B*A

    for (key,value) in commutator.ops
        exp_val += value * PauliOperators.expectation_value(key, ket)
    end
    return exp_val
end


