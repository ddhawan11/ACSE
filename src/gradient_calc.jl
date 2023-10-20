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


"""
    compute_grad2(A::PauliSum{N}, B::PauliSum{N}, ket) where N

Compute gradient of expectation value of `A` w.r.t. a rotation about a series of Pauli axes, `B` evaluated at zero
all the operators in `B` should commute, otherwise its not well defined.

d/dθ < exp(-θB) A exp(θB) >  = <[A,B]> = <AB> - <BA> = <AB> + adj(<AB>) = 2R(<AB>)
"""
function compute_grad2(A::PauliSum{N}, B::PauliSum{N}, ket) where N

    exp_val = zero(ComplexF64) 

    for (oi,coeffi) in A.ops
        for (oj,coeffj) in B.ops
            exp_val += coeffi*coeffj * expectation_value(oi*oj,ket)*get_phase(oi,oj) 
        end
    end            

    return 2*real(exp_val)
end


