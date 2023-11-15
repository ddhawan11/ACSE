using LinearAlgebra

## Calculate Expectation Value
#Calculates the expectation value for H operator, H = \sum_{k} c_{k}P_{k}, where P_{k} is a Pauli
function calc_energy(H, ket)
    energy = zero(ComplexF64)
    for (key,value) in H.ops
        energy += value * PauliOperators.expectation_value(key, ket)
    end
    imag(energy) ≈ 0 || throw("Complex expectation value")
    return energy
end


#Calculates the l2 norm of the operator H (Normalization factor not included)

##There is a need to figure out why the Operator norm and the Matrix norm don't match.
function oper_norm(H::PauliSum{N}) where N
    norm = 0
    for (key,value) in H.ops
        norm += abs(value^2)
    end
    
    norm = sqrt(norm)
    return norm
end



# Calculate the squared l2-norm of a matrix, same as opnorm(H_mat, 2)
function calc_norm(H_mat)
    norm = 0
    for i in size(H_mat)[1]
        for j in size(H_mat)[2]
            norm += H_mat[i,j]
        end
    end
    return norm
end


