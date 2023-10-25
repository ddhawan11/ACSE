## Calculate Expectation Value

function calc_energy(H, ket)
    energy = zero(ComplexF64)
    for (key,value) in H.ops
        energy += value * PauliOperators.expectation_value(key, ket)
    end
    imag(energy) ≈ 0 || throw("Complex expectation value")
    return energy
end


