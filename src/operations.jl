function calc_energy(H, ket)
    energy = zero(ComplexF64)
    for (key,value) in H.ops
        energy += value * PauliOperators.expectation_value(key, ket)
    end
    return energy
end

function reference_state(norb)
    n_spinorb    = norb * 2
    state_binary = 0
    for i in 1:norb
        state_binary += 2^(i-1)
    end
#    println(state_binary)

    KetBitString(n_spinorb, state_binary)

end

