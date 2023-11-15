# Calculates the derivative of shannon entropy of Hamiltonian entropy
function calc_ds(H::PauliSum{N}, A) where N
    deriv_s = 0.0
    for a in A
        for (oi,oi_coeff) in H.ops
            l2_norm, dpk = calc_dpk(H, oi, oi_coeff, a)
            deriv_s -= dpk * log(oi_coeff^2/l2_norm) - (1/log(2))*dpk
        end
        println(deriv_s)
    end
    return deriv_s
end    

# This function calculates the derivative of probability
## WARNING : It assumes you are using qubit pool(Will need to work with multiple generators for the fermionic pool)

function calc_dpk(H::PauliSum{N}, oi::FixedPhasePauli{N}, oi_coeff, a) where N
    l2_norm = oper_norm(H)
    sum_coeff = 0.0
    for (g, g_coeff) in a.ops
        for (hi,pi) in H.ops
            if g*hi == oi
                sum_coeff += 1im*pi*g_coeff                    
            end
        end
    end

#    println("oi_coeff: ", oi_coeff)
#    println("sum_coeff: ", sum_coeff)
    dpk = 2*(1/l2_norm)*sum_coeff*oi_coeff
    return l2_norm, dpk
end


# Calculates entropy for an operator
function calc_entropy(H)
    l2_norm = oper_norm(H)
    s = 0.0
    for (oi,oi_coeff) in H.ops
        pk = oi_coeff^2/l2_norm
        s -= pk*log(pk) 
    end
    return(s)
end
