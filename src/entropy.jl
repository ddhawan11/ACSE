# Derivative of entropy through finite difference
function calc_ds_finite_diff(H, a)
    bfs_thresh = 1e-10
    s_H = calc_entropy(H)
    theta = 0.0
    h = 0.001
    a_n = (theta+h) * a
    a_p = (theta-h) * a
    H_transformed_n = ACSE.BFS(a_n, H, thresh=bfs_thresh)
    H_transformed_p = ACSE.BFS(a_p, H, thresh=bfs_thresh)

    ds = (calc_entropy(H_transformed_n) - calc_entropy(H_transformed_p))/(2*h)
    return ds
end


# Calculates the derivative of shannon entropy of Hamiltonian entropy
function calc_ds(H::PauliSum{N}, A) where N
    for a in A
        deriv_s = 0.0
        for (oi,oi_coeff) in H.ops
            l2_norm, dpk = calc_dpk(H, oi, oi_coeff, a)
            deriv_s -= (dpk * log(oi_coeff^2/l2_norm) - (1/log(2))*dpk)
        end
        ds_fd = calc_ds_finite_diff(H, a)
        @printf("Derivative: %10.8f, From Finite Difference: %10.8f \n", deriv_s, ds_fd)
    end
#    return deriv_s
end    

# This function calculates the derivative of probability
## WARNING : It assumes you are using qubit pool(Will need to work with multiple paulis in a generator for the fermionic pool)

function calc_dpk(H::PauliSum{N}, oi::FixedPhasePauli{N}, oi_coeff, a) where N
    l2_norm = oper_norm(H)
    sum_deriv = 0.0
    for (g, g_coeff) in a.ops
        for (pi,hi) in H.ops
            if pi*g == oi
                sum_deriv -= hi*g_coeff*get_phase(pi,g)
            end
        end
    end

#    println("oi_coeff: ", oi_coeff)
#    println("sum_coeff: ", sum_coeff)
    dpk = 2*(1/l2_norm)*sum_deriv*oi_coeff
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


