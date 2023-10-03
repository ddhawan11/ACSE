using PauliOperators

function acse_residual_pool(n_occ, n_vir)
    norb = n_occ + n_vir
    acse_res = Vector{PauliSum{norb}}([])
    op = PauliSum(norb)
    for i in 1:n_occ
        for j in 1:n_occ
            if i != j
                for k in n_occ+1:norb
                    for l in n_occ+1:norb
                        if k != l
                            oper = jordan_wigner(k,norb) * jordan_wigner(l,norb) * adjoint(jordan_wigner(i,norb)) * adjoint(jordan_wigner(j,norb))
                            oper = clip!(oper - adjoint(oper))
                            if length(oper) != 0
                                push!(acse_res, oper)
                            end
                        end
                    end
                end
            end
        end
    end
    return acse_res
end

