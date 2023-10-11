using PauliOperators

function acse_residual_pool(n_occ, n_vir)
    norb = n_occ + n_vir
    acse_res = Vector{PauliSum{norb}}([])
    fermi_count = 0
    ## Single excitations
    for i in 1:2:norb
        for j in i+2:2:norb
            oper = jordan_wigner(j,norb) * adjoint(jordan_wigner(i,norb))
            oper -= jordan_wigner(i,norb) * adjoint(jordan_wigner(j,norb))
            oper += jordan_wigner(j+1,norb) * adjoint(jordan_wigner(i+1,norb))
            oper -= jordan_wigner(i+1,norb) * adjoint(jordan_wigner(j+1,norb))
            clip!(oper)
            if length(oper) != 0
                push!(acse_res, oper)
            end
            fermi_count+=1
        end
    end

    ## Double excitations a-a or b-b type
    ij = 0
    for i in 1:2:norb
        for j in i:2:norb
            kl = 0
            for k in 1:2:norb
                for l in k:2:norb
                    if ij < kl
                        continue
                    end
                    oper = jordan_wigner(k,norb) * adjoint(jordan_wigner(i,norb)) * jordan_wigner(l,norb) * adjoint(jordan_wigner(j,norb))
                    oper -= jordan_wigner(j,norb) * adjoint(jordan_wigner(l,norb)) * jordan_wigner(i,norb) * adjoint(jordan_wigner(k,norb))
                    oper += jordan_wigner(k+1,norb) * adjoint(jordan_wigner(i+1,norb)) * jordan_wigner(l+1,norb) * adjoint(jordan_wigner(j+1,norb))
                    oper -= jordan_wigner(j+1,norb) * adjoint(jordan_wigner(l+1,norb)) * jordan_wigner(i+1,norb) * adjoint(jordan_wigner(k+1,norb))
                    clip!(oper)
                    fermi_count+=1
                    if length(oper) != 0
                        push!(acse_res, oper)
                    end
                    kl += 1
                end
            end
            ij += 1
        end
    end

    ## Double excitations a-b type
    ij = 0
    for i in 1:2:norb
        for j in 2:2:norb
            kl = 0
            for k in 1:2:norb
                for l in 2:2:norb
                    if ij < kl
                        continue
                    end
                    oper = jordan_wigner(k,norb) * adjoint(jordan_wigner(i,norb)) * jordan_wigner(l,norb) * adjoint(jordan_wigner(j,norb))
                    oper -= jordan_wigner(j,norb) * adjoint(jordan_wigner(l,norb)) * jordan_wigner(i,norb) * adjoint(jordan_wigner(k,norb))                     
                    if i>j
                        continue
                    end
                    oper += jordan_wigner(l-1,norb) * adjoint(jordan_wigner(j-1,norb)) * jordan_wigner(k+1,norb) * adjoint(jordan_wigner(i+1,norb))
                    oper -= jordan_wigner(i+1,norb) * adjoint(jordan_wigner(k+1,norb)) * jordan_wigner(j-1,norb) * adjoint(jordan_wigner(l-1,norb))
                    fermi_count+=1
                    clip!(oper)
                    if length(oper) != 0
                        push!(acse_res, oper)
                    end
                    kl += 1
                end
            end
            ij += 1
        end
    end
    println("Number of fermionic operators:", fermi_count)
    return acse_res
end

function find_generator(A::Vector{PauliSum{N}}, H::PauliSum{N}, ket::KetBitString{N}) where N
    max_grad = 0.0
    curr_oper = PauliSum(N)
    i = 1
    j = 1
    for op in A
        gradient = ACSE.compute_grad(H,op,ket)
#        if abs(gradient) > 1e-8
#            println("finding the operator: ", gradient)
#        end
        gradient = real(gradient)
        if abs(gradient) > abs(max_grad)
            curr_oper = op
            max_grad = gradient
            j = i
        end
        i+=1
    end
    println("Current Index:", j)
    return curr_oper, max_grad
end
