using PauliOperators

# The spin orbitals are labelled as 1a,1b,2a,2b,3a,3b..

## Spin complemented GSD
function acse_residual_pool(N)
    norb = N
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
                    oper = jordan_wigner(k,norb) * (jordan_wigner(l,norb)) * adjoint(jordan_wigner(i,norb)) * adjoint(jordan_wigner(j,norb))
                    oper -= jordan_wigner(j,norb) * (jordan_wigner(i,norb)) * adjoint(jordan_wigner(l,norb)) * adjoint(jordan_wigner(k,norb))
                    oper += jordan_wigner(k+1,norb) * (jordan_wigner(l+1,norb)) * adjoint(jordan_wigner(i+1,norb)) * adjoint(jordan_wigner(j+1,norb))
                    oper -= jordan_wigner(j+1,norb) * (jordan_wigner(i+1,norb)) * adjoint(jordan_wigner(l+1,norb)) * adjoint(jordan_wigner(k+1,norb))
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
                    oper = jordan_wigner(k,norb) * (jordan_wigner(l,norb)) * adjoint(jordan_wigner(i,norb)) * adjoint(jordan_wigner(j,norb))
                    oper -= jordan_wigner(j,norb) * (jordan_wigner(i,norb)) * adjoint(jordan_wigner(l,norb)) * adjoint(jordan_wigner(k,norb))                     
                    if i>j
                        continue
                    end
                    oper += jordan_wigner(l-1,norb) * (jordan_wigner(k+1,norb)) * adjoint(jordan_wigner(j-1,norb)) * adjoint(jordan_wigner(i+1,norb))
                    oper -= jordan_wigner(i+1,norb) * (jordan_wigner(j-1,norb)) * adjoint(jordan_wigner(k+1,norb)) * adjoint(jordan_wigner(l-1,norb))
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


function qubit_pool(N)

    oper_pool = Vector{PauliSum{N}}([])
    oper = PauliSum{N}
    for i in 1:2:N
        for j in i+2:2:N
            oper = PauliSum(Pauli(N, X=[i], Y=[j]))*-1im
            push!(oper_pool,oper)
            oper = PauliSum(Pauli(N, Y=[i], X=[j]))*-1im
            push!(oper_pool,oper)
            oper = PauliSum(Pauli(N, X=[i+1], Y=[j+1]))*-1im
            push!(oper_pool,oper)
            oper = PauliSum(Pauli(N, Y=[i+1], X=[j+1]))*-1im
            push!(oper_pool,oper)
        end
    end

    for i in 1:N
        for j in i+1:N
            for k in j+1:N
                for l in k+1:N
                    if (i+j+k+l) % 2 == 0
                        oper = PauliSum(Pauli(N, Y=[i], X=[j,k,l]))*-1im
                        push!(oper_pool, oper)
                        oper = PauliSum(Pauli(N, Y=[j], X=[i,k,l]))*-1im
                        push!(oper_pool, oper)
                        oper = PauliSum(Pauli(N, Y=[k], X=[i,j,l]))*-1im
                        push!(oper_pool, oper)
                        oper = PauliSum(Pauli(N, Y=[l], X=[i,j,k]))*-1im
                        push!(oper_pool, oper)
                        oper = PauliSum(Pauli(N, X=[i], Y=[j,k,l]))*-1im
                        push!(oper_pool, oper)
                        oper = PauliSum(Pauli(N, X=[j], Y=[i,k,l]))*-1im
                        push!(oper_pool, oper)
                        oper = PauliSum(Pauli(N, X=[k], Y=[i,j,l]))*-1im
                        push!(oper_pool, oper)
                        oper = PauliSum(Pauli(N, X=[l], Y=[i,j,k]))*-1im
                        push!(oper_pool, oper)
                    end
                end
            end
        end
    end
    
    return oper_pool
end


function find_generator(A::Vector{PauliSum{N}}, H::PauliSum{N}, ket::KetBitString{N}) where N
    max_grad = 0.0
    curr_oper = PauliSum(N)
    i = 1
    j = 1
    for op in A
        gradient = ACSE.compute_grad(H,op,ket)

        if abs(gradient) > abs(max_grad)
            curr_oper = op
            max_grad = gradient
            j = i
        end
        i+=1
    end
    imag(max_grad) ≈ 0 || throw("Complex gradients")

    return curr_oper, max_grad, j
end

function find_generator2(A::Vector{PauliSum{N}}, H::PauliSum{N}, ket::KetBitString{N}) where N
    max_grad = 0.0
    curr_oper = PauliSum(N)
    i = 1
    j = 1
    for op in A
        gradient = ACSE.compute_grad2(H,op,ket)

        if abs(gradient) > abs(max_grad)
            curr_oper = op
            max_grad = gradient
            j = i
        end
        i+=1
    end
    imag(max_grad) ≈ 0 || throw("Complex gradients")

    return curr_oper, max_grad, j
end
