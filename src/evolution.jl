using Printf

function BFS(generators::PauliSum{N}, H::PauliSum{N}, ket; thresh=1e-4) where N

    o_transformed = deepcopy(H)
    
    for (g,g_coeff) in generators.ops

        sin_branch = PauliSum(N)
        vcos = cos(g_coeff)
        vsin = sin(g_coeff)

        for (oi,coeff) in o_transformed.ops

            if PauliOperators.commute(oi, g) == false

                # cos branch
                o_transformed[oi] = coeff * vcos 

                # sin branch
                oj = oi * g    # multiply the paulis
                oj = Pauli{N}(PauliOperators.phase(oi, g), oj)
                sum!(sin_branch, oj * vsin * coeff * -1)

            end
        end
        sum!(o_transformed, sin_branch)

    end
    clip!(o_transformed, thresh=thresh)

    return o_transformed
end

function det_evolution(generators::PauliSum{N}, H::PauliSum{N}, ket; thresh=1e-3) where {N}

    #
    # for a single pauli Unitary, U = exp(-i θn Pn/2)
    # U' O U = cos(θ) O + i sin(θ) OP

    opers = H
    H_transformed = PauliSum(N)
    for (g,g_coeff) in generators.ops

        branch_opers = PauliSum(N)

        vcos = cos(g_coeff)
        vsin = sin(g_coeff)
#        sizehint!(branch_opers.ops, 1000)
        for (oi,oi_coeff) in opers.ops
            
            if PauliOperators.commute(oi, g)
                if haskey(branch_opers, oi)
                    branch_opers[oi] += oi_coeff
                else
                    branch_opers[oi] = oi_coeff
                end
            else
                # cos branch
                coeff = vcos * oi_coeff
                           
                if haskey(branch_opers, oi) # Add operator to dictionary if the key doesn't exist
                    branch_opers[oi] += coeff
                else
                    branch_opers[oi] = coeff #Modify the coeff if the key exists already
                end
            
                # sin branch
                coeff = vsin * oi_coeff * -1

                oj = Pauli{N}(0,oi) * Pauli{N}(0,g)    # multiply the pauli's

                if haskey(branch_opers, oj.pauli) # Add operator to dictionary if the key doesn't exist
                    branch_opers[oj.pauli] += coeff * (1im)^oj.θ
                else
                    branch_opers[oj.pauli] = coeff * (1im)^oj.θ  #Modify the coeff if the key exists already
                end
            end
            opers = deepcopy(branch_opers) # Change the list of operators to the next row
        end
    end
    H_transformed = opers
    clip!(H_transformed, thresh=thresh)
    return H_transformed
end

function exact_evolution(generators::PauliSum{N}, H::PauliSum{N}) where N
    
    o_mat = Matrix(H)
    U = Matrix(Pauli(N))

    for (g,g_coeff) in generators.ops
        α = g_coeff
        U = cos(α/2) .* U .- 1*sin(α/2) .* U * Matrix(g)
    end
    evolved_Hamiltonian = U'*o_mat*U

    ei = diag(U'*o_mat*U)
    
    println("ei", ei, "\n", "\n")
end


"""
    evolve_Hamiltonian(A, Hin, ket, bfs_thresh, grad_thresh; max_iter=100, verbose=1, alpha=.1)

TBW
"""
function evolve_Hamiltonian(A, Hin, ket, bfs_thresh, grad_thresh; max_iter=100, verbose=1, alpha=.1)

    H = deepcopy(Hin)
   
    generator_sequence = []
    generators, curr_grad, op_idx = ACSE.find_generator(A, H, ket)

    verbose < 1 || @printf("Number of Paulis in Hamiltonian operator: %i\n", length(H))

    generators *= curr_grad * alpha
    println("curr_grad ", curr_grad)
    gradients = []
    energies = []
    for iter in 1:max_iter
    
        push!(generator_sequence, (op_idx, curr_grad*alpha))
        H_transformed = BFS(generators, H, ket, thresh=bfs_thresh)
        # H_transformed = det_evolution(generators, H, ket, thresh=bfs_thresh)
        # exact_evolution(generators, H)

         
        energy = ACSE.calc_energy(H_transformed, ket)
        verbose < 1 || @printf("Iter: %4i Energy: %10.8f + %10.8fi Gradient: %6.1e #Ops in H: %7i Op IDX: %4i\n", iter, real(energy),imag(energy), abs(curr_grad), length(H_transformed), op_idx)
        H = H_transformed
        # generators, curr_grad, op_idx = ACSE.find_generator(A, H, ket)
        generators, curr_grad, op_idx = ACSE.find_generator2(A, H, ket)
        generators *= curr_grad * alpha 

        push!(energies, energy)
        push!(gradients, curr_grad)
        abs(curr_grad) > grad_thresh || break 
    end
    
    return H, energies, gradients, generator_sequence
end
