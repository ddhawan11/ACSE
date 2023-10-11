using Printf

function BFS(generators::PauliSum{N}, H::PauliSum{N}, ket, grad; thresh=1e-4) where N

    o_transformed = deepcopy(H)
    vcos = cos(grad)
    vsin = sin(grad)
    
    for (g,g_coeff) in generators.ops

        sin_branch = PauliSum(N)
#        vcos = cos(g_coeff)
#        vsin = sin(g_coeff)

        for (oi,coeff) in o_transformed.ops

#            println("g ", g)
#            println("oi ", oi)
            
            if commute(oi, g) == false

#                vcos = cos(g_coeff)
#                vsin = sin(g_coeff)
                
                # cos branch
                o_transformed[oi] = coeff * vcos 

                # sin branch
                oj = oi * g    # multiply the paulis
                oj = Pauli{N}(PauliOperators.phase(oi, g), oj)
                sum!(sin_branch, oj * vsin * coeff * 1im * g_coeff)

            end
        end
        sum!(o_transformed, sin_branch)
#        println("o_trans", o_transformed)
#        exit()
    end
    clip!(o_transformed, thresh=thresh)

    println(length(o_transformed))
    return o_transformed
end

function det_evolution(generators::PauliSum{N}, H::PauliSum{N}, ket, grad ; thresh=1e-3) where {N}

    #
    # for a single pauli Unitary, U = exp(-i θn Pn/2)
    # U' O U = cos(θ) O + i sin(θ) OP
#    vcos = cos(grad)
#    vsin = sin(grad)


    opers = H#PauliSum(o)
    H_transformed = PauliSum(N)
    for (g,g_coeff) in generators.ops

        branch_opers = PauliSum(N)

        vcos = cos(g_coeff)
        vsin = sin(g_coeff)
#        sizehint!(branch_opers.ops, 1000)
        for (oi,oi_coeff) in opers.ops
            
            if commute(oi, g)
                if haskey(branch_opers, oi)
                    branch_opers[oi] += oi_coeff
                else
                    branch_opers[oi] = oi_coeff
                end
                continue
            end
            # cos branch
            coeff = vcos * oi_coeff
                           
            if haskey(branch_opers, oi) # Add operator to dictionary if the key doesn't exist
                branch_opers[oi] += coeff
            else
                branch_opers[oi] = coeff #Modify the coeff if the key exists already
            end
            
            # sin branch
            coeff = vsin * 1im * oi_coeff
            oi = oi * g    # multiply the pauli's

            if haskey(branch_opers, oi) # Add operator to dictionary if the key doesn't exist
                branch_opers[oi] += coeff * (1im)^PauliOperators.phase(oi,g)
            else
                branch_opers[oi] = coeff * (1im)^PauliOperators.phase(oi,g) #Modify the coeff if the key exists already
            end
            opers = deepcopy(branch_opers) # Change the list of operators to the next row
        end
        H_transformed += opers
    end
    return opers
end

function evolve_Hamiltonian(A, H, ket, bfs_thresh, grad_thresh)
    
    generators, curr_grad = ACSE.find_generator(A, H, ket)
    @printf("Number of Paulis in Hamiltonian operator: %i\n", length(H))

    curr_grad = 0.057416530623 ## Hard-coded from adapt for H2
    generators *= curr_grad
    println("curr_grad ", curr_grad)

    while abs(curr_grad) > grad_thresh     

#        H_transformed = BFS(generators, H, ket, curr_grad, thresh=bfs_thresh)
        H_transformed = det_evolution(generators, H, ket, curr_grad, thresh=bfs_thresh)

        energy = ACSE.calc_energy(H_transformed, ket)
        @printf("Energy: %10.8f\n", real(energy))
        exit()
        H = H_transformed
        
        generators, curr_grad = ACSE.find_generator(A, H, ket)
        @printf("Number of Paulis in Hamiltonian operator: %i\n", length(H))
        generators *= curr_grad
        println("curr_grad ", curr_grad)

    end
    
    return H
end
