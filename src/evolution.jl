using Printf

function BFS(generators::PauliSum{N}, H::PauliSum{N}, ket, grad; thresh=1e-4) where N

    o_transformed = deepcopy(H)
    
    for (g,g_coeff) in generators.ops

        sin_branch = PauliSum(N)

        for (oi,coeff) in o_transformed.ops

            println("g ", g)
            println("oi ", oi)
            
            if commute(oi, g) == false

                vcos = cos(g_coeff)
                vsin = sin(g_coeff)
                
                # cos branch
                o_transformed[oi] = coeff * vcos

                # sin branch
                oj = oi * g    # multiply the paulis
                oj = Pauli{N}(PauliOperators.phase(oi, g), oj)
                sum!(sin_branch, oj * vsin * coeff * 1im)

            end
        end
        sum!(o_transformed, sin_branch)
        println("o_trans", o_transformed)
        exit()
    end
    clip!(o_transformed, thresh=thresh)

    println(length(o_transformed))
    return o_transformed
end

function evolve_Hamiltonian(A, H, ket, bfs_thresh, grad_thresh)
    
    generator, curr_grad = ACSE.find_generator(A, H, ket)
    @printf("Number of Paulis in Hamiltonian operator: %i\n", length(H))
    generator *= curr_grad
    println("generator: ", generator)
    H_transformed = ACSE.BFS(generator, H, ket, curr_grad, thresh=bfs_thresh)

    while abs(curr_grad) > grad_thresh     
        println("curr_grad ", curr_grad)
        generator, curr_grad = ACSE.find_generator(A, H_transformed, ket)
        @printf("Number of Paulis in Hamiltonian operator: %i\n", length(H_transformed))

        H_new = ACSE.BFS(generator, H_transformed, ket, curr_grad, thresh=bfs_thresh)

        H_transformed = H_new
    end
    return H_transformed
end
