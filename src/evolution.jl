function BFS(generators::PauliSum{N}, H::PauliSum{N}, ket, grad) where N
    vcos = cos(grad)
    vsin = sin(grad)

    o_transformed = deepcopy(o)
  
    n_ops = zeros(Int,nt)
    
    for (key,value) in generators.ops

        g = key

        sin_branch = PauliSum(N)

        for (oi,coeff) in o_transformed.ops
            
            if commute(oi, g.pauli) == false
                
                # cos branch
                o_transformed[oi] = coeff * vcos

                # sin branch
                oj = g * oi    # multiply the paulis
                sum!(sin_branch, oj * vsin * coeff * 1im)

            end
        end
        sum!(o_transformed, sin_branch) 
        clip!(o_transformed, thresh=thresh)
        n_ops[t] = length(o_transformed)
    end
    return n_ops

end
