function delta(i,j)
    if i==j
        return 1
    else
        return 0
    end
end
    
function build_K_matrix(h1, h2, Ne)

    norb = size(h1)[1]
    K    = zeros(size(h2))
    K   += h2
    for p in 1:norb
        for q in 1:norb
            for s in 1:norb
                for t in 1:norb
                    K[t,s,q,p] += (delta(q,s)*h1[t,p]+delta(t,p)*h1[q,s])/(4*(Ne-1))
                    K[t,s,q,p] -= (delta(p,s)*h1[t,q]+delta(t,q)*h1[p,s])/(4*(Ne-1))
                    
                end
            end
        end
    end
#    display(K)

    return K
end    

function transform_contracted_Hamiltonian(data_dir="H2_0.75", Ne = N_el)

    ints_0 = npzread(data_dir*"/h0.npy");
    ints_1 = npzread(data_dir*"/h1.npy");
    ints_2 = npzread(data_dir*"/h2.npy");
    ints_2 = permutedims(ints_2, (3,2,1,4)) #Change to Physicist Notation
    h1, h2 = spatial_to_spin_orb(ints_1, ints_2) #Convert to Spin Orbitals 

    K = build_K_matrix(h1, h2, Ne)
    
    norb = size(ints_1)[1]
    nuc_rep = PauliSum(Dict(FixedPhasePauli(join(["I" for i in 1:norb*2]))=>complex(ints_0)))
    
    qubit_ham = nuc_rep + ACSE.two_body_transform(K) # Add nuclear repulsion    
    clip!(qubit_ham, thresh=1e-15)

    return qubit_ham

end
