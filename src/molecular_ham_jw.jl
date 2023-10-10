using NPZ
using PauliOperators
using LinearAlgebra

function spatial_to_spin_orb(ints_1, ints_2)
    norb = size(ints_1)[1]
    nqubits = 2*norb

    h1 = zeros(nqubits,nqubits)
    h2 = zeros(nqubits, nqubits, nqubits, nqubits)

    for p in 1:norb
        for q in 1:norb
            h1[2*p-1, 2*q-1] = ints_1[p,q]
            h1[2*p, 2*q] = ints_1[p,q]
        
            for r in 1:norb
                for s in 1:norb
                    h2[2 * p-1, 2 * q, 2 * r, 2 * s-1] = ints_2[p,q,r,s]
                    h2[2 * p, 2 * q-1, 2 * r-1, 2 * s] = ints_2[p,q,r,s]
                    h2[2 * p-1, 2 * q-1, 2 * r-1, 2 * s-1] = ints_2[p,q,r,s]
                    h2[2 * p, 2 * q, 2 * r, 2 * s] = ints_2[p,q,r,s]
                end
            end
        end
    end
    h2 = 0.5 .* h2

    return(h1, h2)
end

function one_body_transform(h1)
    norb = size(h1)[1]
    oper   = PauliSum(norb)
    for i in 1:norb
        for j in 1:norb
            oper += (jordan_wigner(i,norb) * adjoint(jordan_wigner(j,norb))) * h1[i,j]
            clip!(oper)
        end
    end
    return oper
end

function two_body_transform(h2)
    norb = size(h2)[1]
    oper   = PauliSum(norb)
    println(norb)
    for i in 1:norb
        for j in 1:norb
            if i != j
                for k in 1:norb
                    for l in 1:norb
                        if k != l
                            oper += jordan_wigner(i,norb) * jordan_wigner(j,norb) * adjoint(jordan_wigner(k,norb)) * adjoint(jordan_wigner(l,norb)) * h2[i,j,k,l]
                            clip!(oper)
                        end
                    end
                end
            end
        end
    end
    return oper
end

function transform_molecular_Hamiltonian(data_dir="H2_0.75")

    ints_0 = npzread(data_dir*"/h0.npy");
    ints_1 = npzread(data_dir*"/h1.npy");
    ints_2 = npzread(data_dir*"/h2.npy");
    ints_2 = permutedims(ints_2, (3,2,1,4)) #Change to Physicist Notation
    println(ints_0)
    h1, h2 = spatial_to_spin_orb(ints_1, ints_2) #Convert to Spin Orbitals 
#    display(h2)
#    exit()
    qubit_ham = one_body_transform(h1) + two_body_transform(h2)
    clip!(qubit_ham, thresh=1e-15)

    return qubit_ham

end

