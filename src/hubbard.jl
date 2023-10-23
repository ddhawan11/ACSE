using PauliOperators 


"""
    hubbard_jw(graph::Array{T,2}, U, t)

TBW
"""
function hubbard_hamiltonian(graph::Matrix{T}, U, t) where T
    Ni, Nj = size(graph)
    Ni == Nj || throw(DimensionMismatch)
    Norb = Ni
    N = 2*Norb
    H = PauliSum(N)
    for i in 1:Norb
        ia = 2*i-1
        ib = 2*i
        for j in i+1:Norb
            ja = 2*j-1
            jb = 2*j
            abs(graph[i,j]) > 1e-16 || continue

            tij = jordan_wigner(ia, N)*jordan_wigner(ja,N)'
            tij += jordan_wigner(ib, N)*jordan_wigner(jb,N)'
            tij += adjoint(tij)
            clip!(tij)
            sum!(H,t*tij)
        end
            
        ni = jordan_wigner(ia, N)*jordan_wigner(ia,N)'*jordan_wigner(ib, N)*jordan_wigner(ib,N)'
        sum!(H,U*ni)
    end
    return H
end


"""
    hubbard_number(graph::Matrix{T}) where T

TBW
"""
function hubbard_number(graph::Matrix{T}) where T
    Ni, Nj = size(graph)
    Ni == Nj || throw(DimensionMismatch)
    Norb = Ni
    N = 2*Norb
    Na = PauliSum(N)
    Nb = PauliSum(N)
    for i in 1:Norb
        ia = 2*i-1
        ib = 2*i
        sum!(Na, jordan_wigner(ia, N)*jordan_wigner(ia,N)')
        sum!(Nb, jordan_wigner(ib, N)*jordan_wigner(ib,N)')
    end
    return Na, Nb
end

"""
    hubbard_pool(graph::Matrix{T}) where T

TBW
"""
function hubbard_pool(graph::Matrix{T}) where T
    Ni, Nj = size(graph)
    Ni == Nj || throw(DimensionMismatch)
    Norb = Ni
    N = 2*Norb
    G = Vector{PauliSum{N}}([])
    for i in 1:Norb
        ia = 2*i-1
        ib = 2*i
        for j in i+1:Norb
            ja = 2*j-1
            jb = 2*j
            abs(graph[i,j]) > 1e-16 || continue

            tij = jordan_wigner(ia, N)*jordan_wigner(ja,N)'
            tij += jordan_wigner(ib, N)*jordan_wigner(jb,N)'
            tij -= adjoint(tij)
            clip!(tij)
            push!(G, tij)
            
            tij =  jordan_wigner(ia, N)*jordan_wigner(ia, N)*jordan_wigner(ja,N)'*jordan_wigner(ja,N)'
            tij += jordan_wigner(ib, N)*jordan_wigner(ib, N)*jordan_wigner(jb,N)'*jordan_wigner(jb,N)'
            tij -= adjoint(tij)
            clip!(tij)
            push!(G, tij)
            
            #tij = jordan_wigner(ia, N)*jordan_wigner(ja,N)'*jordan_wigner(ib, N)*jordan_wigner(jb,N)'
            #tij -= adjoint(tij)
            #clip!(tij)
            #push!(G, tij)
            
            #tij = jordan_wigner(ia, N)*jordan_wigner(ia,N)'*jordan_wigner(jb, N)*jordan_wigner(jb,N)'
            #tij += jordan_wigner(ib, N)*jordan_wigner(ib,N)'*jordan_wigner(ja, N)*jordan_wigner(ja,N)'
            #tij -= adjoint(tij)
            #clip!(tij)
            #push!(G, tij)
        end
            
        ni = jordan_wigner(ia, N)*jordan_wigner(ia,N)'*jordan_wigner(ib, N)*jordan_wigner(ib,N)'
        push!(G, 1im*ni)
    end
    return G
end


function hubbard_pool2(graph::Matrix{T}) where T
    Ni, Nj = size(graph)
    Ni == Nj || throw(DimensionMismatch)
    Norb = Ni
    N = 2*Norb
    G = Vector{PauliSum{N}}([])
    Gmix = Vector{PauliSum{N}}([])
    for i in 1:Norb
        ia = 2*i-1
        ib = 2*i
        for j in i+1:Norb
            ja = 2*j-1
            jb = 2*j
            abs(graph[i,j]) > 1e-16 || continue

            tij = jordan_wigner(ia, N)*jordan_wigner(ja,N)'
            tij += jordan_wigner(ib, N)*jordan_wigner(jb,N)'
            tij -= adjoint(tij)
            clip!(tij)
            push!(G, tij)
            
        end
            
        ni = jordan_wigner(ia, N)*jordan_wigner(ia,N)'*jordan_wigner(ib, N)*jordan_wigner(ib,N)'
        push!(Gmix, 1im*ni)
    end
    return G,Gmix
end


