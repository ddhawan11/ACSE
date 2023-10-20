using PauliOperators 


"""
    hubbard_jw(graph::Array{T,2}, U, t)

TBW
"""
function hubbard_jw(graph::Matrix{T}, U, t) where T
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