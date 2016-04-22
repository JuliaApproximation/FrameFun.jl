# lowrank.jl


function lowrank_approximation(A; ɛ = 1000eps(eltype(A)), estimated_rank = 2)
    ELT = eltype(A)
    (m,n) = size(A)
    r = min(estimated_rank,n)
    random_matrix = map(ELT, rand(n, r))
    C = apply_multiple(A, random_matrix)
    c = cond(C)
    m = maximum(abs(C))
    while (c < m/ɛ) && (r < n)
        r0 = r
        r = min(2r,n)
        extra_random_matrix = map(ELT, rand(n,r-r0))
        extraC = apply_multiple(A, extra_random_matrix)
        C = [C extraC]
        random_matrix = [random_matrix extra_random_matrix]
        c = cond(C)
        m = maximum(abs(C))
    end
    USV = LAPACK.gesvd!('S','S',C)
    S = USV[2]
    rank = findlast(S.>ɛ)
    u = USV[1][:,1:rank]
    s = S[1:rank]
    sinv = s.^(-1)
    v = USV[3][1:rank,:]
    ut = u'
    vs = v' * diagm(sinv)
    ut,vs,random_matrix
end

function lowranksolve(lr, b)
    ut,vs,m = lr
    m * (vs * (ut * b))
end
