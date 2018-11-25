
apply_multiple(A::AbstractMatrix, B::AbstractMatrix) = A*B

svd_threshold(::Type{T}) where {T} = 10^(T(4)/5*log10(eps(T)))
svd_threshold(::Type{Complex{T}}) where {T} = svd_threshold(T)

function randomizedsvd(A; threshold = svd_threshold(eltype(A)), verbose = false, rankestimate = 16, growth_factor = 2, maxiter = 8)
    verbose && println("Computing randomized SVD of $(typeof(A))")
    T = eltype(A)

    R = min(rankestimate, size(A,2))
    W = rand(T, size(A,2), R)
    C = apply_multiple(A, W)
    sigma = svdvals(C)[end]

    iter = 0
    while (sigma > threshold) && (R < size(A,2)) && (iter < maxiter)
        verbose && println("Smallest singular value : $sigma\t Threshold : $threshold")
        verbose && println("SVD truncated at R = ", R, " out of ", size(A,2))
        R1 = R
        R = min(round(Int,growth_factor*R), size(A,2))
        W_extra = rand(T, size(A,2), R-R1)
        C_extra = apply_multiple(A, W_extra)
        W = [W W_extra]
        C = [C C_extra]
        sigma = svdvals(C)[end]
        iter += 1
    end
    if iter == maxiter
        warning("Maximum number of iteration reached in randomizedsvd")
    end
    USV = svd(C, full=false)
    U = USV.U
    S = USV.S
    Vt = USV.Vt
    rank = findlast(S .> threshold)
    if rank == nothing
        verbose && println("No singular values larger than $threshold found")
        rank = 1
    end
    verbose && println("Smallest singular value : $(S[rank])\t Threshold : $threshold")
    verbose && println("Done: SVD truncated with rank ", rank, " out of ",size(A,2))
    (W, U[:,1:rank], S[1:rank], Vt[1:rank,:])
end
