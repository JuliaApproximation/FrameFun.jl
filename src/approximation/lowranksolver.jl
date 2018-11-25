
"""
A Randomized Truncated SVD solver storing a decomposition of an operator A.
If A is MxN, and of rank R, the cost of constructing this solver is MR^2.

For tensor product operators it returns a decomposition of the linearized system
"""
struct RandomizedSvdSolver{T} <: FE_Solver{T}
    # Keep the original operator
    A          ::  DictionaryOperator
    # Random matrix
    W           ::  MultiplicationOperator
    # Decomposition
    Ut          ::  Array{T,2}
    Sinv        ::  Array{T,1}
    V           ::  Array{T,2}
    # For storing intermediate results when applying
    y           ::  Array{T,1}
    sy          ::  Array{T,1}
    scratch_src ::  Array{T,1}
    scratch_dest::  Array{T,1}

    smallcoefficients   :: Bool
    smalltol            :: Float64

    function RandomizedSvdSolver{T}(A::DictionaryOperator; threshold = default_threshold(A), rankestimate = 5, verbose = false, smallcoefficients = false, smalltol = 10, options...) where {T}
        W, U, S, Vt = randomizedsvd(A, threshold = threshold, rankestimate = rankestimate, verbose = verbose)

        Wsrc = DiscreteVectorDictionary{T}(size(W,2))
        Wdest = src(A)
        WM = MatrixOperator(Wsrc, Wdest, W)
        scratch_src = zeros(T, length(dest(A)))
        scratch_dest = zeros(T, length(src(A)))
        y = zeros(T, size(Vt,2))
        sy = zeros(T, length(S))
        new(A, WM, U', S.^(-1), Vt', y, sy, scratch_src, scratch_dest, smallcoefficients, smalltol)
    end
end

operator(s::RandomizedSvdSolver) = s.A

RandomizedSvdSolver(A::DictionaryOperator; options...) =
    RandomizedSvdSolver{eltype(A)}(A::DictionaryOperator; options...)

apply!(s::RandomizedSvdSolver, coef_dest, coef_src) = _apply!(s, coef_dest, coef_src,
    s.W, s.Ut, s.Sinv, s.V, s.y, s.sy, s.scratch_src, s.scratch_dest,
    s.smallcoefficients, s.smalltol)

# We don't need to (de)linearize coefficients when they are already vectors
function _apply!(s::RandomizedSvdSolver, coef_dest::Vector, coef_src::Vector,
    W, Ut, Sinv, V, y, sy, scratch_src, scratch_dest, smallcoefficients, smalltol)

    # Applying the truncated SVD to P*Rhs
    mul!(sy, Ut, coef_src)
    for i in 1:length(sy)
        sy[i] = sy[i]*Sinv[i]
    end
    if smallcoefficients
        M = maximum(abs.(coef_src))
        sy[abs.(sy) .> smalltol*M] = 0
    end
    mul!(y, V, sy)
    apply!(W, coef_dest, y)
end

function _apply!(s::RandomizedSvdSolver, coef_dest, coef_src,
    W, Ut, Sinv, V, y, sy, scratch_src, scratch_dest, smallcoefficients, smalltol)

    linearize_coefficients!(src(s), scratch_src, coef_src)
    # Call the implementation above for vectors
    _apply!(s, scratch_dest, scratch_src, W, Ut, Sinv, V, y, sy,
        scratch_src, scratch_dest, smallcoefficients, smalltol)
    delinearize_coefficients!(dest(s), coef_dest, scratch_dest)
end


"""
A Truncated SVD solver storing a decomposition of an operator A.
The cost is equal to the computation of the SVD of A.

For tensor product operators it returns a decomposition of the linearized system
"""
struct ExactSvdSolver{T} <: FE_Solver{T}
    # Keep the original operator
    A          ::  DictionaryOperator
    # Decomposition
    Ut          ::  Array{T,2}
    Sinv        ::  Array{T,1}
    V           ::  Array{T,2}
    # For storing intermediate results when applying
    y           ::  Array{T,1}
    sy          ::  Array{T,1}
    scratch_src ::  Array{T,1}

    function ExactSvdSolver{T}(A::DictionaryOperator; threshold = default_threshold(A), verbose = false,options...) where {T}
        C = matrix(A)
        m = maximum(abs.(C))
        USV = LAPACK.gesdd!('S',C)
        S = USV[2]
        maxind = findlast(S.>threshold*m)
        verbose && println("Solver truncated at singular value equal to ", S[maxind], " threshold was ", threshold, " and norm of A ", m)
        Sinv = 1 ./ S[1:maxind]
        y = zeros(T, size(USV[3],1))
        sy = zeros(T, maxind)

        scratch_src = zeros(T, length(dest(A)))
        new(A, USV[1][:,1:maxind]',Sinv,USV[3][1:maxind,:]',y,sy, scratch_src)
    end
end

src(t::ExactSvdSolver) = dest(t.A)
dest(t::ExactSvdSolver) = src(t.A)
inv(t::ExactSvdSolver) = t.A

ExactSvdSolver(A::DictionaryOperator; options...) =
    ExactSvdSolver{eltype(A)}(A::DictionaryOperator; options...)

# We don't need to (de)linearize coefficients when they are already vectors
function apply!(s::ExactSvdSolver, coef_dest::Vector, coef_src::Vector)
    # Applying the truncated SVD to P*Rhs
    mul!(s.sy, s.Ut, coef_src)
    for i =1:length(s.sy)
        s.sy[i]=s.sy[i]*s.Sinv[i]
    end
    mul!(coef_dest, s.V, s.sy)
end

function apply!(s::ExactSvdSolver, coef_dest, coef_src)
    linearize_coefficients!(src(s), s.scratch_src, coef_src)
    mul!(s.sy, s.Ut, s.scratch_src)
    for i =1:length(s.sy)
        s.sy[i]=s.sy[i]*s.Sinv[i]
    end
    mul!(s.y, s.V, s.sy)
    delinearize_coefficients!(dest(s), coef_dest, s.y)
end

# Function with equal functionality, but allocating memory
function truncatedsvd_solve(b::Vector, A::DictionaryOperator; threshold = default_threshold(A), R = 5, growth_factor = 2, verbose = false, smallcoefficients=false,smalltol=10,options...)
    finished=false
    T = eltype(A)
    R = min(R, size(A,2))
    random_matrix = map(T, rand(size(A,2), R))
    C = apply_multiple(A, random_matrix)
    c = cond(C)
    cold = threshold
    while (c < 1/threshold) && (R<size(A,2)) && (c>cold*10)
        verbose && println("c : $c\t cold : $cold\t threshold : $threshold")
        verbose && println("Solver truncated at R = ", R, " dof out of ",size(A,2))
        R0 = R
        R = min(round(Int,growth_factor*R),size(A,2))
        extra_random_matrix = rand(T, size(A,2), R-R0)
        Cextra = apply_multiple(A, extra_random_matrix)
        random_matrix = [random_matrix extra_random_matrix]
        # Extra condition: condition number has to increase by a significant amount each step, otherwise, possibly well conditioned.
        # cold = c
        C = [C Cextra]
        c = cond(C)
    end
    verbose && println("c : $c\t cold : $cold\t threshold : $threshold")
    verbose && println("Solver truncated at R = ", R, " dof out of ",size(A,2))


    y, rr = LAPACK.gelsy!(C,b[:],threshold)
    x1 = random_matrix*y;
    reshape(x1, size(src(A))...)
end
