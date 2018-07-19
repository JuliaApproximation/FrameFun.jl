# lowranksolver.jl

"""
A Randomized Truncated SVD solver storing a decomposition of an operator A.
If A is MxN, and of rank R, the cost of constructing this solver is MR^2.

For tensor product operators it returns a decomposition of the linearized system
"""
struct TruncatedSvdSolver{ELT} <: FE_Solver{ELT}
    # Keep the original operator
    A          ::  DictionaryOperator
    # Random matrix
    W           ::  MultiplicationOperator
    # Decomposition
    Ut          ::  Array{ELT,2}
    Sinv        ::  Array{ELT,1}
    V          ::  Array{ELT,2}
    # For storing intermediate results when applying
    y           ::  Array{ELT,1}
    sy          ::  Array{ELT,1}
    scratch_src ::  Array{ELT,1}
    scratch_dest::  Array{ELT,1}
    smallcoefficients :: Bool
    smalltol :: Float64


    function TruncatedSvdSolver{ELT}(A::DictionaryOperator; cutoff = default_cutoff(A), R = 5, growth_factor = 2, verbose = false, smallcoefficients=false,smalltol=10,options...) where ELT
        finished=false
        USV = ()
        R = min(R, size(A,2))
        random_matrix = map(ELT, rand(size(A,2), R))
        C = apply_multiple(A, random_matrix)
        c = cond(C)
# <<<<<<< HEAD
#         # TODO change things with cold back
#         # cold = 1
#         m = maximum(abs.(C))
#         while (c < 1/cutoff) && (R<size(op,2)) #&& (c/cold > 10)
# =======
        cold = cutoff
        m = maximum(abs.(C))
        while (c < 1/cutoff) && (R<size(A,2)) && (c>cold*10)
            verbose && println("c : $c\t cold : $cold\t cutoff : $cutoff")
            verbose && println("Solver truncated at R = ", R, " dof out of ",size(A,2))
            R0 = R
            R = min(round(Int,growth_factor*R),size(A,2))
            extra_random_matrix = map(ELT, rand(size(A,2), R-R0))
            Cextra = apply_multiple(A, extra_random_matrix)
            random_matrix = [random_matrix extra_random_matrix]
            # Extra condition: condition number has to increase by a significant amount each step, otherwise, possibly well conditioned.
            # cold = c
            C = [C Cextra]
            c = cond(C)

        end
        verbose && println("c : $c\t cold : $cold\t cutoff : $cutoff")
        verbose && println("Solver truncated at R = ", R, " dof out of ",size(A,2))
        USV = LAPACK.gesdd!('S',C)
        S = USV[2]
        maxind = findlast(S.>(maximum(S)*cutoff))
        Sinv = 1 ./ S[1:maxind]
        y = zeros(ELT, size(USV[3],1))
        sy = zeros(ELT, maxind)
        Wsrc = DiscreteVectorDictionary{ELT}(size(random_matrix,2))
        Wdest = src(A)
        W = MatrixOperator(Wsrc, Wdest, random_matrix)

        scratch_src = zeros(ELT, length(dest(A)))
        scratch_dest = zeros(ELT, length(src(A)))
        new(A, W, USV[1][:,1:maxind]',Sinv,USV[3][1:maxind,:]',y,sy, scratch_src, scratch_dest,smallcoefficients,smalltol)
    end
end

src(t::TruncatedSvdSolver) = dest(t.A)
dest(t::TruncatedSvdSolver) = src(t.A)
inv(t::TruncatedSvdSolver) = t.A

TruncatedSvdSolver(A::DictionaryOperator; options...) =
    TruncatedSvdSolver{eltype(A)}(A::DictionaryOperator; options...)

apply!(s::TruncatedSvdSolver, coef_dest, coef_src) = _apply!(s, coef_dest, coef_src,
    s.W, s.Ut, s.Sinv, s.V, s.y, s.sy, s.scratch_src, s.scratch_dest,
    s.smallcoefficients, s.smalltol)

# We don't need to (de)linearize coefficients when they are already vectors
function _apply!(s::TruncatedSvdSolver, coef_dest::Vector, coef_src::Vector,
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

function _apply!(s::TruncatedSvdSolver, coef_dest, coef_src,
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
struct ExactTruncatedSvdSolver{ELT} <: FE_Solver{ELT}
    # Keep the original operator
    A          ::  DictionaryOperator
    # Decomposition
    Ut          ::  Array{ELT,2}
    Sinv        ::  Array{ELT,1}
    V          ::  Array{ELT,2}
    # For storing intermediate results when applying
    y           ::  Array{ELT,1}
    sy          ::  Array{ELT,1}
    scratch_src ::  Array{ELT,1}

    function ExactTruncatedSvdSolver{ELT}(A::DictionaryOperator; cutoff = default_cutoff(A), verbose = false,options...) where ELT
        C = matrix(A)
        m = maximum(abs.(C))
        USV = LAPACK.gesdd!('S',C)
        S = USV[2]
        maxind = findlast(S.>cutoff*m)
        verbose && println("Solver truncated at singular value equal to ", S[maxind], " cutoff was ", cutoff, " and norm of A ", m)
        Sinv = 1 ./ S[1:maxind]
        y = zeros(ELT, size(USV[3],1))
        sy = zeros(ELT, maxind)

        scratch_src = zeros(ELT, length(dest(A)))
        new(A, USV[1][:,1:maxind]',Sinv,USV[3][1:maxind,:]',y,sy, scratch_src)
    end
end

src(t::ExactTruncatedSvdSolver) = dest(t.A)
dest(t::ExactTruncatedSvdSolver) = src(t.A)
inv(t::ExactTruncatedSvdSolver) = t.A

ExactTruncatedSvdSolver(A::DictionaryOperator; options...) =
    ExactTruncatedSvdSolver{eltype(A)}(A::DictionaryOperator; options...)

# We don't need to (de)linearize coefficients when they are already vectors
function apply!(s::ExactTruncatedSvdSolver, coef_dest::Vector, coef_src::Vector)
    # Applying the truncated SVD to P*Rhs
    mul!(s.sy, s.Ut, coef_src)
    for i =1:length(s.sy)
        s.sy[i]=s.sy[i]*s.Sinv[i]
    end
    mul!(coef_dest, s.V, s.sy)
end

function apply!(s::ExactTruncatedSvdSolver, coef_dest, coef_src)
    linearize_coefficients!(src(s), s.scratch_src, coef_src)
    mul!(s.sy, s.Ut, s.scratch_src)
    for i =1:length(s.sy)
        s.sy[i]=s.sy[i]*s.Sinv[i]
    end
    mul!(s.y, s.V, s.sy)
    delinearize_coefficients!(dest(s), coef_dest, s.y)
end

"""
A RestrictionSolver is a solver used specifically for B splines,
taking advantage of the sparse structure.
"""
struct RestrictionSolver{ELT} <: FrameFun.FE_Solver{ELT}
    Aop::DictionaryOperator

    A::Matrix{ELT}
    # Mapping operator of dict to its boundary overlapping prolates
    BE::DictionaryOperator
    # Selection operator for the collocation points in the support of the boundary overlapping prolates
    GR::DictionaryOperator

    cutoff

    scratch_b::Array{ELT}
    scratch_b_linear::Vector{ELT}
    scratch_y1_native::Array{ELT}
    y1::Vector{ELT}
    function RestrictionSolver{ELT}(A::DictionaryOperator, BE::DictionaryOperator, GR::DictionaryOperator; cutoff=FrameFun.default_cutoff(A), options...) where {ELT}
        new(A, truncated_svd(matrix(GR*A*BE),cutoff), BE, GR, cutoff, zeros(ELT,size(dest(GR))...), zeros(ELT, length(dest(GR))), zeros(ELT, size(src(BE))...), zeros(ELT, length(src(BE))))
    end
end

src(t::RestrictionSolver) = dest(t.Aop)
dest(t::RestrictionSolver) = src(t.Aop)
inv(t::RestrictionSolver) = t.Aop

RestrictionSolver(A::DictionaryOperator, BE::DictionaryOperator, GR::DictionaryOperator; options...) =
    RestrictionSolver{eltype(A)}(A, BE, GR; options...)

function BasisFunctions.apply!(S::RestrictionSolver, x, b::Vector)
    apply!(S.GR, S.scratch_b, b)

    mul!(S.y1, S.A, S.scratch_b)
    # y1, rr = LAPACK.gelsy!(S.scratch_A,S.scratch_b, S.cutoff)
    apply!(S.BE, x, S.y1)
end

function truncated_svd(A, cutoff)
    USV = LAPACK.gesdd!('S',A)
    S = USV[2]
    maxind = findlast(S.>(maximum(S)*cutoff))
    Sinv = 1 ./ S[1:maxind]
    USV[3][1:maxind,:]'*(Sinv.*USV[1][:,1:maxind]')
end


# Function with equal functionality, but allocating memory
function truncatedsvd_solve(b::Vector, A::DictionaryOperator; cutoff = default_cutoff(A), R = 5, growth_factor = 2, verbose = false, smallcoefficients=false,smalltol=10,options...)
    finished=false
    ELT = eltype(A)
    R = min(R, size(A,2))
    random_matrix = map(ELT, rand(size(A,2), R))
    C = apply_multiple(A, random_matrix)
    c = cond(C)
    cold = cutoff
    while (c < 1/cutoff) && (R<size(A,2)) && (c>cold*10)
        verbose && println("c : $c\t cold : $cold\t cutoff : $cutoff")
        verbose && println("Solver truncated at R = ", R, " dof out of ",size(A,2))
        R0 = R
        R = min(round(Int,growth_factor*R),size(A,2))
        extra_random_matrix = map(ELT, rand(size(A,2), R-R0))
        Cextra = apply_multiple(A, extra_random_matrix)
        random_matrix = [random_matrix extra_random_matrix]
        # Extra condition: condition number has to increase by a significant amount each step, otherwise, possibly well conditioned.
        # cold = c
        C = [C Cextra]
        c = cond(C)
    end
    verbose && println("c : $c\t cold : $cold\t cutoff : $cutoff")
    verbose && println("Solver truncated at R = ", R, " dof out of ",size(A,2))


    y, rr = LAPACK.gelsy!(C,b[:],cutoff)
    x1 = random_matrix*y;
    reshape(x1, size(src(A))...)
end
