# lowranksolver.jl

"""
A Randomized Truncated SVD solver storing a decomposition of an operator A.
If A is MxN, and of rank R, the cost of constructing this solver is MR^2.

For tensor product operators it returns a decomposition of the linearized system
"""
struct TruncatedSvdSolver{ELT} <: FE_Solver{ELT}
    # Keep the original operator
    op          ::  AbstractOperator
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


    function TruncatedSvdSolver{ELT}(op::AbstractOperator; cutoff = default_cutoff(op), R = 5, growth_factor = 2, verbose = false, smallcoefficients=false,smalltol=10,options...) where ELT
        finished=false
        USV = ()
        R = min(R, size(op,2))
        random_matrix = map(ELT, rand(size(op,2), R))
        C = apply_multiple(op, random_matrix)
        c = cond(C)
# <<<<<<< HEAD
#         # TODO change things with cold back
#         # cold = 1
#         m = maximum(abs.(C))
#         while (c < 1/cutoff) && (R<size(op,2)) #&& (c/cold > 10)
# =======
        cold = cutoff
        m = maximum(abs.(C))
        while (c < 1/cutoff) && (R<size(op,2)) && (c>cold*10)
            verbose && println("c : $c\t cold : $cold\t cutoff : $cutoff")
            verbose && println("Solver truncated at R = ", R, " dof out of ",size(op,2))
            R0 = R
            R = min(round(Int,growth_factor*R),size(op,2))
            extra_random_matrix = map(ELT, rand(size(op,2), R-R0))
            Cextra = apply_multiple(op, extra_random_matrix)
            random_matrix = [random_matrix extra_random_matrix]
            # Extra condition: condition number has to increase by a significant amount each step, otherwise, possibly well conditioned.
            # cold = c
            C = [C Cextra]
            c = cond(C)

        end
        verbose && println("c : $c\t cold : $cold\t cutoff : $cutoff")
        verbose && println("Solver truncated at R = ", R, " dof out of ",size(op,2))
        USV = LAPACK.gesdd!('S',C)
        S = USV[2]
        maxind = findlast(S.>(maximum(S)*cutoff))
        Sinv = 1./S[1:maxind]
        y = zeros(ELT, size(USV[3],1))
        sy = zeros(ELT, maxind)
        Wsrc = Span(DiscreteSet{ELT}(size(random_matrix,2)))
        Wdest = src(op)
        W = MatrixOperator(Wsrc, Wdest, random_matrix)

        scratch_src = zeros(ELT, length(dest(op)))
        scratch_dest = zeros(ELT, length(src(op)))
        new(op, W, USV[1][:,1:maxind]',Sinv,USV[3][1:maxind,:]',y,sy, scratch_src, scratch_dest,smallcoefficients,smalltol)
    end
end

src(t::TruncatedSvdSolver) = dest(t.op)
dest(t::TruncatedSvdSolver) = src(t.op)
inv(t::TruncatedSvdSolver) = t.op

TruncatedSvdSolver(op::AbstractOperator; options...) =
    TruncatedSvdSolver{eltype(op)}(op::AbstractOperator; options...)

# We don't need to (de)linearize coefficients when they are already vectors
function apply!(s::TruncatedSvdSolver, coef_dest::Vector, coef_src::Vector)
    # Applying the truncated SVD to P*Rhs
    A_mul_B!(s.sy, s.Ut, coef_src)
    for i =1:length(s.sy)
        s.sy[i]=s.sy[i]*s.Sinv[i]
    end
    if s.smallcoefficients
        s.sy[abs.(s.sy).>s.smalltol*maximum(abs.(coef_src))]=0
    end
    A_mul_B!(s.y, s.V, s.sy)

    apply!(s.W, coef_dest, s.y)
end

function apply!(s::TruncatedSvdSolver, coef_dest, coef_src)
    linearize_coefficients!(src(s), s.scratch_src, coef_src)
    A_mul_B!(s.sy, s.Ut, s.scratch_src)
    for i =1:length(s.sy)
        s.sy[i]=s.sy[i]*s.Sinv[i]
    end
    if s.smallcoefficients
        s.sy[abs.(s.sy).>s.smalltol*maximum(abs.(coef_src))]=0
    end
    A_mul_B!(s.y, s.V, s.sy)
    apply!(s.W, s.scratch_dest, s.y)
    delinearize_coefficients!(dest(s), coef_dest, s.scratch_dest)
end


"""
A Truncated SVD solver storing a decomposition of an operator A.
The cost is equal to the computation of the SVD of A.

For tensor product operators it returns a decomposition of the linearized system
"""
struct ExactTruncatedSvdSolver{ELT} <: FE_Solver{ELT}
    # Keep the original operator
    op          ::  AbstractOperator
    # Decomposition
    Ut          ::  Array{ELT,2}
    Sinv        ::  Array{ELT,1}
    V          ::  Array{ELT,2}
    # For storing intermediate results when applying
    y           ::  Array{ELT,1}
    sy          ::  Array{ELT,1}
    scratch_src ::  Array{ELT,1}

    function ExactTruncatedSvdSolver{ELT}(op::AbstractOperator; cutoff = default_cutoff(op), verbose = false,options...) where ELT
        C = matrix(op)
        m = maximum(abs.(C))
        USV = LAPACK.gesdd!('S',C)
        S = USV[2]
        maxind = findlast(S.>cutoff*m)
        verbose && println("Solver truncated at singular value equal to ", S[maxind], " cutoff was ", cutoff, " and norm of A ", m)
        Sinv = 1./S[1:maxind]
        y = zeros(ELT, size(USV[3],1))
        sy = zeros(ELT, maxind)

        scratch_src = zeros(ELT, length(dest(op)))
        new(op, USV[1][:,1:maxind]',Sinv,USV[3][1:maxind,:]',y,sy, scratch_src)
    end
end

src(t::ExactTruncatedSvdSolver) = dest(t.op)
dest(t::ExactTruncatedSvdSolver) = src(t.op)
inv(t::ExactTruncatedSvdSolver) = t.op

ExactTruncatedSvdSolver(op::AbstractOperator; options...) =
    ExactTruncatedSvdSolver{eltype(op)}(op::AbstractOperator; options...)

# We don't need to (de)linearize coefficients when they are already vectors
function apply!(s::ExactTruncatedSvdSolver, coef_dest::Vector, coef_src::Vector)
    # Applying the truncated SVD to P*Rhs
    A_mul_B!(s.sy, s.Ut, coef_src)
    for i =1:length(s.sy)
        s.sy[i]=s.sy[i]*s.Sinv[i]
    end
    A_mul_B!(coef_dest, s.V, s.sy)
end

function apply!(s::ExactTruncatedSvdSolver, coef_dest, coef_src)
    linearize_coefficients!(src(s), s.scratch_src, coef_src)
    A_mul_B!(s.sy, s.Ut, s.scratch_src)
    for i =1:length(s.sy)
        s.sy[i]=s.sy[i]*s.Sinv[i]
    end
    A_mul_B!(s.y, s.V, s.sy)
    delinearize_coefficients!(dest(s), coef_dest, s.y)
end
