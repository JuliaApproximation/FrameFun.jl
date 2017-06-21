# lowranksolver.jl

"""
A Randomized Truncated SVD solver storing a decomposition of an operator A.
If A is MxN, and of rank R, the cost of constructing this solver is MR^2.

For tensor product operators it returns a decomposition of the linearized system
"""
immutable TruncatedSvdSolver{ELT} <: AbstractOperator{ELT}
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

    function TruncatedSvdSolver(op::AbstractOperator; cutoff = default_cutoff(problem), R = 5, growth_factor = sqrt(2), verbose = false, smallcoefficients=false,options...)
        finished=false
        USV = ()
        R = min(R, size(op,2))
        random_matrix = map(ELT, rand(size(op,2), R))
        C = apply_multiple(op, random_matrix)
        c = cond(C)
        m = norm(C)
        while (c < 1/cutoff) && (R<size(op,2))
            R0 = R
            R = min(round(Int,growth_factor*R),size(op,2))
            verbose && println("Solver truncated at R = ", R, " dof out of ",size(op,2))
            verbose && println(c," ",m," ",cutoff)            
            extra_random_matrix = map(ELT, rand(size(op,2), R-R0))
            Cextra = apply_multiple(op, extra_random_matrix)
            random_matrix = [random_matrix extra_random_matrix]
            
            C = [C Cextra]
            c = cond(C)
            m = norm(C)
        end
        verbose && println("Solver truncated at R = ", R, " dof out of ",size(op,2))
            verbose && println(c," ",m," ",cutoff)            
        USV = LAPACK.gesdd!('S',C)
        S = USV[2]
        maxind = findlast(S.>cutoff*m)
        Sinv = 1./S[1:maxind]
        y = zeros(ELT, size(USV[3],1))
        sy = zeros(ELT, maxind)
        Wsrc = ELT <: Complex ? Cn{ELT}(size(random_matrix,2)) : Rn{ELT}(size(random_matrix,2))
        Wdest = src(op)
        W = MatrixOperator(Wsrc, Wdest, random_matrix)

        scratch_src = zeros(ELT, length(dest(op)))
        scratch_dest = zeros(ELT, length(src(op)))
        new(op, W, USV[1][:,1:maxind]',Sinv,USV[3][1:maxind,:]',y,sy, scratch_src, scratch_dest,smallcoefficients)
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
        s.sy[abs(s.sy).>100*maximum(abs(coef_src))]=0
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
        s.sy[abs(s.sy).>100*maximum(abs(coef_src))]=0
    end
    A_mul_B!(s.y, s.V, s.sy)
    apply!(s.W, s.scratch_dest, s.y)
    delinearize_coefficients!(dest(s), coef_dest, s.scratch_dest)
end


"""
A Double Randomized Truncated SVD solver storing a decomposition of an operator A.
If A is MxN, and of rank R, the cost of constructing this solver is MR^2.

For tensor product operators it returns a decomposition of the linearized system
"""
immutable DoubleTruncatedSvdSolver{ELT} <: AbstractOperator{ELT}
    # Keep the original operator
    op          ::  AbstractOperator
    # Random matrix
    W           ::  MultiplicationOperator
    Wb          ::  MultiplicationOperator
    # Decomposition
    Ut          ::  Array{ELT,2}
    VS          ::  Array{ELT,2}
    # For storing intermediate results when applying
    y           ::  Array{ELT,1}
    yb          ::  Array{ELT,1}
    sy          ::  Array{ELT,1}
    scratch_src ::  Array{ELT,1}
    scratch_dest::  Array{ELT,1}

    function DoubleTruncatedSvdSolver(op::AbstractOperator; cutoff = default_cutoff(problem), R = 5, growth_factor = sqrt(2), verbose = false, options...)
        finished=false
        USV = ()
        R = min(R, size(op,2))
        random_matrix = map(ELT, rand(size(op,2), R))
        random_matrixb = map(ELT, rand(2*R,size(op,1)))
        C = apply_multiple(op, random_matrix)
        D = random_matrixb * C
        c = cond(C)
        m = maximum(abs(C))
        while (c < 1/cutoff) && (R<size(op,2))
            verbose && println("DoubleSolver truncated at R = ", R, " dof out of ",size(op,2))
            R0 = R
            R = min(round(Int,sqrt(2)*R),size(op,2))
            extra_random_matrix = map(ELT, rand(size(op,2), R-R0))
            random_matrixb = map(ELT, rand(growth_factor*R,size(op,1)))
            Cextra = apply_multiple(op, extra_random_matrix)
            random_matrix = [random_matrix extra_random_matrix]
            C = [C Cextra]
            D = random_matrixb * C
            c = cond(D)
            m = maximum(abs(D))
        end
        verbose && println("DoubleSolver truncated at R = ", R, " dof out of ",size(op,2))

        USV = LAPACK.gesdd!('S',D)
        S = USV[2]
        maxind = findlast(S.>cutoff*m)
        Sinv = 1./S[1:maxind]
        y = zeros(ELT, size(USV[3],1))
        yb = zeros(ELT, size(random_matrixb,1))
        sy = zeros(ELT, maxind)
        Wsrc = ELT <: Complex ? Cn{ELT}(size(random_matrix,2)) : Rn{ELT}(size(random_matrix,2))
        Wdest = src(op)
        W = MatrixOperator(Wsrc, Wdest, random_matrix)
        Wsrcb = dest(op)
        Wdestb = ELT <: Complex ? Cn{ELT}(size(random_matrixb,1)) : Rn{ELT}(size(random_matrixb,1))
        Wb = MatrixOperator(Wsrcb, Wdestb, random_matrixb)
        scratch_src = zeros(ELT, length(dest(op)))
        scratch_dest = zeros(ELT, length(src(op)))
        new(op, W, Wb, USV[1][:,1:maxind]',USV[3][1:maxind,:]'*diagm(Sinv[:]),y,yb,sy, scratch_src, scratch_dest)
    end
end

src(t::DoubleTruncatedSvdSolver) = dest(t.op)
dest(t::DoubleTruncatedSvdSolver) = src(t.op)
inv(t::DoubleTruncatedSvdSolver) = t.op

DoubleTruncatedSvdSolver(op::AbstractOperator; options...) =
    DoubleTruncatedSvdSolver{eltype(op)}(op::AbstractOperator; options...)

# We don't need to (de)linearize coefficients when they are already vectors
function apply!(s::DoubleTruncatedSvdSolver, coef_dest::Vector, coef_src::Vector)
    # Applying the truncated SVD to P*Rhs
    apply!(s.Wb, s.yb, coef_src)
    A_mul_B!(s.sy, s.Ut, s.yb)
    A_mul_B!(s.y, s.VS, s.sy)
    apply!(s.W, coef_dest, s.y)
end

function apply!(s::DoubleTruncatedSvdSolver, coef_dest, coef_src)
    linearize_coefficients!(src(s), s.scratch_src, coef_src)
    apply!(s.Wb, s.yb, s.scratch_src)
    A_mul_B!(s.sy, s.Ut, s.yb)
    A_mul_B!(s.y, s.VS, s.sy)
    apply!(s.W, s.scratch_dest, s.y)
    delinearize_coefficients!(dest(s), coef_dest, s.scratch_dest)
end
