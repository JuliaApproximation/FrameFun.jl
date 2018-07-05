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
        Sinv = 1./S[1:maxind]
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
    A_mul_B!(sy, Ut, coef_src)
    for i in 1:length(sy)
        sy[i] = sy[i]*Sinv[i]
    end
    if smallcoefficients
        M = maximum(abs.(coef_src))
        sy[abs.(sy) .> smalltol*M] = 0
    end
    A_mul_B!(y, V, sy)
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
        Sinv = 1./S[1:maxind]
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

    A_mul_B!(S.y1, S.A, S.scratch_b)
    # y1, rr = LAPACK.gelsy!(S.scratch_A,S.scratch_b, S.cutoff)
    apply!(S.BE, x, S.y1)
end

function truncated_svd(A, cutoff)
    USV = LAPACK.gesdd!('S',A)
    S = USV[2]
    maxind = findlast(S.>(maximum(S)*cutoff))
    Sinv = 1./S[1:maxind]
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

# Function with equal functionality, but allocating memory
restriction_solve(b::Vector, A::DictionaryOperator, BE::IndexExtensionOperator,
        GR::IndexRestrictionOperator; cutoff=FrameFun.default_cutoff(A), options...) =
    BE*LAPACK.gelsy!(matrix(GR*A*BE),GR*b,cutoff)[1]

function decomposition_solve(b::Vector, A::DictionaryOperator, cart_indices, classified_indices; cutoff=default_cutoff(A), verbose=false, info=false, options...)
    bins = unique(classified_indices)
    # assign sequence numbers to each of the bins.
    seq_nr = assign_sequence_nro(bins)

    primal = basis(src(A))
    omega = grid(dest(A))
    g = supergrid(omega)
    x1 = zeros(primal)
    b_ =  copy(b)
    t = similar(x1)

    for d in 1:1<<length(bins[1])
        # first solve all parts with a low sequence number
        bpart = bins[find(seq_nr.==d)]
        # Each of the parts with an equal sequence number can be solved independently
        verbose && println("$(length(bpart)) parts in $(d)th flow")
        for i in 1:size(bpart,1)
            mo = cart_indices[find(classified_indices.==bpart[i,:])]
            xx, yy = FrameFun._azselection_restriction_operators(primal, g, omega, mo)
            verbose && println("\t$(i)\t has size ($(size(yy,1)),$(size(xx,1)))")
            op = yy*A*xx'
            a = matrix(op)
            y = LAPACK.gelsy!(a, yy*b_, cutoff)[1]
            # x1 = x1 + xx'*y
            apply!(xx', t, y)
            x1 .+= t
        end
        # Remove the solved part after all parts with equal seq number are dealt with.
        if d!=1<<length(bins[1])
            # b_ = b-A*x1
            apply!(A, b_, x1)
            b_ .= b .- b_
        end
    end
    x1
end

function decomposition_info(b::Vector, A::DictionaryOperator, cart_indices, classified_indices; cutoff=default_cutoff(A), info=false, options...)
    bins = unique(classified_indices)
    seq_nr = assign_sequence_nro(bins)
    primal = basis(src(A))
    omega = grid(dest(A))
    g = supergrid(omega)
    r = zeros(Int, length(bins), 2)
    ii = 1
    for d in 1:1<<length(bins[1])
        bpart = bins[find(seq_nr.==d)]
        println("$(length(bpart)) parts in $(d)th flow")
        for i in 1:size(bpart,1)
            mo = cart_indices[find(classified_indices.==bpart[i,:])]
            xx, yy = FrameFun._azselection_restriction_operators(primal, g, omega, mo)
            op = yy*A*xx'
            println("\t$(i)\t has size ($(size(yy,1)),$(size(xx,1)))")
            r[ii,1] = size(yy,1)
            r[ii,2] = size(xx,1)
            ii += 1
        end
    end
    r
end

using Plots
function decomposition_plot(A::DictionaryOperator, cart_indices, classified_indices)
    bins = unique(classified_indices)
    seq_nr = assign_sequence_nro(bins)
    # clrs = [:blue,:red,:green,:black]
    clrs = Plots.colormap("blues",1<<length(bins[1]))
    seq_nr = assign_sequence_nro(bins)
    primal = basis(src(A))
    plot()
    for d in 1:1<<length(bins[1])
        bpart = bins[find(seq_nr.==d)]
        for i in 1:size(bpart,1)
            mo = cart_indices[find(classified_indices.==bpart[i,:])]
            m = falses(size(primal))
            m[mo] = true
            scatter!(m,c=clrs[d])
        end
    end
    scatter!()
end

# Function with equal functionality, but allocating memory
restriction_info(b::Vector, A::DictionaryOperator, BE::IndexExtensionOperator,
        GR::IndexRestrictionOperator; cutoff=FrameFun.default_cutoff(A), options...) =
    (println("Selection has size ($(size(GR,1)),$(size(BE,2)))"); [size(GR,1),size(BE,2)]')
