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
        Wsrc = DiscreteVectorSet{ELT}(size(random_matrix,2))
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

# function BasisFunctions.apply!(S::RestrictionSolver, x, b)
#     copy!(S.scratch_A, S.A)
#     apply!(S.GR, S.scratch_b, b)
#     BasisFunctions.linearize_coefficients!(dest(S.GR), S.scratch_b_linear, S.scratch_b)
#     y1, rr = LAPACK.gelsy!(S.scratch_A,S.scratch_b_linear, S.cutoff)
#     BasisFunctions.delinearize_coefficients!(src(S.BE), S.scratch_y1_native, y1)
#     apply!(S.BE, x, S.scratch_y1_native)
# end


"""
A DivideAndConquerSolver is a solver used specifically for B splines,
taking advantage of the sparse structure.
"""
struct DivideAndConquerSolver{ELT} <: FrameFun.FE_Solver{ELT}
    A::DictionaryOperator{ELT}
    # Small problem between A1 and A2
    A0::Matrix{ELT}
    A1::Matrix{ELT}
    A2::Matrix{ELT}
    # Mapping operator of dict to its boundary overlapping prolates
    BE0::DictionaryOperator
    BE1::DictionaryOperator
    BE2::DictionaryOperator
    # Selection operator for the collocation points in the support of the boundary overlapping prolates
    GR0::DictionaryOperator
    GR1::DictionaryOperator
    GR2::DictionaryOperator

    cutoff

    scratch_b0::Array{ELT}
    scratch_b1::Array{ELT}
    scratch_b2::Array{ELT}

    # x0::Array{ELT}
    x1::Array{ELT}
    # x2::Array{ELT}

    y0::Array{ELT}
    y1::Array{ELT}
    y2::Array{ELT}

    bnew::Array{ELT}

    function DivideAndConquerSolver{ELT}(
        A::DictionaryOperator,
        BE0::DictionaryOperator, BE1::DictionaryOperator, BE2::DictionaryOperator,
        GR0::DictionaryOperator, GR1::DictionaryOperator, GR2::DictionaryOperator;
        cutoff=FrameFun.default_cutoff(A), options...) where {ELT}

        new(A, truncated_svd(matrix(GR0*A*BE0), cutoff), truncated_svd(matrix(GR1*A*BE1), cutoff), truncated_svd(matrix(GR2*A*BE2), cutoff),
        BE0, BE1, BE2,
        GR0, GR1, GR2,
        cutoff,
        zeros(ELT, length(dest(GR0))), zeros(ELT, length(dest(GR1))), zeros(ELT, length(dest(GR2))),
        zeros(ELT, size(dest(BE0))...), #zeros(ELT, size(dest(BE1))...), zeros(ELT, size(dest(BE1))...),
        zeros(ELT, length(src(BE0))), zeros(ELT, length(src(BE1))), zeros(ELT, length(src(BE2))),
        zeros(ELT, size(dest(A))...))
    end
end

function truncated_svd(A, cutoff)
    USV = LAPACK.gesdd!('S',A)
    S = USV[2]
    maxind = findlast(S.>(maximum(S)*cutoff))
    Sinv = 1./S[1:maxind]
    USV[3][1:maxind,:]'*(Sinv.*USV[1][:,1:maxind]')
end

src(t::DivideAndConquerSolver) = dest(t.A)
dest(t::DivideAndConquerSolver) = src(t.A)
inv(t::DivideAndConquerSolver) = t.A

DivideAndConquerSolver( A::DictionaryOperator,
                        BE0::DictionaryOperator, BE1::DictionaryOperator, BE2::DictionaryOperator,
                        GR0::DictionaryOperator, GR1::DictionaryOperator, GR2::DictionaryOperator
                        ; options...) =
    DivideAndConquerSolver{eltype(A)}(A, BE0, BE1, BE2, GR0, GR1, GR2; options...)

function BasisFunctions.apply!(S::DivideAndConquerSolver, x, b::Vector)
    apply!(S.GR0, S.scratch_b0, b)

    A_mul_B!(S.y0, S.A0, S.scratch_b0)
    apply!(S.BE0, x, S.y0)

    apply!(S.A, S.bnew, x)
    for i in 1:length(S.bnew)
        S.bnew[i] = b[i] - S.bnew[i]
    end

    apply!(S.GR1, S.scratch_b1, S.bnew)
    A_mul_B!(S.y1, S.A1, S.scratch_b1)
    apply!(S.BE1, S.x1, S.y1)

    for i in 1:length(S.x1)
        x[i] = x[i] + S.x1[i]
    end

    apply!(S.GR2, S.scratch_b2, S.bnew)
    A_mul_B!(S.y2, S.A2, S.scratch_b2)
    apply!(S.BE2, S.x1, S.y2)

    for i in 1:length(S.x1)
        x[i] = x[i] + S.x1[i]
    end
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

# restriction_solve(A, BE, GR, b; cutoff=FrameFun.default_cutoff(A), options...) =
#     BE*reshape(LAPACK.gelsy!(matrix(GR*A*BE),(GR*b)[:],cutoff)[1], size(src(BE)))


function divideandconqer_solve(b::Vector, A::DictionaryOperator, BE0::IndexExtensionOperator, BE1::IndexExtensionOperator, BE2::IndexExtensionOperator,
        GR0::IndexRestrictionOperator, GR1::IndexRestrictionOperator, GR2::IndexRestrictionOperator;  cutoff=FrameFun.default_cutoff(A), options...)
    x0 = BE0*LAPACK.gelsy!(matrix(GR0*A*BE0), GR0*b, cutoff)[1]
    Lb = length(b)::Int
    Lx = length(x0)::Int

    bnew = zeros(size(dest(A)))

    # bnew = A*x0
    apply!(A, bnew, x0)
    bnew .= b .- bnew
    # for i in 1:Lb
    #     bnew[i] = b[i] - bnew[i]
    # end
    x1 = BE1*LAPACK.gelsy!(matrix(GR1*A*BE1), GR1*bnew, cutoff)[1]
    x0 .= x0 .+ x1
    # for i in 1:Lx
    #     x0[i] = x0[i] + x1[i]
    # end

    x1 = BE2*LAPACK.gelsy!(matrix(GR2*A*BE2), GR2*bnew, cutoff)[1]
    x0 .= x0 .+ x1
    # for i in 1:Lx
    #     x0[i] = x0[i] + x1[i]
    # end

    x0
end

function divideandconqerN_solve(b::Vector, A, ops::DictionaryOperator...; divideandconqerN_op_lengths=nothing, options...)
    A0 = Array{CompositeOperator}(divideandconqerN_op_lengths[1])
    G0 = Array{IndexRestrictionOperator}(divideandconqerN_op_lengths[2])
    A1 = Array{CompositeOperator}(divideandconqerN_op_lengths[3])
    G1 = Array{IndexRestrictionOperator}(divideandconqerN_op_lengths[4])
    i = 1
    for ARRAY in (A0, G0, A1, G1)
        for j in 1:length(ARRAY)
            ARRAY[j] = ops[i]
            i += 1
        end
    end
    if length(divideandconqerN_op_lengths) == 4
        divideandconqerN_solve(b, A, A0, G0, A1, G1; options...)
    else
        A2 = Array{CompositeOperator}(divideandconqerN_op_lengths[1])
        G2 = Array{IndexRestrictionOperator}(divideandconqerN_op_lengths[2])
        A3 = Array{CompositeOperator}(divideandconqerN_op_lengths[3])
        G3 = Array{IndexRestrictionOperator}(divideandconqerN_op_lengths[4])
        i = 1
        for ARRAY in (A2, G2, A3, G3)
            for j in 1:length(ARRAY)
                ARRAY[j] = ops[i]
                i += 1
            end
        end
        divideandconqerN_solve(b, A, A0, G0, A1, G1, A2, G2, A3, G3; options...)
    end
end

function divideandconqerN_solve(b::Vector, A, A0::Vector{OP1}, GR0::Vector{OP2},
        A1::Vector{OP3}, GR1::Vector{OP4};
        cutoff=FrameFun.default_cutoff(A), options...) where {OP1<:DictionaryOperator, OP2<:IndexRestrictionOperator, OP3<:DictionaryOperator, OP4<:IndexRestrictionOperator}
    omega = BasisFunctions.grid(dest(A))
    gamma = supergrid(omega)
    omega_restriction = restriction_operator(gridbasis(gamma), gridbasis(omega))
    b_ext = omega_restriction'b
    x0 = zeros(src(A))
    for (gr,a) in zip(GR0,A0)
        y0 = LAPACK.gelsy!(matrix(a), gr*b_ext, cutoff)[1]
        x0[a.operators[1].subindices] .+= y0
    end

    bnew = similar(b)
    apply!(A, bnew, x0)

    for (i_i,i) in enumerate(subindices(omega))
        b_ext[i] -= bnew[i_i]
    end

    for (gr,a) in zip(GR1,A1)
        y0 = LAPACK.gelsy!(matrix(a), gr*b_ext, cutoff)[1]
        x0[a.operators[1].subindices] .+= y0
    end
    x0
end

function divideandconqerN_solve(b::Vector, A,
        A0::Vector{OP1}, GR0::Vector{OP2},
        A1::Vector{OP3}, GR1::Vector{OP4},
        A2::Vector{OP5}, GR2::Vector{OP6},
        A3::Vector{OP7}, GR3::Vector{OP8};
        cutoff=FrameFun.default_cutoff(A), options...) where {OP1<:DictionaryOperator, OP2<:IndexRestrictionOperator, OP3<:DictionaryOperator, OP4<:IndexRestrictionOperator, OP5<:DictionaryOperator, OP6<:IndexRestrictionOperator, OP7<:DictionaryOperator, OP8<:IndexRestrictionOperator}
    omega = BasisFunctions.grid(dest(A))
    gamma = supergrid(omega)
    omega_restriction = restriction_operator(gridbasis(gamma), gridbasis(omega))
    b_ext = omega_restriction'b
    b_ext_copy = copy(b_ext)
    x0 = zeros(src(A))

    bnew = similar(b)

    # split domain of split domain
    for (gr,a) in zip(GR2,A2)
        y0 = LAPACK.gelsy!(matrix(a), gr*b_ext, cutoff)[1]
        x0[a.operators[1].subindices] .+= y0
    end
    apply!(A, bnew, x0)
    for (i_i,i) in enumerate(subindices(omega))
        b_ext[i] -= bnew[i_i]
    end

    # pieces of split domain
    for (gr,a) in zip(GR3,A3)
        y0 = LAPACK.gelsy!(matrix(a), gr*b_ext, cutoff)[1]
        x0[a.operators[1].subindices] .+= y0
    end
    apply!(A, bnew, x0)
    for (i_i,i) in enumerate(subindices(omega))
        b_ext[i] = b_ext_copy[i] - bnew[i_i]
    end

    # pieces of domain
    for (gr,a) in zip(GR1,A1)
        y0 = LAPACK.gelsy!(matrix(a), gr*b_ext, cutoff)[1]
        x0[a.operators[1].subindices] .+= y0
    end
    x0
end

struct DomainDecompositionSolver{ELT} <: FE_Solver{ELT}
    tree    :: AbstractDomainDecompositionNode
    basis   :: Dictionary
    gamma   :: AbstractGrid
    omega   :: AbstractGrid
    domain  :: Domain
    A       :: DictionaryOperator
end

DomainDecompositionSolver(fplatform::BasisFunctions.GenericPlatform, i; options...) =
    DomainDecompositionSolver(primal(fplatform.super_platform, i), grid(sampler(fplatform.super_platform, i)),
        grid(sampler(fplatform, i)), domain(primal(fplatform, i)), A(fplatform, i); options...)


DomainDecompositionSolver(basis, gamma, omega, domain, A; options...) =
    DomainDecompositionSolver{coeftype(basis)}(create_tree(basis, gamma, omega, domain; options...), basis, gamma, omega, domain, A)

src(s::DomainDecompositionSolver) = extensionframe(s.basis, s.domain)
dest(s::DomainDecompositionSolver{ELT}) where {ELT} = gridbasis(s.omega, ELT)

BasisFunctions.apply(s::DomainDecompositionSolver, src) = domaindecomposition_solve(src, s.A, s)

domaindecomposition_solve(b::Vector, A::DictionaryOperator, s::DomainDecompositionSolver; options...) =
    solve(b, A, s.tree, s.basis, s.gamma, s.omega; options...)
