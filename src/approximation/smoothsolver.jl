# smoothsolver.jl

"""
A fast FE solver based on a low-rank approximation of the plunge region. The plunge region
is isolated using a projection operator. This algorithm contains an extra smoothing step

"""
struct AZSmoothSolver{ELT} <: FE_Solver{ELT}
    TS          :: DictionaryOperator
    A           ::  DictionaryOperator
    Zt          ::  DictionaryOperator
    plunge_op   ::  DictionaryOperator    # store the operator because it allocates memory
    b           ::  Array{ELT,1}
    blinear     ::  Array{ELT,1}
    syv         ::  Array{ELT,1}
    x1          ::  Array{ELT}
    x2          ::  Array{ELT}
    x3          ::  Array{ELT}
    Q           ::  Array{ELT,2}
    D           ::  DictionaryOperator
    AD          ::  DictionaryOperator

    function AZSmoothSolver{ELT}(A::DictionaryOperator, Zt::DictionaryOperator, D::DictionaryOperator; cutoff = default_cutoff(A), cutoffv=sqrt(cutoff), R = estimate_plunge_rank(A), verbose=false,  options...) where ELT
        plunge_op = plunge_operator(A, Zt)
        # Create Random matrices
        TS1 = TruncatedSvdSolver(plunge_op*A; cutoff = cutoff, verbose=verbose,R=R,options...)
        TS2 = TruncatedSvdSolver(Zt*plunge_op; cutoff = cutoffv, verbose=verbose, R=R, options...)
        AD = inv(D)
        ADV = (TS2.Ut)'.*diagonal(AD)
        # Orthogonal basis for D^(-1)V_mid
        Q, R = qr(ADV)
        # Pre-allocation
        b = zeros(size(dest(plunge_op)))
        blinear = zeros(ELT, length(dest(A)))
        x1 = zeros(size(src(A)))
        x2 = zeros(size(src(A)))
        x3 = zeros(size(src(A)))
        syv = zeros(size(TS2.Ut,1))
        new(TS1, A, Zt, plunge_op, b, blinear, syv, x1, x2, x3, Q, D, AD)
    end
end

function AZSmoothSolver(A::DictionaryOperator, Zt::DictionaryOperator; options...)
    D = IdxnScalingOperator(src(A); options...)
    AZSmoothSolver{eltype(A)}(A, Zt, D; options...)
end

AZSmoothSolver(A::DictionaryOperator; scaling=nothing, options...) =
        AZSmoothSolver{eltype(A)}(A, 1/scaling*A'; options...)

apply!(s::AZSmoothSolver, coef_dest, coef_src) =
    _apply!(s, coef_dest, coef_src, s.TS, s.A, s.Zt, s.plunge_op, s.b, s.blinear, s.syv, s.x1, s.x2, s.x3, s.Q, s.D, s.AD)

function _apply!(s::AZSmoothSolver, coef_dest, coef_src, TS, A, Zt, plunge_op, b, blinear, syv, x1, x2, x3, Q, D, AD)
    # Apply plunge operator to right hand side
    apply!(plunge_op, b, coef_src)
    BasisFunctions.linearize_coefficients!(dest(A), blinear, b)
    # x2 is solving for the middle singular values
    apply!(TS, x2, blinear)
    # Project Dx2 onto the orthogonal complement of D^(-1)V_mid
    apply!(D, x1, x2)
    ## apply!(MantrixOperator(s.Q'),s.syv,s.x1)
    ## apply!(MatrixOperator(s.Q),s.x3,s.syv)
    ## for i = 1:length(s.x1)
    ##     s.x1[i] = s.x1[i]- s.x3[i]
    ## end
    ## apply!(AD,s.x3,s.x1)
    ## for i = 1:length(s.x2)
    ##     s.x2[i] = s.x2[i] - s.x3[i]
## end
    apply!(MatrixOperator(Q'), syv, x1)
    apply!(MatrixOperator(Q), x3, syv)
    apply!(AD, x2, x3)

    # post smoothing step
    apply!(A, b, x2)
    apply!(Zt, x1, coef_src-b)
    for i in 1:length(coef_dest)
        coef_dest[i] = x1[i] + x2[i]
    end
    coef_dest
end

# For Dictionary's that have a DC component
dc_index(b::ChebyshevBasis) = 1
dc_index(b::FourierBasis) = 1

# An index scaling operator, used to generate weights for the polynomial scaling algorithm.
struct IdxnScalingOperator{ELT} <: DictionaryOperator{ELT}
    src     ::  Dictionary
    order   ::  Int
    scale   ::  Function
end

IdxnScalingOperator(dict::Dictionary; order=1, scale = default_scaling_function) =
    IdxnScalingOperator{BasisFunctions.op_eltype(dict,dict)}(dict, order, scale)

dest(op::IdxnScalingOperator) = src(op)

default_scaling_function(i) = 10.0^-4+(abs(i))+abs(i)^2+abs(i)^3
default_scaling_function(i,j) = 1+(abs(i)^2+abs(j)^2)

is_inplace(::IdxnScalingOperator) = true
is_diagonal(::IdxnScalingOperator) = true

ctranspose(op::IdxnScalingOperator) = DiagonalOperator(src(op), conj(diagonal(op)))

function apply_inplace!(op::IdxnScalingOperator, dest::Dictionary1d, srcdict, coef_srcdest)
    ELT = eltype(op)
    for i in eachindex(dest)
        coef_srcdest[i] *= op.scale(ELT(BasisFunctions.value(native_index(dest,i))))^op.order
    end
    coef_srcdest
end

function apply_inplace!(op::IdxnScalingOperator, dest::Dictionary2d, srcdict, coef_srcdest)
    ELT = eltype(op)
    for i in eachindex(coef_srcdest)
        ni = recursive_native_index(dest,i)
        coef_srcdest[i]*=op.scale(ELT(BasisFunctions.value(ni[1])),ELT(BasisFunctions.value(ni[2])))^op.order
    end
    coef_srcdest
end
inv(op::IdxnScalingOperator) = IdxnScalingOperator(op.src, order=op.order*-1, scale=op.scale)
