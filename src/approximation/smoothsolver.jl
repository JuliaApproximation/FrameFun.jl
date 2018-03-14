# smoothsolver.jl

"""
A fast FE solver based on a low-rank approximation of the plunge region. The plunge region
is isolated using a projection operator. This algorithm contains an extra smoothing step

"""
struct AZSmoothSolver{ELT} <: FE_Solver{ELT}
    TS :: AbstractOperator
    A     ::  AbstractOperator
    Zt    ::  AbstractOperator
    plunge_op   ::  AbstractOperator    # store the operator because it allocates memory
    b     ::  Array{ELT,1}
    blinear     ::  Array{ELT,1}
    syv         ::  Array{ELT,1}
    x2          ::  Array{ELT}
    x1          ::  Array{ELT}
    x3          ::  Array{ELT}
    Q          ::  Array{ELT,2}
    D           ::  AbstractOperator

    function AZSmoothSolver{ELT}(A::AbstractOperator, Zt::AbstractOperator; cutoff = default_cutoff(A), cutoffv=sqrt(cutoff), R = estimate_plunge_rank(A), verbose=false,  options...) where ELT
        plunge_op = plunge_operator(A, Zt)
        # Create Random matrices
        TS1 = TruncatedSvdSolver(plunge_op*A; cutoff = cutoff, verbose=verbose,R=R,options...)
        TS2 = TruncatedSvdSolver(Zt*plunge_op; cutoff = cutoffv, verbose=verbose, R=R, options...)
        # D = Sobolev operator
        D = IdxnScalingOperator(src(A); options...)
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
        new(TS1,A, Zt, plunge_op, b,blinear,syv,x1,x2,x3,Q, D)
    end
end

AZSmoothSolver(A::AbstractOperator, Zt::AbstractOperator; options...) =
    AZSmoothSolver{eltype(A)}(A, Zt; options...)
    
AZSmoothSolver(A::AbstractOperator; scaling=nothing, options...) =
        AZSmoothSolver{eltype(A)}(A, 1/scaling*A'; options...)

function apply!(s::AZSmoothSolver, destarg, src, coef_dest, coef_src)
    A = s.A
    At = s.Zt
    P = s.plunge_op
    # Apply plunge operator to right hand side
    apply!(P,s.b, coef_src)
    BasisFunctions.linearize_coefficients!(dest(A), s.blinear, s.b)
    # x2 is solving for the middle singular values
    apply!(s.TS,s.x2,s.blinear)
    # smoothing x2 step
    D = s.D
    AD = inv(D)
    # Project Dx2 onto the orthogonal complement of D^(-1)V_mid
    apply!(D,s.x1,s.x2)
    ## apply!(MantrixOperator(s.Q'),s.syv,s.x1)
    ## apply!(MatrixOperator(s.Q),s.x3,s.syv)
    ## for i = 1:length(s.x1)
    ##     s.x1[i] = s.x1[i]- s.x3[i]
    ## end
    ## apply!(AD,s.x3,s.x1)
    ## for i = 1:length(s.x2)
    ##     s.x2[i] = s.x2[i] - s.x3[i]
## end
    apply!(MatrixOperator(s.Q'),s.syv,s.x1)
    apply!(MatrixOperator(s.Q),s.x3,s.syv)
    apply!(AD,s.x2,s.x3)

    # post smoothing step
    apply!(A, s.b, s.x2)
    apply!(At, s.x1, coef_src-s.b)
    for i = 1:length(coef_dest)
        coef_dest[i] = s.x1[i] + s.x2[i]
    end

end

# For FunctionSets that have a DC component
dc_index(b::ChebyshevBasis) = 1
dc_index(b::FourierBasis) = 1

# An index scaling operator, used to generate weights for the polynomial scaling algorithm.
struct IdxnScalingOperator{ELT} <: AbstractOperator{ELT}
    src     ::  Span
    order   ::  Int
    scale   ::  Function
end

IdxnScalingOperator(span::Span; order=1, scale = default_scaling_function) =
    IdxnScalingOperator{BasisFunctions.op_eltype(span,span)}(span, order, scale)

dest(op::IdxnScalingOperator) = src(op)

default_scaling_function(i) = 10.0^-4+(abs(i))+abs(i)^2+abs(i)^3
default_scaling_function(i,j) = 1+(abs(i)^2+abs(j)^2)

is_inplace(::IdxnScalingOperator) = true
is_diagonal(::IdxnScalingOperator) = true

ctranspose(op::IdxnScalingOperator) = DiagonalOperator(src(op), conj(diagonal(op)))
function apply_inplace!(op::IdxnScalingOperator, destspan, srcspan, coef_srcdest)
    dest = set(destspan)
    ELT = eltype(op)
    for i in eachindex(dest)
        coef_srcdest[i] *= op.scale(ELT(BasisFunctions.index(native_index(dest,i))))^op.order
    end
    coef_srcdest
end

function apply_inplace!(op::IdxnScalingOperator, destspan::Span{A,SET}, srcspan, coef_srcdest) where {A,TS1,TS2, SET<:TensorProductSet{Tuple{TS1,TS2}}}
    dest = set(destspan)
    ELT = eltype(op)
    for i in eachindex(coef_srcdest)
        ni = native_index(dest,i)
        coef_srcdest[i]*=op.scale(ELT(BasisFunctions.index(ni[1])),ELT(BasisFunctions.index(ni[2])))^op.order
    end
    coef_srcdest
end
inv(op::IdxnScalingOperator) = IdxnScalingOperator(op.src, order=op.order*-1, scale=op.scale)
