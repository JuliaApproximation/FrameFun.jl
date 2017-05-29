# smoothsolver.jl

"""
A fast FE solver based on a low-rank approximation of the plunge region. The plunge region
is isolated using a projection operator. This algorithm contains an extra smoothing step

"""
struct FE_SmoothProjectionSolver{ELT} <: FE_Solver{ELT}
    TS :: TruncatedSvdSolver
        problem     ::  FE_DiscreteProblem
        plunge_op   ::  AbstractOperator    # store the operator because it allocates memory
        b     ::  Array{ELT,1}
        blinear     ::  Array{ELT,1}
        syv         ::  Array{ELT,1}
        x2          ::  Array{ELT}
        x1          ::  Array{ELT}
        x3          ::  Array{ELT}
        Q          ::  Array{ELT,2}
        D           ::  AbstractOperator

    function FE_SmoothProjectionSolver{ELT}(problem::FE_DiscreteProblem; cutoff = default_cutoff(problem), cutoffv=sqrt(cutoff), R = estimate_plunge_rank(problem), options...) where ELT
        plunge_op = plunge_operator(problem)
        # Create Random matrices
        TS1 = TruncatedSvdSolver(plunge_operator(problem)*operator(problem); cutoff = cutoff, options...)
        TS2 = TruncatedSvdSolver(operator_transpose(problem)*plunge_operator(problem); cutoff = cutoff, options...)
        # D = Sobolev operator
        D = IdxnScalingOperator(frequency_basis(problem); options...)
        AD = inv(D)
        ADV = (TS2.Ut)'.*diagonal(AD)
        # Orthogonal basis for D^(-1)V_mid
        Q, R = qr(ADV)
        # Pre-allocation
        b = zeros(size(dest(plunge_op)))
        blinear = zeros(ELT, length(dest(plunge_op)))
        x1 = zeros(size(src(operator(problem))))
        x2 = zeros(size(src(operator(problem))))
        x3 = zeros(size(src(operator(problem))))
        syv = zeros(size(TS2.Ut,1))
        new(TS1,problem, plunge_op, b,blinear,syv,x1,x2,x3,Q, D)
    end
end


FE_SmoothProjectionSolver(problem::FE_DiscreteProblem; options...) =
        FE_SmoothProjectionSolver{eltype(problem)}(problem; options...)

function apply!(s::FE_SmoothProjectionSolver, destarg, src, coef_dest, coef_src)
    A = operator(s)
    At = operator_transpose(s)
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
    apply!(MatrixOperator(s.Q'),s.syv,s.x1)
    apply!(MatrixOperator(s.Q),s.x3,s.syv)
    for i = 1:length(s.x1)
        s.x1[i] = s.x1[i]- s.x3[i]
    end
    apply!(AD,s.x3,s.x1)
    for i = 1:length(s.x2)
        s.x2[i] = s.x2[i] - s.x3[i]
    end

    # post smoothing step
    apply!(A, s.b, s.x2)
    apply!(At, s.x1, coef_src-s.b)
    for i = 1:length(coef_dest)
        coef_dest[i] = s.x1[i] + s.x2[i]
    end

    # Normalization
    apply!(normalization(problem(s)), coef_dest, coef_dest)
end

# For FunctionSets that have a DC component
dc_index(b::ChebyshevBasis) = 1
dc_index(b::FourierBasis) = 1

# An index scaling operator, used to generate weights for the polynomial scaling algorithm.
struct IdxnScalingOperator{ELT} <: AbstractOperator{ELT}
    src     ::  FunctionSet
    order   ::  Int
    scale   ::  Function
end

IdxnScalingOperator(src::FunctionSet; order=1, scale = default_scaling_function) =
    IdxnScalingOperator{eltype(src)}(src, order, scale)

dest(op::IdxnScalingOperator) = src(op)

default_scaling_function(i) = 10.0^-4+(abs(i))+abs(i)^2+abs(i)^3
default_scaling_function(i,j) = 1+(abs(i)^2+abs(j)^2)

is_inplace(::IdxnScalingOperator) = true
is_diagonal(::IdxnScalingOperator) = true

ctranspose(op::IdxnScalingOperator) = DiagonalOperator(src(op), conj(diagonal(op)))
function apply_inplace!(op::IdxnScalingOperator, dest, src, coef_srcdest)
    ELT = eltype(op)
    for i in eachindex(dest)
        coef_srcdest[i] *= op.scale(ELT(BasisFunctions.index(native_index(dest,i))))^op.order
    end
    coef_srcdest
end

function apply_inplace!{TS1,TS2}(op::IdxnScalingOperator, dest::TensorProductSet{Tuple{TS1,TS2}}, src, coef_srcdest)
    ELT = eltype(op)
    for i in eachindex(coef_srcdest)
        ni = native_index(dest,i)
        coef_srcdest[i]*=op.scale(ELT(BasisFunctions.index(ni[1])),ELT(BasisFunctions.index(ni[2])))^op.order
    end
    coef_srcdest
end
inv(op::IdxnScalingOperator) = IdxnScalingOperator(op.src, order=op.order*-1, scale=op.scale)
