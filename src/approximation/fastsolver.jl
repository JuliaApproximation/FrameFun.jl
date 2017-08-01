
# fastsolver.jl

"""
A fast FE solver based on a low-rank approximation of the plunge region. The plunge region
is isolated using a projection operator.
For more details, see the paper 'Fast algorithms for the computation of Fourier extensions of arbitrary length'
http://arxiv.org/abs/1509.00206
"""
struct FE_ProjectionSolver{ELT} <: FE_Solver{ELT}
    TS :: AbstractOperator
    op     ::  AbstractOperator
    plunge_op   ::  AbstractOperator    # store the operator because it allocates memory
    b
    blinear     ::  Array{ELT,1}
    x2
    x1
    scaling
    function FE_ProjectionSolver{ELT}(op::AbstractOperator, scaling; cutoff = default_cutoff(op), trunc = TruncatedSvdSolver, R = estimate_plunge_rank(op), verbose=false,options...) where ELT
        plunge_op = plunge_operator(op, scaling)
        TS = trunc(plunge_op*op; cutoff=cutoff, R=R, verbose=verbose, options...)
        b = zeros(ELT, dest(plunge_op))
        blinear = zeros(ELT, length(dest(plunge_op)))
        x1 = zeros(ELT, src(op))
        x2 = zeros(ELT, src(op))
        new(TS, op, plunge_op, b,blinear,x1,x2,scaling)
    end
end

FE_ProjectionSolver(op::AbstractOperator, scaling; options...) =
    FE_ProjectionSolver{eltype(op)}(op, scaling; options...)

function plunge_operator(op, scaling)
    A = op
    Ap = op'
    I = ScalingOperator(dest(A),scaling)

    A*Ap - I
end

default_cutoff(op::AbstractOperator) = 10^(4/5*log10(eps(real(eltype(op)))))
function estimate_plunge_rank(op::AbstractOperator)
    nml=length(src(op))^2/length(dest(op))
    N = dimension(src(op))
    if N==1
        return min(round(Int, 9*log(nml)),length(src(op)))
    else
        return min(round(Int, 9*log(nml)*nml^((N-1)/N)),length(src(op)))
    end
end

apply!(s::FE_ProjectionSolver, dest, src, coef_dest, coef_src) =
    apply!(s, dest, src, coef_dest, coef_src, s.op, s.op', s.plunge_op, s.x1, s.x2)

function apply!(s::FE_ProjectionSolver, destset, srcset, coef_dest, coef_src, A, At, P, x1, x2)
    # Applying plunge to the right hand side
    apply!(P, s.b, coef_src)
    BasisFunctions.linearize_coefficients!(dest(A), s.blinear, s.b)
    apply!(s.TS,x2,s.blinear)
    # x2 solves the middle guys
    apply!(A, s.b, x2)
    apply!(At, x1, coef_src-s.b)
    # x1 solves the remainder through one application of the operator
    for i in eachindex(x1)
        coef_dest[i] = x1[i]/s.scaling+x2[i]
    end
end
