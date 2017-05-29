
# fastsolver.jl

"""
A fast FE solver based on a low-rank approximation of the plunge region. The plunge region
is isolated using a projection operator.
For more details, see the paper 'Fast algorithms for the computation of Fourier extensions of arbitrary length'
http://arxiv.org/abs/1509.00206
"""
struct FE_ProjectionSolver{ELT} <: FE_Solver{ELT}
    TS :: AbstractOperator
    problem     ::  FE_DiscreteProblem
    plunge_op   ::  AbstractOperator    # store the operator because it allocates memory
    b
    blinear     ::  Array{ELT,1}
    x2
    x1

    function FE_ProjectionSolver{ELT}(problem::FE_DiscreteProblem; cutoff = default_cutoff(problem), trunc = TruncatedSvdSolver, R = estimate_plunge_rank(problem), options...) where ELT
        TS = trunc(plunge_operator(problem)*operator(problem); cutoff=cutoff, R=R, verbose=true, options...)
        plunge_op = plunge_operator(problem)
        b = zeros(ELT, dest(plunge_op))
        blinear = zeros(ELT, length(dest(plunge_op)))
        x1 = zeros(ELT, src(operator(problem)))
        x2 = zeros(ELT, src(operator(problem)))
        new(TS, problem, plunge_op, b,blinear,x1,x2)
    end
end

FE_ProjectionSolver(problem::FE_DiscreteProblem; options...) =
    FE_ProjectionSolver{eltype(problem)}(problem; options...)


function plunge_operator(problem::FE_DiscreteProblem)
    A = operator(problem)
    Ap = operator_transpose(problem)
    I = IdentityOperator(time_basis_restricted(problem))

    A*Ap - I
end

default_cutoff(problem::FE_DiscreteProblem) = 10^(4/5*log10(eps(numtype(frequency_basis(problem)))))
estimate_plunge_rank{N}(problem::FE_DiscreteProblem{N}) = min(round(Int, 9*log(param_N(problem))*(param_M(problem)*param_N(problem)/param_L(problem))^(1-1/N) + 2),param_N(problem))

estimate_plunge_rank(problem::FE_DiscreteProblem{1,BigFloat}) = round(Int, 28*log(param_N(problem)) + 5)

apply!(s::FE_ProjectionSolver, dest, src, coef_dest, coef_src) =
    apply!(s, dest, src, coef_dest, coef_src, operator(s), operator_transpose(s), s.plunge_op, s.x1, s.x2)

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
        x1[i] += x2[i]
    end
    # normalization
    apply!(normalization(problem(s)), coef_dest, x1)
end
