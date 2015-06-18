

immutable FE_ProjectionSolver{ELT} <: FE_Solver
    problem     ::  FE_DiscreteProblem
    plunge_op   ::  AbstractOperator    # store the operator because it allocates memory
    W           ::  AffineMap{ELT}
    Ut          ::  Array{ELT,2}
    V           ::  Array{ELT,2}
    Sinv        ::  Array{ELT,2}


    function FE_ProjectionSolver(problem::FE_DiscreteProblem)
        plunge_op = plunge_operator(problem)
        R = estimate_plunge_rank(problem)
        W = AffineMap( map(ELT, rand(param_N(problem), R)) )
        USV= svd(matrix(plunge_op * operator(problem) * W))
        maxind=maximum(find(USV[2].>1e-6))
        S=USV[2]
        Sinv=[1./S[1:maxind];zeros(length(S)-maxind,1)]
        new(problem, plunge_op, W, USV[1]',USV[3],diagm(Sinv[:]))
    end
end

FE_ProjectionSolver(problem::FE_DiscreteProblem) = FE_ProjectionSolver{eltype(problem)}(problem)

function plunge_operator(problem::FE_DiscreteProblem)
    A = operator(problem)
    Ap = operator_transpose(problem)
    I = IdentityOperator(time_basis_restricted(problem))
    lambda = param_L(problem)

    A*Ap - lambda * I
end

estimate_plunge_rank{N}(problem::FE_DiscreteProblem{N}) = min(round(Int, 9*log(param_N(problem))*(param_M(problem)*param_N(problem)/param_L(problem))^(1-1/N) + 2),param_N(problem))

estimate_plunge_rank(problem::FE_DiscreteProblem{1,BigFloat}) = round(Int, 18*log(param_N(problem)) + 5)

function solve!(s::FE_ProjectionSolver, coef::Array, rhs::Array)
    A = operator(s)
    At = operator_transpose(s)
    P = s.plunge_op
    L = param_L(problem(s))

    b = P*rhs
    y = (s.V*(s.Sinv*(s.Ut*b)))
    x2 = reshape(s.W * y,size(src(A)))
    x1 = At * (rhs - A*x2)/L
    coef[:] = x1+x2
end


