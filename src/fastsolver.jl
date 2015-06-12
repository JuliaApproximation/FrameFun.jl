

immutable FE_ProjectionSolver{ELT} <: FE_Solver
    problem     ::  FE_DiscreteProblem
    plunge_op   ::  AbstractOperator    # store the operator because it allocates memory
    W           ::  AffineMap{ELT}
    m           ::  Array{ELT,2}


    function FE_ProjectionSolver(problem::FE_DiscreteProblem)
        plunge_op = plunge_operator(problem)
        R = estimate_plunge_rank(problem)
        W = AffineMap( map(ELT, rand(param_N(problem), R)) )
        m = matrix(plunge_op * operator(problem) * W)

        new(problem, plunge_op, W, m)
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
    F=svd(s.m)
    smax=maximum(find(F[2].>1e-13))
    iS=diagm(1./F[2][1:smax])
    V=F[:3]
    U=F[:1]
    y=V[:,1:smax]*(iS*(U[:,1:smax]'*b))
    x2 = reshape(s.W * y,size(src(A)))
    x1 = At * (rhs - A*x2)/L
    coef[:] = x1+x2
end


