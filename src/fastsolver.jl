

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

estimate_plunge_rank(problem::FE_DiscreteProblem) = round(Int, 9*log(param_N(problem)) + 2)

estimate_plunge_rank(problem::FE_DiscreteProblem{1,BigFloat}) = round(Int, 18*log(param_N(problem)) + 5)

function solve!(s::FE_ProjectionSolver, coef::Array, rhs::Array)
    A = operator(s)
    At = operator_transpose(s)
    P = s.plunge_op
    L = param_L(problem(s))

    u = P*rhs

    y = s.m \ (P * rhs)
    x2 = s.W * y
    x1 = At * (rhs - A*x2)/L
    coef[:] = x1+x2
end


