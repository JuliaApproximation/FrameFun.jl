# fe_solvers.jl


abstract FE_Solver

frequency_basis(s::FE_Solver) = frequency_basis(operator(s))

time_basis(s::FE_Solver) = time_basis(operator(s))

eltype(s::FE_Solver) = eltype(operator(s))

numtype(s::FE_Solver) = numtype(operator(s))

size(s::FE_Solver) = size(operator(s))

size(s::FE_Solver, j::Int) = size(operator(s), j)


function solve(s::FE_Solver, f::Function, elt = eltype(operator(s)))
    coef = Array(elt, size(s, 2))
    rhs = Array(elt, size(s, 1))
    solve!(s, coef, rhs, f)
    SetExpansion(frequency_basis(s), coef)
end

function solve!{T}(s::FE_Solver, coef::Array{T}, rhs::Array{T}, f::Function)
    rhs!(operator(s), rhs, f)
    solve!(s, coef, rhs)
end





immutable FE_DirectSolver{ELT} <: FE_Solver
    op      ::  FE_DiscreteOperator
    matrix  ::  Array{ELT,2}

    FE_DirectSolver(op::FE_DiscreteOperator) = new(op, matrix(op))
end

FE_DirectSolver(problem::FE_Problem) = FE_DirectSolver{eltype(problem)}(FE_DiscreteOperator(problem))

eltype{ELT}(s::FE_DirectSolver{ELT}) = ELT

operator(s::FE_DirectSolver) = s.op


function solve!{T}(s::FE_DirectSolver, coef::Array{T}, rhs::Array{T})
    coef[:] = s.matrix \ rhs
end


abstract FE_IterativeSolver <: FE_Solver


immutable FE_IterativeSolverLSQR <: FE_IterativeSolver
    op      ::  FE_DiscreteOperator
    opt     ::  OperatorTranspose

    FE_IterativeSolverLSQR(op::FE_DiscreteOperator) = new(op, transpose(op))
end

FE_IterativeSolverLSQR(problem::FE_Problem) = FE_IterativeSolverLSQR(FE_DiscreteOperator(problem))


operator(s::FE_IterativeSolver) = s.op

function solve!{T}(s::FE_IterativeSolverLSQR, coef::Array{T}, rhs::Array{T})
    op = s.op
    opt = s.opt

    my_A_mul_B!(output, x) =  ( apply!(op,  reshape(output, size(dest(op ))), reshape(x, size(src(op )))); output )
    my_Ac_mul_B!(output, y) = ( apply!(opt, reshape(output, size(dest(opt))), reshape(y, size(src(opt)))); output )

    matcfcn = MatrixCFcn{T}(size(op, 1), size(op, 2), my_A_mul_B!, my_Ac_mul_B!)

    coef[:] = 0
    y,ch = my_lsqr!(coef, matcfcn, rhs, maxiter = 100)

    println("Stopped after ", ch.mvps, " iterations with residual ", abs(ch.residuals[end]), ".")
end



