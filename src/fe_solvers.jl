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
end

function FE_DirectSolver(problem::FE_Problem)
    op = FE_DiscreteOperator(problem)
    mat = matrix(op)
    FE_DirectSolver{eltype(problem)}(op, mat)
end

eltype{ELT}(s::FE_DirectSolver{ELT}) = ELT

operator(s::FE_DirectSolver) = s.op


function solve!{T}(s::FE_DirectSolver, coef::Array{T}, rhs::Array{T})
    coef[:] = s.matrix \ rhs
end


immutable FE_IterativeSolver{ELT} <: FE_Solver
    op      ::  FE_DiscreteOperator
end

operator(s::FE_IterativeSolver) = s.op






