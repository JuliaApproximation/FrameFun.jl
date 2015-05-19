# fe_solvers.jl


abstract FE_Solver

frequency_basis(s::FE_Solver) = frequency_basis(operator(s))

time_basis(s::FE_Solver) = time_basis(operator(s))


function solve(s::FE_Solver, f::Function, elt = eltype(operator(s)))
    coef = Array(elt, size(s, 2))
    rhs = Array(elt, size(s, 1))
    solve!(s, coef, rhs, f)
    SetExpansion(frequency_basis(s), coef)
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

numtype(s::FE_DirectSolver) = numtype(s.op)

operator(s::FE_DirectSolver) = s.op

size(s::FE_DirectSolver) = size(s.matrix)

size(s::FE_DirectSolver, j::Int) = size(s.matrix, j)


function solve!{T}(s::FE_DirectSolver, coef::Array{T}, rhs::Array{T})
    @assert length(coef) == size(s, 2)
    @assert length(rhs) == size(s, 1)

    coef[:] = s.matrix \ rhs
end

function solve!{T}(s::FE_DirectSolver, coef::Array{T}, rhs::Array{T}, f::Function)
    rhs!(operator(s), rhs, f)
    solve!(s, coef, rhs)
end


