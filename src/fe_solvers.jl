# fe_solvers.jl


abstract FE_Solver

problem(s::FE_Solver) = s.problem

# Delegation methods
for op in (:numtype,:eltype,:frequency_basis,:frequency_basis_ext, :time_basis, :time_basis_ext,
    :restricted_time_basis, :size, :operator, :operator_transpose)
    @eval $op(s::FE_Solver) = $op(problem(s))
end

size(s::FE_Solver, j::Int) = size(problem(s), j)

function solve(s::FE_Solver, f::Function, elt = eltype(s))
    coef = Array(elt, size(s, 2))
    rhs = Array(elt, size(s, 1))
    solve!(s, coef, rhs, f)
    SetExpansion(frequency_basis(s), coef)
end

function solve!{T}(s::FE_Solver, coef::Array{T}, rhs::Array{T}, f::Function)
    rhs!(problem(s), rhs, f)
    solve!(s, coef, rhs)
end





immutable FE_DirectSolver{ELT} <: FE_Solver
    problem ::  FE_Problem
    matrix  ::  Array{ELT,2}

    FE_DirectSolver(problem::FE_Problem) = new(problem, matrix(operator(problem)))
end

FE_DirectSolver(problem::FE_Problem) = FE_DirectSolver{eltype(problem)}(problem)

eltype{ELT}(s::FE_DirectSolver{ELT}) = ELT


function solve!{T}(s::FE_DirectSolver, coef::Array{T}, rhs::Array{T})
    coef[:] = s.matrix \ rhs
end



abstract FE_IterativeSolver <: FE_Solver


immutable FE_IterativeSolverLSQR <: FE_IterativeSolver
    problem ::  FE_Problem
end


function solve!{T}(s::FE_IterativeSolverLSQR, coef::Array{T}, rhs::Array{T})
    op = operator(s)
    opt = operator_transpose(s)

    my_A_mul_B!(output, x) =  ( apply!(op,  reshape(output, size(dest(op ))), reshape(x, size(src(op )))); output )
    my_Ac_mul_B!(output, y) = ( apply!(opt, reshape(output, size(dest(opt))), reshape(y, size(src(opt)))); output )

    matcfcn = MatrixCFcn{T}(size(op, 1), size(op, 2), my_A_mul_B!, my_Ac_mul_B!)

    coef[:] = 0
    y,ch = lsqr!(coef, matcfcn, rhs, maxiter = 100)

    println("Stopped after ", ch.mvps, " iterations with residual ", abs(ch.residuals[end]), ".")
end



