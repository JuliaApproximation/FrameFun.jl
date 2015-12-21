# fe_solvers.jl


abstract FE_Solver{ELT,SRC,DEST} <: AbstractOperator{SRC,DEST}

eltype{ELT,SRC,DEST}(::Type{FE_Solver{ELT,SRC,DEST}}) = ELT

problem(s::FE_Solver) = s.problem

domain(s::FE_Solver) = domain(problem(s))

# Delegation methods
for op in (:frequency_basis, :frequency_basis_ext, :time_basis, :time_basis_ext,
    :time_basis_restricted, :operator, :operator_transpose)
    @eval $op(s::FE_Solver) = $op(problem(s))
end

size(s::FE_Solver, j::Int) = size(problem(s), j)

src(s::FE_Solver) = time_basis_restricted(s)

dest(s::FE_Solver) = frequency_basis(s)

function solve(s::AbstractOperator, f::Function, p::FE_Problem)
    ELT = eltype(s)
    coef = Array(ELT, size(frequency_basis(p)))
    rhs = Array(ELT, size(time_basis_restricted(p)))
        
    rhs!(p, rhs, f)
    apply!(s, coef, rhs)
    coef = reshape(coef, size(frequency_basis(p)))
    norm = normalization_operator(frequency_basis_ext(p), frequency_basis(p), ELT)
    coef = norm * coef
    frame = DomainFrame(domain(p), frequency_basis(p))
    SetExpansion(frame, coef)
end


function solve!{T}(s::FE_Solver, coef::AbstractArray{T}, rhs::AbstractArray{T}, f::Function)
    rhs!(problem(s), rhs, f)
    apply!(s, coef, rhs)
end




immutable FE_DirectSolver{ELT,SRC,DEST} <: FE_Solver{ELT,SRC,DEST}
    problem ::  FE_Problem
    QR      ::  Factorization{ELT}

    function FE_DirectSolver(problem::FE_Problem)
        new(problem, qrfact(matrix(operator(problem))))
    end
end

function FE_DirectSolver(problem::FE_Problem)
    ELT = eltype(problem)
    SRC = typeof(time_basis_restricted(problem))
    DEST = typeof(frequency_basis(problem))
    FE_DirectSolver{ELT,SRC,DEST}(problem)
end

FE_DirectSolver(p::FE_TensorProductProblem) = TensorProductOperator(map(FE_DirectSolver,p.problems)...)


function solve!{T}(s::FE_DirectSolver, coef::AbstractArray{T}, rhs::AbstractArray{T})
    coef[:] = s.QR \ rhs
end


function apply!(s::FE_DirectSolver, dest, src, coef_dest, coef_src)
    coef_dest[:] = matrix(operator(problem(s))) \ coef_src
end

## immutable FE_DirectSolver{ELT} <: FE_Solver
##     problem ::  FE_Problem

##     function FE_DirectSolver(problem::FE_Problem)
##         new(problem)
##     end
## end

## FE_DirectSolver(problem::FE_Problem) = FE_DirectSolver{eltype(problem)}(problem)

## FE_DirectSolver(p::FE_TensorProductProblem) = TensorProductOperator(map(FE_DirectSolver,p.problems)...)


## @debug function apply!(s::FE_DirectSolver, dest, src, coef_dest, coef_src)
##     @bp
##     coef_dest[:] = matrix(operator(problem(s))) \ coef_src
## end

## abstract FE_IterativeSolver <: FE_Solver


## immutable FE_IterativeSolverLSQR <: FE_IterativeSolver
##     problem ::  FE_Problem
## end


## function solve!{T}(s::FE_IterativeSolverLSQR, coef::AbstractArray{T}, rhs::AbstractArray{T})
##     op = operator(s)
##     opt = operator_transpose(s)

##     my_A_mul_B!(output, x) =  ( apply!(op,  reshape(output, size(dest(op ))), reshape(x, size(src(op )))); output )
##     my_Ac_mul_B!(output, y) = ( apply!(opt, reshape(output, size(dest(opt))), reshape(y, size(src(opt)))); output )

##     matcfcn = MatrixCFcn{T}(size(op, 1), size(op, 2), my_A_mul_B!, my_Ac_mul_B!)

##     coef[:] = 0
##     y,ch = lsqr!(coef, matcfcn, rhs, maxiter = 100)

##     println("Stopped after ", ch.mvps, " iterations with residual ", abs(ch.residuals[end]), ".")
## end



