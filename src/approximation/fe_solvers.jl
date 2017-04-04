# fe_solvers.jl


abstract FE_Solver{ELT} <: AbstractOperator{ELT}

op(s::FE_Solver) = s.op

## # Delegation methods
## for op in (:frequency_basis, :frequency_basis_ext, :time_basis, :time_basis_ext,
##            :time_basis_restricted, :operator, :operator_transpose, :domain)
##     @eval $op(s::FE_Solver) = $op(op(s))
## end

size(s::FE_Solver, j::Int) = size(transpose(op(s)), j)

src(s::FE_Solver) = dest(op(s))

dest(s::FE_Solver) = src(op(s))



immutable FE_DirectSolver{ELT} <: FE_Solver{ELT}
    problem ::  FE_Problem
    QR      ::  Factorization

    function FE_DirectSolver(problem::FE_Problem)
        new(problem, qrfact(matrix(operator(problem)),Val{true}))
    end
end

FE_DirectSolver(problem::FE_Problem; options...) =
    FE_DirectSolver{eltype(problem)}(problem)

function apply!(s::FE_DirectSolver, coef_dest, coef_src)
    coef_dest[:] = s.QR \ coef_src
    apply!(normalization(problem(s)), coef_dest, coef_dest)
end


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
