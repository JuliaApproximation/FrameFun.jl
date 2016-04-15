# fe_solvers.jl


abstract FE_Solver{SRC,DEST} <: AbstractOperator{SRC,DEST}

problem(s::FE_Solver) = s.problem

# Delegation methods
for op in (:frequency_basis, :frequency_basis_ext, :time_basis, :time_basis_ext,
    :time_basis_restricted, :operator, :operator_transpose, :domain)
    @eval $op(s::FE_Solver) = $op(problem(s))
end

size(s::FE_Solver, j::Int) = size(problem(s), j)

src(s::FE_Solver) = time_basis_restricted(s)

dest(s::FE_Solver) = frequency_basis(s)



immutable FE_DirectSolver{SRC,DEST} <: FE_Solver{SRC,DEST}
    problem ::  FE_Problem
    QR      ::  Factorization

    function FE_DirectSolver(problem::FE_Problem)
        new(problem, qrfact(matrix(operator(problem))))
    end
end

function FE_DirectSolver(problem::FE_Problem; options...)
    SRC = typeof(time_basis_restricted(problem))
    DEST = typeof(frequency_basis(problem))
    FE_DirectSolver{SRC,DEST}(problem)
end

function apply!(s::FE_DirectSolver, dest, src, coef_dest, coef_src)
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
