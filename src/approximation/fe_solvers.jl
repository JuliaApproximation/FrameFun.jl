
abstract type FE_Solver{T} <: BasisFunctions.AbstractSolverOperator{T} end

operator(s::FE_Solver) = s.A


struct DirectSolver{ELT} <: FE_Solver{ELT}
    A       ::  DictionaryOperator
    QR      ::  Factorization

    src_scratch :: Vector{ELT}

    function DirectSolver{ELT}(A::DictionaryOperator) where ELT
        new(A, qr(matrix(A),Val(true)), zeros(ELT, length(dest(A))))
    end
end

DirectSolver(A::DictionaryOperator{ELT}; options...) where {ELT} =
    DirectSolver{eltype(A)}(A)

function apply!(s::DirectSolver, coef_dest, coef_src::MultiArray)
    linearize_coefficients!(s.src_scratch, coef_src)
    apply!(s, coef_dest, s.src_scratch)
end

function apply!(s::DirectSolver, coef_dest, coef_src)
    coef_dest[:] = s.QR \ coef_src
end

## abstract FE_IterativeSolver <: FE_Solver


## struct FE_IterativeSolverLSQR <: FE_IterativeSolver
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
