# fe_solvers.jl


abstract type FE_Solver{ELT} <: AbstractOperator{ELT} end

op(s::FE_Solver) = s.op

## # Delegation methods
## for op in (:frequency_basis, :frequency_basis_ext, :time_basis, :time_basis_ext,
##            :time_basis_restricted, :operator, :operator_transpose, :domain)
##     @eval $op(s::FE_Solver) = $op(op(s))
## end

size(s::FE_Solver, j::Int) = size(transpose(op(s)), j)

src(s::FE_Solver) = dest(op(s))

dest(s::FE_Solver) = src(op(s))



struct FE_DirectSolver{ELT} <: FE_Solver{ELT}
    op      ::  AbstractOperator
    QR      ::  Factorization

    function FE_DirectSolver{ELT}(op::AbstractOperator,scaling) where ELT
        new(op, qrfact(matrix(op),Val{true}))
    end
end

FE_DirectSolver{ELT}(op::AbstractOperator{ELT}, scaling; options...) =
    FE_DirectSolver{eltype(op)}(op,scaling)

function apply!(s::FE_DirectSolver, coef_dest, coef_src)
    coef_dest[:] = s.QR \ coef_src
end

immutable ContinuousDirectSolver{T} <: AbstractOperator{T}
  src                     :: FunctionSet
  mixedgramfactorization  :: Factorization
  normalizationofb        :: AbstractOperator
  scratch                 :: Array{T,1}
end

dest(s::ContinuousDirectSolver) = s.src

ContinuousDirectSolver(frame::ExtensionFrame; options...) =
    ContinuousDirectSolver{eltype(frame)}(frame, qrfact(matrix(MixedGram(frame; options...)),Val{true}), DualGram(basis(frame); options...), zeros(eltype(frame),length(frame)))

function apply!(s::ContinuousDirectSolver, coef_dest, coef_src)
  apply!(s.normalizationofb, s.scratch, coef_src)
  coef_dest[:] = s.mixedgramfactorization \ s.scratch
end

immutable ContinuousTruncatedSolver{T} <: AbstractOperator{T}
  src                     :: FunctionSet
  mixedgramsvd            :: AbstractOperator
  normalizationofb        :: AbstractOperator
  scratch                 :: Array{T,1}
end

dest(s::ContinuousTruncatedSolver) = s.src

function ContinuousTruncatedSolver(frame::ExtensionFrame; cutoff=1e-5, fullsvd=false, options...)
    fullsvd? solver = FrameFun.ExactTruncatedSvdSolver: solver = FrameFun.TruncatedSvdSolver
    ContinuousTruncatedSolver{eltype(frame)}(frame, solver(MixedGram(frame; options...); cutoff=cutoff, options...), DualGram(basis(frame); options...), zeros(eltype(frame),length(frame)))
end

function ContinuousTruncatedSvdSolver(frame::ExtensionFrame; cutoff=1e-5, fullsvd=false, options...)
    fullsvd? solver = FrameFun.ExactTruncatedSvdSolver: solver = FrameFun.TruncatedSvdSolver
    ContinuousTruncatedSolver{eltype(frame)}(frame, solver(MixedGram(frame; options...); cutoff=cutoff, options...), DualGram(basis(frame); options...), zeros(eltype(frame),length(frame)))
end

function apply!(s::ContinuousTruncatedSolver, coef_dest, coef_src)
  apply!(s.normalizationofb, s.scratch, coef_src)
  coef_dest[:] = s.mixedgramsvd*s.scratch
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
