# fe_solvers.jl

# TODO: the FE_Solver types should become operator types, so that we can reuse TensorProductOperator

abstract FE_Solver{SRC,DEST} <: AbstractOperator{SRC,DEST}

problem(s::FE_Solver) = s.problem

domain(s::FE_Solver) = domain(problem(s))

# Delegation methods
for op in (:numtype,:eltype,:frequency_basis,:frequency_basis_ext, :time_basis, :time_basis_ext,
    :time_basis_restricted, :operator, :operator_transpose)
    @eval $op(s::FE_Solver) = $op(problem(s))
end

size(s::FE_Solver, j::Int) = size(problem(s), j)

src(s::FE_Solver) = time_basis_restricted(s)

dest(s::FE_Solver) = frequency_basis(s)

function solve(s::AbstractOperator, f::Function, p::FE_Problem, elt = eltype(s))
    coef = Array(elt, size(frequency_basis(p)))
    println("frequency basis ",frequency_basis(p))
    println("frequency basis size",size(frequency_basis(p)))
    rhs = Array(elt, size(time_basis_restricted(p)))
        
    rhs!(p, rhs, f)
    apply!(s, coef, rhs)
    coef = reshape(coef, size(frequency_basis(p)))
    norm = normalization_operator(frequency_basis_ext(p), frequency_basis(p))
    println("sizes : ", size(s),size(coef),size(rhs))
    coef = norm * coef
    frame = DomainFrame(domain(p), frequency_basis(p))
    SetExpansion(frame, coef)
end


function solve!{T}(s::FE_Solver, coef::AbstractArray{T}, rhs::AbstractArray{T}, f::Function)
    rhs!(problem(s), rhs, f)
    apply!(s, coef, rhs)
end

## immutable FE_TensorProductSolver{ELT,N,ID,ND} <: FE_Solver
##     solvers1d :: Array{FE_Solver,1}
##     ns :: Array{Int,1}
##     ms :: Array{Int,1}
##     grid :: TensorProductGrid
##     basis :: TensorProductSet
##     basis_ext :: TensorProductSet
##     function FE_TensorProductSolver(problems::Array{FE_DiscreteProblem,1},solver_type::Tuple)
##         solvers = Array(FE_Solver,length(problems))
##         ns=Int[]
##         ms=Int[]
##         grids=AbstractGrid[]
##         bases=FunctionSet[]
##         bases_ext=FunctionSet[]
##         for i=1:length(problems)
##             solvers[i]=solver_type[i](problems[i])
##             push!(ns,size(frequency_basis(problems[i]))...)
##             push!(ms,size(time_basis_restricted(problems[i]))...)
##             push!(grids,grid(time_basis_restricted(problems[i])))
##             if dim(problems[i])==1
##                 push!(bases,frequency_basis(problems[i]))
##                 push!(bases_ext,frequency_basis_ext(problems[i]))
##             else
##                 for j=1:dim(problems[i])
##                     push!(bases,set(frequency_basis(problems[i]),j))
##                     push!(bases_ext,set(frequency_basis_ext(problems[i]),j))
##                 end
##             end
##         end
##         new(solvers,ns,ms,TensorProductGrid(grids...),TensorProductSet(bases...),TensorProductSet(bases_ext...))
##     end
## end

## function solve(s::FE_TensorProductSolver, f::Function, elt = eltype(s))
##     # Find out why eltype(s) is Type::Any
##     coef = Array(elt, size(s, 2)...)
##     rhs = Array(elt, size(s,1)...)
##     solve!(s, coef, rhs,f)
##     coef = reshape(coef, size(s.basis))
##     norm = normalization_operator(frequency_basis_ext(s), frequency_basis(s))
##     coef = norm * coef
##     frame = DomainFrame(domain(s), frequency_basis(s))
##     SetExpansion(frame, coef)
## end

## size(s::FE_TensorProductSolver,j)= j==1 ? s.ms : s.ns
## grid(s::FE_TensorProductSolver) = s.grid
## frequency_basis(s::FE_TensorProductSolver) = s.basis
## frequency_basis_ext(s::FE_TensorProductSolver) = s.basis_ext
## FE_TensorProductSolver(problems::Array{FE_DiscreteProblem,1},solver_type) = FE_TensorProductSolver{eltype(problem),sum(map(dim,problems)),length(problems),tuple(map(dim,problems)...)}(problems,solver_type)
## # HackyWacky
## problem(s::FE_TensorProductSolver) = problem(s.solvers1d[1])

## # More HackyWacky
## domain(s::FE_TensorProductSolver) = TensorProductDomain([domain(problem(i)) for i in s.solvers1d]...)

## # Ugly starting function, eventually replace by NxNx... coefficients
## function solve!{ELT,N,ID,ND,T}(s::FE_TensorProductSolver{ELT,N,ID,ND}, coef::AbstractArray{T,N}, rhs::AbstractArray{T,ID}, f::Function)
##     solvers = s.solvers1d
##     sr=size(rhs)
##     sc=size(coef)
##     rhs!(grid(s), rhs, f)
##      # The last N-1 dimensions of rhs are solved (recursively)
##     temp=zeros(T,sr[1],sc[1+ND[1]:end]...)
##     _solve!(s,coef,rhs,solvers,temp)
## end

## function _solve!{ELT,T,N,ID,ND}(s::FE_TensorProductSolver{ELT,N,ID,ND}, coef::AbstractArray{T}, rhs::AbstractArray{T}, solvers::Array{FE_Solver},temp::Array{T})
##     # The current recursion level
##     reclevel=ID-ndims(rhs)
##     # The Cumulative Dimension distribution
##     NDC=cumsum([ND...])
##     sr=size(rhs)
##     sc=size(coef)
##     # Indexing
##     indices = fill(Colon(),N)
##     # Allocate space for lower-dimensional rhs
##     tmpa=zeros(T,sr[2:end]...)
##      # The last N-1 dimensions of rhs are solved (recursively)
##     newtemp=zeros(T,sr[2],sc[NDC[reclevel+1]+2:end]...)
##     # Do a recursive solve along the first dimension.
##     for j=1:sr[1]
##         tmp=slice(rhs,j,indices[reclevel+2:ID]...)
##         _solve!(s,slice(temp,j,indices[1:N-ND[reclevel+1]]...) ,tmp+tmpa, solvers[2:end],newtemp)
##     end
##     # do the remaining 1d solves
##     # Allocate space for lower-dimensional rhs
##     tmpa=zeros(T,sr[1])
##     for i=1:prod(sc[ND[reclevel+1]+1:end])
##         tmp=slice(temp,:,i)
##         _solve!(s,slice(coef,indices[1:ND[1+reclevel]]...,i) ,tmpa+tmp, solvers[1:end],tmpa)
##     end    
##  end
## #End of the recursion
## function _solve!{T}(s::FE_TensorProductSolver, coef::AbstractArray{T}, rhs::AbstractArray{T,1}, solver1d::Array{FE_Solver}, temp::Array{T})
##     solve!(solver1d[1], coef, rhs)
## end


immutable FE_DirectSolver{ELT} <: FE_Solver
    problem ::  FE_Problem
    QR :: Base.LinAlg.QRCompactWY{ELT,Array{ELT,2}}

    function FE_DirectSolver(problem::FE_Problem)
        new(problem, qrfact(matrix(operator(problem))))
    end
end

FE_DirectSolver(problem::FE_Problem) = FE_DirectSolver{eltype(problem)}(problem)

FE_DirectSolver(p::FE_TensorProductProblem) = TensorProductOperator(map(FE_DirectSolver,p.problems)...)


function solve!{T}(s::FE_DirectSolver, coef::AbstractArray{T}, rhs::AbstractArray{T})
    coef[:] = s.QR \ rhs
end


function apply!(s::FE_DirectSolver, dest, src, coef_dest, coef_src)
    coef_dest[:] = s.QR \ coef_src
end


abstract FE_IterativeSolver <: FE_Solver


immutable FE_IterativeSolverLSQR <: FE_IterativeSolver
    problem ::  FE_Problem
end


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



