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
    rhs = Array(elt, size(time_basis_restricted(problem(s))))
    solve!(s, coef, rhs,f)
    SetExpansion(frequency_basis(s), coef)
end


function solve!{T}(s::FE_Solver, coef::Array{T}, rhs::Array{T}, f::Function)
    rhs!(problem(s), rhs, f)
    solve!(s, coef, rhs)
end

immutable FE_TensorProductSolver{ELT} <: FE_Solver
    problem :: FE_Problem
    solvers1d :: Array{FE_Solver,1}
    function FE_TensorProductSolver(problem::FE_Problem)
        fbasis=frequency_basis(problem)
        tbasise=time_basis_ext(problem)
        tbasisr=time_basis_restricted(problem)
        solvers=FE_Solver[]
        for i=1:dim(problem)
            # subproblem parameters
            a=left(tbasisr,i)
            b=right(tbasisr,i)
            n=round(Int,(size(fbasis,i)-1)/2)
            T=(size(tbasise,i))/(size(tbasisr,i)-1)
            gamma=(size(tbasisr,i)-1)/(size(fbasis,i)-1)
            subproblem=default_fourier_problem(Interval(a,b),n,T,gamma)
            # Collect all solvers
            push!(solvers,FE_ProjectionSolver(subproblem))
        end
        new(problem,solvers)
    end
    
end

FE_TensorProductSolver(problem::FE_Problem) = FE_TensorProductSolver{eltype(problem)}(problem)

function solve!{T,N}(s::FE_TensorProductSolver, coef::Array{T}, rhs::Array{T,N}, solvers::Array{FE_Solver}=s.solvers1d)
    println("entering recursion, coef = $(size(coef)), rhs = $(size(rhs))")
    sr=size(rhs)
    sc=size(frequency_basis(s.problem))
    # Don't know is rhs is supposed to be preserved
    interold = deepcopy(rhs)
    # Indexing
    indices = fill(Colon(),N)
    for i = 1:2
        println(" i = $i/$N")
        internew = zeros(T,sc[1:i]...,sr[i+1:end]...)
        # Allocate space for lower-dimensional rhs
        tmpa=zeros(T,sr[1:i-1]...,sr[i+1:N]...)
        # Do a recursive solve for each dimensional row (this won't work for three dimensions)
        intertemp=zeros(T,sc[[1:i-1;i+1:N]]...)

        for j=1:size(interold,N-i+1)
            println(" j = $j/$(size(interold,N-i+1))")
            tmp=slice(interold,indices[i+1:N]...,j,indices[1:i-1]...)
            solve!(s,intertemp ,tmpa+tmp, solvers[[1:i-1;i+1:end]])
            slice(internew,indices[i+1:N]...,j,indices[1:i-1]...)[:]=intertemp[:]
        end
        interold=internew
    end
    # remove this in favor of N-dimensional coefficients?
    coef[:] = interold[:]
 end
#End of the recursion
function solve!{T}(s::FE_TensorProductSolver, coef::Array{T}, rhs::Array{T,1}, solver1d::Array{FE_Solver})
    solve!(solver1d[1], coef, rhs)
end


immutable FE_DirectSolver{ELT} <: FE_Solver
    problem ::  FE_Problem
    QR :: Base.LinAlg.QRCompactWY{ELT,Array{ELT,2}}

    function FE_DirectSolver(problem::FE_Problem)
        new(problem, qrfact(matrix(operator(problem))))
    end
end

FE_DirectSolver(problem::FE_Problem) = FE_DirectSolver{eltype(problem)}(problem)

eltype{ELT}(s::FE_DirectSolver{ELT}) = ELT


function solve!{T}(s::FE_DirectSolver, coef::Array{T}, rhs::Array{T})
    coef[:] = s.QR \ rhs
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



