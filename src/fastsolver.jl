

immutable FE_ProjectionSolver{ELT} <: FE_Solver
    problem     ::  FE_DiscreteProblem
    plunge_op   ::  AbstractOperator    # store the operator because it allocates memory
    W           ::  MatrixOperator
    Ut          ::  Array{ELT,2}
    VS           ::  Array{ELT,2}
    b           :: Array{ELT,1}
    y           :: Array{ELT,1}
    sy          :: Array{ELT,1}
    x2          :: Array{ELT}
    x1          :: Array{ELT}

    function FE_ProjectionSolver(problem::FE_DiscreteProblem)
        plunge_op = plunge_operator(problem)
        R = estimate_plunge_rank(problem)
        W = MatrixOperator( map(ELT, rand(param_N(problem), R)) )
        # Added
        USV= svd(matrix(plunge_op * operator(problem) * W))
        maxind=maximum(find(USV[2].>1e-6))
        S=USV[2]
        Sinv=1./S[1:maxind]
        b=zeros(size(dest(plunge_op)))
        y=zeros(size(USV[3],1))
        x1=zeros(size(src(operator(problem))))
        x2=zeros(size(src(operator(problem))))
        sy=zeros(maxind,)
        new(problem, plunge_op, W, USV[1][:,1:maxind]',USV[3][:,1:maxind]*diagm(Sinv[:]),b,y,sy,x1,x2)
    end
end

FE_ProjectionSolver(problem::FE_DiscreteProblem) = FE_ProjectionSolver{eltype(problem)}(problem)

function plunge_operator(problem::FE_DiscreteProblem)
    A = operator(problem)
    Ap = operator_transpose(problem)
    I = IdentityOperator(time_basis_restricted(problem))

    A*Ap - I
end


estimate_plunge_rank{N}(problem::FE_DiscreteProblem{N}) = min(round(Int, 9*log(param_N(problem))*(param_M(problem)*param_N(problem)/param_L(problem))^(1-1/N) + 2),param_N(problem))

estimate_plunge_rank(problem::FE_DiscreteProblem{1,BigFloat}) = round(Int, 28*log(param_N(problem)) + 5)

function solve!{T}(s::FE_ProjectionSolver, coef::AbstractArray{T}, rhs::AbstractArray{T})
    A = operator(s)
    At = operator_transpose(s)
    
    P = s.plunge_op
    apply!(P,s.b,rhs)
    A_mul_B!(s.sy,s.Ut,s.b)
    A_mul_B!(s.y,s.VS,s.sy)
    apply!(s.W,s.x2,s.y)
    #x2 = reshape(s.W * y,size(src(A)))
    apply!(A,s.b,s.x2)
    apply!(At,s.x1,rhs-s.b)
    for i=1:length(coef)
        coef[i]=s.x1[i]+s.x2[i]
    end
    println("solution norm", norm(A*coef-rhs))
    println("original",rhs)
    println("coef",coef)
    
end


