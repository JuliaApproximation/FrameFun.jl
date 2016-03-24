# fastsolver.jl

"""
A fast FE solver based on a low-rank approximation of the plunge region. The plunge region
is isolated using a projection operator.
For more details, see the paper 'Fast algorithms for the computation of Fourier extensions of arbitrary length'
http://arxiv.org/abs/1509.00206
"""
immutable FE_ProjectionSolver{ELT,SRC,DEST} <: FE_Solver{ELT,SRC,DEST}
    problem     ::  FE_DiscreteProblem
    plunge_op   ::  AbstractOperator    # store the operator because it allocates memory
    W           ::  MatrixOperator
    Ut          ::  Array{ELT,2}
    VS          ::  Array{ELT,2}
    b           ::  Array{ELT,1}
    y           ::  Array{ELT,1}
    sy          ::  Array{ELT,1}
    x2          ::  Array{ELT}
    x1          ::  Array{ELT}

    function FE_ProjectionSolver(problem::FE_DiscreteProblem)
        plunge_op = plunge_operator(problem)
        R = estimate_plunge_rank(problem)
        W = MatrixOperator( map(ELT, rand(param_N(problem), R)) )
        ## println("max operator forward",maximum(svd(matrix(operator(problem).op2))[2]))
        ## println("min operator forward",minimum(svd(matrix(operator(problem).op2))[2]))
        ## println("max operator backward",maximum(svd(matrix(operator_transpose(problem).op2))[2]))
        ## println("min operator backward",minimum(svd(matrix(operator_transpose(problem).op2))[2]))
        USV = LAPACK.gesvd!('S','S',matrix(plunge_op * operator(problem) * W))
        S = USV[2]
        limit = 10^(3/4*log10(eps(numtype(frequency_basis(problem)))))
        maxind = findlast(S.>limit)
        Sinv = 1./S[1:maxind]
        b = zeros(size(dest(plunge_op)))
        y = zeros(size(USV[3],1))
        x1 = zeros(size(src(operator(problem))))
        x2 = zeros(size(src(operator(problem))))
        sy = zeros(maxind,)
        new(problem, plunge_op, W, USV[1][:,1:maxind]',USV[3][1:maxind,:]'*diagm(Sinv[:]),b,y,sy,x1,x2)
    end
end

eltype{ELT,SRC,DEST}(::Type{FE_ProjectionSolver{ELT,SRC,DEST}}) = ELT

function FE_ProjectionSolver(problem::FE_DiscreteProblem)
    ELT = eltype(problem)
    SRC = typeof(time_basis_restricted(problem))
    DEST = typeof(frequency_basis(problem))
    FE_ProjectionSolver{ELT,SRC,DEST}(problem)
end


FE_ProjectionSolver(p::FE_TensorProductProblem) = TensorProductOperator(map(FE_ProjectionSolver,p.problems)...)

function plunge_operator(problem::FE_DiscreteProblem)
    A = operator(problem)
    Ap = operator_transpose(problem)
    I = IdentityOperator(time_basis_restricted(problem))

    A*Ap - I
end
 
estimate_plunge_rank{N}(problem::FE_DiscreteProblem{N}) = min(round(Int, 9*log(param_N(problem))*(param_M(problem)*param_N(problem)/param_L(problem))^(1-1/N) + 2),param_N(problem))

estimate_plunge_rank(problem::FE_DiscreteProblem{1,BigFloat}) = round(Int, 28*log(param_N(problem)) + 5)

function apply!(s::FE_ProjectionSolver, dest, src, coef_dest, coef_src)
    A = operator(s)
    At = operator_transpose(s)
    P = s.plunge_op
    apply!(P,s.b, coef_src)
    A_mul_B!(s.sy, s.Ut, s.b)
    A_mul_B!(s.y, s.VS, s.sy)
    apply!(s.W, s.x2, s.y)
    #x2 = reshape(s.W * y,size(src(A)))
    apply!(A, s.b, s.x2)
    apply!(At, s.x1, coef_src-s.b)
    for i = 1:length(coef_dest)
        coef_dest[i] = s.x1[i] + s.x2[i]
    end
    apply!(normalization(problem(s)), coef_dest)
end


