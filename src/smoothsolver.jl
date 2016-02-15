# smoothsolver.jl

"""
A fast FE solver based on a low-rank approximation of the plunge region. The plunge region
is isolated using a projection operator. This algorithm contains an extra smoothing step

"""
immutable FE_SmoothProjectionSolver{ELT,SRC,DEST} <: FE_Solver{ELT,SRC,DEST}
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
        x3          ::  Array{ELT}
        QQ          ::  Array{ELT,2}
        scale       ::  Function

        function FE_SmoothProjectionSolver(problem::FE_DiscreteProblem, scale::Function)
                plunge_op = plunge_operator(problem)
                R = estimate_plunge_rank(problem)
                WU = MatrixOperator( map(ELT, rand(param_N(problem), R)) )
                WV = MatrixOperator( map(ELT, rand(param_M(problem), R)) )
                USVU = LAPACK.gesvd!('S','S',matrix(plunge_op * operator(problem) * WU))
                USVV = LAPACK.gesvd!('S','S',matrix(operator_transpose(problem) * plunge_op * WV))
                SU = USVU[2]
                limit = 10^(4/5*log10(eps(numtype(frequency_basis(problem)))))
                maxind = findlast(SU.>limit)
                Sinv = 1./SU[1:maxind]
                D = IdxnScalingOperator(frequency_basis(problem), scale)
                AD = inv(D)
                limitv = 10^(1/2*log10(eps(numtype(frequency_basis(problem)))))
                maxindv = findlast(SU.>limitv)
                ADV = (USVV[1][:,1:maxindv])./diag(matrix(AD))
                Q, R = qr(ADV)
                b = zeros(size(dest(plunge_op)))
                y = zeros(size(USVU[3],1))
                x1 = zeros(size(src(operator(problem))))
                x2 = zeros(size(src(operator(problem))))
                x3 = zeros(size(src(operator(problem))))
                sy = zeros(maxind,)
                new(problem, plunge_op, WU, USVU[1][:,1:maxind]',USVU[3][1:maxind,:]'*diagm(Sinv[:]),b,y,sy,x1,x2,x3,Q*Q', scale)
        end
end

eltype{ELT,SRC,DEST}(::Type{FE_ProjectionSolver{ELT,SRC,DEST}}) = ELT


function FE_SmoothProjectionSolver(problem::FE_DiscreteProblem; scale=(i->1+i^3))
        ELT = eltype(problem)
        SRC = typeof(time_basis_restricted(problem))
        DEST = typeof(frequency_basis(problem))
        FE_SmoothProjectionSolver{ELT,SRC,DEST}(problem, scale)
end


FE_SmoothProjectionSolver(p::FE_TensorProductProblem) = TensorProductOperator(map(FE_SmoothProjectionSolver,p.problems)...)

@debug function apply!(s::FE_SmoothProjectionSolver, destarg, src, coef_dest, coef_src)
        A = operator(s)
        At = operator_transpose(s)
        P = s.plunge_op
        apply!(P,s.b, coef_src)
        A_mul_B!(s.sy, s.Ut, s.b)
        A_mul_B!(s.y, s.VS, s.sy)
        apply!(s.W, s.x2, s.y)
        # smoothing x2 step
        D = IdxnScalingOperator(frequency_basis(s.problem), s.scale)
        AD = inv(D)
        apply!(D,s.x1,s.x2)
        apply!(MatrixOperator(s.QQ,frequency_basis(s.problem),frequency_basis(s.problem)),s.x3,s.x1)
        for i = 1:length(s.x1)
                s.x1[i] = s.x1[i]- s.x3[i]
        end
        apply!(AD,s.x3,s.x1)
        for i = 1:length(s.x2)
                s.x2[i] = s.x2[i] - s.x3[i]
        end
        
        # post smoothing step
        apply!(A, s.b, s.x2)
        apply!(At, s.x1, coef_src-s.b)
        for i = 1:length(coef_dest)
                coef_dest[i] = s.x1[i] + s.x2[i]
        end
        apply!(normalization(problem(s)), coef_dest)
end

# For FunctionSets that have a DC component
dc_index(b::ChebyshevBasis) = 1
dc_index(b::FourierBasis) = 1
