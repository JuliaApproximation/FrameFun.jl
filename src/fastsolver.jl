# fastsolver.jl

"""
A fast FE solver based on a low-rank approximation of the plunge region. The plunge region
is isolated using a projection operator.
For more details, see the paper 'Fast algorithms for the computation of Fourier extensions of arbitrary length'
http://arxiv.org/abs/1509.00206
"""
immutable FE_ProjectionSolver{ELT} <: FE_Solver{ELT}
    problem     ::  FE_DiscreteProblem
    plunge_op   ::  AbstractOperator    # store the operator because it allocates memory
    W           ::  MultiplicationOperator
    Ut          ::  Array{ELT,2}
    VS          ::  Array{ELT,2}
    b
    blinear     ::  Array{ELT,1}
    y           ::  Array{ELT,1}
    sy          ::  Array{ELT,1}
    x2
    x1

    function FE_ProjectionSolver(problem::FE_DiscreteProblem; cutoff = default_cutoff(problem), R = 5, verbose = false, options...)
        plunge_op = plunge_operator(problem)
        finished=false
        USV = ()
        sold = 0.6
        snew = 0.6
        Rold = 1
        while R<=param_N(problem)
            random_matrix = map(ELT, rand(param_N(problem), R))
            Wsrc = ELT <: Complex ? Cn{ELT}(size(random_matrix,2)) : Rn{ELT}(size(random_matrix,2))
            Wdest = src(operator(problem))
            W = MatrixOperator(Wsrc, Wdest, random_matrix)
        
            USV = LAPACK.gesdd!('S',matrix(plunge_op * operator(problem) * W))
            verbose && println("minimal singular value: ",minimum(USV[2])," cutoff : ",cutoff)
            verbose && println("R/Rest/Rmax: ",R,"/",estimate_plunge_rank(problem),"/",param_N(problem))
            if minimum(USV[2])<cutoff || R==param_N(problem)
                S = USV[2]
                maxind = findlast(S.>cutoff)
                Sinv = 1./S[1:maxind]
                b = zeros(ELT, dest(plunge_op))
                blinear = zeros(ELT, length(dest(plunge_op)))
                y = zeros(ELT, size(USV[3],1))
                x1 = zeros(ELT, src(operator(problem)))
                x2 = zeros(ELT, src(operator(problem)))
                sy = zeros(ELT, maxind)
                return new(problem, plunge_op, W, USV[1][:,1:maxind]',USV[3][1:maxind,:]'*diagm(Sinv[:]),b,blinear,y,sy,x1,x2)
            else
                # Assuming exponential decay, take the optimal step
                snew = minimum(USV[2])
                alpha = (log(snew)-log(sold))/(R-Rold)
                step = (log(cutoff)-log(snew))/alpha
                sold = snew
                Rold = R
                if (sold<10.0^(-3))
                    #started converging
                    R = min(round(Int,R+step),param_N(problem))
                else
                    #double
                    R = min(round(Int,2*R),param_N(problem))
                end
            end
        end
    end
end

FE_ProjectionSolver(problem::FE_DiscreteProblem; options...) =
    FE_ProjectionSolver{eltype(problem)}(problem; options...)


function plunge_operator(problem::FE_DiscreteProblem)
    A = operator(problem)
    Ap = operator_transpose(problem)
    I = IdentityOperator(time_basis_restricted(problem))

    A*Ap - I
end

default_cutoff(problem::FE_DiscreteProblem) = 10^(4/5*log10(eps(numtype(frequency_basis(problem)))))
estimate_plunge_rank{N}(problem::FE_DiscreteProblem{N}) = min(round(Int, 9*log(param_N(problem))*(param_M(problem)*param_N(problem)/param_L(problem))^(1-1/N) + 2),param_N(problem))

estimate_plunge_rank(problem::FE_DiscreteProblem{1,BigFloat}) = round(Int, 28*log(param_N(problem)) + 5)

apply!(s::FE_ProjectionSolver, dest, src, coef_dest, coef_src) =
    apply!(s, dest, src, coef_dest, coef_src, operator(s), operator_transpose(s), s.plunge_op, s.W, s.x1, s.x2)

function apply!(s::FE_ProjectionSolver, destset, srcset, coef_dest, coef_src, A, At, P, W, x1, x2)
    apply!(P, s.b, coef_src)
    BasisFunctions.linearize_coefficients!(dest(A), s.blinear, s.b)
    A_mul_B!(s.sy, s.Ut, s.blinear)
    A_mul_B!(s.y, s.VS, s.sy)
    apply!(W, s.x2, s.y)
    #x2 = reshape(s.W * y,size(src(A)))
    apply!(A, s.b, x2)
    apply!(At, x1, coef_src-s.b)
    for i in eachindex(x1)
        x1[i] += x2[i]
    end
    apply!(normalization(problem(s)), coef_dest, x1)
end
