
"""
A fast FE solver based on a low-rank approximation of the plunge region. The plunge region
is isolated using a projection operator.
For more details, see the paper 'Fast algorithms for the computation of Fourier extensions of arbitrary length'
http://arxiv.org/abs/1509.00206
"""
struct FE_TridiagonalSolver{ELT} <: AbstractSolverOperator{ELT}
    Ut      :: Array{ELT,2}
    V       :: Array{ELT,2}
    Sinv    :: Array{ELT,1}
    A       ::  DictionaryOperator
    b
    blinear ::  Array{ELT,1}
    x2
    x1
    y
    y2
    scaling
    function FE_TridiagonalSolver{ELT}(A::DictionaryOperator, scaling; threshold = default_threshold(A), irange = estimate_plunge_range(A,scaling,estimate_plunge_rank(A)), verbose=false,options...) where ELT
        R = estimate_plunge_rank(A)
        finished = false
        V=[]
        U=[]
        S=[]
        while !finished
            (V,U)=mideigs(size(A,2),size(A,1),scaling,irange)
            S = zeros(ELT,size(V,2))
            for i=1:size(V,2)
                S[i]=U[:,i]'*apply(A,V[:,i])
            end
            if all(minimum(abs.(S)).< threshold) || R>size(A,1)
                finished=true
            else
                R = 2*R
                irange=estimate_plunge_range(A,scaling,R)
            end
        end
        I = abs.(S).>(threshold*maximum(abs.(S)))
        b = zeros(dest(A))
        blinear = zeros(ELT, length(dest(A)))
        x1 = zeros(src(A))
        x2 = zeros(src(A))
        y = zeros(ELT,sum(I))
        y2 = zeros(ELT,sum(I))
    new(U[:,I]',V[:,I],S[I].^(-1), A, b,blinear,x1,x2,y,y2,scaling)
    end
end

FE_TridiagonalSolver(A::DictionaryOperator; scaling = nothing, options...) =
    FE_TridiagonalSolver(A, scaling; options...)

# TODO: remove the scaling factor here
FE_TridiagonalSolver(A::DictionaryOperator, scaling; options...) =
    FE_TridiagonalSolver{eltype(A)}(A, scaling; options...)

operator(op::FE_TridiagonalSolver) = op.A

# Estimate for the rank of (A*Zt-I)*A when computing the low rank decomposition. If check fails, rank estimate is steadily increased.
# TODO: remove this estimate and use adaptivity
estimate_plunge_rank(A::DictionaryOperator) =
    estimate_plunge_rank(src(A), dest(A))

estimate_plunge_rank(src::ExtensionFrame, dest::Dictionary) =
    estimate_plunge_rank(superdict(src), support(src), dest)

estimate_plunge_rank(src::Dictionary, dest::Dictionary) =
    default_estimate_plunge_rank(src, dest)

estimate_plunge_rank(src::Dictionary, domain::Domain, dest::Dictionary) =
    default_estimate_plunge_rank(src, dest)

function default_estimate_plunge_rank(src::Dictionary, dest::Dictionary)
    nml=length(src)^2/length(dest)
    N = dimension(src)
    if N==1
        return max(1,min(round(Int, 9*log(nml)),length(src)))
    else
        return max(1,min(round(Int, 9*log(nml)*nml^((N-1)/N)),length(src)))
    end
end

function estimate_plunge_range(A,L,C)
    mid = size(A,2)-round(Int,size(A,1)*size(A,2)/L)
    mid .+ (max(1-mid,-round(Int,C)):min(size(A,2)-mid,round(Int,C)))
end

apply!(s::FE_TridiagonalSolver, coef_dest, coef_src) =
apply!(s, coef_dest, coef_src, s.A, s.A', s.x1, s.x2,s.y,s.y2)

function apply!(s::FE_TridiagonalSolver, coef_dest, coef_src, A, At, x1, x2,y,y2)
    # Applying plunge to the right hand side
    BasisFunctions.linearize_coefficients!(dest(A), s.blinear, coef_src)
    mul!(y,s.Ut,s.blinear)
    for i =1:length(y)
        y[i]*=s.Sinv[i]
        if abs(y[i])>10
            y[i]=0
        end
    end
    mul!(x2,s.V,y)
    # x2 solves the middle guys
    apply!(A, s.b, x2)
    apply!(At, x1, coef_src-s.b)
    # x1 solves the remainder through one application of the Aerator
    for i in eachindex(x1)
        coef_dest[i] = x1[i]/s.scaling+x2[i]
    end
end

function mideigs(N,M,L,irange)
    J2 = Chamzas(L,M-1,N)
    J1 = Chamzas(L,N-1,M)
    E,V1 = eigen(J1,irange)
    for i = 1:size(V1,1)
        for j = 1:size(V1,2)
            if i%2 == 0
                V1[i,j]*=-1
            end
        end
    end
    V1b = circshift(V1,round(Int,(N-1)/2)+1)
    E,V2 = eigen(J2,irange .+ (M-N))
    (V1b,V2)
end

function Chamzas(N,M,F)
    bs = sin.((pi/N)*(1:M)).*sin.((pi/N)*(M .- (0:(M-1))))
    cs = -cos.((pi/N)*(2*(0:M) .- M))*cos.((pi/N)*(F))
    T = SymTridiagonal(cs,bs)
end
