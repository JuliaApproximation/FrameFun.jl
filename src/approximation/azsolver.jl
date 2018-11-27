
"""
Fast implementation of the AZ Algorithm. Returns an operator that is approximately
`pinv(A)`, based on a suitable `Zt` so that `(A*Zt-I)*A` is low rank.

Steps when applying to right hand side b:

1. (A*Zt-I)*A*x2=(A*Zt-I)*b, solve for x2
This is by default done using a low rank decomposition of (A*Zt-I)*A (RandomizedSvdSolver).

2. x1 = Zt*(b-A*x2)
Can be done fast if Zt and A are fast.

3. x = x1+x2
This is the solution.
"""
struct AZSolver{T} <: AbstractSolverOperator{T}
    A           ::  DictionaryOperator # Store for application in step 2
    Zt          ::  DictionaryOperator # Store for application in step 2
    plunge_op   ::  DictionaryOperator # (A*Zt-I), store because it allocates memory
    psolver     ::  DictionaryOperator # The low rank plunge solver for (A*Zt-I)*A
    b                                # Scratch for right hand size
    blinear     ::  Array{T,1}     # Scratch for linearized right hand side (necessary for svd inproducts)
    x2
    x1

    function AZSolver{T}(A::DictionaryOperator, Zt::DictionaryOperator, plunge_op::DictionaryOperator, psolver::DictionaryOperator) where {T}
        # Allocate scratch space
        b = zeros(src(psolver))
        blinear = zeros(T, length(src(psolver)))
        x1 = zeros(src(A))
        x2 = zeros(src(A))
        new(A, Zt, plunge_op, psolver, b, blinear, x1, x2)
    end
end

# Set type of scratch space based on operator eltype.
AZSolver(A::DictionaryOperator{T}, Zt::DictionaryOperator{T}, plunge_op::DictionaryOperator{T}, psolver::DictionaryOperator{T}) where {T} =
    AZSolver{eltype(A)}(A, Zt, plunge_op, psolver)

# # If no Zt is supplied, Zt=A' (up to scaling) by default.
AZSolver(A::DictionaryOperator{T}; scaling=nothing, options...) where {T} =
    AZSolver(A, T(1)/T(scaling)*A'; options...)

operator(s::AZSolver) = s.A

function plunge_operator(A, Zt)
    I = IdentityOperator(dest(A))
    A*Zt - I
end

default_regularization(A; options...) = RandomizedSvdSolver(A; options...)

function AZSolver(A::DictionaryOperator{T}, Zt::DictionaryOperator{T};
            REG = default_regularization,
            rankestimate = estimate_plunge_rank(A),
            threshold = default_threshold(A),
            options...) where {T}

    plunge_op = plunge_operator(A, Zt)
    psolver = REG(plunge_op*A; threshold = threshold, rankestimate = rankestimate, options...)
    AZSolver(A, Zt, plunge_op, psolver)
end

default_threshold(A::DictionaryOperator) = regularization_threshold(eltype(A))

# Estimate for the rank of (A*Zt-I)*A when computing the low rank decomposition. If check fails, rank estimate is steadily increased.
estimate_plunge_rank(A::DictionaryOperator) =
    estimate_plunge_rank(src(A), dest(A))

estimate_plunge_rank(src::ExtensionFrame, dest::Dictionary) =
    estimate_plunge_rank(superdict(src), domain(src), dest)

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

apply!(s::AZSolver, coef_dest, coef_src) = _apply!(s, coef_dest, coef_src,
        s.plunge_op, s.A, s.Zt, s.b, s.blinear, s.psolver, s.x1, s.x2)

function _apply!(s::AZSolver, coef_dest, coef_src, plunge_op::DictionaryOperator, A, Zt, b, blinear, psolver, x1, x2)
    # Step 1:
    # Compute (A*Zt-I)*b
    apply!(plunge_op, b, coef_src)
    # Solve x2 = ((A*Zt-I)*A)^-1(A*Zt-I)*b
    apply!(psolver, x2, b)

    # Step 2:
    # Store A*x2 in s.b and subtract from the right hand side
    apply!(A, b, x2)
    for i in eachindex(b)
        b[i] = coef_src[i] - b[i]
    end
    # Then compute x1 =  Zt*(b-A*x2)
    apply!(Zt, x1, b)
    # Step 3:
    # x = x1 + x2
    for i in eachindex(x1)
        coef_dest[i] = x1[i] + x2[i]
    end
end

function AZSolver(platform::Platform, i; options...)
    A = matrix_A(platform, i)
    Zt = matrix_Zt(platform, i)
    AZSolver(A, Zt; options...)
end
