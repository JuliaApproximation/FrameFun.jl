
# fastsolver.jl

"""
Fast implementation of the AZ Algorithm. Returns an operator that is approximately
`pinv(A)`, based on a suitable `Zt` so that `(A*Zt-I)*A` is low rank.

Steps when applying to right hand side b:

1. (A*Zt-I)*A*x2=(A*Zt-I)*b, solve for x2
This is by default done using a low rank decomposition of (A*Zt-I)*A (TruncatedSvdSolver).

2. x1 = Zt*(b-A*x2)
Can be done fast if Zt and A are fast.

3. x = x1+x2
This is the solution.
"""
struct AZSolver{ELT} <: AbstractSolverOperator{ELT}
    TS          ::  DictionaryOperator # The low rank decomposition of (A*Zt-I)*A
    A           ::  DictionaryOperator # Store for application in step 2
    Zt          ::  DictionaryOperator # Store for application in step 2
    plunge_op   ::  DictionaryOperator # (A*Zt-I), store because it allocates memory
    b                                # Scratch for right hand size
    blinear     ::  Array{ELT,1}     # Scratch for linearized right hand side (necessary for svd inproducts)
    x2
    x1

    function AZSolver{ELT}(trunc::AbstractSolverOperator, A::DictionaryOperator, Zt::DictionaryOperator, plunge_op::DictionaryOperator;
            options...) where ELT
        # Allocate scratch space
        b = zeros(src(trunc))
        blinear = zeros(ELT, length(src(trunc)))
        x1 = zeros(src(A))
        x2 = zeros(src(A))
        new(trunc, A, Zt, plunge_op, b, blinear, x1, x2)
    end
end

# Set type of scratch space based on operator eltype.
AZSolver(trunc::AbstractSolverOperator{ELT}, A::DictionaryOperator{ELT}, Zt::DictionaryOperator{ELT}, plunge_op::DictionaryOperator;
        options...) where {ELT} =
    AZSolver{eltype(A)}(trunc, A, Zt, plunge_op; options...)

# If no Zt is supplied, Zt=A' (up to an optional scaling) by default.
function AZSolver(A::DictionaryOperator{ELT}; scaling=nothing, options...) where {ELT}
    if scaling == nothing
        AZSolver(A, A'; options...)
    else
        AZSolver(A, ELT(1)/ELT(scaling)*A'; options...)
    end
end

operator(op::AZSolver) = op.A

function plunge_operator(A, Zt)
    I = IdentityOperator(dest(A))
    A*Zt - I
end

function AZSolver(A::DictionaryOperator{ELT}, Zt::DictionaryOperator{ELT};
        cutoff = default_cutoff(A), TRUNC=TruncatedSvdSolver, R = estimate_plunge_rank(A), verbose=false, options...) where {ELT}
    # Calculate (A*Zt-I)
    plunge_op = plunge_operator(A, Zt)
    TS = TRUNC(plunge_op*A; cutoff=cutoff, R=R, verbose=verbose, options...)
    AZSolver(TS, A, Zt, plunge_op; options...)
end

default_cutoff(A::DictionaryOperator) = 10^(4/5*log10(eps(real(eltype(A)))))

# Estimate for the rank of (A*Zt-I)*A when computing the low rank decomposition. If check fails, rank estimate is steadily increased.
@inline estimate_plunge_rank(A::DictionaryOperator) =
    estimate_plunge_rank(src(A), dest(A))

@inline estimate_plunge_rank(src::ExtensionFrame, dest::Dictionary) =
    estimate_plunge_rank(superdict(src), domain(src), dest)

@inline estimate_plunge_rank(src::Dictionary, dest::Dictionary) =
    default_estimate_plunge_rank(src, dest)

@inline estimate_plunge_rank(src::Dictionary, domain::Domain, dest::Dictionary) =
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
        s.plunge_op, s.A, s.Zt, s.b, s.blinear, s.TS, s.x1, s.x2)

function _apply!(s::AZSolver{ELT}, coef_dest, coef_src, plunge_op::DictionaryOperator, A, Zt, b, blinear, TS, x1, x2) where {ELT}
    # Step 1:
    # Compute (A*Zt-I)*b
    apply!(plunge_op, b, coef_src)
    # Solve x2 = ((A*Zt-I)*A)^-1(A*Zt-I)*b
    apply!(TS,x2,b)

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

function AZSolver(platform::BasisFunctions.Platform, i; options...)
    a = A(platform, i)
    zt = Zt(platform, i)
    AZSolver(a, zt; options...)
end
