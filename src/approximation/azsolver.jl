
# fastsolver.jl

"""
Fast implementation of the AZ Algorithm. Returns an operator that is approximately pinv(A), based on a suitable Zt so that (A*Zt-I)*A is low rank.

Steps when applying to right hand side b:

1. (A*Zt-I)*A*x2=(A*Zt-I)*b, solve for x2
This is by default done using a low rank decomposition of (A*Zt-I)*A (TruncatedSvdSolver).

2. x1 = Zt*(b-A*x2)
Can be done fast if Zt and A are fast.

3. x = x1+x2
This is the solution.

"""
struct AZSolver{ELT} <: FE_Solver{ELT}
    TS          ::  AbstractOperator # The low rank decomposition of (A*Zt-I)*A
    A           ::  AbstractOperator # Store for application in step 2
    Zt          ::  AbstractOperator # Store for application in step 2
    plunge_op   ::  AbstractOperator # (A*Zt-I), store because it allocates memory
    b                                # Scratch for right hand size
    blinear     ::  Array{ELT,1}     # Scratch for linearized right hand side (necessary for svd inproducts)
    x2
    x1
    function AZSolver{ELT}(A::AbstractOperator, Zt::AbstractOperator;
            cutoff = default_cutoff(A), trunc = TruncatedSvdSolver, R = estimate_plunge_rank(A), verbose=false,options...) where ELT
        # Calculate (A*Zt-I)
        plunge_op = plunge_operator(A, Zt)
        # Calculate low rank INVERSE of (A*Zt-I)*A
        TS = trunc(plunge_op*A; cutoff=cutoff, R=R, verbose=verbose, options...)
        # Allocate scratch space
        b = zeros(dest(plunge_op))
        blinear = zeros(ELT, length(dest(plunge_op)))
        x1 = zeros(src(A))
        x2 = zeros(src(A))
        # Construct AZSolver object
        new(TS, A, Zt, plunge_op, b,blinear,x1,x2)
    end

    function AZSolver{ELT}(A::AbstractOperator, Zt::AbstractOperator, RD::AbstractOperator,  SB::IndexRestrictionOperator;
            cutoff = default_cutoff(A), trunc = RestrictionSolver, verbose=false, options...) where ELT
        # Calculate (A*Zt-I)
        plunge_op = plunge_operator(A, Zt)
        # Calculate low rank INVERSE of (A*Zt-I)*A
        TS = trunc(plunge_op*A, RD, SB; cutoff=cutoff, verbose=verbose, options...)
        # Allocate scratch space
        b = zeros(dest(plunge_op))
        blinear = zeros(ELT, length(dest(plunge_op)))
        x1 = zeros(src(A))
        x2 = zeros(src(A))
        # Construct AZSolver object
        new(TS, A, Zt, plunge_op, b,blinear,x1,x2)
    end
end

# Set type of scratch space based on operator eltype.
AZSolver(A::AbstractOperator, Zt::AbstractOperator; options...) =
    AZSolver{eltype(A)}(A, Zt; options...)

AZSolver(A::AbstractOperator, Zt::AbstractOperator, RD::AbstractOperator, SB::AbstractOperator; options...) =
    AZSolver{eltype(A)}(A, Zt, RD, SB; options...)

# If no Zt is supplied, Zt=A' (up to scaling) by default.
AZSolver(A::AbstractOperator; scaling=nothing, options...) =
    AZSolver{eltype(A)}(A, 1/scaling*A'; options...)

function plunge_operator(A, Zt)
    I = IdentityOperator(dest(A))
    A*Zt - I
end

default_cutoff(A::AbstractOperator) = 10^(4/5*log10(eps(real(eltype(A)))))

# Estimate for the rank of (A*Zt-I)*A when computing the low rank decomposition. If check fails, rank estimate is steadily increased.
@inline estimate_plunge_rank(A::AbstractOperator) =
    estimate_plunge_rank(dictionary(src(A)), dictionary(dest(A)))

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
        return min(round(Int, 9*log(nml)),length(src))
    else
        return min(round(Int, 9*log(nml)*nml^((N-1)/N)),length(src))
    end
end

function apply!(s::AZSolver, destset, srcset, coef_dest, coef_src)
    # Step 1:
    # Consruct (A*Zt-I)*b
    apply!(s.plunge_op, s.b, coef_src)
    BasisFunctions.linearize_coefficients!(dest(s.A), s.blinear, s.b)
    # Solve x2 = ((A*Zt-I)*A)^-1(A*Zt-I)*b
    apply!(s.TS,s.x2,s.blinear)
    # Step 2:
    # Store A*x2 in b
    apply!(s.A, s.b, s.x2)
    # Compute x1 =  Zt*(b-A*x2)
    apply!(s.Zt, s.x1, coef_src-s.b)
    # Step 3:
    # x = x1 + x2
    for i in eachindex(s.x1)
        coef_dest[i] = s.x1[i]+s.x2[i]
    end
end


# Function with equal functionality, but allocating memory
function az_solve(A::AbstractOperator, Zt::AbstractOperator, b;
        cutoff = default_cutoff(A), trunc = truncatedsvd_solve, verbose=false, options...)
    P = plunge_operator(A, Zt)
    x2 = trunc(P*A, P*b; cutoff=cutoff, verbose=verbose, options...)
    x1 = Zt*(b-A*x2)
    x1 + x2
end

# Function with equal functionality, but allocating memory
function az_solve(A::AbstractOperator, Zt::AbstractOperator, RD::AbstractOperator, SB::AbstractOperator, b;
        cutoff = default_cutoff(A), trunc = restriction_solve, verbose=false, options...)
    P = plunge_operator(A, Zt)
    x2 = trunc(P*A, RD, SB, P*b; cutoff=cutoff, verbose=verbose, options...)
    x1 = Zt*(b-A*x2)
    x1 + x2
end
