
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
    plunge_op   ::  Union{Void,AbstractOperator} # (A*Zt-I), store because it allocates memory
    b           ::  Array{ELT}                  # Scratch for right hand size
    blinear     ::  Array{ELT,1}     # Scratch for linearized right hand side (necessary for svd inproducts)
    x2          ::  Array{ELT}
    x1          ::  Array{ELT}

    function AZSolver{ELT}(trunc::FE_Solver, A::AbstractOperator, Zt::AbstractOperator; plunge_op = nothing, options...) where ELT
        # Allocate scratch space
        b = zeros(src(trunc))
        blinear = zeros(ELT, length(src(trunc)))
        x1 = zeros(src(A))
        x2 = zeros(src(A))
        # Construct AZSolver object
        new(trunc, A, Zt, plunge_op, b, blinear, x1, x2)
    end
end

function AZSolver(A::AbstractOperator{ELT}, Zt::AbstractOperator{ELT}, util_operators::AbstractOperator{ELT}...;
        use_plunge=true, options...) where {ELT}
    plunge_op = nothing
    OP = A
    if use_plunge
        # Calculate (A*Zt-I)
        plunge_op = plunge_operator(A, Zt)
        OP = plunge_op*A
    end
    _AZSolver(A::AbstractOperator, Zt::AbstractOperator, OP, util_operators...; plunge_op=plunge_op, options...)
end

function _AZSolver(A::AbstractOperator{ELT}, Zt::AbstractOperator{ELT}, OP::AbstractOperator{ELT}, util_operators::AbstractOperator{ELT}...;
        cutoff = default_cutoff(A), TRUNC=TruncatedSvdSolver, R = estimate_plunge_rank(A), verbose=false, options...) where {ELT}
    # Calculate low rank INVERSE of (A*Zt-I)*A
    TS = TRUNC(OP, util_operators...; cutoff=cutoff, R=R, verbose=verbose, options...)
    AZSolver(TS, A, Zt; options...)
end

# Set type of scratch space based on operator eltype.
AZSolver(trunc::FE_Solver{ELT}, A::AbstractOperator{ELT}, Zt::AbstractOperator{ELT}; options...) where {ELT} =
    AZSolver{eltype(A)}(trunc, A, Zt; options...)

AZSSolver(A::AbstractOperator, Zt::AbstractOperator, RD::AbstractOperator, EF::AbstractOperator;
        TRUNC = RestrictionSolver, use_plunge=false, options...) =
    AZSolver(A, Zt, RD, EF; use_plunge=use_plunge, TRUNC=TRUNC, options...)

AZSDCSolver(A::AbstractOperator, Zt::AbstractOperator,
        RD0::AbstractOperator, RD1::AbstractOperator, RD2::AbstractOperator,
        EF0::AbstractOperator, EF1::AbstractOperator, EF2::AbstractOperator;
        TRUNC = DivideAndConquerSolver, use_plunge=false, options...) =
    AZSolver(A, Zt, RD0, RD1, RD2, EF0, EF1, EF2; use_plunge=use_plunge, TRUNC=TRUNC, options...)

# If no Zt is supplied, Zt=A' (up to scaling) by default.
AZSolver(A::AbstractOperator{ELT}; scaling=nothing, options...) where {ELT} =
    AZSolver(A, ELT(1)/ELT(scaling)*A'; options...)

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
    if typeof(s.plunge_op) <: Void
        # Solve x2 = A*x=b
        apply!(s.TS,s.x2,coef_src)
    else
        # Consruct (A*Zt-I)*b
        apply!(s.plunge_op, s.b, coef_src)
        # Solve x2 = ((A*Zt-I)*A)^-1(A*Zt-I)*b
        apply!(s.TS,s.x2,s.b)
    end
    # Step 2:
    # Store A*x2 in s.b
    apply!(s.A, s.b, s.x2)
    # Store b-A*x2 in s.b
    for i in 1:length(s.b)
        s.b[i] = coef_src[i] - s.b[i]
    end
    # s.b .= coef_src .- s.b
    # Compute x1 =  Zt*(b-A*x2)
    apply!(s.Zt, s.x1, s.b)
    # Step 3:
    # x = x1 + x2
    for i in 1:length(s.x1)
        coef_dest[i] = s.x1[i] + s.x2[i]
    end
    # coef_dest .= s.x1 .+ s.x2
end

function AZSolver(platform::BasisFunctions.Platform, i; options...)
    a = A(platform, i)
    zt = Zt(platform, i)
    AZSolver(a, zt; options...)
end

function AZSSolver(platform::BasisFunctions.Platform, i; options...)
    frame_restriction, grid_restriction = FrameFun.spline_util_restriction_operators(platform, i)
    AZSSolver(A(platform, i), Zt(platform, i), frame_restriction', grid_restriction; options...)
end


function AZSDCSolver(fplatform::BasisFunctions.Platform, i, dim::Int, range; options...)
    platform = fplatform.super_platform
    # The grid on Gamma
    gamma = grid(sampler(platform, i))
    # The grid on Omega
    omega = grid(sampler(fplatform, i))

    dom = domain(primal(fplatform, i))
    p = primal(platform, i)

    ops = FrameFun.divide_and_conquer_restriction_operators(omega, gamma, p, dom, dim, range)

    if length(ops) == 2
        frame_restriction, grid_restriction = ops
        warn("going over to AZSSolver")
        AZSSolver(A(fplatform, i), Zt(fplatform, i), frame_restriction', grid_restriction; options...)
    else
        frame_res0, frame_res1, frame_res2, grid_res0, grid_res1, grid_res2 = ops
        FrameFun.AZSDCSolver(A(fplatform, i), Zt(fplatform, i), frame_res0', frame_res1', frame_res2', grid_res0, grid_res1, grid_res2)
    end
end

# Function with equal functionality, but allocating memory
function az_solve(b, A::AbstractOperator, Zt::AbstractOperator, util_operators::AbstractOperator...; use_plunge=true,
        cutoff = default_cutoff(A), trunc = truncatedsvd_solve, verbose=false, options...)
    op = A
    if use_plunge
        P = plunge_operator(A, Zt)
        x2 = trunc(P*b, P*A, util_operators...; cutoff=cutoff, verbose=verbose, options...)
    else
        x2 = trunc(b, A, util_operators... ; cutoff=cutoff, verbose=verbose, options...)
    end
    x1 = Zt*(b-A*x2)
    x1 + x2
end

# Function with equal functionality, but allocating memory
azs_solve(b, A::AbstractOperator, Zt::AbstractOperator, RD::AbstractOperator, SB::AbstractOperator;
        trunc = restriction_solve, use_plunge=false, options...) =
    az_solve(b, A, Zt, RD, SB; trunc=trunc, use_plunge=use_plunge, options...)

# Function with equal functionality, but allocating memory
azsdc_solve(b, A::AbstractOperator, Zt::AbstractOperator,
        RD0::AbstractOperator, RD1::AbstractOperator, RD2::AbstractOperator,
        EF0::AbstractOperator, EF1::AbstractOperator, EF2::AbstractOperator;
        trunc = divideandconqer_solve, use_plunge=false, options...) =
    az_solve(b, A, Zt, RD0, RD1, RD2, EF0, EF1, EF2; trunc=trunc, use_plunge=use_plunge, options...)


azsdcN_solve(b, A::AbstractOperator, Zt::AbstractOperator,
        A0::Vector{OP1}, GR0::Vector{OP2},
        A1::Vector{OP3}, GR1::Vector{OP4};
        trunc = divideandconqerN_solve, use_plunge=false, options...) where {OP1<:AbstractOperator, OP2<:AbstractOperator, OP3<:AbstractOperator, OP4<:AbstractOperator}=
    az_solve(b, A, Zt, A0..., GR0..., A1..., GR1...;
    divideandconqerN_op_lengths=[length(A0), length(GR0), length(A1), length(GR1)],
    trunc=trunc, use_plunge=use_plunge, options...)

az_tree_solve(b, A::AbstractOperator, Zt::AbstractOperator, solver::DomainDecompositionSolver;
        trunc=domaindecomposition_solve, use_plunge=false, options...) =
    az_solve(b, A, Zt, solver; use_plunge=use_plunge ,trunc=trunc)


azsdcN_solve(b, A::AbstractOperator, Zt::AbstractOperator,
        A0::Vector{OP1}, GR0::Vector{OP2},
        A1::Vector{OP3}, GR1::Vector{OP4},
        A2::Vector{OP5}, GR2::Vector{OP6},
        A3::Vector{OP7}, GR3::Vector{OP8};
        trunc = divideandconqerN_solve, use_plunge=false, options...) where {OP1<:AbstractOperator, OP2<:AbstractOperator, OP3<:AbstractOperator, OP4<:AbstractOperator,  OP5<:AbstractOperator, OP6<:AbstractOperator, OP7<:AbstractOperator, OP8<:AbstractOperator}=
    az_solve(b, A, Zt, A0..., GR0..., A1..., GR1..., A2..., GR2..., A3..., GR3...;
    divideandconqerN_op_lengths=[length(A0), length(GR0), length(A1), length(GR1), length(A2), length(GR2), length(A3), length(GR3)],
    trunc=trunc, use_plunge=use_plunge, options...)


function az_solve(platform::BasisFunctions.Platform, i, f::Function; R=0, options...)
    a = A(platform, i)
    zt = Zt(platform, i)
    s = sampler(platform, i)
    (R == 0) && (R=estimate_plunge_rank(a))
    az_solve(s*f, a, zt; R=R, options...)
end

function azs_solve(platform::BasisFunctions.Platform, i, f::Function; options...)
    a = A(platform, i)
    zt = Zt(platform, i)
    s = sampler(platform, i)
    rd,sb = spline_util_restriction_operators(platform, i)
    azs_solve(s*f, a, zt, rd', sb; options...)
end

function azsdc_solve(fplatform::BasisFunctions.Platform, i, f::Function, dim::Int, range; options...)
    platform = fplatform.super_platform
    a = A(fplatform, i)
    zt = Zt(fplatform, i)
    S = sampler(fplatform, i)

    # The grid on Gamma
    gamma = grid(sampler(platform, i))
    # The grid on Omega
    omega = grid(S)

    dom = domain(primal(fplatform, i))
    p = primal(platform, i)

    ops = FrameFun.divide_and_conquer_restriction_operators(omega, gamma, p, dom, dim, range; options...)

    if length(ops) == 2
        frame_res, grid_res = ops
        warn("going over to azs_solve")
        azs_solve(S*f, a, zt, frame_res', grid_res; options...)
    else
        frame_res0, frame_res1, frame_res2, grid_res0, grid_res1, grid_res2 = ops
        azsdc_solve(S*f, a, zt, frame_res0', frame_res1', frame_res2', grid_res0, grid_res1, grid_res2; options...)
    end
end

function azsdcN_solve(fplatform::BasisFunctions.Platform, i, f::Function; recur=nothing, options...)
    platform = fplatform.super_platform
    a = A(fplatform, i)
    zt = Zt(fplatform, i)
    S = sampler(fplatform, i)
    D = dimension(src(a))
    (recur==nothing) && (recur = D==3? 2:1)

    ops = FrameFun.divide_and_conquer_N_util_operators(fplatform, i; recur=recur, options...)
    azsdcN_solve(S*f, a, zt, ops...; options...)
end

function az_tree_solve(fplatform::BasisFunctions.Platform, i, f::Function;
        options...)
    platform = fplatform.super_platform
    a = A(fplatform, i)
    zt = Zt(fplatform, i)
    S = sampler(fplatform, i)

    # The grid on Gamma
    gamma = grid(sampler(platform, i))
    # The grid on Omega
    omega = grid(S)

    dom = domain(primal(fplatform, i))
    basis = primal(platform, i)

    solver = DomainDecompositionSolver(basis, gamma, omega, dom; options...)
    az_tree_solve(S*f, a, zt, solver; options...)
end

# max_alloc(platform, i, fun) = *(max_system_size(platform, i, fun)...)*sizeof(Float64)/(1024)^3
# function max_system_size(platform, i, fun::typeof(FrameFun.azs_solve))
#     rd,sb = FrameFun.spline_util_restriction_operators(platform, i)
#     ((size(sb,1), size(rd,1)))
# end
# function max_system_size(platform, i, fun::typeof(FrameFun.az_tree_solve))
#     tree = FrameFun.DomainDecompositionSolver(platform, i).tree
#     cs = FrameFun.leave_containers(tree)
#     a = FrameFun.DomainDecompositionLeaf[]
#     for c in cs
#         for i in c.container
#             push!(a, i)
#         end
#     end
#     maximum(map(size, map(FrameFun.DMZ,a)))
# end
