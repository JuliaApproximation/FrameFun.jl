
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
struct AZSolver{ELT,AFIRST<:Union{Val{false},Val{true}}} <: FE_Solver{ELT}
    TS          ::  DictionaryOperator # The low rank decomposition of (A*Zt-I)*A
    A           ::  DictionaryOperator # Store for application in step 2
    Zt          ::  DictionaryOperator # Store for application in step 2
    plunge_op   ::  Union{Void,DictionaryOperator} # (A*Zt-I), store because it allocates memory
    b                                # Scratch for right hand size
    blinear     ::  Array{ELT,1}     # Scratch for linearized right hand side (necessary for svd inproducts)
    x2
    x1

    function AZSolver{ELT,AFIRST}(trunc::FE_Solver, A::DictionaryOperator, Zt::DictionaryOperator, plunge_op::Union{Void,DictionaryOperator};
            options...) where {ELT,AFIRST}
        # Allocate scratch space
        b = zeros(src(trunc))
        blinear = zeros(ELT, length(src(trunc)))
        x1 = zeros(src(A))
        x2 = zeros(src(A))
        new(trunc, A, Zt, plunge_op, b, blinear, x1, x2)
    end
end

# Set type of scratch space based on operator eltype.
AZSolver(trunc::FE_Solver{ELT}, A::DictionaryOperator{ELT}, Zt::DictionaryOperator{ELT}, plunge_op::Union{Void,DictionaryOperator}=nothing;
        afirst=true, options...) where {ELT} =
    AZSolver{eltype(A),Val{afirst}}(trunc, A, Zt, plunge_op; options...)

# If no Zt is supplied, Zt=A' (up to scaling) by default.
AZSolver(A::DictionaryOperator{ELT}; scaling=nothing, options...) where {ELT} =
    AZSolver(A, ELT(1)/ELT(scaling)*A'; options...)

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

function AZSSolver(A::DictionaryOperator, Zt::DictionaryOperator, RD::DictionaryOperator, EF::DictionaryOperator; options...)
    TS = RestrictionSolver(A, RD, EF; options...)
    AZSolver(TS, A, Zt; options...)
end

function AZSDCSolver(A::DictionaryOperator, Zt::DictionaryOperator,
        RD0::DictionaryOperator, RD1::DictionaryOperator, RD2::DictionaryOperator,
        EF0::DictionaryOperator, EF1::DictionaryOperator, EF2::DictionaryOperator;
        options...)
    TS = DivideAndConquerSolver(A, RD0, RD1, RD2, EF0, EF1, EF2; options...)
    AZSolver(TS, A, Zt; options...)
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
        return min(round(Int, 9*log(nml)),length(src))
    else
        return min(round(Int, 9*log(nml)*nml^((N-1)/N)),length(src))
    end
end

apply!(s::AZSolver, coef_dest, coef_src) = _apply!(s, coef_dest, coef_src,
        s.plunge_op, s.A, s.Zt, s.b, s.blinear, s.TS, s.x1, s.x2)

function _apply!(s::AZSolver{ELT,Val{true}}, coef_dest, coef_src, plunge_op::DictionaryOperator, A, Zt, b, blinear, TS, x1, x2) where {ELT}
    # Step 1:
    # Consruct (A*Zt-I)*b
    apply!(plunge_op, b, coef_src)
    # Solve x2 = ((A*Zt-I)*A)^-1(A*Zt-I)*b
    apply!(TS,x2,b)
    _apply2!(s, coef_dest, coef_src, A, Zt, b, x1, x2)
end

function _apply!(s::AZSolver{ELT,Val{true}}, coef_dest, coef_src, plunge_op::Void, A, Zt, b, blinear, TS, x1, x2) where {ELT}
    # Step 1:
    # Solve x2 = A*x=b
    apply!(TS,x2,coef_src)
    _apply2!(s, coef_dest, coef_src, A, Zt, b, x1, x2)
end

function _apply2!(s::AZSolver{ELT,Val{true}}, coef_dest, coef_src, A, Zt, b, x1, x2) where {ELT}
    # Step 2:
    # Store A*x2 in s.b
    apply!(A, b, x2)
    # Store b-A*x2 in s.b
    # Compute x1 =  Zt*(b-A*x2)
    # b .= coef_src .- b
    for i in eachindex(b)
        b[i] = coef_src[i] - b[i]
    end
    apply!(Zt, x1, b)
    # Step 3:
    # x = x1 + x2
    # coef_dest .= s.x1 .+ s.x2
    for i in eachindex(x1)
        coef_dest[i] = x1[i] + x2[i]
    end
end

function _apply!(s::AZSolver{ELT,Val{false}}, coef_dest, coef_src, plunge_op, A, Zt, b, blinear, TS, x1, x2) where {ELT}
    # Step 1:
    # Compute x2 =  Zt*b
    apply!(Zt, x2, coef_src)
    _apply2!(s, coef_dest, coef_src, plunge_op, A, Zt, b, blinear, TS, x1, x2)
end

function _apply2!(s::AZSolver{ELT,Val{false}}, coef_dest, coef_src, plunge_op::Void, A, Zt, b, blinear, TS, x1, x2) where {ELT}
    # Step 2:
    # Store A*x2 in s.b
    apply!(A, b, x2)
    # Store b-A*x2 in s.b
    for i in eachindex(b)
        b[i] = coef_src[i] - b[i]
    end
    # Solve x2 = A*x=(b-A*x2)
    apply!(TS,x1,b)
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

function AZSSolver(fplatform::BasisFunctions.Platform, i; options...)
    a = BasisFunctions.A(fplatform, i; options...)
    zt = BasisFunctions.Zt(fplatform, i; options...)
    platform = fplatform.super_platform
    s = BasisFunctions.sampler(fplatform, i)
    omega = grid(s)
    gamma = supergrid(omega)
    domain = FrameFun.domain(src(a))
    frame_restriction, grid_restriction = azselection_restriction_operators(primal(platform, i), gamma, omega, domain)
    # frame_restriction, grid_restriction = FrameFun.spline_util_restriction_operators(platform, i)
    AZSSolver(a, zt, frame_restriction', grid_restriction; options...)
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
function az_solve(b, A::DictionaryOperator, Zt::DictionaryOperator, util_operators::DictionaryOperator...; use_plunge=true,
        cutoff = default_cutoff(A), trunc = truncatedsvd_solve, verbose=false, options...)
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
azs_solve(b, A::DictionaryOperator, Zt::DictionaryOperator, RD::DictionaryOperator, SB::DictionaryOperator;
        trunc = restriction_solve, use_plunge=false, options...) =
    az_solve(b, A, Zt, RD, SB; trunc=trunc, use_plunge=use_plunge, options...)

# Function with equal functionality, but allocating memory
azsdc_solve(b, A::DictionaryOperator, Zt::DictionaryOperator,
        RD0::DictionaryOperator, RD1::DictionaryOperator, RD2::DictionaryOperator,
        EF0::DictionaryOperator, EF1::DictionaryOperator, EF2::DictionaryOperator;
        trunc = divideandconqer_solve, use_plunge=false, options...) =
    az_solve(b, A, Zt, RD0, RD1, RD2, EF0, EF1, EF2; trunc=trunc, use_plunge=use_plunge, options...)


azsdcN_solve(b, A::DictionaryOperator, Zt::DictionaryOperator,
        A0::Vector{OP1}, GR0::Vector{OP2},
        A1::Vector{OP3}, GR1::Vector{OP4};
        trunc = divideandconqerN_solve, use_plunge=false, options...) where {OP1<:DictionaryOperator, OP2<:DictionaryOperator, OP3<:DictionaryOperator, OP4<:DictionaryOperator}=
    az_solve(b, A, Zt, A0..., GR0..., A1..., GR1...;
    divideandconqerN_op_lengths=[length(A0), length(GR0), length(A1), length(GR1)],
    trunc=trunc, use_plunge=use_plunge, options...)

az_tree_solve(b, A::DictionaryOperator, Zt::DictionaryOperator, solver::DomainDecompositionSolver;
        trunc=domaindecomposition_solve, use_plunge=false, options...) =
    az_solve(b, A, Zt, solver; use_plunge=use_plunge ,trunc=trunc)


azsdcN_solve(b, A::DictionaryOperator, Zt::DictionaryOperator,
        A0::Vector{OP1}, GR0::Vector{OP2},
        A1::Vector{OP3}, GR1::Vector{OP4},
        A2::Vector{OP5}, GR2::Vector{OP6},
        A3::Vector{OP7}, GR3::Vector{OP8};
        trunc = divideandconqerN_solve, use_plunge=false, options...) where {OP1<:DictionaryOperator, OP2<:DictionaryOperator, OP3<:DictionaryOperator, OP4<:DictionaryOperator,  OP5<:DictionaryOperator, OP6<:DictionaryOperator, OP7<:DictionaryOperator, OP8<:DictionaryOperator}=
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

    solver = DomainDecompositionSolver(basis, gamma, omega, dom, a; options...)
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
