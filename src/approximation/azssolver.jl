
"""

"""
struct AZSSolver{ELT,AFIRST<:Union{Val{false},Val{true}}} <: FE_Solver{ELT}
    TS          ::  DictionaryOperator # The low rank decomposition of (A*Zt-I)*A
    A           ::  DictionaryOperator # Store for application in step 2
    Zt          ::  DictionaryOperator # Store for application in step 2
    b                                # Scratch for right hand size
    blinear     ::  Array{ELT,1}     # Scratch for linearized right hand side (necessary for svd inproducts)
    x2
    x1

    function AZSSolver{ELT,AFIRST}(trunc::FE_Solver, A::DictionaryOperator, Zt::DictionaryOperator;
            options...) where {ELT,AFIRST}
        # Allocate scratch space
        b = zeros(src(trunc))
        blinear = zeros(ELT, length(src(trunc)))
        x1 = zeros(src(A))
        x2 = zeros(src(A))
        new(trunc, A, Zt, b, blinear, x1, x2)
    end
end

# Set type of scratch space based on operator eltype.
AZSSolver(trunc::FE_Solver{ELT}, A::DictionaryOperator{ELT}, Zt::DictionaryOperator{ELT};
        afirst=true,options...) where {ELT} =
    AZSSolver{eltype(A),Val{afirst}}(trunc, A, Zt; options...)


function AZSSolver(A::DictionaryOperator, Zt::DictionaryOperator, RD::DictionaryOperator, EF::DictionaryOperator; options...)
    TS = RestrictionSolver(A, RD, EF; options...)
    AZSSolver(TS, A, Zt; options...)
end

apply!(s::AZSSolver, coef_dest, coef_src) = _apply!(s, coef_dest, coef_src, s.A, s.Zt, s.b, s.blinear, s.TS, s.x1, s.x2)

function _apply!(s::AZSSolver, coef_dest, coef_src, A, Zt, b, blinear, TS, x1, x2)
    # Step 1:
    # Solve x2 = A^-1*b
    apply!(TS,x2,coef_src)
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

function _apply!(s::AZSSolver{ELT,Val{false}}, coef_dest, coef_src, A, Zt, b, blinear, TS, x1, x2) where {ELT}
    # Step 1:
    # Compute x2 =  Zt*b
    apply!(Zt, x2, coef_src)
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
