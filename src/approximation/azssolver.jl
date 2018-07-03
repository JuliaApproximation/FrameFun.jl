
# fastsolver.jl

"""

"""
struct AZSSolver{ELT} <: FE_Solver{ELT}
    TS          ::  DictionaryOperator # The low rank decomposition of (A*Zt-I)*A
    A           ::  DictionaryOperator # Store for application in step 2
    Zt          ::  DictionaryOperator # Store for application in step 2
    b                                # Scratch for right hand size
    blinear     ::  Array{ELT,1}     # Scratch for linearized right hand side (necessary for svd inproducts)
    x2
    x1

    function AZSSolver{ELT}(trunc::FE_Solver, A::DictionaryOperator, Zt::DictionaryOperator;
            options...) where ELT
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
        options...) where {ELT} =
    AZSSolver{eltype(A)}(trunc, A, Zt; options...)


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

function AZSSolver(platform::BasisFunctions.Platform, i; options...)
    frame_restriction, grid_restriction = FrameFun.spline_util_restriction_operators(platform, i)
    AZSSolver(A(platform, i), Zt(platform, i), frame_restriction', grid_restriction; options...)
end
