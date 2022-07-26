using BasisFunctions: DictionarySolverOperator
import BasisFunctions: operator, apply!

export AZSolver
"""
Fast implementation of the AZ Algorithm. Returns an operator that is approximately
`pinv(A)`, based on a suitable `Zt` so that `(A*Zt-I)*A` is low rank.

Steps when applying to right hand side b:

1. (I-A*Zt)*A*x2=(I-A*Zt)*b, solve for x2
This is by default done using a low rank decomposition of (I-A*Zt)*A (RandomizedSvdSolver).

2. x1 = Zt*(b-A*x2)
Can be done fast if Zt and A are fast.

3. x = x1+x2
This is the solution.
"""
struct AZSolver{T} <: DictionarySolverOperator{T}
    A           ::  DictionaryOperator # Store for application in step 2
    Zt          ::  DictionaryOperator # Store for application in step 2
    plunge_op   ::  DictionaryOperator # (A*Zt-I), store because it allocates memory
    psolver     ::  DictionaryOperator # The low rank plunge solver for (A*Zt-I)*A
    Pb                                 # Scratch space for the plunge opreator applied to the right hand size
    x1                                 # Scratch space for the intermediate result vector x1

    function AZSolver{T}(A::DictionaryOperator, Zt::DictionaryOperator, plunge_op::DictionaryOperator,
            psolver::DictionaryOperator) where {T}
        # Allocate scratch space
        Pb = zeros(T, src(psolver))
        x1 = zeros(T, src(A))
        new(A, Zt, plunge_op, psolver, Pb, x1)
    end
end

AZSolver(A::DictionaryOperator{T}, Zt::DictionaryOperator{T}, plunge_op::DictionaryOperator{T}, psolver::DictionaryOperator{T}) where {T} =
    AZSolver{eltype(A)}(A, Zt, plunge_op, psolver)

operator(s::AZSolver) = s.A

default_regularization(A; options...) = pSVD_solver(A; options...)

function default_regularization(A::DictionaryOperator{BigFloat}; options...)
    @warn "Our own implementation of a randomized SVD solver is used with `BigFloat`s.
        This implementation is not as robust as the partial SVD solver in `LowRankApprox.jl`.
        Concider using `DoubleFloats.jl` instead."
    rSVD_solver(A; options...)
end
function default_regularization(A::DictionaryOperator{Complex{BigFloat}}; options...)
    @warn "Our own implementation of a randomized SVD solver is used with `BigFloat`s.
        This implementation is not as robust as the partial SVD solver in `LowRankApprox.jl`.
        Concider using `DoubleFloats.jl` instead."
    rSVD_solver(A; options...)
end

function AZSolver(A::DictionaryOperator, Zt::DictionaryOperator;
            REG = default_regularization,
            rankestimate = 40,
            threshold = default_threshold(A),
            verbose = false,
            weightedAZ = false,
            options...)

    plunge_op = I-A*Zt
    if weightedAZ
        verbose && @info "Weighted AZ"
        AZ_Cweight = haskey(options,:AZ_Cweight) ? options[:AZ_Cweight] : error("No options `AZ_Cweight`")
        psolver = AZ_Cweight*REG(plunge_op*A*AZ_Cweight; verbose, threshold, rankestimate, options...)
    else
        psolver = REG(plunge_op*A; verbose, threshold, rankestimate, options...)
    end
    AZSolver(A, Zt, plunge_op, psolver)
end

default_threshold(A::DictionaryOperator) = regularization_threshold(eltype(A))

apply!(s::AZSolver, coef_dest, coef_src) = _apply!(s, coef_dest, coef_src,
        s.plunge_op, s.A, s.Zt, s.psolver, s.Pb, s.x1)

function _apply!(s::AZSolver, coef_dest, coef_src, plunge_op::DictionaryOperator, A, Zt, psolver, Pb, x1)
    # Step 1: Solve (I-A*Zt)*A * x1 = (I-A*Zt)*b
    apply!(plunge_op, Pb, coef_src)
    apply!(psolver, x1, Pb)

    # Step 2: Compute x2 =  Zt*(b-A*x1)
    apply!(A, Pb, x1)
    Pb .= coef_src .- Pb
    apply!(Zt, coef_dest, Pb)

    # Step 3: x = x1 + x2
    coef_dest .+= x1
end
