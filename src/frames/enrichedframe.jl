# enrichedframe.jl

"""
An EnrichedFrame is the concatenation of a basis with a bounded number of extra
functions.
"""
immutable EnrichedFrame{N,T} <: DerivedSet{N,T}
    superset    ::  FunctionSet{N,T}

    function EnrichedFrame(set)
        @assert is_enrichedframe(set)
        new(set)
    end
end

EnrichedFrame{N,T}(set::FunctionSet{N,T}) = EnrichedFrame{N,T}(set)

is_enrichedframe(set::FunctionSet) = false

is_enrichedframe(set::MultiSet) = is_enrichedframe(set, element(set, 1), element(set, 2:composite_length(set)))
is_enrichedframe(set::MultiSet, set1::FunctionSet, sets::FunctionSet) =
    has_transform(set1)


similar_set(f::EnrichedFrame, set::FunctionSet) = EnrichedFrame(set)

# The basis of the enriched frame is the first element
basis(f::EnrichedFrame) = element(superset(f),1)

name(f::EnrichedFrame) = "The basis " * name(element(f,1)) * " enriched with " * name(element(f,2))


############
# Extension
############

# We only extend the basis, not the other set(s)
extension_size(f::EnrichedFrame) = extension_size(basis(f))

resize(f::EnrichedFrame, n) =
    similar_set(f, MultiSet([resize(basis(f), n); elements(f)[2:end]]))


############
# Solver
############

"""
Solver to approximate in an enriched frame.
"""
immutable EnrichedFrameSolver{T} <: AbstractOperator{T}
    src     ::  FunctionSet
    dest    ::  FunctionSet

    Tpre
    Trans
    Tpost
    E
    R
    A
    Aprime
    B
    P
    C
    lr
end


function EnrichedFrameSolver(frame)
    set1 = element(frame, 1)
    extended_set = resize(set1, 2*length(frame))
    extended_grid = grid(extended_set)
    dgs = DiscreteGridSpace(extended_set)
    Tpre, Trans, Tpost = transform_operators(extended_set, dgs)
    E = extension_operator(set1, extended_set)
    R = restriction_operator(extended_set, set1)
    A = Trans*E
    Aprime = A'
    B = evaluation_operator(frame, extended_grid)
    P = IdentityOperator(dgs) - A*A'
    C = P*B
    lr = lowrank_approximation(C)

    T = eltype(frame)
    EnrichedFrameSolver{T}(dgs, frame, Tpre, Trans, Tpost, E, R, A, Aprime, B, P, C, lr)
end


function apply!(op::EnrichedFrameSolver, coef_dest, coef_src)
    B = coef_src
    B1 = inv(op.Tpost) * B
    PB = op.P * B1
    x2 = lowranksolve(op.lr, PB)
    B2 = B1 - inv(op.Tpost) * (op.B * x2)
    x1 = op.Aprime * B2
    x1 = (op.R*inv(op.Tpre)*op.E) * x1
    x2[1:length(src(op.E))] += x1
    coef_dest[:] = x2
end

approximation_operator(frame::EnrichedFrame; options...) = EnrichedFrameSolver(frame)
