# weightedsumframe.jl

"""
A WeightedSumFrame is the concatenation of multiple sets, possibly weighted
except for the first one.
"""
immutable WeightedSumFrame{N,T} <: DerivedSet{N,T}
    superset    ::  FunctionSet{N,T}

    function WeightedSumFrame(set)
        @assert is_weightedsumset(set)
        new(set)
    end
end

WeightedSumFrame{N,T}(set::FunctionSet{N,T}) = WeightedSumFrame{N,T}(set)

is_weightedsumset(set::FunctionSet) = false

is_weightedsumset(set::MultiSet) = _is_weightedsumset(set, elements(set)...)

_is_weightedsumset(set::MultiSet, sets::FunctionSet...) = false
_is_weightedsumset(set::MultiSet, sets::WeightedSet...) = false
_is_weightedsumset(set::MultiSet, set1::FunctionSet, set2::WeightedSet) = true

similar_set(f::WeightedSumFrame, set::FunctionSet) = WeightedSumFrame(set)

sumframe(sets::FunctionSet...) = WeightedSumFrame(multiset(sets...))

name(f::WeightedSumFrame) = "A weighted sum frame"

"""
Solver to approximate functions of the form f(x) = p(x) + w(x) q(x), where w is a given
weight function, and p and q are well approximated by set1 and set2 respectively.
"""
immutable ConcatSolver{T} <: AbstractOperator{T}
    src     ::  FunctionSet
    dest    ::  FunctionSet
    set1    ::  FunctionSet
    set2    ::  FunctionSet
    fun

    A1
    A2
    D
    P
    C
    A1prime
    A2prime
    normalization1
    normalization2
    lr
end


function ConcatSolver(concatset, set1, set2, weightfunction)
    ELT = eltype(concatset)
    extended_set = resize(set1, 2*(length(set1)+length(set2)))
    ls_grid = grid(extended_set)
    dgs = DiscreteGridSpace(extended_set)
    T = transform_operator(extended_set, dgs)
    N = transform_operator_post(dgs, extended_set)
    E1 = extension_operator(set1, extended_set)
    E2 = extension_operator(set2, extended_set)
    R1 = restriction_operator(extended_set, set1)
    R2 = restriction_operator(extended_set, set2)
    D = DiagonalOperator(dgs, ELT[weightfunction(x) for x in ls_grid])
    A1 = T*E1
    A2 = T*E2
    P = IdentityOperator(dgs) - A1*A1'
    C = P*D*A2
    lr = lowrank_approximation(C)
    normalization1 = matrix(R1*N*E1)
    normalization2 = matrix(R2*N*E2)
    ConcatSolver{ELT}(dgs, concatset, set1, set2, weightfunction,
        A1, A2, D, P, C, A1', A2', normalization1, normalization2, lr)
end


function apply!(op::ConcatSolver, coef_dest, coef_src)
    B = coef_src
    PB = op.P * B
    x2 = lowranksolve(op.lr, PB)
    B2 = B - op.D * (op.A2 * x2)
    x1 = op.A1prime * B2
    x1b = op.normalization1*x1
    x2b = op.normalization2*x2
    coef_dest[:] = [x1b; x2b]
end

function approximation_operator(set::WeightedSumFrame; options...)
    approximation_operator_concat(set, element(set, 1), element(set, 2); options...)
end

approximation_operator_concat(set::WeightedSumFrame, set1::FunctionSet, set2::FunctionSet) =
    default_approximation_operator(set)

function approximation_operator_concat(set::WeightedSumFrame, set1::FunctionSet, set2::WeightedSet)
    if has_transform(set1) && has_transform(superset(set2))
        ConcatSolver(set, set1, superset(set2), weightfunction(set2))
    else
        default_approximation_operator(set)
    end
end

default_threshold{T <: Number}(::Type{T}) = 1000eps(T)

default_threshold{T <: Number}(::Type{Complex{T}}) = 1000eps(T)

function lowrank_approximation(A; ɛ = default_threshold(eltype(A)), estimated_rank = 2)
    T = eltype(A)
    (m,n) = size(A)
    r = min(estimated_rank,n)
    random_matrix = map(T, rand(n, r))
    C = apply_multiple(A, random_matrix)
    c = cond(C)
    m = maximum(abs(C))
    while (c < m/ɛ) && (r < n)
        r0 = r
        r = min(2r,n)
        extra_random_matrix = map(T, rand(n,r-r0))
        extraC = apply_multiple(A, extra_random_matrix)
        C = [C extraC]
        random_matrix = [random_matrix extra_random_matrix]
        c = cond(C)
        m = maximum(abs(C))
    end
    USV = LAPACK.gesvd!('S','S',C)
    S = USV[2]
    rank = findlast(S.>ɛ)
    u = USV[1][:,1:rank]
    s = S[1:rank]
    sinv = s.^(-1)
    v = USV[3][1:rank,:]
    ut = u'
    vs = v' * diagm(sinv)
    ut,vs,random_matrix
end

function lowranksolve(lr, b)
    ut,vs,m = lr
    m * (vs * (ut * b))
end
