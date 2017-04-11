# SumFrame.jl

"""
A SumFrame is the concatenation of multiple sets, where the first set is a basis
with an associated unitary transform (such as Fourier or Chebyshev).
"""
immutable SumFrame{N,T} <: DerivedSet{N,T}
    superset    ::  FunctionSet{N,T}

    function SumFrame(set)
        @assert is_sumset(set)
        new(set)
    end
end

SumFrame{N,T}(set::FunctionSet{N,T}) = SumFrame{N,T}(set)

is_sumset(set::FunctionSet) = false

is_sumset(set::MultiSet) = has_transform(element(set,1))

similar_set(f::SumFrame, set::FunctionSet) = SumFrame(set)

sumframe(sets::FunctionSet...) = SumFrame(multiset(sets...))

name(f::SumFrame) = "A sum frame"



default_threshold{T <: Number}(::Type{T}) = 1000eps(T)

default_threshold{T <: Number}(::Type{Complex{T}}) = 1000eps(T)


"Solver to approximate functions using a SumFrame."
immutable SumFrameSolver{T} <: AbstractOperator{T}
    src     ::  FunctionSet
    dest    ::  FunctionSet
    set1    ::  FunctionSet
    set2    ::  FunctionSet

    A1
    P
    A1prime
    normalization1
    lr
    L
end


function SumFrameSolver(sumset)
    ELT = eltype(sumset)
    set1 = element(sumset, 1)
    set2 = tail(sumset)
    extended_set = resize(set1, 2*length(sumset))
    dgs = DiscreteGridSpace(extended_set)
    T = transform_operator(extended_set, dgs)
    Npre = transform_operator_pre(extended_set, dgs)
    Npost = transform_operator_post(extended_set, dgs)
    E1 = extension_operator(set1, extended_set)
    R1 = restriction_operator(extended_set, set1)
    A1 = T*E1
    P = IdentityOperator(dgs) - A1*A1'
    L = evaluation_operator(set2, dgs)
    lr = TruncatedSvdSolver(P*L, cutoff = default_threshold(ELT), growth_factor = 2)
    normalization1 = inv(DiagonalOperator(set1, set1, diagonal(R1*Npre*E1)))
    SumFrameSolver{ELT}(dgs, sumset, set1, set2, A1, P, A1', normalization1, lr, L)
end


function apply!(op::SumFrameSolver, coef_dest, coef_src)
    B = coef_src
    PB = op.P * B
    x2 = op.lr * PB
    B2 = B - op.L * x2
    x1 = op.normalization1 * (op.A1prime * B2)
    coef_dest[:] = [x1[:]; x2[:]]
end

approximation_operator(set::SumFrame; options...) = SumFrameSolver(set)

# The SumFrameSolver is also the default for a MultiSet
function approximation_operator(set::MultiSet; options...)
    if is_sumset(set)
        SumFrameSolver(set)
    else
        default_approximation_operator(set; options...)
    end
end
