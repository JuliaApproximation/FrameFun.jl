# extensionframe.jl

"""
An ExtensionFrame is the restriction of a basis to a subset of its domain. This
results in a frame that implicitly represents extensions of functions on the
smaller set to the larger set.
"""
struct ExtensionFrame{T} <: DerivedSet{T}
    domain      ::  Domain
    basis       ::  FunctionSet{T}

    function ExtensionFrame{T}(domain::Domain, basis::FunctionSet) where {T}
        # @assert is_basis(basis)
        new(domain, basis)
    end
end

const ExtensionSpan{A,F <: ExtensionFrame} = Span{A,F}

"The span of the basis of the given extension frame span."
basisspan(s::ExtensionSpan) = Span(basis(s), coeftype(s))
domain(s::ExtensionSpan) = domain(set(s))
basis(s::ExtensionSpan) = basis(set(s))

ExtensionFrame{T}(domain::Domain, basis::FunctionSet{T}) =
    ExtensionFrame{T}(domain, basis)

# superset is the function for DerivedSet's to obtain the underlying set
superset(f::ExtensionFrame) = f.basis

basis(f::ExtensionFrame) = f.basis
domain(f::ExtensionFrame) = f.domain

similar_set(f::ExtensionFrame, set::FunctionSet) = ExtensionFrame(domain(f), set)

is_frame(f::ExtensionFrame) = true
is_biorthogonal(f::ExtensionFrame) = false
is_orthogonal(f::ExtensionFrame) = false
is_orthonormal(f::ExtensionFrame) = false

# The following properties do not hold for extension frames
# - there is no interpolation grid
has_grid(f::ExtensionFrame) = false
# - there is no unitary transform
has_transform(f::ExtensionFrame) = false
has_transform(f::ExtensionFrame, dgs) = false
# - there is no antiderivative (in general)
has_antiderivative(f::ExtensionFrame) = false

name(f::ExtensionFrame) = "An extension frame of " * name(f.basis)

in_support(f::ExtensionFrame, idx::Int, x) = x ∈ domain(f)

is_compatible(d1::ExtensionFrame, d2::ExtensionFrame) = is_compatible(basis(d1),basis(d2))

function (*)(d1::ExtensionFrame, d2::ExtensionFrame, args...)
    @assert is_compatible(d1,d2)
    (mset, mcoef) = (*)(basis(d1),basis(d2),args...)
    df = ExtensionFrame(domain(d1) ∩ domain(d2),mset)
    (df, mcoef)
end

# Should we check whether x lies in the domain?
call_set_expansion(e::SetExpansion, s::ExtensionFrame, coef, x...) = eval_expansion(basis(s), coef, x...)

eval_element(s::ExtensionFrame, idx::Int, x) = eval_element(basis(s), idx, x)

grid(f::ExtensionFrame) = subgrid(grid(basis(f)),domain(f))

BasisFunctions.default_oversampling(f::ExtensionFrame) = length(subgrid(BasisFunctions.oversampled_grid(basis(f), BasisFunctions.default_oversampling(basis(f))), domain(f)))/length(basis(f))


"""
Make an ExtensionFrame, but match tensor product domains with tensor product sets
in a suitable way.

For example: an interval ⊗ a disk (= a cylinder) combined with a 3D Fourier series, leads to a
tensor product of a Fourier series on the interval ⊗ a 2D Fourier series on the disk.
"""
function extensionframe(domain::ProductDomain, basis::TensorProductSet)
    ExtensionFrames = FunctionSet[]
    dc = 1
    for i = 1:nb_elements(domain)
        el = element(domain, i)
        range = dc:dc+dimension(el)-1
        push!(ExtensionFrames, ExtensionFrame(el, element(basis, range)))
        dc += dimension(el)
    end
    tensorproduct(ExtensionFrames...)
end

extensionframe(domain::Domain, basis::FunctionSet) = ExtensionFrame(domain, basis)
extensionframe(basis::FunctionSet, domain::Domain) = extensionframe(domain, basis)
extensionspan(span::Span, domain::Domain) = Span(extensionframe(domain, set(span)))

left(d::ExtensionFrame, x...) = leftendpoint(domain(d))
right(d::ExtensionFrame, x...) = rightendpoint(domain(d))

DualGram(f::ExtensionSpan; options...) = wrap_operator(f, f, DualGram(basisspan(f); options...)*Gram(f; options...)*DualGram(basisspan(f); options...))

MixedGram(f::ExtensionSpan; options...) = wrap_operator(f, f, DualGram(basisspan(f); options...)*Gram(f; options...))

DiscreteDualGram(f::ExtensionSpan; oversampling=BasisFunctions.default_oversampling(set(f))) = wrap_operator(f, f, DiscreteDualGram(basisspan(f); oversampling=BasisFunctions.basis_oversampling(set(f), oversampling))*DiscreteGram(f; oversampling=oversampling)*DiscreteDualGram(basisspan(f); oversampling=BasisFunctions.basis_oversampling(set(f), oversampling)))

DiscreteMixedGram(f::ExtensionSpan; oversampling=BasisFunctions.default_oversampling(set(f))) = wrap_operator(f, f, DiscreteDualGram(basisspan(f); oversampling=BasisFunctions.basis_oversampling(set(f), oversampling))*DiscreteGram(f; oversampling=oversampling))

BasisFunctions.discrete_dual_evaluation_operator(span::ExtensionSpan; oversampling=1, options...) =
    BasisFunctions.grid_evaluation_operator(span, gridspace(span, BasisFunctions.oversampled_grid(set(span), oversampling)), BasisFunctions.oversampled_grid(set(span), oversampling); options...)*DiscreteDualGram(basisspan(span); oversampling=BasisFunctions.basis_oversampling(set(span), oversampling))


import BasisFunctions: dot
import BasisFunctions: native_nodes

native_nodes(basis::FunctionSet, domain::Interval) =
  BasisFunctions.clip_and_cut(native_nodes(basis), leftendpoint(domain), rightendpoint(domain))

dot(span::ExtensionSpan, f1::Function, f2::Function; options...)  =
    dot(basisspan(span), domain(span), f1::Function, f2::Function; options...)

dot(span::ExtensionSpan, f1::Int, f2::Function; options...) =
    dot(span, x->eval_element(basis(span), f1, x), f2; options...)

dot(span::ExtensionSpan, f1::Int, f2::Int; options...) =
    dot(span, f1 ,x->eval_element(basis(span), f2, x); options...)

dot(span::Span, domain::Interval, f1::Function, f2::Function; options...) =
    dot(span, f1, f2, native_nodes(set(span), domain); options...)
# # TODO now we assume that domainunion contains sections that do not overlap
# dot(set::FunctionSet, domain::DomainUnion, f1::Function, f2::Function; options...) =
#     dot(set, firstelement(domain), f1, f2; options...) +
#     dot(set, secondelement(domain), f1, f2; options...)
# 
# continuous_approximation_operator(span::ExtensionSpan; solver = ContinuousDirectSolver, options...) = solver(span; options...)

#################
## Gram operators extended
#################
dual(span::ExtensionSpan; options...) = extensionspan(dual(basisspan(span); options...), domain(span))

dot(frame1::ExtensionSpan, frame2::ExtensionSpan, f1::Int, f2::Int; options...) =
    dot(frame1, frame2, x->eval_element(set(frame1), f1, x), x->eval_element(set(frame2), f2, x); options...)

function dot(frame1::ExtensionFrame, frame2::ExtensionFrame, f1::Function, f2::Function; options...)
    d1 = domain(frame1); d2 = domain(frame2)
    @assert d1 == d2
    dot(basisspan(frame1), basisspan(frame2), d1, f1, f2; options...)
end

dot(set1::ExtensionSpan, set2::ExtensionSpan, domain::Interval, f1::Function, f2::Function; options...) =
    dot(set1, set2, f1, f2, native_nodes(set(set1), set(set2), domain); options...)

function native_nodes(set1::FunctionSet, set2::FunctionSet, domain::Interval)
    @assert left(set1) ≈ left(set2)
    @assert left(set2) ≈ left(set2)
    native_nodes(set1, domain)
end
