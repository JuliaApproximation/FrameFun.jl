# extensionframe.jl

"""
An ExtensionFrame is the restriction of a basis to a subset of its domain. This
results in a frame that implicitly represents extensions of functions on the
smaller set to the larger set.
"""
struct ExtensionFrame{S,T} <: DerivedDict{S,T}
    domain      ::  Domain
    basis       ::  Dictionary{S,T}

    function ExtensionFrame{S,T}(domain::Domain, basis::Dictionary) where {S,T}
        # @assert is_basis(basis)
        new(domain, basis)
    end
end


ExtensionFrame{S,T}(domain::Domain, basis::Dictionary{S,T}) =
    ExtensionFrame{S,T}(domain, basis)

# superdict is the function for DerivedDict's to obtain the underlying set
superdict(f::ExtensionFrame) = f.basis

basis(f::ExtensionFrame) = f.basis
domain(f::ExtensionFrame) = f.domain

"The span of the basis of the given extension frame span."


similar_dictionary(f::ExtensionFrame, dict::Dictionary) = ExtensionFrame(domain(f), dict)

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

unsafe_eval_element(s::ExtensionFrame, idx::Int, x) =
    unsafe_eval_element(basis(s), idx, x)

grid(f::ExtensionFrame) = subgrid(grid(basis(f)),domain(f))

"""
Make an ExtensionFrame, but match tensor product domains with tensor product sets
in a suitable way.

For example: an interval ⊗ a disk (= a cylinder) combined with a 3D Fourier series, leads to a
tensor product of a Fourier series on the interval ⊗ a 2D Fourier series on the disk.
"""
function extensionframe(domain::ProductDomain, basis::TensorProductDict)
    ExtensionFrames = Dictionary[]
    dc = 1
    for i = 1:nb_elements(domain)
        el = element(domain, i)
        range = dc:dc+dimension(el)-1
        push!(ExtensionFrames, ExtensionFrame(el, element(basis, range)))
        dc += dimension(el)
    end
    tensorproduct(ExtensionFrames...)
end

extensionframe(domain::Domain, basis::Dictionary) = ExtensionFrame(domain, basis)
extensionframe(basis::Dictionary, domain::Domain) = extensionframe(domain, basis)
extensionspan(span::Span, domain::Domain) = Span(extensionframe(domain, dictionary(span)))

left(d::ExtensionFrame, x...) = leftendpoint(domain(d))
right(d::ExtensionFrame, x...) = rightendpoint(domain(d))

DualGram(f::ExtensionFrame; options...) = wrap_operator(f, f, DualGram(basis(f); options...)*Gram(f; options...)*DualGram(basis(f); options...))

MixedGram(f::ExtensionFrame; options...) = wrap_operator(f, f, DualGram(basis(f); options...)*Gram(f; options...))

DiscreteDualGram(f::ExtensionFrame; oversampling=BasisFunctions.default_oversampling(f)) = wrap_operator(f, f, DiscreteDualGram(basis(f); oversampling=BasisFunctions.basis_oversampling(f, oversampling))*DiscreteGram(f; oversampling=oversampling)*DiscreteDualGram(basis(f); oversampling=BasisFunctions.basis_oversampling(f, oversampling)))

DiscreteMixedGram(f::ExtensionFrame; oversampling=BasisFunctions.default_oversampling(f)) = wrap_operator(f, f, DiscreteDualGram(basis(f); oversampling=BasisFunctions.basis_oversampling(f, oversampling))*DiscreteGram(f; oversampling=oversampling))

BasisFunctions.discrete_dual_evaluation_operator(dict::ExtensionFrame; oversampling=1, options...) =
    BasisFunctions.grid_evaluation_operator(dict, gridbasis(dict, BasisFunctions.oversampled_grid(dict, oversampling)), BasisFunctions.oversampled_grid(dict, oversampling); options...)*DiscreteDualGram(basis(dict); oversampling=BasisFunctions.basis_oversampling(dict, oversampling))


import BasisFunctions: dot
import BasisFunctions: native_nodes

native_nodes(basis::Dictionary, domain::Domains.AbstractInterval) =
    BasisFunctions.clip_and_cut(native_nodes(basis), infimum(domain), supremum(domain))

dot(span::ExtensionFrame, f1::Function, f2::Function; options...)  =
    dot(basis(span), domain(span), f1::Function, f2::Function; options...)

dot(span::ExtensionFrame, f1::Int, f2::Function; options...) =
    dot(span, x->eval_element(basis(span), f1, x), f2; options...)

dot(span::ExtensionFrame, f1::Int, f2::Int; options...) =
    dot(span, f1 ,x->eval_element(basis(span), f2, x); options...)

dot(span::Dictionary, domain::Domains.AbstractInterval, f1::Function, f2::Function; options...) =
    dot(span, f1, f2, native_nodes(span, domain); options...)
# # TODO now we assume that domainunion contains sections that do not overlap
# dot(dict::Dictionary, domain::DomainUnion, f1::Function, f2::Function; options...) =
#     dot(dict, firstelement(domain), f1, f2; options...) +
#     dot(dict, secondelement(domain), f1, f2; options...)

#continuous_approximation_operator(span::ExtensionFrame; solver = ContinuousDirectSolver, options...) = solver(span; options...)

#################
## Gram operators extended
#################
dual(span::ExtensionFrame; options...) = extensionframe(dual(basis(span); options...), domain(span))

dot(frame1::ExtensionFrame, frame2::ExtensionFrame, f1::Int, f2::Int; options...) =
    dot(frame1, frame2, x->eval_element(frame1, f1, x), x->eval_element(frame2, f2, x); options...)

function dot(frame1::ExtensionFrame, frame2::ExtensionFrame, f1::Function, f2::Function; options...)
    d1 = domain(frame1); d2 = domain(frame2)
    @assert d1 == d2
    dot(basis(frame1), basis(frame2), d1, f1, f2; options...)
end

dot(set1::Dictionary, set2::Dictionary, domain::Domains.AbstractInterval, f1::Function, f2::Function; options...) =
    dot(set1, set2, f1, f2, native_nodes(set1, set2, domain); options...)

function native_nodes(set1::Dictionary, set2::Dictionary, domain::Domains.AbstractInterval)
    @assert infimum(support(set1) )≈ infimum(support((set2)))
    @assert supremum(support((set1))) ≈ supremum(support((set2)))
    native_nodes(set1, domain)
end
