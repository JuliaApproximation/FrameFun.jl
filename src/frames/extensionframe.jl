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

const ExtensionSpan{A,S,T,D <: ExtensionFrame} = Span{A,S,T,D}

ExtensionFrame{S,T}(domain::Domain, basis::Dictionary{S,T}) =
    ExtensionFrame{S,T}(domain, basis)

# superdict is the function for DerivedDict's to obtain the underlying set
superdict(f::ExtensionFrame) = f.basis
superspan(s::ExtensionSpan) = Span(superdict(dictionary(s)), coeftype(s))

basis(f::ExtensionFrame) = f.basis
domain(f::ExtensionFrame) = f.domain

"The span of the basis of the given extension frame span."
basisspan(s::ExtensionSpan) = Span(basis(s), coeftype(s))
domain(s::ExtensionSpan) = domain(dictionary(s))
basis(s::ExtensionSpan) = basis(dictionary(s))


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

const TensorProductExtensionFrameDict{N,N1,S,T} = TensorProductDict{N,NTuple{N1,DT},S,T} where {N,N1,DT<:ExtensionFrame,S,T}
const TensorProductExtensionFrameSpan{A,S,T,D <: TensorProductExtensionFrameDict} = Span{A,S,T,D}
# TODO remove need of this function
function flatten(dict::TensorProductExtensionFrameDict)
    basis = tensorproduct([superdict(dicti) for dicti in elements(dict)]...)
    domain = Domains.ProductDomain([FrameFun.domain(dicti) for dicti in elements(dict)]...)
    ExtensionFrame(domain, basis)
end
flatten(s::TensorProductExtensionFrameSpan) = Span(flatten(dictionary(s)))

extensionframe(domain::Domain, basis::Dictionary) = ExtensionFrame(domain, basis)
extensionframe(basis::Dictionary, domain::Domain) = extensionframe(domain, basis)
extensionspan(span::Span, domain::Domain) = Span(extensionframe(domain, dictionary(span)))

left(d::ExtensionFrame, x...) = leftendpoint(domain(d))
right(d::ExtensionFrame, x...) = rightendpoint(domain(d))

DualGram(f::ExtensionSpan; options...) = wrap_operator(f, f, DualGram(basisspan(f); options...)*Gram(f; options...)*DualGram(basisspan(f); options...))

MixedGram(f::ExtensionSpan; options...) = wrap_operator(f, f, DualGram(basisspan(f); options...)*Gram(f; options...))

DiscreteDualGram(f::ExtensionSpan; oversampling=BasisFunctions.default_oversampling(dictionary(f))) = wrap_operator(f, f, DiscreteDualGram(basisspan(f); oversampling=BasisFunctions.basis_oversampling(dictionary(f), oversampling))*DiscreteGram(f; oversampling=oversampling)*DiscreteDualGram(basisspan(f); oversampling=BasisFunctions.basis_oversampling(dictionary(f), oversampling)))

DiscreteMixedGram(f::ExtensionSpan; oversampling=BasisFunctions.default_oversampling(dictionary(f))) = wrap_operator(f, f, DiscreteDualGram(basisspan(f); oversampling=BasisFunctions.basis_oversampling(dictionary(f), oversampling))*DiscreteGram(f; oversampling=oversampling))

BasisFunctions.discrete_dual_evaluation_operator(span::ExtensionSpan; oversampling=1, options...) =
    BasisFunctions.grid_evaluation_operator(span, gridspace(span, BasisFunctions.oversampled_grid(dictionary(span), oversampling)), BasisFunctions.oversampled_grid(dictionary(span), oversampling); options...)*DiscreteDualGram(basisspan(span); oversampling=BasisFunctions.basis_oversampling(dictionary(span), oversampling))


import BasisFunctions: dot
import BasisFunctions: native_nodes

native_nodes(basis::Dictionary, domain::Interval) =
  BasisFunctions.clip_and_cut(native_nodes(basis), leftendpoint(domain), rightendpoint(domain))

dot(span::ExtensionSpan, f1::Function, f2::Function; options...)  =
    dot(basisspan(span), domain(span), f1::Function, f2::Function; options...)

dot(span::ExtensionSpan, f1::Int, f2::Function; options...) =
    dot(span, x->eval_element(basis(span), f1, x), f2; options...)

dot(span::ExtensionSpan, f1::Int, f2::Int; options...) =
    dot(span, f1 ,x->eval_element(basis(span), f2, x); options...)

dot(span::Span, domain::Interval, f1::Function, f2::Function; options...) =
    dot(span, f1, f2, native_nodes(dictionary(span), domain); options...)
# # TODO now we assume that domainunion contains sections that do not overlap
# dot(dict::Dictionary, domain::DomainUnion, f1::Function, f2::Function; options...) =
#     dot(dict, firstelement(domain), f1, f2; options...) +
#     dot(dict, secondelement(domain), f1, f2; options...)

#continuous_approximation_operator(span::ExtensionSpan; solver = ContinuousDirectSolver, options...) = solver(span; options...)

#################
## Gram operators extended
#################
dual(span::ExtensionSpan; options...) = extensionspan(dual(basisspan(span); options...), domain(span))

dot(frame1::ExtensionSpan, frame2::ExtensionSpan, f1::Int, f2::Int; options...) =
    dot(frame1, frame2, x->eval_element(dictionary(frame1), f1, x), x->eval_element(dictionary(frame2), f2, x); options...)

function dot(frame1::ExtensionFrame, frame2::ExtensionFrame, f1::Function, f2::Function; options...)
    d1 = domain(frame1); d2 = domain(frame2)
    @assert d1 == d2
    dot(basisspan(frame1), basisspan(frame2), d1, f1, f2; options...)
end

dot(set1::ExtensionSpan, set2::ExtensionSpan, domain::Interval, f1::Function, f2::Function; options...) =
    dot(set1, set2, f1, f2, native_nodes(dictionary(set1), dictionary(set2), domain); options...)

function native_nodes(set1::Dictionary, set2::Dictionary, domain::Interval)
    @assert left(set1) ≈ left(set2)
    @assert left(set2) ≈ left(set2)
    native_nodes(set1, domain)
end


grid_evaluation_operator(s::TensorProductExtensionFrameSpan, dgs::DiscreteGridSpace, grid::AbstractGrid; options...) =
    grid_evaluation_operator(flatten(s), dgs, grid; options...)
grid_evaluation_operator(s::TensorProductExtensionFrameSpan, dgs::DiscreteGridSpace, grid::AbstractSubGrid; options...) =
    grid_evaluation_operator(flatten(s), dgs, grid; options...)

function BasisFunctions.DWTSamplingOperator(span::FrameFun.ExtensionSpan, oversampling::Int=1, recursion::Int=0)
    S = BasisFunctions.GridSamplingOperator(gridspace(BasisFunctions.dwt_oversampled_grid(dictionary(span), oversampling, recursion), coeftype(span)))
    new_oversampling = Int(length(supergrid(grid(S)))/length(span))>>recursion
    E = extension_operator(gridspace(S), gridspace(supergrid(grid(S)), coeftype(span)))
    W = BasisFunctions.WeightOperator(FrameFun.basisspan(span), new_oversampling, recursion)
    BasisFunctions.DWTSamplingOperator(S, W*E)
end

##################
# platform
##################

BasisFunctions.sampler(platform::Platform, sampler::GridSamplingOperator, domain::Domain) =
    GridSamplingOperator(gridspace(sampler), grid(sampler), domain)
BasisFunctions.GridSamplingOperator(dgs::DiscreteGridSpace, grid::AbstractGrid, domain::Domain) =
    GridSamplingOperator(gridspace(FrameFun.subgrid(grid, domain), coeftype(dgs)))

function BasisFunctions.sampler(platform::BasisFunctions.GenericPlatform, sampler::BasisFunctions.DWTSamplingOperator, domain::Domain)
    S = BasisFunctions.sampler(platform, sampler.sampler, domain)
    E = extension_operator(gridspace(S), gridspace(supergrid(grid(S)), coeftype(primal(platform, 1))))
    BasisFunctions.DWTSamplingOperator(S,sampler.weight*E)
end

extension_frame_sampler(platform::Platform, domain::Domain) = n->sampler(platform, platform.sampler_generator(n), domain)
dual_extension_frame_sampler(platform::Platform, domain::Domain) = n->sampler(platform, platform.dual_sampler_generator(n), domain)

function extension_frame_platform(platform::BasisFunctions.GenericPlatform, domain::Domain)
    primal = n->extensionframe(platform.primal_generator(n), domain)
    dual = n->platform.dual_generator(n)
    sampler = extension_frame_sampler(platform, domain)
    dual_sampler = dual_extension_frame_sampler(platform, domain)
    BasisFunctions.GenericPlatform(super_platform=platform, primal = primal, dual = dual, sampler = sampler,dual_sampler = dual_sampler,
        params = platform.parameter_sequence, name = "extension frame of " * platform.name)
end
