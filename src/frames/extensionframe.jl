
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

ExtensionFrame(domain::Domain, basis::Dictionary{S,T}) where {S,T} =
    ExtensionFrame{S,T}(domain, basis)

# superdict is the function for DerivedDict's to obtain the underlying set
superdict(f::ExtensionFrame) = f.basis

basis(f::ExtensionFrame) = f.basis
domain(f::ExtensionFrame) = f.domain

similar_dictionary(f::ExtensionFrame, dict::Dictionary) = ExtensionFrame(domain(f), dict)

is_basis(f::ExtensionFrame) = false
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

dict_in_support(f::ExtensionFrame, x) = x ∈ domain(f)
dict_in_support(f::ExtensionFrame, idx, x) = x ∈ domain(f) && in_support(basis(f), idx, x)

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
    for i = 1:numelements(domain)
        el = element(domain, i)
        range = dc:dc+dimension(el)-1
        push!(ExtensionFrames, ExtensionFrame(el, element(basis, range)))
        dc += dimension(el)
    end
    tensorproduct(ExtensionFrames...)
end

const TensorProductExtensionFrameDict{N,N1,S,T} = TensorProductDict{N,NTuple{N1,DT},S,T} where {N,N1,DT<:ExtensionFrame,S,T}

# TODO remove need of this function
function flatten(dict::TensorProductExtensionFrameDict)
    basis = tensorproduct([superdict(dicti) for dicti in elements(dict)]...)
    domain = DomainSets.ProductDomain([FrameFun.domain(dicti) for dicti in elements(dict)]...)
    ExtensionFrame(domain, basis)
end

extensionframe(domain::Domain, basis::Dictionary) = ExtensionFrame(domain, basis)
extensionframe(basis::Dictionary, domain::Domain) = extensionframe(domain, basis)

left(d::ExtensionFrame, x...) = leftendpoint(domain(d))
right(d::ExtensionFrame, x...) = rightendpoint(domain(d))

DualGram(f::ExtensionFrame; options...) = wrap_operator(f, f, DualGram(basis(f); options...)*Gram(f; options...)*DualGram(basis(f); options...))

MixedGram(f::ExtensionFrame; options...) = wrap_operator(f, f, DualGram(basis(f); options...)*Gram(f; options...))

DiscreteDualGram(f::ExtensionFrame; oversampling=BasisFunctions.default_oversampling(f)) = wrap_operator(f, f, DiscreteDualGram(basis(f); oversampling=BasisFunctions.basis_oversampling(f, oversampling))*DiscreteGram(f; oversampling=oversampling)*DiscreteDualGram(basis(f); oversampling=BasisFunctions.basis_oversampling(f, oversampling)))

DiscreteMixedGram(f::ExtensionFrame; oversampling=BasisFunctions.default_oversampling(f)) = wrap_operator(f, f, DiscreteDualGram(basis(f); oversampling=BasisFunctions.basis_oversampling(f, oversampling))*DiscreteGram(f; oversampling=oversampling))

BasisFunctions.discrete_dual_evaluation_operator(dict::ExtensionFrame; oversampling=1, options...) =
    BasisFunctions.grid_evaluation_operator(dict, gridbasis(dict, BasisFunctions.oversampled_grid(dict, oversampling)), BasisFunctions.oversampled_grid(dict, oversampling); options...)*DiscreteDualGram(basis(dict); oversampling=BasisFunctions.basis_oversampling(dict, oversampling))

import BasisFunctions: measure, restrict
import BasisFunctions: innerproduct

measure(f::ExtensionFrame) = restrict(measure(basis(f)), domain(f))
innerproduct(f1::ExtensionFrame, i, f2::ExtensionFrame, j, measure) =
    innerproduct(basis(f1), i, basis(f2), j, measure)

import BasisFunctions: dot
import BasisFunctions: native_nodes

native_nodes(basis::Dictionary, domain::DomainSets.AbstractInterval) =
    BasisFunctions.clip_and_cut(native_nodes(basis), infimum(domain), supremum(domain))

dot(dict::ExtensionFrame, f1::Function, f2::Function; options...)  =
    dot(basis(dict), domain(dict), f1::Function, f2::Function; options...)

dot(dict::ExtensionFrame, f1::Int, f2::Function; options...) =
    dot(dict, x->eval_element(basis(dict), f1, x), f2; options...)

dot(dict::ExtensionFrame, f1::Int, f2::Int; options...) =
    dot(dict, f1 ,x->eval_element(basis(dict), f2, x); options...)

dot(dict::Dictionary, domain::DomainSets.AbstractInterval, f1::Function, f2::Function; options...) =
    dot(dict, f1, f2, native_nodes(dict, domain); options...)
# # TODO now we assume that domainunion contains sections that do not overlap
# dot(dict::Dictionary, domain::DomainUnion, f1::Function, f2::Function; options...) =
#     dot(dict, firstelement(domain), f1, f2; options...) +
#     dot(dict, secondelement(domain), f1, f2; options...)

#continuous_approximation_operator(dict::ExtensionFrame; solver = ContinuousDirectSolver, options...) = solver(dict; options...)

#################
## Gram operators extended
#################
dual(dict::ExtensionFrame; options...) = extensionframe(dual(basis(dict); options...), domain(dict))

dot(frame1::ExtensionFrame, frame2::ExtensionFrame, f1::Int, f2::Int; options...) =
    dot(frame1, frame2, x->eval_element(frame1, f1, x), x->eval_element(frame2, f2, x); options...)

function dot(frame1::ExtensionFrame, frame2::ExtensionFrame, f1::Function, f2::Function; options...)
    d1 = domain(frame1); d2 = domain(frame2)
    @assert d1 == d2
    dot(basis(frame1), basis(frame2), d1, f1, f2; options...)
end

dot(set1::Dictionary, set2::Dictionary, domain::DomainSets.AbstractInterval, f1::Function, f2::Function; options...) =
    dot(set1, set2, f1, f2, native_nodes(set1, set2, domain); options...)

function native_nodes(set1::Dictionary, set2::Dictionary, domain::DomainSets.AbstractInterval)
    @assert infimum(support(set1) )≈ infimum(support((set2)))
    @assert supremum(support((set1))) ≈ supremum(support((set2)))
    native_nodes(set1, domain)
end

grid_evaluation_operator(s::TensorProductExtensionFrameDict, dgs::GridBasis, grid::AbstractGrid; options...) =
    grid_evaluation_operator(flatten(s), dgs, grid; options...)
grid_evaluation_operator(s::TensorProductExtensionFrameDict, dgs::GridBasis, grid::AbstractSubGrid; options...) =
    grid_evaluation_operator(flatten(s), dgs, grid; options...)



##################
# platform
##################

sampler(platform::Platform, sampler::GridSamplingOperator, domain::Domain) =
    GridSamplingOperator(gridbasis(sampler), grid(sampler), domain, sampler.scaling)
GridSamplingOperator(dgs::GridBasis, grid::AbstractGrid, domain::Domain, scaling) =
    GridSamplingOperator(gridbasis(FrameFun.subgrid(grid, domain), coeftype(dgs)), scaling=scaling)

extension_frame_sampler(platform::Platform, domain::Domain) = n->sampler(platform, platform.sampler_generator(n), domain)
dual_extension_frame_sampler(platform::Platform, domain::Domain) = n->sampler(platform, platform.dual_sampler_generator(n), domain)

function extension_frame_platform(platform::GenericPlatform, domain::Domain)
    primal = n->extensionframe(platform.primal_generator(n), domain)
    dual = n->platform.dual_generator(n)
    sampler = extension_frame_sampler(platform, domain)
    dual_sampler = dual_extension_frame_sampler(platform, domain)
    GenericPlatform(super_platform=platform, primal = primal, dual = dual, sampler = sampler,dual_sampler = dual_sampler,
        params = platform.parameter_sequence, name = "extension frame of " * platform.name)
end

struct ExtensionFramePlatform <: FramePlatform
    basisplatform   ::  Platform
    domain          ::  Domain
end

dictionary(p::ExtensionFramePlatform, param) = extensionframe(dictionary(p.basisplatform, param), p.domain)
