
"""
An ExtensionFrame is the restriction of a basis to a subset of its domain. This
results in a frame that implicitly represents extensions of functions on the
smaller set to the larger set.
"""
struct ExtensionFrame{S,T} <: DerivedDict{S,T}
    domain      ::  Domain
    basis       ::  Dictionary{S,T}

    function ExtensionFrame{S,T}(domain::Domain, basis::Dictionary) where {S,T}
        # @assert isbasis(basis)
        new(domain, basis)
    end
end

ExtensionFrame(domain::Domain, basis::Dictionary{S,T}) where {S,T} =
    ExtensionFrame{S,T}(domain, basis)

# superdict is the function for DerivedDict's to obtain the underlying set
superdict(f::ExtensionFrame) = f.basis

basis(f::ExtensionFrame) = f.basis
support(f::ExtensionFrame) = f.domain

similar_dictionary(f::ExtensionFrame, dict::Dictionary) = ExtensionFrame(support(f), dict)

isbasis(f::ExtensionFrame) = false
isframe(f::ExtensionFrame) = true
isbiorthogonal(f::ExtensionFrame) = false
isorthogonal(f::ExtensionFrame) = false
isorthonormal(f::ExtensionFrame) = false

# The following properties do not hold for extension frames
# - there is no interpolation grid
hasinterpolationgrid(f::ExtensionFrame) = false
# - there is no unitary transform
hastransform(f::ExtensionFrame) = false
hastransform(f::ExtensionFrame, dgs) = false
# - there is no antiderivative (in general)
hasantiderivative(f::ExtensionFrame) = false

name(f::ExtensionFrame) = "An extension frame of " * name(f.basis)

dict_in_support(f::ExtensionFrame, x) = x ∈ support(f)
dict_in_support(f::ExtensionFrame, idx, x) = x ∈ support(f) && in_support(basis(f), idx, x)

iscompatible(d1::ExtensionFrame, d2::ExtensionFrame) = iscompatible(basis(d1),basis(d2))

function (*)(d1::ExtensionFrame, d2::ExtensionFrame, args...)
    @assert iscompatible(d1,d2)
    (mset, mcoef) = (*)(basis(d1),basis(d2),args...)
    df = ExtensionFrame(support(d1) ∩ support(d2),mset)
    (df, mcoef)
end

unsafe_eval_element(s::ExtensionFrame, idx::Int, x) =
    unsafe_eval_element(basis(s), idx, x)

function interpolation_grid(f::ExtensionFrame)
    @warn "interpolation_grid called on extensionframe $f"
    subgrid(interpolation_grid(basis(f)),support(f))
end

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
    domain = DomainSets.ProductDomain([FrameFun.support(dicti) for dicti in elements(dict)]...)
    ExtensionFrame(domain, basis)
end

extensionframe(domain::Domain, basis::Dictionary) = ExtensionFrame(domain, basis)
extensionframe(basis::Dictionary, domain::Domain) = extensionframe(domain, basis)

leftendpoint(d::ExtensionFrame, args...) = leftendpoint(support(d))
rightendpoint(d::ExtensionFrame, args...) = rightendpoint(support(d))


import BasisFunctions: measure, restrict
import BasisFunctions: innerproduct

measure(f::ExtensionFrame) = restrict(measure(basis(f)), support(f))

innerproduct(f1::ExtensionFrame, i, f2::ExtensionFrame, j, measure; options...) =
    innerproduct(basis(f1), i, basis(f2), j, measure; options...)

grid_evaluation_operator(s::TensorProductExtensionFrameDict, dgs::GridBasis, grid::AbstractGrid; options...) =
    grid_evaluation_operator(flatten(s), dgs, grid; options...)



##################
# platform
##################

SolverStyle(dict::ExtensionFrame, ::OversamplingStyle) = AZStyle()

GridSampling(dgs::GridBasis, grid::AbstractGrid, domain::Domain, scaling) =
    GridSampling(GridBasis{coefficienttype(dgs)}(subgrid(grid, domain)), scaling=scaling)
