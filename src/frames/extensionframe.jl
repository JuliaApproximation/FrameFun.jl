
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

basis(f::ExtensionFrame) = superdict(f)
support(f::ExtensionFrame) = f.domain

similar_dictionary(f::ExtensionFrame, dict::Dictionary) = ExtensionFrame(support(f), dict)

isbasis(f::ExtensionFrame) = false
isframe(f::ExtensionFrame) = true
isbiorthogonal(f::ExtensionFrame, measure::Measure) = false
isorthogonal(f::ExtensionFrame, measure::Measure) = false
isorthonormal(f::ExtensionFrame, measure::Measure) = false

# The following properties do not hold for extension frames
# - there is no interpolation grid
hasinterpolationgrid(f::ExtensionFrame) = false
# - there is no unitary transform
hastransform(f::ExtensionFrame) = false
hastransform(f::ExtensionFrame, dgs) = false
# - there is no antiderivative (in general)
hasantiderivative(f::ExtensionFrame) = false

name(f::ExtensionFrame) = "Extension frame"


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

const ExtensionFrameTensor = TensorProductDict{D,NTuple{D,DICT}} where D where {DICT<:ExtensionFrame}

extensionframe(domain::Domain, basis::Dictionary) = ExtensionFrame(domain, basis)
extensionframe(basis::Dictionary, domain::Domain) = extensionframe(domain, basis)


import BasisFunctions: measure, restrict
import BasisFunctions: innerproduct, innerproduct_native

hasmeasure(::ExtensionFrame) = true
measure(f::ExtensionFrame) = BasisFunctions.submeasure(measure(basis(f)), support(f))

innerproduct_native(f1::ExtensionFrame, i, f2::ExtensionFrame, j, measure::SubMeasure; options...) =
    innerproduct(superdict(f1), i, superdict(f2), j, measure; options...)

# This routine will be called if we have two mapped dictionaries, where the measure is a subset of a
# mapped measure. We check for this case and undo the mapping.
innerproduct_native(dict1::MappedDict, i, dict2::MappedDict, j, measure::SubMeasure; options...) =
    _innerproduct_native(support(measure), supermeasure(measure), dict1, i, dict2, j, measure; options...)

function _innerproduct_native(domain::AbstractInterval, superμ::MappedMeasure, dict1::MappedDict, i, dict2::MappedDict, j, measure::SubMeasure; options...)
    if iscompatible(dict1, dict2) && iscompatible(mapping(dict1), mapping(superμ))
        supermap = mapping(superμ)
        newdomain = Interval(applymap(inv(supermap), infimum(domain)), applymap(inv(supermap), supremum(domain)))
        innerproduct_native(superdict(dict1), i, superdict(dict2), j, SubMeasure(supermeasure(superμ), newdomain))
    else
        dict_default_innerproduct(dict1, i, dict2, j, measure; options...)
    end
end

function _innerproduct_native(domain::UnionDomain, superμ::MappedMeasure, dict1::MappedDict, i, dict2::MappedDict, j, measure::SubMeasure; options...)
    z = zero(coefficienttype(dict1))
    for d in elements(domain)
        z += _innerproduct_native(d, superμ, dict1, i, dict2, j, measure; options...)
    end
    z
end

_innerproduct_native(domain, superμ, dict1::MappedDict, i, dict2::MappedDict, j, measure::SubMeasure; options...) =
    default_dict_innerproduct(dict1, i, dict2, j, measure; options...)

extensiondual(dict::ExtensionFrame, measure; options...) =
    extensionframe(support(dict), gramdual(superdict(dict), supermeasure(measure); options...),)

function BasisFunctions.gramdual(dict::ExtensionFrame, measure::Measure; options...)
    @debug "Are you sure you want `dualtype=gramdual` and not `extensiondual`"
    @warn "Changing to extensiondual"
    extensiondual(dict, measure; options...)
    # BasisFunctions.default_gramdual(dict, measure; options...)
end

superdict(dict::ExtensionFrameTensor) = TensorProductDict(map(superdict, elements(dict))...)
support(dict::ExtensionFrameTensor) = ProductDomain((map(support, elements(dict)))...)

## Printing

string(f::ExtensionFrame) = name(f) * " of " * name(f.basis)

modifiersymbol(dict::ExtensionFrame) = PrettyPrintSymbol{:𝔼}(dict)

string(s::PrettyPrintSymbol{:𝔼}) = _string(s, s.object)
_string(s::PrettyPrintSymbol{:𝔼}, dict::ExtensionFrame) =
    "Extension frame, from $(support(dict)) to $(support(superdict(dict)))"



##################
# platform
##################

SolverStyle(dict::ExtensionFrame, ::OversamplingStyle) = AZStyle()

GridSampling(dgs::GridBasis, grid::AbstractGrid, domain::Domain, scaling) =
    GridSampling(GridBasis{coefficienttype(dgs)}(subgrid(grid, domain)), scaling=scaling)
