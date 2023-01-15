
export ExtensionFrame,
    extensionframe,
    ⇐, →


using BasisFunctions: default_in_support, unsafe_eval_element,
    default_dict_innerproduct, Measure

import BasisFunctions: superdict, support, similardictionary, isbasis,
    isbiorthogonal, isorthogonal, isorthonormal, hasinterpolationgrid, hastransform,
    hasantiderivative, name, string, strings, dict_in_support, iscompatible, *,
    unsafe_eval_element1, interpolation_grid, measure, restrict, dict_innerproduct,
    dict_innerproduct_native, gramdual, modifiersymbol, GridSampling, in_support,
    points, support

include("submeasure.jl")

export ExtensionFrame
"""
    struct ExtensionFrame{S,T} <: DerivedDict{S,T}

An ExtensionFrame is the restriction of a basis to a subset of its domain. This
results in a frame that implicitly represents extensions of functions on the
smaller set to the larger set.

Created by the function `extensionframe`.

See also [`extensionframe`](@ref)
```jldocs
julia> extensionframe(Fourier(10),0.0..0.5)
```
"""
struct ExtensionFrame{S,T,DICT<:Dictionary{S,T},DOM<:Domain} <: DerivedDict{S,T}
    domain      ::  DOM
    basis       ::  DICT

    function ExtensionFrame{S,T,DICT,DOM}(domain::DOM, basis::DICT) where {S,T,DICT<:Dictionary{S,T},DOM<:Domain}
        # @assert isbasis(basis)
        new(domain, basis)
    end
end

ExtensionFrame(domain::Domain, basis::Dictionary{S,T}) where {S,T} =
    ExtensionFrame{S,T,typeof(basis),typeof(domain)}(domain, basis)
ExtensionFrame(basis::Dictionary, domain::Domain) =
    ExtensionFrame(domain, basis)

superdict(f::ExtensionFrame) = f.basis

export basis
"""
    basis(f::ExtensionFrame)

Return the dictionary on the bouding box.
"""
basis(f::ExtensionFrame) = superdict(f)
support(f::ExtensionFrame) = f.domain
support(f::ExtensionFrame, i) = f.domain

similardictionary(f::ExtensionFrame, dict::Dictionary) = ExtensionFrame(support(f), dict)

isbasis(f::ExtensionFrame) = false
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

iscompatible(d1::ExtensionFrame, d2::ExtensionFrame) = iscompatible(basis(d1),basis(d2))

function BasisFunctions.expansion_multiply(d1::ExtensionFrame, d2::ExtensionFrame, args...)
    @assert iscompatible(d1,d2)
    (mset, mcoef) = BasisFunctions.expansion_multiply(basis(d1),basis(d2),args...)
    df = ExtensionFrame(support(d1) ∩ support(d2),mset)
    (df, mcoef)
end

unsafe_eval_element1(dict::ExtensionFrame, idx::Int, x) =
    in_support(dict, idx, x) ? unsafe_eval_element(basis(dict), idx, x) : zero(codomaintype(dict))
unsafe_eval_element1(dict::ExtensionFrame, idx, x) =
    in_support(dict, idx, x) ? unsafe_eval_element(basis(dict), idx, x) : zero(codomaintype(dict))
# unsafe_eval_element1(dict::ExtensionFrame, idx, grid::AbstractGrid) =
    # _default_unsafe_eval_element_in_grid(dict, idx, grid)
in_support(dict::ExtensionFrame, idx::Int, x) =
    default_in_support(dict, idx, x)
in_support(dict::ExtensionFrame, idx, x) =
    default_in_support(dict, idx, x)

function interpolation_grid(f::ExtensionFrame)
    @warn "interpolation_grid called on extensionframe $f"
    subgrid(interpolation_grid(basis(f)),support(f))
end

export extensionframe
"""
    extensionframe(basis::Dictionary, domain::Domain)

    extensionframe(domain::Domain, basis::Dictionary)

Make an ExtensionFrame, but match tensor product domains with tensor product sets
in a suitable way.

# examples
For example: an interval ⊗ a disk (= a cylinder) combined with a 3D Fourier series, leads to a
tensor product of a Fourier series on the interval ⊗ a 2D Fourier series on the disk.
```jldocs
julia> extensionframe(Fourier(10)^2,Disk(.4, SVector(0.5,0.5)))
Dictionary 𝔼(F ⊗ F)

𝔼   :   Extension frame, from A mapped domain based on the 2-dimensional unit ball to 0.0..1.0 (Unit) x 0.0..1.0 (Unit)
F   :   Fourier series
            ↳ length = 10
            ↳ Float64 -> Complex{Float64}
            ↳ support = 0.0..1.0 (Unit)


julia> extensionframe(Fourier(10)^2,(0.0..0.5)^2)
Dictionary 𝔼(F) ⊗ 𝔼(F)

𝔼   :   Extension frame, from 0.0..0.5 to 0.0..1.0 (Unit)
F   :   Fourier series
            ↳ length = 10
            ↳ Float64 -> Complex{Float64}
            ↳ support = 0.0..1.0 (Unit)

```
"""
extensionframe(domain::Domain, basis::Dictionary) = ExtensionFrame(domain, basis)
extensionframe(basis::Dictionary, domain::Domain) = extensionframe(domain, basis)

"Make an extension frame (symbol ⇐ is \\Leftarrow)"
⇐(basis::Dictionary, domain::Domain) = extensionframe(basis, domain)

import BasisFunctions: →
→(domain::Domain, basis::Dictionary) = extensionframe(basis, domain)

function extensionframe(domain::ProductDomain, basis::TensorProductDict)
    ExtensionFrames = Dictionary[]
    dc = 1
    for i = 1:ncomponents(domain)
        el = component(domain, i)
        range = dc:dc+dimension(el)-1
        push!(ExtensionFrames, ExtensionFrame(el, tensorproduct(component(basis, j) for j in range)...))
        dc += dimension(el)
    end
    tensorproduct(ExtensionFrames...)
end

using BasisFunctions: TupleProductDict

const ExtensionFrameTensor = Union{TupleProductDict{D,NTuple{D,DICT}} where D where {DICT<:ExtensionFrame},
    TupleProductDict{2,Tuple{DICT1,DICT2}} where {DICT1<:ExtensionFrame,DICT2<:ExtensionFrame},
    TupleProductDict{3,Tuple{DICT1,DICT2,DICT3}} where {DICT1<:ExtensionFrame,DICT2<:ExtensionFrame,DICT3<:ExtensionFrame},
    TupleProductDict{4,Tuple{DICT1,DICT2,DICT3,DICT4}} where {DICT1<:ExtensionFrame,DICT2<:ExtensionFrame,DICT3<:ExtensionFrame,DICT4<:ExtensionFrame}}
const ExtensionFrameSuper = Union{<:ExtensionFrame,<:ExtensionFrameTensor}
export ExtensionFrameTensor, ExtensionFrameSuper



hasmeasure(::ExtensionFrame) = true
measure(f::ExtensionFrame) = submeasure(measure(basis(f)), support(f))

dict_innerproduct_native(f1::ExtensionFrame, i, f2::ExtensionFrame, j, measure; options...) =
    dict_innerproduct(superdict(f1), i, superdict(f2), j, measure; options...)

# This routine will be called if we have two mapped dictionaries, where the measure is a subset of a
# mapped measure. We check for this case and undo the mapping.
dict_innerproduct_native(dict1::MappedDict, i, dict2::MappedDict, j, measure; options...) =
    _dict_innerproduct_native(support(measure), supermeasure(measure), dict1, i, dict2, j, measure; options...)

function _dict_innerproduct_native(domain::AbstractInterval, superμ, dict1::MappedDict, i, dict2::MappedDict, j, measure; options...)
    if BasisFunctions.ismappedmeasure(superμ) && iscompatible(dict1, dict2) && iscompatible(forward_map(dict1), forward_map(superμ))
        supermap = forward_map(superμ)
        newdomain = Interval(applymap(inv(supermap), infimum(domain)), applymap(inv(supermap), supremum(domain)))
        dict_innerproduct_native(superdict(dict1), i, superdict(dict2), j, submeasure(supermeasure(superμ), newdomain))
    else
        default_dict_innerproduct(dict1, i, dict2, j, measure; options...)
    end
end

function _dict_innerproduct_native(domain::UnionDomain, superμ, dict1::MappedDict, i, dict2::MappedDict, j, measure; options...)
    z = zero(coefficienttype(dict1))
    for d in components(domain)
        z += _dict_innerproduct_native(d, superμ, dict1, i, dict2, j, measure; options...)
    end
    z
end

_dict_innerproduct_native(domain, superμ, dict1::MappedDict, i, dict2::MappedDict, j, measure; options...) =
    BasisFunctions.default_dict_innerproduct(dict1, i, dict2, j, measure; options...)

export extensiondual
"""
    extensiondual(dict::ExtensionFrame, measure::Measure; options...)

Return the extension frame of the dual of the dictionary on the boundingbox.
"""
extensiondual(dict::ExtensionFrame, measure::Measure; options...) =
    extensionframe(support(dict), gramdual(superdict(dict), supermeasure(measure); options...),)
extensiondual(dict::ExtensionFrameTensor, measure::Measure; options...) =
    TensorProductDict(map((dicti,measurei)->extensionframe(support(dicti), gramdual(superdict(dicti), supermeasure(measurei); options...)),
        components(dict), components(measure))...)

function gramdual(dict::ExtensionFrame, measure::Measure; verbose=false, options...)
    verbose && println("WARN: Are you sure you want `dualtype=gramdual` and not `extensiondual`?")
    default_gramdual(dict, measure; verbose, options...)
end

superdict(dict::ExtensionFrameTensor) = TensorProductDict(map(superdict, components(dict))...)
basis(dict::ExtensionFrameTensor) = TensorProductDict(map(basis, components(dict))...)
support(dict::ExtensionFrameTensor) = ProductDomain((map(support, components(dict)))...)

## Printing

Display.object_parentheses(d::ExtensionFrame) = true
Display.stencil_parentheses(d::ExtensionFrame) = true
Display.displaystencil(d::ExtensionFrame) = _stencil(d, superdict(d), support(d))
_stencil(d::ExtensionFrame, dict, domain) = [domain, " → ", dict]

GridSampling(dgs::GridBasis, grid::AbstractGrid, domain::Domain, scaling) =
    GridSampling(GridBasis{coefficienttype(dgs)}(subgrid(grid, domain)), scaling=scaling)
