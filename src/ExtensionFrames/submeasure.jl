
using BasisFunctions, GridArrays, DomainSets
using BasisFunctions: FFreq, DiscreteMeasure, DiscreteProductMeasure,
    innerproduct_fourier_part, default_dict_innerproduct, AbstractMeasure,
	OuterProductArray
import BasisFunctions: quadweights, stencilarray, name, element, elements,
    iscomposite, weights, grid, subindices, isnormalized, supermeasure,
    support, unsafe_weight, strings, innerproduct_native, restrict, discretemeasure
export submeasure, DiscreteTensorSubMeasure, SubMeasure, DiscreteSubMeasure

"""
    struct SubMeasure{M,D,T} <: Measure{T}

Submeasure of a Measure.
"""
struct SubMeasure{M<:Measure,D<:Domain,T} <: Measure{T}
    measure     ::  M
    domain      ::  D

    function SubMeasure(measure::Measure, domain::Domain)
        @assert eltype(domain) == eltype(support(measure))
        new{typeof(measure), typeof(domain), eltype(domain)}(measure, domain)
    end
end

"""
    submeasure(measure::Measure, domain::Domain)

Restrict a measure to a subset of its support
"""
submeasure(measure::Measure, domain::Domain) = SubMeasure(measure, domain)
name(m::SubMeasure) = "Restriction of a measure"
supermeasure(measure::SubMeasure) = measure.measure
support(measure::SubMeasure) = measure.domain
unsafe_weight(m::SubMeasure, x) = unsafe_weight(supermeasure(m), x)

"""
    restrict(measure::Measure, domain::Domain)

Restrict a measure to a subset of its support
"""
restrict(measure::Measure, domain::Domain) = submeasure(measure, domain)
strings(m::SubMeasure) = (name(m), (string(support(m)),), strings(supermeasure(m)))
submeasure(measure::ProductMeasure, domain::ProductDomain) = ProductMeasure(map(submeasure, elements(measure), elements(domain))...)


quadweights(grid::AbstractSubGrid, measure::SubMeasure) =
    quadweights(supergrid(grid), supermeasure(measure))[subindices(grid)]

function innerproduct_native(b1::Fourier, i::FFreq, b2::Fourier, j::FFreq, m::SubMeasure{<:FourierMeasure}; options...)
	d = support(m)
	if d isa AbstractInterval
		innerproduct_fourier_part(b1, i, b2, j, infimum(d), supremum(d))
	else
		default_dict_innerproduct(b1, i, b2, j, m; options...)
	end
end

"""
    struct DiscreteSubMeasure{T,M<:DiscreteMeasure{T},G<:AbstractGrid} <: DiscreteMeasure{T}

Submeasure of a discrete measure.
"""
struct DiscreteSubMeasure{T,M<:DiscreteMeasure{T},G<:AbstractGrid} <: DiscreteMeasure{T}
    supermeasure   :: M
    subgrid        :: G
    DiscreteSubMeasure(measure::DiscreteMeasure{T}, grid::AbstractGrid) where T = new{T,typeof(measure),typeof(grid)}(measure, grid)
    DiscreteSubMeasure(measure::DiscreteSubMeasure{T}, grid::AbstractGrid) where T = error()
end

submeasure(measure::DiscreteMeasure, domain::Domain) = _discretesubmeasure(subgrid(points(measure), domain), weights(measure))
subindices(measure::DiscreteSubMeasure) = subindices(measure.subgrid)
supermeasure(measure::DiscreteSubMeasure) = measure.supermeasure
discretemeasure(grid::Union{AbstractSubGrid,TensorSubGrid}) = DiscreteSubMeasure(discretemeasure(supergrid(grid)), grid)
_discretesubmeasure(grid::Union{AbstractSubGrid,TensorSubGrid},weights) = DiscreteSubMeasure(discretemeasure(supergrid(grid),weights), grid)
restrict(measure::DiscreteMeasure, domain::Domain) = DiscreteSubMeasure(measure, subgrid(points(measure), domain))
weights(measure::DiscreteSubMeasure) = subweights(measure, subindices(measure), weights(supermeasure(measure)))
subweights(_, subindices, w) = w[subindices]
points(measure::DiscreteSubMeasure) = measure.subgrid
isnormalized(::DiscreteSubMeasure) = false
name(m::DiscreteSubMeasure) = "Restriction of a "*name(supermeasure(m))

"""
    const DiscreteTensorSubMeasure{T,G,W} = DiscreteSubMeasure{T,M,G} where {T,M<:DiscreteProductMeasure,G<:TensorSubGrid}

A tensor product of discrete submeasures.
"""
const DiscreteTensorSubMeasure{T,G,W} = DiscreteSubMeasure{T,M,G} where {T,M<:DiscreteProductMeasure,G<:TensorSubGrid}
name(m::DiscreteTensorSubMeasure) = "Tensor of submeasures (supermeasure:"*name(supermeasure(m))
elements(m::DiscreteTensorSubMeasure) = map(_discretesubmeasure,elements(points(m)),elements(weights(supermeasure(m))))
element(m::DiscreteTensorSubMeasure, i) = _discretesubmeasure(element(points(m),i),element(weights(supermeasure(m)),i))
iscomposite(m::DiscreteTensorSubMeasure) = true
weights(m::DiscreteTensorSubMeasure) = OuterProductArray(map(weights, elements(m))...)


function stencilarray(m::DiscreteTensorSubMeasure)
    A = Any[]
    push!(A, element(m,1))
    for i = 2:length(elements(m))
        push!(A," ⊗ ")
        push!(A, element(m,i))
    end
    A
end
