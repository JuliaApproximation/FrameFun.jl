
using BasisFunctions: FFreq, DiscreteWeight, DiscreteProductWeight,
    innerproduct_fourier_part, default_dict_innerproduct, Measure, Weight,
	OuterProductArray
import BasisFunctions: quadweights, stencilarray, name, weights,
	grid, subindices, isnormalized, supermeasure,
    support, unsafe_weightfun, strings, restrict, discretemeasure
export submeasure, DiscreteTensorSubWeight, SubWeight, DiscreteSubWeight

"""
    struct SubWeight{M,D,T} <: Weight{T}

Submeasure of a Weight.
"""
struct SubWeight{M<:Weight,D<:Domain,T} <: Weight{T}
    measure     ::  M
    domain      ::  D

    function SubWeight(measure::Weight, domain::Domain)
        @assert eltype(domain) == eltype(support(measure))
        new{typeof(measure), typeof(domain), eltype(domain)}(measure, domain)
    end
end

"""
    submeasure(measure::Weight, domain::Domain)

Restrict a measure to a subset of its support
"""
submeasure(measure::Weight, domain::Domain) = SubWeight(measure, domain)
name(m::SubWeight) = "Restriction of a measure"
supermeasure(measure::SubWeight) = measure.measure
support(measure::SubWeight) = measure.domain
unsafe_weightfun(m::SubWeight, x) = unsafe_weightfun(supermeasure(m), x)

"""
    restrict(measure::Weight, domain::Domain)

Restrict a measure to a subset of its support
"""
restrict(measure::Weight, domain::Domain) = submeasure(measure, domain)
strings(m::SubWeight) = (name(m), (string(support(m)),), strings(supermeasure(m)))
submeasure(measure::ProductWeight, domain::ProductDomain) = ProductWeight(map(submeasure, components(measure), components(domain))...)


quadweights(grid::SubGrid, measure::SubWeight) =
    quadweights(supergrid(grid), supermeasure(measure))[subindices(grid)]

function dict_innerproduct_native(b1::Fourier, i::FFreq, b2::Fourier, j::FFreq, m::SubWeight{<:FourierWeight}; options...)
	d = support(m)
	if d isa AbstractInterval
		innerproduct_fourier_part(b1, i, b2, j, infimum(d), supremum(d))
	else
		default_dict_innerproduct(b1, i, b2, j, m; options...)
	end
end

"""
    struct DiscreteSubWeight{T,M<:DiscreteWeight{T},G<:AbstractGrid} <: DiscreteWeight{T}

Submeasure of a discrete measure.
"""
struct DiscreteSubWeight{T,M<:DiscreteWeight{T},G<:AbstractGrid} <: DiscreteWeight{T}
    supermeasure   :: M
    subgrid        :: G
    DiscreteSubWeight(measure::DiscreteWeight{T}, grid::AbstractGrid) where T = new{T,typeof(measure),typeof(grid)}(measure, grid)
    DiscreteSubWeight(measure::DiscreteSubWeight{T}, grid::AbstractGrid) where T = error()
end

submeasure(measure::DiscreteWeight, domain::Domain) = _discretesubmeasure(subgrid(points(measure), domain), weights(measure))
subindices(measure::DiscreteSubWeight) = subindices(measure.subgrid)
supermeasure(measure::DiscreteSubWeight) = measure.supermeasure
discretemeasure(grid::Union{SubGrid,ProductSubGrid}) = DiscreteSubWeight(discretemeasure(supergrid(grid)), grid)
_discretesubmeasure(grid::Union{SubGrid,ProductSubGrid},weights) = DiscreteSubWeight(discretemeasure(supergrid(grid),weights), grid)
restrict(measure::DiscreteWeight, domain::Domain) = DiscreteSubWeight(measure, subgrid(points(measure), domain))
weights(measure::DiscreteSubWeight) = subweights(measure, subindices(measure), weights(supermeasure(measure)))
subweights(_, subindices, w) = w[subindices]
points(measure::DiscreteSubWeight) = measure.subgrid
isnormalized(::DiscreteSubWeight) = false
name(m::DiscreteSubWeight) = "Restriction of a "*name(supermeasure(m))

function discretemeasure(grid::Union{SubGrid,ProductSubGrid}, weights::Ones)
    @assert size(grid) == size(weights)
    DiscreteSubWeight(discretemeasure(supergrid(grid)), grid)
end

"""
    const DiscreteTensorSubWeight{T,G,W} = DiscreteSubWeight{T,M,G} where {T,M<:DiscreteProductWeight,G<:ProductSubGrid}

A tensor product of discrete submeasures.
"""
const DiscreteTensorSubWeight{T,G,W} = DiscreteSubWeight{T,M,G} where {T,M<:DiscreteProductWeight,G<:ProductSubGrid}
name(m::DiscreteTensorSubWeight) = "Tensor of submeasures (supermeasure:"*name(supermeasure(m))
components(m::DiscreteTensorSubWeight) = map(_discretesubmeasure,components(points(m)),components(weights(supermeasure(m))))
component(m::DiscreteTensorSubWeight, i) = _discretesubmeasure(component(points(m),i),component(weights(supermeasure(m)),i))
weights(m::DiscreteTensorSubWeight) = OuterProductArray(map(weights, components(m))...)


function stencilarray(m::DiscreteTensorSubWeight)
    A = Any[]
    push!(A, component(m,1))
    for i = 2:length(components(m))
        push!(A," ⊗ ")
        push!(A, component(m,i))
    end
    A
end
