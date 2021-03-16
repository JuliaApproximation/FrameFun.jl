
using BasisFunctions, GridArrays, DomainSets
using BasisFunctions: FFreq, DiscreteWeight, DiscreteProductWeight,
    innerproduct_fourier_part, default_dict_innerproduct, Measure, Weight,
	OuterProductArray
import BasisFunctions: quadweights, stencilarray, name, element, elements,
    iscomposite, weights, grid, subindices, isnormalized, supermeasure,
    support, unsafe_weightfun, strings, innerproduct_native, restrict, discretemeasure
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
submeasure(measure::ProductWeight, domain::ProductDomain) = ProductWeight(map(submeasure, elements(measure), elements(domain))...)


quadweights(grid::AbstractSubGrid, measure::SubWeight) =
    quadweights(supergrid(grid), supermeasure(measure))[subindices(grid)]

function innerproduct_native(b1::Fourier, i::FFreq, b2::Fourier, j::FFreq, m::SubWeight{<:FourierWeight}; options...)
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
discretemeasure(grid::Union{AbstractSubGrid,TensorSubGrid}) = DiscreteSubWeight(discretemeasure(supergrid(grid)), grid)
_discretesubmeasure(grid::Union{AbstractSubGrid,TensorSubGrid},weights) = DiscreteSubWeight(discretemeasure(supergrid(grid),weights), grid)
restrict(measure::DiscreteWeight, domain::Domain) = DiscreteSubWeight(measure, subgrid(points(measure), domain))
weights(measure::DiscreteSubWeight) = subweights(measure, subindices(measure), weights(supermeasure(measure)))
subweights(_, subindices, w) = w[subindices]
points(measure::DiscreteSubWeight) = measure.subgrid
isnormalized(::DiscreteSubWeight) = false
name(m::DiscreteSubWeight) = "Restriction of a "*name(supermeasure(m))

"""
    const DiscreteTensorSubWeight{T,G,W} = DiscreteSubWeight{T,M,G} where {T,M<:DiscreteProductWeight,G<:TensorSubGrid}

A tensor product of discrete submeasures.
"""
const DiscreteTensorSubWeight{T,G,W} = DiscreteSubWeight{T,M,G} where {T,M<:DiscreteProductWeight,G<:TensorSubGrid}
name(m::DiscreteTensorSubWeight) = "Tensor of submeasures (supermeasure:"*name(supermeasure(m))
elements(m::DiscreteTensorSubWeight) = map(_discretesubmeasure,elements(points(m)),elements(weights(supermeasure(m))))
element(m::DiscreteTensorSubWeight, i) = _discretesubmeasure(element(points(m),i),element(weights(supermeasure(m)),i))
iscomposite(m::DiscreteTensorSubWeight) = true
weights(m::DiscreteTensorSubWeight) = OuterProductArray(map(weights, elements(m))...)


function stencilarray(m::DiscreteTensorSubWeight)
    A = Any[]
    push!(A, element(m,1))
    for i = 2:length(elements(m))
        push!(A," ⊗ ")
        push!(A, element(m,i))
    end
    A
end
