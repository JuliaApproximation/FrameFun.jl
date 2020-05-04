## Platforms for certain bases and frames
module FrameFunPlatforms
using ..Platforms
import ..Platforms: dictionary, SolverStyle, measure, dualdictionary,
    correctparamformat, unsafe_dictionary
using ..FrameFunInterface
import ..FrameFunInterface: discrete_gram_measure
using BasisFunctions, ..ExtensionFrames, DomainSets
using BasisFunctions: AbstractMeasure

export FourierPlatform
"""
    struct FourierPlatform{T} <: BasisPlatform

A platform of Fourier series on [0,1]
"""
struct FourierPlatform{T} <: BasisPlatform
end

correctparamformat(::FourierPlatform, ::Int) = true
correctparamformat(::FourierPlatform, _) = false

FourierPlatform() = FourierPlatform{Float64}()

unsafe_dictionary(p::FourierPlatform{T}, n) where {T} = Fourier{T}(n)

SolverStyle(p::FourierPlatform, ::OversamplingStyle) = DualStyle()

measure(p::FourierPlatform{T}) where {T} = FourierMeasure{T}()

discrete_gram_measure(p::FourierPlatform, n, L;
            dict = dictionary(p, n), normalizedsampling = false, options...) =
    normalizedsampling ? BasisFunctions.NormalizedDiracComb(oversampling_grid(dict, L)) :
                        discretemeasure(oversampling_grid(dict, L))

export ChebyshevPlatform
"""
    struct ChebyshevPlatform{T} <: BasisPlatform

A platform of ChebyshevT dictionaries on [-1,1]
"""
struct ChebyshevPlatform{T} <: BasisPlatform
end

correctparamformat(::ChebyshevPlatform, ::Int) = true
correctparamformat(::ChebyshevPlatform, _) = false

ChebyshevPlatform() = ChebyshevPlatform{Float64}()

unsafe_dictionary(p::ChebyshevPlatform{T}, n) where {T} = ChebyshevT{T}(n)

measure(platform::ChebyshevPlatform{T}) where {T} = ChebyshevMeasure{T}()

export FourierExtensionPlatform
"""
    struct FourierExtensionPlatform <: FramePlatform

A platform of fourier extensions with [0,1] as bouding box.
"""
struct FourierExtensionPlatform <: FramePlatform
    basisplatform   ::  FourierPlatform
    domain          ::  Domain

    function FourierExtensionPlatform(basisplatform, domain::Domain{T}) where {T}
        # Only supported for intervals
        # @assert issubset(domain, UnitInterval{T}())
        new(basisplatform, domain)
    end
end

correctparamformat(platform::FourierExtensionPlatform, param) =
    correctparamformat(platform.basisplatform, param)

function FourierExtensionPlatform(domain::Domain{T}) where {T}
    basisplatform = FourierPlatform{T}()
    FourierExtensionPlatform(basisplatform, domain)
end

unsafe_dictionary(p::FourierExtensionPlatform, n) =
    extensionframe(p.domain, unsafe_dictionary(p.basisplatform, n))

dualdictionary(platform::FourierExtensionPlatform, param, measure::AbstractMeasure; options...) =
   extensionframe(dualdictionary(platform.basisplatform, param, supermeasure(measure); options...), platform.domain)

measure(platform::FourierExtensionPlatform) = measure(dictionary(platform, 1))

discrete_gram_measure(p::FourierExtensionPlatform, n, L;
            dict = dictionary(p, n), normalizedsampling = false, options...) =
    normalizedsampling ?
    restrict(BasisFunctions.NormalizedDiracComb(oversampling_grid(dictionary(p.basisplatform, n), L)), p.domain) :
    restrict(discretemeasure(oversampling_grid(dictionary(p.basisplatform, n), L)), p.domain)

export OPSExtensionFramePlatform
"""
    struct OPSExtensionFramePlatform <: FramePlatform

A platform of classical orthogonal polynomial extensions. The bouding box is [-1,1]
"""
struct OPSExtensionFramePlatform <: FramePlatform
    basisplatform   ::Platform
    domain          ::Domain
    function OPSExtensionFramePlatform(dict::BasisFunctions.OPS, domain::Domain1d)
        new(platform(dict),domain)
    end
end

correctparamformat(p::OPSExtensionFramePlatform, param) =
    correctparamformat(p.basisplatform, param)


unsafe_dictionary(p::OPSExtensionFramePlatform, n) =
    ExtensionFrame(p.domain, unsafe_dictionary(p.basisplatform, n))

measure(platform::OPSExtensionFramePlatform) = restrict(measure(platform.basisplatform), platform.domain)

dualdictionary(platform::OPSExtensionFramePlatform, param, measure::AbstractMeasure; options...) =
   extensionframe(dualdictionary(platform.basisplatform, param, supermeasure(measure); options...), platform.domain)

discretemeasure(::DiscreteStyle, p::OPSExtensionFramePlatform, n, L;
            dict = dictionary(p, n), options...) =
    _opsnodesmeasure(resize(dict,L))
_opsnodesmeasure(frame::ExtensionFrame) = restrict(gauss_rule(basis(frame)), support(frame))
_opsnodesmeasure(ops::BasisFunctions.OrthogonalPolynomials) = gauss_rule(ops)
end
