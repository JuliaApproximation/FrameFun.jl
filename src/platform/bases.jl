## Platforms for certain bases and frames


struct FourierPlatform{T} <: BasisPlatform
end

FourierPlatform() = FourierPlatform{Float64}()

dictionary(p::FourierPlatform{T}, n) where {T} = Fourier{T}(n)

first_parameters(p::FourierPlatform) = (8,8)

SolverStyle(p::FourierPlatform, ::OversamplingStyle) = DualStyle()

measure(p::FourierPlatform{T}) where {T} = FourierMeasure{T}()

discrete_gram_measure(p::FourierPlatform, n, L;
            dict = dictionary(p, n), normalizedsampling = false, options...) =
    normalizedsampling ? BasisFunctions.DiracCombProbabilityMeasure(oversampling_grid(dict, L)) :
                        discretemeasure(oversampling_grid(dict, L))


struct ChebyshevPlatform{T} <: BasisPlatform
end

ChebyshevPlatform() = ChebyshevPlatform{Float64}()

dictionary(p::ChebyshevPlatform{T}, n) where {T} = ChebyshevT{T}(n)

first_parameters(p::ChebyshevPlatform) = (8,8)

measure(platform::ChebyshevPlatform{T}) where {T} = ChebyshevMeasure{T}()



struct FourierExtensionPlatform <: FramePlatform
    basisplatform   ::  FourierPlatform
    domain          ::  Domain

    function FourierExtensionPlatform(basisplatform, domain::Domain{T}) where {T}
        @assert issubset(domain, UnitInterval{T}())
        new(basisplatform, domain)
    end
end

function FourierExtensionPlatform(domain::Domain{T}) where {T}
    basisplatform = FourierPlatform{T}()
    FourierExtensionPlatform(basisplatform, domain)
end

dictionary(p::FourierExtensionPlatform, n) =
    extensionframe(p.domain, dictionary(p.basisplatform, n))

measure(platform::FourierExtensionPlatform) = measure(dictionary(platform, 1))

discrete_gram_measure(p::FourierExtensionPlatform, n, L;
            dict = dictionary(p, n), normalizedsampling = false, options...) =
    normalizedsampling ?
    restrict(BasisFunctions.DiracCombProbabilityMeasure(oversampling_grid(dictionary(p.basisplatform, n), L)), p.domain) :
    restrict(discretemeasure(oversampling_grid(dictionary(p.basisplatform, n), L)), p.domain)


struct ExtensionFramePlatform <: FramePlatform
    basisplatform   ::  Platform
    domain          ::  Domain
end

dictionary(p::ExtensionFramePlatform, n) =
    ExtensionFrame(p.domain, dictionary(p.basisplatform, n))

measure(platform::ExtensionFramePlatform) = restrict(measure(platform.basisplatform), platform.domain)

struct OPSExtensionFramePlatform <: FramePlatform
    basisplatform   ::Platform
    domain          ::Domain
    function OPSExtensionFramePlatform(dict::BasisFunctions.OPS, domain::Domain1d)
        new(ModelPlatform(dict),domain)
    end
end


dictionary(p::OPSExtensionFramePlatform, n) =
    ExtensionFrame(p.domain, dictionary(p.basisplatform, n))

measure(platform::OPSExtensionFramePlatform) = restrict(measure(platform.basisplatform), platform.domain)


discretemeasure(::DiscreteStyle, p::OPSExtensionFramePlatform, n, L;
            dict = dictionary(p, n), options...) =
    _opsnodesmeasure(resize(dict,L))
_opsnodesmeasure(frame::ExtensionFrame) = restrict(gauss_rule(basis(frame)), support(frame))
_opsnodesmeasure(ops::BasisFunctions.OrthogonalPolynomials) = gauss_rule(ops)
