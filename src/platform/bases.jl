## Platforms for certain bases and frames


struct FourierPlatform{T} <: BasisPlatform
end

FourierPlatform() = FourierPlatform{Float64}()

dictionary(p::FourierPlatform{T}, n) where {T} = Fourier{T}(n)

first_parameters(p::FourierPlatform) = (8,8)

SolverStyle(p::FourierPlatform, ::OversamplingStyle) = DualStyle()

measure(p::FourierPlatform{T}) where {T} = FourierMeasure{T}()

# discrete_gram_measure(p::FourierPlatform, n, L; dict = dictionary(p, n), options...) =
#     BasisFunctions.DiracCombProbabilityMeasure(oversampling_grid(dict, L))


struct ChebyshevPlatform{T} <: BasisPlatform
end

ChebyshevPlatform() = ChebyshevPlatform{Float64}()

dictionary(p::ChebyshevPlatform{T}, n) where {T} = ChebyshevT{T}(n)

first_parameters(p::ChebyshevPlatform) = (8,8)

measure(platform::ChebyshevPlatform{T}) where {T} = ChebyshevMeasure{T}()


"""
A `ModelPlatform` is a platform based on a model dictionary. The platform is
defined by resizing the dictionary, using its own implementation of `resize`.
All other operations are the defaults for the model dictionary.

This platform is convenient to compute adaptive approximations based on an
example of a dictionary from the desired family.
"""
struct ModelPlatform <: Platform
    model   ::  Dictionary
end

model(p::ModelPlatform) = p.model

dictionary(p::ModelPlatform, n) = resize(model(p), n)

first_sizeparameter(p::ModelPlatform) = typeof(model(p)) <: Dictionary1d ? length(model(p)) : size(model(p))

SamplingStyle(p::ModelPlatform) = SamplingStyle(dictionary(p, first_sizeparameter(p)))

SolverStyle(p::ModelPlatform, dstyle::SamplingStyle) = SolverStyle(model(p), dstyle)

measure(platform::ModelPlatform) = measure(model(platform))


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

# discrete_gram_measure(p::FourierExtensionPlatform, n, L; dict = dictionary(p, n), options...) =
#     restrict(BasisFunctions.DiracCombProbabilityMeasure(oversampling_grid(dictionary(p.basisplatform, n), L)), p.domain)


struct ExtensionFramePlatform <: FramePlatform
    basisplatform   ::  Platform
    domain          ::  Domain
end

dictionary(p::ExtensionFramePlatform, n) =
    ExtensionFrame(p.domain, dictionary(p.basisplatform, n))

measure(platform::ExtensionFramePlatform) = restrict(measure(platform.basisplatform), platform.domain)
