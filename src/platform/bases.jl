## Platforms for certain bases and frames


struct FourierPlatform{T} <: BasisPlatform
end

FourierPlatform() = FourierPlatform{Float64}()

dictionary(p::FourierPlatform{T}, n) where {T} = Fourier{T}(n)

first_parameters(p::FourierPlatform) = (8,8)

SolverStyle(p::FourierPlatform, ::OversamplingStyle) = DualStyle()

discrete_normalization(p::FourierPlatform, n, L; S = samplingoperator(p, n, L)) =
    quadraturenormalization(S)

functionspace(p::FourierPlatform{T}) where {T} = BasisFunctions.FourierSpace{T}()


struct ChebyshevPlatform{T} <: BasisPlatform
end

ChebyshevPlatform() = ChebyshevPlatform{Float64}()

dictionary(p::ChebyshevPlatform{T}, n) where {T} = ChebyshevT{T}(n)

first_parameters(p::ChebyshevPlatform) = (8,8)

discrete_normalization(p::ChebyshevPlatform, n, L; S = samplingoperator(p, n, L)) =
    quadraturenormalization(S)

functionspace(p::ChebyshevPlatform{T}) where {T} = BasisFunctions.ChebyshevSpace{T}()


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

discrete_normalization(p::ModelPlatform, n, L, S) = quadraturenormalization(S)


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
    ExtensionFrame(p.domain, dictionary(p.basisplatform, n))

function discrete_normalization(p::FourierExtensionPlatform, n, L; S = samplingoperator(p, n, L))
    grid1 = grid(dest(S))
    grid2 = supergrid(grid1)
    val = scalar(quadraturenormalization(coefficienttype(dest(S)), grid2))
    ScalingOperator(dest(S), val)
end



struct ExtensionFramePlatform <: FramePlatform
    basisplatform   ::  Platform
    domain          ::  Domain
end

dictionary(p::ExtensionFramePlatform, n) =
    ExtensionFrame(p.domain, dictionary(p.basisplatform, n))

function discrete_normalization(p::ExtensionFramePlatform, n, L; S = samplingoperator(p, n, L))
    grid1 = grid(dest(S))
    grid2 = supergrid(grid1)
    # TODO: make this more generic
    # For now, assume a ScalingOperator from which we can extract the scalar value
    val = scalar(quadraturenormalization(coefficienttype(dest(S)), grid2))
    ScalingOperator(dest(S), val)
end
