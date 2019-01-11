## Platforms for certain bases and frames

struct FourierPlatform{T} <: BasisPlatform
end

FourierPlatform() = FourierPlatform{Float64}()

dictionary(p::FourierPlatform{T}, n) where {T} = Fourier{T}(n)

dualdictionary(p::FourierPlatform, n; dict = dictionary(p, n)) = dict

SolverStyle(p::FourierPlatform, ::OversamplingStyle) = DualStyle()

dualsamplingoperator(p::FourierPlatform, n, m; S = samplingoperator(p, n; M=m)) =
    quadraturenormalization(S) * S

functionspace(p::FourierPlatform{T}) where {T} =
    L2Space{Complex{T}}(zero(T)..one(T))


struct ChebyshevPlatform{T} <: BasisPlatform
end

ChebyshevPlatform() = ChebyshevPlatform{Float64}()

dictionary(p::ChebyshevPlatform{T}, n) where {T} = ChebyshevT{T}(n)

function dualdictionary(p::ChebyshevPlatform{T}, n; dict = dictionary(p, n)) where {T}
    scaling = ScalingOperator(dict, 2/convert(T, pi)) *
        BasisFunctions.CoefficientScalingOperator(dict, 1, one(T)/2)
    scaling * dict
end

dualsamplingoperator(p::ChebyshevPlatform, n, m; S = samplingoperator(p, n; M=m)) =
    quadraturenormalization(S) * S




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

dualdictionary(p::FourierExtensionPlatform, n; dict = dictionary(p, n)) = dict

function dualsamplingoperator(p::FourierExtensionPlatform, n, m; S = samplingoperator(p, n; M=m))
    grid1 = grid(dest(S))
    grid2 = supergrid(grid1)
    val = scalar(quadraturenormalization(coefficienttype(dest(S)), grid2))
    val * S
end



struct ExtensionFramePlatform <: FramePlatform
    basisplatform   ::  Platform
    domain          ::  Domain
end

dictionary(p::ExtensionFramePlatform, n) =
    ExtensionFrame(p.domain, dictionary(p.basisplatform, n))

dualdictionary(p::ExtensionFramePlatform, n; dict = dictionary(p, n)) =
    dualdictionary(p.basisplatform, n; dict=dict)

function dualsamplingoperator(p::ExtensionFramePlatform, n, m; S = samplingoperator(p, n; M=m))
    grid1 = grid(dest(S))
    grid2 = supergrid(grid1)
    # TODO: make this more generic
    # We can not just apply the quadrature normalization of the grid of S to S,
    # we need the normalization of the underlying basis first, and then we need to
    # restrict that to a subdomain
    # For now, assume a ScalingOperator from which we can extract the scalar value
    val = scalar(quadraturenormalization(coefficienttype(dest(S)), grid2))
    val * S
end
