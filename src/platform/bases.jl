## Platforms for certain bases

struct FourierPlatform{T} <: BasisPlatform
end

FourierPlatform() = FourierPlatform{Float64}()

Dictionary(p::FourierPlatform{T}, n; options...) where {T} = FourierBasis{T}(n)

DiscretizationStyle(p::FourierPlatform) = InterpolationStyle()
SolverStyle(p::FourierPlatform, ::DiscretizationStyle) = TransformStyle()


struct ChebyshevPlatform{T} <: BasisPlatform
end

ChebyshevPlatform() = ChebyshevPlatform{Float64}()

Dictionary(p::ChebyshevPlatform{T}, n; options...) where {T} = ChebyshevBasis{T}(n)

DiscretizationStyle(p::ChebyshevPlatform) = InterpolationStyle()
SolverStyle(p::ChebyshevPlatform, ::DiscretizationStyle) = TransformStyle()

function interpolation_grid(dict::ChebyshevBasis; secondkind = false, options...)
    if secondkind
        secondgrid(dict)
    else
        grid(dict)
    end
end
