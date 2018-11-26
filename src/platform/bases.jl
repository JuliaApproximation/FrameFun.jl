## Platforms for Fourier series

struct FourierPlatform{T} <: BasisPlatform
end

FourierPlatform() = FourierPlatform{Float64}()

dictionary(p::FourierPlatform{T}, n; options...) where {T} = FourierBasis{T}(n)


struct ChebyshevPlatform{T} <: BasisPlatform
end

ChebyshevPlatform() = ChebyshevPlatform{Float64}()

dictionary(p::ChebyshevPlatform{T}, n; options...) where {T} = ChebyshevBasis{T}(n)

function interpolation_grid(dict::ChebyshevBasis; secondkind = false, options...)
    if secondkind
        secondgrid(dict)
    else
        grid(dict)
    end
end
