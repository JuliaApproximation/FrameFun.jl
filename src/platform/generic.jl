
"""
A `GenericPlatform` stores a primal and dual dictionary generator, along with
a sequence of parameter values.
"""
struct GenericPlatform <: Platform
    super_platform
    primal_generator
    dual_generator
    sampler_generator
    dual_sampler_generator
    parameter_sequence
    name
end

GenericPlatform(; super_platform=nothing, primal = None, dual = primal, sampler = None, dual_sampler=sampler, params = None,
    name = "Generic Platform") = GenericPlatform(super_platform, primal, dual, sampler, dual_sampler, params, name)

function primal(platform::GenericPlatform, i)
    param = platform.parameter_sequence[i]
    platform.primal_generator(param)
end

function dual(platform::GenericPlatform, i)
    param = platform.parameter_sequence[i]
    platform.dual_generator(param)
end

function sampler(platform::GenericPlatform, i)
    param = platform.parameter_sequence[i]
    platform.sampler_generator(param)
end

function dual_sampler(platform::GenericPlatform, i)
    param = platform.parameter_sequence[i]
    platform.dual_sampler_generator(param)
end

name(platform::GenericPlatform) = platform.name

matrix_A(platform::Platform, i; options...) = apply(sampler(platform, i), primal(platform, i); options...)

matrix_Zt(platform::Platform, i; options...) = matrix_Zt(dual(platform, i), dual_sampler(platform, i); options...)

matrix_Zt(dual, dual_sampler; options...) = apply(dual_sampler, dual; options...)'


"""
Initalized with a series of generators, it generates tensorproduct dictionaries
given a series of lengths.
"""
struct TensorGenerator{T}
    fun
end
(TG::TensorGenerator)(n::Int...) = TG(collect(n))
(TG::TensorGenerator)(n::AbstractVector{Int}) = tensorproduct(TG.fun(n))

tensor_generator(::Type{T}, generators...) where {T} = TensorGenerator{T}( n ->([gi(ni)  for (ni, gi) in  zip(collect(n), collect(generators))]))



"A doubling sequence that produces odd values, i.e., the value after `n` is `2n+1`."
struct OddDoublingSequence <: DimensionSequence
    initial ::  Int

	# Default constructor to guarantee that the initial value is odd
	OddDoublingSequence(initial::Int) = (@assert isodd(initial); new(initial))
end

initial(s::OddDoublingSequence) = s.initial

OddDoublingSequence() = OddDoublingSequence(1)

getindex(s::OddDoublingSequence, idx::Int) = initial(s) * (2<<(idx-1)) - 1


fourier_platform(;options...) = fourier_platform(Float64; options...)

fourier_platform(n::Int; options...) = fourier_platform(Float64, n; options...)

fourier_platform(::Type{T};options...) where {T} = fourier_platform(T, 1; options...)

function fourier_platform(::Type{T}, n::Int; oversampling=1) where {T}
	primal = Fourier{T}
	dual = Fourier{T}
        sampler = n -> GridSampling(GridBasis{T}(PeriodicEquispacedGrid(round(Int,oversampling*n), UnitInterval{T}())))
        dual_sampler = n->(1/length(dest(sampler(n))))*sampler(n)
	params = isodd(n) ? OddDoublingSequence(n) : DoublingSequence(n)
	GenericPlatform(primal = primal, dual = dual, sampler = sampler, dual_sampler=dual_sampler,
		params = params, name = "Fourier series")
end
