##########################
# Approximation problems
##########################

## Types
#  Approximation problem types simply group the arguments a user supplies
#  to the Fun constructor.
#  The implement an interface, but all functionality is delegated to the
#  underlying dictionary or platform.

abstract type ApproximationProblem end

struct DictionaryApproximation <: ApproximationProblem
    dict    ::  Dictionary
end

dictionary(ap::DictionaryApproximation) = ap.dict

coefficienttype(ap::DictionaryApproximation) = coefficienttype(ap.dict)

struct PlatformApproximation <: ApproximationProblem
    platform    ::  Platform
    param
    dict        ::  Dictionary

    PlatformApproximation(platform, param) = new(platform, param, dictionary(platform, param))
end

dictionary(ap::PlatformApproximation) = ap.dict

coefficienttype(ap::PlatformApproximation) = coefficienttype(ap.dict)


struct AdaptiveApproximation <: ApproximationProblem
    platform    ::  Platform
end

dictionary(ap::AdaptiveApproximation, param) = ap.platform[param]



## Construct an approximation problem
# - from a dictionary
approximationproblem(dict::Dictionary) = approximationproblem(coefficienttype(dict), dict)
approximationproblem(dict::Dictionary, domain::Domain) =
    approximationproblem(promote_type(coefficienttype(dict),eltype(domain)), dict, domain)

approximationproblem(::Type{T}, dict::Dictionary) where {T} =
    DictionaryApproximation(promote_coefficienttype(dict, T))

function approximationproblem(::Type{T}, dict::Dictionary, domain::Domain) where {T}
    if domain == support(dict)
        approximationproblem(T, dict)
    else
        approximationproblem(T, ExtensionFrame(domain, promote_coefficienttype(dict, T)))
    end
end

# - from a platform
approximationproblem(platform, param) = PlatformApproximation(platform, param)
approximationproblem(platform) = AdaptiveApproximation(platform)

approximationproblem(ap::AdaptiveApproximation, param) = approximationproblem(ap.platform, param)


########################
# Some helper functions
########################

dictionarylength(ap::DictionaryApproximation) = length(ap.dict)
dictionarylength(ap::PlatformApproximation) = length(ap.dict)


########################
# Interface delegation
########################

## SamplingStyle and SolverStyle

SamplingStyle(ap::DictionaryApproximation) = SamplingStyle(dictionary(ap))
SolverStyle(ap::DictionaryApproximation, dstyle::SamplingStyle) = SolverStyle(dictionary(ap), dstyle)

SamplingStyle(ap::PlatformApproximation) = SamplingStyle(ap.platform)
SolverStyle(ap::PlatformApproximation, dstyle::SamplingStyle) = SolverStyle(ap.platform, dstyle)

SamplingStyle(ap::AdaptiveApproximation) = SamplingStyle(ap.platform)
SolverStyle(ap::AdaptiveApproximation, dstyle::SamplingStyle) = SolverStyle(ap.platform, dstyle)

## Sampling operator

samplingoperator(dict::Dictionary, args...; options...) = samplingoperator(approximationproblem(dict, args...); options...)
samplingoperator(platform::Platform, args...; options...) = samplingoperator(approximationproblem(platform, args...); options...)

function samplingoperator(ap::ApproximationProblem;
        M = nothing,
        samplingstyle = SamplingStyle(ap, M),
        options...)
    if M == nothing
        samplingoperator(samplingstyle, ap; options...)
    else
        samplingoperator(samplingstyle, ap; M=M, options...)
    end
end

samplinglength(S::AbstractOperator) = length(dest(S))

oversamplingoperator(platform::Platform, n, m) = samplingoperator(platform, n;
    samplingstyle=OversamplingStyle(), M=m)

function samplingoperator(samplingstyle::DiscreteStyle, ap::ApproximationProblem; T = coefficienttype(ap), options...)
    grid = sampling_grid(samplingstyle, ap; options...)
    GridSampling(grid, T)
end

samplingoperator(::GenericSamplingStyle, ap::PlatformApproximation; options...) =
    samplingoperator(ap.platform, ap.param; dict = ap.dict, options...)

samplingoperator(samplingstyle::GramStyle, ap::ApproximationProblem;
        measure = measure(dictionary(ap)), options...) =
    ProjectionSampling(dictionary(ap), measure)


## Discrete sampling grid

# - interpolation: we ask the dictionary
sampling_grid(::InterpolationStyle, ap; options...) =
    interpolation_grid(dictionary(ap))
# - generic grid: it should be passed as an argument
sampling_grid(::GridStyle, ap; grid, options...) = grid

# - oversampling: we compute an oversampled grid. Its size can optionally be given
# using the M keyword arguments. Here, we compute sensible defaults based on
# oversampling by a factor of 2
oversamplingsize(ap::ApproximationProblem) = oversamplingsize(dictionary(ap))
oversamplingsize(dict::Dictionary) = 2length(dictionary)

sampling_grid(::OversamplingStyle, ap; options...) = oversampledgrid(ap; options...)

oversampledsize(dict::Dictionary, factor) = round(Int, factor*length(dict))

oversampledgrid(ap::DictionaryApproximation;
        oversamplingfactor = 2,
        M = oversampledsize(dictionary(ap), oversamplingfactor), options...) =
    oversampledgrid(dictionary(ap), M)

oversampledsize(ap::PlatformApproximation, factor) = oversampledsize(ap.platform, ap.param, ap.dict, factor)
oversampledsize(p::Platform, param, dict, factor) = oversampledsize(dict, factor)

oversampledgrid(ap::PlatformApproximation;
        oversamplingfactor = 2,
        M = oversampledsize(ap, oversamplingfactor), options...) =
    oversampledgrid(ap.platform, ap.param, ap.dict, M)

oversampledgrid(platform::Platform, param, dict, M) =
    oversampledgrid(dict, M)

##  Solver

for op in (:solver, :AZ_Zt)
    @eval $op(dict::Dictionary, args...; options...) = $op(approximationproblem(dict, args...); options...)
    @eval $op(platform::Platform, args...; options...) = $op(approximationproblem(platform, args...); options...)
end

function solver(ap::ApproximationProblem;
            samplingstyle = SamplingStyle(ap),
            solverstyle = SolverStyle(ap, samplingstyle), options...)
    A, S = discretization(ap; options...)
    solver(solverstyle, ap, A; S = S, options...)
end

AZ_Zt(ap::ApproximationProblem; options...) = AZ_Zt(ap, samplingoperator(ap; options...))


## Dual solvers

dualdictionary(ap::PlatformApproximation) =
    dualdictionary(ap.platform, ap.param; dict = ap.dict)

dualsamplingoperator(ap::PlatformApproximation, S; options...) =
    dualsamplingoperator(ap.platform, ap.param, samplinglength(S); S=S, options...)

dualdictionary(ap::DictionaryApproximation) = dualdictionary(dictionary(ap))

dualsamplingoperator(ap::DictionaryApproximation, S) =
    dualsamplingoperator(dictionary(ap), S)


# # TODO: clean up these scaling factors. They have to do with a normalization of the
# # sampling operators.
Zt_scaling_factor(S::Dictionary, A) = length(supergrid(grid(dest(A))))
Zt_scaling_factor(S::DerivedDict, A) = Zt_scaling_factor(superdict(S), A)
Zt_scaling_factor(S::ChebyshevT, A) = length(supergrid(grid(dest(A))))/2


function AZ_Zt(ap::ApproximationProblem, S)
    dtilde = dualdictionary(ap)
    Stilde = dualsamplingoperator(ap, S)
    apply(Stilde, dtilde)'
end


# The discretization and solver

for op in (:discretization, :dualdiscretization)
    @eval $op(dict::Dictionary, args...; options...) = $op(approximationproblem(dict, args...); options...)
    @eval $op(platform::Platform, args...; options...) = $op(approximationproblem(platform, args...); options...)
end

function discretization(ap::ApproximationProblem; options...)
    d = dictionary(ap)
    S = samplingoperator(ap; options...)
    A = apply(S, d)
    A, S
end

function dualdiscretization(ap::ApproximationProblem; options...)
    A, S = discretization(ap; options...)
    dtilde = dualdictionary(ap)
    Stilde = dualsamplingoperator(ap, S)
    Z = apply(Stilde, dtilde)
    A, Z, S, Stilde
end


for op in (:plungeoperator, :smoothingoperator)
    @eval $op(dict::Dictionary, args...; options...) = $op(approximationproblem(dict, args...); options...)
    @eval $op(platform::Platform, args...; options...) = $op(approximationproblem(platform, args...); options...)
end

function plungeoperator(ap::ApproximationProblem; options...)
    A,Z,S,Stilde = dualdiscretization(ap; options...)
    I = IdentityOperator(dest(A))
    I - A * Z'
end

smoothingoperator(ap::ApproximationProblem; options...) =
    WeightedSmoothingOperator(dictionary(ap); options...)
