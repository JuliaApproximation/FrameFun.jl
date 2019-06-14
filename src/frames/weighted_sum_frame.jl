
"""
A WeightedSumPlatform is the union of a finite number of copies of a single frame,
each weighted by a function.
"""
struct WeightedSumPlatform <: FramePlatform
    P       :: Platform
    weights :: NTuple{N,Function} where N
end

WeightedSumPlatform(platform::Platform, weights::Function...) where N = WeightedSumPlatform(platform, weights)

SolverStyle(platform::WeightedSumPlatform, ::SamplingStyle) = AZStyle()

weight(platform::WeightedSumPlatform, i) = platform.weights[i]

dictionary(platform::WeightedSumPlatform, i) =
    MultiDict([weight(platform, j) * dictionary(platform.P, i) for j in 1:length(platform.weights)])
measure(platform::WeightedSumPlatform, param) = measure(platform.P, param)

function azdual_dict(sstyle::SamplingStyle, platform::WeightedSumPlatform, param, L, measure::Measure; options...)
   denom = (x...)->sum(map(w->abs(w(x...))^2, platform.weights))
   # TODO: discuss, what is the relation between param, L of a platform and platform.P
   MultiDict([((x...)->(platform.weights[j](x...)/denom(x...))) * azdual_dict(sstyle, platform.P, param, L, measure; options...) for j=1:length(platform.weights)])
end

function BasisFunctions.gramdual(dict::MultiDict, measure::Measure; options...)
    @debug "Are you sure you want `dualtype=gramdual` and not `weightedsumdual`"
    BasisFunctions.default_gramdual(dict, measure; options...)
end


oversampling_grid(p::WeightedSumPlatform, param, L; dict, options...) =
    oversampling_grid(superdict(element(dict,1)), L)

discrete_normalization(platform::WeightedSumPlatform, param, L; options...) =
    discrete_normalization(platform.P, param, L; options...)

measure(platform::WeightedSumPlatform) = measure(platform.P)

param_next(p::WeightedSumPlatform, param) = 2param
