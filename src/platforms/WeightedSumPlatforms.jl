module WeightedSumPlatforms
using ..Platforms
import ..Platforms: SolverStyle, dictionary, measure, dualdictionary, param_first, param_next
using BasisFunctions: Dictionary, Measure, MultiDict

export WeightedSumPlatform
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


function dualdictionary(platform::WeightedSumPlatform, param, measure::Measure; options...)
    denom = (x...)->sum(map(w->abs(w(x...))^2, platform.weights))
    # TODO: discuss, what is the relation between param, L of a platform and platform.P
    MultiDict([((x...)->(platform.weights[j](x...)/denom(x...))) * dualdictionary(platform.P, param, measure; options...) for j=1:length(platform.weights)])
end

measure(platform::WeightedSumPlatform) = measure(platform.P)

param_first(platform::WeightedSumPlatform) = param_first(platform.P)
param_next(platform::WeightedSumPlatform, param) = param_next(platform.P, param)

end
