module WeightedSumPlatforms
using ..Platforms
import ..Platforms: SolverStyle, dictionary, measure, dualdictionary, param_first,
    param_next, correctparamformat, unsafe_dictionary
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

correctparamformat(platform::WeightedSumPlatform, param::NTuple) =
    all(map(parami->correctparamformat(platform.P, parami), param))

correctparamformat(::WeightedSumPlatform, _) = false

WeightedSumPlatform(platform::Platform, weights::Function...) where N = WeightedSumPlatform(platform, weights)

SolverStyle(platform::WeightedSumPlatform, ::SamplingStyle) = AZStyle()

weight(platform::WeightedSumPlatform, i) = platform.weights[i]

# dictionary(platform::WeightedSumPlatform, i) =
#     MultiDict([weight(platform, j) * dictionary(platform.P, i) for j in 1:length(platform.weights)])

unsafe_dictionary(platform::WeightedSumPlatform, i::NTuple) =
    MultiDict([weight(platform, j) * unsafe_dictionary(platform.P, i[j]) for j in 1:length(platform.weights)])


function dualdictionary(platform::WeightedSumPlatform, param::NTuple, measure::Measure; options...)
    denom = (x...)->sum(map(w->abs(w(x...))^2, platform.weights))
    # TODO: discuss, what is the relation between param, L of a platform and platform.P
    MultiDict([((x...)->(platform.weights[j](x...)/denom(x...))) * dualdictionary(platform.P, param[j], measure; options...) for j=1:length(platform.weights)])
end

measure(platform::WeightedSumPlatform) = measure(platform.P)

param_first(platform::WeightedSumPlatform) = ntuple(k->param_first(platform.P),Val(length(platform.weights)))
param_next(platform::WeightedSumPlatform, param::NTuple) =
    ntuple(k->param_next(platform.P, param[k]), Val(length(platform.weights)))

end
