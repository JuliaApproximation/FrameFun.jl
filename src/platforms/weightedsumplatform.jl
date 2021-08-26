
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

unsafe_dictionary(platform::WeightedSumPlatform, i::NTuple) =
    MultiDict([weight(platform, j) * unsafe_dictionary(platform.P, i[j]) for j in 1:length(platform.weights)])

WeightedSumPlatform(platform::Platform, weights::Function...) where N = WeightedSumPlatform(platform, weights)

SolverStyle(platform::WeightedSumPlatform, ::SamplingStyle) = AZStyle()

weight(platform::WeightedSumPlatform, i) = platform.weights[i]

# dictionary(platform::WeightedSumPlatform, i) =
#     MultiDict([weight(platform, j) * dictionary(platform.P, i) for j in 1:length(platform.weights)])



function dualdictionary(platform::WeightedSumPlatform, param::NTuple, measure::Measure; options...)
    denom = (x...)->sum(map(w->abs(w(x...))^2, platform.weights))
    # TODO: discuss, what is the relation between param, L of a platform and platform.P
    MultiDict([((x...)->(platform.weights[j](x...)/denom(x...))) * dualdictionary(platform.P, param[j], measure; options...) for j=1:length(platform.weights)])
end

measure(platform::WeightedSumPlatform) = measure(platform.P)

param_first(platform::WeightedSumPlatform) = ntuple(k->param_first(platform.P),Val(length(platform.weights)))
param_double(platform::WeightedSumPlatform, param::NTuple) =
    ntuple(k->param_double(platform.P, param[k]), Val(length(platform.weights)))
param_inbetween(platform::WeightedSumPlatform, param1::NTuple, param2::NTuple) =
    ntuple(k->param_inbetween(platform.P, param1[k], param2[k]), Val(length(platform.weights)))

function default_param_path(platform::WeightedSumPlatform)
    firstpath = default_param_path(platform.P)
    secondpath = NTupleParameterPath{length(platform.weights)}(first(firstpath))
    HierarchyPath(firstpath, secondpath)
end
