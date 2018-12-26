
"""
A WeightedSumFrame is the union of a finite number of copies of a single frame,
each weighted by a function.
"""
struct WeightedSumPlatform <: FramePlatform
    P       :: Platform
    weights :: Vector{Function}
end

WeightedSumPlatform(platform::Platform, weights) = WeightedSumPlatform(platform, weights)

SolverStyle(platform::WeightedSumPlatform, ::SamplingStyle) = AZStyle()

weight(platform::WeightedSumPlatform, i) = platform.weights[i]

Dictionary(platform::WeightedSumPlatform, i) =
    MultiDict([weight(platform, j) * Dictionary(platform.P,i) for j in 1:length(platform.weights)])

function dualdictionary(platform::WeightedSumPlatform, i; dict = Dictionary(platform, i))
    denom = (x...)->sum(map(w->abs(w(x...))^2, platform.weights))
    MultiDict([((x...)->(platform.weights[j](x...)/denom(x...))) * dualdictionary(platform.P,i) for j=1:length(platform.weights)])
end

oversampledgrid(p::WeightedSumPlatform, param, dict, M) = oversampledgrid(element(dict,1), M)

dualsamplingoperator(platform::WeightedSumPlatform, n, m; S = samplingoperator(platform, n; M=m)) =
    dualsamplingoperator(platform.P, n, m, S=S)
