
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

dictionary(platform::WeightedSumPlatform, i) =
    MultiDict([weight(platform, j) * dictionary(platform.P, i) for j in 1:length(platform.weights)])

function dualplatformdictionary(sstyle::DiscreteStyle, platform::WeightedSumPlatform, i; options...)
    denom = (x...)->sum(map(w->abs(w(x...))^2, platform.weights))
    MultiDict([((x...)->(platform.weights[j](x...)/denom(x...))) * dualplatformdictionary(sstyle, platform.P, i; options...) for j=1:length(platform.weights)])
end

oversampling_grid(p::WeightedSumPlatform, param, L; dict, options...) =
    oversampling_grid(superdict(element(dict,1)), L)

discrete_normalization(platform::WeightedSumPlatform, n, L; options...) =
    discrete_normalization(platform.P, n, L; options...)

measure(platform::WeightedSumPlatform) = measure(platform.P)

nextsize(p::WeightedSumPlatform, n) = 2n
