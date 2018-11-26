
"""
A WeightedSumFrame is the union of a finite number of copies of a single frame,
each weighted by function.
"""
struct WeightedSumPlatform <: FramePlatform
    weights :: Vector{Function}
    P       :: Platform
end

weight(P::WeightedSumPlatform, i) = weights[i]

primal(P::WeightedSumPlatform, i) = MultiDict([P.weights[j]*primal(P.P,i) for j=1:length(P.weights)])

function dual(P::WeightedSumPlatform, i)
    denom = (x...)->sum(map(w->abs(w(x...))^2,P.weights))
    MultiDict([((x...)->(P.weights[j](x...)/denom(x...)))*dual(P.P,i) for j=1:length(P.weights)])
end
sampler(P::WeightedSumPlatform,i) = sampler(P.P,i)
dual_sampler(P::WeightedSumPlatform,i) = dual_sampler(P.P,i)
