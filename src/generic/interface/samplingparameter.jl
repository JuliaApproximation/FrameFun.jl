
# Functionality for computing a suitable sampling parameter of an approximation problem

function samplingparameter(ss::SamplingStyle, platform::Platform, plt_par;
        verbose=false, L = -1, smpl_par = -1, options...)
    # Was the parameter specified as an option?
    if (L > 0) || (smpl_par > 0)
        smpl_par = max(L, smpl_par)
    else
        # It wasn't. Try to deduce its value from the options given.
        smpl_par = deduce_samplingparameter(ss, platform, plt_par; verbose=verbose, options...)
    end
    verbose && println("Sampling parameter: using smpl_par = $smpl_par")
    return smpl_par
end

samplingparameter(ss::GenericSamplingStyle, platform::Platform, plt_par; options...) = nothing

samplingparameter(ss::WeightedDiscreteStyle, platform::Platform, plt_par; options...) =
    samplingparameter(unweighted(ss), platform, plt_par; options...)

samplingparameter(ss::DiscreteGramStyle, platform::Platform, plt_par; options...) =
    samplingparameter(OversamplingStyle(), platform, plt_par; options...)

deduce_samplingparameter(ss::SamplingStyle, platform::Platform, plt_par; options...) =
    deduce_samplingparameter(ss, dictionary(platform, plt_par); options...)

function deduce_oversampling_parameter(ss::OversamplingStyle, dict; verbose = false, oversamplingfactor = 2, options...)
    # In the absence of smpl_par, we deduce M and then find the best matching smpl_par
    # M is either supplied, or we compute it based on the (default) oversamplingfactor
    M = haskey(options, :M) ? options[:M] : round(Int, oversamplingfactor * length(dict))
    verbose && println("Sampling parameter: oversamplingfactor=$oversamplingfactor, length=$(length(dict)), M=$M")
    smpl_par = match_and_correct_sampling_parameter(dict, M; samplingstyle=ss, options...)
    verbose && println("Sampling parameter: best match for M = $M is smpl_par = $smpl_par")
    smpl_par
end

deduce_samplingparameter(ss::OversamplingStyle, dict::Dictionary; options...) =
    deduce_oversampling_parameter(ss, dict; options...)

deduce_samplingparameter(ss::InterpolationStyle, dict::Dictionary; options...) =
    match_and_correct_sampling_parameter(dict, length(dict); samplingstyle=ss, options...)

# TODO: implement this one better (more general)
deduce_samplingparameter(::GramStyle, dict::Dictionary; options...) = length(dict)

deduce_samplingparameter(::RectangularGramStyle, dict::Dictionary; projectiondict, options...) =
    length(projectiondict)

deduce_samplingparameter(::GridStyle, dict::Dictionary; options...) = platformparameter(dict)


## Product sampling

oversamplingfactor_to_tuple(of::Real, dim::Int) = ntuple(t->of^(1/dim), Val(dim))
oversamplingfactor_to_tuple(of::Tuple, dim::Int) = (@assert length(of)==dim; of)


deduce_samplingparameter(ss::ProductSamplingStyle, platform::ProductPlatform, plt_par; options...) =
    map((style,pl,ppar)->samplingparameter(style, pl, ppar; options...),
        factors(ss), factors(platform), productparameter(platform, plt_par))

function deduce_samplingparameter(ss::ProductSamplingStyle{NTuple{N,OversamplingStyle}},
        platform::ProductPlatform, plt_par;
        oversamplingfactor=2, M = -1, verbose=false, options...) where {N}
    dict = dictionary(platform, plt_par)
    if !(M == -1)
        if M isa Integer
            oversamplingfactor = M / length(dict)
        else
            @assert length(M) == N
            oversamplingfactor = M ./ map(length, factors(dict))
        end
    end
    oversamplingfactor = oversamplingfactor_to_tuple(oversamplingfactor, N)
    smpl_par = tuple(map((pl,ppar,of)->deduce_samplingparameter(OversamplingStyle(), pl, ppar;
        oversamplingfactor=of, verbose=verbose, options...),
            factors(platform), productparameter(platform, plt_par), oversamplingfactor)...)
    verbose && println("Sampling parameter: best match for M = $M is smpl_par = $smpl_par")
    smpl_par
end
