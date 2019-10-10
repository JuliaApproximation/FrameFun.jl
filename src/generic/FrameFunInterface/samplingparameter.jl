
function samplingparameter(samplingstyle::SamplingStyle, ap::ApproximationProblem; verbose=false, forget=false, options...)
    if haskey(options, :L)
        # The user specified L as an option
        L = options[:L]
        setsamplingparam!(ap, L)
    elseif samplingparam(ap) != nothing && !forget
        # Has it been computed before
        L = samplingparam(ap)
    else
        # It wasn't. We deduce its value from the options given and store the outcome.
        L = deduce_samplingparameter(samplingstyle, ap; verbose=verbose, options...)
        setsamplingparam!(ap, L)
    end
    verbose && println("Sampling parameter: using L = $L")
    return L
end

samplingparameter(samplingstyle::DiscreteGramStyle, ap::ApproximationProblem; options...) =
    samplingparameter(OversamplingStyle(), ap; options...)

deduce_samplingparameter(ss::SamplingStyle, ap::ApproximationProblem; options...) =
    deduce_samplingparameter(ss, platform(ap), parameter(ap); options...)

deduce_samplingparameter(ss::ProductSamplingStyle, ap::ApproximationProblem; options...) =
    map((x,style)->samplingparameter(style, x; options...), ApproximationProblems.unsafe_elements(ap), ss.styles)

function deduce_samplingparameter(ss::ProductSamplingStyle{NTuple{N,OversamplingStyle}}, ap::ApproximationProblem;
        dict=dictionary(ap), oversamplingfactor=2, verbose=false, options...) where N
    oversamplingfactor = oversamplingfactor^(1/N)
    if haskey(options, :M)
        M = options[:M]
        if M isa Integer
            M = round.(Int,collect(size(dict))/length(dict)*M)
        end
    else
        M = round.(Int,oversamplingfactor * collect(size(dict)))
    end
    L = tuple(map((api, Mi)->deduce_samplingparameter(OversamplingStyle(), api;
        oversamplingfactor=oversamplingfactor, M=Mi, verbose=verbose, options...), ApproximationProblems.unsafe_elements(ap), M)...)
    verbose && println("Sampling parameter: best match for M = $M is L = $L")
    L
end

function deduce_oversampling_parameter(ss::OversamplingStyle, args...; dict=dictionary(args...), verbose = false, oversamplingfactor = 2, options...)
    # In the absence of L, we deduce M and then find the best matching L
    # M is either supplied, or we compute it based on the (default) oversamplingfactor
    M = haskey(options, :M) ? options[:M] : round(Int, oversamplingfactor * length(dict))
    verbose && println("Sampling parameter: oversamplingfactor=$oversamplingfactor, length=$(length(dict)), M=$M")
    L = match_and_correct_sampling_parameter(args..., M; samplingstyle=ss, options...)
    verbose && println("Sampling parameter: best match for M = $M is L = $L")
    L
end

deduce_samplingparameter(ss::SamplingStyle, platform::Platform, param; dict = dictionary(platform, param), options...) =
    deduce_samplingparameter(ss, dict; options...)
deduce_samplingparameter(ss::OversamplingStyle, dict::Dictionary; options...) =
    deduce_oversampling_parameter(ss, dict; dict=dict, options...)
deduce_samplingparameter(ss::OversamplingStyle, platform::Platform, param; dict=dictionary(platform, param), options...) =
    deduce_oversampling_parameter(ss, platform, param; dict=dict, options...)
deduce_samplingparameter(ss::InterpolationStyle, dict::Dictionary; options...) =
    match_and_correct_sampling_parameter(dict, length(dict); samplingstyle=ss, options...)
# TODO: implement this one better (more general)
deduce_samplingparameter(::GramStyle, dict::Dictionary; options...) = length(dict)

deduce_samplingparameter(::RectangularGramStyle, dict::Dictionary; projectiondict, options...) = length(projectiondict)

deduce_samplingparameter(::GridStyle, dict::Dictionary; options...) = param(dict)
