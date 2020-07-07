
samplingparameter(samplingstyle::SamplingStyle, ap::ApproximationProblem; samplingstrategy=DefaultSamplingStrategy(), options...)=
    samplingparameter(samplingstyle, samplingstrategy, ap; options...)

samplingparameter(samplingstyle::GenericSamplingStyle, ap::ApproximationProblem; options...) = -1

function samplingparameter(samplingstyle::SamplingStyle, samplingstrategy::SamplingStrategy, ap::ApproximationProblem;
        verbose=false, forget=false, options...)
    if haskey(options, :L)
        # The user specified L as an option
        L = options[:L]
        setsamplingparam!(ap, L)
    elseif samplingparam(ap) != nothing && !forget
        # Has it been computed before
        L = samplingparam(ap)
    else
        # It wasn't. We deduce its value from the options given and store the outcome.
        L = deduce_samplingparameter(samplingstyle, samplingstrategy, ap; verbose=verbose, options...)
        setsamplingparam!(ap, L)
    end
    verbose && println("Sampling parameter: using L = $L")
    return L
end

samplingparameter(samplingstyle::DiscreteGramStyle, samplingstrategy::SamplingStrategy, ap::ApproximationProblem; options...) =
    samplingparameter(OversamplingStyle(), samplingstrategy, ap; options...)

deduce_samplingparameter(ss::SamplingStyle, samplingstrategy::SamplingStrategy, ap::ApproximationProblem; options...) =
    deduce_samplingparameter(ss, samplingstrategy, platform(ap), parameter(ap); options...)

deduce_samplingparameter(ss::ProductSamplingStyle, samplingstrategy::SamplingStrategy, ap::ProductPlatformApproximation; options...) =
    map((x,style)->samplingparameter(style, samplingstrategy, x; options...), ApproximationProblems.unsafe_elements(ap), ss.styles)

function deduce_samplingparameter(ss::ProductSamplingStyle{NTuple{N,OversamplingStyle}}, samplingstrategy::SamplingStrategy, ap::ProductPlatformApproximation;
        dict=dictionary(ap), oversamplingfactor=2, verbose=false, options...) where N
    oversamplingfactor = oversamplingfactor_float2tuple(samplingstrategy, oversamplingfactor, dict)
    if haskey(options, :M)
        M = options[:M]
        if M isa Integer
            M = M_int2tuple(samplingstrategy, M, dict)
        end
    else
        M = oversamplingM(samplingstrategy, oversamplingfactor, dict)
    end
    L = tuple(map((api, Mi)->deduce_samplingparameter(OversamplingStyle(), samplingstrategy, api;
        oversamplingfactor=oversamplingfactor, M=Mi, verbose=verbose, options...), ApproximationProblems.unsafe_elements(ap), M)...)
    verbose && println("Sampling parameter: best match for M = $M is L = $L")
    L
end

M_int2tuple(::DefaultSamplingStrategy, M::Int, dict::TensorProductDict) =
    ntuple(k->round(Int,M*size(dict, k)/length(dict)))
oversamplingM(::DefaultSamplingStrategy, oversamplingfactor::Real, dict::TensorProductDict) =
    ntuple(k->round(Int,oversamplingfactor*size(dict, k)),Val(dimension(dict)))
oversamplingfactor_float2tuple(::DefaultSamplingStrategy, oversamplingfactor::Real, dict::TensorProductDict) =
    oversamplingfactor^(1/dimension(dict))

function deduce_oversampling_parameter(ss::OversamplingStyle, samplingstrategy::SamplingStrategy, args...; dict=dictionary(args...), verbose = false, oversamplingfactor = 2, options...)
    # In the absence of L, we deduce M and then find the best matching L
    # M is either supplied, or we compute it based on the (default) oversamplingfactor
    M = haskey(options, :M) ? options[:M] : round(Int, oversamplingfactor * length(dict))
    verbose && println("Sampling parameter: oversamplingfactor=$oversamplingfactor, length=$(length(dict)), M=$M")
    L = match_and_correct_sampling_parameter(samplingstrategy, args..., M; samplingstyle=ss, options...)
    verbose && println("Sampling parameter: best match for M = $M is L = $L")
    L
end

deduce_samplingparameter(ss::SamplingStyle, samplingstrategy::SamplingStrategy, platform::Platform, param; dict = dictionary(platform, param), options...) =
    deduce_samplingparameter(ss, samplingstrategy, dict; options...)
deduce_samplingparameter(ss::OversamplingStyle, samplingstrategy::SamplingStrategy, dict::Dictionary; options...) =
    deduce_oversampling_parameter(ss, samplingstrategy, dict; dict=dict, options...)
deduce_samplingparameter(ss::OversamplingStyle, samplingstrategy::SamplingStrategy, platform::Platform, param; dict=dictionary(platform, param), options...) =
    deduce_oversampling_parameter(ss, samplingstrategy, platform, param; dict=dict, options...)
deduce_samplingparameter(ss::InterpolationStyle, samplingstrategy::SamplingStrategy, dict::Dictionary; options...) =
    match_and_correct_sampling_parameter(samplingstrategy, dict, length(dict); samplingstyle=ss, options...)
# TODO: implement this one better (more general)
deduce_samplingparameter(::GramStyle, samplingstrategy::SamplingStrategy, dict::Dictionary; options...) = length(dict)

deduce_samplingparameter(::RectangularGramStyle, samplingstrategy::SamplingStrategy, dict::Dictionary; projectiondict, options...) = length(projectiondict)

deduce_samplingparameter(::GridStyle, samplingstrategy::SamplingStrategy, dict::Dictionary; options...) = param(dict)
