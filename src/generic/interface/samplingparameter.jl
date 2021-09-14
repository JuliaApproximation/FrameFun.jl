
samplingparameter(samplingstyle::GenericSamplingStyle, ap::ApproximationProblem; options...) = -1

function samplingparameter(samplingstyle::SamplingStyle, ap::ApproximationProblem;
        verbose=false, forget=false, options...)
    # Was it stored?
    if _samplingparameter(ap) != nothing && !forget
        smpl_par = _samplingparameter(ap)
    else
        # It wasn't stored. Perhaps the user specified the parameter as an option?
        if haskey(options, :L)
            smpl_par = options[:L]
        elseif haskey(options, :smpl_par)
            smpl_par = options[:smpl_par]
        else
            # The user didn't. Try to deduce its value from the options given.
            smpl_par = deduce_samplingparameter(samplingstyle, ap; verbose=verbose, options...)
        end
        # store the outcome for later
        _samplingparameter!(ap, smpl_par)
    end
    verbose && println("Sampling parameter: using smpl_par = $smpl_par")
    return smpl_par
end

samplingparameter(samplingstyle::DiscreteGramStyle, ap::ApproximationProblem; options...) =
    samplingparameter(OversamplingStyle(), ap; options...)

deduce_samplingparameter(ss::SamplingStyle, ap::ApproximationProblem; options...) =
    deduce_samplingparameter(ss, platform(ap), platformparameter(ap); options...)

deduce_samplingparameter(ss::ProductSamplingStyle, ap::ProductPlatformApproximation; options...) =
    map((x,style)->samplingparameter(style, x; options...), FrameFun.unsafe_components(ap), ss.styles)

function deduce_oversampling_parameter(ss::OversamplingStyle, args...; dict=dictionary(args...), verbose = false, oversamplingfactor = 2, options...)
    # In the absence of smpl_par, we deduce M and then find the best matching smpl_par
    # M is either supplied, or we compute it based on the (default) oversamplingfactor
    M = haskey(options, :M) ? options[:M] : round(Int, oversamplingfactor * length(dict))
    verbose && println("Sampling parameter: oversamplingfactor=$oversamplingfactor, length=$(length(dict)), M=$M")
    smpl_par = match_and_correct_sampling_parameter(args..., M; samplingstyle=ss, options...)
    verbose && println("Sampling parameter: best match for M = $M is smpl_par = $smpl_par")
    smpl_par
end

deduce_samplingparameter(ss::SamplingStyle, platform::Platform, plt_par; dict = dictionary(platform, plt_par), options...) =
    deduce_samplingparameter(ss, dict; options...)
deduce_samplingparameter(ss::OversamplingStyle, dict::Dictionary; options...) =
    deduce_oversampling_parameter(ss, dict; dict=dict, options...)
deduce_samplingparameter(ss::OversamplingStyle, platform::Platform, param; dict=dictionary(platform, param), options...) =
    deduce_oversampling_parameter(ss, platform, param; dict=dict, options...)
deduce_samplingparameter(ss::InterpolationStyle, dict::Dictionary; options...) =
    match_and_correct_sampling_parameter(dict, length(dict); samplingstyle=ss, options...)

# TODO: implement this one better (more general)
deduce_samplingparameter(::GramStyle, dict::Dictionary; options...) = length(dict)

deduce_samplingparameter(::RectangularGramStyle, dict::Dictionary; projectiondict, options...) =
    length(projectiondict)

deduce_samplingparameter(::GridStyle, dict::Dictionary; options...) = platformparameter(dict)


M_int2tuple(M::Int, dict::TensorProductDict) =
    ntuple(k->round(Int,M*size(dict, k)/length(dict)))
oversamplingM(oversamplingfactor::Real, dict::TensorProductDict) =
    ntuple(k->round(Int,oversamplingfactor*size(dict, k)),Val(dimension(dict)))
oversamplingfactor_float2tuple(oversamplingfactor::Real, dict::TensorProductDict) =
    oversamplingfactor^(1/dimension(dict))

function deduce_samplingparameter(ss::ProductSamplingStyle{NTuple{N,OversamplingStyle}}, ap::ProductPlatformApproximation;
        dict=dictionary(ap), oversamplingfactor=2, verbose=false, options...) where N
    oversamplingfactor = oversamplingfactor_float2tuple(oversamplingfactor, dict)
    if haskey(options, :M)
        M = options[:M]
        if M isa Integer
            M = M_int2tuple(M, dict)
        end
    else
        M = oversamplingM(oversamplingfactor, dict)
    end
    smpl_par = tuple(map((api, Mi)->deduce_samplingparameter(OversamplingStyle(), api;
        oversamplingfactor=oversamplingfactor, M=Mi, verbose=verbose, options...), FrameFun.unsafe_components(ap), M)...)
    verbose && println("Sampling parameter: best match for M = $M is smpl_par = $smpl_par")
    smpl_par
end
