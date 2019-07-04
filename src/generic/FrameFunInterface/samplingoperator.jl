function samplingoperator(samplingstyle::DiscreteStyle, ap::ApproximationProblem;
            T = coefficienttype(ap), options...)
    dmeasure = discretemeasure(samplingstyle, ap; options...)
    g  = grid(dmeasure)
    w = BasisFunctions.weights(dmeasure)
    GS = GridSampling(g, T)
    D = DiagonalOperator(dest(GS), dest(GS), sqrt.(w))
    D*GS
end

samplingoperator(samplingstyle::ProductSamplingStyle, ap::ApproximationProblem; options...) =
    tensorproduct( elements(samplingoperator, samplingstyle, ap; options...)... )

samplingoperator(samplingstyle::DiscreteGramStyle, ap::ApproximationProblem; options...) =
    ProjectionSampling(dictionary(ap), measure(samplingstyle, ap; options...))

samplingoperator(::GenericSamplingStyle, ap::PlatformApproximation; options...) =
    genericsamplingoperator(ap.platform, ap.param; dict = ap.dict, options...)

samplingoperator(samplingstyle::GramStyle, ap::ApproximationProblem; options...) =
    ProjectionSampling(dictionary(ap), measure(samplingstyle, ap; options...))

dualsamplingoperator(samplingstyle::GramStyle, ap::ApproximationProblem; options...) =
    ProjectionSampling(
        haskey(options, :dualdict) ? options[:dualdict] : azdual_dict(samplingstyle, ap; options...),
            measure(samplingstyle, ap; options...))

dualsamplingoperator(samplingstyle::DiscreteGramStyle, ap::ApproximationProblem; options...) =
    ProjectionSampling(
        haskey(options, :dualdict) ? options[:dualdict] : azdual_dict(samplingstyle, ap; options...),
            measure(samplingstyle, ap; options...))

samplingoperator(samplingstyle::RectangularGramStyle, ap::ApproximationProblem;
            projectiondict, options...) =
    ProjectionSampling(projectiondict, measure(samplingstyle, ap; options...))



# First, if necessary, we recompute the samplingoperator using the given options.
function dualsamplingoperator(samplingstyle::SamplingStyle, ap::ApproximationProblem; options...)
    if haskey(options, :S)
        S = options[:S]
    else
        S = samplingoperator(samplingstyle, ap; options...)
    end
    S
end
