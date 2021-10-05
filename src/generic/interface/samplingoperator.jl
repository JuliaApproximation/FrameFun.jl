
## First we treat weighted sampling

export sampling_weights
@ap_interface sampling_weights
@ap_sampling_property sampling_weights

sampling_weights(sstyle::UnweightedDiscreteStyle, ap; options...) =
    BasisFunctions.uniformweights(sampling_grid(sstyle, ap; options...))

sampling_weights(sstyle::NormalizedSampling, ap; options...) =
    BasisFunctions.normalizing_weights(sampling_grid(sstyle, ap; options...), measure(ap))

sampling_weights(sstyle::WeightedSampling, ap; options...) =
    sampling_weights(sstyle, platform(ap), platformparameter(ap), samplingparameter(ap); options...)


## Now we can make a sampling operator

export sampling_operator, dual_sampling_operator
@ap_interface sampling_operator
@ap_sampling_property sampling_operator

@ap_interface dual_sampling_operator
@ap_sampling_property dual_sampling_operator

sampling_operator(ss::UnweightedDiscreteStyle, ap::ApproximationProblem; options...) =
    GridSampling(sampling_grid(ss, ap; options...), coefficienttype(ap))

@deprecate samplingoperator(args...; options...) sampling_operator(args...; options...)
@deprecate dualsamplingoperator(args...; options...) dual_sampling_operator(args...; options...)

function sampling_operator(ss::WeightedDiscreteStyle, ap::ApproximationProblem; options...)
    grid = sampling_grid(ss, ap; options...)
    weights = sampling_weights(ss, ap; options...)
    T = coefficienttype(ap)
    S = GridSampling(grid, T)
    D = DiagonalOperator(dest(S), dest(S), weights)
    D*S
end

sampling_operator(ss::DiscreteGramStyle, ap::ApproximationProblem; options...) =
    ProjectionSampling(dictionary(ap), discretemeasure(ap))

sampling_operator(ss::GenericSamplingStyle, ap::PlatformApproximation; options...) =
    sampling_operator(ss, platform(ap), platformparameter(ap), samplingparameter(ap); options...)

sampling_operator(ss::GramStyle, ap::ApproximationProblem; options...) =
    ProjectionSampling(dictionary(ap), measure(ap))

sampling_operator(samplingstyle::RectangularGramStyle, ap::ApproximationProblem;
            projectiondict, options...) =
    ProjectionSampling(projectiondict, measure(ap))

# TODO: fix, remove reference to azdual
dual_sampling_operator(samplingstyle::GramStyle, ap::ApproximationProblem; options...) =
    ProjectionSampling(
        haskey(options, :dualdict) ? options[:dualdict] : azdual(samplingstyle, ap; options...),
            measure(ap))

# TODO: fix, remove reference to azdual
dual_sampling_operator(samplingstyle::DiscreteGramStyle, ap::ApproximationProblem; options...) =
    ProjectionSampling(
        haskey(options, :dualdict) ? options[:dualdict] : azdual(samplingstyle, ap; options...),
            measure(ap))

function dual_sampling_operator(ss::SamplingStyle, ap::ApproximationProblem; options...)
    S = sampling_operator(ap; options...)
    todual(S)
end

todual(S::BasisFunctions.GridSampling) = S
todual(S::BasisFunctions.ProjectionSampling) = S
todual(S::BasisFunctions.GenericCompositeOperator) = compose(S.operators[1], map(inv, S.operators[2:end])...)
