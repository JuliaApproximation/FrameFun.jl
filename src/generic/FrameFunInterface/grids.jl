sampling_grid(::SamplingStyleSuper{InterpolationStyle}, ap::ApproximationProblem; options...) =
    interpolation_grid(platform(ap), parameter(ap); dict = ap.dict, options...)

# - generic grid: we invoke platform_grid on the platform
sampling_grid(ss::SamplingStyleSuper{GridStyle}, ap::ApproximationProblem; options...) =
    platform_grid(ss, ap; options...)

sampling_grid(ss::SamplingStyleSuper{OversamplingStyle}, ap::ApproximationProblem; options...) =
    oversampling_grid(ss, ap; options...)

sampling_grid(ss::SamplingStyleSuper{DiscreteGramStyle}, ap::ApproximationProblem; options...) =
    oversampling_grid(ss, ap; options...)

platform_grid(samplingstyle::SamplingStyle, ap::ApproximationProblem; options...) =
    haskey(options, :grid) ? options[:grid] : error("optional parameter `grid` should be provided using `GridStyle` and platform_grid")

platform_grid(ss::ProductSamplingStyle, ap::ApproximationProblem; options...) =
    ProductGrid(elements(platform_grid, ss, ap; options...)...)
