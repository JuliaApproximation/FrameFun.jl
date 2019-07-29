for f in (:sampling_grid, :interpolation_grid, :platform_grid)
    @eval begin
        $f(ss::ProductSamplingStyle, ap::ApproximationProblem; options...) =
            ProductGrid(elements($f, ss, ap; options...)...)
    end
end

sampling_grid(ss::SamplingStyle, ap::ApproximationProblem; options...) =
    sampling_grid(ss, platform(ap), parameter(ap), samplingparameter(ss, ap; options...); options...)

sampling_grid(ss::InterpolationStyle, ap::ApproximationProblem; options...) =
    interpolation_grid(ss, platform(ap), parameter(ap); dict = ap.dict, options...)

# - generic grid: we invoke platform_grid on the platform
sampling_grid(ss::GridStyle, ap::ApproximationProblem; options...) =
    platform_grid(ss, platform(ap), parameter(ap); options...)

sampling_grid(ss::DiscreteGramStyle, ap::ApproximationProblem; options...) =
    oversampling_grid(ss, platform(ap), parameter(ap), samplingparameter(ss,ap;options...); options...)

sampling_grid(ss::OversamplingStyle, ap::ApproximationProblem; options...) =
    oversampling_grid(ss, platform(ap) ,parameter(ap), samplingparameter(ss,ap;options...);options...)

platform_grid(samplingstyle::SamplingStyle, platform::Platform, param; options...) =
    haskey(options, :grid) ? options[:grid] : error("optional parameter `grid` should be provided using `GridStyle` and platform_grid")
