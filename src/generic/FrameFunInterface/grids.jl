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
sampling_grid(samplingstyle::GridStyle, ap::ApproximationProblem; options...) =
    haskey(options,:grid) ? options[:grid] : platform_grid(samplingstyle, platform(ap), parameter(ap); options...)
platform_grid(samplingstyle::GridStyle, ap::ApproximationProblem, args...; options...) =
    haskey(options,:grid) ? options[:grid] : platform_grid(samplingstyle, platform(ap), parameter(ap); options...)
platform_grid(samplingstyle::SamplingStyle, ap::ApproximationProblem, args...; options...) =
    haskey(options,:grid) ? options[:grid] : platform_grid(samplingstyle, platform(ap), parameter(ap); options...)
sampling_grid(ss::DiscreteGramStyle, ap::ApproximationProblem; options...) =
    oversampling_grid(ss, platform(ap), parameter(ap), samplingparameter(ss,ap;options...); options...)


sampling_grid(ss::OversamplingStyle, ap::ApproximationProblem; options...) =
    oversampling_grid(ss, platform(ap) ,parameter(ap), samplingparameter(ss,ap;options...);options...)

sampling_grid(samplingstyle::GridStyle, platform::Platform, param, args...; options...) =
    haskey(options,:grid) ?  options[:grid] : platform_grid(samplingstyle, platform(ap), parameter(ap); options...)
sampling_grid(samplingstyle::OversamplingStyle, platform::Platform, param, L, args...; options...) =
    oversampling_grid(samplingstyle, platform, param, L; options...)
