
# For product sampling we make a product grid
for f in (:sampling_grid, :interpolation_grid, :platform_grid)
    @eval begin
        $f(ss::ProductSamplingStyle, ap::ApproximationProblem; options...) =
            ProductGrid(ap_components($f, ss, ap; options...)...)
    end
end

sampling_grid(ss::SamplingStyle, ap::ApproximationProblem; options...) =
    sampling_grid(ss, platform(ap), platformparameter(ap), samplingparameter(ss, ap; options...); options...)

sampling_grid(ss::InterpolationStyle, ap::ApproximationProblem; options...) =
    interpolation_grid(ss, platform(ap), platformparameter(ap); dict = ap.dict, options...)

# - generic grid: we invoke platform_grid on the platform
sampling_grid(samplingstyle::GridStyle, ap::ApproximationProblem; options...) =
    haskey(options,:grid) ? options[:grid] : platform_grid(samplingstyle, platform(ap), platformparameter(ap); options...)
platform_grid(samplingstyle::GridStyle, ap::ApproximationProblem, args...; options...) =
    haskey(options,:grid) ? options[:grid] : platform_grid(samplingstyle, platform(ap), platformparameter(ap); options...)
sampling_grid(ss::DiscreteGramStyle, ap::ApproximationProblem; options...) =
    oversampling_grid(ss, platform(ap), platformparameter(ap), samplingparameter(ss,ap;options...); options...)


sampling_grid(ss::OversamplingStyle, ap::ApproximationProblem; options...) =
    oversampling_grid(ss, platform(ap), platformparameter(ap), samplingparameter(ss,ap;options...);options...)

sampling_grid(samplingstyle::GridStyle, platform::Platform, param, args...; options...) =
    haskey(options,:grid) ?  options[:grid] : platform_grid(samplingstyle, platform(ap), platformparameter(ap); options...)
sampling_grid(samplingstyle::OversamplingStyle, platform::Platform, param, L, args...; options...) =
    oversampling_grid(samplingstyle, platform, param, L; options...)
