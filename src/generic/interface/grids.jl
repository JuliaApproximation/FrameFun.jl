import BasisFunctions: interpolation_grid
export interpolation_grid
@ap_interface interpolation_grid
@ap_sampling_property_to_platform interpolation_grid

interpolation_grid(ss::SamplingStyle, platform::Platform, plt_par, smpl_par; options...) =
    interpolation_grid(dictionary(platform, plt_par))

export oversampling_grid
@ap_interface oversampling_grid
@ap_sampling_property_to_platform oversampling_grid

oversampling_grid(ss::SamplingStyle, platform::Platform, plt_par, smpl_par; options...) =
    oversampling_grid(dictionary(platform, plt_par), smpl_par)
include("oversampling_grid.jl")

export platform_grid
@ap_interface platform_grid
@ap_sampling_property_to_platform platform_grid

# there is no default, signal a hopefully informative error
platform_grid(ss::SamplingStyle, platform::Platform, plt_par, smpl_par; grid = nothing, options...) =
    grid == nothing ? error("This platform does not implement a `platform_grid`.") : grid


export sampling_grid
@ap_interface sampling_grid
@ap_sampling_property sampling_grid

"Return the sampling grid that would be used for an approximation problem."
sampling_grid(ss::InterpolationStyle, ap::ApproximationProblem; options...) =
    interpolation_grid(ap; options...)
sampling_grid(ss::GridStyle, ap::ApproximationProblem; grid = nothing, options...) =
    grid == nothing ? platform_grid(ap; options...) : grid
sampling_grid(ss::OversamplingStyle, ap::ApproximationProblem; options...) =
    oversampling_grid(ap; options...)
sampling_grid(ss::WeightedDiscreteStyle, ap::ApproximationProblem; options...) =
    sampling_grid(unweighted(ss), ap; options...)

to_product(grids::AbstractGrid...) = productgrid(grids...)
# function sampling_grid(ss::ProductSamplingStyle, ap::ApproximationProblem; options...)
#     @assert nfactors(ss)==nfactors(platform(ap))==length(platformparameter(ap))==length(samplingparameter(ap))
#     productgrid( map((st,pl,p_par,s_par)->sampling_grid(pl, p_par, s_par; samplingstyle=st, options...),
#         factors(ss), factors(platform(ap)), productparameter(ap), samplingparameter(ap))...)
# end

# sampling_grid(ss::ProductSamplingStyle, ap::ApproximationProblem; options...) =
#     sampling_grid(ss, platform(ap), productparameter(ap), samplingparameter(ap); options...)
#
# function sampling_grid(ss::ProductSamplingStyle, platform::ProductPlatform, plt_par, smpl_par; options...)
#     @assert nfactors(ss)==nfactors(platform)==length(plt_par)==length(smpl_par)
#     productgrid(map(
#         (pl,pl_par,s_par, st) -> sampling_grid(pl, pl_par, s_par; samplingstyle=st, options...),
#         factors(platform), plt_par, smpl_par, factors(ss))...)
# end
