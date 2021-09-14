
import BasisFunctions: discretemeasure, measure, interpolation_grid

# INTERFACEFUNCTIONS = (:samplingoperator, :dualsamplingoperator, :samplingparameter,
#     :sampling_grid, :interpolation_grid, :platform_grid, :discretization, :dualdiscretization,
#     :azdual_dict, :measure, :discretemeasure, :solver, :AZ_A, :AZ_Z, :AZ_Zt, :plungeoperator,
#     :plungematrix, :plungerank, :smoothingoperator)


"""
    @approximation_interface `function`

Generate an interface for the given function in terms of an approximation
problem and a sampling style. These are generated from the arguments of the
function.
"""
macro approximation_interface(ex)
    def = Meta.parse("default_"*string(ex))
    ret = quote
        # generate an approximation problem
        $(ex)(dict::Dictionary, args...; options...) = $(ex)(approximationproblem(dict, args...); options...)
        $(ex)(platform::Platform, args...; options...) = $(ex)(approximationproblem(platform, args...); options...)

        # add a sampling style
        function $(ex)(ap::ApproximationProblem; samplingstyle=SamplingStyle(ap), options...)
            # compute the sampling parameter with the given options to make sure it is cached
            samplingparameter(samplingstyle, ap; options...)
            $(ex)(samplingstyle, ap; options...)
        end

        # translate the approximation problem to a dictionary via a platform
        # these calls can be intercepted by the platform/dictionary to generate custom behaviour
        $(ex)(ss::SamplingStyle, ap::ApproximationProblem, args...; options...) =
            $(ex)(ss, platform(ap), platformparameter(ap), args...; options...)
        $(ex)(ss::SamplingStyle, platform::Platform, plt_par, args...; dict = dictionary(platform, plt_par), options...) =
            $(ex)(ss, dict, args...; options...)
        $(ex)(ss::SamplingStyle, platform::ParametrizedPlatform, plt_par, args...; options...) =
            $(ex)(ss, platform.platform, platform.path[plt_par], args...; options...)
        $(ex)(ss::SamplingStyle, dict::Dictionary, args...; options...) =
            $(def)(ss, dict, args...; options...)
    end
    esc(ret)
end

"""
    @ap_interface `function`

Generate an interface for the given function in terms of an approximation
problem. The approximation problem is constructed from the arguments of the
function.
"""
macro ap_interface(ex)
    ret = quote
        # generate an approximation problem
        $(ex)(dict::Dictionary, args...; options...) = $(ex)(approximationproblem(dict, args...); options...)
        $(ex)(platform::Platform, args...; options...) = $(ex)(approximationproblem(platform, args...); options...)
    end
    esc(ret)
end

# macro add_samplingparameter_interface(ex)
#     intermediate = Meta.parse("_"*string(ex))
#     ret = quote
#         $(ex)(ss::SamplingStyle, ap::ApproximationProblem; options...) =
#             $(ex)(ss, platform(ap), parameter(ap), samplingparameter(ss, ap; options...); options...)
#     end
#     esc(ret)
# end

macro add_samplingstyle_interface(ex)
    ret = quote
        function $(ex)(ap::ApproximationProblem; samplingstyle = SamplingStyle(ap), options...)
            # compute the sampling parameter once to make sure it is cached
            samplingparameter(samplingstyle, ap; options...)
            $(ex)(samplingstyle, ap; options...)
        end
    end
    esc(ret)
end


# Utility functions to create tensor FrameFunInterface functions

function ap_components(fun, ss::ProductSamplingStyle, ap::ApproximationProblem, args...; options...)
    samplingparameter(ss, ap; options...)
    @assert length(components(ss))==length(productcomponents(ap))
    map((ssi, api, argi)->fun(ssi, api, argi...; options...), components(ss), productcomponents(ap), zip(args...))
end

function ap_components(fun, ss::ProductSamplingStyle, ap::ApproximationProblem; options...)
    samplingparameter(ss, ap; options...)
    @assert length(components(ss))==length(productcomponents(ap))
    map((ssi, api)->fun(ssi, api; options...), components(ss), productcomponents(ap))
end


default_aznormalization(a...) = false





export interpolation_grid
@approximation_interface interpolation_grid
default_interpolation_grid(ss::InterpolationStyle, dict::Dictionary; options...) =
    interpolation_grid(dict)
function default_interpolation_grid(ss::SamplingStyle, dict::Dictionary; oversamplingfactor=1, options...)
    oversamplingfactor != 1 && @warn "interpolation does not support option `oversamplingfactor=$oversamplingfactor`, use `OversamplingStyle()` instead."
    interpolation_grid(dict)
end


export oversampling_grid
@approximation_interface oversampling_grid
# add a sampling parameter to the arguments
oversampling_grid(ss::SamplingStyle, ap::ApproximationProblem; options...) =
    oversampling_grid(ss, platform(ap), platformparameter(ap), samplingparameter(ss, ap; options...); options...)
oversampling_grid(ss::SamplingStyle, platform::Platform, param, L; dict = dictionary(platform, param), options...) =
    oversampling_grid(dict, L)
include("oversampling_grid.jl")


export samplingparameter
@ap_interface samplingparameter
@add_samplingstyle_interface samplingparameter
include("samplingparameter.jl")


export samplingoperator, dualsamplingoperator
@ap_interface samplingoperator
@ap_interface dualsamplingoperator
@add_samplingstyle_interface samplingoperator
@add_samplingstyle_interface dualsamplingoperator
include("samplingoperator.jl")



export sampling_grid, platform_grid
@ap_interface sampling_grid
@add_samplingstyle_interface sampling_grid
@approximation_interface platform_grid
include("grids.jl")
default_sampling_grid(ss::SamplingStyle, dict::Dictionary, L; options...) = oversampling_grid(dict, L)
default_platform_grid(ss::SamplingStyle, dict::Dictionary, args...; grid, options...) = grid


export discretemeasure
@approximation_interface discretemeasure
discretemeasure(ss::ProductSamplingStyle, ap::ApproximationProblem; options...) =
    productmeasure(ap_components(discretemeasure, ss, ap; options...)...)
discretemeasure(ss::SamplingStyle, ap::ApproximationProblem; options...) =
    discretemeasure(ss, platform(ap), platformparameter(ap), ap; options...)
discretemeasure(ss::SamplingStyle, platform::Platform, param, ap; options...) =
    discretemeasure(sampling_grid(ss, ap; options...))

export measure
@approximation_interface measure
default_measure(::SamplingStyle, dict::Dictionary; options...) =
    measure(dict)
measure(ss::DiscreteGramStyle, ap::ApproximationProblem; options...) =
    discrete_gram_measure(ss, ap; options...)
@approximation_interface discrete_gram_measure
discrete_gram_measure(ss::SamplingStyleSuper{DiscreteGramStyle}, ap::ApproximationProblem; options...) =
    discretemeasure(ss, ap; options...)
discrete_gram_measure(ss::SamplingStyle, ap::ApproximationProblem; options...) =
    error("Only with `DiscreteGramStyle`")


export azdual_dict
@approximation_interface azdual_dict
include("azdual_dict.jl")


export discretization, dualdiscretization
@approximation_interface discretization
@approximation_interface dualdiscretization
include("discretization.jl")


export normalized_discretization, normalized_dualdiscretization
@approximation_interface normalized_discretization
@approximation_interface normalized_dualdiscretization
include("normalized_discretization.jl")


export solver
@ap_interface solver
include("solver.jl")


## The AZ algorithm
export AZ_A
@ap_interface AZ_A
AZ_A(ap::ApproximationProblem; options...) =
    normalized_discretization(ap; options...)


export AZ_Z
@ap_interface AZ_Z
AZ_Z(ap::ApproximationProblem; options...) =
    normalized_dualdiscretization(ap; options...)


export AZ_Zt
@ap_interface AZ_Zt
AZ_Zt(ap::ApproximationProblem; options...) = AZ_Z(ap; options...)'


export plungeoperator
@ap_interface plungeoperator
function plungeoperator(ap::ApproximationProblem; options...)
    A = AZ_A(ap; options...)
    Zt = AZ_Zt(ap; options...)
    I - A*Zt
end


export plungematrix, firstAZstepoperator
@ap_interface plungematrix
function plungematrix(ap::ApproximationProblem; options...)
    A = AZ_A(ap; options...)
    P = plungeoperator(ap; options...)
    P * A
end
const firstAZstepoperator = plungematrix


export plungerank
@ap_interface plungerank
function plungerank(ap::ApproximationProblem;
            rankestimate = 40,
            options...)
    C = plungematrix(ap; options...)
    Q = rSVD_solver(C; rankestimate = rankestimate, options...)
    length(Q.Sinv)
end

export smoothingoperator
smoothingoperator(ap::ApproximationProblem; options...) =
    WeightedSmoothingOperator(dictionary(ap); options...)


export normalizationoperator
normalizationoperator(ap::ApproximationProblem; sstyle=SamplingStyle(ap), opts...) =
    normalizationoperator(sstyle, ap; opts...)

function normalizationoperator(sstyle::DiscreteStyle, ap::ApproximationProblem; T=coefficienttype(ap), options...)
    grid = sampling_grid(sstyle, ap; options...)
    sampling_normalization(GridBasis{T}(grid),
        measure(sstyle, ap; options...))
end


include("approximate.jl")
