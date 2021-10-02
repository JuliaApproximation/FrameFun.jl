
import BasisFunctions: discretemeasure, measure


# """
#     @approximation_interface `function`
#
# Generate an interface for the given function in terms of an approximation
# problem and a sampling style. These are generated from the arguments of the
# function.
# """
# macro approximation_interface(ex)
#     opt = Meta.parse("options_"*string(ex))
#     def = Meta.parse("default_"*string(ex))
#     ret = quote
#         # generate an approximation problem and pass to options parsing
#         $(ex)(dict::Dictionary, args...; options...) =
#             $(opt)(approximationproblem(dict, args...); options...)
#         $(ex)(platform::Platform, args...; options...) =
#             $(opt)(approximationproblem(platform, args...); options...)
#         # parse a few options and invoke the function again
#         $(opt)(ap::ApproximationProblem; options...) =
#             $(ex)(ap; ap_process_options(ap; options...)...)
#
#         # add a sampling style
#         $(ex)(ap::ApproximationProblem; options...) =
#             $(ex)(SamplingStyle(ap), ap; options...)
#
#         # translate the approximation problem to a dictionary via a platform
#         # these calls can be intercepted by the platform/dictionary to generate custom behaviour
#         $(ex)(ss::SamplingStyle, ap::ApproximationProblem, args...; options...) =
#             $(ex)(ss, platform(ap), platformparameter(ap), args...; options...)
#         $(ex)(ss::SamplingStyle, platform::Platform, plt_par, args...; dict = dictionary(platform, plt_par), options...) =
#             $(ex)(ss, dict, args...; options...)
#         $(ex)(ss::SamplingStyle, platform::ParametrizedPlatform, plt_par, args...; options...) =
#             $(ex)(ss, platform.platform, platform.path[plt_par], args...; options...)
#         $(ex)(ss::SamplingStyle, dict::Dictionary, args...; options...) =
#             $(def)(ss, dict, args...; options...)
#     end
#     esc(ret)
# end

"Compute a property and store it in the ap cache."
macro ap_property(property, default)
    ret = quote
        $(property)(ap::ApproximationProblem; options...) = ap_hasproperty(ap, $(property)) ?
            cache(ap, $(property)) : cache!(ap, $(property), $(default))
    end
    esc(ret)
end
macro ap_property(property)
    compute = Meta.parse("compute_"*string(property))
    ret = quote
        $(property)(ap::ApproximationProblem; options...) = ap_hasproperty(ap, $(property)) ?
            cache(ap, $(property)) : cache!(ap, $(property), $(compute)(ap; options...))
    end
    esc(ret)
end



"""
Return the cached property of an approximation problem, where `property` is a
function.

If the result is not cached yet, compute it by invoking:
property(ap; options...) =
    property(SamplingStyle(ap), platform(ap), platformparameter(ap), samplingparameter(ap); options...)

Then cache the result.
"""
macro ap_sampling_property_to_platform(property)
    ret = quote
        $(property)(ap::ApproximationProblem; options...) = ap_hasproperty(ap, $(property)) ?
            cache(ap, $(property)) : cache!(ap, $(property), $(property)(SamplingStyle(ap), ap; options...))
        $(property)(sstyle::SamplingStyle, ap::ApproximationProblem; options...) =
            $(property)(sstyle, platform(ap), platformparameter(ap), samplingparameter(ap); options...)
    end
    esc(ret)
end

macro ap_sampling_property(property)
    ret = quote
        $(property)(ap::ApproximationProblem; options...) = ap_hasproperty(ap, $(property)) ?
            cache(ap, $(property)) : cache!(ap, $(property), $(property)(SamplingStyle(ap), ap; options...))
        function $(property)(ss::ProductSamplingStyle, ap::ApproximationProblem; options...)
            @assert nfactors(ss) == nfactors(platform(ap))
            to_product( map((st,pl,p_par,s_par)->$(property)(pl, p_par, s_par; samplingstyle=st, options...),
                factors(ss), factors(platform(ap)), productparameter(ap), samplingparameter(ap))...)
        end
    end
    esc(ret)
end

to_product(args...) = tensorproduct(args...)


"""
    @ap_interface `function`

Generate an interface for the given function in terms of an approximation
problem. The approximation problem is constructed from the arguments of the
function.
"""
macro ap_interface(ex)
    opt = Meta.parse("options_"*string(ex))
    ret = quote
        # generate an approximation problem and pass to options parsing
        $(ex)(dict::Dictionary, args...; options...) =
            $(opt)(approximationproblem(dict, args...); options...)
        $(ex)(platform::Platform, args...; options...) =
            $(opt)(approximationproblem(platform, args...); options...)
        # Accept a function as the first argument
        $(ex)(fun, dict::Dictionary, args...; options...) =
            $(opt)(approximationproblem(fun, dict, args...); options...)
        $(ex)(fun, platform::Platform, args...; options...) =
            $(opt)(approximationproblem(fun, platform, args...); options...)
        # parse a few options and invoke the function again
        $(opt)(ap::ApproximationProblem; options...) =
            $(ex)(ap; ap_process_options(ap; options...)...)
    end
    esc(ret)
end

"A generic routine to process a few arguments given to the Fun interface."
function ap_generic_options(ap::ApproximationProblem;
        samplingstyle = -1,
        solverstyle = -1,
        normalizedsampling = false,
        options...)

    if samplingstyle == -1
        samplingstyle = SamplingStyle(platform(ap))
    end
    if normalizedsampling && (samplingstyle isa DiscreteStyle)
        samplingstyle = NormalizedSampling(samplingstyle)
    end
    if solverstyle == -1
        solverstyle = SolverStyle(platform(ap), samplingstyle)
    end
    cache!(ap, SamplingStyle, samplingstyle)
    cache!(ap, SolverStyle, solverstyle)
    options
end

ap_process_options(ap::AdaptiveApproximation; options...) =
    ap_generic_options(ap; options...)

function ap_process_options(ap::PlatformApproximation; options...)
    opt = ap_generic_options(ap; options...)
    # compute and cache the sampling parameter
    samplingparameter(ap; opt...)
    opt
end


# macro add_samplingstyle_interface(ex)
#     ret = quote
#         $(ex)(ap::ApproximationProblem; options...) =
#             $(ex)(SamplingStyle(ap), ap; options...)
#     end
#     esc(ret)
# end


# Utility functions to create tensor FrameFunInterface functions

function ap_components(fun, ss::ProductSamplingStyle, ap::ApproximationProblem, args...; options...)
    @assert length(components(ss))==length(factors(ap))
    map((ssi, api, argi)->fun(ssi, api, argi...; options...), components(ss), factors(ap), zip(args...))
end

function ap_components(fun, ss::ProductSamplingStyle, ap::ApproximationProblem; options...)
    @assert length(factors(ss))==length(factors(ap))
    map((ssi, api)->fun(ssi, api; options...), factors(ss), factors(ap))
end


default_aznormalization(a...) = false


## Basic properties

# - SamplingStyle and SolverStyle
@ap_property SamplingStyle SamplingStyle(platform(ap))
@ap_property SolverStyle SolverStyle(platform(ap), SamplingStyle(ap))

# - sampling parameter
export samplingparameter
@ap_interface samplingparameter

function samplingparameter(ap::ApproximationProblem; options...)
    ap_hasproperty(ap, samplingparameter) ? cache(ap, samplingparameter) : begin
        smpl_par = samplingparameter(SamplingStyle(ap), platform(ap), platformparameter(ap); options...)
        if platform(ap) isa ProductPlatform
            @assert length(smpl_par) == length(productparameter(ap))
        end
        cache!(ap, samplingparameter, smpl_par)
    end
end

## Other properties are implemented in separate files
include("samplingparameter.jl")
include("grids.jl")
include("samplingoperator.jl")


export discretemeasure
@ap_interface discretemeasure
discretemeasure(ap; options...) = discretemeasure(sampling_grid(ap; options...), sampling_weights(ap; options...))


export azdual
@ap_interface azdual
@ap_sampling_property azdual
include("azdual.jl")


export discretization, dualdiscretization
@ap_interface discretization
@ap_property discretization apply(sampling_operator(ap), dictionary(ap); options...)

@ap_interface dualdiscretization
@ap_property dualdiscretization apply(dual_sampling_operator(ap), dualdictionary(ap); options...)

export sample_data
@ap_interface sample_data
@ap_property sample_data apply(sampling_operator(ap), ap_fun(ap); options...)

export full_discretization
@ap_interface full_discretization
@ap_property full_discretization (discretization(ap; options...), sample_data(ap; options...))

# include("discretization.jl")

#
# export normalized_discretization, normalized_dualdiscretization
# @ap_interface normalized_discretization
# @ap_interface normalized_dualdiscretization
# include("normalized_discretization.jl")


export solver
@ap_interface solver
include("solver.jl")


## The AZ algorithm
export AZ_A
@ap_interface AZ_A
AZ_A(ap::ApproximationProblem; options...) = discretization(ap; options...)


export AZ_Z
@ap_interface AZ_Z
AZ_Z(ap::ApproximationProblem; options...) = dualdiscretization(ap; options...)


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



include("approximate.jl")
