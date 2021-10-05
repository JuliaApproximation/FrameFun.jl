
import BasisFunctions: discretemeasure, measure


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
@ap_property discretemeasure discretemeasure(sampling_grid(ap; options...), sampling_weights(ap; options...))

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


export solver
@ap_interface solver
include("solver.jl")

include("az.jl")
include("approximate.jl")

export smoothingoperator
smoothingoperator(ap::ApproximationProblem; options...) =
    WeightedSmoothingOperator(dictionary(ap); options...)
