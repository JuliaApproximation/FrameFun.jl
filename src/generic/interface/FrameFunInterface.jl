
import BasisFunctions: discretemeasure, measure, interpolation_grid, components

export SamplingStrategy
abstract type SamplingStrategy end

export DefaultSamplingStrategy
"""
The default strategy tries to keep the oversamplingfactor correct by distributing equally over all degrees of freedom.

"""
struct DefaultSamplingStrategy <: SamplingStrategy end

default_aznormalization(a...) = false

export smoothingoperator
smoothingoperator(ap::ApproximationProblem; options...) =
    WeightedSmoothingOperator(dictionary(ap); options...)

export normalizationoperator
normalizationoperator(::DictionaryOperatorStyle, ap::ApproximationProblem; sstyle=SamplingStyle(ap), opts...) =
    normalizationoperator(sstyle, ap; opts...)

function normalizationoperator(sstyle::DiscreteStyle, ap::ApproximationProblem; T=coefficienttype(ap), options...)
    grid = sampling_grid(sstyle, ap; options...)
    sampling_normalization(GridBasis{T}(grid),
        measure(sstyle, ap; options...))
end

# Utility functions to create tensor FrameFunInterface functions

function components(fun, ss::ProductSamplingStyle, ap::ProductPlatformApproximation, args...; options...)
    samplingparameter(ss, ap; options...)
    @assert length(components(ss))==length(components(ap))
    map((ssi, api, argi)->fun(ssi, api, argi...; options...), components(ss), components(ap), zip(args...))
end

function components(fun, ss::ProductSamplingStyle, ap::ProductPlatformApproximation; options...)
    samplingparameter(ss, ap; options...)
    @assert length(components(ss))==length(components(ap))
    map((ssi, api)->fun(ssi, api; options...), components(ss), components(ap))
end


INTERFACEFUNCTIONS = (:samplingoperator, :dualsamplingoperator, :samplingparameter,
    :sampling_grid, :interpolation_grid, :platform_grid, :discretization, :dualdiscretization,
    :azdual_dict, :measure, :discretemeasure, :solver, :AZ_A, :AZ_Z, :AZ_Zt, :plungeoperator,
    :plungematrix, :plungerank, :smoothingoperator)


# parameter                                     # rename to
# platform: the platform to approximate with    #
# dict: the dict to approximate with            #
# M: grid size of problem                       # M
# L: grid size of bounding box                  # sparam
# param: index of platform                      # pparam
# N: (length(dict))

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
        function $(ex)(ap::ApproximationProblem, args...; samplingstyle=SamplingStyle(ap), options...)
            # compute the sampling parameter to make sure it is cached
            samplingparameter(samplingstyle, ap; options...)
            $(ex)(samplingstyle, ap, args...; options...)
        end

        # translate the approximation problem to a dictionary via a platform
        # these calls can be intercepted by the platform to generate custom behaviour
        $(ex)(ss::SamplingStyle, ap::ApproximationProblem, args...; options...) =
            $(ex)(ss, platform(ap), parameter(ap), args...; options...)
        $(ex)(ss::SamplingStyle, platform::Platform, param, args...; dict = dictionary(platform, param), options...) =
            $(ex)(ss, dict, args...; options...)
        $(ex)(ss::SamplingStyle, platform::ParametrizedPlatform, param, args...; options...) =
            $(ex)(ss, platform.platform, platform.path[param], args...; options...)
        $(ex)(ss::SamplingStyle, dict::Dictionary, args...; options...) =
            $(def)(ss, dict, args...; options...)
    end
    esc(ret)
end

"""
    @basic_interface `function`

Generate an interface for the given function in terms of an approximation
problem. The approximation problem is constructed from the arguments of the
function.
"""
macro basic_interface(ex)
    ret = quote
        # generate an approximation problem
        $(ex)(dict::Dictionary, args...; options...) = $(ex)(approximationproblem(dict, args...); options...)
        $(ex)(platform::Platform, args...; options...) = $(ex)(approximationproblem(platform, args...); options...)

        function $(ex)(ap::ApproximationProblem, args...; samplingstyle=SamplingStyle(ap),options...)
            # compute the sampling parameter once to make sure it is cached
            samplingparameter(samplingstyle, ap; options...)
            $(ex)(samplingstyle, ap, args...; options...)
        end
    end
    esc(ret)
end

macro add_samplingparameter_interface(ex)
    intermediate = Meta.parse("_"*string(ex))
    ret = quote
        $(ex)(ss::SamplingStyle, ap::ApproximationProblem, args...; options...) =
            $(intermediate)(ss, ap, samplingparameter(ss, ap; options...), args...; options...)
        $(intermediate)(ss::SamplingStyle, ap::ApproximationProblem, L, args...; options...) =
            $(ex)(ss, platform(ap), parameter(ap), L, args...; options...)
    end
    esc(ret)
end

macro add_samplingstyle_interface(ex)
    ret = quote
        $(ex)(ap::ApproximationProblem; samplingstyle = SamplingStyle(ap), options...) =
            $(ex)(samplingstyle, ap; options...)
    end
    esc(ret)
end

macro addproblemstyle(ex)
    ret = quote
        $(ex)(ap::ApproximationProblem; problemstyle = ProblemStyle(ap), options...) =
            $(ex)(problemstyle, ap; options...)
    end
    esc(ret)
end



# We delegate some functionality from the platform to the dictionary corresponding
# to a value of the platform parameter param.
# Platforms may choose to override these to enable different behaviour.
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
    oversampling_grid(ss, platform(ap), parameter(ap), samplingparameter(ss, ap; options...); options...)
oversampling_grid(ss::SamplingStyle, platform::Platform, param, L; dict = dictionary(platform, param), options...) =
    oversampling_grid(dict, L)
include("oversampling_grid.jl")


export samplingparameter
@basic_interface samplingparameter
@add_samplingstyle_interface samplingparameter
include("samplingparameter.jl")


export samplingoperator, dualsamplingoperator
@basic_interface samplingoperator
@basic_interface dualsamplingoperator
include("samplingoperator.jl")



export sampling_grid, platform_grid
@basic_interface sampling_grid
@approximation_interface platform_grid
include("grids.jl")
default_sampling_grid(ss::SamplingStyle, dict::Dictionary, L; options...) = oversampling_grid(dict, L)
default_platform_grid(ss::SamplingStyle, dict::Dictionary, args...; grid, options...) = grid


export discretemeasure
@approximation_interface discretemeasure
discretemeasure(ss::ProductSamplingStyle, ap::ApproximationProblem; options...) =
    productmeasure(components(discretemeasure, ss, ap; options...)...)
discretemeasure(ss::SamplingStyle, ap::ApproximationProblem; options...) =
    discretemeasure(ss, platform(ap), parameter(ap), ap; options...)
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
@basic_interface solver
include("solver.jl")


## The AZ algorithm
export AZ_A
@basic_interface AZ_A
@addproblemstyle AZ_A
AZ_A(pstyle::ProblemStyle, ap::ApproximationProblem; options...) =
    normalized_discretization(ap; options...)


AZ_A(::GenericOperatorStyle, ap::ApproximationProblem; options...) =
    SynthesisOperator(dictionary(ap), measure(ap; options...))


export AZ_Z
@basic_interface AZ_Z
@addproblemstyle AZ_Z
AZ_Z(pstyle::ProblemStyle, ap::ApproximationProblem; options...) =
    normalized_dualdiscretization(ap; options...)
AZ_Z(::GenericOperatorStyle, ap::ApproximationProblem; options...) =
    SynthesisOperator(azdual_dict(ap; options...), measure(ap; options...))


export AZ_Zt
@basic_interface AZ_Zt
@addproblemstyle AZ_Zt
AZ_Zt(pstyle::ProblemStyle, ap::ApproximationProblem; options...) = AZ_Z(pstyle, ap; options...)'
samplingoperator(pstyle::GenericOperatorStyle, ap; options...) = AZ_Zt(pstyle, ap; options...)
samplingoperator(pstyle::DictionaryOperatorStyle, ap; options...) = samplingoperator(ap; options...)


export plungeoperator
@basic_interface plungeoperator
@addproblemstyle plungeoperator
function plungeoperator(problemstyle::ProblemStyle, ap::ApproximationProblem; options...)
    A = AZ_A(problemstyle, ap; options...)
    Zt = AZ_Zt(problemstyle, ap; options...)
    I - A*Zt
end


export plungematrix, firstAZstepoperator
@basic_interface plungematrix
@addproblemstyle plungematrix
function plungematrix(pstyle::ProblemStyle, ap::ApproximationProblem; options...)
    A = AZ_A(pstyle, ap; options...)
    P = plungeoperator(pstyle, ap; options...)
    P * A
end
function plungematrix(pstyle::GenericOperatorStyle, ap::ApproximationProblem; options...)
    S = samplingoperator(pstyle, ap; options...)
    A = AZ_A(pstyle, ap; options...)
    Zt = AZ_Zt(pstyle, ap; options...)
    tmp1 = apply(Zt, A; options...)
    tmp2 = apply(S, A; options...)
    tmp2 - tmp2*tmp1
end
const firstAZstepoperator = plungematrix


export plungerank
@basic_interface plungerank
function plungerank(ap::ApproximationProblem;
            rankestimate = 40,
            options...)
    C = plungematrix(ap; options...)
    Q = rSVD_solver(C; rankestimate = rankestimate, options...)
    length(Q.Sinv)
end


import BasisFunctions: approximate
export Fun
include("approximate.jl")
