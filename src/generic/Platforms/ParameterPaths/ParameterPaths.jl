module ParameterPaths

using ..Platforms, BasisFunctions
import Base: getindex, unsafe_getindex, first
import ..Platforms: unsafe_dictionary, correctparamformat, param_first, param_double,
    param_inbetween, param_increment, dualdictionary, measure, SamplingStyle, SolverStyle, ProblemStyle

# For the moment only size is in parameter space. (no spline degree, extensionparameter T, ...)
export ParameterPath
"""
    abstract type ParameterPath{OUT} end

An indexable type that returns something of type `OUT`.
The index can have a lot of complexity, e.g.:
A dictionary expansion and a function together determine a platform parameter
that leads to a expansion with a smaller error, but more degrees of freedom.
"""
abstract type ParameterPath{OUT} end

correctparamformat(::ParameterPath, param) = false

export outtype
outtype(::ParameterPath{OUT}) where OUT = OUT

function getindex(path::ParameterPath, index)
    if !correctparamformat(path, index)
        return ArgumentError("Index $index not compatible with parameterpath $path")
    end
    unsafe_getindex(path, index)
end

export default_param_path
default_param_path(platform::Platform) =
    throw(ArgumentError("`default_param_path` not implemented for platform of type $(typeof(platform))"))

first(path::ParameterPath) = unsafe_getindex(path,  param_first(path))


export LinearParameterPath1D
"""
    struct LinearParameterPath1D{T<:Real} <: ParameterPath{Int}

"""
struct LinearParameterPath1D{T<:Real} <: ParameterPath{Int}
    constant    ::  T
end

LinearParameterPath1D() = LinearParameterPath1D(1)
param_first(path::LinearParameterPath1D) = 1
param_double(path::LinearParameterPath1D, index::Real) =
    2index
param_inbetween(path::LinearParameterPath1D, index1::Real, index2::Real) =
    (index1+index2)/2

param_increment(path::LinearParameterPath1D, param::Real) =
    param + 1/path.constant

unsafe_getindex(path::LinearParameterPath1D, index::Real) =
    round(Int, index*path.constant)

correctparamformat(::LinearParameterPath1D, ::Real) = true

export CartesianParameterPath
"""
    struct CartesianParameterPath{TNTuple{N,Int}} <: ParameterPath{NTuple{N,Int}}

"""
struct CartesianParameterPath{N,T<:NTuple{N,<:Real}} <: ParameterPath{NTuple{N,Int}}
    constant    ::  T
end

CartesianParameterPath{N}() where N = CartesianParameterPath(CartesianIndex{N}(1).I)

param_first(path::CartesianParameterPath) = 1
param_double(path::CartesianParameterPath, param::Real) =
    2param
param_inbetween(path::CartesianParameterPath{N}, param1::Real, param2::Real) where N =
    max(param1, min((param1+param2)/2, minimum(((((param2).^(1/N)).*path.constant.-1)./path.constant).^N)))
param_double(path::CartesianParameterPath, param::NTuple) =
    (2^(1/N)).*param
param_inbetween(path::CartesianParameterPath{N}, param1::NTuple, param2::NTuple) where N =
    max.(param1, min.((param1 .+ param2)./2, (((((param2).^(1/N)).*path.constant.-1)./path.constant).^N)))
param_increment(path::CartesianParameterPath{N}, param::Real) where N =
    minimum((round.(Int,((param).^(1/N)).*path.constant.+1)./path.constant).^N)
param_increment(path::CartesianParameterPath{N}, param::NTuple{N}) where N =
    ((round.(Int,((param).^(1/N)).*path.constant.+1)./path.constant).^N)


unsafe_getindex(path::CartesianParameterPath{N}, index::Real) where N =
    round.(Int, (index^(1/N)).*path.constant)
unsafe_getindex(path::CartesianParameterPath{N}, index::Vararg{N,<:Real}) where N =
    unsafe_getindex(path, index)
unsafe_getindex(path::CartesianParameterPath{N}, index::CartesianIndex{N}) where N =
    unsafe_getindex(path, index.I)
unsafe_getindex(path::CartesianParameterPath{N}, index::NTuple{N,<:Real}) where N =
    round.(Int, (index).*path.constant)
correctparamformat(::CartesianParameterPath, ::Real) = true
correctparamformat(::CartesianParameterPath{N}, ::NTuple{N,<:Real}) where N = true
correctparamformat(::CartesianParameterPath{N}, ::CartesianIndex{N}) where N = true
correctparamformat(::CartesianParameterPath{N}, ::Vararg{N,<:Real}) where N = true


export IncrementalCartesianParameterPath
"""
    struct IncrementalCartesianParameterPath{TNTuple{N,Int}} <: ParameterPath{NTuple{N,Int}}

"""
struct IncrementalCartesianParameterPath{N} <: ParameterPath{NTuple{N,Int}}
end

param_first(path::IncrementalCartesianParameterPath) = 1
param_double(path::IncrementalCartesianParameterPath, param::Real) =
    2param
param_inbetween(path::IncrementalCartesianParameterPath{N}, param1::Real, param2::Real) where N =
    max(param1, min((param1+param2)/2, param2-1))
param_increment(path::IncrementalCartesianParameterPath{N}, param::Real) where N =
    param + 1


function unsafe_getindex(path::IncrementalCartesianParameterPath{N}, index::Real) where N
    q, r = divrem(round(Int, index-1), N)
    ntuple(k-> k<=r ? q + 2 : q +1, Val(N))
end

correctparamformat(::IncrementalCartesianParameterPath, ::Real) = true




export NTupleParameterPath
struct NTupleParameterPath{N,OUT} <: ParameterPath{NTuple{N,OUT}}
    first::OUT
    NTupleParameterPath{N}(out::OUT) where {N,OUT} = new{N,OUT}(out)
end

correctparamformat(::NTupleParameterPath{N,OUT}, param) where {N,OUT} = param isa OUT

unsafe_getindex(path::NTupleParameterPath{N}, param) where N =
    ntuple(k->param, Val(N))

export HierarchyPath
"""
"""
struct HierarchyPath{OUT,PATH1<:ParameterPath,PATH2<:ParameterPath{OUT}} <: ParameterPath{OUT}
    pathfirst :: PATH1
    pathsecond :: PATH2
    function HierarchyPath(pathfirst::ParameterPath, pathsecond::ParameterPath)
        if !correctparamformat(pathsecond,first(pathfirst))
            error("parameters produced by the first path are not compatible with the input of the second path.")
        end
        new{outtype(pathsecond),typeof(pathfirst),typeof(pathsecond)}(pathfirst,pathsecond)
    end
end

correctparamformat(path::HierarchyPath, param) =
    correctparamformat(path.pathfirst, param)

unsafe_getindex(path::HierarchyPath, param) =
    unsafe_getindex(path.pathsecond,path.pathfirst[param])

param_first(path::HierarchyPath) =
    param_first(path.pathfirst)

param_double(path::HierarchyPath, param) =
    param_double(path.pathfirst, param)

param_inbetween(path::HierarchyPath, param1, param2) =
    param_inbetween(path.pathfirst, param1, param2)

param_increment(path::HierarchyPath, param1) =
    param_increment(path.pathfirst, param1)

export ProductPath
"""
"""
struct ProductPath{N,PATHS,OUT} <: ParameterPath{OUT}
    paths :: PATHS
    function ProductPath(paths::ParameterPath...)
        new{length(paths),typeof(paths),typeof(map(first, paths))}(paths)
    end
end

correctparamformat(path::ProductPath, param) =
    all(map(correctparamformat,path.paths, param))

unsafe_getindex(path::ProductPath, param) =
    map(unsafe_getindex, path.paths, param)

param_first(path::ProductPath) =
    map(param_first, path.paths)

param_double(path::ProductPath, param) =
    map(param_double, path.paths, param)

param_inbetween(path::ProductPath, param1, param2) =
    map(param_inbetween, path.paths, param1, param2)

default_param_path(platform::ProductPlatform) where N =
    ProductPath(map(default_param_path, elements(platform))...)

struct ParametrizedPlatform{PLATFORM<:Platform, PATH<:ParameterPath} <: Platform
    platform :: PLATFORM
    path     :: PATH

    function ParametrizedPlatform(platform::Platform, path::ParameterPath)
        if !correctparamformat(platform, first(path))
            throw(ArgumentError("The given path returns no parameters that are compatible with the given platform"))
        end
        new{typeof(platform),typeof(path)}(platform, path)
    end
end

export parametrizedplatform
parametrizedplatform(platform::Platform) =
    parametrizedplatform(platform, default_param_path(platform))

parametrizedplatform(platform::Platform, path::ParameterPath) =
    ParametrizedPlatform(platform, path)

correctparamformat(platform::ParametrizedPlatform, param) =
    correctparamformat(platform.path, param)

unsafe_dictionary(platform::ParametrizedPlatform, param) =
    unsafe_dictionary(platform.platform, platform.path[param])

default_param_path(platform::ParametrizedPlatform) = platform.path

param_first(platform::ParametrizedPlatform) =
    param_first(platform.path)

param_double(platform::ParametrizedPlatform, param) =
    param_double(platform.path, param)

param_inbetween(platform::ParametrizedPlatform, param1, param2) =
    param_inbetween(platform.path, param1, param2)

param_increment(platform::ParametrizedPlatform, param1) =
    param_increment(platform.path, param1)

dualdictionary(platform::ParametrizedPlatform, param, measure::Measure; options...) =
    dualdictionary(platform.platform, param, measure; options...)

measure(platform::ParametrizedPlatform) =
    measure(platform.platform)

SamplingStyle(platform::ParametrizedPlatform) =
    SamplingStyle(platform.platform)

SolverStyle(platform::ParametrizedPlatform, ss::SamplingStyle) =
    SolverStyle(platform.platform, ss)

ProblemStyle(platform::ParametrizedPlatform) =
    ProblemStyle(platform.platform)

end # module
