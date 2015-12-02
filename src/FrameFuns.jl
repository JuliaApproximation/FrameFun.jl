# FrameFuns.jl

module FrameFuns

using FixedSizeArrays
using BasisFunctions
using Debug

using Base.Cartesian

import Base: +, *, /, ==, |, &, -, \

import Base: length, eltype, size, push!

import Base: eachindex, start, next, done, getindex, in

import Base: show, showcompact, call, convert

import BasisFunctions: src, dest, matrix, matrix!, apply!, numtype

import BasisFunctions: dim, index_dim, grid, left, right

import BasisFunctions: operator, coefficients, set, is_basis, is_frame, normalization_operator

import BasisFunctions: call_set, call_set!, call_expansion, call_expansion!

import BasisFunctions: True, False


# from domains.jl
export Interval, Circle, Square, Cube, Sphere, Cylinder, atomium

export numtype

# from funs.jl
export ExpFun, ChebyFun, FrameFun

# from domainframe.jl
export DomainFrame, basis, call_set, call_set!


include("domains.jl")

include("subgrid.jl")

include("domainframe.jl")

include("fe_problem.jl")

include("fe_solvers.jl")

include("fastsolver.jl")

include("funs.jl")

include("fe_fourier.jl")

# TODO: try out Plots.jl
FrameFuns.include("plots.jl")

end # module


