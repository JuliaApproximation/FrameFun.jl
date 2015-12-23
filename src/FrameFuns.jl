# FrameFuns.jl

module FrameFuns

using FixedSizeArrays
using BasisFunctions
using Debug
#using ApproxFun

using Base.Cartesian

import Base: +, *, /, ==, |, &, -, \

import Base: length, eltype, size, push!

import Base: eachindex, start, next, done, getindex, in

import Base: show, showcompact, call, convert

import BasisFunctions: ⊗

import BasisFunctions: src, dest, matrix, matrix!, apply!, numtype

import BasisFunctions: dim, index_dim, grid, left, right

import BasisFunctions: operator, coefficients, set, is_basis, is_frame, transform_normalization_operator

import BasisFunctions: call_set, call_set!, call_expansion, call_expansion!

import BasisFunctions: True, False

import BasisFunctions: tp_length, left, right

import BasisFunctions: show_setexpansion

# from box.jl
export BBox, left, right

# from domains.jl
export Interval, Circle, Square, Cube, Sphere, Cylinder, atomium, ⊗, box

export numtype

# from funs.jl
export ExpFun, ChebyFun, Fun, FrameFun

# from domainframe.jl
export DomainFrame, basis, call_set, call_set!

# from fractal.jl
export Mandelbrot, JuliaSet


include("box.jl")

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


