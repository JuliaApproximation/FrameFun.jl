# FrameFun.jl
__precompile__(true)
module FrameFun

using Base.Cartesian
using StaticArrays
using RecipesBase

using Domains
using BasisFunctions
using Combinatorics: permutations

if VERSION < v"0.7-"
    import BasisFunctions.adjoint
    mul! = A_mul_B!
    eigen  = eig
else
    using LinearAlgebra
    using Printf
    import LinearAlgebra: adjoint
end

###############################
## Exhaustive list of imports
###############################

## - Imports from Base
import Base: +, *, /, ^, ==, |, &, -, \, <, <=, >, >=

import Base: intersect, union, isapprox, setdiff, in

# Arrays
import Base: length, eltype, size, push!, similar
import Base: inv

# Iteration and indexing
import Base: eachindex, start, next, done, getindex, unsafe_getindex,
    checkbounds

# Display
import Base: show, showcompact

import Base: promote, promote_rule, convert, promote_eltype

import Base: broadcast

import Base: minimum, maximum

## - Imports from Domains
import Domains: indomain, dimension
# - for mapped domains
import Domains: domain
# - for composite structures
import Domains: element, elements, numelements
# - for cartesian products
import Domains: cartesianproduct, ×


## - Imports from BasisFunctions
import BasisFunctions: src, dest, matrix, matrix!, apply!, apply_inplace!,
    dimension, codomaintype, domaintype

import BasisFunctions: tensorproduct, ⊗

import BasisFunctions: grid, left, right, stepsize, sample

import BasisFunctions: is_basis, is_frame, is_orthogonal, is_orthonormal, is_biorthogonal,
    has_transform, has_grid, has_derivative,
    has_antiderivative, has_extension, has_grid_transform

import BasisFunctions: operator, matrix, is_diagonal, is_inplace, ⊕

import BasisFunctions: coefficients, dictionary,
    transform_operator_pre, transform_operator_post, evaluation_operator, interpolation_operator,
    differentiation_operator, antidifferentiation_operator, approximation_operator,
    extend, extension_size, extension_operator, restriction_operator,
    default_approximation_operator, has_extension, wrap_operator, grid_evaluation_operator

import BasisFunctions: superdict, similar_dictionary,
    promote_domaintype, promote_domainsubtype

import BasisFunctions: eval_element, eval_element_derivative, eval_expansion,
    unsafe_eval_element, unsafe_eval_element_derivative,
    name, in_support, dict_in_support, domain

import BasisFunctions: differentiate, ∂x, ∂y, ∂z, ∫∂x, ∫∂y, ∫∂z, ∫, is_compatible

import BasisFunctions: True, False, resize, promote_eltype

import BasisFunctions: show_setexpansion

import BasisFunctions: postprocess, plotgrid

import BasisFunctions: flatten

import BasisFunctions: Span

import BasisFunctions: coefficient_type, coeftype

# about subgrids
import BasisFunctions: AbstractSubGrid, IndexSubGrid, is_subindex, supergrid,
    similar_subgrid, mask, subindices

import BasisFunctions: Gram, DualGram, MixedGram, DiscreteGram, DiscreteDualGram,
    DiscreteMixedGram, gram_entry
import BasisFunctions: dual, primal, sampler, dual_sampler

import BasisFunctions: discrete_approximation_operator, continuous_approximation_operator

import BasisFunctions: primal, dual, Zt, A, sampler, dual_sampler

###############################
## Exhaustive list of exports
###############################
# from funs.jl
export ExpFun, ChebyFun, Fun, DictFun, sampling_grid, domain, abserror, maxerror, L2error, expansion

# from subgrid.jl
export MaskedGrid, boundary_grid

# from domains/boundingbox.jl
export BoundingBox, BBox, BBox1, BBox2, BBox3, BBox4
export boundingbox

# from domains/extensions.jl
export dist, normal

# from frames/extensionframe.jl
export ExtensionFrame, basis, domain, extensionframe
export Gram, DualGram, MixedGram
export extension_frame_platform

# from frames/weighted_sum_frame.jl
export WeightedSumFrame, weightfunctions
export hassuperdict

# from DiffEquation.jl
export DirichletBC, NeumannBC, DiffEquation, solve

# from constructors.jl
export FunConstructor

# from space.jl
export FourierSpace, ChebyshevSpace, ⊕, add, construct
# from recipes.jl

# from randomgrid.jl
export randomgrid, randompoint

# from domains/fractals.jl
export mandelbrot, juliaset

# from domains/atomium.jl
export atomium

# from domains/polardomain.jl
export polardomain

# from domains/characteristic.jl
export characteristic

export FeFun, FeFunNd

# from diffequation.jl
export operator
export FECollocationOperator

# from approximation
export DirectSolver, AZSolver, AZSmoothSolver, TridiagonalSolver, AZSSolver, AZSDCSolver

include("sampling/subgrid.jl")


#include("domains/boundingbox.jl")
include("domains/extensions.jl")

include("frames/extensionframe.jl")
include("frames/weighted_sum_frame.jl")

include("fun/basisdomains.jl")
include("fun/funs.jl")


include("approximation/fe_solvers.jl")
include("approximation/continuous_solver.jl")
include("approximation/lowranksolver.jl")
include("approximation/azsolver.jl")
include("approximation/tridiagonalsolver.jl")
include("approximation/smoothsolver.jl")
include("approximation/fe_approx.jl")

include("approximation/oversampling.jl")

include("recipes.jl")

include("diffequation.jl")

include("approximation/space.jl")
include("approximation/constructors.jl")

include("domains/fourierdomains.jl")
include("domains/fractals.jl")
include("domains/atomium.jl")
include("domains/characteristic.jl")
include("domains/polardomain.jl")

include("sampling/randomgrid.jl")
include("sampling/oversampledgrid.jl")

include("fun/fefun.jl")


end # module
