
###################
# Helper functions
###################

determine_return_type(fun, ::Type{S}) where {S} = Base.Core.Compiler.return_type(fun, (S,))
determine_return_type(fun, S::Type{<:Tuple}) = Base.Core.Compiler.return_type(fun, S)

interpolation_grid(dict::Dictionary; options...) = grid(dict)
interpolation_grid(platform::Platform, param, dict::Dictionary; options...) = interpolation_grid(dict; options...)

function directsolver(A; directsolver = :svd, verbose = false, options...)
    if directsolver == :svd
        verbose && println("Using direct solver: SVD")
        SVD_solver(A)
    elseif directsolver == :qr
        verbose && println("Using direct solver: QR")
        QR_solver(A)
    else
        DS = directsolver(A)
        verbose && println("Using direct solver: $DS")
        DS
    end
end

# TODO: clean up these scaling factors. They have to do with a normalization of the
# sampling operators.
scaling_factor(S::Dictionary, A) = length(supergrid(grid(dest(A))))
scaling_factor(S::DerivedDict, A) = scaling_factor(superdict(S), A)
scaling_factor(S::ChebyshevBasis, A) = length(supergrid(grid(dest(A))))/2

default_Zt(dict::Dictionary, A) = (one(eltype(A))/convert(eltype(A), scaling_factor(dict, A))) * A'

DiscretizationStyle(dict::Dictionary) = is_basis(dict) ? InterpolationStyle() : OversamplingStyle()

SolverStyle(dict::Dictionary, dstyle::InterpolationStyle) = has_transform(dict) ? TransformStyle() : DirectStyle()
SolverStyle(dict::Dictionary, dstyle::OversamplingStyle) = AZStyle()


##############################
# Dictionaries from platforms
##############################

"""
A platformdictionary groups a dictionary with the platform and parameter value
it originated from.
"""
struct PlatformDictionary
    platform    ::  Platform
    param
    dict        ::  Dictionary
end

domaintype(dict::PlatformDictionary) = domaintype(dict.dict)

DiscretizationStyle(dict::PlatformDictionary) = DiscretizationStyle(dict.platform)
SolverStyle(dict::PlatformDictionary, ::DiscretizationStyle) = SolverStyle(dict.platform)

for op in (:approximation_grid, :interpolation_grid, :oversamplingfactor)
    @eval $op(pdict::PlatformDictionary; options...) = $op(pdict.platform, pdict.param, pdict.dict; options...)
end

const AnyDictionary = Union{Dictionary,PlatformDictionary}

Base.getindex(platform::Platform, param) = PlatformDictionary(platform, param, Dictionary(platform, param))


####################
# The Fun interface
####################


function Fun(dict::Dictionary, fun;
            discretizationstyle = DiscretizationStyle(dict),
            solverstyle = SolverStyle(dict, discretizationstyle),
            coefficienttype = promote_type(codomaintype(dict), determine_return_type(fun, domaintype(dict))),
            options...)
    approximate(discretizationstyle, solverstyle, promote_coefficient_type(dict, coefficienttype), fun; options...)
end

function Fun(platform::Platform, param, fun; verbose = false, options...)
    pdict = platform[param]
    verbose && println("Platform: using dictionary $(pdict.dict)")

    Fun(pdict, fun; verbose = verbose, options...)
end

function Fun(dict::PlatformDictionary, fun;
            discretizationstyle = DiscretizationStyle(dict),
            solverstyle = SolverStyle(dict, discretizationstyle),
            options...)
    approximate(discretizationstyle, solverstyle, dict, fun; options...)
end


# Construct the approximation problem and solve it
function approximate(discretizationstyle::DiscretizationStyle, solverstyle::SolverStyle, dict::AnyDictionary, fun;
            verbose = false, options...)
    verbose && println("Approximate: discretizing with style $discretizationstyle")
    A, B = discretization(discretizationstyle, dict, fun; verbose = verbose, options...)

    verbose && println("Approximate: solving with style $solverstyle")
    C = solve(solverstyle, dict, A, B; verbose = verbose, options...)
    combine(dict, C)
end

combine(dict::Dictionary, C) = DictFun(dict, C)
combine(pdict::PlatformDictionary, C) = DictFun(pdict.dict, C)


discretization(dict::Dictionary, fun; discretizationstyle = DiscretizationStyle(dict), options...) =
    discretization(discretizationstyle, dict, fun; options...)

discretization(::InterpolationStyle, dict::Dictionary, fun; options...) =
    grid_discretization(dict, fun, interpolation_grid(dict; options...))

discretization(::InterpolationStyle, pdict::PlatformDictionary, fun; options...) =
    grid_discretization(pdict.dict, fun, interpolation_grid(pdict; options...))

discretization(::GridStyle, dict::Dictionary, fun; grid, options...) =
    grid_discretization(dict, fun, grid)

discretization(::GridStyle, pdict::PlatformDictionary, fun;
            grid = approximation_grid(pdict), options...) =
    grid_discretization(pdict.dict, fun, grid)

function discretization(::OversamplingStyle, dict::Dictionary, fun;
            oversamplingfactor = 2, options...)
    grid = oversampled_grid(dict; oversamplingfactor = oversamplingfactor, options...)[1]
    grid_discretization(dict, fun, grid)
end

function discretization(::OversamplingStyle, pdict::PlatformDictionary, fun;
        oversamplingfactor = oversamplingfactor(pdict), verbose = false, options...)
    grid = oversampled_grid(pdict.dict; oversamplingfactor = oversamplingfactor, options...)[1]
    verbose && println("Discretization: oversampled grid of length $(length(grid))")
    grid_discretization(dict, fun, grid)
end

function grid_discretization(dict, fun, grid)
    A = evaluation_operator(dict, grid)
    B = sample(grid, fun, coefficient_type(dict))
    A, B
end


solve(style::SolverStyle, dict, A, B; options...) = solver(style, dict, A; options...) * B

solver(::InverseStyle, dict, A; options...) = inv(A)
solver(::DirectStyle, dict, A; options...) = directsolver(A; options...)
solver(::AZStyle, dict, A; Zt = default_Zt(dict, A), options...) = AZSolver(A, Zt; options...)
# TODO: Select a transform based operator here
solver(::TransformStyle, dict, A; options...) = inv(A)


solver(::TridiagonalProlateStyle, dict, A; scaling = scaling_factor(dict, A), options...) =
    FE_TridiagonalSolver(A, scaling; options...)
