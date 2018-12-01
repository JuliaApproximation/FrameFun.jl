
###################
# Helper functions
###################

determine_return_type(fun, ::Type{S}) where {S} = Base.Core.Compiler.return_type(fun, (S,))
determine_return_type(fun, S::Type{<:Tuple}) = Base.Core.Compiler.return_type(fun, S)

function promote_dictionary(dict, fun)
    T = promote_type(codomaintype(dict), determine_return_type(fun, domaintype(dict)))
    promote_coefficienttype(dict, T)
end

interpolation_grid(dict::Dictionary; options...) = grid(dict)
interpolation_grid(platform::Platform, param, dict::Dictionary; options...) = interpolation_grid(dict; options...)
oversampled_grid(platform::Platform, param, dict::Dictionary; options...) = oversampled_grid(dict; options...)


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

DiscretizationStyle(dict::Dictionary) = is_basis(dict) ? InterpolationStyle() : OversamplingStyle()

SolverStyle(dict::Dictionary, dstyle::InterpolationStyle) = has_transform(dict) ? TransformStyle() : DirectStyle()
SolverStyle(dict::Dictionary, dstyle::OversamplingStyle) = AZStyle()



##########################
# Approximation problems
##########################

abstract type ApproximationProblem end

struct DictionaryApproximation <: ApproximationProblem
    dict    ::  Dictionary
end

Dictionary(ap::DictionaryApproximation) = ap.dict

DiscretizationStyle(ap::DictionaryApproximation) = DiscretizationStyle(Dictionary(ap))
SolverStyle(ap::DictionaryApproximation, dstyle::DiscretizationStyle) = SolverStyle(Dictionary(ap), dstyle)

for op in (:interpolation_grid, :oversampled_grid, :approximation_grid)
    @eval $op(ap::DictionaryApproximation; options...) = $op(Dictionary(ap); options...)
end

# TODO: clean up these scaling factors. They have to do with a normalization of the
# sampling operators.
Zt_scaling_factor(S::Dictionary, A) = length(supergrid(grid(dest(A))))
Zt_scaling_factor(S::DerivedDict, A) = Zt_scaling_factor(superdict(S), A)
Zt_scaling_factor(S::ChebyshevBasis, A) = length(supergrid(grid(dest(A))))/2

AZ_Zt(dict::Dictionary, A) = (one(eltype(A))/convert(eltype(A), Zt_scaling_factor(dict, A))) * A'

AZ_Zt(ap::DictionaryApproximation, A) = AZ_Zt(Dictionary(ap), A)



struct PlatformApproximation <: ApproximationProblem
    platform    ::  Platform
    param
    dict        ::  Dictionary
end

# Compute the dictionary if it was not provided
PlatformApproximation(platform, param) = PlatformApproximation(platform, param, Dictionary(platform, param))

Dictionary(ap::PlatformApproximation) = ap.dict

DiscretizationStyle(ap::PlatformApproximation) = DiscretizationStyle(ap.platform)
SolverStyle(dict::PlatformApproximation, dstyle::DiscretizationStyle) = SolverStyle(ap.platform, dstyle)

for op in (:interpolation_grid, :oversampled_grid, :approximation_grid)
    @eval $op(ap::PlatformApproximation; options...) = $op(ap.platform, ap.param, ap.dict; options...)
end

AZ_Zt(ap::PlatformApproximation, A) = AZ_Zt(Dictionary(ap), A)

Base.getindex(platform::Platform, param) = PlatformApproximation(platform, param)


approximationproblem(dict::Dictionary) = DictionaryApproximation(dict)
approximationproblem(platform, param) = PlatformApproximation(platform, param)


####################
# The Fun interface
####################

# The `Fun` interface first turns the approximation problem into a subtype of
# ApproximationProblem. Next, it calls `approximate`. Finally, it discards all
# the operators that were computed and simply returns the function.
# Users wanting to access the operators can call `approximate` directly.

function Fun(dict::Dictionary, domain::Domain, fun;
        coefficienttype = promote_type(codomaintype(dict), determine_return_type(fun, domaintype(dict))),
        options...)
    promoted_dict = promote_coefficienttype(dict, coefficienttype)
    Fun(ExtensionFrame(domain, promoted_dict), fun; options...)
end

function Fun(dict::Dictionary, fun;
            coefficienttype = promote_type(codomaintype(dict), determine_return_type(fun, domaintype(dict))),
            options...)
    promoted_dict = promote_coefficienttype(dict, coefficienttype)
    ap = DictionaryApproximation(promoted_dict)
    Fun(ap, fun; options...)
end

function Fun(platform::Platform, param, fun; verbose = false, options...)
    ap = platform[param]
    verbose && println("Platform: using dictionary $(Dictionary(ap))")
    Fun(ap, fun; verbose = verbose, options...)
end

function Fun(ap::ApproximationProblem, fun; options...)
    A, B, C, F = approximate(ap, fun; options...)
    F
end

# The difference between Fun and approximate is that approximate returns all the
# operators it constructed.
function approximate(ap::ApproximationProblem, fun;
            discretizationstyle = DiscretizationStyle(ap),
            solverstyle = SolverStyle(ap, discretizationstyle),
            verbose = false,
            options...)

    if verbose
        println("Fun: discretizing with style $discretizationstyle")
        println("Fun: solving with style $solverstyle")
    end
    approximate(discretizationstyle, solverstyle, ap, fun; verbose=verbose, options...)
end

# Construct the approximation problem and solve it
function approximate(discretizationstyle::DiscretizationStyle, solverstyle::SolverStyle,
            ap::ApproximationProblem, fun; options...)
    A, B = discretization(discretizationstyle, ap, fun; options...)
    C = solve(solverstyle, ap, A, B; options...)
    A, B, C, DictFun(Dictionary(ap), C)
end

# Convenience function that one can call directly in order to obtain the discretization
# without going through Fun (which does not return the operators)
discretization(dict::Dictionary, fun; discretizationstyle = DiscretizationStyle(dict), options...) =
    discretization(discretizationstyle, DictionaryApproximation(dict), fun; options...)

discretization(::InterpolationStyle, ap::ApproximationProblem, fun; options...) =
    grid_discretization(Dictionary(ap), fun, interpolation_grid(ap; options...))

discretization(::GridStyle, ap::ApproximationProblem, fun;
            grid = approximation_grid(ap), options...) =
    grid_discretization(Dictionary(ap), fun, grid)


function discretization(::OversamplingStyle, ap::ApproximationProblem, fun;
            oversamplingfactor = 2, verbose = false, options...)
    grid = oversampled_grid(ap; oversamplingfactor = oversamplingfactor, options...)[1]
    verbose && println("Discretization: oversampled grid of length $(length(grid))")
    discretization(GridStyle(), ap, fun; grid = grid, verbose = verbose, options...)
end

function grid_discretization(dict, fun, grid)
    A = evaluation_operator(dict, grid)
    B = sample(grid, fun, coefficienttype(dict))
    A, B
end


solve(style::SolverStyle, ap, A, B; options...) = solver(style, ap, A; options...) * B

solver(::InverseStyle, ap, A; options...) = inv(A)
solver(::DirectStyle, ap, A; options...) = directsolver(A; options...)
solver(::AZStyle, ap, A; Zt = AZ_Zt(ap, A), options...) = AZSolver(A, Zt; options...)
solver(::AZSmoothStyle, ap, A; Zt = AZ_Zt(ap, A), options...) =
    AZSmoothSolver(A, Zt; options...)

# TODO: Select a transform based operator here
solver(::TransformStyle, ap, A; options...) = inv(A)


solver(::TridiagonalProlateStyle, ap, A; scaling = Zt_scaling_factor(Dictionary(ap), A), options...) =
    FE_TridiagonalSolver(A, scaling; options...)
