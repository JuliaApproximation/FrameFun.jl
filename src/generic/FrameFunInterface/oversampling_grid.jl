using DomainSets
import BasisFunctions: resize
## Functions to provide oversampled grids for FrameFun approximations.

resize(grid::ProductGrid, n) = ProductGrid(map(resize, elements(grid), n)...)

equispacedgrid(domain::AbstractInterval, L::Int) =
    EquispacedGrid(L, infimum(domain), supremum(domain))

equispacedgrid(domain::ProductDomain, L::Tuple) =
    cartesianproduct(map(equispacedgrid, elements(domain), L))

equispacedgrid(domain::Domain, L::Tuple) = subgrid(equispacedgrid(boundingbox(domain), L), domain)

function oversampling_grid(dict::Dictionary, L)
    if hasinterpolationgrid(dict)
        try
            return resize(interpolation_grid(dict), L)
        catch
            return interpolation_grid(resize(dict, L))
        end
    else
        # If the dictionary does not have an associated grid, we compute an
        # equispaced grid on its support.
        equispacedgrid(support(dict), L)
    end
end


oversampling_grid(dict::TensorProductDict, L) = ProductGrid(map(oversampling_grid, elements(dict), L)...)
oversampling_grid(dict::TensorProductDict, L::CartesianIndex) = oversampling_grid(dict, L.I)
oversampling_grid(dict::TensorProductDict, L::Int) = error("Expects a tuple or a CartesianIndex")

oversampling_grid(dict::BasisFunctions.CompositeDict, L) = oversampling_grid(element(dict,1), L)

oversampling_grid(dict::WeightedDict, L) = oversampling_grid(superdict(dict), L)

oversampling_grid(dict::MappedDict, L) = apply_map(oversampling_grid(superdict(dict), L), mapping(dict))

oversampling_grid(dict::OperatedDict, L) = oversampling_grid(superdict(dict), L)

# Make an initial guess for the sampling parameter L. The problem is that we have to
# determine the type of L. We try to infer this from the dimension of the dictionary and/or
# of its support.
initialguess(ap, M) = initialguess(dictionary(ap), M)
initialguess(dict::Dictionary1d, M::Int) = M
initialguess(dict::BasisFunctions.CompositeDict, M) = initialguess(element(dict,1), M)
initialguess(dict::BasisFunctions.CompositeDict, M::Int) = initialguess(element(dict,1), M)
initialguess(dict::Dictionary, M) = _initialguess(dict, M, support(dict), size(dict))
_initialguess(dict::Dictionary, M, domain::Domain, size) = size
_initialguess(dict::Dictionary, M, domain::Domain1d, size::Tuple{Int}) = M
function _initialguess(dict::Dictionary, M, domain::Domain2d, size::Tuple{Int})
    n = round(Int,sqrt(size[1]))
    (n,n)
end

# TODO Should everything be implemented for Dictionaries and Platforms or do we shift to platforms only
match_and_correct_sampling_parameter(strategy::SamplingStrategy, platform::Platform, param, M; dict=dictionary(platform, param), options...) =
    match_and_correct_sampling_parameter(strategy, platform, param, M, initialguess(dict, M); dict=dict, options...)
match_and_correct_sampling_parameter(::SamplingStrategy, dict::Dictionary, M; options...) =
    match_and_correct_sampling_parameter(dict, M, initialguess(dict, M); options...)

function match_and_correct_sampling_parameter(strategy::SamplingStrategy, platform::Platform, param, M, L_init; dict=dictionary(platform, param), options...)
    L_trial = match_sampling_parameter(dict, M, L_init)
    # The above does not always return a suitable L since it does not account for the platform/Dictionary
    # correct the result to a suitable one if necesarry
    correct_sampling_parameter(strategy, platform, param, L_trial; dict=dict, options...)
end

function match_and_correct_sampling_parameter(dict::Dictionary, M, L_init; options...)
    L_trial = match_sampling_parameter(dict, M, L_init)
    # The above does not always return a suitable L since it does not account for the platform/Dictionary
    # correct the result to a suitable one if necesarry
    correct_sampling_parameter(dict, L_trial; options...)
end

correct_sampling_parameter(dict::Dictionary, L_trial; options...) = L_trial
correct_sampling_parameter(strategy::SamplingStrategy, platform::Platform, param, L_trial; options...) = L_trial
correct_sampling_parameter(p::ProductPlatform, param, L; options...) =
    tuple(map(x->correct_sampling_parameter(x...; options...), zip(elements(p), param, L))...)

function match_sampling_parameter(samplingobject::Dictionary, M, L_init)
    objective(L) = length(oversampling_grid(samplingobject, L))-M
    monotonic_bisect(objective, L_init)
end

struct BisectionType end
Base.broadcastable(bt::BisectionType) = Ref(bt)

bisect_smaller(::BisectionType, L::Int) = L>>1
bisect_smaller(::BisectionType, L::Number) = L/2
bisect_larger(::BisectionType, L::Number) = 2L

bisect_average(bt::BisectionType, L1, L2) = bisect_average(bt, promote(L1,L2)...)
bisect_average(::BisectionType, L1::Int, L2::Int) = ceil(Int, (L1+L2)/2)
bisect_average(::BisectionType, L1::T, L2::T) where {T} = convert(T, (L1+L2)/2)

bisect_smaller(bt::BisectionType, L::Tuple) = bisect_smaller.(bt, L)
bisect_larger(bt::BisectionType, L::Tuple) = bisect_larger.(bt, L)
bisect_average(bt::BisectionType, L1::Tuple, L2::Tuple) = bisect_average.(bt, L1, L2)

bisect_small(l1) = bisect_smaller(BisectionType(), l1)
bisect_larger(l1) = bisect_larger(BisectionType(), l2)
bisect_average(l1, l2) = bisect_average(BisectionType(), l1, l2)


"""
Find the root of a monotonically increasing function with a given starting value
`L_init` and using bisection.

The value is increased or decreased first (using `bisect_larger` and
`bisect_smaller`) until an interval is found for which the objective function
switches sign. Next, the smallest value of `L` for which the objective function
becomes positive is located using bisection (with `bisect_average` to average out
two values of `L`).

The user can optionally supply an argument `bt`, which will be passed as the
first argument to the `bisect_X` routines. Thus, the variable `L` can have any
type as long as making it larger/smaller or computing averages can be meaningfully
defined.

The bisection is run until convergence. This ensures that if `L` has a discrete
nature, the optimal solution is always found exactly, in finite time.
"""
function monotonic_bisect(objective, L_init; maxiterations = 20, bt = BisectionType())
    local L1
    Z = objective(L_init)
    if Z ≈ 0
        return L_init
    end
    iterations = 0
    if Z > 0
        while (Z > 0) && (iterations < maxiterations)
            iterations += 1
            L1 = L_init
            L_init = bisect_smaller(bt, L1)
            Z = objective(L_init)
        end
        iterations == maxiterations && error("Maximal number of iterations reached in monotonic_bisect.")
        monotonic_bisect(objective, L_init, L1, bt)
    else
        while (Z < 0) && (iterations <  maxiterations)
            iterations += 1
            L1 = L_init
            L_init = bisect_larger(bt, L1)
            Z = objective(L_init)
        end
        iterations == maxiterations && error("Maximal number of iterations reached in monotonic_bisect.")
        monotonic_bisect(objective, L1, L_init, bt)
    end
end

function monotonic_bisect(objective, Lmin, Lmax, bt)
    Lmin, Lmax, bt
    Lmid = bisect_average(bt, Lmin, Lmax)
    if (Lmid == Lmax) || (Lmid == Lmin)
        return Lmid
    end
    Z = objective(Lmid)
    if Z ≈ 0
        return Lmid
    elseif Z > 0
        monotonic_bisect(objective, Lmin, Lmid, bt)
    else
        monotonic_bisect(objective, Lmid, Lmax, bt)
    end
end
