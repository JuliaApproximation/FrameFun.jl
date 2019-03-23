
## Functions to provide oversampled grids for FrameFun approximations.

resize(grid::ProductGrid, n) = ProductGrid(map(resize, elements(grid), n)...)

equispacedgrid(domain::AbstractInterval, L::Int) =
    EquispacedGrid(L, infimum(domain), supremum(domain))

equispacedgrid(domain::ProductDomain, L::Tuple) =
    cartesianproduct(map(equispacedgrid, elements(domain), L))

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

oversampling_grid(dict::ExtensionFrame, L) = subgrid(oversampling_grid(superdict(dict),L), support(dict))

oversampling_grid(dict::WeightedDict, L) = oversampling_grid(superdict(dict), L)

oversampling_grid(dict::MappedDict, L) = apply_map(oversampling_grid(superdict(dict), L), mapping(dict))

oversampling_grid(dict::OperatedDict, L) = oversampling_grid(superdict(dict), L)

initialguess(dict::Dictionary1d, M::Int) = M
initialguess(dict::Dictionary, M) = size(dict)
initialguess(dict::BasisFunctions.CompositeDict, M) = initialguess(element(dict,1), M)
initialguess(ap, M) = initialguess(dictionary(ap), M)

match_sampling_parameter(dict, M::Int) = match_sampling_parameter(dict, M, initialguess(dict, M))

function match_sampling_parameter(samplingobject, M::Int, L0)
    objective(L) = length(oversampling_grid(samplingobject, L))-M
    monotonic_bisect(objective, L0)
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
`L0` and using bisection.

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
function monotonic_bisect(objective, L0; maxiterations = 20, bt = BisectionType())
    local L1
    Z = objective(L0)
    if Z ≈ 0
        return L0
    end
    iterations = 0
    if Z > 0
        while (Z > 0) && (iterations < maxiterations)
            iterations += 1
            L1 = L0
            L0 = bisect_smaller(bt, L1)
            Z = objective(L0)
        end
        iterations == maxiterations && error("Maximal number of iterations reached in monotonic_bisect.")
        monotonic_bisect(objective, L0, L1, bt)
    else
        while (Z < 0) && (iterations <  maxiterations)
            iterations += 1
            L1 = L0
            L0 = bisect_larger(bt, L1)
            Z = objective(L0)
        end
        iterations == maxiterations && error("Maximal number of iterations reached in monotonic_bisect.")
        monotonic_bisect(objective, L1, L0, bt)
    end
end

function monotonic_bisect(objective, Lmin, Lmax, bt)
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

function oversampled_grid(dict::Dictionary, M)
    L = match_sampling_parameter(dict, M)
    oversampling_grid(dict, L)
end
