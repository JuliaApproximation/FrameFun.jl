
## Functions to provide oversampled grids for FrameFun approximations.

equispacedgrid(domain::AbstractInterval, smpl_par::Int) =
    EquispacedGrid(smpl_par, infimum(domain), supremum(domain))

equispacedgrid(domain::ProductDomain, smpl_par::Tuple) =
    cartesianproduct(map(equispacedgrid, components(domain), smpl_par))

equispacedgrid(domain::Domain, smpl_par::Tuple) = subgrid(equispacedgrid(boundingbox(domain), smpl_par), domain)

function oversampling_grid(dict::Dictionary, smpl_par)
    if hasinterpolationgrid(dict)
        try
            return resize(interpolation_grid(dict), smpl_par)
        catch
            return interpolation_grid(resize(dict, smpl_par))
        end
    else
        # If the dictionary does not have an associated grid, we compute an
        # equispaced grid on its support.
        equispacedgrid(support(dict), smpl_par)
    end
end


oversampling_grid(dict::TensorProductDict, smpl_par) = ProductGrid(map(oversampling_grid, components(dict), smpl_par)...)
oversampling_grid(dict::TensorProductDict, smpl_par::CartesianIndex) = oversampling_grid(dict, smpl_par.I)
oversampling_grid(dict::TensorProductDict, smpl_par::Int) = error("Expects a tuple or a CartesianIndex")

oversampling_grid(dict::BasisFunctions.CompositeDict, smpl_par) = oversampling_grid(component(dict,1), smpl_par)

oversampling_grid(dict::WeightedDict, smpl_par) = oversampling_grid(superdict(dict), smpl_par)

oversampling_grid(dict::MappedDict, smpl_par) = apply_map(oversampling_grid(superdict(dict), smpl_par), forward_map(dict))

oversampling_grid(dict::OperatedDict, smpl_par) = oversampling_grid(superdict(dict), smpl_par)

# Make an initial guess for the sampling parameter smpl_par. The problem is that we have to
# determine the type of smpl_par. We try to infer this from the dimension of the dictionary and/or
# of its support.
initialguess(ap, M) = initialguess(dictionary(ap), M)
initialguess(dict::Dictionary1d, M::Int) = M
initialguess(dict::BasisFunctions.CompositeDict, M) = initialguess(component(dict,1), M)
initialguess(dict::BasisFunctions.CompositeDict, M::Int) = initialguess(component(dict,1), M)
initialguess(dict::Dictionary, M) = _initialguess(dict, M, support(dict), size(dict))
_initialguess(dict::Dictionary, M, domain::Domain, size) = size
_initialguess(dict::Dictionary, M, domain::Domain1d, size::Tuple{Int}) = M
function _initialguess(dict::Dictionary, M, domain::Domain2d, size::Tuple{Int})
    n = round(Int,sqrt(size[1]))
    (n,n)
end

match_and_correct_sampling_parameter(platform::Platform, plt_par, M; dict=dictionary(platform, plt_par), options...) =
    match_and_correct_sampling_parameter(platform, plt_par, M, initialguess(dict, M); dict=dict, options...)
match_and_correct_sampling_parameter(dict::Dictionary, M; options...) =
    match_and_correct_sampling_parameter(dict, M, initialguess(dict, M); options...)

function match_and_correct_sampling_parameter(platform::Platform, plt_par, M, smpl_par_init; dict=dictionary(platform, plt_par), options...)
    smpl_par_trial = match_sampling_parameter(dict, M, smpl_par_init)
    # The above does not always return a suitable smpl_par since it does not account for the platform/Dictionary
    # correct the result to a suitable one if necesarry
    correct_sampling_parameter(platform, plt_par, smpl_par_trial; dict=dict, options...)
end

function match_and_correct_sampling_parameter(dict::Dictionary, M, smpl_par_init; options...)
    smpl_par_trial = match_sampling_parameter(dict, M, smpl_par_init)
    # The above does not always return a suitable smpl_par since it does not account for the platform/Dictionary
    # correct the result to a suitable one if necesarry
    correct_sampling_parameter(dict, smpl_par_trial; options...)
end

correct_sampling_parameter(dict::Dictionary, smpl_par_trial; options...) = smpl_par_trial
correct_sampling_parameter(platform::Platform, plt_par, smpl_par_trial; options...) = smpl_par_trial
correct_sampling_parameter(p::ProductPlatform, plt_par, smpl_par; options...) =
    tuple(map(x->correct_sampling_parameter(x...; options...), zip(components(p), plt_par, smpl_par))...)

function match_sampling_parameter(samplingobject::Dictionary, M, smpl_par_init)
    objective(smpl_par) = length(oversampling_grid(samplingobject, smpl_par))-M
    monotonic_bisect(objective, smpl_par_init)
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
