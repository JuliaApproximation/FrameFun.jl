
## Functions to provide oversampled grids for FrameFun approximations.

function oversampledgrid(b::Dictionary, M)
    if has_interpolationgrid(b)
        interpolation_grid(resize(b, M))
    else
        error("Don't know how to compute an oversampled grid")
    end
end

oversampledgrid(b::BasisFunctions.DerivedDict, M) =
    oversampledgrid(superdict(b), M)

oversampledgrid(b::MappedDict, M) =
    mapped_grid(oversampledgrid(superdict(b), M), mapping(b))

oversampledgrid(dict::ExtensionFrame, M) =
    oversampledgrid(support(dict), basis(dict), M)

function oversampledgrid(domain::AbstractInterval, basis::FourierBasis, M)
    area = supremum(domain)-infimum(domain)
    L = ceil(Int,  M/area)-1
    subgrid(oversampledgrid(basis, L), domain)
end

function oversampledgrid(domain::Domain, basis::Dictionary, M)
    grid1 = oversampledgrid(basis, M)
    grid2 = subgrid(grid1, domain)

    # Initial guess : This could be way off if the original size was small.
    ratio = max(1,length(grid2)) / length(grid1)
    L = ceil(Int, M/ratio)
    grid3 = oversampledgrid(basis, L)
    grid4 = subgrid(grid3, domain)

    # If the number of sampling points is correct, return
    if length(grid4) == M
        return grid4
    end
    # Else let L grow
    Lmin = L
    Lmax = L
    while length(grid4) < M
        Lmin = L
        L = 2L
        grid3 = oversampledgrid(basis, L)
        grid4 = subgrid(grid3, domain)
        Lmax = L
    end

    # Now apply bisection to find the right L
    iterations = 0
    while (Lmax-Lmin) > 1 && iterations < 40
        Lmid = (Lmin+Lmax) >> 1
        grid3 = oversampledgrid(basis, Lmid)
        grid4 = subgrid(grid3, domain)
        length(grid4) < M ? Lmin=Lmid : Lmax=Lmid
        iterations += 1
    end
    grid3 = oversampledgrid(basis, Lmax)
    grid4 = subgrid(grid3, domain)
    grid4
end
