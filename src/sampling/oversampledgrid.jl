
## Functions to provide oversampled grids for FrameFun approximations.
## TODO: make more generic, right now this relies on resizing bases and only works with the default grid

# Sidenote: very slow when testing with no inlining.

oversampled_grid(set::ExtensionFrame; options...) =
    oversampled_grid(domain(set), basis(set); options...)

function oversampled_grid(domain, basis::Dictionary; oversamplingfactor, options...)
    N = dimension(basis)
    n_goal = length(basis) * oversamplingfactor^N
    grid1 = BasisFunctions.grid(basis)
    grid2 = FrameFun.subgrid(grid1, domain)
    ratio = max(1,length(grid2)) / length(grid1)
    # Initial guess : This could be way off if the original size was small.
    newsize = ceil(Int,n_goal/ratio)
    n = BasisFunctions.approx_length(basis, newsize)
    large_basis = resize(basis, n)
    grid3 = BasisFunctions.grid(large_basis)
    grid4 = FrameFun.subgrid(grid3, domain)
    # If the number of sampling points is correct, return
    if length(grid4)==n_goal
        return grid4, large_basis
    end
    maxN = newsize
    #
    while length(grid4)<n_goal
        newsize = 2*newsize
        n = BasisFunctions.approx_length(basis, newsize)
        large_basis = resize(basis, n)
        grid3 = BasisFunctions.grid(large_basis)
        grid4 = FrameFun.subgrid(grid3, domain)
        maxN = newsize
    end
    minN = newsize>>>1
    its = 0
    while (maxN-minN) >1 && its < 40
        midpoint = (minN+maxN) >>> 1
        n = BasisFunctions.approx_length(basis,  midpoint)
        large_basis = resize(basis, n)
        grid3 = BasisFunctions.grid(large_basis)
        grid4 = FrameFun.subgrid(grid3, domain)
        length(grid4)<n_goal ? minN=midpoint : maxN=midpoint
        its += 1
    end
    n = BasisFunctions.approx_length(basis,  maxN)
    large_basis = resize(basis, n)
    grid3 = BasisFunctions.grid(large_basis)
    grid4 = FrameFun.subgrid(grid3, domain)
    grid4, large_basis
end
