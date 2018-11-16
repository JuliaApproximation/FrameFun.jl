## Series


## Plot full DictFun as the underlying expansion
@recipe function f(F::DictFun; plot_ext = false)
    # plot_ext || (title --> "DictFun")
    plot_ext ? Expansion(basis(F),coefficients(F)) : expansion(F)
end

# When supplying a function along with the DictFun, plot the error
@recipe function f(F::DictFun, target::Function; plot_ext = true)
    plot_ext ? (Expansion(basis(F),coefficients(F)), target) : (expansion(F), target)
end

# Postprocessing when the set is a frame: set values outside the domain to NaN
function postprocess(D::Domain, grid, vals, value=NaN)
    mgrid = subgrid(grid, D)
    for i in eachindex(grid)
        if ~is_subindex(i, mgrid)
            vals[i] = value
        end
    end
    vals
end

postprocess(B::ExtensionFrame, args...) = postprocess(domain(B), args...)

# Plotgrids are determined by the underlying set
plotgrid(B::ExtensionFrame, n) = plotgrid(basis(B),n)

# Plot a domain

@recipe function f(dom::Domain2d; n=300, distance=false, xlim=[-1,1], ylim=[-1,1])
    seriescolor --> :tempo
     seriestype --> :heatmap
    aspect_ratio --> 1
    cbar --> false
    # xrange = linspace(xlim[1],xlim[2],n)
    xrange = EquispacedGrid(n,xlim[1],xlim[2])
    # yrange = linspace(ylim[1],ylim[2],n)
    yrange = EquispacedGrid(n,ylim[1],ylim[2])
    # grid = [SVector(i,j) for i in xrange , j in yrange]
    grid = xrange×yrange

    plotdata = distance ? dist.(grid, Ref(dom)) : in.(grid, Ref(dom))
    collect(xrange),collect(yrange),plotdata
end

# @recipe function f(dom::Domain2d; n=300)
#     seriescolor --> :blues
#     seriestype --> :heatmap
#     B = boundingbox(dom)
#     grid = equispaced_aspect_grid(B,n)
#     Z = indomain_broadcast(grid, dom)
#     grid, 1./Z
# end
