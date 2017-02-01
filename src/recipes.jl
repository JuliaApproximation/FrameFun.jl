## Series


## Plot full SetFun as the underlying expansion
@recipe function f(F::SetFun; plot_ext = false)
    plot_ext || (title --> "SetFun")
    plot_ext ? SetExpansion(basis(F),coefficients(F)) : expansion(F)
end

# When supplying a function along with the SetFun, plot the error
@recipe function f(F::SetFun, target::Function; plot_ext = true)
    plot_ext ? (SetExpansion(basis(F),coefficients(F)), target) : (expansion(F), target)
end

@recipe function f(dom::AbstractDomain2d; n=300)
    seriescolor --> :blues
    seriestype --> :heatmap
    B = boundingbox(dom)
    grid = equispaced_aspect_grid(B,n)
    Z = indomain_grid(grid, dom)
    grid, 1./Z
end

# Postprocessing when the set is a frame: set values outside the domain to NaN
function postprocess(D::AbstractDomain, grid, vals, value=NaN)
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
