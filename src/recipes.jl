## Series


## Plot full FrameFuns as the underlying expansion
@recipe function f(F::FrameFun; plot_ext = false)
    plot_ext || title --> "FrameFun"
    plot_ext ? SetExpansion(basis(F),coefficients(F)) : expansion(F)
end

# When supplying a function along with the FrameFun, plot the error
@recipe function f(F::FrameFun, target::Function; plot_ext = true)
    plot_ext ? (SetExpansion(basis(F),coefficients(F)), target) : (expansion(F), target)
end

@recipe function f(dom::AbstractDomain2d; n=300)
    seriescolor --> :blues
    seriestype --> :heatmap
    B = boundingbox(dom)
    grid = equispaced_aspect_grid(B,n)
    Z = evalgrid(grid, dom)
    grid, 1./Z
end

# Postprocessing when the set is a frame: set values outside the domain to NaN
function postprocess(B::DomainFrame, grid, vals)
    mgrid = subgrid(grid, domain(B))
    for i in eachindex(grid)
        if ~in(i, mgrid)
            vals[i] = NaN
        end
    end
    vals
end
        
# Plotgrids are determined by the underlying set
plotgrid(B::DomainFrame, n) = plotgrid(basis(B),n)

# Plot a domain
