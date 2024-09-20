
import BasisFunctions:
    plotgrid,
    plot_postprocess

# Postprocessing when the set is a frame: set values outside the domain to NaN
function plot_postprocess(D::Domain, grid, vals, value=NaN)
    mgrid = GridArrays.MaskedGrid(grid, D)
    vals[.!(mask(mgrid))] .= value
    vals
end
plot_postprocess(B::ExtensionFrame, args...) = plot_postprocess(support(B), args...)
# Plotgrids are determined by the underlying set
plotgrid(B::ExtensionFrame, n) = plotgrid(basis(B),n)
