## Series
import BasisFunctions: plotgrid, postprocess
# Postprocessing when the set is a frame: set values outside the domain to NaN
function postprocess(D::Domain, grid, vals, value=NaN)
    mgrid = subgrid(grid, D)
    for i in eachindex(grid)
        if ~issubindex(i, mgrid)
            vals[i] = value
        end
    end
    vals
end
postprocess(B::ExtensionFrame, args...) = postprocess(support(B), args...)
# Plotgrids are determined by the underlying set
plotgrid(B::ExtensionFrame, n) = plotgrid(basis(B),n)
