## Series
import BasisFunctions: plotgrid, postprocess
using GridArrays: MaskedGrid
# Postprocessing when the set is a frame: set values outside the domain to NaN
function postprocess(D::Domain, grid, vals, value=NaN)
    mgrid = MaskedGrid(grid, D)
    vals[.!(mask(mgrid))] .= value
    vals
end
postprocess(B::ExtensionFrame, args...) = postprocess(support(B), args...)
# Plotgrids are determined by the underlying set
plotgrid(B::ExtensionFrame, n) = plotgrid(basis(B),n)
