"""
A grid that contains the points of `omega_grid` that are not evaluated to zero by the elements that overlap with boundary_grid.
"""
boundary_support_grid(basis, boundary_grid::Union{MaskedGrid,IndexSubGrid}, omega_grid::Union{MaskedGrid,IndexSubGrid}) =
    boundary_support_grid(basis, basis, boundary_grid, omega_grid)

function boundary_support_grid(basis::Dictionary, dual::Dictionary, boundary_grid::Union{MaskedGrid,IndexSubGrid}, omega_grid::Union{MaskedGrid,IndexSubGrid})
    boundary_element_m = BasisFunctions.coefficient_index_mask_of_overlapping_elements(dual, boundary_grid)
    gamma = supergrid(omega_grid)
    m = BasisFunctions.grid_index_mask_in_element_support(basis, gamma, boundary_element_m)
    m .= m .& FrameFun.mask(omega_grid)
    MaskedGrid(gamma,m)
end


spline_util_restriction_operators(platform::BasisFunctions.GenericPlatform, i) =
    spline_util_restriction_operators(primal(platform, i), sampler(platform, i))

spline_util_restriction_operators(dict::Dictionary, sampler::GridSamplingOperator) =
    spline_util_restriction_operators(dict, grid(sampler))

spline_util_restriction_operators(frame::ExtensionFrame, grid::Union{IndexSubGrid,MaskedGrid}) =
    spline_util_restriction_operators(superdict(frame), grid, supergrid(grid), Domains.domain(frame))

spline_util_restriction_operators(dict::TensorProductDict{N,NTuple{N1,DT},S,T}, grid::AbstractGrid) where {N,N1,DT<:ExtensionFrame,S,T} =
    spline_util_restriction_operators(FrameFun.flatten(dict), grid)

spline_util_restriction_operators(dict::Dictionary, omega::AbstractGrid, grid::AbstractGrid, domain::Domain) =
    spline_util_restriction_operators(dict, omega, boundary_grid(grid, domain))

"""
Frame restriction operator and grid restriction operator.

The former restricts `dict` to the elements that overlap with the boundary and
the latter restricts `omega` to the points in the span of the dict elements that
overlap with the boundary.
"""
spline_util_restriction_operators(dict::Dictionary, omega::AbstractGrid, boundary::AbstractGrid) =
    _spline_util_restriction_operators(dict, omega, boundary_support_grid(dict, boundary, omega))


function BasisFunctions.grid_restriction_operator(src::Dictionary, dest::Dictionary, src_grid::Union{IndexSubGrid,MaskedGrid}, dest_grid::MaskedGrid)
    @assert supergrid(src_grid) == supergrid(dest_grid)
    IndexRestrictionOperator(src, dest, FrameFun.relative_indices(dest_grid,src_grid))
end

function _spline_util_restriction_operators(dict::Dictionary, grid::AbstractGrid, DMZ::AbstractGrid)
    boundary_indices = BasisFunctions.coefficient_indices_of_overlapping_elements(dict, DMZ)
    frame_restriction = IndexRestrictionOperator(dict, dict[boundary_indices], boundary_indices)
    grid_restriction = BasisFunctions.restriction_operator(gridbasis(grid), gridbasis(DMZ))
    frame_restriction, grid_restriction
end

"""
Frame restriction operator and grid restriction operator.

The former restricts `dict` to the elements that overlap with the boundary and
the latter restricts `grid` to the points in the span of the dict elements that
overlap with a region defined as the span of the dict elements that overlap with the boundary.
"""
function _spline_util_restriction_operators(dict::Dictionary, grid::AbstractGrid, boundary::AbstractGrid, DMZ::AbstractGrid)
    boundary_indices = BasisFunctions.coefficient_indices_of_overlapping_elements(dict, boundary)
    frame_restriction = IndexRestrictionOperator(dict, dict[boundary_indices], boundary_indices)
    grid_restriction = BasisFunctions.restriction_operator(gridbasis(grid), gridbasis(DMZ))
    frame_restriction, grid_restriction
end

estimate_plunge_rank(src::BSplineTranslatesBasis, domain::Domain, dest::GridBasis) =
    estimate_plunge_rank(src, domain, grid(dest))
estimate_plunge_rank(src::BSplineTranslatesBasis, domain::Domain, grid::Union{IndexSubGrid,MaskedGrid}) =
    estimate_plunge_rank_bspline(src, domain, supergrid(grid))

estimate_plunge_rank(src::TensorProductDict{N,DT,S,T}, domain::Domain, dest::GridBasis) where {N,DT<:NTuple{N1,BasisFunctions.BSplineTranslatesBasis} where {N1},S,T} =
    estimate_plunge_rank(src, domain, grid(dest))

estimate_plunge_rank(src::TensorProductDict{N,DT,S,T}, domain::Domain, grid::MaskedGrid) where {N,DT<:NTuple{N1,BasisFunctions.BSplineTranslatesBasis} where {N1},S,T} =
    estimate_plunge_rank_bspline(src, domain, supergrid(grid))

function estimate_plunge_rank_bspline(src, domain::Domain, grid::AbstractGrid)
    boundary = boundary_grid(grid, domain)
    sum(BasisFunctions.coefficient_index_mask_of_overlapping_elements(src, boundary))
end
