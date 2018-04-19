
# Some of the util functions situated here can be located to more suitable locations.

"""
Index of elements of `B` that overlap with `boundary`.
"""
function boundary_element_indices(B, boundary::MaskedGrid)
    s = Set{Int}()
    for x in boundary
        push!(s,overlapping_elements(B,x)...)
    end
    collect(s)
end


"""
A grid that contains the points of `omega_grid` that are not evaluated to zero by the elements that overlap with boundary_grid.
"""
function boundary_support_grid(B, boundary_grid::MaskedGrid, omega_grid::MaskedGrid)
    boundary_indices = boundary_element_indices(B,boundary_grid)
    s = Set{Int}()
    for i in boundary_indices
        push!(s,BasisFunctions.support_indices(B,supergrid(omega_grid),i)...)
    end
    a = collect(s)
    mask = zeros(Bool,size(omega_grid.mask))
    mask[a] = true
    mask .= mask .& omega_grid.mask
    MaskedGrid(supergrid(omega_grid),mask)
end

az_selection_util_operators(platform::BasisFunctions.GenericPlatform, i) =
    az_selection_util_operators(primal(platform, i), sampler(platform, i))

az_selection_util_operators(frame::ExtensionFrame, sampler::GridSamplingOperator) =
    az_selection_util_operators(superdict(frame), BasisFunctions.grid(sampler), supergrid(BasisFunctions.grid(sampler)), Domains.domain(frame))

function az_selection_util_operators(dict::Dictionary, omega::AbstractGrid, grid::AbstractGrid, domain::Domain)
    boundary = boundary_grid(grid, domain)
    boundary_indices = boundary_element_indices(dict, boundary)
    boundary_support = boundary_support_grid(dict, boundary, omega)
    IndexRestrictionOperator(Span(dict), Span(dict[boundary_indices]), boundary_indices), restriction_operator(omega, boundary_support)
end

estimate_plunge_rank(src::BSplineTranslatesBasis, domain::Domain, dest::GridBasis) =
    estimate_plunge_rank(src, domain, grid(dest))
estimate_plunge_rank(src::BSplineTranslatesBasis, domain::Domain, grid::MaskedGrid) =
    estimate_plunge_rank_bspline(src, domain, supergrid(grid))

estimate_plunge_rank(src::TensorProductDict{N,DT,S,T}, domain::Domain, dest::GridBasis) where {N,DT<:NTuple{N1,BasisFunctions.BSplineTranslatesBasis} where {N1},S,T} =
    estimate_plunge_rank(src, domain, grid(dest))

estimate_plunge_rank(src::TensorProductDict{N,DT,S,T}, domain::Domain, grid::MaskedGrid) where {N,DT<:NTuple{N1,BasisFunctions.BSplineTranslatesBasis} where {N1},S,T} =
    estimate_plunge_rank_bspline(src, domain, supergrid(grid))

function estimate_plunge_rank_bspline(src, domain::Domain, grid::AbstractGrid)
    boundary = boundary_grid(grid, domain)
    length(boundary_element_indices(src, boundary))
end
