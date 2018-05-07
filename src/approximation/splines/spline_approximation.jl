
# Some of the util functions situated here can be located to more suitable locations.

"""
Index of elements of `B` that overlap with `boundary`.
"""
function boundary_element_indices(B, boundary::AbstractGrid)
    s = Set{CartesianIndex{dimension(B)}}()
    for x in boundary
        push!(s,BasisFunctions.overlapping_elements(B,x)...)
    end
    collect(s)
end


"""
A grid that contains the points of `omega_grid` that are not evaluated to zero by the elements that overlap with boundary_grid.
"""
function boundary_support_grid(B, boundary_grid::Union{MaskedGrid}, omega_grid::Union{MaskedGrid,IndexSubGrid})
    boundary_indices = FrameFun.boundary_element_indices(B,boundary_grid)
    s = Set{CartesianIndex{dimension(B)}}()
    for i in boundary_indices
        push!(s,BasisFunctions.support_indices(B,supergrid(omega_grid),i)...)
    end
    a = collect(s)
    m = BitArray(size(FrameFun.mask(omega_grid)))#zeros(Bool,)
    m[:] = 0
    m[a] = true
    m .= m .& FrameFun.mask(omega_grid)
    MaskedGrid(supergrid(omega_grid),m)
end


function boundary_support_grid(B, boundary_grid::IndexSubGrid, omega_grid::MaskedGrid)
    boundary_indices = FrameFun.boundary_element_indices(B,boundary_grid)
    s = Set{CartesianIndex{dimension(B)}}()
    for i in boundary_indices
        push!(s,BasisFunctions.support_indices(B,supergrid(omega_grid),i)...)
    end
    a = collect(s)
    b = CartesianIndex{dimension(B)}[]
    for ai in a
        if is_subindex(ai, omega_grid)
            push!(b,ai)
        end
    end
    IndexSubGrid(supergrid(boundary_grid), b)
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


function BasisFunctions.grid_restriction_operator(src::Span, dest::Span, src_grid::Union{IndexSubGrid,MaskedGrid}, dest_grid::MaskedGrid)
    @assert supergrid(src_grid) == supergrid(dest_grid)
    IndexRestrictionOperator(src, dest, FrameFun.relative_indices(dest_grid,src_grid))
end

function _spline_util_restriction_operators(dict::Dictionary, grid::AbstractGrid, DMZ::AbstractGrid)
    boundary_indices = FrameFun.boundary_element_indices(dict, DMZ)
    frame_restriction = IndexRestrictionOperator(Span(dict), Span(dict[boundary_indices]), boundary_indices)
    grid_restriction = BasisFunctions.restriction_operator(gridspace(grid), gridspace(DMZ))
    frame_restriction, grid_restriction
end

"""
Frame restriction operator and grid restriction operator.

The former restricts `dict` to the elements that overlap with the boundary and
the latter restricts `grid` to the points in the span of the dict elements that
overlap with a region defined as the span of the dict elements that overlap with the boundary.
"""
function _spline_util_restriction_operators(dict::Dictionary, grid::AbstractGrid, boundary::AbstractGrid, DMZ::AbstractGrid)
    boundary_indices = FrameFun.boundary_element_indices(dict, boundary)
    frame_restriction = IndexRestrictionOperator(Span(dict), Span(dict[boundary_indices]), boundary_indices)
    grid_restriction = BasisFunctions.restriction_operator(gridspace(grid), gridspace(DMZ))
    frame_restriction, grid_restriction
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


function divide_and_conquer_restriction_operators(fplatform::BasisFunctions.Platform, i, dim, range)
    platform = fplatform.super_platform
    basis = primal(platform, i)
    S = sampler(fplatform, i)
    s = sampler(platform, i)
    domain = FrameFun.domain(primal(fplatform, i))
    gamma = BasisFunctions.grid(s)
    omega = BasisFunctions.grid(S)
    FrameFun.divide_and_conquer_restriction_operators(omega, gamma, basis, domain, dim, range)
end

# The grid on the boundary of omega
divide_and_conquer_restriction_operators(omega::AbstractGrid, gamma::AbstractGrid, basis::Dictionary, domain::Domains.Domain, dim::Int, range::AbstractVector) =
    divide_and_conquer_restriction_operators(omega::AbstractGrid, gamma::AbstractGrid, basis::Dictionary, FrameFun.boundary_support_grid(basis, boundary_grid(gamma, domain), omega)::AbstractGrid, dim::Int, range::AbstractVector)

function intersect_boundary(DMZ::AbstractGrid, gamma::AbstractGrid, dim::Int, range)
    a = leftendpoint(gamma); b = rightendpoint(gamma)
    dx = BasisFunctions.stepsize(elements(gamma)[dim])
    mid = mean(range)
    D = dimension(gamma)
    # Some domain restriction one dimension of the bounding box
    if dim ==1
        split_domain = interval(mid-dx/2,mid+dx/2)×Domains.ProductDomain([interval(a[i],b[i]) for i in 2:D]...)
    elseif dim==D
        split_domain = ProductDomain([interval(a[i],b[i]) for i in 1:D-1]...)×interval(mid-dx/2,mid+dx/2)
    else
        split_domain = ProductDomain([interval(a[i],b[i]) for i in 1:dim-1]...)×interval(mid-dx/2,mid+dx/2)×ProductDomain([interval(a[i],b[i]) for i in dim+1:D]...)
    end

    # The grid on the intersection of split_domain and the boundary
    split_grid = FrameFun.subgrid(DMZ, split_domain)
    if length(split_grid) == 0
        error("check your ranges for splitting the domain.")
    end
    split_grid
end

function divide_and_conquer_restriction_operators(omega::AbstractGrid, gamma::AbstractGrid, basis::Dictionary, DMZ::AbstractGrid, dim::Int, range::AbstractVector)
    split_grid = FrameFun.intersect_boundary(DMZ, gamma, dim, range)
    #  The points on the support of the functions overlapping with split_grid
    mid = FrameFun.boundary_support_grid(basis, split_grid, DMZ)
    # Look at the grid away from the intersection and split in two disconnected parts.
    left_over = DMZ-mid
    if length(left_over) != 0
        left, right = split(left_over)
        leftDMZ = FrameFun.boundary_support_grid(basis, left, DMZ)
        rightDMZ = FrameFun.boundary_support_grid(basis, right, DMZ)

        if !(.5 < length(left)/length(right) < 2)
            warn("ratio between the separated grids is not distributed evenly: $(length(left)) and $(length(right))")
        end

        BR0, DMZ_R0 = _spline_util_restriction_operators(basis, omega, mid)
        BR1, DMZ_R1 = _spline_util_restriction_operators(basis, omega, left, leftDMZ)
        BR2, DMZ_R2 = _spline_util_restriction_operators(basis, omega, right, rightDMZ)

        BR0, BR1, BR2, DMZ_R0, DMZ_R1, DMZ_R2
    else
        warn("No grid splitting possible")
        BR, DMZ_R = _spline_util_restriction_operators(basis, omega, DMZ)
    end
end


function split_DMZ(DMZ::AbstractGrid, basis, gamma, ranges; options...)
    split_grid = intersect_DMZ(DMZ::AbstractGrid, basis, gamma::AbstractGrid, ranges; options...)
    split_grid_DMZ = FrameFun.boundary_support_grid(basis, split_grid, DMZ)
    left_over = DMZ - split_grid_DMZ
    if length(left_over) == 0
        return [DMZ], AbstractGrid[]
    end
    FrameFun.split_in_IndexGrids(split_grid_DMZ), FrameFun.split_in_IndexGrids(left_over)
end


function intersect_DMZ(DMZ::AbstractGrid, basis, gamma::AbstractGrid, ranges; verbose=false, factor=3, options...)
    domain_grid = domain_grid_Nd(basis, ranges, factor; options...)
    split_domain = split_domain_Nd(gamma, domain_grid)
    verbose && println(split_domain)

    # The grid on the intersection of split_domain and the boundary
    split_grid = FrameFun.subgrid(DMZ, split_domain)
    if length(split_grid) == 0
        error("check your mids for splitting the domain.")
    end
    split_grid
end



domain_grid_Nd(basis, ranges, factor; options...) = [domain_grid_1d(d, basis, ranges, factor; options...) for d in 1:dimension(basis)]

function domain_grid_1d(d, basis, ranges, factor; shift=false)
    basis_coarseness = BasisFunctions.support_length_of_compact_function(element(basis, d))
    mid = mean(ranges[d])
    shift && (mid = mid + basis_coarseness*factor/2)
    e = mid; es = [e]
    e = mid - factor*basis_coarseness
    while e>minimum(ranges[d]) + .5basis_coarseness
        push!(es, e)
        e = e - factor*basis_coarseness
    end
    e = mid + factor*basis_coarseness
    while e<maximum(ranges[d]) - .5basis_coarseness
        push!(es, e)
        e = e + factor*basis_coarseness
    end

    tuple(es...)
end



function split_domain_1d(gamma, dim, mid)
    a = leftendpoint(gamma); b = rightendpoint(gamma)
    dx = BasisFunctions.stepsize(elements(gamma)[dim])
    D = length(dx)
    if dim ==1
        split_domain = interval(mid-dx/2,mid+dx/2)×Domains.ProductDomain([interval(a[i],b[i]) for i in 2:length(a)]...)
    elseif dim==D
        split_domain = ProductDomain([interval(a[i],b[i]) for i in 1:length(a)-1]...)×interval(mid-dx/2,mid+dx/2)
    else
        split_domain = ProductDomain([interval(a[i],b[i]) for i in 1:dim-1]...)×interval(mid-dx/2,mid+dx/2)×ProductDomain([interval(a[i],b[i]) for i in dim+1:length(a)]...)
    end
    split_domain
end

split_domain_Nd(gamma, mids::Vector{ELT}) where {ELT<:Real} = UnionDomain([split_domain_1d(gamma, d, mids[d]) for d in 1:length(mids)]...)

split_domain_Nd(gamma, mids::Vector{NTuple{N,ELT}}) where {N,ELT<:Real} = UnionDomain([       UnionDomain([split_domain_1d(gamma, d, x) for x in mids[d]]...) for d in 1:length(mids)]...)


function divide_and_conquer_N_util_operators(fplatform::BasisFunctions.Platform, i, ranges; options...)
    platform = fplatform.super_platform
    basis = primal(platform, i)
    S = sampler(fplatform, i)
    s = sampler(platform, i)
    domain = FrameFun.domain(primal(fplatform, i))
    gamma = BasisFunctions.grid(s)
    omega = BasisFunctions.grid(S)
    FrameFun.divide_and_conquer_N_util_operators(omega, gamma, basis, domain, ranges; options...)
end


# The grid on the boundary of omega
divide_and_conquer_N_util_operators(omega::AbstractGrid, gamma::AbstractGrid, basis::Dictionary, domain::Domains.Domain, ranges::AbstractVector; options...) =
    divide_and_conquer_N_util_operators(omega::AbstractGrid, gamma::AbstractGrid, basis::Dictionary, FrameFun.boundary_support_grid(basis, boundary_grid(gamma, domain), omega)::AbstractGrid, ranges::AbstractVector; options...)

function FrameFun.divide_and_conquer_N_util_operators(omega::AbstractGrid, gamma::AbstractGrid, basis::Dictionary, DMZ::AbstractGrid, ranges::AbstractVector; recur=nothing, options...)
    DMZs, GRs = FrameFun.split_DMZ(DMZ, basis, gamma, ranges; options...)
    OP = GridSamplingOperator(gridspace(gamma))*basis
    ops = []
    (recur == nothing) && (recur = dimension(basis)==3 ? 2 :1 )
    if recur >= 1
        ops = push!(ops, FrameFun.util_operators(OP, DMZ, DMZs, GRs, gamma, basis)...)
    end
    if recur >= 2
        @assert length(DMZs)==1
        try
        SplitDMZs, SplitGRs = FrameFun.split_DMZ(DMZs[1], basis,  gamma, ranges; shift=true)
        ops = push!(ops, FrameFun.util_operators(OP, DMZs[1], SplitDMZs, SplitGRs, gamma, basis)...)
        catch
            warn("no extra splitting was possible")
        end
    end
    tuple(ops...)
end


function util_operators(OP, DMZ, DMZs, GRs, gamma, basis)
    A0 = Array{CompositeOperator}(length(DMZs))
    GR0 = Array{IndexRestrictionOperator}(length(DMZs))
    # FE0 = Array{IndexExtensionOperator}(length(DMZs))
    for (i,dmz) in enumerate(DMZs)
        frame_restriction, grid_restriction = FrameFun._spline_util_restriction_operators(basis, gamma, dmz)
        A0[i] = grid_restriction*OP*frame_restriction'
        GR0[i] = grid_restriction
        # FE0[i] = frame_restriction'
    end

    A1 = Vector{CompositeOperator}(length(GRs))
    GR1 = Vector{IndexRestrictionOperator}(length(GRs))
    # FE1 = Vector{IndexExtensionOperator}(length(GRs))
    for (i,gr) in enumerate(GRs)
        dmz = FrameFun.boundary_support_grid(basis, gr, DMZ)
        frame_restriction, grid_restriction = FrameFun._spline_util_restriction_operators(basis, gamma, gr, dmz)
        A1[i] = grid_restriction*OP*frame_restriction'
        GR1[i] = grid_restriction
        # FE1[i] = frame_restriction'
    end
    A0, GR0, A1, GR1
end
